#!/usr/bin/env python3
"""Audit the six-study QC release after Snakemake production finishes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


TIERS = {
    "primary_pass": "primary_membership.tsv.gz",
    "sensitivity_only": "sensitivity_membership.tsv.gz",
    "technical_exclude": "technical_membership.tsv.gz",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def membership_ids(path: Path) -> set[str]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        values = [row["global_cell_id"] for row in rows]
    if len(values) != len(set(values)):
        raise ValueError(f"duplicate global_cell_id in {path}")
    return set(values)


def set_hash(values: set[str]) -> str:
    return hashlib.sha256(("\n".join(sorted(values)) + "\n").encode()).hexdigest()


def add_check(checks: list[dict[str, str]], name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def audit_study(root: Path, study: dict[str, str], review_status: str) -> dict:
    run = root / "runs" / study["output_dir"]
    expected_libraries = int(study["expected_libraries"])
    expected_cells = int(study["expected_called_cells"])
    checks: list[dict[str, str]] = []
    required = [
        run / "run_complete.ok",
        run / "independent_validation.tsv",
        run / "tier_counts_and_hashes.tsv",
        run / "figures/qc_contact_sheets/contact_sheet_index.tsv",
        *[run / filename for filename in TIERS.values()],
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        return {
            "scope": study["scope"],
            "output_dir": study["output_dir"],
            "status": "PENDING",
            "checks": [{"check": "required_outputs", "status": "PENDING", "detail": ";".join(missing)}],
        }

    complete = json.loads((run / "run_complete.ok").read_text(encoding="utf-8"))
    add_check(checks, "run_complete_status", complete.get("status") == "PASS", str(complete.get("status")))
    add_check(checks, "library_count", int(complete.get("libraries", -1)) == expected_libraries,
              f"{complete.get('libraries')}/{expected_libraries}")
    add_check(checks, "called_cell_count", int(complete.get("called_cells", -1)) == expected_cells,
              f"{complete.get('called_cells')}/{expected_cells}")
    for key in ("scrublet_hard_exclusions", "scrublet_sample_gate_changes", "scrublet_tier_changes"):
        add_check(checks, key, int(complete.get(key, -1)) == 0, str(complete.get(key)))

    statuses = sorted(run.glob("per_library/*/status.json"))
    parsed_statuses = [json.loads(path.read_text(encoding="utf-8")) for path in statuses]
    add_check(checks, "per_library_status_count", len(statuses) == expected_libraries,
              f"{len(statuses)}/{expected_libraries}")
    add_check(checks, "per_library_all_completed",
              len(statuses) == expected_libraries and all(row.get("status") == "completed" for row in parsed_statuses),
              f"completed={sum(row.get('status') == 'completed' for row in parsed_statuses)}")

    validation = read_tsv(run / "independent_validation.tsv")
    add_check(checks, "independent_validation",
              bool(validation) and all(row.get("status") == "PASS" for row in validation),
              f"pass={sum(row.get('status') == 'PASS' for row in validation)}/{len(validation)}")

    tier_rows = read_tsv(run / "tier_counts_and_hashes.tsv")
    tier_lookup = {row["tier"]: row for row in tier_rows}
    add_check(checks, "tier_inventory", set(tier_lookup) == set(TIERS), ",".join(sorted(tier_lookup)))
    observed_sets: dict[str, set[str]] = {}
    for tier, filename in TIERS.items():
        path = run / filename
        values = membership_ids(path)
        observed_sets[tier] = values
        expected = tier_lookup.get(tier, {})
        count_ok = expected and len(values) == int(expected["n_cells"])
        hash_ok = expected and set_hash(values) == expected["membership_sha256"]
        add_check(checks, f"membership_count:{tier}", bool(count_ok),
                  f"{len(values)}/{expected.get('n_cells', 'missing')}")
        add_check(checks, f"membership_hash:{tier}", bool(hash_ok),
                  f"{set_hash(values)}/{expected.get('membership_sha256', 'missing')}")
    union = set().union(*observed_sets.values())
    pairwise_disjoint = all(
        observed_sets[left].isdisjoint(observed_sets[right])
        for index, left in enumerate(TIERS)
        for right in list(TIERS)[index + 1:]
    )
    add_check(checks, "tier_disjointness", pairwise_disjoint, str(pairwise_disjoint))
    add_check(checks, "tier_union_called_cells", len(union) == expected_cells, f"{len(union)}/{expected_cells}")

    contact_rows = read_tsv(run / "figures/qc_contact_sheets/contact_sheet_index.tsv")
    contact_count = sum(int(row["n_libraries"]) for row in contact_rows)
    contact_paths = [
        path if path.is_absolute() else root / path
        for path in (Path(row["path"]) for row in contact_rows)
    ]
    contact_paths_exist = bool(contact_rows) and all(path.is_file() for path in contact_paths)
    add_check(checks, "contact_sheet_library_count", contact_count == expected_libraries,
              f"{contact_count}/{expected_libraries}")
    add_check(checks, "contact_sheet_paths", contact_paths_exist, f"pages={len(contact_rows)}")
    add_check(checks, "manual_visual_review", review_status == "PASS", review_status)

    failed = [row for row in checks if row["status"] == "FAIL"]
    only_visual = failed and all(row["check"] == "manual_visual_review" for row in failed)
    status = "VISUAL_PENDING" if only_visual else "FAIL" if failed else "PASS"
    return {"scope": study["scope"], "output_dir": study["output_dir"], "status": status, "checks": checks}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    studies = read_tsv(root / "production/studies.tsv")
    reviews = {row["scope"]: row["status"] for row in read_tsv(root / "production/manual_visual_review.tsv")}
    results = [audit_study(root, row, reviews.get(row["scope"], "MISSING")) for row in studies]
    statuses = [row["status"] for row in results]
    if statuses and all(value == "PASS" for value in statuses):
        overall = "PASS"
    elif any(value == "FAIL" for value in statuses):
        overall = "FAIL"
    elif any(value == "PENDING" for value in statuses):
        overall = "PENDING"
    elif statuses and all(value in {"PASS", "VISUAL_PENDING"} for value in statuses):
        overall = "VISUAL_PENDING"
    else:
        overall = "FAIL"
    report = {
        "status": overall,
        "n_studies": len(studies),
        "expected_libraries": sum(int(row["expected_libraries"]) for row in studies),
        "expected_called_cells": sum(int(row["expected_called_cells"]) for row in studies),
        "studies": results,
    }
    state = root / "runs/_production_acceptance"
    state.mkdir(parents=True, exist_ok=True)
    (state / "acceptance.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (state / "acceptance.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scope", "output_dir", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in results)
    print(json.dumps({key: report[key] for key in ("status", "n_studies", "expected_libraries", "expected_called_cells")}, sort_keys=True))
    if args.strict and overall != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
