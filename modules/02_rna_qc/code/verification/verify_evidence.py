"""Reconcile the accepted seven-study RNA QC ledgers with manuscript Table S2."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re


TIERS = {"primary_pass", "sensitivity_only", "technical_exclude"}


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_run(run_dir, expected, gse_acceptance=None):
    summary = json.loads((run_dir / "run_summary.json").read_text())
    assert summary["status"] == "PASS", run_dir
    assert summary["libraries"] == expected["libraries"], run_dir
    assert summary["called_cells"] == expected["cells"], run_dir
    assert all(summary[key] == 0 for key in (
        "scrublet_hard_exclusions", "scrublet_sample_gate_changes", "scrublet_tier_changes"
    )), run_dir
    rows = read_tsv(run_dir / "tier_counts_and_hashes.tsv")
    assert {row["tier"] for row in rows} == TIERS and len(rows) == 3, run_dir
    counts = {}
    for row in rows:
        assert row["study_id"] == summary["study_id"], run_dir
        assert re.fullmatch(r"[0-9a-f]{64}", row["membership_sha256"]), run_dir
        count = int(row["n_cells"])
        assert count == summary["tier_counts"][row["tier"]], run_dir
        if gse_acceptance is not None:
            assert count == gse_acceptance["tiers"][row["tier"]]["n_cells"], run_dir
            assert row["membership_sha256"] == gse_acceptance["tiers"][row["tier"]]["cell_set_sha256"], run_dir
        counts[row["tier"]] = count
    assert sum(counts.values()) == expected["cells"], run_dir
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s2", required=True, type=Path)
    parser.add_argument("--evidence-dir", required=True, type=Path)
    parser.add_argument("--gse-acceptance", required=True, type=Path)
    args = parser.parse_args()

    inventory = read_tsv(args.s2)
    expected = {}
    for row in inventory:
        item = expected.setdefault(row["study"], {"libraries": 0, "cells": 0})
        item["libraries"] += 1
        item["cells"] += int(row["called_cells_before_downstream_qc"])
    assert len(inventory) == 210 and len(expected) == 7
    assert sum(item["cells"] for item in expected.values()) == 843893

    evidence = args.evidence_dir / "runs"
    accepted = json.loads((evidence / "_production_acceptance/acceptance.json").read_text())
    assert accepted["status"] == "PASS" and accepted["n_studies"] == 6
    assert accepted["expected_libraries"] == 193
    assert accepted["expected_called_cells"] == 778207
    totals = Counter()
    seen = set()
    for study in accepted["studies"]:
        assert study["status"] == "PASS"
        assert all(check["status"] == "PASS" for check in study["checks"])
        run_dir = evidence / study["output_dir"]
        summary = json.loads((run_dir / "run_summary.json").read_text())
        study_id = summary["study_id"]
        assert study_id in expected and study_id not in seen
        seen.add(study_id)
        totals.update(check_run(run_dir, expected[study_id]))
    assert len(seen) == 6 and "GSE221739" not in seen

    gse = json.loads(args.gse_acceptance.read_text())
    assert gse["status"] == "PASS" and gse["independent_all_pass"] is True
    assert gse["tiers_disjoint_and_complete"] is True
    assert gse["libraries"] == 17 and gse["called_cells"] == 65686
    totals.update(check_run(evidence / "gse221739_production", expected["GSE221739"], gse))
    assert sum(totals.values()) == 843893
    print("PASS: 7 studies, 210 libraries, 843893 called cells; " +
          ", ".join(f"{tier}={totals[tier]}" for tier in sorted(TIERS)))


if __name__ == "__main__":
    main()
