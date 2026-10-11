#!/usr/bin/env python3
"""Materialize the accepted six-study cluster annotation input contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import tarfile
from collections import Counter, defaultdict
from pathlib import Path


ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
HARVEST_ROOT = (
    PROJECT_ROOT
    / "outputs/multica-dispatch/2026-08-29-six-study-annotation-v2/harvest"
)
R9_ROOT = HARVEST_ROOT / "production_r9"
WOT_R10B_ARCHIVE = HARVEST_ROOT / "DOC-629_R10B/R10B_WOT_RELEASE.tar.gz"

STUDIES = {
    "GSE103221": ("Mus_musculus", "ANN05_gse103221_r9", "R9"),
    "GSE118258": ("Homo_sapiens", "ANN06_gse118258_r9", "R9"),
    "GSE150311": ("Homo_sapiens", "ANN07_gse150311_r9", "R9"),
    "GSE171127": ("Mus_musculus", "ANN08_gse171127_r9", "R9"),
    "GSE178325": ("Homo_sapiens", "ANN09_gse178325_r9", "R9"),
    "WOT_GSE115943": ("Mus_musculus", "ANN10_wot_r9", "R10B"),
}

EXPECTED_CLUSTER_COUNTS = {
    "GSE103221": 15,
    "GSE118258": 12,
    "GSE150311": 27,
    "GSE171127": 11,
    "GSE178325": 18,
    "WOT_GSE115943": 11,
}
EXPECTED_PRIMARY_CELLS = 606_186
EXPECTED_SENSITIVITY_CELLS = 76_435
EXPECTED_STRICT_ANCHORS = {
    ("GSE171127", "C006"),
    ("WOT_GSE115943", "C007"),
    ("WOT_GSE115943", "C009"),
    ("WOT_GSE115943", "C010"),
    ("WOT_GSE115943", "C011"),
}


def bool_text(value: str) -> str:
    value = value.strip().lower()
    if value not in {"true", "false"}:
        raise ValueError(f"Expected boolean text, got {value!r}")
    return value


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_tsv(data: bytes, source: str) -> tuple[list[str], list[dict[str, str]]]:
    text = data.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames:
        raise ValueError(f"Missing TSV header: {source}")
    rows = list(reader)
    return list(reader.fieldnames), rows


def read_regular_tsv(path: Path) -> tuple[bytes, list[str], list[dict[str, str]]]:
    data = path.read_bytes()
    header, rows = parse_tsv(data, str(path))
    return data, header, rows


def read_tar_tsv(
    archive_path: Path, member: str
) -> tuple[bytes, list[str], list[dict[str, str]]]:
    with tarfile.open(archive_path, "r:gz") as archive:
        extracted = archive.extractfile(member)
        if extracted is None:
            raise FileNotFoundError(f"Archive member missing: {member}")
        data = extracted.read()
    header, rows = parse_tsv(data, f"{archive_path}::{member}")
    return data, header, rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_policy() -> dict[tuple[str, str], dict[str, str]]:
    path = ANALYSIS_ROOT / "config/anchor_candidate_policy.tsv"
    _, _, rows = read_regular_tsv(path)
    result = {}
    for row in rows:
        key = (row["study_id"], row["local_cluster_id"])
        if key in result:
            raise ValueError(f"Duplicate policy key: {key}")
        result[key] = row
    return result


def classify_anchor(
    master: dict[str, str], policy: dict[str, str] | None
) -> tuple[str, str, str, str]:
    hold = master["manual_hold"] == "true"
    context = master["library_context_triggered"] == "true"
    strict_identity = (
        master["hard_identity"] not in {"", "NA"}
        and master["public_state_eligible"] == "true"
        and master["cross_study_vote_allowed"] == "true"
        and not hold
        and not context
    )
    role = policy["nominated_role"] if policy else "not_nominated"

    if strict_identity:
        return "A_confirmed", "accepted_anchor", "true", "accepted_hard_identity"
    if hold:
        return "C_held", "excluded_hold", "false", "manual_hold"
    if role == "reference_only":
        return "R_reference_only", "reference_only", "false", "external_or_endpoint_reference"
    if role == "validation_only":
        return "V_validation_only", "validation_only", "false", "validation_pool_not_seed"
    if context:
        return "C_context", "excluded_context", "false", "library_context_triggered"
    if role == "candidate":
        return "B_review_required", "review_required", "false", "candidate_requires_adjudication"
    if role == "local_only":
        return "C_local_only", "local_only", "false", "local_semantics_not_transferable"
    return "C_not_nominated", "not_nominated", "false", "no_anchor_nomination"


def materialize(output_dir: Path) -> list[tuple[str, bool, str]]:
    policy = load_policy()
    input_manifest = []
    all_master = []

    for expected_study, (species, release_dir, authority_release) in STUDIES.items():
        source_specs = []
        if authority_release == "R10B":
            for table_name in ("cluster_annotation.tsv", "cluster_working_annotation.tsv"):
                member = f"R10B_WOT_RELEASE/release/{table_name}"
                data, header, rows = read_tar_tsv(WOT_R10B_ARCHIVE, member)
                source_specs.append((table_name, f"{WOT_R10B_ARCHIVE}::{member}", data, header, rows))
        else:
            release = R9_ROOT / release_dir / "release"
            for table_name in ("cluster_annotation.tsv", "cluster_working_annotation.tsv"):
                path = release / table_name
                data, header, rows = read_regular_tsv(path)
                source_specs.append((table_name, str(path), data, header, rows))

        source_tables = {}
        for table_name, locator, data, header, rows in source_specs:
            source_tables[table_name] = rows
            input_manifest.append(
                {
                    "study_id": expected_study,
                    "authority_release": authority_release,
                    "table_name": table_name,
                    "source_locator": locator,
                    "sha256": sha256_bytes(data),
                    "row_count": len(rows),
                    "column_count": len(header),
                }
            )

        ann_rows = {
            (row["study_id"], row["local_cluster_id"]): row
            for row in source_tables["cluster_annotation.tsv"]
        }
        working_rows = {
            (row["study_id"], row["local_state_id"]): row
            for row in source_tables["cluster_working_annotation.tsv"]
        }
        if set(ann_rows) != set(working_rows):
            raise ValueError(f"Annotation/working key mismatch for {expected_study}")
        if {key[0] for key in ann_rows} != {expected_study}:
            raise ValueError(f"Unexpected study ID in authority input for {expected_study}")

        for key in sorted(ann_rows):
            ann = ann_rows[key]
            work = working_rows[key]
            for field in ("manual_hold", "public_state_eligible", "cross_study_vote_allowed"):
                if bool_text(ann[field]) != bool_text(work[field]):
                    raise ValueError(f"Authority disagreement for {key} field {field}")
            if ann["hard_identity"] != work["hard_identity"]:
                raise ValueError(f"Hard identity disagreement for {key}")

            all_master.append(
                {
                    "authority_release": authority_release,
                    "species": species,
                    "study_id": ann["study_id"],
                    "scope_id": ann["scope_id"],
                    "local_cluster_id": ann["local_cluster_id"],
                    "n_primary_cells": ann["n_primary_cells"],
                    "n_sensitivity_cells": ann["n_sensitivity_cells"],
                    "partition_stability": ann["partition_stability"],
                    "local_annotation_label": ann["annotation_label"],
                    "local_working_state": work["working_semantic_annotation"],
                    "confirmatory_annotation": work["confirmatory_annotation"],
                    "hard_identity": ann["hard_identity"],
                    "transition_mode": ann["transition_mode"],
                    "confidence_class": ann["confidence_class"],
                    "annotation_confidence": ann["annotation_confidence"],
                    "abstention_reason": ann["abstention_reason"],
                    "active_identity_programs": ann["active_identity_programs"],
                    "active_process_programs": ann["active_process_programs"],
                    "contributing_programs": work["contributing_programs"],
                    "strongest_family": work["strongest_family"],
                    "strongest_family_evidence": work["strongest_family_evidence"],
                    "second_strongest_family": work["second_strongest_family"],
                    "second_strongest_family_evidence": work["second_strongest_family_evidence"],
                    "family_margin": work["family_margin"],
                    "route_relation": work["route_relation"],
                    "library_context_triggered": bool_text(work["library_context_triggered"]),
                    "dominant_library_id": work["dominant_library_id"],
                    "dominant_library_fraction": work["dominant_library_fraction"],
                    "observed_day_or_stage": work["observed_day_or_stage"],
                    "observed_route": work["observed_route"],
                    "observed_donor": work["observed_donor"],
                    "observed_capture_or_well": work["observed_capture_or_well"],
                    "observed_context_summary": work["observed_context_summary"],
                    "biological_restriction_id": ann["biological_restriction_id"],
                    "manual_hold": bool_text(work["manual_hold"]),
                    "public_state_eligible": bool_text(work["public_state_eligible"]),
                    "cross_study_vote_allowed": bool_text(work["cross_study_vote_allowed"]),
                    "identity_evidence_replication_status": work[
                        "identity_evidence_replication_status"
                    ],
                    "single_library_identity_evidence_only": bool_text(
                        work["single_library_identity_evidence_only"]
                    ),
                    "prohibited_claim": work["prohibited_claim"],
                }
            )

    all_master.sort(key=lambda row: (row["study_id"], row["local_cluster_id"]))
    master_fields = list(all_master[0])
    write_tsv(output_dir / "cluster_annotation_master.tsv", master_fields, all_master)

    restriction_rows = []
    anchor_rows = []
    for row in all_master:
        key = (row["study_id"], row["local_cluster_id"])
        candidate = policy.get(key)
        tier, status, eligible, tier_reason = classify_anchor(row, candidate)
        restrictions = []
        if row["manual_hold"] == "true":
            restrictions.append("manual_hold")
        if row["library_context_triggered"] == "true":
            restrictions.append("library_context")
        if row["public_state_eligible"] == "false":
            restrictions.append("public_state_ineligible")
        if row["cross_study_vote_allowed"] == "false":
            restrictions.append("cross_study_vote_ineligible")
        if row["partition_stability"] == "rare_or_unresolved":
            restrictions.append("rare_or_unresolved_partition")
        if row["confidence_class"] in {"unresolved", "mixed_identity", "process_dominant"}:
            restrictions.append(row["confidence_class"])

        proposed_parent = candidate["proposed_parent_family"] if candidate else "unassigned"
        proposed_subtype = candidate["proposed_subtype"] if candidate else "NA"
        secondary_parent = candidate["secondary_parent_family"] if candidate else "NA"
        policy_note = candidate["policy_note"] if candidate else "not_nominated"
        required_evidence = (
            candidate["required_next_evidence"] if candidate else "No action unless nominated"
        )

        anchor_rows.append(
            {
                "study_id": row["study_id"],
                "species": row["species"],
                "local_cluster_id": row["local_cluster_id"],
                "n_primary_cells": row["n_primary_cells"],
                "local_working_state": row["local_working_state"],
                "hard_identity": row["hard_identity"],
                "proposed_parent_family": proposed_parent,
                "proposed_subtype": proposed_subtype,
                "secondary_parent_family": secondary_parent,
                "anchor_tier_initial": tier,
                "adjudication_status": status,
                "alignment_anchor_eligible": eligible,
                "tier_reason": tier_reason,
                "policy_note": policy_note,
                "required_next_evidence": required_evidence,
                "manual_hold": row["manual_hold"],
                "library_context_triggered": row["library_context_triggered"],
                "public_state_eligible": row["public_state_eligible"],
                "cross_study_vote_allowed": row["cross_study_vote_allowed"],
                "observed_day_or_stage": row["observed_day_or_stage"],
                "observed_route": row["observed_route"],
            }
        )
        restriction_rows.append(
            {
                "study_id": row["study_id"],
                "local_cluster_id": row["local_cluster_id"],
                "restriction_classes": ";".join(restrictions) if restrictions else "none",
                "biological_restriction_id": row["biological_restriction_id"],
                "manual_hold": row["manual_hold"],
                "library_context_triggered": row["library_context_triggered"],
                "public_state_eligible": row["public_state_eligible"],
                "cross_study_vote_allowed": row["cross_study_vote_allowed"],
                "anchor_tier_initial": tier,
                "alignment_anchor_eligible": eligible,
                "prohibited_claim": row["prohibited_claim"],
            }
        )

    write_tsv(output_dir / "anchor_adjudication.tsv", list(anchor_rows[0]), anchor_rows)
    write_tsv(
        output_dir / "restriction_registry.tsv", list(restriction_rows[0]), restriction_rows
    )

    state_registry_source = ANALYSIS_ROOT / "config/state_family_registry.tsv"
    state_data, state_header, state_rows = read_regular_tsv(state_registry_source)
    write_tsv(output_dir / "state_family_registry.tsv", state_header, state_rows)
    input_manifest.extend(
        [
            {
                "study_id": "ALL",
                "authority_release": "PLAN_V1",
                "table_name": "anchor_candidate_policy.tsv",
                "source_locator": str(ANALYSIS_ROOT / "config/anchor_candidate_policy.tsv"),
                "sha256": sha256_file(ANALYSIS_ROOT / "config/anchor_candidate_policy.tsv"),
                "row_count": len(policy),
                "column_count": 8,
            },
            {
                "study_id": "ALL",
                "authority_release": "PLAN_V1",
                "table_name": "state_family_registry.tsv",
                "source_locator": str(state_registry_source),
                "sha256": sha256_bytes(state_data),
                "row_count": len(state_rows),
                "column_count": len(state_header),
            },
        ]
    )
    write_tsv(output_dir / "input_manifest.tsv", list(input_manifest[0]), input_manifest)

    study_summary = []
    by_study = defaultdict(list)
    for row, anchor in zip(all_master, anchor_rows):
        by_study[row["study_id"]].append((row, anchor))
    for study_id in sorted(by_study):
        pairs = by_study[study_id]
        tiers = Counter(anchor["anchor_tier_initial"] for _, anchor in pairs)
        study_summary.append(
            {
                "study_id": study_id,
                "species": pairs[0][0]["species"],
                "n_clusters": len(pairs),
                "n_primary_cells": sum(int(row["n_primary_cells"]) for row, _ in pairs),
                "n_sensitivity_cells": sum(
                    int(row["n_sensitivity_cells"]) for row, _ in pairs
                ),
                "n_A_confirmed": tiers["A_confirmed"],
                "n_B_review_required": tiers["B_review_required"],
                "n_C_or_held": sum(
                    count for tier, count in tiers.items() if tier.startswith("C_")
                ),
                "n_validation_only": tiers["V_validation_only"],
                "n_reference_only": tiers["R_reference_only"],
            }
        )
    write_tsv(output_dir / "study_annotation_summary.tsv", list(study_summary[0]), study_summary)

    observed_cluster_counts = Counter(row["study_id"] for row in all_master)
    strict_anchor_keys = {
        (row["study_id"], row["local_cluster_id"])
        for row in anchor_rows
        if row["anchor_tier_initial"] == "A_confirmed"
    }
    wot_c001 = next(
        row
        for row in all_master
        if row["study_id"] == "WOT_GSE115943" and row["local_cluster_id"] == "C001"
    )
    checks = [
        ("cluster_total", len(all_master) == 94, f"observed={len(all_master)} expected=94"),
        (
            "cluster_counts_by_study",
            dict(observed_cluster_counts) == EXPECTED_CLUSTER_COUNTS,
            f"observed={dict(observed_cluster_counts)} expected={EXPECTED_CLUSTER_COUNTS}",
        ),
        (
            "primary_cell_closure",
            sum(int(row["n_primary_cells"]) for row in all_master)
            == EXPECTED_PRIMARY_CELLS,
            f"observed={sum(int(row['n_primary_cells']) for row in all_master)} expected={EXPECTED_PRIMARY_CELLS}",
        ),
        (
            "sensitivity_cell_closure",
            sum(int(row["n_sensitivity_cells"]) for row in all_master)
            == EXPECTED_SENSITIVITY_CELLS,
            f"observed={sum(int(row['n_sensitivity_cells']) for row in all_master)} expected={EXPECTED_SENSITIVITY_CELLS}",
        ),
        (
            "strict_anchor_set",
            strict_anchor_keys == EXPECTED_STRICT_ANCHORS,
            f"observed={sorted(strict_anchor_keys)} expected={sorted(EXPECTED_STRICT_ANCHORS)}",
        ),
        (
            "wot_c001_r10b_semantics",
            wot_c001["authority_release"] == "R10B"
            and wot_c001["local_annotation_label"] == "early_source_context"
            and wot_c001["local_working_state"] == "source_like"
            and wot_c001["transition_mode"] == "stable_source"
            and wot_c001["route_relation"] == "not_applicable"
            and wot_c001["hard_identity"] == "NA"
            and wot_c001["manual_hold"] == "true"
            and wot_c001["public_state_eligible"] == "false"
            and wot_c001["cross_study_vote_allowed"] == "false",
            "R10B source-like hold fields checked",
        ),
        (
            "no_eligible_anchor_with_hold_or_context",
            all(
                anchor["manual_hold"] == "false"
                and anchor["library_context_triggered"] == "false"
                for anchor in anchor_rows
                if anchor["alignment_anchor_eligible"] == "true"
            ),
            "eligible anchors must be free of hold and library-context gates",
        ),
        (
            "policy_keys_exist",
            set(policy).issubset(
                {(row["study_id"], row["local_cluster_id"]) for row in all_master}
            ),
            "every candidate policy key must exist in authority inputs",
        ),
        (
            "state_family_count",
            len(state_rows) == 5,
            f"observed={len(state_rows)} expected=5",
        ),
    ]
    validation_rows = [
        {"check_id": check_id, "status": "PASS" if passed else "FAIL", "details": details}
        for check_id, passed, details in checks
    ]
    write_tsv(output_dir / "validation_summary.tsv", list(validation_rows[0]), validation_rows)

    tier_counts = Counter(row["anchor_tier_initial"] for row in anchor_rows)
    report = [
        "# Unified annotation input materialization report",
        "",
        "Decision: `PASS`" if all(passed for _, passed, _ in checks) else "Decision: `FAIL`",
        "",
        "## Scope",
        "",
        "- Authority: five R9 releases plus WOT R10B.",
        f"- Frozen clusters: {len(all_master)}.",
        f"- Primary cells represented: {sum(int(row['n_primary_cells']) for row in all_master):,}.",
        f"- Sensitivity cells represented: {sum(int(row['n_sensitivity_cells']) for row in all_master):,}.",
        "- No cluster membership or authority annotation field was modified.",
        "",
        "## Initial anchor tiers",
        "",
    ]
    for tier in sorted(tier_counts):
        report.append(f"- `{tier}`: {tier_counts[tier]} clusters.")
    report.extend(
        [
            "",
            "Only A_confirmed clusters are currently alignment-anchor eligible. B_review_required clusters require independent adjudication; C, V and R rows cannot train the first model.",
            "",
            "## Stop condition",
            "",
            "Do not train scVI/scANVI until mouse and human each have an accepted minimum anchor set spanning source, transition and target-like families. This release does not authorize integration, trajectory, APA or cross-species alignment.",
            "",
        ]
    )
    (output_dir / "report.md").write_text("\n".join(report), encoding="utf-8")

    manifest_rows = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "release_manifest.tsv":
            manifest_rows.append(
                {
                    "relative_path": path.name,
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
            )
    write_tsv(output_dir / "release_manifest.tsv", list(manifest_rows[0]), manifest_rows)

    if not all(passed for _, passed, _ in checks):
        failed = [check_id for check_id, passed, _ in checks if not passed]
        raise RuntimeError(f"Validation failed: {failed}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ANALYSIS_ROOT / "release")
    args = parser.parse_args()
    materialize(args.output.resolve())


if __name__ == "__main__":
    main()

