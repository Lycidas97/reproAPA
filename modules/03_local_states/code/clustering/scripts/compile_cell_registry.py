#!/usr/bin/env python3
"""Compile a lossless primary and sensitivity clustering registry."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from cluster_utils import write_json, write_marker


IDENTITY_COLUMNS = [
    "study_id",
    "scope_id",
    "cell_id",
    "library_id",
    "barcode",
    "species",
    "sample_accession",
    "sample_title",
    "experimental_time",
    "condition",
]


def require_unique(frame: pd.DataFrame, label: str) -> None:
    if frame["cell_id"].isna().any() or frame["cell_id"].duplicated().any():
        raise ValueError(f"{label}: missing or duplicate cell_id")


def compile_dataset(output_root: Path, dataset: str) -> tuple[pd.DataFrame, dict]:
    scope = output_root / dataset / "scopes" / "studywide"
    release = json.loads((scope / "release/scope_release.json").read_text(encoding="utf-8"))
    if not release.get("all_input_cells_retained") or not release.get("all_sensitivity_cells_retained"):
        raise ValueError(f"{dataset}: scope release does not retain all primary and sensitivity cells")

    primary_cells = pd.read_csv(scope / "prepared/cells.tsv.gz", sep="\t", dtype=str)
    primary_clusters = pd.read_csv(
        scope / "frozen/cell_clusters.tsv.gz", sep="\t", dtype=str,
        usecols=["cell_id", "cluster_id", "selected_resolution"],
    )
    program_cells = pd.read_csv(scope / "program_audit/cell_program_scores.tsv.gz", sep="\t", dtype=str)
    cluster_transfer = pd.read_csv(scope / "full_partition_audit/cluster_transfer_audit.tsv", sep="\t", dtype=str)
    cluster_programs = pd.read_csv(scope / "program_audit/cluster_program_scores.tsv", sep="\t", dtype=str)
    for label, frame in {
        "primary cells": primary_cells,
        "primary clusters": primary_clusters,
        "primary program scores": program_cells,
    }.items():
        require_unique(frame, f"{dataset} {label}")

    cluster_attributes = cluster_transfer[[
        "cluster_id", "cluster_release_class", "review_reasons"
    ]].rename(columns={"review_reasons": "cluster_stability_review_reasons"}).merge(
        cluster_programs[["cluster_id", "requires_program_review", "program_review_reasons"]],
        on="cluster_id", how="outer", validate="one_to_one",
    )
    if cluster_attributes.isna().any().any():
        raise ValueError(f"{dataset}: incomplete cluster stability or program attributes")

    primary = primary_cells[IDENTITY_COLUMNS].merge(
        primary_clusters, on="cell_id", how="left", validate="one_to_one"
    ).merge(
        program_cells[[
            "cell_id", "s_phase", "g2m_phase", "immediate_early_stress",
            "interferon_response", "cell_cycle_phase",
        ]],
        on="cell_id", how="left", validate="one_to_one",
    ).merge(cluster_attributes, on="cluster_id", how="left", validate="many_to_one")
    primary["qc_tier"] = "primary_pass"
    primary["cluster_assignment_probability"] = "1.0"
    primary["normalized_cluster_entropy"] = "0.0"
    primary["projection_status"] = "direct_fit"

    sensitivity_cells = pd.read_csv(scope / "sensitivity/prepared/cells.tsv.gz", sep="\t", dtype=str)
    sensitivity_projection = pd.read_csv(scope / "sensitivity/projection/sensitivity_projection.tsv.gz", sep="\t", dtype=str)
    require_unique(sensitivity_cells, f"{dataset} sensitivity cells")
    require_unique(sensitivity_projection, f"{dataset} sensitivity projection")
    sensitivity = sensitivity_cells[IDENTITY_COLUMNS].merge(
        sensitivity_projection[[
            "cell_id", "assigned_cluster", "maximum_cluster_probability",
            "normalized_cluster_entropy", "projection_status",
        ]],
        on="cell_id", how="left", validate="one_to_one",
    ).rename(columns={
        "assigned_cluster": "cluster_id",
        "maximum_cluster_probability": "cluster_assignment_probability",
    }).merge(cluster_attributes, on="cluster_id", how="left", validate="many_to_one")
    sensitivity["selected_resolution"] = str(release["selected_resolution"])
    sensitivity["qc_tier"] = "sensitivity_only"
    for column in ["s_phase", "g2m_phase", "immediate_early_stress", "interferon_response", "cell_cycle_phase"]:
        sensitivity[column] = pd.NA

    ordered = IDENTITY_COLUMNS + [
        "qc_tier", "cluster_id", "selected_resolution", "cluster_assignment_probability",
        "normalized_cluster_entropy", "projection_status", "cluster_release_class",
        "cluster_stability_review_reasons", "requires_program_review", "program_review_reasons",
        "s_phase", "g2m_phase", "immediate_early_stress", "interferon_response", "cell_cycle_phase",
    ]
    combined = pd.concat([primary[ordered], sensitivity[ordered]], ignore_index=True)
    combined.insert(0, "dataset", dataset)
    require_unique(combined, f"{dataset} combined registry")
    if combined[["cluster_id", "cluster_release_class", "projection_status"]].isna().any().any():
        raise ValueError(f"{dataset}: unassigned retained cell or missing cluster release class")
    if len(primary) != int(release["n_cells"]):
        raise ValueError(f"{dataset}: primary count does not match scope release")
    if len(sensitivity) != int(release["n_sensitivity_cells_projected"]):
        raise ValueError(f"{dataset}: sensitivity count does not match scope release")
    summary = {
        "dataset": dataset,
        "study_id": str(release["study_id"]),
        "primary_cells": len(primary),
        "sensitivity_cells": len(sensitivity),
        "retained_cells": len(combined),
        "clusters": int(combined["cluster_id"].nunique()),
        "stable_candidate_cells": int(combined["cluster_release_class"].eq("stable_candidate").sum()),
        "rare_or_unresolved_cells": int(combined["cluster_release_class"].eq("rare_or_unresolved").sum()),
    }
    return combined, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--datasets", nargs="+", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output_root = Path(args.output_root)
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)

    frames = []
    summaries = []
    for dataset in sorted(args.datasets):
        frame, summary = compile_dataset(output_root, dataset)
        frames.append(frame)
        summaries.append(summary)
    registry = pd.concat(frames, ignore_index=True)
    require_unique(registry, "project registry")
    registry.to_csv(
        target / "clustering_cell_registry.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(summaries).to_csv(target / "dataset_summary.tsv", sep="\t", index=False)
    summary = {
        "status": "PASS",
        "datasets": len(summaries),
        "primary_cells": int(sum(row["primary_cells"] for row in summaries)),
        "sensitivity_cells": int(sum(row["sensitivity_cells"] for row in summaries)),
        "retained_cells": len(registry),
        "duplicate_cell_ids": int(registry["cell_id"].duplicated().sum()),
        "technical_exclude_cells_in_registry": int(registry["qc_tier"].eq("technical_exclude").sum()),
        "all_retained_cells_assigned": bool(registry["cluster_id"].notna().all()),
    }
    write_json(target / "summary.json", summary)
    write_marker(target / ".ok", **summary)


if __name__ == "__main__":
    main()
