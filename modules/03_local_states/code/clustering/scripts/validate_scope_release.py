#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from cluster_utils import ordered_sha256, write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope-root", required=True)
    args = parser.parse_args()
    root = Path(args.scope_root)
    required = [
        root / "selection" / "selection.tsv", root / "selection" / "resolution_selection.png",
        root / "tree" / "resolution_tree.png", root / "tree" / "resolution_tree_selected_window.png",
        root / "frozen" / "cell_clusters.tsv.gz",
        root / "markers" / "top_markers.tsv", root / "technical" / "summary.json",
        root / "classifier" / "summary.json", root / "classifier" / "classifier_confusion.png",
        root / "visual" / "partition_umap.png", root / "visual" / "umap_coordinates.tsv.gz",
        root / "full_partition_audit" / "summary.json",
        root / "full_partition_audit" / "cluster_transfer_audit.tsv",
        root / "full_partition_audit" / "evaluation_to_full_transfer.png",
        root / "program_audit" / "summary.json", root / "program_audit" / "cluster_program_scores.tsv",
        root / "program_audit" / "program_overlay.png",
        root / "sensitivity" / "prepared" / "cells.tsv.gz",
        root / "sensitivity" / "projection" / "summary.json",
        root / "sensitivity" / "projection" / "sensitivity_projection.tsv.gz",
        root / "sensitivity" / "projection" / "sensitivity_projection.png",
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise ValueError(f"missing or empty scope release artifacts: {missing}")
    cells = pd.read_csv(root / "prepared" / "cells.tsv.gz", sep="\t", dtype=str)
    partition = pd.read_csv(root / "frozen" / "cell_clusters.tsv.gz", sep="\t", dtype=str)
    if len(cells) != len(partition) or partition["cell_id"].duplicated().any():
        raise ValueError("final partition does not retain each input cell exactly once")
    if ordered_sha256(cells["cell_id"]) != ordered_sha256(partition["cell_id"]):
        raise ValueError("final partition cell order drift")
    transfer = pd.read_csv(root / "full_partition_audit" / "cluster_transfer_audit.tsv", sep="\t", dtype=str)
    if set(transfer["cluster_id"]) != set(partition["cluster_id"]):
        raise ValueError("full-partition audit does not cover every frozen cluster")
    if transfer["full_n_cells"].astype(int).sum() != len(partition):
        raise ValueError("full-partition audit cluster sizes do not close")
    sensitivity_cells = pd.read_csv(root / "sensitivity" / "prepared" / "cells.tsv.gz", sep="\t", dtype=str)
    sensitivity_projection = pd.read_csv(
        root / "sensitivity" / "projection" / "sensitivity_projection.tsv.gz", sep="\t", dtype=str
    )
    if len(sensitivity_cells) != len(sensitivity_projection):
        raise ValueError("sensitivity projection does not retain every sensitivity-tier cell")
    if len(sensitivity_cells) and ordered_sha256(sensitivity_cells["cell_id"]) != ordered_sha256(sensitivity_projection["cell_id"]):
        raise ValueError("sensitivity projection cell order drift")
    decision = json.loads((root / "selection" / "selected_resolution.json").read_text(encoding="utf-8"))
    transfer_summary = json.loads((root / "full_partition_audit" / "summary.json").read_text(encoding="utf-8"))
    program_summary = json.loads((root / "program_audit" / "summary.json").read_text(encoding="utf-8"))
    sensitivity_summary = json.loads((root / "sensitivity" / "projection" / "summary.json").read_text(encoding="utf-8"))
    statuses = [decision["status"], transfer_summary["status"], program_summary["status"], sensitivity_summary["status"]]
    status = "PASS" if all(item == "PASS" for item in statuses) else "REVIEW"
    write_json(root / "release" / "scope_release.json", {"scope_id": cells["scope_id"].iloc[0],
                                                           "study_id": cells["study_id"].iloc[0],
                                                           "n_cells": len(cells),
                                                           "n_clusters": partition["cluster_id"].nunique(),
                                                           "n_stable_candidate_clusters": int(transfer["cluster_release_class"].eq("stable_candidate").sum()),
                                                           "n_rare_or_unresolved_clusters": int(transfer["cluster_release_class"].eq("rare_or_unresolved").sum()),
                                                           "n_sensitivity_cells_projected": len(sensitivity_projection),
                                                           "selected_resolution": decision["selected_resolution"],
                                                           "selection_status": decision["status"],
                                                           "evaluation_to_full_status": transfer_summary["status"],
                                                           "evaluation_to_full_adjusted_rand_index": transfer_summary["adjusted_rand_index"],
                                                           "program_overlay_status": program_summary["status"],
                                                           "n_program_review_clusters": program_summary["n_program_review_clusters"],
                                                           "program_review_clusters": program_summary["program_review_clusters"],
                                                           "sensitivity_projection_status": sensitivity_summary["status"],
                                                           "all_input_cells_retained": True,
                                                           "all_sensitivity_cells_retained": True,
                                                           "manual_biological_visual_review_required": True,
                                                           "status": status})
    write_marker(root / "release" / ".ok", status=status, n_cells=len(cells),
                 n_sensitivity_cells=len(sensitivity_projection))


if __name__ == "__main__":
    main()
