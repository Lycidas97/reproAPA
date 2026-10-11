#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml
from scipy import sparse

from cluster_utils import leiden, ordered_sha256, write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-config", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--scope-root", required=True)
    args = parser.parse_args()
    study = yaml.safe_load(Path(args.study_config).read_text(encoding="utf-8"))
    project = yaml.safe_load(Path(args.project_config).read_text(encoding="utf-8"))
    pipe = study["scopes"][args.scope]["pipeline"]
    root = Path(args.scope_root)
    target = root / "frozen"
    target.mkdir(parents=True, exist_ok=True)
    decision = json.loads((root / "selection" / "selected_resolution.json").read_text(encoding="utf-8"))
    graph = sparse.load_npz(root / "embedding" / "knn_graph.npz")
    cells = pd.read_csv(root / "prepared" / "cells.tsv.gz", sep="\t", dtype=str)
    implementation = str(pipe["leiden_implementation"])
    labels = leiden(graph, float(decision["selected_resolution"]), int(project["selection"]["seed"]),
                    int(pipe["leiden_iterations"]), implementation)
    metadata_columns = [
        column for column in ["species", "sample_accession", "sample_title", "experimental_time", "condition"]
        if column in cells
    ]
    partition = cells[["study_id", "scope_id", "cell_id", "library_id", *metadata_columns]].copy()
    partition["cluster_id"] = [f"C{int(value) + 1:03d}" for value in labels]
    partition["selected_resolution"] = float(decision["selected_resolution"])
    partition.to_csv(target / "cell_clusters.tsv.gz", sep="\t", index=False,
                     compression={"method": "gzip", "mtime": 0})
    sizes = partition.groupby("cluster_id", sort=True).size().rename("n_cells").reset_index()
    sizes["fraction"] = sizes["n_cells"] / len(partition)
    sizes.to_csv(target / "cluster_sizes.tsv", sep="\t", index=False)
    write_json(target / "metadata.json", {"study_id": study["study_id"], "scope_id": args.scope,
                                           "n_cells": len(partition), "n_clusters": len(sizes),
                                           "selected_resolution": decision["selected_resolution"],
                                           "ordered_cell_id_sha256": ordered_sha256(partition["cell_id"]),
                                           "all_input_cells_retained": True})
    write_marker(target / ".ok", status="PASS", n_cells=len(partition), n_clusters=len(sizes))


if __name__ == "__main__":
    main()
