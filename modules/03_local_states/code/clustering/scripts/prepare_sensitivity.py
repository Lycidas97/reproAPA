#!/usr/bin/env python3
"""Prepare sensitivity-tier cells in the frozen primary feature universe."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from cluster_utils import ordered_sha256, write_marker
from prepare_scope import load_10x, load_h5ad, load_membership, normalize_log1p, set_sha256


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--scope-root", required=True)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    scope_cfg = config["scopes"][args.scope]
    sensitivity_cfg = copy.deepcopy(config)
    sensitivity_cfg["membership_path"] = config["sensitivity_membership_path"]
    membership = load_membership(sensitivity_cfg, args.scope, scope_cfg)
    expected_n = int(scope_cfg["sensitivity_n_cells"])
    if len(membership) != expected_n:
        raise ValueError("sensitivity membership count drift")
    if ordered_sha256(membership["cell_id"]) != scope_cfg["sensitivity_source_membership_order_sha256"]:
        raise ValueError("sensitivity membership order hash drift")
    if set_sha256(membership["cell_id"]) != scope_cfg["sensitivity_cell_set_sha256"]:
        raise ValueError("sensitivity membership cell-set hash drift")

    root = Path(args.scope_root)
    target = root / "sensitivity" / "prepared"
    target.mkdir(parents=True, exist_ok=True)
    primary_features = pd.read_csv(root / "prepared" / "features.tsv.gz", sep="\t", dtype=str)
    target_ids = primary_features["gene_id"].astype(str).tolist()
    chunks = []
    positions = []
    assigned = np.zeros(len(membership), dtype=np.int8)
    for source in scope_cfg["expression_sources"]:
        selector = source.get("selector", "source_scope_id")
        selected_positions = np.flatnonzero(
            membership[selector].astype(str).eq(str(source["value"])).to_numpy()
        )
        if not len(selected_positions):
            continue
        assigned[selected_positions] += 1
        selected_membership = membership.iloc[selected_positions]
        if source["mode"] == "10x_h5":
            matrix, features = load_10x(source["path"], selected_membership)
        elif source["mode"] == "h5ad":
            matrix, features = load_h5ad(source["path"], selected_membership)
        else:
            raise ValueError(f"unsupported expression source mode: {source['mode']}")
        indexer = pd.Index(features["gene_id"].astype(str)).get_indexer(target_ids)
        if np.any(indexer < 0):
            raise ValueError("primary feature universe is absent from a sensitivity expression source")
        chunks.append(matrix[:, indexer].tocsr())
        positions.append(selected_positions)
    if np.any(assigned != 1):
        raise ValueError("each sensitivity cell must map to exactly one expression source")
    if chunks:
        stacked = sparse.vstack(chunks, format="csr")
        stacked_positions = np.concatenate(positions)
        reorder = np.argsort(stacked_positions, kind="stable")
        if not np.array_equal(stacked_positions[reorder], np.arange(len(membership))):
            raise ValueError("sensitivity expression source row contract is incomplete")
        raw = stacked[reorder].tocsr()
        matrix = normalize_log1p(raw, float(scope_cfg["pipeline"]["normalization_target_sum"]))
    else:
        matrix = sparse.csr_matrix((0, len(primary_features)), dtype=np.float32)
    sparse.save_npz(target / "log1p_gex.npz", matrix, compressed=True)
    membership.to_csv(
        target / "cells.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    metadata = {
        "study_id": config["study_id"],
        "scope_id": args.scope,
        "n_cells": len(membership),
        "n_features": len(primary_features),
        "membership_cell_set_sha256": set_sha256(membership["cell_id"]),
        "primary_feature_universe_sha256": ordered_sha256(primary_features["gene_id"]),
        "used_to_fit_hvg_pca_neighbors_or_clusters": False,
        "status": "PASS",
    }
    (target / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_marker(target / ".ok", status="PASS", n_cells=len(membership), n_features=len(primary_features))


if __name__ == "__main__":
    main()
