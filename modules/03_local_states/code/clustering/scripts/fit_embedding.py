#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from threadpoolctl import threadpool_limits

from cluster_utils import compute_pca, graph_from_coordinates, write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    pipe = config["scopes"][args.scope]["pipeline"]
    root = Path(args.output_root) / "scopes" / args.scope
    prepared = root / "prepared"
    target = root / "embedding"
    target.mkdir(parents=True, exist_ok=True)
    matrix = sparse.load_npz(prepared / "log1p_gex.npz").tocsr()
    features = pd.read_csv(prepared / "features.tsv.gz", sep="\t", dtype=str)
    with threadpool_limits(limits=args.threads):
        coordinates, selected, variance, hvg_audit, model_state = compute_pca(matrix, features, pipe)
        used = coordinates[:, :int(pipe["pca_n_components_used"])]
        graph = graph_from_coordinates(used, int(pipe["neighbor_n_neighbors"]), int(pipe["pca_seed"]), args.threads, str(pipe["neighbor_backend"]))
    np.save(target / "pca.npy", coordinates, allow_pickle=False)
    np.savez_compressed(target / "pca_model.npz", **model_state)
    np.save(target / "explained_variance_ratio.npy", variance, allow_pickle=False)
    sparse.save_npz(target / "knn_graph.npz", graph, compressed=True)
    pd.DataFrame({"feature_index": selected, "gene_id": features.iloc[selected]["gene_id"].to_numpy(),
                  "gene_symbol": features.iloc[selected]["gene_symbol"].to_numpy()}).to_csv(target / "selected_features.tsv", sep="\t", index=False)
    write_json(target / "metadata.json", {
        "study_id": config["study_id"], "scope_id": args.scope, "n_cells": matrix.shape[0],
        "n_features": matrix.shape[1], "n_hvg": len(selected), "pca_fitted": coordinates.shape[1],
        "pca_used": int(pipe["pca_n_components_used"]), "neighbors": int(pipe["neighbor_n_neighbors"]),
        "backend": pipe["neighbor_backend"], "hvg_audit": hvg_audit,
    })
    write_marker(target / ".ok", status="PASS", scope_id=args.scope, n_cells=matrix.shape[0])


if __name__ == "__main__":
    main()
