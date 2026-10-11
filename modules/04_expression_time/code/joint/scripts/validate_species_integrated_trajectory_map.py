#!/usr/bin/env python3
"""Independent structural checks for the species-balanced display coembedding."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    root = args.results
    summary = json.loads((root / "summary.json").read_text())
    nodes = pd.read_csv(root / "common_metacell_nodes.tsv.gz", sep="\t", low_memory=False)
    arrows = pd.read_csv(root / "within_species_study_trajectory_arrows.tsv", sep="\t")
    local = sparse.load_npz(root / "local_expression_kernel.npz")
    temporal = sparse.load_npz(root / "within_species_time_kernel.npz")
    cross = sparse.load_npz(root / "cross_study_expression_mnn_kernel.npz")
    graph = sparse.load_npz(root / "time_constrained_graph.npz")
    checks = pd.read_csv(root / "cross_study_mnn_checks.tsv", sep="\t")
    xy = nodes[["UMAP1", "UMAP2"]].to_numpy(np.float64)
    species = nodes.species_key.to_numpy(str)
    if len(nodes) != 11104 or not nodes.metacell_id.is_unique:
        raise ValueError("Frozen metacell node closure failed")
    if not np.isfinite(xy).all() or graph.shape != (len(nodes), len(nodes)):
        raise ValueError("Coordinates or graph shape is invalid")
    for name, matrix in (("local", local), ("temporal", temporal), ("cross", cross), ("graph", graph)):
        if not np.isfinite(matrix.data).all():
            raise ValueError(f"{name} graph is not finite")
    if (graph - graph.T).nnz:
        raise ValueError("Final UMAP graph is not symmetric")
    time_edges = temporal.tocoo()
    if np.any(species[time_edges.row] != species[time_edges.col]):
        raise ValueError("A temporal edge crosses species")
    arrow_species = arrows.source_group.str.extract(r"^(.*) \|", expand=False)
    if not arrows.route.eq(arrow_species).all() or not arrows.route.isin(nodes.route.astype(str)).all():
        raise ValueError("Displayed trajectory arrow route assignment failed")
    if arrows[["source_x", "source_y", "target_x", "target_y"]].isna().any().any():
        raise ValueError("Displayed arrows are missing coordinates")
    if checks[["recall_a_to_b", "recall_b_to_a"]].to_numpy().min() < .99:
        raise ValueError("Exact MNN neighbor check failed")
    queried = cKDTree(xy).query(xy, k=31, workers=-1)[1][:, 1:]
    mixing = float(np.mean(species[queried] != species[:, None]))
    if not np.isclose(mixing, summary["cross_species_30nn_fraction"], rtol=1e-6, atol=1e-8):
        raise ValueError("Reported cross-species mixing is not reproducible from coordinates")
    overview_png = root / "figures/human_mouse_species_integrated_trajectories.png"
    single_map_png = root / "figures/human_mouse_species_integrated_study_trajectory_map.png"
    width, height = Image.open(overview_png).size
    single_width, single_height = Image.open(single_map_png).size
    result = {
        "status": "PASS",
        "nodes": int(len(nodes)),
        "time_graph_components": int(connected_components(graph, directed=False)[0]),
        "cross_species_temporal_edges": 0,
        "displayed_study_trajectories": int(arrows.route.nunique()),
        "displayed_trajectory_arrows": int(len(arrows)),
        "mnn_pair_checks": int(len(checks)),
        "mnn_min_exact_recall": float(checks[["recall_a_to_b", "recall_b_to_a"]].to_numpy().min()),
        "cross_species_30nn_fraction": mixing,
        "overview_figure_pixels": {"width": int(width), "height": int(height)},
        "single_map_figure_pixels": {"width": int(single_width), "height": int(single_height)},
        "interpretation": "Species integration is a display alignment; temporal OT and arrows remain study/species internal.",
    }
    (root / "independent_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
