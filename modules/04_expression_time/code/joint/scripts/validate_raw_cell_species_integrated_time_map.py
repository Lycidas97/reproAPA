#!/usr/bin/env python3
"""Independent validation for the raw-cell human/mouse integrated time map."""

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output
    summary = json.loads((output / "summary.json").read_text())
    preflight = json.loads((output / "preflight.json").read_text())
    nodes = pd.read_csv(output / "raw_cell_nodes.tsv.gz", sep="\t", dtype={"cell_id": str})
    arrows = pd.read_csv(output / "raw_cell_trajectory_arrows.tsv", sep="\t")

    require(summary["status"] == "PASS" and preflight["status"] == "PASS", "Run status is not PASS")
    require(summary["analysis_unit"] == "raw primary cell", "Analysis unit is not raw primary cell")
    require(summary["metacells_used"] is False, "Metacells were used")
    require(len(nodes) == 551201 and nodes.cell_id.is_unique, "Raw-cell registry does not contain 551,201 unique cells")
    require(nodes.species_key.value_counts().to_dict() == {"mouse": 336087, "human": 215114},
            "Raw-cell species counts do not close")
    require({"library_key", "study_id", "cell_id", "UMAP1", "UMAP2"}.issubset(nodes.columns),
            "Raw-cell provenance or coordinates are missing")
    xy = nodes[["UMAP1", "UMAP2"]].to_numpy()
    spatial_neighbors = cKDTree(xy).query(xy, k=31, workers=-1)[1][:, 1:]
    species = nodes.species_key.to_numpy()
    spatial_cross_fraction = (species[spatial_neighbors] != species[:, None]).mean(axis=1)
    require(len(arrows) == 61 and arrows.route.nunique() == 12, "Direct raw-cell arrow count does not close")
    require(arrows.groupby("species").size().to_dict() == {"human": 35, "mouse": 26},
            "Human/mouse direct arrow counts do not close")
    require((arrows.source_cells > 0).all() and (arrows.target_cells > 0).all(),
            "A trajectory endpoint has no raw cells")
    require(arrows.evidence.str.contains("raw-cell").all(), "An arrow lacks raw-cell OT evidence")
    late = arrows.loc[arrows.route.eq("gse171127_late_pooled")]
    require(len(late) == 2 and (late.source_dataset != late.target_dataset).all(),
            "Frozen D6-to-D12 direct mouse trajectory was not retained")

    cross = sparse.load_npz(output / "cross_species_expression_mnn.npz").tocsr()
    graph = sparse.load_npz(output / "integrated_time_graph.npz").tocsr()
    require(cross.shape == (len(nodes), len(nodes)) and graph.shape == cross.shape,
            "Integrated graph node count does not match raw-cell registry")
    require((cross - cross.T).nnz == 0 and (graph - graph.T).nnz == 0,
            "Cross-species or integrated graph is not symmetric")
    human = 215114
    require(cross[:human, :human].nnz == 0 and cross[human:, human:].nnz == 0,
            "Expression MNN graph contains within-species edges")
    require(cross[:human, human:].nnz == summary["cross_species_mnn"]["retained_undirected_edges"],
            "Cross-species MNN edge count does not close")
    require(connected_components(graph, directed=False)[0] == 1, "Integrated graph is disconnected")
    recall = summary["cross_species_mnn"]
    require(min(recall["sampled_exact_neighbor_recall_human_to_mouse"],
                recall["sampled_exact_neighbor_recall_mouse_to_human"]) >= 0.90,
            "FAISS sampled exact-neighbor recall is below the configured floor")
    figure = output / "figures" / "human_mouse_raw_cell_integrated_trajectories.png"
    require(figure.is_file() and figure.stat().st_size > 1000000, "Main raw-cell figure is missing or implausibly small")

    result = {
        "status": "PASS",
        "raw_primary_cells": int(len(nodes)),
        "human_cells": int((nodes.species_key == "human").sum()),
        "mouse_cells": int((nodes.species_key == "mouse").sum()),
        "metacells_used": False,
        "mean_cross_species_30NN_fraction": float(spatial_cross_fraction.mean()),
        "cells_with_any_cross_species_30NN_fraction": float((spatial_cross_fraction > 0).mean()),
        "direct_raw_cell_ot_arrows": int(len(arrows)),
        "direct_raw_cell_ot_routes": int(arrows.route.nunique()),
        "faiss_sampled_exact_recall_human_to_mouse": recall["sampled_exact_neighbor_recall_human_to_mouse"],
        "faiss_sampled_exact_recall_mouse_to_human": recall["sampled_exact_neighbor_recall_mouse_to_human"],
        "integrated_graph_components": 1,
        "cross_species_ot_or_common_clock_created": False,
    }
    (output / "independent_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
