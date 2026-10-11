#!/usr/bin/env python3
"""Project sensitivity-tier cells into a primary-only PCA/cluster reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from pynndescent import NNDescent
from scipy import sparse
from threadpoolctl import threadpool_limits

from cluster_utils import stratified_indices, write_json, write_marker


def cluster_color_map(cluster_levels: list[str]) -> dict[str, tuple]:
    palette = np.vstack([
        np.asarray(plt.get_cmap(name).colors)
        for name in ("tab20", "tab20b", "tab20c")
    ])
    if len(cluster_levels) > len(palette):
        palette = plt.get_cmap("turbo")(np.linspace(0.02, 0.98, len(cluster_levels)))
    return {cluster_id: tuple(palette[index]) for index, cluster_id in enumerate(cluster_levels)}


def transform(matrix: sparse.csr_matrix, model_path: str) -> np.ndarray:
    with np.load(model_path, allow_pickle=False) as model:
        selected = model["selected_feature_indices"].astype(int)
        mean = model["feature_mean"].astype(np.float32)
        scale = model["feature_scale"].astype(np.float32)
        components = model["components"].astype(np.float32)
        clip = float(model["scale_max_value"])
    dense = matrix[:, selected].toarray().astype(np.float32, copy=False)
    dense -= mean
    dense /= scale
    np.clip(dense, -clip, clip, out=dense)
    return np.asarray(dense @ components.T, dtype=np.float32)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sensitivity-matrix", required=True)
    parser.add_argument("--sensitivity-cells", required=True)
    parser.add_argument("--pca-model", required=True)
    parser.add_argument("--primary-pca", required=True)
    parser.add_argument("--primary-clusters", required=True)
    parser.add_argument("--primary-umap", required=True)
    parser.add_argument("--neighbors", type=int, required=True)
    parser.add_argument("--support-quantile", type=float, required=True)
    parser.add_argument("--pca-components-used", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--max-plot-cells", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    sensitivity_matrix = sparse.load_npz(args.sensitivity_matrix).tocsr()
    sensitivity_cells = pd.read_csv(args.sensitivity_cells, sep="\t", dtype=str)
    primary_pca = np.load(args.primary_pca, mmap_mode="r")[:, :args.pca_components_used]
    primary_clusters = pd.read_csv(args.primary_clusters, sep="\t", dtype=str)
    primary_umap = pd.read_csv(args.primary_umap, sep="\t", dtype={"cell_id": str})
    if len(primary_pca) != len(primary_clusters) or not primary_clusters["cell_id"].equals(primary_umap["cell_id"]):
        raise ValueError("primary projection reference has inconsistent cell order")
    if sensitivity_matrix.shape[0] != len(sensitivity_cells):
        raise ValueError("sensitivity matrix/cell count mismatch")
    cluster_levels = sorted(primary_clusters["cluster_id"].unique())
    if not len(sensitivity_cells):
        empty = sensitivity_cells.copy()
        empty["assigned_cluster"] = pd.Series(dtype=str)
        empty["projection_status"] = pd.Series(dtype=str)
        empty.to_csv(target / "sensitivity_projection.tsv.gz", sep="\t", index=False,
                     compression={"method": "gzip", "mtime": 0})
        write_json(target / "summary.json", {"n_sensitivity_cells": 0, "status": "PASS"})
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        ax.text(0.5, 0.5, "No sensitivity-tier cells", ha="center", va="center")
        ax.axis("off")
        fig.savefig(target / "sensitivity_projection.png", dpi=220)
        plt.close(fig)
        write_marker(target / ".ok", status="PASS", n_sensitivity_cells=0)
        return

    sensitivity_pca = transform(sensitivity_matrix, args.pca_model)[:, :args.pca_components_used]
    np.save(target / "sensitivity_pca.npy", sensitivity_pca, allow_pickle=False)
    k = min(args.neighbors, len(primary_clusters) - 1)
    if k < 2:
        raise ValueError("primary reference is too small for sensitivity projection")
    reference = np.array(primary_pca, dtype=np.float32, copy=True, order="C")
    with threadpool_limits(limits=args.threads):
        index = NNDescent(
            reference,
            n_neighbors=min(max(k + 1, 51), len(primary_clusters) - 1),
            metric="euclidean",
            random_state=args.seed,
            n_jobs=args.threads,
            low_memory=True,
        )
        own_indices, own_distances = index.neighbor_graph
        rows = np.arange(len(primary_clusters))[:, None]
        own_work = own_distances.copy()
        own_work[own_indices == rows] = np.inf
        own_sorted = np.sort(own_work, axis=1)
        primary_k_distance = own_sorted[:, min(k - 1, own_sorted.shape[1] - 1)]
        neighbor_indices, neighbor_distances = index.query(sensitivity_pca, k=k)
    finite_primary = primary_k_distance[np.isfinite(primary_k_distance)]
    if not len(finite_primary):
        raise ValueError("primary neighbor support distribution is empty")
    support_threshold = float(np.quantile(finite_primary, args.support_quantile))
    distance_scale = max(float(np.median(finite_primary)), 1e-6)
    weights = np.exp(-neighbor_distances / distance_scale).astype(np.float32)
    weights /= np.maximum(weights.sum(axis=1, keepdims=True), 1e-12)
    primary_label = primary_clusters["cluster_id"].to_numpy()
    probabilities = np.zeros((len(sensitivity_cells), len(cluster_levels)), dtype=np.float32)
    for index_value, cluster_id in enumerate(cluster_levels):
        probabilities[:, index_value] = np.sum(
            weights * (primary_label[neighbor_indices] == cluster_id), axis=1
        )
    maximum = probabilities.max(axis=1)
    assigned = np.asarray(cluster_levels)[probabilities.argmax(axis=1)]
    entropy = -np.sum(probabilities * np.log(np.maximum(probabilities, 1e-12)), axis=1)
    entropy /= np.log(max(len(cluster_levels), 2))
    kth_distance = neighbor_distances[:, -1]
    in_support = kth_distance <= support_threshold
    projection_status = np.where(~in_support, "out_of_primary_support", np.where(maximum >= 0.6, "high_confidence", "ambiguous"))
    primary_xy = primary_umap[["umap_1", "umap_2"]].to_numpy(float)
    projected_xy = np.sum(weights[:, :, None] * primary_xy[neighbor_indices], axis=1)

    output_columns = [
        column for column in ["study_id", "scope_id", "cell_id", "library_id", "species", "sample_accession", "sample_title", "experimental_time", "condition"]
        if column in sensitivity_cells
    ]
    output = sensitivity_cells[output_columns].copy()
    output["assigned_cluster"] = assigned
    output["maximum_cluster_probability"] = maximum
    output["normalized_cluster_entropy"] = entropy
    output["primary_k_neighbor_distance"] = kth_distance
    output["primary_support_threshold"] = support_threshold
    output["projection_status"] = projection_status
    output["projected_umap_1"] = projected_xy[:, 0]
    output["projected_umap_2"] = projected_xy[:, 1]
    for index_value, cluster_id in enumerate(cluster_levels):
        output[f"probability_{cluster_id}"] = probabilities[:, index_value]
    output.to_csv(target / "sensitivity_projection.tsv.gz", sep="\t", index=False,
                  compression={"method": "gzip", "mtime": 0})
    composition = pd.crosstab(output["assigned_cluster"], output["projection_status"], dropna=False)
    composition.to_csv(target / "projection_cluster_composition.tsv", sep="\t")
    summary = {
        "n_sensitivity_cells": len(output),
        "n_primary_cells": len(primary_clusters),
        "neighbors": k,
        "support_quantile": args.support_quantile,
        "support_threshold": support_threshold,
        "high_confidence_fraction": float(np.mean(projection_status == "high_confidence")),
        "ambiguous_fraction": float(np.mean(projection_status == "ambiguous")),
        "out_of_primary_support_fraction": float(np.mean(projection_status == "out_of_primary_support")),
        "used_to_fit_hvg_pca_neighbors_or_clusters": False,
        "all_sensitivity_cells_retained": True,
        "status": "PASS",
    }
    write_json(target / "summary.json", summary)

    primary_selected = stratified_indices(
        primary_clusters["cluster_id"], min(args.max_plot_cells, len(primary_clusters)), args.seed
    )
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.8), constrained_layout=True)
    for ax in axes:
        ax.scatter(primary_xy[primary_selected, 0], primary_xy[primary_selected, 1], s=0.5,
                   linewidths=0, color="#c7c7c7", alpha=0.35, rasterized=True)
        ax.set(xticks=[], yticks=[], xlabel="Projected UMAP 1", ylabel="Projected UMAP 2")
    cluster_colors = cluster_color_map(cluster_levels)
    axes[0].scatter(projected_xy[:, 0], projected_xy[:, 1],
                    c=[cluster_colors[value] for value in assigned], s=1.2,
                    linewidths=0, rasterized=True)
    axes[0].set_title("Sensitivity cells: assigned primary cluster")
    axes[0].legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markersize=5,
                   markerfacecolor=cluster_colors[cluster_id], markeredgewidth=0,
                   label=cluster_id)
            for cluster_id in cluster_levels
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=min(6, len(cluster_levels)),
        fontsize=6,
        columnspacing=0.8,
        handletextpad=0.3,
        frameon=False,
    )
    confidence = axes[1].scatter(projected_xy[:, 0], projected_xy[:, 1], c=maximum, cmap="viridis",
                                 vmin=0, vmax=1, s=1.2, linewidths=0, rasterized=True)
    axes[1].set_title("Maximum assignment probability")
    fig.colorbar(confidence, ax=axes[1], shrink=0.75)
    support_levels = ["high_confidence", "ambiguous", "out_of_primary_support"]
    support_colors = {
        "high_confidence": "#2ca25f",
        "ambiguous": "#f0a202",
        "out_of_primary_support": "#d73027",
    }
    axes[2].scatter(projected_xy[:, 0], projected_xy[:, 1],
                    c=[support_colors[value] for value in projection_status], s=1.2,
                    linewidths=0, rasterized=True)
    axes[2].set_title("Projection support class")
    axes[2].legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markersize=5,
                   markerfacecolor=support_colors[level], markeredgewidth=0,
                   label=level)
            for level in support_levels
        ],
        loc="best",
        frameon=False,
    )
    fig.savefig(target / "sensitivity_projection.png", dpi=220)
    plt.close(fig)
    write_marker(target / ".ok", status="PASS", n_sensitivity_cells=len(output))


if __name__ == "__main__":
    main()
