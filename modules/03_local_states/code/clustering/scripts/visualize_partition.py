#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap

from cluster_utils import stratified_indices, write_marker


def categorical_scatter(ax, coordinates, values, title, legend_limit=20):
    labels, encoded = np.unique(values, return_inverse=True)
    palette = plt.get_cmap("tab20", max(len(labels), 2))
    ax.scatter(coordinates[:, 0], coordinates[:, 1], c=encoded, cmap=palette, s=1.0, linewidths=0, rasterized=True)
    ax.set(title=title, xlabel="UMAP 1", ylabel="UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    if len(labels) <= legend_limit:
        handles = [plt.Line2D([], [], marker="o", linestyle="", color=palette(index), label=label, markersize=4)
                   for index, label in enumerate(labels)]
        ax.legend(handles=handles, frameon=False, fontsize=7, markerscale=1, loc="best")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pca", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--max-plot-cells", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    coordinates = np.load(args.pca, mmap_mode="r")
    frame = pd.read_csv(args.clusters, sep="\t", dtype=str)
    reducer = umap.UMAP(n_neighbors=30, min_dist=0.3, metric="euclidean", random_state=args.seed,
                        n_components=2, low_memory=True)
    embedding = reducer.fit_transform(np.asarray(coordinates, dtype=np.float32)).astype(np.float32)
    metadata_columns = [
        column for column in ["experimental_time", "condition", "sample_accession"] if column in frame
    ]
    output = frame[["cell_id", "cluster_id", "library_id", *metadata_columns]].copy()
    output["umap_1"] = embedding[:, 0]
    output["umap_2"] = embedding[:, 1]
    output.to_csv(target / "umap_coordinates.tsv.gz", sep="\t", index=False,
                  compression={"method": "gzip", "mtime": 0})
    selected = stratified_indices(frame["cluster_id"], min(args.max_plot_cells, len(frame)), args.seed)
    panels = [("cluster_id", "Selected Leiden clusters"), ("library_id", "Library / sampling stratum")]
    if "experimental_time" in frame:
        panels.append(("experimental_time", "Experimental sampling time"))
    if "condition" in frame:
        panels.append(("condition", "Experimental condition"))
    ncols = 2
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 5.5 * nrows), constrained_layout=True)
    axes = np.asarray(axes).reshape(-1)
    for ax, (column, title) in zip(axes, panels):
        categorical_scatter(ax, embedding[selected], frame.iloc[selected][column].to_numpy(), title)
    for ax in axes[len(panels):]:
        ax.axis("off")
    fig.savefig(target / "partition_umap.png", dpi=220)
    plt.close(fig)
    write_marker(target / ".ok", status="PASS", n_cells=len(frame), plotted_n=len(selected))


if __name__ == "__main__":
    main()
