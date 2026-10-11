#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cluster_utils import write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--minimum-cluster-cells", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(args.metrics, sep="\t")
    valid = metrics.loc[metrics["valid_for_selection"].astype(str).str.lower().eq("true")].copy()
    if valid.empty:
        raise ValueError("no resolution has a valid multi-cluster partition")
    threshold = float(valid["ci95_low"].max())
    eligible = valid.loc[valid["median_cluster_silhouette"] >= threshold].copy()
    if eligible.empty:
        raise ValueError("chooseR selection produced no eligible resolution")
    chosen = eligible.sort_values(["n_clusters", "resolution"], kind="mergesort").iloc[-1]
    boundary = bool(chosen["resolution"] in {metrics["resolution"].min(), metrics["resolution"].max()})
    small_cluster = bool(chosen["minimum_cluster_cells"] < args.minimum_cluster_cells)
    status = "REVIEW" if boundary or small_cluster else "PASS"
    reasons = []
    if boundary:
        reasons.append("selected_resolution_on_grid_boundary")
    if small_cluster:
        reasons.append("selected_partition_contains_small_cluster")
    decision = {
        "selected_resolution": float(chosen["resolution"]),
        "selection_threshold": threshold,
        "selected_n_clusters": int(chosen["n_clusters"]),
        "selected_median_cluster_silhouette": float(chosen["median_cluster_silhouette"]),
        "selected_ci95_low": float(chosen["ci95_low"]),
        "selected_ci95_high": float(chosen["ci95_high"]),
        "selected_minimum_cluster_cells": int(chosen["minimum_cluster_cells"]),
        "selection_rule": "max n_clusters among resolutions with median_cluster_silhouette >= max(ci95_low)",
        "status": status,
        "review_reasons": reasons,
    }
    write_json(target / "selected_resolution.json", decision)
    pd.DataFrame([decision | {"review_reasons": ";".join(reasons) or "NONE"}]).to_csv(
        target / "selection.tsv", sep="\t", index=False,
    )
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    ax.errorbar(metrics["resolution"], metrics["median_cluster_silhouette"],
                yerr=[metrics["median_cluster_silhouette"] - metrics["ci95_low"],
                      metrics["ci95_high"] - metrics["median_cluster_silhouette"]],
                fmt="o-", color="#2f5d7c", ecolor="#9aa7b0", capsize=2, linewidth=1)
    ax.axhline(threshold, color="#bf3f35", linestyle="--", linewidth=1, label="chooseR threshold")
    ax.axvline(float(chosen["resolution"]), color="#252525", linestyle=":", linewidth=1.2, label="selected")
    ax.set(xlabel="Leiden resolution", ylabel="Median cluster silhouette", ylim=(-1, 1))
    ax.legend(frameon=False)
    fig.savefig(target / "resolution_selection.png", dpi=220)
    plt.close(fig)
    write_marker(target / ".ok", status=status, selected_resolution=decision["selected_resolution"])


if __name__ == "__main__":
    main()
