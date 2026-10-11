#!/usr/bin/env python3
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap
from sklearn.metrics import silhouette_samples

from cluster_utils import write_json, write_marker


def confidence_interval(values: np.ndarray, n_resamples: int, seed: int):
    values = np.asarray(values, dtype=float)
    if len(values) < 2 or np.allclose(values, values[0]):
        return float(np.median(values)), float(np.median(values))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = bootstrap((values,), np.median, confidence_level=0.95, n_resamples=n_resamples,
                           method="BCa", random_state=np.random.default_rng(seed))
    low, high = result.confidence_interval
    if not np.isfinite([low, high]).all():
        rng = np.random.default_rng(seed)
        medians = np.median(rng.choice(values, size=(n_resamples, len(values)), replace=True), axis=1)
        low, high = np.quantile(medians, [0.025, 0.975])
    return float(low), float(high)


def coclustering_counts(labels: np.ndarray) -> np.ndarray:
    n_cells = labels.shape[1]
    counts = np.zeros((n_cells, n_cells), dtype=np.uint8)
    for row in labels:
        present = row >= 0
        for cluster in np.unique(row[present]):
            members = np.flatnonzero(row == cluster)
            counts[np.ix_(members, members)] += 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--partitions", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--ci-resamples", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    cohort = pd.read_csv(args.cells, sep="\t", dtype=str)
    with np.load(args.partitions, allow_pickle=False) as archive:
        resolutions = archive["resolutions"].astype(float)
        full = archive["full_labels"].astype(np.int16)
        sampled = archive["sampled_labels"].astype(np.int16)
    presence = (sampled[0] >= 0).astype(np.float32)
    observed_together = presence.T @ presence
    if np.any(observed_together <= 0):
        raise ValueError("some cell pairs were never jointly sampled")
    cell_rows = []
    cluster_rows = []
    resolution_rows = []
    cocluster_rows = []
    for r_index, resolution in enumerate(resolutions):
        labels = full[r_index]
        counts = coclustering_counts(sampled[r_index])
        frequency = counts.astype(np.float32) / observed_together
        distance = 1.0 - frequency
        np.fill_diagonal(distance, 0.0)
        unique, sizes = np.unique(labels, return_counts=True)
        valid = 1 < len(unique) < len(labels)
        silhouettes = silhouette_samples(distance, labels, metric="precomputed") if valid else np.full(len(labels), -1.0)
        cluster_scores = []
        for cluster, size in zip(unique, sizes):
            mask = labels == cluster
            score = float(np.mean(silhouettes[mask]))
            cluster_scores.append(score)
            cluster_rows.append({"resolution": resolution, "cluster_id": f"C{int(cluster) + 1:03d}",
                                 "n_cells": int(size), "mean_silhouette": score,
                                 "median_cell_silhouette": float(np.median(silhouettes[mask]))})
        low, high = confidence_interval(np.asarray(cluster_scores), args.ci_resamples, args.seed + r_index)
        resolution_rows.append({"resolution": resolution, "n_clusters": len(unique),
                                "valid_for_selection": valid,
                                "median_cluster_silhouette": float(np.median(cluster_scores)),
                                "ci95_low": low, "ci95_high": high,
                                "minimum_cluster_cells": int(sizes.min()),
                                "minimum_cluster_fraction": float(sizes.min() / len(labels)),
                                "negative_cell_silhouette_fraction": float(np.mean(silhouettes < 0))})
        cell_rows.append(pd.DataFrame({"cell_id": cohort["cell_id"], "resolution": resolution,
                                      "cluster_id": [f"C{int(value) + 1:03d}" for value in labels],
                                      "silhouette": silhouettes}))
        for source in unique:
            source_mask = labels == source
            for destination in unique:
                destination_mask = labels == destination
                values = frequency[np.ix_(source_mask, destination_mask)]
                if source == destination and values.shape[0] > 1:
                    values = values[~np.eye(values.shape[0], dtype=bool)]
                cocluster_rows.append({"resolution": resolution, "cluster_1": f"C{int(source) + 1:03d}",
                                       "cluster_2": f"C{int(destination) + 1:03d}",
                                       "mean_coclustering_frequency": float(values.mean()) if values.size else 1.0})
        del counts, frequency, distance, silhouettes
    pd.concat(cell_rows, ignore_index=True).to_csv(target / "cell_silhouette.tsv.gz", sep="\t", index=False,
                                                    compression={"method": "gzip", "mtime": 0})
    pd.DataFrame(cluster_rows).to_csv(target / "cluster_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(resolution_rows).to_csv(target / "resolution_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(cocluster_rows).to_csv(target / "cluster_coclustering.tsv.gz", sep="\t", index=False,
                                       compression={"method": "gzip", "mtime": 0})
    write_json(target / "metadata.json", {"n_cells": len(cohort), "n_resolutions": len(resolutions),
                                           "ci_method": "BCa bootstrap of median cluster silhouette",
                                           "distance": "1 - pairwise co-clustering frequency"})
    write_marker(target / ".ok", status="PASS", n_cells=len(cohort), n_resolutions=len(resolutions))


if __name__ == "__main__":
    main()
