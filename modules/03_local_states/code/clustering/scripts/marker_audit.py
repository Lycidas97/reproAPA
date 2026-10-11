#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from cluster_utils import write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--top-n", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    matrix = sparse.load_npz(args.matrix).tocsr()
    features = pd.read_csv(args.features, sep="\t", dtype=str)
    partition = pd.read_csv(args.clusters, sep="\t", dtype=str)
    labels, encoded = np.unique(partition["cluster_id"], return_inverse=True)
    indicator = sparse.csr_matrix((np.ones(len(encoded), dtype=np.float32), (encoded, np.arange(len(encoded)))),
                                  shape=(len(labels), len(encoded)))
    sums = np.asarray((indicator @ matrix).toarray(), dtype=np.float64)
    detected = np.asarray((indicator @ (matrix > 0).astype(np.float32)).toarray(), dtype=np.float64)
    sizes = np.bincount(encoded).astype(float)
    total_sum = sums.sum(axis=0)
    total_detected = detected.sum(axis=0)
    rows = []
    for index, label in enumerate(labels):
        outside_n = len(encoded) - sizes[index]
        mean_in = sums[index] / sizes[index]
        mean_out = (total_sum - sums[index]) / max(outside_n, 1)
        pct_in = detected[index] / sizes[index]
        pct_out = (total_detected - detected[index]) / max(outside_n, 1)
        score = (mean_in - mean_out) * np.maximum(pct_in - pct_out, 0)
        order = np.argsort(-score, kind="stable")[:args.top_n]
        for rank, feature_index in enumerate(order, start=1):
            rows.append({"cluster_id": label, "rank": rank, "gene_id": features.iloc[feature_index]["gene_id"],
                         "gene_symbol": features.iloc[feature_index]["gene_symbol"],
                         "mean_log1p_in": mean_in[feature_index], "mean_log1p_out": mean_out[feature_index],
                         "detected_fraction_in": pct_in[feature_index], "detected_fraction_out": pct_out[feature_index],
                         "marker_score": score[feature_index]})
    pd.DataFrame(rows).to_csv(target / "top_markers.tsv", sep="\t", index=False)
    write_marker(target / ".ok", status="PASS", n_clusters=len(labels), top_n=args.top_n)


if __name__ == "__main__":
    main()

