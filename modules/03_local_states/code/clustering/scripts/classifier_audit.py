#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix
from sklearn.model_selection import train_test_split

from cluster_utils import stratified_indices, write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pca", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--max-cells", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    coordinates = np.load(args.pca, mmap_mode="r")
    frame = pd.read_csv(args.clusters, sep="\t", dtype=str)
    selected = stratified_indices(frame["cluster_id"], min(args.max_cells, len(frame)), args.seed)
    x = np.asarray(coordinates[selected], dtype=np.float32)
    y = frame.iloc[selected]["cluster_id"].to_numpy()
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=args.seed, stratify=y)
    model = LogisticRegression(max_iter=1000, solver="lbfgs", class_weight="balanced")
    model.fit(x_train, y_train)
    predicted = model.predict(x_test)
    labels = np.unique(y)
    matrix = confusion_matrix(y_test, predicted, labels=labels, normalize="true")
    score = float(balanced_accuracy_score(y_test, predicted))
    pd.DataFrame(matrix, index=labels, columns=labels).to_csv(target / "normalized_confusion_matrix.tsv", sep="\t")
    per_cluster = pd.DataFrame({"cluster_id": labels, "held_out_recall": np.diag(matrix)})
    per_cluster.to_csv(target / "per_cluster_recall.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(max(5, len(labels) * 0.35), max(4, len(labels) * 0.35)), constrained_layout=True)
    image = ax.imshow(matrix, vmin=0, vmax=1, cmap="viridis")
    ax.set(xticks=np.arange(len(labels)), yticks=np.arange(len(labels)), xticklabels=labels, yticklabels=labels,
           xlabel="Predicted cluster", ylabel="Observed cluster", title=f"Held-out balanced accuracy = {score:.3f}")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.tick_params(axis="y", labelsize=7)
    fig.colorbar(image, ax=ax, label="Row-normalized fraction")
    fig.savefig(target / "classifier_confusion.png", dpi=220)
    plt.close(fig)
    write_json(target / "summary.json", {"balanced_accuracy": score, "evaluation_n": len(selected),
                                         "classifier": "multinomial logistic regression on frozen PCA",
                                         "used_for_resolution_selection": False})
    write_marker(target / ".ok", status="PASS", balanced_accuracy=score)


if __name__ == "__main__":
    main()

