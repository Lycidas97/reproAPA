#!/usr/bin/env python3
"""Audit transfer from the chooseR evaluation cohort to the frozen full partition."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

from cluster_utils import write_json, write_marker


def resolution_column(frame: pd.DataFrame, resolution: float) -> str:
    expected = f"res_{resolution:g}"
    if expected in frame:
        return expected
    candidates = [column for column in frame if column.startswith("res_")]
    matches = [column for column in candidates if math.isclose(float(column[4:]), resolution, abs_tol=1e-8)]
    if len(matches) != 1:
        raise ValueError(f"cannot identify selected resolution column for {resolution}")
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluation-partitions", required=True)
    parser.add_argument("--selected-json", required=True)
    parser.add_argument("--frozen-clusters", required=True)
    parser.add_argument("--minimum-cluster-cells", type=int, required=True)
    parser.add_argument("--minimum-evaluation-cluster-cells", type=int, required=True)
    parser.add_argument("--minimum-transfer-ari", type=float, required=True)
    parser.add_argument("--minimum-cluster-jaccard", type=float, required=True)
    parser.add_argument("--minimum-cluster-purity", type=float, required=True)
    parser.add_argument("--minimum-cluster-recall", type=float, required=True)
    parser.add_argument("--maximum-cluster-count-drift-fraction", type=float, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)

    decision = json.loads(Path(args.selected_json).read_text(encoding="utf-8"))
    resolution = float(decision["selected_resolution"])
    evaluation = pd.read_csv(args.evaluation_partitions, sep="\t", dtype=str)
    column = resolution_column(evaluation, resolution)
    evaluation = evaluation[["cell_id", column]].rename(columns={column: "evaluation_cluster"})
    frozen = pd.read_csv(args.frozen_clusters, sep="\t", dtype=str)
    if frozen["cell_id"].duplicated().any() or evaluation["cell_id"].duplicated().any():
        raise ValueError("duplicate cell IDs in partition transfer audit")
    joined = evaluation.merge(
        frozen[["cell_id", "cluster_id"]], on="cell_id", how="left", validate="one_to_one"
    )
    if joined["cluster_id"].isna().any():
        raise ValueError("evaluation cohort is not a subset of the frozen full partition")

    table = pd.crosstab(joined["cluster_id"], joined["evaluation_cluster"], dropna=False).sort_index()
    full_sizes = frozen.groupby("cluster_id", sort=True).size()
    row_ind, col_ind = linear_sum_assignment(-table.to_numpy())
    matches = {table.index[row]: table.columns[col] for row, col in zip(row_ind, col_ind)}
    rows = []
    for cluster_id in table.index:
        matched = matches.get(cluster_id)
        eval_n = int(table.loc[cluster_id].sum())
        full_n = int(full_sizes.loc[cluster_id])
        if matched is None:
            overlap = matched_total = 0
            purity = recall = jaccard = 0.0
        else:
            overlap = int(table.loc[cluster_id, matched])
            matched_total = int(table[matched].sum())
            purity = overlap / eval_n if eval_n else 0.0
            recall = overlap / matched_total if matched_total else 0.0
            union = eval_n + matched_total - overlap
            jaccard = overlap / union if union else 0.0
        reasons = []
        if full_n < args.minimum_cluster_cells:
            reasons.append("rare_full_cluster")
        if eval_n < args.minimum_evaluation_cluster_cells:
            reasons.append("underrepresented_in_evaluation_cohort")
        if matched is None:
            reasons.append("no_one_to_one_evaluation_match")
        elif (
            jaccard < args.minimum_cluster_jaccard
            and (purity < args.minimum_cluster_purity or recall < args.minimum_cluster_recall)
        ):
            reasons.append("unstable_evaluation_to_full_transfer")
        rows.append({
            "cluster_id": cluster_id,
            "full_n_cells": full_n,
            "evaluation_n_cells": eval_n,
            "evaluation_fraction_of_full": eval_n / full_n,
            "matched_evaluation_cluster": matched or "UNMATCHED",
            "overlap_n_cells": overlap,
            "matched_evaluation_cluster_n_cells": matched_total,
            "transfer_purity": purity,
            "transfer_recall": recall,
            "transfer_jaccard": jaccard,
            "cluster_release_class": "stable_candidate" if not reasons else "rare_or_unresolved",
            "review_reasons": ";".join(reasons) or "NONE",
        })
    cluster_audit = pd.DataFrame(rows)
    cluster_audit.to_csv(target / "cluster_transfer_audit.tsv", sep="\t", index=False)
    table.to_csv(target / "evaluation_to_full_contingency.tsv", sep="\t")

    ari = float(adjusted_rand_score(joined["evaluation_cluster"], joined["cluster_id"]))
    ami = float(adjusted_mutual_info_score(joined["evaluation_cluster"], joined["cluster_id"]))
    evaluation_n_clusters = int(evaluation["evaluation_cluster"].nunique())
    full_n_clusters = int(frozen["cluster_id"].nunique())
    drift = abs(full_n_clusters - evaluation_n_clusters)
    drift_limit = max(2, int(math.ceil(args.maximum_cluster_count_drift_fraction * evaluation_n_clusters)))
    overall_reasons = []
    if ari < args.minimum_transfer_ari:
        overall_reasons.append("evaluation_to_full_ari_below_threshold")
    if drift > drift_limit:
        overall_reasons.append("cluster_count_drift_above_threshold")
    unresolved = int(cluster_audit["cluster_release_class"].eq("rare_or_unresolved").sum())
    if unresolved:
        overall_reasons.append("one_or_more_clusters_require_review")
    status = "PASS" if not overall_reasons else "REVIEW"
    summary = {
        "selected_resolution": resolution,
        "evaluation_n_cells": len(evaluation),
        "full_n_cells": len(frozen),
        "evaluation_n_clusters": evaluation_n_clusters,
        "full_n_clusters": full_n_clusters,
        "cluster_count_absolute_drift": drift,
        "cluster_count_drift_limit": drift_limit,
        "adjusted_rand_index": ari,
        "adjusted_mutual_information": ami,
        "unresolved_cluster_count": unresolved,
        "all_full_cells_retained": True,
        "status": status,
        "review_reasons": overall_reasons,
    }
    write_json(target / "summary.json", summary)

    normalized = table.div(table.sum(axis=1), axis=0).fillna(0.0)
    width = max(7.0, 0.45 * normalized.shape[1] + 3.0)
    height = max(5.0, 0.36 * normalized.shape[0] + 2.0)
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    image = ax.imshow(normalized.to_numpy(), aspect="auto", interpolation="nearest", vmin=0, vmax=1, cmap="Blues")
    ax.set_xticks(np.arange(normalized.shape[1]), normalized.columns, rotation=90)
    ax.set_yticks(np.arange(normalized.shape[0]), normalized.index)
    ax.set(xlabel="Evaluation-cohort cluster", ylabel="Frozen full-data cluster",
           title=f"Evaluation-to-full transfer (ARI={ari:.3f}, status={status})")
    fig.colorbar(image, ax=ax, label="Fraction within full-data cluster")
    fig.savefig(target / "evaluation_to_full_transfer.png", dpi=220)
    plt.close(fig)
    write_marker(target / ".ok", status=status, ari=ari, unresolved_clusters=unresolved)


if __name__ == "__main__":
    main()
