#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from sklearn.metrics import normalized_mutual_info_score

from cluster_utils import write_json, write_marker


def contingency_outputs(frame: pd.DataFrame, field: str, target: Path, label: str | None = None) -> dict:
    label = label or field
    counts = pd.crosstab(frame["cluster_id"], frame[field], dropna=False)
    counts.to_csv(target / f"cluster_by_{label}_counts.tsv", sep="\t")
    fractions = counts.div(counts.sum(axis=1), axis=0)
    fractions.to_csv(target / f"cluster_by_{label}_fractions.tsv", sep="\t")
    nmi = 0.0 if frame[field].nunique() == 1 else float(
        normalized_mutual_info_score(frame["cluster_id"], frame[field])
    )
    return {
        f"n_{label}_levels": int(frame[field].nunique()),
        f"cluster_{label}_normalized_mutual_information": nmi,
        f"maximum_cluster_dominant_{label}_fraction": float(fractions.max(axis=1).max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.clusters, sep="\t", dtype=str)
    counts = pd.crosstab(frame["cluster_id"], frame["library_id"], dropna=False)
    fractions = counts.div(counts.sum(axis=1), axis=0)
    library_summary = contingency_outputs(frame, "library_id", target, "library")
    per_cluster = pd.DataFrame({"cluster_id": counts.index, "n_cells": counts.sum(axis=1).to_numpy(),
                                "dominant_library_fraction": fractions.max(axis=1).to_numpy(),
                                "n_libraries_observed": (counts > 0).sum(axis=1).to_numpy()})
    per_cluster.to_csv(target / "cluster_technical_metrics.tsv", sep="\t", index=False)
    summary = {
        "n_cells": len(frame),
        "n_clusters": int(frame["cluster_id"].nunique()),
        "n_libraries": int(frame["library_id"].nunique()),
        "cluster_library_normalized_mutual_information": library_summary["cluster_library_normalized_mutual_information"],
        "maximum_cluster_dominant_library_fraction": library_summary["maximum_cluster_dominant_library_fraction"],
        "used_for_resolution_selection": False,
    }
    for field, label in [("experimental_time", "time"), ("condition", "condition")]:
        if field in frame:
            summary.update(contingency_outputs(frame, field, target, label))
            libraries_per_level = frame.groupby(field, dropna=False)["library_id"].nunique()
            summary[f"minimum_libraries_per_{label}"] = int(libraries_per_level.min())
            summary[f"fraction_{label}_levels_with_one_library"] = float((libraries_per_level == 1).mean())
    if "experimental_time" in frame:
        summary["time_library_confounding_warning"] = bool(
            summary["fraction_time_levels_with_one_library"] > 0
        )
    write_json(target / "summary.json", summary)
    write_marker(target / ".ok", status="PASS", nmi=summary["cluster_library_normalized_mutual_information"])


if __name__ == "__main__":
    main()
