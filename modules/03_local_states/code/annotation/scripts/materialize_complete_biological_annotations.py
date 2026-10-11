#!/usr/bin/env python3
"""Build a complete cluster annotation table without conflating labels and anchors."""

import argparse
import csv
from pathlib import Path


STUDY_SLUG = {
    "GSE103221": "gse103221", "GSE118258": "gse118258",
    "GSE150311": "gse150311", "GSE171127": "gse171127",
    "GSE178325": "gse178325", "WOT_GSE115943": "wot",
}
ANCHOR_LABEL = {
    "source_fibroblast": "fibroblast_source",
    "epithelializing": "epithelializing_transition",
    "mouse_ips_naive": "mouse_naive_pluripotent_target",
    "human_naive": "human_naive_pluripotent_target",
    "human_primed": "human_primed_pluripotent_target",
    "trophoblast_like": "trophoblast_like_alternative_fate",
}


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--adjudication", type=Path, required=True)
    parser.add_argument("--non-anchor-labels", type=Path, required=True)
    parser.add_argument("--marker-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    master = read_tsv(args.master)
    adjudication = {(r["study_id"], r["local_cluster_id"]): r for r in read_tsv(args.adjudication)}
    manual_rows = read_tsv(args.non_anchor_labels)
    manual = {(r["study_id"], r["local_cluster_id"]): r for r in manual_rows}
    if len(manual) != len(manual_rows):
        raise ValueError("Duplicate non-anchor annotation key")
    expected_non_anchor = {
        key for key, row in adjudication.items()
        if row["alignment_training_anchor_eligible"] != "true"
    }
    if set(manual) != expected_non_anchor:
        missing = sorted(expected_non_anchor - set(manual))
        extra = sorted(set(manual) - expected_non_anchor)
        raise ValueError(f"Non-anchor annotation closure failed; missing={missing}, extra={extra}")

    markers = {}
    for study, slug in STUDY_SLUG.items():
        path = args.marker_root / slug / "scopes/studywide/markers/top_markers.tsv"
        by_cluster = {}
        for row in read_tsv(path):
            if int(row["rank"]) <= 10:
                by_cluster.setdefault(row["cluster_id"], []).append(row["gene_symbol"])
        markers.update({(study, cluster): ";".join(genes) for cluster, genes in by_cluster.items()})

    output_rows = []
    for row in master:
        key = (row["study_id"], row["local_cluster_id"])
        decision = adjudication[key]
        is_anchor = decision["alignment_training_anchor_eligible"] == "true"
        if is_anchor:
            subtype = decision["proposed_subtype"]
            if subtype not in ANCHOR_LABEL:
                raise ValueError(f"Unmapped anchor subtype: {subtype}")
            label, confidence = ANCHOR_LABEL[subtype], "high"
            rationale = decision["review_evidence_basis"]
            scope = "cross_study_training_anchor"
        else:
            label = manual[key]["descriptive_biological_label"]
            confidence = manual[key]["descriptive_confidence"]
            rationale = manual[key]["biological_rationale"]
            scope = "study_descriptive_only"
        output_rows.append({
            "study_id": row["study_id"], "species": row["species"],
            "local_cluster_id": row["local_cluster_id"], "n_primary_cells": row["n_primary_cells"],
            "descriptive_biological_label": label, "descriptive_confidence": confidence,
            "annotation_scope": scope, "alignment_training_anchor_eligible": str(is_anchor).lower(),
            "local_working_state": row["local_working_state"],
            "observed_day_or_stage": row["observed_day_or_stage"], "observed_route": row["observed_route"],
            "top_markers": markers[key], "biological_rationale": rationale,
            "manual_hold": row["manual_hold"], "library_context_triggered": row["library_context_triggered"],
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_rows[0].keys(), delimiter="\t")
        writer.writeheader(); writer.writerows(output_rows)


if __name__ == "__main__":
    main()
