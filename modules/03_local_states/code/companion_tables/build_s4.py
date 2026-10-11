"""Rebuild manuscript Table S4 from frozen partition and sensitivity records."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
from statistics import median


SOURCE_LABEL = "notes/2026-09-13-methods-integration-annotation/sources/remote-release-details.json"
WOT_SOURCE_LABEL = "notes/2026-09-13-methods-integration-annotation/sources/wot-discovery.json"
FIELDS = ("stability_runs", "median_ari", "minimum_ari",
          "median_broad_family_precision", "minimum_broad_family_precision")


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build(release_path, annotations_path, runs_path, overlaps_path):
    release = json.loads(release_path.read_text())
    partitions = []
    for unit, path in release["scopes"].items():
        if unit == "wot":
            continue

        def record(suffix):
            return json.loads(release["files"][path + "/" + suffix]["text"])

        frozen = record("frozen/metadata.json")
        selected = record("selection/selected_resolution.json")
        embedding = record("embedding/metadata.json")
        cohort = record("cohort/metadata.json")
        sweep = record("sweep/metadata.json")
        partitions.append(dict(
            analysis_unit=unit, n_cells=frozen["n_cells"],
            n_local_states=frozen["n_clusters"], n_hvg=embedding["n_hvg"],
            n_pcs=embedding["pca_used"], n_neighbors=embedding["neighbors"],
            selected_resolution=round(selected["selected_resolution"], 3),
            selection_method="chooseR-style subsampling", evaluation_cells=cohort["evaluation_n"],
            evaluation_clusters=selected["selected_n_clusters"], seed=cohort["seed"],
            resolution_grid=";".join(str(round(value, 3)) for value in sweep["resolutions"]),
            selection_status=selected["status"],
            review_reasons=";".join(selected["review_reasons"]), source_partition=SOURCE_LABEL,
        ))
    partitions.append(dict(
        analysis_unit="wot_common_to_2i", n_cells=158631, n_local_states=12,
        n_hvg=2500, n_pcs=30, n_neighbors=30, selected_resolution=0.45,
        selection_method="previously frozen route-specific partition",
        evaluation_cells="not applicable", evaluation_clusters="not applicable",
        seed=218, resolution_grid="0.45", selection_status="frozen partition reused",
        review_reasons="not selected by chooseR; separate source", source_partition=WOT_SOURCE_LABEL,
    ))
    assert len(partitions) == 11
    assert sum(int(row["n_cells"]) for row in partitions) == 551201
    assert sum(int(row["n_local_states"]) for row in partitions) == 111

    families = {row["original_state_id"]: row["display_group"] for row in read_tsv(annotations_path)}
    assert len(families) == 111
    runs = read_tsv(runs_path)
    ari_by_run = {}
    for row in runs:
        if row["role"] == "sensitivity":
            key = row["unit"], row["resolution"], row["seed"]
            assert key not in ari_by_run
            ari_by_run[key] = float(row["ARI"])
    assert len(runs) == 110 and len(ari_by_run) == 99

    counts_by_run = defaultdict(lambda: defaultdict(int))
    with gzip.open(overlaps_path, "rt", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["role"] != "sensitivity":
                continue
            key = row["unit"], row["resolution"], row["seed"]
            assert key in ari_by_run and row["original_state_id"] in families
            counts_by_run[key][row["alternative_cluster"], families[row["original_state_id"]]] += int(row["cells"])
    assert set(counts_by_run) == set(ari_by_run)

    cells_by_unit = {row["analysis_unit"]: int(row["n_cells"]) for row in partitions}
    metrics_by_unit = defaultdict(list)
    for key, counts in counts_by_run.items():
        cluster_totals = Counter()
        for (cluster, _), count in counts.items():
            cluster_totals[cluster] += count
        cells = sum(cluster_totals.values())
        assert cells == cells_by_unit[key[0]]
        precision = sum(count * count / cluster_totals[cluster]
                        for (cluster, _), count in counts.items()) / cells
        metrics_by_unit[key[0]].append((ari_by_run[key], precision))
    assert set(metrics_by_unit) == set(cells_by_unit)

    for row in partitions:
        values = metrics_by_unit[row["analysis_unit"]]
        assert len(values) == 9
        ari = [value[0] for value in values]
        precision = [value[1] for value in values]
        row.update(stability_runs=9, median_ari=round(median(ari), 6),
                   minimum_ari=round(min(ari), 6),
                   median_broad_family_precision=round(median(precision), 6),
                   minimum_broad_family_precision=round(min(precision), 6))
    return partitions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-details", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--partition-runs", required=True, type=Path)
    parser.add_argument("--partition-overlaps", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = build(args.release_details, args.annotations, args.partition_runs, args.partition_overlaps)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} partitions; {sum(row['n_cells'] for row in rows)} cells")


if __name__ == "__main__":
    main()
