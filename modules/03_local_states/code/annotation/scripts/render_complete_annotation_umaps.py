#!/usr/bin/env python3
"""Render frozen UMAPs with a complete cluster-to-biological-label key."""

import argparse
import csv
import hashlib
import textwrap
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw

from render_annotated_cluster_umaps import (
    BROAD_COLORS, BROAD_STATE, coordinate_rows, draw_panel, font,
)


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--complete-annotations", type=Path, required=True)
    parser.add_argument("--r9-root", type=Path, required=True)
    parser.add_argument("--wot-archive", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    master = {(r["study_id"], r["local_cluster_id"]): r for r in read_tsv(args.master)}
    annotations = {(r["study_id"], r["local_cluster_id"]): r for r in read_tsv(args.complete_annotations)}
    if set(master) != set(annotations):
        raise ValueError("Complete annotation keys differ from annotation master")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index = []
    for study in sorted({key[0] for key in master}):
        rows, authority, source_hash = coordinate_rows(study, args.r9_root, args.wot_archive)
        points = defaultdict(list)
        for row in rows:
            if row["record_type"] == "cell":
                points[row["local_cluster_id"]].append((float(row["x"]), float(row["y"])))
        expected = {key[1] for key in master if key[0] == study}
        if set(points) != expected:
            raise ValueError(f"{study}: UMAP cluster closure failed")
        broad = {cluster: BROAD_STATE[master[(study, cluster)]["local_working_state"]] for cluster in expected}
        image = Image.new("RGB", (3000, 1540), "white")
        draw = ImageDraw.Draw(image)
        draw.text((55, 12), f"{study} | complete biological annotation | original frozen UMAP", fill="#111111", font=font(30, True))
        draw_panel(image, (55, 65, 1840, 1510), "Color: unified broad state; labels: frozen cluster IDs", points, broad, BROAD_COLORS)
        draw.text((1900, 75), "Cluster biological labels", fill="#111111", font=font(25, True))
        draw.text((1900, 112), "H/M/L = descriptive confidence; * = training anchor", fill="#444444", font=font(16))
        y = 150
        for cluster in sorted(expected):
            row = annotations[(study, cluster)]
            confidence = row["descriptive_confidence"][0].upper()
            anchor = "*" if row["alignment_training_anchor_eligible"] == "true" else " "
            draw.rectangle((1900, y + 4, 1915, y + 19), fill=BROAD_COLORS[broad[cluster]])
            draw.text((1925, y), f"{cluster} {confidence}{anchor}", fill="#111111", font=font(17, True))
            label = row["descriptive_biological_label"].replace("_", " ")
            lines = textwrap.wrap(label, width=48) or [label]
            draw.multiline_text((2020, y), "\n".join(lines), fill="#222222", font=font(16), spacing=1)
            y += max(42, 19 * len(lines) + 7)
        output = args.output_dir / f"{study.lower()}_complete_biological_annotation_umap.png"
        image.save(output)
        index.append({"study_id": study, "figure": str(output), "n_cells": sum(map(len, points.values())),
                      "n_clusters": len(expected), "coordinate_authority": authority,
                      "coordinate_sha256": source_hash,
                      "figure_sha256": hashlib.sha256(output.read_bytes()).hexdigest(), "validation": "PASS"})
    with open(args.output_dir / "figure_index.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index[0].keys(), delimiter="\t")
        writer.writeheader(); writer.writerows(index)


if __name__ == "__main__":
    main()
