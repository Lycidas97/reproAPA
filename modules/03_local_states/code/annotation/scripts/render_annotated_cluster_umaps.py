#!/usr/bin/env python3
"""Project frozen cluster annotations onto the original study UMAPs."""

import argparse
import csv
import gzip
import hashlib
import io
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


STUDY_SOURCES = {
    "GSE103221": "ANN05_gse103221_r9/release/figures/data/frozen_clusters.tsv.gz",
    "GSE118258": "ANN06_gse118258_r9/release/figures/data/frozen_clusters.tsv.gz",
    "GSE150311": "ANN07_gse150311_r9/release/figures/data/frozen_clusters.tsv.gz",
    "GSE171127": "ANN08_gse171127_r9/release/figures/data/frozen_clusters.tsv.gz",
    "GSE178325": "ANN09_gse178325_r9/release/figures/data/frozen_clusters.tsv.gz",
}

BROAD_STATE = {
    "source_mesenchymal_myeloid_program": "local_special_source",
    "source_like": "source_retaining",
    "epithelializing_intermediate_candidate": "reprogramming_transition",
    "source_target_mixed_transition": "reprogramming_transition",
    "mixed_identity_candidate": "ambiguous_mixed",
    "mouse_ips_naive_candidate": "pluripotency_target_like",
    "mouse_ips_naive_target_stable": "pluripotency_target_like",
    "human_naive_candidate": "pluripotency_target_like",
    "human_primed_candidate": "pluripotency_target_like",
    "target_program_high": "pluripotency_target_like",
    "trophoblast_trophectoderm_like": "alternative_fate",
    "neural_alternative": "alternative_fate",
    "process_dominant": "process_dominant",
    "unresolved": "unresolved",
}

LOCAL_COLORS = {
    "source_like": "#4E79A7", "source_mesenchymal_myeloid_program": "#8C564B",
    "epithelializing_intermediate_candidate": "#F28E2B", "source_target_mixed_transition": "#FF9D4D",
    "mixed_identity_candidate": "#B07AA1", "mouse_ips_naive_candidate": "#59A14F",
    "mouse_ips_naive_target_stable": "#2E8B57", "human_naive_candidate": "#76B7B2",
    "human_primed_candidate": "#1F9E89", "target_program_high": "#277DA1",
    "trophoblast_trophectoderm_like": "#E15759", "neural_alternative": "#AF7AA1",
    "process_dominant": "#EDC948", "unresolved": "#BAB0AC",
}
BROAD_COLORS = {
    "source_retaining": "#4E79A7", "reprogramming_transition": "#F28E2B",
    "pluripotency_target_like": "#59A14F", "alternative_fate": "#E15759",
    "ambiguous_mixed": "#B07AA1", "process_dominant": "#EDC948",
    "local_special_source": "#8C564B", "unresolved": "#BAB0AC",
}
STATUS_COLORS = {
    "training_anchor": "#2CA02C", "validation_reserve": "#1F77B4",
    "study_specific_identity": "#9467BD", "reference_only": "#17BECF",
    "hold_or_context": "#D62728", "non_anchor": "#BDBDBD",
}


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def anchor_status(row):
    decision = row["review_decision"]
    if row["alignment_training_anchor_eligible"] == "true":
        return "training_anchor"
    if decision == "reserve_validation":
        return "validation_reserve"
    if decision == "accepted_study_specific_identity":
        return "study_specific_identity"
    if row["anchor_tier_initial"] == "R_reference_only":
        return "reference_only"
    if row["manual_hold"] == "true" or row["library_context_triggered"] == "true":
        return "hold_or_context"
    return "non_anchor"


def coordinate_rows(study, r9_root, wot_archive):
    if study != "WOT_GSE115943":
        path = r9_root / STUDY_SOURCES[study]
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t")), str(path), hashlib.sha256(path.read_bytes()).hexdigest()
    member = "R10B_WOT_RELEASE/release/figures/data/frozen_clusters.tsv.gz"
    with tarfile.open(wot_archive, "r:gz") as archive:
        payload = archive.extractfile(member).read()
    with gzip.open(io.BytesIO(payload), "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows, f"{wot_archive}::{member}", hashlib.sha256(payload).hexdigest()


def font(size, bold=False):
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    try:
        return ImageFont.truetype(name, size)
    except OSError:
        return ImageFont.load_default()


def transform(points, box):
    x0, y0, x1, y1 = box
    xs = [p[0] for p in points]; ys = [p[1] for p in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    dx = max_x - min_x or 1; dy = max_y - min_y or 1
    scale = min((x1 - x0) / dx, (y1 - y0) / dy)
    ox = x0 + ((x1 - x0) - dx * scale) / 2
    oy = y0 + ((y1 - y0) - dy * scale) / 2
    return [(int(ox + (x - min_x) * scale), int(y1 - ((y - min_y) * scale + (oy - y0)))) for x, y in points]


def draw_panel(image, box, title, cluster_points, labels, colors):
    draw = ImageDraw.Draw(image)
    title_font, label_font, legend_font = font(25, True), font(15, True), font(14)
    draw.text((box[0], 48), title, fill="#111111", font=title_font)
    all_raw = [point for values in cluster_points.values() for point in values]
    all_px = transform(all_raw, (box[0], box[1] + 55, box[2], box[3] - 115))
    px_by_cluster, offset = {}, 0
    for cluster, values in cluster_points.items():
        px_by_cluster[cluster] = all_px[offset:offset + len(values)]; offset += len(values)
    for cluster in sorted(px_by_cluster, key=lambda key: len(px_by_cluster[key]), reverse=True):
        draw.point(px_by_cluster[cluster], fill=colors[labels[cluster]])
    occupied = []
    candidates = [(0, 0), (0, -26), (0, 26), (-35, 0), (35, 0), (-35, -26), (35, -26), (-35, 26), (35, 26)]
    for cluster in sorted(px_by_cluster, key=lambda key: len(px_by_cluster[key]), reverse=True):
        points = px_by_cluster[cluster]
        cx = int(sum(p[0] for p in points) / len(points)); cy = int(sum(p[1] for p in points) / len(points))
        text = cluster
        placement = None
        for dx, dy in candidates:
            tx, ty = cx + dx, cy + dy
            bbox = draw.textbbox((tx, ty), text, font=label_font, anchor="mm")
            padded = (bbox[0]-3, bbox[1]-2, bbox[2]+3, bbox[3]+2)
            if padded[0] < box[0] or padded[2] > box[2] or padded[1] < box[1]+50 or padded[3] > box[3]-115:
                continue
            if not any(not (padded[2] < old[0] or padded[0] > old[2] or padded[3] < old[1] or padded[1] > old[3]) for old in occupied):
                placement = (tx, ty, padded)
                break
        if placement is None:
            tx, ty = cx, cy
            bbox = draw.textbbox((tx, ty), text, font=label_font, anchor="mm")
            placement = (tx, ty, (bbox[0]-3, bbox[1]-2, bbox[2]+3, bbox[3]+2))
        tx, ty, bbox = placement
        occupied.append(bbox)
        if (tx, ty) != (cx, cy):
            draw.line((cx, cy, tx, ty), fill="#666666", width=1)
        draw.rounded_rectangle((bbox[0]-3, bbox[1]-2, bbox[2]+3, bbox[3]+2), radius=3, fill="#FFFFFFDD", outline="#555555")
        draw.text((tx, ty), text, fill="#111111", font=label_font, anchor="mm")
    legend_y = box[3] - 90
    x = box[0]
    for value in sorted(set(labels.values())):
        item_width = 28 + int(draw.textlength(value, font=legend_font))
        if x + item_width > box[2]:
            x = box[0]; legend_y += 23
        draw.rectangle((x, legend_y, x + 14, legend_y + 14), fill=colors[value])
        draw.text((x + 20, legend_y - 2), value, fill="#222222", font=legend_font)
        x += item_width


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--adjudication", type=Path, required=True)
    parser.add_argument("--r9-root", type=Path, required=True)
    parser.add_argument("--wot-archive", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    master = read_tsv(args.master); adjudication = read_tsv(args.adjudication)
    master_by = {(r["study_id"], r["local_cluster_id"]): r for r in master}
    adj_by = {(r["study_id"], r["local_cluster_id"]): r for r in adjudication}
    if set(master_by) != set(adj_by):
        raise ValueError("Master and adjudication cluster keys differ")
    unknown = sorted({r["local_working_state"] for r in master} - set(BROAD_STATE))
    if unknown:
        raise ValueError(f"Unmapped local states: {unknown}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index_rows = []
    for study in sorted({r["study_id"] for r in master}):
        coord_rows, authority, source_hash = coordinate_rows(study, args.r9_root, args.wot_archive)
        coord_rows = [r for r in coord_rows if r["record_type"] == "cell"]
        cluster_points = defaultdict(list)
        for row in coord_rows:
            cluster_points[row["local_cluster_id"]].append((float(row["x"]), float(row["y"])))
        expected = {k[1]: int(v["n_primary_cells"]) for k, v in master_by.items() if k[0] == study}
        observed = {k: len(v) for k, v in cluster_points.items()}
        if observed != expected:
            raise ValueError(f"{study}: coordinate counts differ from annotation master")
        local = {c: master_by[(study, c)]["local_working_state"] for c in expected}
        broad = {c: BROAD_STATE[local[c]] for c in expected}
        status = {c: anchor_status(adj_by[(study, c)]) for c in expected}
        image = Image.new("RGB", (3300, 1120), "white")
        draw = ImageDraw.Draw(image)
        draw.text((55, 12), f"{study} | original frozen UMAP | cluster membership unchanged", fill="#111111", font=font(30, True))
        boxes = [(55, 65, 1075, 1090), (1140, 65, 2160, 1090), (2225, 65, 3245, 1090)]
        draw_panel(image, boxes[0], "Detailed working annotation", cluster_points, local, LOCAL_COLORS)
        draw_panel(image, boxes[1], "Unified broad state", cluster_points, broad, BROAD_COLORS)
        draw_panel(image, boxes[2], "Alignment evidence role", cluster_points, status, STATUS_COLORS)
        output = args.output_dir / f"{study.lower()}_annotated_cluster_umap.png"
        image.save(output, optimize=True)
        index_rows.append({"study_id": study, "figure": str(output), "n_cells": len(coord_rows),
                           "n_clusters": len(expected), "coordinate_authority": authority,
                           "coordinate_sha256": source_hash, "figure_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
                           "projection_validation": "PASS"})
    with open(args.output_dir / "figure_index.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index_rows[0].keys(), delimiter="\t")
        writer.writeheader(); writer.writerows(index_rows)


if __name__ == "__main__":
    main()
