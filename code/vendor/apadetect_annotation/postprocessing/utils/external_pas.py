"""External PAS resource loading and strand-aware overlap support."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


EXTERNAL_PAS_COLUMNS = [
    "pas_id",
    "source_name",
    "source_species",
    "source_version",
    "source_available_for_species",
    "overlap_status",
    "distance_to_source_pas",
    "distance_bin",
    "source_pas_id",
    "source_region_class",
]


def compute_external_pas_overlap(
    pas_map: pd.DataFrame,
    resources: Iterable[dict | str] | None,
    annotation_species: str,
    default_window_bp: int = 50,
) -> pd.DataFrame:
    """Return one PAS-by-source evidence row for every configured resource.

    Available BED/TSV resources are searched on normalized chromosome and the
    same strand. The nearest source PAS is recorded even when it falls outside
    the positive-support window, making ``absent`` distinct from unavailable.
    """
    resources = list(resources or [])
    if not resources:
        return _unavailable_rows(pas_map, {}, annotation_species)

    frames = []
    for raw_resource in resources:
        resource = _normalize_resource(raw_resource)
        source_species = str(resource.get("species", annotation_species))
        configured_status = str(resource.get("status", "available")).lower()
        species_matches = not source_species or source_species.lower() == str(annotation_species).lower()
        if configured_status == "unavailable" or not species_matches:
            frames.append(_unavailable_rows(pas_map, resource, source_species or annotation_species))
            continue

        table = load_external_pas_resource(resource)
        window_bp = int(resource.get("window_bp", default_window_bp))
        frames.append(
            _match_resource(
                pas_map,
                table,
                resource,
                source_species or annotation_species,
                window_bp=max(0, window_bp),
                sparse=configured_status == "sparse",
            )
        )
    return pd.concat(frames, ignore_index=True)[EXTERNAL_PAS_COLUMNS]


def load_external_pas_resource(resource: dict) -> pd.DataFrame:
    """Load a headered TSV/CSV or headerless BED into a normalized PAS table."""
    path = Path(str(resource.get("path", "")))
    if not path.exists():
        raise FileNotFoundError(f"External PAS resource not found: {path}")

    resource_format = str(resource.get("format", "auto")).lower()
    if resource_format == "auto":
        suffixes = "".join(path.suffixes).lower()
        resource_format = "bed" if ".bed" in suffixes else "tsv"

    if resource_format == "bed":
        raw = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str).fillna("")
        if raw.shape[1] < 3:
            raise ValueError(f"External BED requires at least 3 columns: {path}")
        default_names = ["chrom", "start", "end", "source_pas_id", "score", "strand"]
        raw.columns = default_names[: raw.shape[1]] + [f"extra_{idx}" for idx in range(raw.shape[1] - 6)]
    elif resource_format in {"tsv", "csv"}:
        separator = "," if resource_format == "csv" else "\t"
        raw = pd.read_csv(path, sep=separator, comment="#", dtype=str).fillna("")
    else:
        raise ValueError(f"Unsupported external PAS resource format: {resource_format}")

    chrom_col = _column(resource, "chrom_column", raw, ["chrom", "chr", "seqname"])
    strand_col = _column(resource, "strand_column", raw, ["strand"])
    cut_col = _column(resource, "cut_site_column", raw, ["cut_site", "position", "pas_position"], required=False)
    start_col = _column(resource, "start_column", raw, ["start", "chromStart"], required=cut_col is None)
    end_col = _column(resource, "end_column", raw, ["end", "chromEnd"], required=cut_col is None)
    id_col = _column(resource, "id_column", raw, ["source_pas_id", "pas_id", "name", "id"], required=False)
    region_col = _column(
        resource,
        "region_class_column",
        raw,
        ["source_region_class", "region_class", "region"],
        required=False,
    )

    if cut_col is not None:
        cut_site = pd.to_numeric(raw[cut_col], errors="coerce")
        if bool(resource.get("one_based", False)):
            cut_site = cut_site - 1
    else:
        start = pd.to_numeric(raw[start_col], errors="coerce")
        end = pd.to_numeric(raw[end_col], errors="coerce")
        cut_site = np.floor((start + end - 1) / 2)

    normalized = pd.DataFrame(
        {
            "chrom": raw[chrom_col].astype(str),
            "chrom_key": raw[chrom_col].map(_chrom_key),
            "cut_site": cut_site,
            "strand": raw[strand_col].astype(str),
            "source_pas_id": raw[id_col].astype(str) if id_col is not None else "",
            "source_region_class": raw[region_col].astype(str) if region_col is not None else "",
        }
    )
    normalized = normalized.dropna(subset=["cut_site"])
    normalized = normalized[normalized["strand"].isin(["+", "-"])].copy()
    normalized["cut_site"] = normalized["cut_site"].astype(int)
    if normalized.empty:
        raise ValueError(f"External PAS resource has no valid chrom/cut_site/strand rows: {path}")
    return normalized.sort_values(["chrom_key", "strand", "cut_site", "source_pas_id"]).reset_index(drop=True)


def summarize_external_status(values: pd.Series) -> str:
    """Aggregate source-specific evidence without treating absence as support."""
    statuses = [str(value) for value in values if str(value)]
    if not statuses:
        return "unavailable"
    if "supported" in statuses:
        return "supported"
    if "absent" in statuses:
        return "absent"
    if "unavailable_or_sparse" in statuses:
        return "unavailable_or_sparse"
    return ";".join(sorted(set(statuses)))


def _normalize_resource(resource: dict | str) -> dict:
    if isinstance(resource, dict):
        return dict(resource)
    return {"name": Path(str(resource)).stem, "path": str(resource), "status": "available"}


def _match_resource(
    pas_map: pd.DataFrame,
    source: pd.DataFrame,
    resource: dict,
    source_species: str,
    window_bp: int,
    sparse: bool,
) -> pd.DataFrame:
    lookup = {
        (str(chrom), str(strand)): group.reset_index(drop=True)
        for (chrom, strand), group in source.groupby(["chrom_key", "strand"], sort=False)
    }
    rows = []
    for pas in pas_map.itertuples(index=False):
        group = lookup.get((_chrom_key(pas.chrom), str(pas.strand)))
        nearest = _nearest_source(group, int(pas.cut_site))
        if nearest is None:
            status = "unavailable_or_sparse" if sparse else "absent"
            distance = ""
            source_pas_id = ""
            source_region = ""
        else:
            distance = abs(int(nearest.cut_site) - int(pas.cut_site))
            status = "supported" if distance <= window_bp else ("unavailable_or_sparse" if sparse else "absent")
            source_pas_id = nearest.source_pas_id
            source_region = nearest.source_region_class
        rows.append(
            {
                "pas_id": pas.pas_id,
                "source_name": str(resource.get("name", Path(str(resource.get("path", ""))).stem)),
                "source_species": source_species,
                "source_version": str(resource.get("version", "")),
                "source_available_for_species": "true",
                "overlap_status": status,
                "distance_to_source_pas": distance,
                "distance_bin": _distance_bin(distance, status),
                "source_pas_id": source_pas_id,
                "source_region_class": source_region,
            }
        )
    return pd.DataFrame(rows, columns=EXTERNAL_PAS_COLUMNS)


def _nearest_source(group: pd.DataFrame | None, cut_site: int):
    if group is None or group.empty:
        return None
    positions = group["cut_site"].to_numpy(dtype=int)
    index = int(np.searchsorted(positions, cut_site, side="left"))
    candidates = [candidate for candidate in [index - 1, index] if 0 <= candidate < len(group)]
    best = min(candidates, key=lambda candidate: (abs(int(positions[candidate]) - cut_site), int(positions[candidate])))
    return group.iloc[best]


def _unavailable_rows(pas_map: pd.DataFrame, resource: dict, source_species: str) -> pd.DataFrame:
    source_name = str(resource.get("name", ""))
    source_version = str(resource.get("version", ""))
    return pd.DataFrame(
        [
            {
                "pas_id": row.pas_id,
                "source_name": source_name,
                "source_species": source_species,
                "source_version": source_version,
                "source_available_for_species": "false",
                "overlap_status": "unavailable",
                "distance_to_source_pas": "",
                "distance_bin": "unavailable",
                "source_pas_id": "",
                "source_region_class": "",
            }
            for row in pas_map.itertuples(index=False)
        ],
        columns=EXTERNAL_PAS_COLUMNS,
    )


def _column(resource: dict, key: str, table: pd.DataFrame, defaults: list[str], required: bool = True):
    configured = resource.get(key)
    if configured is not None:
        if isinstance(configured, int):
            try:
                return table.columns[int(configured)]
            except IndexError as exc:
                raise ValueError(f"External PAS {key} index out of range: {configured}") from exc
        if str(configured) not in table.columns:
            raise ValueError(f"External PAS column not found for {key}: {configured}")
        return str(configured)
    for candidate in defaults:
        if candidate in table.columns:
            return candidate
    if required:
        raise ValueError(f"External PAS resource is missing a required column for {key}: {defaults}")
    return None


def _chrom_key(value) -> str:
    chrom = str(value).strip()
    lowered = chrom.lower()
    if lowered.startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in {"M", "MT"}:
        return "MT"
    return chrom


def _distance_bin(distance, status: str) -> str:
    if status in {"unavailable", "unavailable_or_sparse"} and distance == "":
        return status
    if distance == "":
        return "unavailable"
    distance = int(distance)
    if distance == 0:
        return "0 bp"
    if distance <= 10:
        return "1-10 bp"
    if distance <= 25:
        return "11-25 bp"
    if distance <= 50:
        return "26-50 bp"
    if distance <= 100:
        return "51-100 bp"
    if distance <= 250:
        return "101-250 bp"
    return ">250 bp"
