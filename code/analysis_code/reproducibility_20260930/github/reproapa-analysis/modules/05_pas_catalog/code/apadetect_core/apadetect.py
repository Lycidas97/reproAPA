#!/usr/bin/env python3
"""Integrated apaDetect PAS calling, localization, and hard quantification."""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import os
import sys
import time
from array import array
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pysam
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

from endpoint_convolution import fit_kernel  # noqa: E402
from barcode_order_selection import (  # noqa: E402
    DEFAULT_MULTI_SELECTION_MODE,
    select_candidates as select_barcode_candidates,
    select_coverage_candidates,
)
from hard_interval_tool import (  # noqa: E402
    CompactMolecules,
    quantify_sweep,
    write_matrix,
)
from barcode_whitelist import load_barcode_whitelist  # noqa: E402
from candidate_filters import CandidateFilter, apply_candidate_filters
from multi_endpoint_mixture import (  # noqa: E402
    fit_mixture,
    oriented_coordinate,
    predicted_genomic,
    use_mixture_mle,
)
from conformity_decision_pipeline import select_order  # noqa: E402
from truncated_convolution import bin_masses  # noqa: E402
from junction_filter import (  # noqa: E402
    JunctionEvidence,
    junction_frame,
    junction_overlap_mask,
)
from native_event_backend import (  # noqa: E402
    collect_coverage_events_native,
    native_binary_sha256,
    resolve_event_backend,
)


BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]
NO_CANDIDATES_EXIT_CODE = 3


class NoCandidatesError(RuntimeError):
    """The chromosome has no PAS candidates that can enter quantification."""


@dataclass(frozen=True)
class CoverageComponent:
    start: int
    end: int
    strand: str
    candidate_indices: tuple[int, ...]


def read_strand(read: pysam.AlignedSegment, mode: str) -> str:
    strand = "-" if read.flag & 0x10 else "+"
    if mode == "opposite":
        return "+" if strand == "-" else "-"
    return strand


def usable_read(read: pysam.AlignedSegment) -> bool:
    return not read.flag & 0xD04


def estimate_read_length(bam_path: Path, limit: int = 10_000) -> int:
    lengths = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not usable_read(read):
                continue
            # Match the original apaDetect estimator: soft-clipped poly(A/T)
            # sequence is part of the sequenced read length.
            length = read.query_length or read.query_alignment_length or read.reference_length
            if length:
                lengths[int(length)] += 1
            if sum(lengths.values()) >= limit:
                break
    if not lengths:
        raise ValueError("No usable aligned reads were found for read-length estimation")
    return int(lengths.most_common(1)[0][0])


def triangular_smooth_reference(coverage: np.ndarray, read_length: int) -> np.ndarray:
    half = read_length // 2
    if half < 1:
        return coverage.astype(float, copy=True)
    kernel = np.concatenate((np.arange(1, half + 1), np.arange(half, 0, -1)))
    kernel = kernel / (half + 1)
    kernel = kernel / kernel.sum()
    return np.convolve(coverage, kernel, mode="same")


def _range_sum(prefix: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    values = np.zeros(len(starts), dtype=float)
    valid = ends >= starts
    values[valid] = prefix[ends[valid] + 1] - prefix[starts[valid]]
    return values


def triangular_smooth(coverage: np.ndarray, read_length: int) -> np.ndarray:
    """Exact O(n) triangular smoothing for large kernels, reference otherwise."""
    half = read_length // 2
    kernel_size = 2 * half
    if half < 1:
        return coverage.astype(float, copy=True)
    if coverage.size < kernel_size or kernel_size <= 512:
        return triangular_smooth_reference(coverage, read_length)

    values = coverage.astype(float, copy=False)
    positions = np.arange(coverage.size, dtype=float)
    prefix = np.r_[0.0, np.cumsum(values)]
    weighted = np.r_[0.0, np.cumsum(values * positions)]
    indices = np.arange(coverage.size, dtype=np.int64)
    conv_index = indices + (kernel_size - 1) // 2
    left_start = np.maximum(0, conv_index - kernel_size + 1)
    left_end = np.minimum(coverage.size - 1, conv_index - half)
    right_start = np.maximum(0, conv_index - half + 1)
    right_end = np.minimum(coverage.size - 1, conv_index)
    left = (
        (kernel_size - conv_index) * _range_sum(prefix, left_start, left_end)
        + _range_sum(weighted, left_start, left_end)
    )
    right = (
        (conv_index + 1) * _range_sum(prefix, right_start, right_end)
        - _range_sum(weighted, right_start, right_end)
    )
    return (left + right) / (half * (half + 1))


def collect_coverage_events(
    bam_path: Path,
    chromosome: str,
    strand_mode: str,
) -> tuple[
    int,
    dict[str, tuple[np.ndarray, np.ndarray]],
    JunctionEvidence,
]:
    positions = {"+": [], "-": []}
    deltas = {"+": [], "-": []}
    junctions = JunctionEvidence()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chrom_length = bam.get_reference_length(chromosome)
        for read in bam.fetch(chromosome):
            if not usable_read(read):
                continue
            strand = read_strand(read, strand_mode)
            blocks = read.get_blocks()
            if len(blocks) > 1:
                junctions.add_alignment(
                    read.cigartuples,
                    read.reference_start,
                    strand,
                )
            for start, end in blocks:
                start = max(0, int(start))
                end = min(chrom_length, int(end))
                if start < end:
                    positions[strand].extend((start, end))
                    deltas[strand].extend((1, -1))
    result = {}
    for strand in ("+", "-"):
        if not positions[strand]:
            result[strand] = (
                np.array([], dtype=np.int64),
                np.array([], dtype=np.int32),
            )
            continue
        pos = np.asarray(positions[strand], dtype=np.int64)
        delta = np.asarray(deltas[strand], dtype=np.int32)
        order = np.argsort(pos, kind="mergesort")
        pos, delta = pos[order], delta[order]
        unique, starts = np.unique(pos, return_index=True)
        result[strand] = (
            unique,
            np.add.reduceat(delta, starts).astype(np.int32, copy=False),
        )
    return chrom_length, result, junctions


def covered_intervals(
    positions: np.ndarray,
    deltas: np.ndarray,
    chrom_length: int,
    merge_gap: int = 0,
) -> list[tuple[int, int]]:
    intervals = []
    coverage = 0
    start = None
    for position, delta in zip(positions, deltas):
        previous = coverage
        coverage += int(delta)
        position = int(np.clip(position, 0, chrom_length))
        if previous == 0 and coverage > 0:
            start = position
        elif previous > 0 and coverage == 0 and start is not None:
            if start < position:
                intervals.append((start, position))
            start = None
    if start is not None:
        intervals.append((start, chrom_length))
    merged = []
    for start, end in intervals:
        if merged and start - merged[-1][1] <= merge_gap:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def coverage_from_events(
    positions: np.ndarray,
    deltas: np.ndarray,
    start: int,
    end: int,
) -> np.ndarray:
    diff = np.zeros(end - start + 1, dtype=np.int32)
    mask = (positions >= start) & (positions <= end)
    np.add.at(diff, positions[mask] - start, deltas[mask])
    coverage = diff[:-1]
    np.cumsum(coverage, out=coverage)
    return coverage


def call_peaks_from_coverage(
    coverage: np.ndarray,
    read_length: int,
    strand: str,
    chromosome: str,
    min_coverage: int,
    offset: int = 0,
) -> pd.DataFrame:
    if coverage.size < read_length or not coverage.size or coverage.max() < min_coverage:
        return pd.DataFrame(columns=BED_COLUMNS)
    oriented = coverage if strand == "+" else coverage[::-1]
    smoothed = triangular_smooth(oriented, read_length)
    derivative = np.diff(smoothed)
    valleys = np.where(np.diff(np.sign(derivative)) > 0)[0] + 1
    if not len(valleys):
        return pd.DataFrame(columns=BED_COLUMNS)
    mask = np.concatenate(
        (
            [False],
            (
                np.maximum.reduceat(smoothed, valleys[:-1])
                - smoothed[valleys[1:]]
            )
            >= min_coverage,
        )
    )
    downstream = valleys[mask]
    if not len(downstream):
        return pd.DataFrame(columns=BED_COLUMNS)
    if strand == "-":
        downstream = coverage.size - downstream - 1
        valleys = coverage.size - valleys - 1
    downstream.sort()
    valleys.sort()
    if strand == "+":
        upstream = valleys[np.searchsorted(valleys, downstream, side="left") - 1]
        starts, ends = upstream, downstream
    else:
        upstream = valleys[np.searchsorted(valleys, downstream, side="right")]
        starts, ends = downstream, upstream
    frame = pd.DataFrame(
        {
            "chrom": chromosome,
            "start": starts + offset,
            "end": ends + offset,
            "score": downstream + offset,
            "strand": strand,
        }
    )
    frame["name"] = (
        frame["chrom"].astype(str)
        + ":"
        + frame["start"].astype(str)
        + "-"
        + frame["end"].astype(str)
        + ":"
        + frame["strand"]
    )
    return frame[BED_COLUMNS]


def call_candidates(
    chromosome: str,
    chrom_length: int,
    events: dict[str, tuple[np.ndarray, np.ndarray]],
    read_length: int,
    min_coverage: int,
) -> tuple[pd.DataFrame, dict[str, list[tuple[int, int]]]]:
    frames = []
    components = {}
    padding = 2 * max(1, read_length // 2)
    for strand in ("+", "-"):
        positions, deltas = events[strand]
        components[strand] = covered_intervals(positions, deltas, chrom_length)
        islands = covered_intervals(
            positions, deltas, chrom_length, merge_gap=2 * padding
        )
        for core_start, core_end in islands:
            start = max(0, core_start - padding)
            end = min(chrom_length, core_end + padding)
            coverage = coverage_from_events(positions, deltas, start, end)
            frame = call_peaks_from_coverage(
                coverage,
                read_length,
                strand,
                chromosome,
                min_coverage,
                offset=start,
            )
            if not frame.empty:
                frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=BED_COLUMNS), components
    candidates = (
        pd.concat(frames, ignore_index=True)
        .drop_duplicates(["chrom", "start", "end", "strand"])
        .sort_values(["chrom", "strand", "start", "end"])
        .reset_index(drop=True)
    )
    return candidates[BED_COLUMNS], components


def fixed_coordinates(candidates: pd.DataFrame, read_length: int) -> pd.Series:
    half = read_length // 2
    return candidates["score"].astype(int) + np.where(
        candidates["strand"].eq("+"), -half, half
    )


def canonical_chromosome(chromosome: str) -> str:
    name = str(chromosome)
    return name if name.startswith("chr") else f"chr{name}"


def load_asnr_intervals(
    asnr_bed: Path,
) -> dict[tuple[str, str], tuple[list[int], list[int]]]:
    """Load and merge ASNR intervals for exact same-strand overlap queries."""
    if not asnr_bed.exists():
        raise FileNotFoundError(f"ASNR BED does not exist: {asnr_bed}")
    grouped = defaultdict(list)
    with asnr_bed.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"Invalid ASNR BED record at line {line_number}")
            chrom, start, end, _, _, strand = fields[:6]
            if strand not in {"+", "-"}:
                raise ValueError(f"Invalid ASNR strand at line {line_number}: {strand}")
            grouped[(canonical_chromosome(chrom), strand)].append(
                (int(start), int(end))
            )

    result = {}
    for key, values in grouped.items():
        merged = []
        for start, end in sorted(values):
            if start >= end:
                raise ValueError(f"ASNR BED contains a non-positive interval: {key} {start} {end}")
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        result[key] = (
            [start for start, _ in merged],
            [end for _, end in merged],
        )
    return result


def asnr_overlap_mask(candidates: pd.DataFrame, asnr_bed: Path) -> pd.Series:
    """Match ``bedtools subtract -s -A`` removal decisions."""
    intervals = load_asnr_intervals(asnr_bed)
    removed = []
    for row in candidates.itertuples():
        starts, ends = intervals.get(
            (canonical_chromosome(row.chrom), row.strand),
            ((), ()),
        )
        index = bisect.bisect_left(starts, int(row.end)) - 1
        removed.append(index >= 0 and ends[index] > int(row.start))
    return pd.Series(removed, index=candidates.index, dtype=bool)


def assign_components(
    candidates: pd.DataFrame,
    component_intervals: dict[str, list[tuple[int, int]]],
) -> list[CoverageComponent]:
    result = []
    for strand in ("+", "-"):
        intervals = component_intervals[strand]
        groups = [[] for _ in intervals]
        starts = [start for start, _ in intervals]
        ends = [end for _, end in intervals]
        current = candidates[
            candidates["strand"].eq(strand)
            & ~candidates["asnr_filtered"]
            & ~candidates["junction_filtered"]
        ]
        for row in current.itertuples():
            left, right = int(row.start), int(row.end)
            interval_index = bisect.bisect_right(ends, left)
            best_index = None
            best_overlap = 0
            while interval_index < len(intervals) and starts[interval_index] < right:
                start, end = intervals[interval_index]
                overlap = max(0, min(right, end) - max(left, start))
                if overlap > best_overlap:
                    best_index = interval_index
                    best_overlap = overlap
                interval_index += 1
            if best_index is not None:
                groups[best_index].append(int(row.Index))
        result.extend(
            CoverageComponent(start, end, strand, tuple(groups[index]))
            for index, (start, end) in enumerate(intervals)
            if groups[index]
        )
    return result


def collect_component_data(
    bam_path: Path,
    chromosome: str,
    components: list[CoverageComponent],
    strand_mode: str,
    barcode_tag: str,
    umi_tag: str,
    valid_barcodes: set[str] | None = None,
) -> tuple[
    dict[int, np.ndarray],
    dict[int, list[tuple[str, int]]],
    dict[str, int],
]:
    by_strand = {
        strand: [(index, item) for index, item in enumerate(components) if item.strand == strand]
        for strand in ("+", "-")
    }
    starts = {
        strand: [component.start for _, component in values]
        for strand, values in by_strand.items()
    }
    endpoint_molecules = defaultdict(dict)
    barcode_molecules = defaultdict(dict)
    stats = Counter()
    read_counter = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chromosome):
            if not usable_read(read):
                continue
            strand = read_strand(read, strand_mode)
            values = by_strand[strand]
            if not values:
                continue
            endpoint = int(read.reference_end if strand == "+" else read.reference_start)
            quantification_endpoint = int(
                read.reference_end - 1 if strand == "+" else read.reference_start
            )
            local = bisect.bisect_right(starts[strand], endpoint) - 1
            for neighbor in (local - 1, local, local + 1):
                if not 0 <= neighbor < len(values):
                    continue
                component_index, component = values[neighbor]
                if not (
                    component.start <= endpoint <= component.end
                    and read.reference_start < component.end
                    and read.reference_end > component.start
                ):
                    continue
                try:
                    barcode = str(read.get_tag(barcode_tag))
                    umi = str(read.get_tag(umi_tag))
                    if (
                        valid_barcodes is not None
                        and barcode not in valid_barcodes
                    ):
                        stats["outside_barcode_whitelist_reads"] += 1
                        continue
                    molecule = (barcode, umi)
                    tagged = True
                except KeyError:
                    if valid_barcodes is not None:
                        stats["missing_tag_reads"] += 1
                        continue
                    molecule = (str(read.query_name), str(read_counter))
                    tagged = False
                oriented = oriented_coordinate(endpoint, strand)
                previous = endpoint_molecules[component_index].get(molecule)
                if previous is None or oriented > previous:
                    endpoint_molecules[component_index][molecule] = oriented
                if tagged:
                    oriented_quantification = oriented_coordinate(
                        quantification_endpoint, strand
                    )
                    previous = barcode_molecules[component_index].get(molecule)
                    if previous is None or oriented_quantification > previous:
                        barcode_molecules[component_index][molecule] = (
                            oriented_quantification
                        )
                    stats["tagged_reads"] += 1
                else:
                    stats["missing_tag_reads"] += 1
            read_counter += 1
    endpoints = {
        index: np.asarray(list(values.values()), dtype=int)
        for index, values in endpoint_molecules.items()
    }
    tagged = {
        index: [(barcode, endpoint) for (barcode, _), endpoint in values.items()]
        for index, values in barcode_molecules.items()
    }
    stats["tagged_molecules"] = sum(len(values) for values in tagged.values())
    return endpoints, tagged, dict(stats)


def collect_component_endpoints(
    bam_path: Path,
    chromosome: str,
    components: list[CoverageComponent],
    strand_mode: str,
    barcode_tag: str,
    umi_tag: str,
) -> dict[int, np.ndarray]:
    endpoints, _, _ = collect_component_data(
        bam_path,
        chromosome,
        components,
        strand_mode,
        barcode_tag,
        umi_tag,
    )
    return endpoints


def terminal_endpoint(
    component: CoverageComponent,
    endpoints: np.ndarray,
    cluster_width: int,
    min_molecules: int,
) -> tuple[float | None, int]:
    genomic_boundary = component.end if component.strand == "+" else component.start
    oriented_boundary = oriented_coordinate(genomic_boundary, component.strand)
    support = int(
        np.sum(
            (endpoints >= oriented_boundary - cluster_width)
            & (endpoints <= oriented_boundary)
        )
    )
    return (
        (float(genomic_boundary), support)
        if support >= min_molecules
        else (None, support)
    )


def train_kernel(
    candidates: pd.DataFrame,
    components: list[CoverageComponent],
    endpoints_by_component: dict[int, np.ndarray],
    read_length: int,
    min_molecules: int,
) -> dict | None:
    lower, upper = -read_length, read_length // 2
    profiles = []
    for component_index, component in enumerate(components):
        if len(component.candidate_indices) != 1:
            continue
        endpoints = endpoints_by_component.get(component_index, np.array([], dtype=int))
        if len(endpoints) < min_molecules:
            continue
        index = component.candidate_indices[0]
        fixed = oriented_coordinate(
            int(candidates.loc[index, "fixed_coordinate"]), component.strand
        )
        offsets = endpoints - fixed
        offsets = offsets[(offsets >= lower) & (offsets <= upper)]
        if len(offsets) < min_molecules:
            continue
        counts = np.bincount(offsets - lower, minlength=upper - lower + 1)
        profiles.append(counts / counts.sum())
    if len(profiles) < 20:
        return None
    coordinates = np.arange(lower, upper + 1)
    fitted = fit_kernel(
        coordinates,
        np.mean(profiles, axis=0),
        read_length,
        bin_masses,
    )
    fitted.pop("probabilities")
    return {
        "model": "TN[-L,0](mu,sigma1^2)+N(0,sigma2^2)",
        "training_role": "sample_internal_clean_single_components",
        "read_length": read_length,
        "window_lower": lower,
        "window_upper": upper,
        "training_sites": len(profiles),
        "endpoint_location_calibration": 0.0,
        "fixed_offset": 0.0,
        **fitted,
    }


def localize_components(
    candidates: pd.DataFrame,
    components: list[CoverageComponent],
    endpoints_by_component: dict[int, np.ndarray],
    kernel: dict | None,
    bic_threshold: float,
    min_model_molecules: int,
    barcode_molecules_by_component: dict[int, list[tuple[str, int]]],
    barcode_bic_threshold: float,
    barcode_min_informative: int,
    multi_selection_mode: str,
    endpoint_cluster_width: int,
    endpoint_min_molecules: int,
    coverage_by_component: dict[int, tuple[np.ndarray, np.ndarray]] | None = None,
    read_length: int = 1,
    coverage_max_gap_read_lengths: float = 2.0,
    coverage_valley_ratio: float = 0.5,
) -> pd.DataFrame:
    rows = []
    for component_index, component in enumerate(components):
        original = list(component.candidate_indices)
        endpoints = endpoints_by_component.get(component_index, np.array([], dtype=int))
        fitted_coordinates = {}
        bic_history = []
        barcode_history = []
        if len(original) == 1:
            selected = original
            index = original[0]
            fixed = float(candidates.loc[index, "fixed_coordinate"])
            endpoint, endpoint_support = terminal_endpoint(
                component, endpoints, endpoint_cluster_width, endpoint_min_molecules
            )
            fixed_oriented = oriented_coordinate(int(round(fixed)), component.strand)
            endpoint_oriented = (
                oriented_coordinate(int(round(endpoint)), component.strand)
                if endpoint is not None
                else None
            )
            decision_coordinate = fixed
            decision = "fixed_single"
            if endpoint_oriented is not None and endpoint_oriented > fixed_oriented:
                decision_coordinate = float(endpoint)
                decision = "downstream_endpoint_single"
        else:
            endpoint_support = 0
            if multi_selection_mode in {"barcode", "coverage"}:
                oriented_fixed = np.asarray(
                    [
                        oriented_coordinate(
                            int(candidates.loc[index, "fixed_coordinate"]),
                            component.strand,
                        )
                        for index in original
                    ],
                    dtype=float,
                )
                barcode_selected, barcode_history = select_barcode_candidates(
                    original,
                    oriented_fixed,
                    barcode_molecules_by_component.get(component_index, []),
                    barcode_bic_threshold,
                    barcode_min_informative,
                )
                if multi_selection_mode == "barcode":
                    selected = barcode_selected
                else:
                    if coverage_by_component is None:
                        raise ValueError("coverage selection requires component coverage")
                    positions, smoothed_coverage = coverage_by_component[
                        component_index
                    ]
                    selected, bic_history = select_coverage_candidates(
                        original,
                        oriented_fixed,
                        positions,
                        smoothed_coverage,
                        read_length,
                        coverage_max_gap_read_lengths,
                        coverage_valley_ratio,
                        barcode_selected,
                    )
            elif kernel is None or len(endpoints) < min_model_molecules * len(original):
                selected = original
                decision = "fixed_multi_low_support"
            else:
                selected, bic_history = select_order(
                    endpoints,
                    original,
                    candidates.rename(
                        columns={"fixed_coordinate": "corrected_coordinate"}
                    ),
                    component.strand,
                    kernel,
                    bic_threshold,
                )
                decision = (
                    "mixture_selected"
                    if len(selected) > 1
                    else "single_selected_by_bic"
                )
                p0 = np.asarray(
                    [
                        oriented_coordinate(
                            int(candidates.loc[index, "fixed_coordinate"]),
                            component.strand,
                        )
                        for index in selected
                    ]
                )
                final_fit = (
                    fit_mixture(endpoints, p0, kernel)
                    if len(selected) > 1 and use_mixture_mle(p0, 300.0)
                    else None
                )
                if final_fit is not None:
                    taus = final_fit["raw_taus"] - float(
                        kernel["endpoint_location_calibration"]
                    )
                    for position, index in enumerate(selected):
                        fixed = int(candidates.loc[index, "fixed_coordinate"])
                        fitted_coordinates[index] = predicted_genomic(
                            fixed, taus[position], component.strand
                        )
        evidence_by_candidate = {
            int(item["tested_candidate_index"]): item
            for item in barcode_history
        }
        downstream_index = max(
            original,
            key=lambda index: oriented_coordinate(
                int(candidates.loc[index, "fixed_coordinate"]),
                component.strand,
            ),
        )
        for index in original:
            fixed = float(candidates.loc[index, "fixed_coordinate"])
            if len(original) > 1:
                decision_coordinate = (
                    fitted_coordinates.get(index, fixed)
                    if index in selected
                    else np.nan
                )
                if multi_selection_mode in {"barcode", "coverage"} and index == downstream_index:
                    decision = "downstream_anchor"
                elif multi_selection_mode == "coverage":
                    decision = next(
                        item["reason"] for item in bic_history
                        if item["tested_candidate_index"] == index
                    )
                elif multi_selection_mode == "barcode" and index in selected:
                    decision = "barcode_confirmed"
                elif multi_selection_mode == "barcode":
                    decision = "ambiguous_merged"
            evidence = evidence_by_candidate.get(index, {})
            rows.append(
                {
                    "candidate_index": index,
                    "component_index": component_index,
                    "component_start": component.start,
                    "component_end": component.end,
                    "component_candidates": len(original),
                    "component_molecules": len(endpoints),
                    "selected_order": len(selected),
                    "selected": index in selected,
                    "decision": decision,
                    "endpoint_terminal_support": endpoint_support,
                    "fixed_coordinate": fixed,
                    "decision_coordinate": decision_coordinate,
                    "bic_history_json": json.dumps(
                        bic_history,
                        sort_keys=True,
                        default=lambda value: value.item(),
                    ),
                    "barcode_history_json": json.dumps(
                        barcode_history,
                        sort_keys=True,
                        default=lambda value: value.item(),
                    ),
                    "barcode_delta_bic": evidence.get(
                        "beta_binomial_delta_bic", np.nan
                    ),
                    "barcode_rho": evidence.get("beta_binomial_rho", np.nan),
                    "barcode_count": evidence.get("barcode_count", np.nan),
                    "barcode_informative_count": evidence.get(
                        "informative_barcode_count", np.nan
                    ),
                    "barcode_evidence_accepted": bool(
                        evidence.get("accepted", False)
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_quantification_intervals(
    calls: pd.DataFrame,
    chromosome: str,
) -> pd.DataFrame:
    selected = calls[calls["selected"]].copy()
    rows = []
    pas_number = 1
    for (_, strand), group in selected.groupby(["component_index", "strand"], sort=False):
        component_start = int(group["component_start"].iloc[0])
        component_end = int(group["component_end"].iloc[0])
        ordered = group.sort_values("decision_coordinate")
        coordinates = ordered["decision_coordinate"].round().astype(int).to_numpy()
        boundaries = (
            ((coordinates[:-1] + coordinates[1:]) // 2).astype(int)
            if len(coordinates) > 1
            else np.array([], dtype=int)
        )
        starts = np.r_[component_start, boundaries]
        ends = np.r_[boundaries, component_end]
        for row, start, end in zip(ordered.itertuples(), starts, ends):
            if start >= end:
                raise ValueError("PAS localization produced a non-positive quantification interval")
            rows.append(
                {
                    "chrom": chromosome,
                    "start": int(start),
                    "end": int(end),
                    "name": f"PAS|{chromosome}|{pas_number:08d}",
                    "score": int(round(row.decision_coordinate)),
                    "strand": strand,
                    "candidate_index": int(row.candidate_index),
                }
            )
            pas_number += 1
    return pd.DataFrame(rows).sort_values(
        ["chrom", "strand", "start", "end"]
    ).reset_index(drop=True)


def write_matrix_optimized(
    counts,
    output_dir: Path,
    *,
    write_matrix_file: bool = True,
    barcode_order: list[str] | None = None,
) -> dict:
    """Write the same logical outputs without sorting every UMI string."""
    if not isinstance(counts, CompactMolecules):
        if not write_matrix_file:
            raise TypeError("Deferred matrix output requires CompactMolecules")
        return write_matrix(counts, output_dir)
    pair_counts = Counter(
        (peak_id, barcode_id) for peak_id, barcode_id, _ in counts.molecules
    )
    used_peak_ids = {peak_id for peak_id, _ in pair_counts}
    used_barcode_ids = {barcode_id for _, barcode_id in pair_counts}
    peak_ids = sorted(used_peak_ids, key=lambda value: counts.peaks[value])
    if barcode_order is None:
        barcode_ids = sorted(
            used_barcode_ids, key=lambda value: counts.barcodes[value]
        )
        barcode_names = [counts.barcodes[value] for value in barcode_ids]
    else:
        allowed_barcodes = set(barcode_order)
        unexpected = {
            counts.barcodes[value]
            for value in used_barcode_ids
            if counts.barcodes[value] not in allowed_barcodes
        }
        if unexpected:
            raise ValueError(
                "Counts contain barcodes outside the supplied whitelist"
            )
        barcode_ids = [
            counts.barcode_ids[barcode]
            for barcode in barcode_order
            if barcode in counts.barcode_ids
        ]
        barcode_names = list(barcode_order)
    peak_rank = {value: index for index, value in enumerate(peak_ids)}
    barcode_rank = (
        {value: index for index, value in enumerate(barcode_ids)}
        if barcode_order is None
        else {
            counts.barcode_ids[barcode]: index
            for index, barcode in enumerate(barcode_order)
            if barcode in counts.barcode_ids
        }
    )
    ordered_pairs = sorted(
        pair_counts,
        key=lambda value: (peak_rank[value[0]], barcode_rank[value[1]]),
    )
    rows, columns, values = array("i"), array("i"), array("i")
    with gzip.open(output_dir / "counts.tsv.gz", "wt", compresslevel=1) as handle:
        handle.write("pas_id\tbarcode\tumi_count\n")
        for peak_id, barcode_id in ordered_pairs:
            count = pair_counts[(peak_id, barcode_id)]
            if write_matrix_file:
                rows.append(peak_rank[peak_id])
                columns.append(barcode_rank[barcode_id])
                values.append(count)
            handle.write(
                f"{counts.peaks[peak_id]}\t{counts.barcodes[barcode_id]}"
                f"\t{count}\n"
            )
    matrix_path = output_dir / "matrix.mtx"
    if write_matrix_file:
        matrix = coo_matrix(
            (
                np.frombuffer(values, dtype=np.int32),
                (
                    np.frombuffer(rows, dtype=np.int32),
                    np.frombuffer(columns, dtype=np.int32),
                ),
            ),
            shape=(len(peak_ids), len(barcode_names)),
            dtype=np.int32,
        )
        mmwrite(matrix_path, matrix)
    else:
        matrix_path.unlink(missing_ok=True)
    peak_names = [counts.peaks[value] for value in peak_ids]
    (output_dir / "pas.tsv").write_text(
        "\n".join(peak_names) + ("\n" if peak_names else "")
    )
    (output_dir / "barcodes.tsv").write_text(
        "\n".join(barcode_names) + ("\n" if barcode_names else "")
    )
    return {
        "quantified_pas": len(peak_ids),
        "barcodes": len(barcode_names),
        "nonzero_entries": len(ordered_pairs),
        "unique_molecules": len(counts.molecules),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--asnr-bed", type=Path, required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--read-length", type=int)
    parser.add_argument("--min-coverage", type=int, default=20)
    parser.add_argument(
        "--event-backend",
        choices=["auto", "native", "python"],
        default="auto",
    )
    parser.add_argument("--native-events-path", type=Path)
    parser.add_argument("--min-junction-count", type=int, default=20)
    parser.add_argument("--junction-min-anchor", type=int, default=8)
    parser.add_argument("--junction-window-left", type=int, default=10)
    parser.add_argument("--junction-window-right", type=int, default=100)
    parser.add_argument("--strand-mode", choices=["same", "opposite"], default="same")
    parser.add_argument("--barcode-tag", default="CB")
    parser.add_argument("--umi-tag", default="UB")
    parser.add_argument(
        "--barcode-whitelist",
        type=Path,
        help="One barcode per line in the exact expression-matrix column order",
    )
    parser.add_argument("--endpoint-cluster-width", type=int, default=5)
    parser.add_argument("--endpoint-min-molecules", type=int, default=3)
    parser.add_argument("--model-min-molecules", type=int, default=20)
    parser.add_argument("--bic-threshold", type=float, default=1500.2874445228445)
    parser.add_argument("--barcode-bic-threshold", type=float, default=10.0)
    parser.add_argument("--barcode-min-informative", type=int, default=10)
    parser.add_argument("--coverage-max-gap-read-lengths", type=float, default=2.0)
    parser.add_argument("--coverage-valley-ratio", type=float, default=0.5)
    parser.add_argument(
        "--multi-selection-mode",
        choices=["coverage", "mixture", "barcode"],
        default=DEFAULT_MULTI_SELECTION_MODE,
        help="Coverage-led conservative selection (default), legacy barcode selection, or experimental mixture fitting",
    )
    parser.add_argument(
        "--skip-matrix",
        action="store_true",
        help="Write count-table shards but defer Matrix Market output",
    )
    parser.add_argument(
        "--kernel",
        type=Path,
        help="Reuse one sample-internal kernel learned by a genome-level run",
    )
    args = parser.parse_args()

    if not args.bam.exists() or not args.fasta.exists() or not args.asnr_bed.exists():
        raise FileNotFoundError("BAM, FASTA, and ASNR BED inputs must exist")
    if args.barcode_whitelist and not args.barcode_whitelist.exists():
        raise FileNotFoundError("Barcode whitelist does not exist")
    barcode_order = (
        load_barcode_whitelist(args.barcode_whitelist)
        if args.barcode_whitelist
        else None
    )
    valid_barcodes = set(barcode_order) if barcode_order is not None else None
    if args.kernel and args.multi_selection_mode != "mixture":
        parser.error(
            "--kernel is only valid with --multi-selection-mode mixture"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    phase_started = started
    timings = {}
    read_length = args.read_length or estimate_read_length(args.bam)
    timings["read_length_seconds"] = time.perf_counter() - phase_started

    phase_started = time.perf_counter()
    event_backend, native_events_path = resolve_event_backend(
        args.event_backend,
        args.native_events_path,
    )
    if event_backend == "native":
        chrom_length, events, junction_evidence = collect_coverage_events_native(
            native_events_path,
            args.bam,
            args.chromosome,
            args.strand_mode,
        )
    else:
        chrom_length, events, junction_evidence = collect_coverage_events(
            args.bam, args.chromosome, args.strand_mode
        )
    candidates, component_intervals = call_candidates(
        args.chromosome,
        chrom_length,
        events,
        read_length,
        args.min_coverage,
    )
    timings["candidate_generation_seconds"] = time.perf_counter() - phase_started
    if candidates.empty:
        raise NoCandidatesError(
            "No derivative PAS candidates passed the coverage threshold"
        )
    candidates[BED_COLUMNS].to_csv(
        args.output_dir / "raw_candidates.bed", sep="\t", header=False, index=False
    )
    phase_started = time.perf_counter()
    fixed_coordinate = fixed_coordinates(candidates, read_length)
    timings["asnr_filter_seconds"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    junction_points = junction_evidence.points(
        args.min_junction_count,
        args.junction_min_anchor,
    )
    candidates, retained_mask = apply_candidate_filters(
        candidates,
        [
            CandidateFilter(
                "asnr",
                lambda frame: asnr_overlap_mask(frame, args.asnr_bed),
            ),
            CandidateFilter(
                "junction",
                lambda frame: junction_overlap_mask(
                    frame,
                    junction_points,
                    args.junction_window_left,
                    args.junction_window_right,
                ),
            ),
        ],
    )
    junction_filtered = candidates.pop("junction_filtered")
    candidates["fixed_coordinate"] = fixed_coordinate
    candidates["junction_filtered"] = junction_filtered
    timings["junction_filter_seconds"] = time.perf_counter() - phase_started
    junction_frame(args.chromosome, junction_points).to_csv(
        args.output_dir / "junctions.bed",
        sep="\t",
        header=False,
        index=False,
    )
    candidates.reset_index(drop=True).to_csv(
        args.output_dir / "candidate_qc.tsv", sep="\t", index_label="candidate_index"
    )
    candidates.loc[~candidates["asnr_filtered"], BED_COLUMNS].to_csv(
        args.output_dir / "asnr_filtered_candidates.bed",
        sep="\t",
        header=False,
        index=False,
    )
    candidates.loc[retained_mask, BED_COLUMNS].to_csv(
        args.output_dir / "filtered_candidates.bed",
        sep="\t",
        header=False,
        index=False,
    )
    components = assign_components(candidates, component_intervals)
    if not components:
        raise NoCandidatesError(
            "All derivative candidates were removed by ASNR/junction filtering"
        )
    phase_started = time.perf_counter()
    endpoints, barcode_molecules, barcode_collection_stats = collect_component_data(
        args.bam,
        args.chromosome,
        components,
        args.strand_mode,
        args.barcode_tag,
        args.umi_tag,
        valid_barcodes,
    )
    timings["endpoint_collection_seconds"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    if args.multi_selection_mode != "mixture":
        kernel = None
        kernel_source = "not_used_candidate_selection"
    elif args.kernel:
        kernel = json.loads(args.kernel.read_text())
        if int(kernel["read_length"]) != read_length:
            raise ValueError(
                f"Supplied kernel read length {kernel['read_length']} "
                f"does not match sample read length {read_length}"
            )
        kernel_source = "supplied_sample_kernel"
    else:
        kernel = train_kernel(
            candidates,
            components,
            endpoints,
            read_length,
            args.model_min_molecules,
        )
        kernel_source = (
            "trained_on_current_chromosome"
            if kernel is not None
            else "insufficient_clean_single_components"
        )
    timings["sample_kernel_seconds"] = time.perf_counter() - phase_started
    if kernel is not None:
        (args.output_dir / "sample_kernel.json").write_text(
            json.dumps(kernel, indent=2) + "\n"
        )
    phase_started = time.perf_counter()
    coverage_by_component = {}
    for component_index, component in enumerate(components):
        positions, deltas = events[component.strand]
        coverage = coverage_from_events(
            positions, deltas, component.start, component.end
        )
        genomic = np.arange(component.start, component.end, dtype=float)
        if component.strand == "+":
            oriented_positions, oriented_coverage = genomic, coverage
        else:
            oriented_positions, oriented_coverage = -genomic[::-1], coverage[::-1]
        coverage_by_component[component_index] = (
            oriented_positions,
            triangular_smooth(oriented_coverage, read_length),
        )
    calls = localize_components(
        candidates,
        components,
        endpoints,
        kernel,
        args.bic_threshold,
        args.model_min_molecules,
        barcode_molecules,
        args.barcode_bic_threshold,
        args.barcode_min_informative,
        args.multi_selection_mode,
        args.endpoint_cluster_width,
        args.endpoint_min_molecules,
        coverage_by_component,
        read_length,
        args.coverage_max_gap_read_lengths,
        args.coverage_valley_ratio,
    ).merge(
        candidates.drop(columns=["fixed_coordinate"]).reset_index(
            names="candidate_index"
        ),
        on="candidate_index",
        how="left",
    )
    timings["localization_seconds"] = time.perf_counter() - phase_started
    calls.to_csv(args.output_dir / "pas_calls.tsv", sep="\t", index=False)
    intervals = build_quantification_intervals(calls, args.chromosome)
    intervals.to_csv(
        args.output_dir / "pas_intervals.tsv", sep="\t", index=False
    )
    intervals[BED_COLUMNS].to_csv(
        args.output_dir / "pas_intervals.bed", sep="\t", header=False, index=False
    )
    phase_started = time.perf_counter()
    counts, scan_stats = quantify_sweep(
        args.bam,
        intervals[BED_COLUMNS],
        args.barcode_tag,
        args.umi_tag,
        compact=True,
        valid_barcodes=valid_barcodes,
    )
    timings["quantification_scan_seconds"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    matrix_stats = write_matrix_optimized(
        counts,
        args.output_dir,
        write_matrix_file=not args.skip_matrix,
        barcode_order=barcode_order,
    )
    timings["matrix_write_seconds"] = time.perf_counter() - phase_started
    summary = {
        "bam": str(args.bam),
        "fasta": str(args.fasta),
        "asnr_bed": str(args.asnr_bed),
        "chromosome": args.chromosome,
        "read_length": read_length,
        "event_backend": event_backend,
        "native_events_path": (
            str(native_events_path) if event_backend == "native" else None
        ),
        "native_events_sha256": (
            native_binary_sha256(native_events_path)
            if event_backend == "native"
            else None
        ),
        "raw_candidates": int(len(candidates)),
        "asnr_removed": int(candidates["asnr_filtered"].sum()),
        "junctions_retained": sum(map(len, junction_points.values())),
        "junction_removed": int(
            ((~candidates["asnr_filtered"]) & candidates["junction_filtered"]).sum()
        ),
        "filtered_candidates": int(retained_mask.sum()),
        "coverage_components": len(components),
        "sample_kernel_status": (
            kernel_source
        ),
        "sample_kernel_training_sites": int(kernel["training_sites"]) if kernel else 0,
        "selected_pas": int(calls["selected"].sum()),
        "barcode_merged": int(calls["decision"].eq("ambiguous_merged").sum()),
        "barcode_confirmed": int(calls["decision"].eq("barcode_confirmed").sum()),
        "coverage_shoulder_merged": int(
            calls["decision"].eq("shoulder_merged").sum()
        ),
        "coverage_distance_retained": int(
            calls["decision"].eq("distance_retained").sum()
        ),
        "coverage_valley_retained": int(
            calls["decision"].eq("valley_retained").sum()
        ),
        "coverage_barcode_rescued": int(
            calls["decision"].eq("barcode_rescued").sum()
        ),
        "coverage_max_gap_read_lengths": args.coverage_max_gap_read_lengths,
        "coverage_valley_ratio": args.coverage_valley_ratio,
        "barcode_bic_threshold": args.barcode_bic_threshold,
        "barcode_min_informative": args.barcode_min_informative,
        "multi_selection_mode": args.multi_selection_mode,
        "barcode_collection_stats": barcode_collection_stats,
        "barcode_whitelist": (
            str(args.barcode_whitelist) if args.barcode_whitelist else None
        ),
        "barcode_whitelist_size": (
            len(barcode_order) if barcode_order is not None else None
        ),
        "endpoint_single_calls": int(
            calls["decision"].eq("downstream_endpoint_single").sum()
        ),
        "timings": timings,
        "scan_stats": scan_stats,
        "matrix_written": not args.skip_matrix,
        **matrix_stats,
        "elapsed_seconds": time.perf_counter() - started,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except NoCandidatesError as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(NO_CANDIDATES_EXIT_CODE)
