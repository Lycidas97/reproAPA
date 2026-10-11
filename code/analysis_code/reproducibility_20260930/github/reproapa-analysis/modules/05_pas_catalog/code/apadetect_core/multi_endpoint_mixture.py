#!/usr/bin/env python3
"""Fixed-candidate truncated-convolution mixture for overlapping PAS peaks."""

from __future__ import annotations

import argparse
import bisect
import heapq
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.optimize import linear_sum_assignment
from scipy.special import logsumexp

from coverage_components import CoverageComponent, assign_candidates  # noqa: E402
from endpoint_convolution import (  # noqa: E402
    Kernel,
    endpoint_probabilities,
    estimate_tau,
)
from truncated_convolution import bin_masses  # noqa: E402


def frozen_calibration(values: np.ndarray, frozen: dict) -> np.ndarray:
    x = np.asarray(frozen["x_thresholds"], dtype=float)
    y = np.asarray(frozen["y_thresholds"], dtype=float)
    isotonic = np.interp(values, x, y, left=y[0], right=y[-1])
    target = float(frozen["gate_target"])
    threshold = float(frozen["update_threshold"])
    return np.where(np.abs(isotonic - target) >= threshold, isotonic, target)


def model_fixed_offset(kernel: dict, frozen: dict) -> float:
    """Return the frozen sample offset and reject inconsistent model files."""
    target = float(frozen["gate_target"])
    stored = float(kernel.get("fixed_offset", target))
    if not np.isclose(stored, target):
        raise ValueError(
            f"kernel fixed_offset {stored} differs from frozen gate_target {target}"
        )
    return target


def use_mixture_mle(p0_values: np.ndarray, max_median_gap: float) -> bool:
    """Use mixture refinement only for candidate groups shown resolvable in development."""
    gaps = np.diff(np.sort(np.unique(np.asarray(p0_values, dtype=float))))
    return len(gaps) > 0 and float(np.median(gaps)) <= max_median_gap


def oriented_coordinate(genomic: int, strand: str) -> int:
    return int(genomic if strand == "+" else -genomic)


def predicted_genomic(p0: int, tau: float, strand: str) -> float:
    return float(p0 + tau if strand == "+" else p0 - tau)


def component_endpoints(
    bam: pysam.AlignmentFile,
    chromosome: str,
    start: int,
    end: int,
    strand: str,
) -> np.ndarray:
    endpoints = []
    for read in bam.fetch(chromosome, start, end):
        read_strand = "-" if read.is_reverse else "+"
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or read_strand != strand
            or not read.cigarstring
            or re.search(r"\d+N", read.cigarstring)
        ):
            continue
        endpoint = read.reference_end if strand == "+" else read.reference_start
        if start <= endpoint <= end:
            endpoints.append(oriented_coordinate(endpoint, strand))
    return np.asarray(endpoints, dtype=int)


def stream_positive_intervals(
    bam_path: str, chromosome: str
) -> tuple[int, dict[str, list[tuple[int, int]]]]:
    """Return exact positive-coverage unions without retaining per-read events."""
    heaps = {"+": [], "-": []}
    merged = {"+": [], "-": []}
    current = {"+": None, "-": None}

    def consume(strand: str, start: int, end: int) -> None:
        active = current[strand]
        if active is None:
            current[strand] = [start, end]
        elif start <= active[1]:
            active[1] = max(active[1], end)
        else:
            merged[strand].append((active[0], active[1]))
            current[strand] = [start, end]

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        chromosome_length = bam.get_reference_length(chromosome)
        for read in bam.fetch(chromosome):
            read_start = int(read.reference_start)
            for strand in ("+", "-"):
                while heaps[strand] and heaps[strand][0][0] <= read_start:
                    consume(strand, *heapq.heappop(heaps[strand]))
            if read.is_unmapped:
                continue
            if read.is_reverse:
                strand = "-"
            elif read.is_secondary:
                continue
            else:
                strand = "+"
            for start, end in read.get_blocks():
                start = max(0, int(start))
                end = min(int(end), chromosome_length)
                if start < end:
                    heapq.heappush(heaps[strand], (start, end))
        for strand in ("+", "-"):
            while heaps[strand]:
                consume(strand, *heapq.heappop(heaps[strand]))
            if current[strand] is not None:
                merged[strand].append(tuple(current[strand]))
    return chromosome_length, merged


def assign_candidates_indexed(
    candidates: pd.DataFrame,
    intervals: list[tuple[int, int]],
    strand: str,
) -> list[CoverageComponent]:
    """Preserve maximum-overlap assignment using a disjoint-interval index."""
    groups = [[] for _ in intervals]
    starts = [start for start, _ in intervals]
    ends = [end for _, end in intervals]
    current = candidates[candidates["strand"].eq(strand)]
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
    return [
        CoverageComponent(start, end, strand, tuple(groups[index]))
        for index, (start, end) in enumerate(intervals)
    ]


def iter_component_endpoint_batches(
    bam: pysam.AlignmentFile,
    chromosome: str,
    components_by_strand: dict[str, list[CoverageComponent]],
    candidates: pd.DataFrame,
    kernel: dict,
    max_gap: int = 1_000,
    max_span: int = 250_000,
    max_components: int = 8,
):
    """Yield bounded endpoint batches so fitted components can be released promptly."""
    for strand in ("+", "-"):
        selected = [
            (index, component)
            for index, component in enumerate(components_by_strand[strand])
            if component.candidate_indices
        ]
        chunks = []
        current = []
        chunk_start = chunk_end = None
        for item in selected:
            component = item[1]
            if current and (
                component.start - chunk_end > max_gap
                or component.end - chunk_start > max_span
                or len(current) >= max_components
            ):
                chunks.append(current)
                current = []
            if not current:
                chunk_start = component.start
                chunk_end = component.end
            else:
                chunk_end = max(chunk_end, component.end)
            current.append(item)
        if current:
            chunks.append(current)

        for chunk in chunks:
            starts = [component.start for _, component in chunk]
            values = {index: [] for index, _ in chunk}
            windows = {}
            for component_index, component in chunk:
                p0_values = [
                    oriented_coordinate(
                        candidates.loc[index, "corrected_coordinate"], strand
                    )
                    for index in component.candidate_indices
                ]
                if len(p0_values) == 1:
                    lower, upper = (
                        int(kernel["window_lower"]),
                        int(kernel["window_upper"]),
                    )
                else:
                    lower, upper = -220, 140
                windows[component_index] = (p0_values, lower, upper)
            fetch_start = chunk[0][1].start
            fetch_end = max(component.end for _, component in chunk)
            for read in bam.fetch(chromosome, fetch_start, fetch_end):
                read_strand = "-" if read.is_reverse else "+"
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read_strand != strand
                    or not read.cigarstring
                    or re.search(r"\d+N", read.cigarstring)
                ):
                    continue
                endpoint = read.reference_end if strand == "+" else read.reference_start
                local_index = bisect.bisect_right(starts, endpoint) - 1
                for neighbor in (local_index - 1, local_index, local_index + 1):
                    if not 0 <= neighbor < len(chunk):
                        continue
                    component_index, component = chunk[neighbor]
                    if (
                        component.start <= endpoint <= component.end
                        and read.reference_start < component.end
                        and read.reference_end > component.start
                    ):
                        oriented_endpoint = oriented_coordinate(endpoint, strand)
                        p0_values, lower, upper = windows[component_index]
                        if any(
                            lower <= oriented_endpoint - p0 <= upper
                            for p0 in p0_values
                        ):
                            values[component_index].append(oriented_endpoint)
            yield strand, {
                component_index: np.asarray(endpoints, dtype=int)
                for component_index, endpoints in values.items()
            }


def probability_banks(
    endpoints: np.ndarray,
    p0_values: np.ndarray,
    taus: np.ndarray,
    kernel: dict,
) -> np.ndarray:
    banks = np.empty((len(p0_values), len(taus), len(endpoints)), dtype=float)
    minimum_offset = int(endpoints.min() - p0_values.max() - taus.max())
    maximum_offset = int(endpoints.max() - p0_values.min() - taus.min())
    offset_grid = np.arange(minimum_offset, maximum_offset + 1, dtype=int)
    mass_lookup = np.clip(
        bin_masses(
            offset_grid,
            float(kernel["mu"]),
            float(kernel["sigma1"]),
            float(kernel["sigma2"]),
            -float(kernel["read_length"]),
            0.0,
            quadrature_order=12,
        ),
        1e-300,
        None,
    )
    for component, p0 in enumerate(p0_values):
        for tau_index, tau in enumerate(taus):
            offsets = endpoints - int(p0) - int(tau)
            banks[component, tau_index] = mass_lookup[
                offsets - minimum_offset
            ]
    return banks


def fit_mixture(
    endpoints: np.ndarray,
    p0_values: np.ndarray,
    kernel: dict,
    search_radius: int = 60,
    max_iterations: int = 80,
    compress_endpoints: bool = True,
) -> dict:
    if compress_endpoints:
        endpoints, endpoint_counts = np.unique(
            endpoints, return_counts=True
        )
        endpoint_counts = endpoint_counts.astype(float)
    else:
        endpoint_counts = np.ones(len(endpoints), dtype=float)
    taus = np.arange(-search_radius, search_radius + 1, dtype=int)
    banks = probability_banks(endpoints, p0_values, taus, kernel)
    log_banks = np.log(banks)
    fixed_offset = float(kernel.get("fixed_offset", -3.0))
    initial_raw_tau = int(
        round(fixed_offset + float(kernel["endpoint_location_calibration"]))
    )
    initial_index = int(np.argmin(np.abs(taus - initial_raw_tau)))
    starts = [np.full(len(p0_values), initial_index, dtype=int)]
    nearest = np.argmin(np.abs(endpoints[:, None] - p0_values[None, :]), axis=1)
    independent = []
    for component in range(len(p0_values)):
        assigned = nearest == component
        if endpoint_counts[assigned].sum() < 5:
            independent.append(initial_index)
        else:
            independent.append(
                int(
                    np.argmax(
                        log_banks[component][:, assigned]
                        @ endpoint_counts[assigned]
                    )
                )
            )
    starts.append(np.asarray(independent, dtype=int))
    fits = []
    for start_indices in starts:
        tau_indices = start_indices.copy()
        weights = np.full(len(p0_values), 1.0 / len(p0_values))
        for iteration in range(max_iterations):
            log_joint = np.vstack(
                [
                    np.log(max(weights[component], 1e-12))
                    + log_banks[component, tau_indices[component]]
                    for component in range(len(p0_values))
                ]
            )
            responsibilities = np.exp(
                log_joint - logsumexp(log_joint, axis=0, keepdims=True)
            )
            weighted_responsibilities = (
                responsibilities * endpoint_counts[None, :]
            )
            new_weights = (weighted_responsibilities.sum(axis=1) + 1.0) / (
                endpoint_counts.sum() + len(p0_values)
            )
            new_indices = np.asarray(
                [
                    int(
                        np.argmax(
                            log_banks[component]
                            @ weighted_responsibilities[component]
                        )
                    )
                    for component in range(len(p0_values))
                ]
            )
            if np.array_equal(new_indices, tau_indices) and np.max(
                np.abs(new_weights - weights)
            ) < 1e-7:
                tau_indices = new_indices
                weights = new_weights
                break
            tau_indices = new_indices
            weights = new_weights
        log_joint = np.vstack(
            [
                np.log(max(weights[component], 1e-12))
                + log_banks[component, tau_indices[component]]
                for component in range(len(p0_values))
            ]
        )
        fits.append(
            {
                "raw_taus": taus[tau_indices],
                "weights": weights,
                "log_likelihood": float(
                    (logsumexp(log_joint, axis=0) * endpoint_counts).sum()
                ),
                "iterations": iteration + 1,
            }
        )
    return max(fits, key=lambda fit: fit["log_likelihood"])


def fit_single(
    endpoints: np.ndarray,
    oriented_p0: int,
    kernel_values: dict,
    search_radius: int,
    probability_cache: dict[int, np.ndarray] | None = None,
) -> tuple[float, int]:
    kernel = Kernel(
        read_length=int(kernel_values["read_length"]),
        mu=float(kernel_values["mu"]),
        sigma1=float(kernel_values["sigma1"]),
        sigma2=float(kernel_values["sigma2"]),
        window_lower=int(kernel_values["window_lower"]),
        window_upper=int(kernel_values["window_upper"]),
    )
    coordinates = np.arange(kernel.window_lower, kernel.window_upper + 1)
    offsets = endpoints - int(oriented_p0)
    kept = offsets[
        (offsets >= kernel.window_lower) & (offsets <= kernel.window_upper)
    ]
    if not len(kept):
        return float("nan"), 0
    counts = np.bincount(
        kept - kernel.window_lower, minlength=len(coordinates)
    )
    probabilities = probability_cache or {
        tau: endpoint_probabilities(coordinates, tau, kernel, bin_masses)
        for tau in range(-search_radius, search_radius + 1)
    }
    raw_tau, _, _ = estimate_tau(
        counts, probabilities.__getitem__, search_radius
    )
    calibrated = raw_tau - float(kernel_values["endpoint_location_calibration"])
    return float(calibrated), len(kept)


def single_probability_cache(
    kernel_values: dict, search_radius: int
) -> dict[int, np.ndarray]:
    kernel = Kernel(
        read_length=int(kernel_values["read_length"]),
        mu=float(kernel_values["mu"]),
        sigma1=float(kernel_values["sigma1"]),
        sigma2=float(kernel_values["sigma2"]),
        window_lower=int(kernel_values["window_lower"]),
        window_upper=int(kernel_values["window_upper"]),
    )
    coordinates = np.arange(kernel.window_lower, kernel.window_upper + 1)
    return {
        tau: endpoint_probabilities(coordinates, tau, kernel, bin_masses)
        for tau in range(-search_radius, search_radius + 1)
    }


def optimal_errors(truth: np.ndarray, calls: np.ndarray) -> np.ndarray:
    cost = np.abs(truth[:, None] - calls[None, :])
    rows, columns = linear_sum_assignment(cost)
    return cost[rows, columns]


def fit_component_call_rows(
    component_index,
    component,
    endpoints,
    strand,
    candidates,
    kernel,
    frozen,
    fixed_offset,
    args,
    cached_single_probabilities,
):
    candidate_indices = list(component.candidate_indices)
    current = candidates.loc[candidate_indices].copy()
    current["oriented_p0"] = [
        oriented_coordinate(value, strand)
        for value in current["corrected_coordinate"]
    ]
    current = current.sort_values("oriented_p0")
    p0_values = current["oriented_p0"].to_numpy(dtype=int)
    if len(current) == 1:
        calibrated, used_reads = fit_single(
            endpoints,
            int(p0_values[0]),
            kernel,
            args.search_radius,
            cached_single_probabilities,
        )
        if used_reads < args.min_reads:
            return []
        endpoints = endpoints[
            (endpoints - int(p0_values[0]) >= kernel["window_lower"])
            & (endpoints - int(p0_values[0]) <= kernel["window_upper"])
        ]
        calibrated_mle = np.asarray([calibrated])
        fitted = {
            "raw_taus": np.asarray(
                [
                    int(
                        round(
                            calibrated
                            + float(kernel["endpoint_location_calibration"])
                        )
                    )
                ]
            ),
            "weights": np.asarray([1.0]),
        }
        component_model = "single"
    else:
        retained = np.any(
            (endpoints[:, None] - p0_values[None, :] >= -220)
            & (endpoints[:, None] - p0_values[None, :] <= 140),
            axis=1,
        )
        endpoints = endpoints[retained]
        if len(endpoints) < args.min_reads * len(current):
            return []
        fitted = fit_mixture(
            endpoints,
            p0_values,
            kernel,
            search_radius=args.search_radius,
        )
        calibrated_mle = (
            fitted["raw_taus"].astype(float)
            - float(kernel["endpoint_location_calibration"])
        )
        component_model = "mixture"
    if component_model == "mixture":
        mixture_refined = use_mixture_mle(
            p0_values, args.mixture_mle_max_median_gap
        )
        gated_tau = (
            calibrated_mle
            if mixture_refined
            else np.full(len(current), fixed_offset)
        )
    else:
        mixture_refined = False
        gated_tau = frozen_calibration(calibrated_mle, frozen)
    rows = []
    for offset, (candidate_index, candidate) in enumerate(current.iterrows()):
        p0 = int(candidate["corrected_coordinate"])
        rows.append(
            {
                "candidate_index": int(candidate_index),
                "component_index": component_index,
                "component_candidates": len(current),
                "component_endpoint_reads": len(endpoints),
                "component_model": component_model,
                "mixture_refined": mixture_refined,
                "strand": strand,
                "p0": p0,
                "mixture_weight": float(fitted["weights"][offset]),
                "raw_tau": int(fitted["raw_taus"][offset]),
                "mle_tau": float(calibrated_mle[offset]),
                "gated_tau": float(gated_tau[offset]),
                "fixed_coordinate": predicted_genomic(
                    p0, fixed_offset, strand
                ),
                "mle_coordinate": predicted_genomic(
                    p0, calibrated_mle[offset], strand
                ),
                "model_coordinate": predicted_genomic(
                    p0, gated_tau[offset], strand
                ),
            }
        )
    return rows


def inference(args: argparse.Namespace) -> pd.DataFrame:
    kernel = json.loads(args.kernel.read_text())
    frozen = json.loads(args.frozen.read_text())
    fixed_offset = model_fixed_offset(kernel, frozen)
    kernel["fixed_offset"] = fixed_offset
    cached_single_probabilities = single_probability_cache(
        kernel, args.search_radius
    )
    raw_candidates = pd.read_csv(args.candidates, sep="\t")
    candidates = raw_candidates[
        ["chrom", "start", "end", "strand", "corrected_coordinate"]
    ].copy()
    chromosome_length, intervals = stream_positive_intervals(
        args.bam, args.chromosome
    )
    components_by_strand = {
        strand: assign_candidates_indexed(candidates, intervals[strand], strand)
        for strand in ("+", "-")
    }
    rows = []
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for strand, endpoint_batch in iter_component_endpoint_batches(
            bam,
            args.chromosome,
            components_by_strand,
            candidates,
            kernel,
        ):
            for component_index, endpoints in endpoint_batch.items():
                rows.extend(
                    fit_component_call_rows(
                        component_index,
                        components_by_strand[strand][component_index],
                        endpoints,
                        strand,
                        candidates,
                        kernel,
                        frozen,
                        fixed_offset,
                        args,
                        cached_single_probabilities,
                    )
                )
    calls = pd.DataFrame(rows)
    if calls.empty:
        return calls
    calls["te_id"] = raw_candidates.loc[calls["candidate_index"], "te_id"].to_numpy()
    return calls


def evaluate(calls: pd.DataFrame, truth: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for te_id, te_truth in truth.groupby("te_id"):
        order = len(te_truth)
        if order < 2:
            continue
        te_calls = calls[calls["te_id"].eq(te_id)]
        base = {
            "te_id": te_id,
            "truth_order": order,
            "candidate_count": len(te_calls),
            "candidate_sufficient": len(te_calls) >= order,
            "exact_candidate_count": len(te_calls) == order,
        }
        if len(te_calls) < order:
            rows.append(base)
            continue
        truth_coordinates = te_truth["truth_coordinate"].to_numpy(dtype=float)
        p0_errors = optimal_errors(
            truth_coordinates, te_calls["p0"].to_numpy(dtype=float)
        )
        base["candidate_grid_supported"] = bool(
            len(p0_errors) == order and np.all(p0_errors <= 66)
        )
        base["candidate_grid_max_error"] = float(p0_errors.max())
        for label, column in [
            ("fixed", "fixed_coordinate"),
            ("mle", "mle_coordinate"),
            ("model", "model_coordinate"),
        ]:
            errors = optimal_errors(
                truth_coordinates, te_calls[column].to_numpy(dtype=float)
            )
            base[f"{label}_mean_error"] = float(errors.mean())
            base[f"{label}_max_error"] = float(errors.max())
            base[f"{label}_complete_10bp"] = bool(np.all(errors <= 10))
            base[f"{label}_errors"] = ",".join(f"{value:.6g}" for value in errors)
        rows.append(base)
    detail = pd.DataFrame(rows)
    metric_rows = []
    for order, group in [("all", detail), *detail.groupby("truth_order", sort=True)]:
        eligible = group[group["candidate_grid_supported"].eq(True)].copy()
        for method in ("fixed", "mle", "model"):
            errors = np.asarray(
                [
                    float(value)
                    for text in eligible[f"{method}_errors"].dropna()
                    for value in text.split(",")
                ]
            )
            metric_rows.append(
                {
                    "truth_order": order,
                    "method": method,
                    "total_te": len(group),
                    "eligible_te": len(eligible),
                    "candidate_sufficient_rate": float(
                        group["candidate_sufficient"].mean()
                    ),
                    "candidate_grid_support_rate": float(
                        group["candidate_grid_supported"].eq(True).mean()
                    ),
                    "site_mean_error": float(errors.mean()) if len(errors) else np.nan,
                    "site_median_error": float(np.median(errors)) if len(errors) else np.nan,
                    "site_p90_error": float(np.quantile(errors, 0.9)) if len(errors) else np.nan,
                    "complete_10bp_eligible": float(
                        eligible[f"{method}_complete_10bp"].mean()
                    )
                    if len(eligible)
                    else np.nan,
                    "complete_10bp_all": float(
                        group[f"{method}_complete_10bp"].eq(True).mean()
                    ),
                }
            )
    return detail, pd.DataFrame(metric_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--search-radius", type=int, default=60)
    parser.add_argument("--min-reads", type=int, default=20)
    parser.add_argument("--mixture-mle-max-median-gap", type=float, default=300.0)
    args = parser.parse_args()
    calls = inference(args)
    truth = pd.read_csv(args.truth, sep="\t")
    detail, metrics = evaluate(calls, truth)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    calls.to_csv(args.output_dir / "component_calls.tsv", sep="\t", index=False)
    detail.to_csv(args.output_dir / "te_evaluation.tsv", sep="\t", index=False)
    metrics.to_csv(args.output_dir / "metrics.tsv", sep="\t", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
