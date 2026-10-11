"""Conservative multi-candidate selection from barcode/UMI heterogeneity."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import betaln, gammaln


DEFAULT_MULTI_SELECTION_MODE = "coverage"


def coverage_boundary_evidence(
    upstream: float,
    downstream: float,
    positions: np.ndarray,
    smoothed_coverage: np.ndarray,
    read_length: int,
) -> dict[str, float]:
    """Measure separation of two adjacent candidates on oriented coverage."""
    upstream_peak = (positions >= upstream - read_length) & (positions <= upstream + 10)
    downstream_peak = (
        (positions >= downstream - read_length) & (positions <= downstream + 10)
    )
    valley = (positions > upstream) & (positions < downstream)
    upstream_height = float(smoothed_coverage[upstream_peak].max(initial=0))
    downstream_height = float(smoothed_coverage[downstream_peak].max(initial=0))
    valley_height = float(smoothed_coverage[valley].min(initial=np.inf))
    if not np.isfinite(valley_height):
        valley_height = min(upstream_height, downstream_height)
    denominator = max(min(upstream_height, downstream_height), 1.0)
    return {
        "gap_bp": float(downstream - upstream),
        "gap_over_read_length": float((downstream - upstream) / read_length),
        "upstream_peak_coverage": upstream_height,
        "downstream_peak_coverage": downstream_height,
        "valley_coverage": valley_height,
        "valley_ratio": float(np.clip(valley_height / denominator, 0.0, 1.0)),
    }


def select_coverage_candidates(
    candidate_indices: list[int],
    oriented_coordinates: np.ndarray,
    positions: np.ndarray,
    smoothed_coverage: np.ndarray,
    read_length: int,
    maximum_gap_read_lengths: float,
    valley_ratio_threshold: float,
    barcode_selected: list[int] | None = None,
) -> tuple[list[int], list[dict]]:
    """Merge only close shoulders; distance, a deep valley, or barcode evidence rescues."""
    order = np.argsort(np.asarray(oriented_coordinates, dtype=float))
    indices = [candidate_indices[int(index)] for index in order]
    coordinates = np.asarray(oriented_coordinates, dtype=float)[order]
    barcode_selected = set(barcode_selected or [])
    selected = [indices[-1]]
    history = []
    for position in range(len(indices) - 2, -1, -1):
        evidence = coverage_boundary_evidence(
            coordinates[position],
            coordinates[position + 1],
            positions,
            smoothed_coverage,
            read_length,
        )
        if evidence["gap_over_read_length"] > maximum_gap_read_lengths:
            reason = "distance_retained"
        elif evidence["valley_ratio"] <= valley_ratio_threshold:
            reason = "valley_retained"
        elif indices[position] in barcode_selected:
            reason = "barcode_rescued"
        else:
            reason = "shoulder_merged"
        retained = reason != "shoulder_merged"
        if retained:
            selected.append(indices[position])
        history.append(
            {
                "tested_candidate_index": indices[position],
                "retained": retained,
                "reason": reason,
                "maximum_gap_read_lengths": maximum_gap_read_lengths,
                "valley_ratio_threshold": valley_ratio_threshold,
                **evidence,
            }
        )
    selected_set = set(selected)
    return [index for index in indices if index in selected_set], history


def beta_binomial_evidence(
    successes: np.ndarray,
    totals: np.ndarray,
) -> dict[str, float]:
    """Compare a shared usage ratio with beta-binomial heterogeneity."""
    successes = np.asarray(successes, dtype=float)
    totals = np.asarray(totals, dtype=float)
    valid = totals > 0
    successes, totals = successes[valid], totals[valid]
    informative = int(np.sum(totals >= 2))
    result = {
        "barcode_count": float(len(totals)),
        "informative_barcode_count": float(informative),
        "beta_binomial_delta_bic": 0.0,
        "beta_binomial_rho": 0.0,
    }
    if len(totals) < 2 or totals.sum() == 0:
        return result
    combinations = (
        gammaln(totals + 1)
        - gammaln(successes + 1)
        - gammaln(totals - successes + 1)
    )
    mean = float(np.clip(successes.sum() / totals.sum(), 1e-4, 1 - 1e-4))

    def objective(log_concentration: float) -> float:
        concentration = np.exp(log_concentration)
        alpha = mean * concentration
        beta = (1.0 - mean) * concentration
        return -float(
            np.sum(
                combinations
                + betaln(successes + alpha, totals - successes + beta)
                - betaln(alpha, beta)
            )
        )

    fitted = minimize_scalar(
        objective,
        method="bounded",
        bounds=(-5.0, 15.0),
        options={"xatol": 1e-5},
    )
    if not fitted.success or not np.isfinite(fitted.fun):
        raise RuntimeError(f"beta-binomial fit failed: {fitted.message}")
    concentration = float(np.exp(fitted.x))
    probability = float(np.clip(successes.sum() / totals.sum(), 1e-9, 1 - 1e-9))
    binomial_log_likelihood = float(
        np.sum(
            combinations
            + successes * np.log(probability)
            + (totals - successes) * np.log1p(-probability)
        )
    )
    observations = max(len(totals), 2)
    beta_binomial_log_likelihood = -float(fitted.fun)
    result.update(
        {
            "beta_binomial_delta_bic": float(
                (-2.0 * binomial_log_likelihood + np.log(observations))
                - (-2.0 * beta_binomial_log_likelihood + 2.0 * np.log(observations))
            ),
            "beta_binomial_rho": float(1.0 / (concentration + 1.0)),
            "beta_binomial_mean": mean,
        }
    )
    return result


def select_candidates(
    candidate_indices: list[int],
    oriented_coordinates: np.ndarray,
    molecules: list[tuple[str, int]],
    bic_threshold: float,
    minimum_informative_barcodes: int,
) -> tuple[list[int], list[dict]]:
    """Keep an upstream boundary only when barcode evidence is very strong."""
    if len(candidate_indices) < 2:
        return list(candidate_indices), []
    order = np.argsort(np.asarray(oriented_coordinates, dtype=float))
    ordered_indices = [candidate_indices[int(index)] for index in order]
    coordinates = np.asarray(oriented_coordinates, dtype=float)[order]
    boundaries = (coordinates[:-1] + coordinates[1:]) / 2.0
    by_candidate = [defaultdict(int) for _ in ordered_indices]
    for barcode, endpoint in molecules:
        assignment = int(np.searchsorted(boundaries, endpoint, side="right"))
        by_candidate[assignment][barcode] += 1

    selected = [ordered_indices[-1]]
    downstream_group = {len(ordered_indices) - 1}
    history = []
    for position in range(len(ordered_indices) - 2, -1, -1):
        barcodes = set(by_candidate[position])
        for member in downstream_group:
            barcodes.update(by_candidate[member])
        successes = []
        totals = []
        upstream_molecules = 0
        downstream_molecules = 0
        for barcode in barcodes:
            upstream = by_candidate[position].get(barcode, 0)
            downstream = sum(
                by_candidate[member].get(barcode, 0)
                for member in downstream_group
            )
            upstream_molecules += upstream
            downstream_molecules += downstream
            successes.append(upstream)
            totals.append(upstream + downstream)
        evidence = beta_binomial_evidence(
            np.asarray(successes, dtype=int),
            np.asarray(totals, dtype=int),
        )
        accepted = (
            upstream_molecules > 0
            and downstream_molecules > 0
            and evidence["informative_barcode_count"]
            >= minimum_informative_barcodes
            and evidence["beta_binomial_delta_bic"] >= bic_threshold
        )
        history.append(
            {
                "tested_candidate_index": ordered_indices[position],
                "downstream_group_candidate_indices": [
                    ordered_indices[index] for index in sorted(downstream_group)
                ],
                "upstream_molecules": upstream_molecules,
                "downstream_molecules": downstream_molecules,
                "threshold": bic_threshold,
                "accepted": accepted,
                **evidence,
            }
        )
        if accepted:
            selected.append(ordered_indices[position])
            downstream_group = {position}
        else:
            downstream_group.add(position)
    selected_set = set(selected)
    return (
        [index for index in ordered_indices if index in selected_set],
        history,
    )
