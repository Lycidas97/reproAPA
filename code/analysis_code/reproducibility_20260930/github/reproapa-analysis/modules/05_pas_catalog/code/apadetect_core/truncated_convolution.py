"""Numerically stable truncated-normal/Gaussian convolution utilities.

The functions in this module are a clean-room reconstruction from the
DOC-176 mathematical contract.  ``original_center_loglik`` is explicitly
named compatibility behavior for the located apabenchmark implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, log_ndtr, logit


LOG_2PI = float(np.log(2.0 * np.pi))


def _logdiffexp(log_x: np.ndarray, log_y: np.ndarray) -> np.ndarray:
    """Return log(exp(log_x) - exp(log_y)) for log_x >= log_y."""

    ratio = np.exp(np.minimum(log_y - log_x, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        return log_x + np.log1p(-ratio)


def log_cdf_difference(upper: np.ndarray, lower: np.ndarray) -> np.ndarray:
    """Stable log(Phi(upper) - Phi(lower))."""

    upper_arr, lower_arr = np.broadcast_arrays(
        np.asarray(upper, dtype=float), np.asarray(lower, dtype=float)
    )
    if np.any(upper_arr <= lower_arr):
        raise ValueError("upper must be greater than lower")
    use_survival = lower_arr > 0.0
    cdf_result = _logdiffexp(log_ndtr(upper_arr), log_ndtr(lower_arr))
    survival_result = _logdiffexp(
        log_ndtr(-lower_arr), log_ndtr(-upper_arr)
    )
    return np.where(use_survival, survival_result, cdf_result)


def log_density(
    z: np.ndarray | float,
    mu: float,
    sigma1: float,
    sigma2: float,
    a: float,
    b: float,
) -> np.ndarray:
    """Log density of TN[a,b](mu,sigma1^2) + N(0,sigma2^2)."""

    if sigma1 <= 0.0 or sigma2 <= 0.0 or not a < b:
        raise ValueError("require sigma1 > 0, sigma2 > 0, and a < b")
    values = np.asarray(z, dtype=float)
    total_var = sigma1 * sigma1 + sigma2 * sigma2
    total_sd = np.sqrt(total_var)
    conditional_mean = (
        sigma2 * sigma2 * mu + sigma1 * sigma1 * values
    ) / total_var
    conditional_sd = sigma1 * sigma2 / total_sd
    gaussian = (
        -0.5 * ((values - mu) / total_sd) ** 2
        - np.log(total_sd)
        - 0.5 * LOG_2PI
    )
    conditional = log_cdf_difference(
        (b - conditional_mean) / conditional_sd,
        (a - conditional_mean) / conditional_sd,
    )
    truncation = log_cdf_difference(
        np.asarray((b - mu) / sigma1),
        np.asarray((a - mu) / sigma1),
    )
    return gaussian + conditional - truncation


def density(
    z: np.ndarray | float,
    mu: float,
    sigma1: float,
    sigma2: float,
    a: float,
    b: float,
) -> np.ndarray:
    return np.exp(log_density(z, mu, sigma1, sigma2, a, b))


@lru_cache(maxsize=None)
def _legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    return np.polynomial.legendre.leggauss(order)


def bin_masses(
    coordinates: Sequence[int] | np.ndarray,
    mu: float,
    sigma1: float,
    sigma2: float,
    a: float,
    b: float,
    *,
    quadrature_order: int = 16,
) -> np.ndarray:
    """Integrate the density over one-base bins centered at coordinates."""

    centers = np.asarray(coordinates, dtype=float)
    nodes, weights = _legendre_rule(quadrature_order)
    evaluation_points = centers[:, None] + 0.5 * nodes[None, :]
    values = density(evaluation_points, mu, sigma1, sigma2, a, b)
    return 0.5 * np.sum(values * weights[None, :], axis=1)


def window_probabilities(
    coordinates: Sequence[int] | np.ndarray,
    mu: float,
    sigma1: float,
    sigma2: float,
    a: float,
    b: float,
    *,
    quadrature_order: int = 16,
) -> np.ndarray:
    masses = bin_masses(
        coordinates,
        mu,
        sigma1,
        sigma2,
        a,
        b,
        quadrature_order=quadrature_order,
    )
    total = float(masses.sum())
    if not np.isfinite(total) or total <= 0.0:
        raise FloatingPointError("window mass is not finite and positive")
    return masses / total


def center_density_probabilities(
    coordinates: Sequence[int] | np.ndarray,
    mu: float,
    sigma1: float,
    sigma2: float,
    a: float,
    b: float,
) -> np.ndarray:
    values = density(np.asarray(coordinates), mu, sigma1, sigma2, a, b)
    return values / values.sum()


def multinomial_loglik(counts: np.ndarray, probabilities: np.ndarray) -> float:
    counts_arr = np.asarray(counts, dtype=float)
    probs_arr = np.asarray(probabilities, dtype=float)
    if counts_arr.shape != probs_arr.shape:
        raise ValueError("counts and probabilities must have the same shape")
    return float(np.sum(counts_arr * np.log(np.clip(probs_arr, 1e-300, 1.0))))


def original_center_loglik(
    counts: np.ndarray,
    coordinates: np.ndarray,
    mu: float,
    sigma1: float,
    sigma2: float,
    a: float,
    b: float,
) -> float:
    """Located-code compatibility mode: unconditioned center densities."""

    values = np.clip(
        density(coordinates, mu, sigma1, sigma2, a, b), 1e-20, None
    )
    positive = np.asarray(counts) > 0
    return float(np.sum(np.asarray(counts)[positive] * np.log(values[positive])))


def anchored_probabilities(
    coordinates: np.ndarray,
    tau: int,
    read_length: int,
    mu: float,
    sigma1: float,
    sigma2: float,
    *,
    quadrature_order: int = 12,
) -> np.ndarray:
    if not -read_length < mu < 0:
        raise ValueError("anchored mu must be strictly inside [-L,0]")
    return window_probabilities(
        np.asarray(coordinates) - int(tau),
        mu,
        sigma1,
        sigma2,
        -float(read_length),
        0.0,
        quadrature_order=quadrature_order,
    )


@dataclass(frozen=True)
class AnchoredFit:
    tau: int
    mu: float
    sigma1: float
    sigma2: float
    log_likelihood: float
    second_best_gap: float
    boundary_hit: bool
    converged: bool
    across_start_range: float


def fit_anchored_profile(
    counts: np.ndarray,
    coordinates: np.ndarray,
    candidate_taus: Iterable[int],
    read_length: int,
    *,
    starts: Sequence[tuple[float, float, float]] | None = None,
    quadrature_order: int = 8,
    maxiter: int = 180,
) -> AnchoredFit:
    """Fit shape parameters at every prespecified integer tau candidate."""

    counts_arr = np.asarray(counts, dtype=float)
    coordinates_arr = np.asarray(coordinates, dtype=int)
    if starts is None:
        starts = (
            (-0.25 * read_length, 0.20 * read_length, 4.0),
            (-0.50 * read_length, 0.35 * read_length, 10.0),
            (-0.80 * read_length, 0.15 * read_length, 20.0),
        )
    eta_bounds = (-12.0, 12.0)
    log_sd_bounds = (np.log(0.25), np.log(2.0 * read_length))
    all_candidates: list[dict[str, object]] = []

    for tau in candidate_taus:
        start_results = []

        def objective(parameters: np.ndarray) -> float:
            eta_mu, log_sigma1, log_sigma2 = parameters
            mu_value = -read_length + read_length * expit(eta_mu)
            try:
                probabilities = anchored_probabilities(
                    coordinates_arr,
                    int(tau),
                    read_length,
                    float(mu_value),
                    float(np.exp(log_sigma1)),
                    float(np.exp(log_sigma2)),
                    quadrature_order=quadrature_order,
                )
            except (ValueError, FloatingPointError):
                return 1e300
            return -multinomial_loglik(counts_arr, probabilities)

        for start_mu, start_sigma1, start_sigma2 in starts:
            initial = np.array(
                [
                    logit(np.clip((start_mu + read_length) / read_length, 1e-6, 1.0 - 1e-6)),
                    np.log(np.clip(start_sigma1, 0.25, 2.0 * read_length)),
                    np.log(np.clip(start_sigma2, 0.25, 2.0 * read_length)),
                ],
                dtype=float,
            )
            result = minimize(
                objective,
                initial,
                method="L-BFGS-B",
                bounds=[eta_bounds, log_sd_bounds, log_sd_bounds],
                options={"maxiter": maxiter, "ftol": 1e-11, "gtol": 1e-7},
            )
            start_results.append(result)
        best = min(start_results, key=lambda result: float(result.fun))
        logliks = np.array([-float(result.fun) for result in start_results])
        params = np.asarray(best.x, dtype=float)
        tolerance = 1e-3
        fitted_mu = -read_length + read_length * expit(params[0])
        boundary_hit = bool(
            abs(params[0] - eta_bounds[0]) < tolerance
            or abs(params[0] - eta_bounds[1]) < tolerance
            or abs(params[1] - log_sd_bounds[0]) < tolerance
            or abs(params[1] - log_sd_bounds[1]) < tolerance
            or abs(params[2] - log_sd_bounds[0]) < tolerance
            or abs(params[2] - log_sd_bounds[1]) < tolerance
        )
        all_candidates.append(
            {
                "tau": int(tau),
                "mu": float(fitted_mu),
                "sigma1": float(np.exp(params[1])),
                "sigma2": float(np.exp(params[2])),
                "loglik": -float(best.fun),
                "converged": bool(best.success),
                "boundary_hit": boundary_hit,
                "start_range": float(logliks.max() - logliks.min()),
            }
        )

    ordered = sorted(all_candidates, key=lambda row: float(row["loglik"]), reverse=True)
    winner = ordered[0]
    gap = (
        float(winner["loglik"]) - float(ordered[1]["loglik"])
        if len(ordered) > 1
        else float("nan")
    )
    return AnchoredFit(
        tau=int(winner["tau"]),
        mu=float(winner["mu"]),
        sigma1=float(winner["sigma1"]),
        sigma2=float(winner["sigma2"]),
        log_likelihood=float(winner["loglik"]),
        second_best_gap=gap,
        boundary_hit=bool(winner["boundary_hit"]),
        converged=bool(winner["converged"]),
        across_start_range=float(winner["start_range"]),
    )


def located_density(
    coordinates: np.ndarray,
    tau: float,
    mu: float,
    sigma1: float,
    sigma2: float,
    a: float,
    b: float,
) -> np.ndarray:
    return density(np.asarray(coordinates) - tau, mu, sigma1, sigma2, a, b)
