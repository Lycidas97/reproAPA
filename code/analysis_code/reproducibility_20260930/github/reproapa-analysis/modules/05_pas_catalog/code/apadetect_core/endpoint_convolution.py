"""Frozen truncated-convolution kernels for PAS localization."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit


@dataclass(frozen=True)
class Kernel:
    read_length: int
    mu: float
    sigma1: float
    sigma2: float
    window_lower: int
    window_upper: int


def equal_peak_histogram(frame, lower: int, upper: int, min_reads: int = 50):
    """Average conditional endpoint histograms so large peaks do not dominate."""
    coordinates = np.arange(lower, upper + 1)
    profiles = []
    for _, group in frame.groupby("pas", sort=False):
        values = group["distance_to_pas"].to_numpy(dtype=int)
        values = values[(values >= lower) & (values <= upper)]
        if len(values) < min_reads:
            continue
        counts = np.bincount(values - lower, minlength=len(coordinates))
        profiles.append(counts / counts.sum())
    if len(profiles) < 20:
        raise ValueError("at least 20 training peaks are required")
    return coordinates, np.mean(profiles, axis=0), len(profiles)


def fit_kernel(coordinates, target, read_length, bin_masses):
    """Fit one sample-level TN[-L,0] + Normal endpoint kernel."""
    coordinates = np.asarray(coordinates, dtype=int)
    target = np.asarray(target, dtype=float)
    eta_bounds = (-10.0, 10.0)
    log_sd_bounds = (np.log(0.25), np.log(2.0 * read_length))

    def decode(parameters):
        return (
            -read_length + read_length * expit(parameters[0]),
            float(np.exp(parameters[1])),
            float(np.exp(parameters[2])),
        )

    def probabilities(parameters):
        mu, sigma1, sigma2 = decode(parameters)
        masses = bin_masses(
            coordinates,
            mu,
            sigma1,
            sigma2,
            -float(read_length),
            0.0,
            quadrature_order=12,
        )
        return masses / masses.sum()

    def objective(parameters):
        probs = np.clip(probabilities(parameters), 1e-300, 1.0)
        return -float(np.sum(target * np.log(probs)))

    starts = (
        (-0.20 * read_length, 0.20 * read_length, 4.0),
        (-0.50 * read_length, 0.35 * read_length, 10.0),
        (-0.80 * read_length, 0.20 * read_length, 20.0),
    )
    results = []
    for mu, sigma1, sigma2 in starts:
        initial = np.array(
            [
                logit((mu + read_length) / read_length),
                np.log(sigma1),
                np.log(sigma2),
            ]
        )
        results.append(
            minimize(
                objective,
                initial,
                method="L-BFGS-B",
                bounds=[eta_bounds, log_sd_bounds, log_sd_bounds],
                options={"maxiter": 300, "ftol": 1e-12, "gtol": 1e-8},
            )
        )
    best = min(results, key=lambda result: float(result.fun))
    mu, sigma1, sigma2 = decode(best.x)
    return {
        "mu": mu,
        "sigma1": sigma1,
        "sigma2": sigma2,
        "loss": float(best.fun),
        "converged": bool(best.success),
        "start_loss_range": float(
            max(result.fun for result in results)
            - min(result.fun for result in results)
        ),
        "probabilities": probabilities(best.x),
    }


def endpoint_probabilities(coordinates, tau, kernel, bin_masses):
    """Endpoint probabilities conditional on the fixed candidate-relative window."""
    masses = bin_masses(
        np.asarray(coordinates, dtype=int) - int(tau),
        kernel.mu,
        kernel.sigma1,
        kernel.sigma2,
        -float(kernel.read_length),
        0.0,
        quadrature_order=12,
    )
    total = float(masses.sum())
    if not np.isfinite(total) or total <= 0:
        raise FloatingPointError("endpoint window has no finite probability mass")
    return masses / total


def coverage_probabilities(
    coordinates,
    tau,
    kernel,
    length_values,
    length_probabilities,
    bin_masses,
):
    """Expected coverage from endpoint mass integrated over aligned read lengths."""
    coordinates = np.asarray(coordinates, dtype=int)
    lengths = np.asarray(length_values, dtype=int)
    weights = np.asarray(length_probabilities, dtype=float)
    if len(lengths) != len(weights) or np.any(lengths <= 0):
        raise ValueError("invalid aligned-length distribution")
    lower = int(coordinates.min() + 1 - tau)
    upper = int(coordinates.max() + lengths.max() - tau)
    endpoint_grid = np.arange(lower, upper + 1)
    masses = bin_masses(
        endpoint_grid,
        kernel.mu,
        kernel.sigma1,
        kernel.sigma2,
        -float(kernel.read_length),
        0.0,
        quadrature_order=12,
    )
    cumulative = np.r_[0.0, np.cumsum(masses)]
    expected = np.zeros(len(coordinates), dtype=float)
    grid_lower = int(endpoint_grid[0])
    for length, weight in zip(lengths, weights):
        first = coordinates + 1 - int(tau) - grid_lower
        last = coordinates + int(length) - int(tau) - grid_lower + 1
        expected += weight * (cumulative[last] - cumulative[first])
    total = float(expected.sum())
    if not np.isfinite(total) or total <= 0:
        raise FloatingPointError("coverage window has no finite probability mass")
    return expected / total


def estimate_tau(counts, probability_function, search_radius=60):
    """Return the integer shift maximizing a frozen-kernel likelihood."""
    observed = np.asarray(counts, dtype=float)
    if observed.sum() <= 0:
        raise ValueError("counts must have positive mass")
    scores = []
    taus = np.arange(-int(search_radius), int(search_radius) + 1)
    for tau in taus:
        probabilities = np.clip(probability_function(int(tau)), 1e-300, 1.0)
        scores.append(float(np.sum(observed * np.log(probabilities))))
    order = np.argsort(scores)[::-1]
    gap = float(scores[order[0]] - scores[order[1]])
    return int(taus[order[0]]), gap, bool(abs(int(taus[order[0]])) == search_radius)
