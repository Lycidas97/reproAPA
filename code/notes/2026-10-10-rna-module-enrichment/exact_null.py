"""Fixed-stratum over-representation and conditional detectable effects."""
import math
from functools import lru_cache
import numpy as np


@lru_cache(maxsize=200000)
def component(N, K, n, odds=1.0):
    assert 0 <= K <= N and 0 <= n <= N and odds > 0
    lo, hi = max(0, n + K - N), min(n, K)
    x = np.arange(lo, hi + 1)
    base = np.array([math.comb(K, int(k)) * math.comb(N - K, n - int(k))
                     / math.comb(N, n) for k in x])
    if odds == 1:
        return lo, base / base.sum()
    logs = np.log(base) + x * np.log(odds)
    weights = np.exp(logs - logs.max())
    return lo, weights / weights.sum()


def distribution(params, odds=1.0):
    offset, values = 0, np.array([1.0])
    for N, K, n in params:
        if not n or not K:
            continue
        if K == N:
            offset += n
            continue
        lo, part = component(int(N), int(K), int(n), odds)
        offset += lo
        values = np.convolve(values, part)
    assert abs(values.sum() - 1) < 1e-12
    return offset, values / values.sum()


def upper_tail(dist, hits):
    lo, p = dist
    if hits <= lo:
        return 1.0
    return min(1.0, float(p[hits - lo:].sum()))


def rejection_probability(test, generated, alpha):
    lo, p = generated
    return float(sum(weight for j, weight in enumerate(p)
                     if upper_tail(test, lo + j) <= alpha))


def detectable_effect(params, alpha, target_power=0.8, maximum_odds=1000000):
    null = distribution(params)
    lo, p = null
    critical = next((lo + j for j in range(len(p))
                     if upper_tail(null, lo + j) <= alpha), None)
    if critical is None:
        return dict(critical_hits=None, minimum_odds=None, minimum_expected_hits=None,
                    minimum_expected_fold=None, attainable=False, unattainable_reason='no_rejection_support')
    power = lambda odds: rejection_probability(null, distribution(params, odds), alpha)
    if power(maximum_odds) < target_power:
        return dict(critical_hits=critical, minimum_odds=None, minimum_expected_hits=None,
                    minimum_expected_fold=None, attainable=False, unattainable_reason='above_maximum_searched_odds')
    low, high = 1.0, maximum_odds
    for _ in range(48):
        mid = np.sqrt(low * high)
        if power(mid) < target_power:
            low = mid
        else:
            high = mid
    altlo, alt = distribution(params, high)
    mean = float(np.arange(altlo, altlo + len(alt)) @ alt)
    nullmean = float(np.arange(lo, lo + len(p)) @ p)
    return dict(critical_hits=critical, minimum_odds=float(high), minimum_expected_hits=mean,
                minimum_expected_fold=mean / nullmean if nullmean else None, attainable=True, unattainable_reason='')


def bh(values):
    p = np.asarray(values, float)
    order = np.argsort(p, kind='stable')
    out = np.empty(len(p))
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    out[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return out
