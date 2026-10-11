"""Resolve Monte Carlo uncertainty using the exact joint hypergeometric law."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def adjusted(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    result = np.empty(len(values))
    result[order] = np.minimum(1, np.minimum.accumulate(
        (values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1])
    return result


def verify(root, out):
    out.mkdir(parents=True, exist_ok=True)
    read = lambda name: pd.read_csv(root / 'original-nulls' / name, sep='\t', float_precision='round_trip')
    results = read('overlap_threshold_layer_sensitivity.tsv')
    diagnostics = read('stratum_diagnostics.tsv.gz')
    rows, archive = [], {}
    for row in results.itertuples(index=False):
        cells = diagnostics[diagnostics.contrast.eq(row.contrast) & diagnostics.layer.eq(row.layer) &
                            diagnostics.threshold.eq(row.threshold) & diagnostics.side.eq('joint')]
        probabilities = np.array([1.])
        offset = 0
        for cell in cells.itertuples(index=False):
            n, nx, ny = int(cell.n), int(cell.human_events), int(cell.mouse_events)
            lower, upper = max(0, nx + ny - n), min(nx, ny)
            offset += lower
            if upper == lower:
                continue
            distribution = hypergeom.pmf(np.arange(lower, upper + 1), n, nx, ny)
            assert abs(distribution.sum() - 1) < 1e-12
            probabilities = np.convolve(probabilities, distribution)
        assert abs(probabilities.sum() - 1) < 1e-11
        support = offset + np.arange(len(probabilities))
        mean = float(probabilities @ support)
        variance = float(probabilities @ (support - mean) ** 2)
        assert abs(mean - row.analytic_null_mean) < 1e-10
        assert abs(np.sqrt(variance) - row.analytic_null_sd) < 1e-10
        p = min(1., float(probabilities[support >= row.observed].sum()))
        rows.append(dict(contrast=row.contrast, layer=row.layer, threshold=row.threshold,
                         observed=row.observed, analytic_expected=mean, exact_p=p,
                         monte_carlo_p=row.empirical_p, monte_carlo_q16=row.q_16,
                         monte_carlo_permutations=row.permutations, minimum_overlap=offset))
        archive[f'{row.contrast}|{row.layer}|{row.threshold}|probability'] = probabilities
        archive[f'{row.contrast}|{row.layer}|{row.threshold}|support'] = support
    frame = pd.DataFrame(rows)
    frame['exact_q16'] = adjusted(frame.exact_p)
    frame['exact_q2_primary'] = np.nan
    primary = frame.layer.eq('supported50') & frame.threshold.eq(.02)
    frame.loc[primary, 'exact_q2_primary'] = adjusted(frame.loc[primary, 'exact_p'])
    frame.to_csv(out / 'joint-exact-results.tsv', sep='\t', index=False, float_format='%.17g')
    np.savez_compressed(out / 'joint-exact-distributions.npz', **archive)
    result = dict(status='PASS', settings=len(frame), tests_skipped=0,
                  exact_distribution_normalization_mean_and_sd_verified=True,
                  corrected_PM_significant_settings_exact_q16=int((frame.contrast.eq('PM') & frame.exact_q16.lt(.05)).sum()))
    (out / 'exact-verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(frame[frame.contrast.eq('PM')].to_string(index=False))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    verify(args.root, args.out)
