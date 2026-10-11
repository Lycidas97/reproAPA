"""Replay the joint sensitivity grid from its recorded stratum margins."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def verify(root, config, out):
    out.mkdir(parents=True, exist_ok=True)
    source = Path(json.loads(config.read_text())['analysis_root']) / 'cross_overlap'
    read = lambda path: pd.read_csv(path, sep='\t', float_precision='round_trip')
    before = read(source / 'new_shared_members.tsv.gz')
    members = read(root / 'original-nulls/new_shared_members.tsv.gz')
    pd.testing.assert_frame_equal(members[before.columns], before, check_dtype=False, rtol=0, atol=1e-12)
    results = read(root / 'original-nulls/overlap_threshold_layer_sensitivity.tsv')
    old_results = read(source / 'overlap_threshold_layer_sensitivity.tsv')
    diagnostic = read(root / 'original-nulls/stratum_diagnostics.tsv.gz')
    archive = {}
    assert len(results) == 16
    for row in results.itertuples(index=False):
        strata = diagnostic[diagnostic.contrast.eq(row.contrast) & diagnostic.layer.eq(row.layer) &
                            diagnostic.threshold.eq(row.threshold) & diagnostic.side.eq('joint')]
        assert len(strata) == row.strata
        random = np.random.default_rng(row.seed)
        null = np.zeros(row.permutations, dtype=int)
        expected = 0.
        for cell in strata.itertuples(index=False):
            human_count, mouse_count = map(int, cell.stratum.split('|')[-2:])
            if human_count < 2:
                assert cell.human_events == 0
            if mouse_count < 2:
                assert cell.mouse_events == 0
            n, nx, ny = int(cell.n), int(cell.human_events), int(cell.mouse_events)
            assert (n, nx, ny) == (cell.n, cell.human_events, cell.mouse_events)
            expected += nx * ny / n
            null += hypergeom.rvs(n, ny, nx,
                                 size=row.permutations, random_state=random)
        selected = members[members.contrast.eq(row.contrast) & members.layer.eq(row.layer) &
                           members.threshold.eq(row.threshold)]
        assert len(selected) == row.observed
        assert abs(expected - row.analytic_null_mean) < 1e-10
        assert abs(null.mean() - row.null_mean) < 1e-12
        assert abs((1 + np.count_nonzero(null >= row.observed)) / (len(null) + 1) - row.empirical_p) < 1e-14
        for probability, column in [(.025, 'null_q025'), (.975, 'null_q975')]:
            assert np.quantile(null, probability) == getattr(row, column)
        archive[f'{row.contrast}|{row.layer}|{row.threshold}'] = null
    values = results.empirical_p.to_numpy()
    order = np.argsort(values)
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum(1, np.minimum.accumulate(
        (values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1])
    assert np.allclose(adjusted, results.q_16, rtol=0, atol=1e-14)
    np.savez_compressed(out / 'replayed-threshold-draws.npz', **archive)
    result = dict(status='PASS', settings=16, tests_skipped=0,
                  all_sensitivity_candidate_members_unchanged=True,
                  zero_events_outside_actual_measurement_eligibility=True,
                  joint_draws_replayed_from_stratum_margins=True, q_16_reconstructed=True,
                  old_PM_significant_settings=int((old_results.contrast.eq('PM') & old_results.q_16.lt(.05)).sum()),
                  corrected_PM_significant_settings=int((results.contrast.eq('PM') & results.q_16.lt(.05)).sum()))
    (out / 'sensitivity-verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    verify(args.root, args.config, args.out)
