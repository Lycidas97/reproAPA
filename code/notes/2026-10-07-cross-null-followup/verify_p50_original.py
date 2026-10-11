"""Verify corrected p50 primary tests against frozen masks and archived null draws."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent / 'ladder-results/apadetect-p50'
read = lambda path: pd.read_csv(path, sep='\t', float_precision='round_trip')
f = read(root / 'original-nulls/ortholog_matched_features.tsv.gz')
pd.testing.assert_frame_equal(f, read(root / 'analysis/cross_overlap/ortholog_matched_features.tsv.gz'),
                              check_dtype=False, rtol=0, atol=1e-12)
table = read(root / 'original-nulls/overlap_null_results.tsv')
draws = dict(np.load(root / 'original-nulls/overlap_null_draws.npz'))
fields = {
    'joint_measurement_PAS': None,
    'independent_q5_reference_capacity': ['effect_n', 'pas_q3', 'prox_capacity', 'ref_q5'],
    'independent_range10_reference_capacity': ['effect_n', 'pas_q3', 'prox_capacity', 'ref_min10', 'ref_max10'],
}
summary = []
for row in table.itertuples(index=False):
    co, model = row.contrast, row.model
    prob, legacy = {}, {}
    for side in ['h', 'm']:
        label = f'{side}_{co}'
        old_columns = ['h_n', 'm_n', 'h_pas', 'm_pas'] if fields[model] is None else [label + '_' + field for field in fields[model]]
        columns = old_columns + ([f'h_{co}_actual_effect_n', f'm_{co}_actual_effect_n'] if fields[model] is None else [label + '_actual_effect_n'])
        prob[side] = f.groupby(columns, dropna=False)[label].transform('mean').to_numpy()
        legacy[side] = f.groupby(old_columns, dropna=False)[label].transform('mean').to_numpy()
        ineligible = f[label + '_actual_effect_n'].lt(2).to_numpy()
        assert not prob[side][ineligible].any() and not f.loc[ineligible, label].any()
    expected = float(prob['h'] @ prob['m'])
    assert abs(expected - row.analytic_null_mean) < 1e-10
    null = draws[co + '|' + model]
    assert len(null) == row.permutations
    assert abs(null.mean() - row.null_mean) < 1e-12
    assert abs(null.std() - row.movable_null_sd) < 1e-12
    assert abs(expected - null.mean()) <= max(.01, 6 * null.std(ddof=1) / np.sqrt(len(null)))
    observed = int((f['h_' + co] & f['m_' + co]).sum())
    assert observed == row.observed
    p = (1 + (null >= observed).sum()) / (len(null) + 1)
    assert abs(p - row.empirical_p) < 1e-14
    assert np.array_equal(np.quantile(null, [.025, .975]), [row.null_q025, row.null_q975])
    summary.append(dict(contrast=co, model=model, observed=observed, legacy_expected=float(legacy['h'] @ legacy['m']),
                        expected=expected, enrichment=observed / expected, null_mean=row.null_mean,
                        p=p, q_two_contrasts=row.q_within_model_two_contrasts))
for _, part in table.groupby('model'):
    p = part.empirical_p.to_numpy()
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    assert np.allclose(part.q_within_model_two_contrasts, adjusted, rtol=0, atol=1e-14)
result = json.loads((root / 'original-nulls/result.json').read_text())
assert result['status'] == 'complete' and all(result['checks'].values())
for artifact in result['artifacts']:
    assert hashlib.sha256((root / 'original-nulls' / artifact['path']).read_bytes()).hexdigest() == artifact['sha256']
pd.DataFrame(summary).to_csv(root.parent.parent / 'p50-original-results.tsv', sep='\t', index=False, float_format='%.17g')
receipt = dict(status='PASS', tests_skipped=0, null_tests=6,
    features_and_actual_effect_counts_agree_with_independently_verified_ladder=True,
    original_event_members_and_universe_unchanged=True, event_probability_outside_eligibility_is_zero=True,
    analytic_means_archived_tail_probabilities_and_BH_verified=True, downloaded_artifact_hashes_verified=True)
(root.parent.parent / 'p50-original-verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
