"""Independent replay of saved allocation statistics and matching signatures."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze import add_bh, inference

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
inputs = json.loads((HERE/'input_sha256.json').read_text())
assert all(sha(ROOT/p) == digest for p, digest in inputs.items())
frozen = pd.read_csv(ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/recurrence_controls/recurrence_frozen_covariates.tsv.gz', sep='\t')
draws = np.load(HERE/'matched_component_null_draws.npz')
tests = pd.read_csv(HERE/'matched_component_tests.tsv', sep='\t', float_precision='round_trip')
interactions = pd.read_csv(HERE/'matched_interaction_tests.tsv', sep='\t', float_precision='round_trip')
assert len(draws.files) == 24 and len(tests) == 48 and len(interactions) == 24
for key in draws.files:
    sp, co, raw_threshold, level = key.split('|')
    threshold = float(raw_threshold)
    file = HERE/f'background-{sp}-{co}-{raw_threshold.replace(".", "p")}-{level}.tsv.gz'
    frame = pd.read_csv(file, sep='\t', float_precision='round_trip')
    native = frozen[frozen.species.eq(sp)&frozen.contrast.eq(co)&frozen.threshold.eq(threshold)]
    signature = native.groupby('gene_id').apply(lambda q: ';'.join(
        f'{r.unit}:{r.series_pattern}' for r in q.sort_values('unit').itertuples()), include_groups=False)
    assert frame.availability.eq(frame.gene_id.map(signature)).all()
    values, labels = frame.A.to_numpy(), frame.label.to_numpy()
    selected = tests[tests.species.eq(sp)&tests.selection_contrast.eq(co)&tests.threshold.eq(threshold)&tests.null_model.eq(level)]
    for i, (direction, sign) in enumerate([('proximal', -1), ('distal', 1)]):
        row = selected[selected.direction.eq(direction)].iloc[0]
        count = int((labels == sign).sum())
        observed = values[labels == sign].mean()
        assert count == row.n and np.isclose(observed, row.observed, atol=1e-12, rtol=0)
        expected = sum((q.label == sign).sum()*q.A.mean() for _, q in frame.groupby('stratum'))/count
        assert np.isclose(expected, row.analytic_null_mean, atol=1e-12, rtol=0)
        result = inference(observed, draws[key][:, i])
        for column in ['p_upper', 'null_mean', 'null_sd', 'null_p025', 'null_p975', 'tail_count']:
            assert np.isclose(result[column], row[column], atol=1e-12, rtol=0, equal_nan=True), (key, column)
    row = interactions[interactions.species.eq(sp)&interactions.selection_contrast.eq(co)&interactions.threshold.eq(threshold)&interactions.null_model.eq(level)].iloc[0]
    observed = values[labels == 1].mean()+values[labels == -1].mean()
    result = inference(observed, draws[key].sum(1))
    assert np.isclose(observed, row.observed, atol=1e-12, rtol=0)
    assert np.isclose(result['p_upper'], row.p_upper, atol=1e-12, rtol=0, equal_nan=True)
for frame in [tests, interactions]:
    recovered = add_bh(frame.drop(columns=['q_BH_per_threshold', 'q_BH_joint_thresholds',
        'BH_family_size_per_threshold', 'BH_family_size_joint_thresholds']))
    for column in ['q_BH_per_threshold', 'q_BH_joint_thresholds']:
        np.testing.assert_allclose(frame[column], recovered[column], atol=1e-12, rtol=0, equal_nan=True)
effects = pd.read_csv(HERE/'fixed_cohort_component_effects.tsv.gz', sep='\t')
assert np.allclose(effects.DWUI_PS, effects.DWUI_PM+effects.DWUI_MS, atol=1e-12, rtol=0)
assert np.allclose(effects.aligned_allocation, effects.direction_sign*(effects.DWUI_MS-effects.DWUI_PM), atol=1e-12, rtol=0)
verification = pd.read_csv(HERE/'verification.tsv', sep='\t')
assert len(verification) == 7 and verification.max_abs_error.lt(1e-12).all()
assert 'Ran 5 tests' in (HERE/'tests.log').read_text() and (HERE/'tests.log').read_text().rstrip().endswith('OK')
result = dict(status='PASS', input_hashes_verified=len(inputs), null_contexts_verified=24,
    actual_series_identity_preserved=True, class_tests_replayed=48, interaction_tests_replayed=24,
    analytic_expectations_and_BH_replayed=True, matched_effect_addition_and_alignment_verified=True,
    upstream_vector_and_weight_checks=7, tests_run=5, tests_skipped=0)
(HERE/'independent_validation.json').write_text(json.dumps(result, indent=2)+'\n')
outputs = sorted(p for p in HERE.iterdir() if p.is_file() and p.name not in ['manifest.json', 'verify.log'])
(HERE/'manifest.json').write_text(json.dumps(dict(status='PASS', input_sha256=inputs,
    output_sha256={p.name: sha(p) for p in outputs}, independent_checks=result,
    post_run_method_fix='Availability signature records actual series identity, not only the number of series; all contexts regenerated with unchanged scientific definitions.'), indent=2)+'\n')
print(json.dumps(result))
