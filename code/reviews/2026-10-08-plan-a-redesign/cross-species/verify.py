"""Independent receipt, probability, BH, and power checks for the A3 redesign."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read(name):
    return pd.read_csv(HERE/name, sep='\t', float_precision='round_trip')


def main():
    receipt = json.loads((HERE/'result.json').read_text())
    assert receipt['status'] == 'PASS' and receipt['tests'] == 8 and receipt['tests_skipped'] == 0
    assert all(sha(p) == value for p, value in receipt['input_sha256'].items())
    assert all(sha(HERE/name) == value for name, value in receipt['outputs_sha256'].items())
    tests, strata, pmfs, powers = [read(n) for n in ['corrected-tests.tsv', 'stratum-diagnostics.tsv', 'direction-exact-pmf.tsv', 'power-curves.tsv']]
    assert len(tests) == 8 and tests.test_id.is_unique
    assert set(tests.family) == {'A3-redesign-cross-species'}
    np.testing.assert_allclose(tests.q_BH_8, false_discovery_control(tests.raw_p, method='bh'), atol=1e-12)
    np.testing.assert_array_equal(tests.reject_BH_005, tests.q_BH_8 <= .05)
    arrays = np.load(HERE/'simulation-draws.npz')
    old_diagnostic = pd.read_csv(ROOT/'reviews/2026-10-08-plan-a-statistical-strictness/cross-species/direction-exact-pmf-diagnostic.tsv', sep='\t', float_precision='round_trip')
    frozen_p = tests.raw_p.to_numpy()
    checked_power_rows = 0
    for (source, model), q in tests.groupby(['pas_source', 'model']):
        f = pd.read_csv(ROOT/'reviews/2026-10-08-plan-a-analysis/cross-species/results'/(source+'-universe.tsv.gz'), sep='\t', float_precision='round_trip')
        assert len(f) == q.universe.iloc[0]
        for side in ['h', 'm']:
            rows = strata[(strata.pas_source == source)&(strata.model == model)&(strata.side == side)]
            assert rows.background_n.sum() == len(f)
            assert rows.members.sum() == int(f[side+'_label'].ne(0).sum())
            assert rows.proximal.sum() == int(f[side+'_label'].eq(-1).sum())
            assert rows.distal.sum() == int(f[side+'_label'].eq(1).sum())
            assert (rows.members+rows.nonmember == rows.background_n).all()
            assert (rows.proximal+rows.distal == rows.members).all()
            assert rows.shared.sum() == q.overlap.iloc[0]
            # Sensitivity must add only PAS Q3, not route-count fields.
            assert all(len(s.split('|')) == (2 if model == 'D1_logical_capacity' else 3) for s in rows.signature)
        shared = read(source+'-shared-members.tsv')
        assert len(shared) == int(q.overlap.iloc[0]) and shared.human_gene_id.is_unique
        assert shared.same_direction.sum() == q.same_direction.iloc[0]
        observed_overlap = int(q.overlap.iloc[0])
        null = arrays[f'{source}|{model}|null_overlap']
        assert null.shape == (20000,) and null.min() >= 0
        row = q[q.test.eq('overlap_enrichment')].iloc[0]
        np.testing.assert_allclose(row.raw_p, (1+(null >= observed_overlap).sum())/20001, atol=1e-15)
        exact = pmfs[(pmfs.pas_source == source)&(pmfs.model == model)].sort_values('same_direction')
        np.testing.assert_array_equal(exact.same_direction, np.arange(observed_overlap+1))
        p0 = exact.null_probability.to_numpy();tails = np.clip(np.cumsum(p0[::-1])[::-1], 0., 1.)
        np.testing.assert_allclose(p0.sum(), 1., atol=1e-10)
        np.testing.assert_allclose(tails, exact.exact_upper_tail, atol=1e-12)
        direction = q[q.test.eq('direction_concordance')].iloc[0]
        np.testing.assert_allclose(direction.raw_p, tails[int(direction.observed)], atol=1e-12)
        support = np.flatnonzero(p0 > 0)
        assert direction.support_min == support.min() and direction.support_max == support.max()
        np.testing.assert_allclose(direction.minimum_attainable_p, p0[support.max()], atol=1e-12)
        if model == 'D1_logical_capacity':
            old = old_diagnostic[old_diagnostic.pas_source.eq(source)&old_diagnostic.model.eq('D1_logical_capacity')]
            np.testing.assert_allclose(p0, old.sort_values('same_direction').null_probability, atol=1e-10)
        for test in ['overlap', 'direction']:
            position = int(np.flatnonzero(tests.test_id.eq(f'cross|{source}|{model}|{test}'))[0])
            for odds in [1, 1.5, 2, 3, 4, 6, 8, 12, 16, 32, 64]:
                selected = powers[(powers.pas_source == source)&(powers.model == model)&(powers.test == test)&(powers.odds_ratio == odds)]
                assert len(selected) == 3
                if test == 'direction':
                    t = np.arange(len(p0));w = p0*np.exp((t-t.max())*np.log(odds));w /= w.sum()
                    p_values = tails
                else:
                    draws = arrays[f'{source}|{model}|power_overlap_OR{odds}']
                    assert draws.shape == (3000,)
                    t, counts = np.unique(draws, return_counts=True);w = counts/3000
                    p_values = (1+20000-np.searchsorted(np.sort(null), t, side='left'))/20001
                conditional = []
                for value in p_values:
                    pv = frozen_p.copy();pv[position] = value
                    conditional.append(false_discovery_control(pv, method='bh')[position] <= .05)
                masks = {'nominal_005_upper_bound':p_values <= .05,
                         'bonferroni8_000625_lower_bound':p_values <= .00625,
                         'conditional_BH_fixed_other_7_observed_p':np.asarray(conditional)}
                for r in selected.itertuples():
                    np.testing.assert_allclose(r.power, w[masks[r.rule]].sum(), atol=1e-10)
                    np.testing.assert_allclose(r.expected_statistic, t @ w, atol=1e-10)
                    if test == 'direction':
                        assert r.simulations == 0 and r.exact
                    else:
                        assert r.simulations == 3000 and not r.exact
                        assert r.power_ci_low <= r.power <= r.power_ci_high
                    checked_power_rows += 1
                values = selected.set_index('rule').power
                assert values['bonferroni8_000625_lower_bound'] <= values['conditional_BH_fixed_other_7_observed_p']+1e-12
                assert values['conditional_BH_fixed_other_7_observed_p'] <= values['nominal_005_upper_bound']+1e-12
    assert checked_power_rows == 264 and len(powers) == 264
    result = dict(status='PASS', tests_checked=8, power_rows_checked=264, checks_skipped=0,
                  original_inputs_unchanged=True, output_hashes_match=True,
                  BH_independent_scipy_match=True, D1_existing_exact_diagnostic_match=True,
                  joint_family_power_not_claimed=True, verification_code_sha256=sha(__file__))
    (HERE/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
