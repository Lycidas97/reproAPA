"""Frozen A3 redesign: two fixed models, exact directions, marginal power."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq
from scipy.special import gammaln
from scipy.stats import beta, hypergeom, nchypergeom_fisher

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT/'reviews/2026-10-08-plan-a-analysis/cross-species/results'
FAMILY = 'A3-redesign-cross-species'
SOURCES = ['polyapipe', 'apadetect-p50']
MODELS = ['D1_logical_capacity', 'D1_PAS3_measurement']
ODDS = [1, 1.5, 2, 3, 4, 6, 8, 12, 16, 32, 64]
B, POWER_N, SEED = 20000, 3000, 20261008
INPUTS = {}


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read(p):
    INPUTS[str(p)] = sha(p)
    return pd.read_csv(p, sep='\t', float_precision='round_trip')


def save(rows, name):
    pd.DataFrame(rows).to_csv(HERE/name, sep='\t', index=False,
                            na_rep='NA', float_format='%.17g')


def rng_for(key):
    seed = int.from_bytes(hashlib.sha256(f'{SEED}|redesign|{key}'.encode()).digest()[:8], 'little')
    return np.random.default_rng(seed), seed


def groups(f, side, model):
    n_routes = 1 if side == 'h' else 2
    fields = pd.DataFrame({
        'can_proximal': f[side+'_S_prox_capacity'].eq(n_routes) & f[side+'_M_prox_capacity'].eq(n_routes),
        'can_distal': f[side+'_S_dist_capacity'].eq(n_routes) & f[side+'_M_dist_capacity'].eq(n_routes),
    })
    if model == 'D1_PAS3_measurement':
        fields['PAS_Q3'] = f[side+'_pas_q3']
    return [(np.asarray(ii, int), '|'.join(map(str, key if isinstance(key, tuple) else (key,))))
            for key, ii in fields.groupby(list(fields), sort=True).indices.items()]


def validate_universe(f, source):
    assert len(f) == (678 if source == 'polyapipe' else 525)
    assert f.human_gene_id.is_unique and f.mouse_gene_id.is_unique
    assert f.homology_type.eq('ortholog_one2one').all()
    assert f.human_partner_count.eq(1).all() and f.mouse_partner_count.eq(1).all()
    for side in ['h', 'm']:
        n_routes = 1 if side == 'h' else 2
        assert set(f[side+'_label']) <= {-1, 0, 1}
        for name in ['S', 'M']:
            assert f.loc[f[side+'_label'].eq(-1), side+'_'+name+'_prox_capacity'].eq(n_routes).all()
            assert f.loc[f[side+'_label'].eq(1), side+'_'+name+'_dist_capacity'].eq(n_routes).all()
        # Reconstruct the frozen PAS quantiles, without changing their definition.
        q = pd.qcut(f[side+'_PAS_median'], 3, labels=False, duplicates='drop').fillna(0).astype(int)
        np.testing.assert_array_equal(q, f[side+'_pas_q3'])
    refs = [c for c in f if c.endswith(('_S_DWUI', '_M_DWUI', '_PAS_number'))]
    assert np.isfinite(f[refs]).all().all()


def draw_labels(labels, strata, n, rng):
    result = np.zeros((n, len(labels)), np.int8)
    for ii, _ in strata:
        values = labels[ii]
        if np.unique(values).size == 1:
            result[:, ii] = values[0]
        else:
            result[:, ii] = values[np.argsort(rng.random((n, len(ii))), axis=1)]
        for value in [-1, 0, 1]:
            assert (result[:, ii] == value).sum(axis=1).min() == int((values == value).sum())
            assert (result[:, ii] == value).sum(axis=1).max() == int((values == value).sum())
    return result


def overlap_null(x, y, hg, mg, rng):
    draws = []
    for start in range(0, B, 500):
        n = min(500, B-start)
        h, m = draw_labels(x, hg, n, rng), draw_labels(y, mg, n, rng)
        draws.append(((h != 0) & (m != 0)).sum(axis=1))
    return np.concatenate(draws)


def mc_p(null, statistic):
    return (1+len(null)-np.searchsorted(np.sort(null), statistic, side='left'))/(1+len(null))


def binomial(k, n):
    return dict(power=k/n,
                power_ci_low=float(beta.ppf(.025, k, n-k+1)) if k else 0.,
                power_ci_high=float(beta.ppf(.975, k+1, n-k)) if k<n else 1.)


def configurations(f, side, strata, hits):
    labels = f[side+'_label'].to_numpy(int)
    patterns = np.arange(2**len(hits), dtype=np.uint64)
    logprob = np.zeros(len(patterns))
    possible = np.ones(len(patterns), bool)
    probability_membership = np.zeros(len(f))
    diagnostics = []
    movable = 0
    for ii, signature in strata:
        selected = ii[labels[ii] != 0]
        n, k = len(selected), int((labels[selected] == -1).sum())
        probability_membership[ii] = n/len(ii)
        positions = [j for j, index in enumerate(hits) if index in set(selected)]
        r = len(positions)
        if r:
            bits = np.uint64(sum(1 << j for j in positions))
            a = np.bitwise_count(patterns & bits).astype(int)
            left = k-a
            valid = (left >= 0) & (left <= n-r)
            possible &= valid
            logprob[valid] += (gammaln(n-r+1)-gammaln(left[valid]+1)-gammaln(n-r-left[valid]+1)
                               -gammaln(n+1)+gammaln(k+1)+gammaln(n-k+1))
        movable += r if 0<k<n else 0
        diagnostics.append(dict(side=side, signature=signature, background_n=len(ii), members=n,
                                proximal=k, distal=n-k, nonmember=len(ii)-n, shared=r,
                                membership_movable=0<n<len(ii), direction_movable=0<k<n))
    probability = np.where(possible, np.exp(logprob), 0.)
    np.testing.assert_allclose(probability.sum(), 1., rtol=0, atol=1e-10)
    keep = probability > 0
    for j, index in enumerate(hits):
        prox = ((patterns[keep] >> np.uint64(j)) & np.uint64(1)) != 0
        n_routes = 1 if side == 'h' else 2
        if prox.any():
            assert f.loc[index, side+'_S_prox_capacity'] == n_routes
            assert f.loc[index, side+'_M_prox_capacity'] == n_routes
        if (~prox).any():
            assert f.loc[index, side+'_S_dist_capacity'] == n_routes
            assert f.loc[index, side+'_M_dist_capacity'] == n_routes
    return patterns[keep], probability[keep], probability_membership, diagnostics, movable


def exact_direction(h, hp, m, mp, shared):
    pmf = np.zeros(shared+1)
    for start in range(0, len(h), 128):
        t = shared-np.bitwise_count(h[start:start+128, None] ^ m[None, :]).astype(int)
        weights = hp[start:start+128, None]*mp[None, :]
        pmf += np.bincount(t.ravel(), weights=weights.ravel(), minlength=shared+1)
    np.testing.assert_allclose(pmf.sum(), 1., rtol=0, atol=1e-10)
    return pmf


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum.accumulate((p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p));result[order] = np.minimum(adjusted, 1.)
    return result


def conditional_bh_mask(values, observed, position):
    result = []
    for value in values:
        p = observed.copy();p[position] = value
        result.append(bh(p)[position] <= .05)
    return np.asarray(result, bool)


def tilt(pmf, odds):
    t = np.arange(len(pmf))
    weights = pmf*np.exp((t-t.max())*np.log(odds))
    return weights/weights.sum()


def overlap_alternative(x, y, hg, mg, odds, rng):
    result = []
    for start in range(0, POWER_N, 500):
        n = min(500, POWER_N-start)
        h = draw_labels(x, hg, n, rng) != 0
        overlap = np.zeros(n, int)
        for ii, _ in mg:
            N, selected = len(ii), int((y[ii] != 0).sum())
            K = h[:, ii].sum(axis=1)
            for val in np.unique(K):
                use = K == val
                low, high = max(0, selected+int(val)-N), min(selected, int(val))
                if low == high:
                    overlap[use] += low
                else:
                    dist = hypergeom if odds == 1 else nchypergeom_fisher
                    args = (N, int(val), selected) if odds == 1 else (N, int(val), selected, odds)
                    overlap[use] += dist.rvs(*args, size=int(use.sum()), random_state=rng)
        result.append(overlap)
    return np.concatenate(result)


def main():
    parser = argparse.ArgumentParser();parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if not args.execute:
        parser.error('Use --execute only after root authorizes the frozen total protocol.')
    INPUTS[str(HERE/'PROTOCOL.zh.md')] = sha(HERE/'PROTOCOL.zh.md')
    total_protocol = ROOT/'reviews/2026-10-08-plan-a-redesign/PROTOCOL.md'
    INPUTS[str(total_protocol)] = sha(total_protocol)
    old_tests = read(OLD/'raw-pvalues.tsv')
    old_exact = read(OLD/'direction-exact-support.tsv')
    universes = {source:read(OLD/(source+'-universe.tsv.gz')) for source in SOURCES}
    for source, f in universes.items():
        validate_universe(f, source)
    (HERE/'input-audit.json').write_text(json.dumps(dict(status='frozen_before_new_inference',
        frozen_utc=datetime.now(timezone.utc).isoformat(), root_execution_authorized=True,
        input_sha256=INPUTS.copy(), code_sha256=sha(__file__), verification_code_sha256=sha(HERE/'verify.py'), master_seed=SEED,
        null_draws=B, overlap_power_draws=POWER_N, direction_power='exact',
        python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
        scipy=scipy.__version__), indent=2)+'\n')
    tests, strata_rows, pmf_rows, seeds, contexts, archives = [], [], [], [], [], {}
    for source, f in universes.items():
        x, y = f.h_label.to_numpy(np.int8), f.m_label.to_numpy(np.int8)
        hits = np.flatnonzero((x != 0) & (y != 0))
        shared = f.iloc[hits].copy();shared['same_direction'] = shared.h_label.eq(shared.m_label)
        save(shared, source+'-shared-members.tsv')
        for model in MODELS:
            hg, mg = groups(f, 'h', model), groups(f, 'm', model)
            h, hp, hprob, hs, hm = configurations(f, 'h', hg, hits)
            m, mp, mprob, ms, mm = configurations(f, 'm', mg, hits)
            strata_rows.extend(dict(pas_source=source, model=model, **r) for r in hs+ms)
            rng, seed = rng_for(f'{source}|{model}|overlap|null')
            seeds.append(dict(pas_source=source, model=model, purpose='overlap_null', odds_ratio=1, seed=seed))
            null = overlap_null(x, y, hg, mg, rng)
            expected = float(hprob @ mprob)
            assert abs(null.mean()-expected) <= max(.02, 6*null.std(ddof=1)/np.sqrt(B))
            pmf = exact_direction(h, hp, m, mp, len(hits));tails = np.clip(np.cumsum(pmf[::-1])[::-1], 0., 1.)
            support = np.flatnonzero(pmf > 0)
            same = int((x[hits] == y[hits]).sum())
            assert pmf[same] > 0
            base = dict(family=FAMILY, pas_source=source, model=model, universe=len(f),
                        human_members=int((x != 0).sum()), mouse_members=int((y != 0).sum()),
                        overlap=len(hits), same_direction=same, h_shared_direction_movable=hm,
                        m_shared_direction_movable=mm, status='complete')
            k = int((null >= len(hits)).sum())
            interval = binomial(k, B)
            tests.append(dict(**base, test='overlap_enrichment',
                test_id=f'cross|{source}|{model}|overlap', observed=len(hits), expected=expected,
                raw_p=float(mc_p(null, len(hits))), method='independent_stratified_permutation',
                null_draws=B, mc_tail_ci_low=interval['power_ci_low'], mc_tail_ci_high=interval['power_ci_high'],
                null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                sampled_null_min=int(null.min()), sampled_null_max=int(null.max()), minimum_attainable_p=np.nan))
            tests.append(dict(**base, test='direction_concordance',
                test_id=f'cross|{source}|{model}|direction', observed=same,
                expected=float(np.arange(len(pmf)) @ pmf), raw_p=float(tails[same]),
                method='exact_fixed_members_signed_margins', null_draws=0,
                null_mean=float(np.arange(len(pmf)) @ pmf), null_sd=float(np.sqrt(((np.arange(len(pmf))-np.arange(len(pmf))@pmf)**2)@pmf)),
                support_min=int(support.min()), support_max=int(support.max()), minimum_attainable_p=float(pmf[support.max()]),
                human_shared_configurations=len(h), mouse_shared_configurations=len(m)))
            pmf_rows.extend(dict(pas_source=source, model=model, same_direction=t, null_probability=p,
                                 exact_upper_tail=tails[t]) for t, p in enumerate(pmf))
            archives[f'{source}|{model}|null_overlap'] = null
            contexts.append(dict(source=source, model=model, x=x, y=y, hg=hg, mg=mg,
                                 null=null, pmf=pmf, tails=tails))
    tests = pd.DataFrame(tests);tests['q_BH_8'] = bh(tests.raw_p);tests['reject_BH_005'] = tests.q_BH_8.le(.05)
    assert len(tests) == 8 and tests.test_id.is_unique and tests.raw_p.notna().all()
    save(tests, 'corrected-tests.tsv');save(strata_rows, 'stratum-diagnostics.tsv');save(pmf_rows, 'direction-exact-pmf.tsv')
    print(tests.to_string(index=False), flush=True)
    frozen_p = tests.raw_p.to_numpy().copy()
    powers, detectable = [], []
    for c in contexts:
        source, model = c['source'], c['model']
        for test in ['overlap', 'direction']:
            position = int(np.flatnonzero(tests.test_id.eq(f'cross|{source}|{model}|{test}'))[0])
            p_values = c['tails'] if test == 'direction' else mc_p(c['null'], np.arange(len(c['x'])+1))
            rejection = {'nominal_005_upper_bound':p_values <= .05,
                         'bonferroni8_000625_lower_bound':p_values <= .00625,
                         'conditional_BH_fixed_other_7_observed_p':conditional_bh_mask(p_values, frozen_p, position)}
            assert np.all(~rejection['bonferroni8_000625_lower_bound'] | rejection['conditional_BH_fixed_other_7_observed_p'])
            assert np.all(~rejection['conditional_BH_fixed_other_7_observed_p'] | rejection['nominal_005_upper_bound'])
            for odds in ODDS:
                if test == 'direction':
                    distribution = tilt(c['pmf'], odds)
                    if odds == 1:
                        np.testing.assert_allclose(distribution, c['pmf'], atol=1e-12)
                    expected = float(np.arange(len(distribution)) @ distribution)
                    for rule, mask in rejection.items():
                        powers.append(dict(pas_source=source, model=model, test=test, odds_ratio=odds, rule=rule,
                            generator='both_species_joint_OR_to_T', expected_statistic=expected,
                            power=float(distribution[mask].sum()), simulations=0, exact=True))
                else:
                    rng, seed = rng_for(f'{source}|{model}|overlap|power|{odds}')
                    seeds.append(dict(pas_source=source, model=model, purpose='overlap_power', odds_ratio=odds, seed=seed))
                    draws = overlap_alternative(c['x'], c['y'], c['hg'], c['mg'], odds, rng)
                    archives[f'{source}|{model}|power_overlap_OR{odds}'] = draws
                    if odds == 1:
                        assert abs(draws.mean()-c['null'].mean()) <= 6*np.sqrt(draws.var(ddof=1)/POWER_N+c['null'].var(ddof=1)/B)
                    for rule, mask in rejection.items():
                        k = int(mask[draws].sum())
                        powers.append(dict(pas_source=source, model=model, test=test, odds_ratio=odds, rule=rule,
                            generator='human_null_then_mouse_conditional_Fisher_noncentral',
                            expected_statistic=float(draws.mean()), rejections=k, simulations=POWER_N,
                            exact=False, **binomial(k, POWER_N)))
            for rule, mask in rejection.items():
                rows = [r for r in powers if r['pas_source']==source and r['model']==model and r['test']==test and r['rule']==rule]
                passing = [r for r in rows if r['power'] >= .8]
                reachable = bool((mask & (c['pmf'] > 0)).any()) if test == 'direction' else bool(mask.any())
                exact_or = np.nan
                if test == 'direction' and reachable:
                    fun = lambda logodds: float(tilt(c['pmf'], np.exp(logodds))[mask].sum())-.8
                    if fun(np.log(1e6)) >= 0:
                        exact_or = float(np.exp(brentq(fun, 0., np.log(1e6))))
                detectable.append(dict(pas_source=source, model=model, test=test, rule=rule,
                    first_grid_OR_80=passing[0]['odds_ratio'] if passing else np.nan,
                    exact_OR_80=exact_or, exact_OR_search_upper=1e6 if test=='direction' else np.nan,
                    grid_max_OR=64, power_at_grid_max=rows[-1]['power'],
                    attainable_rejection=reachable if test=='direction' else 'not_computed_exact_overlap_support'))
    save(powers, 'power-curves.tsv');save(detectable, 'detectable-effects.tsv');save(seeds, 'seeds.tsv')
    comparison = old_tests.copy();comparison['comparison_scope'] = 'D5_frozen_previous_analysis_not_in_new_family'
    comparison['q_BH_previous4'] = bh(comparison.raw_p)
    previous_exact = old_exact.set_index('pas_source').exact_upper_tail_at_observed
    comparison['previous_exact_direction_p'] = comparison.pas_source.map(previous_exact).where(comparison.test.eq('direction_concordance'))
    save(comparison, 'D5-old-results-comparison.tsv')
    # Include old PMFs/support untouched even if their legacy column naming differs.
    save(old_exact, 'D5-old-exact-direction-support.tsv')
    np.savez_compressed(HERE/'simulation-draws.npz', **archives)
    assert all(sha(p) == fingerprint for p, fingerprint in INPUTS.items())
    output_names = ['direction-exact-pmf.tsv', 'power-curves.tsv', 'corrected-tests.tsv',
                    'stratum-diagnostics.tsv', 'detectable-effects.tsv', 'seeds.tsv',
                    'D5-old-results-comparison.tsv', 'D5-old-exact-direction-support.tsv',
                    'simulation-draws.npz', 'input-audit.json', 'PROTOCOL.zh.md',
                    'run_cross_species.py', 'verify.py'] + [s+'-shared-members.tsv' for s in SOURCES]
    outputs = {name:sha(HERE/name) for name in output_names}
    receipt = dict(status='PASS', tests=8, tests_skipped=0, family=FAMILY, input_sha256=INPUTS,
                   original_inputs_unchanged=True, code_sha256=sha(__file__), outputs_sha256=outputs,
                   actual_joint_BH_power_computed=False, conditional_BH_other_p_frozen=True)
    (HERE/'result.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
