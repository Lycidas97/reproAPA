#!/usr/bin/env python3
"""Frozen-source A1 background and A7 third-reference follow-up."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, hypergeom
from scipy.special import gammaln

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASES = {
    'polyapipe': ROOT / 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables',
    'apadetect-p50': ROOT / 'notes/2026-10-06-dual-main-plot-data/apadetect-p50/supplement/tables',
}
ROUTES = {'human': ['GSE118258'], 'mouse': ['CIP', 'WOT']}
TAU = .02
B = 20000
INPUTS = {}


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read(base, relative, **kwargs):
    path = base / relative
    INPUTS[str(path)] = {'sha256': digest(path), 'bytes': path.stat().st_size}
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, relative):
    kwargs = {'compression': {'method': 'gzip', 'mtime': 0}} if relative.endswith('.gz') else {}
    frame.to_csv(HERE / relative, sep='\t', index=False, na_rep='NA', float_format='%.17g', **kwargs)


def qbin(values, q=3):
    # Identical rank convention to the frozen recurrence null.
    assert values.notna().all()
    if values.nunique() < 2:
        return pd.Series('0', index=values.index)
    return (np.minimum(np.ceil(values.rank(method='average', pct=True) * q).astype(int), q) - 1).astype(str)


def label(ps, pm):
    return np.where((ps >= TAU) & (pm >= TAU), 1,
                    np.where((ps <= -TAU) & (pm <= -TAU), -1, 0))


def recurrence(base):
    r = read(base, 'recurrence/recurrence_members.tsv.gz')
    r = r[r.metric.eq('WPP') & r.threshold.eq(TAU) & r.series_matching & r.contrast.eq('P_M')].copy()
    assert not r.duplicated(['species', 'gene_id']).any()
    r['recurrent_sign'] = r.direction.map({'negative': -1, 'positive': 1})
    return r


def background(pas, base):
    m = read(base, 'matched_pms/matched_state_vectors.tsv.gz')
    m = m[m.setting.eq('supported50_p08')].copy()
    r = recurrence(base)
    dual = read(base, 'candidate_screen/dual_reference_candidates.tsv')
    ladder = read(base, 'candidate_screen/dual_reference_ladder.tsv')
    raw = read(base, 'state_vectors/route_effects.primary.tsv.gz')
    cov = raw[raw.contrast.eq('P_M')][['species', 'route', 'gene_id', 'mean_GEX', 'APA_coverage']]
    assert not cov.duplicated(['species', 'route', 'gene_id']).any()
    m = m.merge(cov, on=['species', 'route', 'gene_id'], validate='one_to_one')
    assert np.max(np.abs(m.DWUI_PS - m.DWUI_PM - m.DWUI_MS)) < 1e-10
    out = []
    for sp, routes in ROUTES.items():
        z = m[m.species.eq(sp) & m.route.isin(routes)].copy()
        available = z.groupby('gene_id').route.agg(set)
        ids = available[available.map(lambda a: a == set(routes))].index
        z = z[z.gene_id.isin(ids)].copy()
        assert not z.duplicated(['gene_id', 'route']).any()
        parts = []
        for rt in routes:
            a = z[z.route.eq(rt)].copy().set_index('gene_id').sort_index()
            a['route_alignment_sign'] = label(a.DWUI_PS, a.DWUI_PM)
            fields = []
            for c in ['mean_GEX', 'APA_coverage', 'PAS_number']:
                a[c + '_q3'] = qbin(a[c]); fields.append(c + '_q3')
            a['M_baseline_q5'] = qbin(a.M_DWUI, 5); fields.append('M_baseline_q5')
            for state in ['M', 'S']:
                a[state + '_can_negative'] = a[state + '_DWUI'].ge(TAU - 1e-9)
                a[state + '_can_positive'] = a[state + '_DWUI'].le(1 - TAU + 1e-9)
                fields.extend([state + '_can_negative', state + '_can_positive'])
            a['route_stratum'] = a[fields].astype(str).agg('|'.join, axis=1)
            cols = ['DWUI_PS', 'DWUI_PM', 'DWUI_MS', 'P_DWUI', 'M_DWUI', 'S_DWUI',
                    'mean_GEX', 'APA_coverage', 'PAS_number', 'n_triplets', 'n_target_libraries',
                    'n_source_libraries', 'n_target_dates', 'route_alignment_sign', 'route_stratum']
            parts.append(a[cols].add_prefix(rt + '__'))
        f = pd.concat(parts, axis=1).reset_index()
        f['species'] = sp; f['pas_source'] = pas
        f['gene_symbol'] = f.gene_id.map(z.drop_duplicates('gene_id').set_index('gene_id').gene_symbol)
        signs = f[[rt + '__route_alignment_sign' for rt in routes]].to_numpy()
        f['alignment_sign'] = np.where((signs == 1).all(1), 1, np.where((signs == -1).all(1), -1, 0))
        ps = f[[rt + '__DWUI_PS' for rt in routes]].to_numpy()
        pm = f[[rt + '__DWUI_PM' for rt in routes]].to_numpy()
        f['PS_only_alignment_sign'] = np.where((ps >= TAU).all(1), 1, np.where((ps <= -TAU).all(1), -1, 0))
        f['PM_only_alignment_sign'] = np.where((pm >= TAU).all(1), 1, np.where((pm <= -TAU).all(1), -1, 0))
        f['aligned_any_direction'] = f.alignment_sign.ne(0)
        rr = r[r.species.eq(sp)].set_index('gene_id')
        f['is_pm_recurrent'] = f.gene_id.isin(rr.index)
        f['recurrent_sign'] = f.gene_id.map(rr.recurrent_sign).fillna(0).astype(int)
        f['is_dual_reference'] = f.is_pm_recurrent & f.alignment_sign.eq(f.recurrent_sign) & f.alignment_sign.ne(0)
        f['opposite_recurrent_alignment'] = f.is_pm_recurrent & f.alignment_sign.eq(-f.recurrent_sign) & f.alignment_sign.ne(0)
        f['stratum'] = f[[rt + '__route_stratum' for rt in routes]].astype(str).agg('||'.join, axis=1)
        expected = ladder[ladder.species.eq(sp)].iloc[0]
        assert int(f.is_pm_recurrent.sum()) == expected.matched_measurable
        assert int(f.is_dual_reference.sum()) == expected.dual_reference
        assert set(f.loc[f.is_dual_reference, 'gene_id']) == set(dual.loc[dual.species.eq(sp), 'gene_id'])
        assert int(r.species.eq(sp).sum()) == expected.recurrent_candidates
        out.append(f)
    return pd.concat(out, ignore_index=True), m, r


def fisher_row(name, family, pas, sp, a, b, c, d, **extra):
    odds, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    return dict(test_id=f'{pas}|{sp}|{name}', family=family, pas_source=pas, species=sp,
                test=name, method='Fisher exact one-sided enrichment', p_value=float(p),
                odds_ratio=float(odds), selected_positive=int(a), selected_negative=int(b),
                comparator_positive=int(c), comparator_negative=int(d), **extra)


def logcomb(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def signed_pmf(npos, nneg, nzero, kpos, kneg):
    """Exact signed-candidate allocation; opposite events cannot count as support."""
    n = npos + nneg + nzero
    pmf = np.zeros(kpos + kneg + 1)
    for xp in range(min(kpos, npos) + 1):
        for xm in range(min(kpos - xp, nneg) + 1):
            x0 = kpos - xp - xm
            if not 0 <= x0 <= nzero:
                continue
            prob = float(np.exp(logcomb(npos, xp) + logcomb(nneg, xm) + logcomb(nzero, x0) - logcomb(n, kpos)))
            ym = np.arange(kneg + 1)
            if n == kpos:
                assert kneg == 0
                negative = np.array([1.])
            else:
                negative = hypergeom.pmf(ym, n - kpos, nneg - xm, kneg)
            pmf[xp:xp + kneg + 1] += prob * negative
    assert abs(pmf.sum() - 1) < 1e-9
    return pmf


def matched_null(f, pas, sp, model='reference_capacity'):
    distribution = np.array([1.])
    observed = int(f.is_dual_reference.sum())
    seed = int.from_bytes(hashlib.sha256(f'A1-background|{pas}|{sp}|{model}|20261008'.encode()).digest()[:8], 'big')
    rng = np.random.default_rng(seed)
    draws = np.zeros(B, int)
    audit = []
    layers = f.groupby('stratum', sort=True) if model != 'unmatched_signed_labels' else [('all', f)]
    for key, a in layers:
        conditional = model == 'conditional_PS_reference_capacity'
        candidate_sign = a.recurrent_sign.where(a.PM_only_alignment_sign.eq(a.recurrent_sign), 0) if conditional else a.recurrent_sign
        outcome = a.PS_only_alignment_sign if conditional else a.alignment_sign
        n = len(a); k = int(candidate_sign.ne(0).sum()); success = int(outcome.ne(0).sum())
        kpos = int(candidate_sign.eq(1).sum()); kneg = int(candidate_sign.eq(-1).sum())
        npos = int(outcome.eq(1).sum()); nneg = int(outcome.eq(-1).sum()); nzero = n - npos - nneg
        pmf = signed_pmf(npos, nneg, nzero, kpos, kneg)
        assert abs(pmf.sum() - 1) < 1e-10
        distribution = np.convolve(distribution, pmf)
        pos = rng.multivariate_hypergeometric([npos, nneg, nzero], kpos, size=B)
        if n == kpos:
            negative = np.zeros(B, int)
        else:
            negative = rng.hypergeometric(nneg - pos[:, 1], n - kpos - nneg + pos[:, 1], kneg, B)
        draws += pos[:, 0] + negative
        audit.append(dict(pas_source=pas, species=sp, stratum=key, genes=n, candidates=k,
                          aligned=success, positive_candidates=kpos, negative_candidates=kneg,
                          positive_aligned=npos, negative_aligned=nneg, candidate_aligned=int(a.is_dual_reference.sum()),
                          membership_movable=(sum(v > 0 for v in [kpos, kneg, n - k]) > 1),
                          outcome_movable=(sum(v > 0 for v in [npos, nneg, nzero]) > 1),
                          informative=(float((np.arange(len(pmf)) - np.arange(len(pmf)) @ pmf) ** 2 @ pmf) > 1e-12)))
    distribution /= distribution.sum()
    value = np.arange(len(distribution))
    expected = float(value @ distribution)
    variance = float((value - expected) ** 2 @ distribution)
    assert abs(draws.mean() - expected) <= max(.1, 5 * np.sqrt(variance / B))
    qs = [int(np.searchsorted(np.cumsum(distribution), q)) for q in [.025, .5, .975]]
    p = float(distribution[observed:].sum())
    diagnostic = pd.DataFrame(audit)
    np.savez_compressed(HERE / f'null-{pas}-{sp}-{model}.npz', exact_probability=distribution,
                        draws=draws, observed=observed, seed=np.array(str(seed)))
    row = dict(test_id=f'{pas}|{sp}|{model}', family='A1-background', pas_source=pas,
               species=sp, test=model, method='Exact signed-label multivariate hypergeometric convolution',
               p_value=p, observed=observed, candidates=int(f.is_pm_recurrent.sum()), background=len(f),
               null_mean=expected, null_p025=qs[0], null_median=qs[1], null_p975=qs[2],
               mc_p=(1 + (draws >= observed).sum()) / (B + 1), seed=str(seed), draws=B,
               strata=len(audit), membership_movable_genes=int(diagnostic.loc[diagnostic.membership_movable, 'genes'].sum()),
               informative_genes=int(diagnostic.loc[diagnostic.informative, 'genes'].sum()),
               interpretation='post-selection association; not independent evidence')
    diagnostic['model'] = model
    return row, diagnostic


def a1(f, r, pas):
    tests = []; counts = []; strata = []
    for sp, a in f.groupby('species'):
        cases = a[a.is_pm_recurrent]; controls = a[~a.is_pm_recurrent]
        tests.append(fisher_row('unmatched_all_measurable_any_direction', 'A1-background', pas, sp,
                               int(cases.aligned_any_direction.sum()), int((~cases.aligned_any_direction).sum()),
                               int(controls.aligned_any_direction.sum()), int((~controls.aligned_any_direction).sum()),
                               interpretation='post-selection association; not independent evidence'))
        for model in ['unmatched_signed_labels', 'reference_capacity', 'conditional_PS_reference_capacity']:
            row, diagnostics = matched_null(a, pas, sp, model); tests.append(row); strata.append(diagnostics)
        dual = cases[cases.is_dual_reference]
        for scope, comparison in [('full_pm_recurrence', r[r.species.eq(sp)]),
                                  ('measurable_pm_recurrence', cases)]:
            remainder = comparison[~comparison.gene_id.isin(dual.gene_id)]
            tests.append(fisher_row('proximal_' + scope, 'A1-background', pas, sp,
                                   int(dual.recurrent_sign.eq(-1).sum()), int(dual.recurrent_sign.eq(1).sum()),
                                   int(remainder.recurrent_sign.eq(-1).sum()), int(remainder.recurrent_sign.eq(1).sum()),
                                   interpretation='extra direction enrichment beyond inherited PM selection'))
        for name, q, flag in [('all_measurable', a, 'aligned_any_direction'), ('pm_recurrent_measurable', cases, 'is_dual_reference'),
                              ('non_pm_recurrent_measurable', controls, 'aligned_any_direction')]:
            counts.append(dict(pas_source=pas, species=sp, scope=name, measurable=len(q), aligned=int(q[flag].sum()),
                               fraction=float(q[flag].mean()), proximal=int((q[flag] & q.alignment_sign.eq(-1)).sum()),
                               distal=int((q[flag] & q.alignment_sign.eq(1)).sum()),
                               aligned_any_direction=int(q.aligned_any_direction.sum()),
                               opposite_recurrent_alignment=int(q.opposite_recurrent_alignment.sum())))
    return tests, counts, strata


def weights(frame):
    levels = [['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'],
              ['time_a', 'sample_unit_a', 'time_b'], ['time_a', 'sample_unit_a'], ['time_a'], []]
    q = frame.reset_index(drop=True).copy()
    result = 1 / q.groupby(levels[0], dropna=False).time_a.transform('size').to_numpy(float)
    for child, parent in zip(levels[:-1], levels[1:]):
        distinct = q[child].drop_duplicates()
        if not parent:
            result /= len(distinct)
        else:
            n = distinct.groupby(parent, dropna=False).size().rename('_n').reset_index()
            result /= q[parent].merge(n, on=parent, how='left', validate='many_to_one')._n.to_numpy()
    assert abs(result.sum() - 1) < 1e-12
    return result


def a7(base, pas, f, m):
    states = read(base, 'intermediate/state_means.tsv.gz')
    states = states[states.species.eq('mouse') & states.setting.eq('supported50_p08') & states.state.eq('E')]
    gates = read(base, 'intermediate/global_contrast_eligibility.tsv.gz')
    gates = gates[gates.species.eq('mouse') & gates.setting.eq('supported50_p08') & gates.contrast.eq('P_E') & gates.n_joint_pas.ge(2)]
    e = states.merge(gates[['route', 'gene_id', 'n_joint_pas']], on=['route', 'gene_id'], validate='one_to_one')
    e = e.rename(columns={'value': 'E_DWUI'})
    q = m[m.species.eq('mouse') & m.route.isin(ROUTES['mouse'])].merge(
        e[['route', 'gene_id', 'E_DWUI', 'n_dates', 'n_units', 'n_joint_pas']], on=['route', 'gene_id'], validate='one_to_one')
    roster = f[f.species.eq('mouse')][['gene_id', 'is_pm_recurrent', 'is_dual_reference', 'recurrent_sign', 'alignment_sign']]
    q = q.merge(roster, on='gene_id', validate='many_to_one')
    q['DWUI_PE'] = q.P_DWUI - q.E_DWUI
    q['E_support_recurrent_direction'] = q.recurrent_sign * q.DWUI_PE >= TAU
    q['triple_alignment_sign'] = np.where((q[['DWUI_PS', 'DWUI_PM', 'DWUI_PE']] >= TAU).all(1), 1,
                                         np.where((q[['DWUI_PS', 'DWUI_PM', 'DWUI_PE']] <= -TAU).all(1), -1, 0))
    q['pas_source'] = pas
    tests = []; summaries = []
    for route, a in q.groupby('route'):
        c = a[a.is_pm_recurrent]; selected = c[c.is_dual_reference]; rest = c[~c.is_dual_reference]
        tests.append(fisher_row(route + '_extra_E_support', 'A7-third-reference', pas, 'mouse',
                               int(selected.E_support_recurrent_direction.sum()), int((~selected.E_support_recurrent_direction).sum()),
                               int(rest.E_support_recurrent_direction.sum()), int((~rest.E_support_recurrent_direction).sum()),
                               route=route, interpretation='same-route stage reference, not same-library or cycle-matched'))
        summaries.append(dict(pas_source=pas, design='same_route_stage_E', route=route, measurable_background=len(a),
                              pm_candidates=len(c), dual_measurable=len(selected),
                              dual_retained=int(selected.E_support_recurrent_direction.sum()),
                              other_candidates=len(rest), other_E_supported=int(rest.E_support_recurrent_direction.sum())))
    available = q.groupby('gene_id').route.agg(set)
    complete = q[q.gene_id.isin(available[available.map(lambda a: a == set(ROUTES['mouse']))].index)]
    paired = complete.groupby('gene_id').agg(is_dual_reference=('is_dual_reference', 'first'),
                    is_pm_recurrent=('is_pm_recurrent', 'first'), recurrent_sign=('recurrent_sign', 'first'),
                    alignment_min=('triple_alignment_sign', 'min'), alignment_max=('triple_alignment_sign', 'max'),
                    E_support=('E_support_recurrent_direction', 'all'))
    summaries.append(dict(pas_source=pas, design='same_route_stage_E', route='CIP;WOT', measurable_background=len(paired),
                          pm_candidates=int(paired.is_pm_recurrent.sum()), dual_measurable=int(paired.is_dual_reference.sum()),
                          dual_retained=int((paired.is_dual_reference & paired.E_support).sum()),
                          other_candidates=int((paired.is_pm_recurrent & ~paired.is_dual_reference).sum()),
                          other_E_supported=int((paired.is_pm_recurrent & ~paired.is_dual_reference & paired.E_support).sum())))
    # Extended same-library quadruplets: join PE support into the frozen target P and source S design.
    t = read(base, 'matched_pms/weighted_triplets.tsv.gz')
    t = t[t.setting.eq('supported50_p08') & t.species.eq('mouse')]
    ep = read(base, 'intermediate/paired_library_effects.tsv.gz')
    ep = ep[ep.species.eq('mouse') & ep.setting.eq('supported50_p08') & ep.contrast.eq('P_E')]
    ep = ep[['route', 'gene_id', 'library_key', 'delta']].rename(columns={'library_key': 'library_key_a', 'delta': 'DWUI_PE'})
    pm = read(base, 'state_vectors/pair_effects.primary.tsv.gz', usecols=['route', 'gene_id', 'library_key_a', 'contrast', 'delta_DWUI'])
    pm = pm[pm.contrast.eq('P_M')][['route', 'gene_id', 'library_key_a', 'delta_DWUI']].rename(columns={'delta_DWUI': 'DWUI_PM'})
    assert not pm.duplicated(['route', 'gene_id', 'library_key_a']).any()
    t = t.merge(ep, on=['route', 'gene_id', 'library_key_a'], validate='many_to_one').merge(
        pm, on=['route', 'gene_id', 'library_key_a'], validate='many_to_one')
    strict_rows = []
    for (route, gene), a in t.groupby(['route', 'gene_id']):
        w = weights(a)
        strict_rows.append(dict(pas_source=pas, route=route, gene_id=gene, gene_symbol=a.gene_symbol.iloc[0],
                                n_triplets=len(a), n_target_libraries=a.library_key_a.nunique(),
                                DWUI_PS=float(w @ a.delta_DWUI), DWUI_PM=float(w @ a.DWUI_PM), DWUI_PE=float(w @ a.DWUI_PE)))
    strict = pd.DataFrame(strict_rows).merge(roster, on='gene_id', validate='many_to_one')
    strict['triple_alignment_sign'] = np.where((strict[['DWUI_PS', 'DWUI_PM', 'DWUI_PE']] >= TAU).all(1), 1,
                                              np.where((strict[['DWUI_PS', 'DWUI_PM', 'DWUI_PE']] <= -TAU).all(1), -1, 0))
    strict['retained_original_dual'] = strict.is_dual_reference & strict.triple_alignment_sign.eq(strict.recurrent_sign)
    for route in ROUTES['mouse']:
        a = strict[strict.route.eq(route)]
        summaries.append(dict(pas_source=pas, design='same_library_P_S_M_E', route=route, measurable_background=len(a),
                              pm_candidates=int(a.is_pm_recurrent.sum()), dual_measurable=int(a.is_dual_reference.sum()),
                              dual_retained=int(a.retained_original_dual.sum()), status='available' if len(a) else 'no_same_library_P_E_in_frozen_endpoint_design'))
    return tests, summaries, q, strict


def main():
    started = datetime.now(timezone.utc).isoformat()
    all_tests = []; counts = []; backgrounds = []; strata = []; es = []; stricts = []; a7_counts = []
    for pas, base in BASES.items():
        print('Reading', pas, flush=True)
        f, m, r = background(pas, base)
        tests, cnt, st = a1(f, r, pas)
        all_tests.extend(tests); counts.extend(cnt); strata.extend(st); backgrounds.append(f)
        print('A1 finished', pas, flush=True)
        tests, cnt, e, strict = a7(base, pas, f, m)
        all_tests.extend(tests); a7_counts.extend(cnt); es.append(e); stricts.append(strict)
        print('A7 finished', pas, flush=True)
    save(pd.DataFrame(all_tests), 'tests.tsv')
    save(pd.DataFrame(counts), 'background_summary.tsv')
    save(pd.concat(backgrounds, ignore_index=True), 'gene_background.tsv.gz')
    save(pd.concat(strata, ignore_index=True), 'strata.tsv.gz')
    save(pd.DataFrame(a7_counts), 'E_summary.tsv')
    save(pd.concat(es, ignore_index=True), 'E_route_gene_effects.tsv.gz')
    save(pd.concat(stricts, ignore_index=True), 'E_same_library_gene_effects.tsv.gz')
    (HERE / 'provenance.json').write_text(json.dumps(dict(started=started, finished=datetime.now(timezone.utc).isoformat(),
        status='PASS', script_sha256=digest(Path(__file__)), protocol_sha256=digest(HERE / 'PROTOCOL.md'),
        python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        input_files=INPUTS, tau=TAU, monte_carlo_draws=B, source_roles={'polyapipe': 'primary', 'apadetect-p50': 'validation'}), indent=2))
    print(pd.DataFrame(counts).to_string(index=False), flush=True)
    print(pd.DataFrame(a7_counts).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
