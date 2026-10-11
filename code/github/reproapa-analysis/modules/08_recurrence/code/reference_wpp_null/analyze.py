from pathlib import Path
import hashlib
import importlib.util
import json
import time

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent
M = P.parents[1]
O = M / 'outputs' / P.name
O.mkdir(parents=True, exist_ok=True)
OLD = M / 'outputs/figure3-study-recurrence-null-2026-09-20'
SOURCE = M / 'outputs/figure3-apa-gex-state-analysis-2026-09-19'
RAW = M / 'analysis/figure3-apa-gex-state-analysis-2026-09-19/raw'
original_path = M / 'analysis/figure3-study-recurrence-null-2026-09-20/analyze.py'
spec = importlib.util.spec_from_file_location('original_recurrence', original_path)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
N = 10_000
SEED = 20260923
MODELS = ['original', 'baseline_q5_capacity', 'baseline_q3_capacity', 'baseline_range10_capacity']
source_files = [original_path]


def read(path, **kwargs):
    source_files.append(path)
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, filename):
    frame.to_csv(O / filename, sep='\t', index=False)


def reconstruct():
    pairs = read(SOURCE / 'matched_library_pair_effects.tsv.gz')
    pairs = pairs[pairs.setting.eq('supported50_p08') & (
        (pairs.contrast.eq('P_S') & pairs.design.eq('true_source')) |
        (pairs.contrast.eq('P_M') & pairs.design.eq('same_library')))].copy()
    state_cols = ['species', 'study_id', 'route', 'library_key', 'time', 'setting', 'state', 'gene_id', 'APA']
    states = []
    for file in sorted(RAW.glob('gse*/apa_gex_states.tsv.gz')):
        d = read(file, usecols=state_cols)
        states.append(d[d.setting.eq('supported50_p08') & d.state.isin(['P', 'S', 'M'])])
    states = pd.concat(states, ignore_index=True)
    join = ['species', 'study_id', 'route', 'setting', 'gene_id', 'library_key', 'time', 'state']
    assert not states.duplicated(join).any()
    for suffix, name in [('a', 'target_WPP'), ('b', 'reference_WPP')]:
        pairs[f'state_{suffix}'] = 'P' if suffix == 'a' else pairs.contrast.str[-1]
        rename = {x: x + '_' + suffix for x in ['library_key', 'time', 'state']}
        rename['APA'] = name
        on = [rename.get(x, x) for x in join]
        pairs = pairs.merge(states.rename(columns=rename), on=on, how='left', validate='many_to_one')
    assert pairs[['target_WPP', 'reference_WPP']].notna().all().all()
    assert pairs[['target_WPP', 'reference_WPP']].ge(-1e-8).all().all()
    assert pairs[['target_WPP', 'reference_WPP']].le(1 + 1e-8).all().all()
    pair_error = (pairs.target_WPP - pairs.reference_WPP - pairs.APA_delta).abs().max()
    assert pair_error < 1e-8
    keys = ['species', 'study_id', 'route', 'setting', 'gene_id', 'gene_symbol', 'contrast', 'design']
    vals = ['target_WPP', 'reference_WPP', 'APA_delta']
    u = pairs.groupby(keys + ['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'], dropna=False)[vals].mean().reset_index()
    u = u.groupby(keys + ['time_a', 'sample_unit_a', 'time_b'])[vals].mean().groupby(keys + ['time_a', 'sample_unit_a'])[vals].mean().reset_index()
    dates = u.groupby(keys + ['time_a'])[vals].mean().reset_index()
    routes = dates.groupby(keys)[vals].mean().reset_index()
    reported = read(SOURCE / 'matched_study_effects.tsv.gz')
    reported = reported[reported.setting.eq('supported50_p08') & (
        (reported.contrast.eq('P_S') & reported.design.eq('true_source')) |
        (reported.contrast.eq('P_M') & reported.design.eq('same_library')))]
    q = routes.merge(reported[keys + ['APA_delta']], on=keys, suffixes=('', '_reported'), validate='one_to_one')
    assert len(q) == len(routes) == len(reported)
    q['delta_reconstruction_error'] = (q.target_WPP - q.reference_WPP - q.APA_delta_reported).abs()
    assert q.delta_reconstruction_error.max() < 1e-8
    save(q, 'route_reference_wpp.tsv')
    save(dates, 'date_reference_wpp.tsv.gz')
    (O / 'baseline_reconstruction.json').write_text(json.dumps({
        'qualified_library_pairs': len(pairs), 'route_effects': len(q),
        'max_pair_error': float(pair_error), 'max_route_error': float(q.delta_reconstruction_error.max()),
        'reference_missing': 0,
    }, indent=2))
    print('Baselines reconstructed:', len(q), 'route effects', flush=True)
    return q


def stratify(frame, model):
    f = old.add_strata(frame, include_pau=False)
    f['original_stratum'] = f.stratum
    if model == 'original':
        return f
    threshold = f.threshold.iloc[0]
    f['can_negative'] = f.reference_WPP_min.ge(threshold - 1e-9)
    f['can_positive'] = f.reference_WPP_max.le(1 - threshold + 1e-9)
    assert f.loc[f.event.eq(-1), 'can_negative'].all()
    assert f.loc[f.event.eq(1), 'can_positive'].all()
    if model.startswith('baseline_q'):
        nbin = 5 if 'q5' in model else 3
        f['baseline_bin'] = f.groupby('unit').reference_WPP.transform(lambda x: old.qbin(x, q=nbin))
        parts = ['original_stratum', 'baseline_bin', 'can_negative', 'can_positive']
    else:
        for c in ['reference_WPP_min', 'reference_WPP_max']:
            f[c + '_bin'] = np.minimum(np.floor(f[c].clip(0, 1) * 10).astype(int), 9)
        parts = ['original_stratum', 'reference_WPP_min_bin', 'reference_WPP_max_bin', 'can_negative', 'can_positive']
    f['stratum'] = f[parts].astype(str).agg('|'.join, axis=1)
    return f


def diagnostics(frame):
    rows = []
    for (unit, stratum), f in frame.groupby(['unit', 'stratum']):
        n = len(f)
        npos, nneg, nzero = [(f.event == v).sum() for v in [1, -1, 0]]
        probs = np.array([npos, nneg, nzero]) / n
        rows.append(dict(unit=unit, stratum=stratum, n_genes=n, positive=int(npos), negative=int(nneg), non_event=int(nzero),
                         mixed_labels=f.event.nunique() > 1,
                         expected_labels_changed=n * (1 - np.sum(probs**2)),
                         baseline_min=f.reference_WPP.min(), baseline_max=f.reference_WPP.max(),
                         baseline_range=f.reference_WPP.max()-f.reference_WPP.min(),
                         expected_squared_baseline_difference=2*np.var(f.reference_WPP, ddof=0)))
    return pd.DataFrame(rows)


def checked_draws(frame, rng):
    # Meaningful invariants: every draw preserves each study/stratum's label counts
    # and every gene's eligibility; labels cannot become geometrically impossible.
    for _ in range(5):
        for _, f in frame.groupby(['unit', 'stratum']):
            label = rng.permutation(f.event.to_numpy())
            assert np.array_equal(np.sort(label), np.sort(f.event))
            if 'can_positive' in f:
                assert not ((label == 1) & ~f.can_positive.to_numpy()).any()
                assert not ((label == -1) & ~f.can_negative.to_numpy()).any()


def run():
    start = time.monotonic()
    route = reconstruct()
    baseline = route.groupby(['species', 'contrast', 'study_id', 'gene_id']).agg(
        reference_WPP=('reference_WPP', 'median'), reference_WPP_min=('reference_WPP', 'min'),
        reference_WPP_max=('reference_WPP', 'max'), baseline_routes=('route', 'nunique')).reset_index().rename(columns={'study_id': 'unit'})
    events = read(OLD / 'event_observability_matrix.tsv.gz')
    events = events[events.level.eq('study')].copy()
    n = len(events)
    events = events.merge(baseline, on=['species', 'contrast', 'unit', 'gene_id'], how='left', validate='many_to_one')
    assert len(events) == n and events.reference_WPP.notna().all()
    assert events.n_routes.eq(events.baseline_routes).all()
    # Recompute frozen event labels directly from route effects, including series unanimity.
    for threshold, f in events.groupby('threshold'):
        direct = route.groupby(['species', 'contrast', 'study_id', 'gene_id']).APA_delta.agg(['min', 'max']).reset_index().rename(columns={'study_id': 'unit'})
        direct['direct_event'] = np.where(direct['min'] >= threshold, 1, np.where(direct['max'] <= -threshold, -1, 0))
        z = f.merge(direct, on=['species', 'contrast', 'unit', 'gene_id'], validate='one_to_one')
        assert z.event.eq(z.direct_event).all()
    save(events, 'event_baseline_matrix.tsv.gz')
    sums, groupstats, diags, strata, archive, sets = [], [], [], [], {}, []
    original_summary = read(OLD / 'recurrence_null_summary.tsv')
    original_summary = original_summary[original_summary.level.eq('study') & original_summary.matching.eq('main')]
    for (species, contrast, threshold), f in events.groupby(['species', 'contrast', 'threshold']):
        runs = [('none', f)]
        if f.unit.nunique() >= 3 and threshold == .02:
            for unit in sorted(f.unit.unique()):
                reduced = f[f.unit.ne(unit)].copy()
                reduced['eligible_units'] = reduced.groupby('gene_id').unit.transform('nunique')
                runs.append((unit, reduced))
        for omitted, base in runs:
            models = MODELS if omitted == 'none' else ['original', 'baseline_q5_capacity']
            for model in models:
                prefix = dict(species=species, contrast=contrast, threshold=threshold, model=model, omitted_study=omitted,
                              independent_studies=base.unit.nunique(), permutations=N)
                key = f'{species}__{contrast}__{threshold:.2f}__{omitted}__{model}'
                seed = SEED + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
                rng = np.random.default_rng(seed)
                frame = stratify(base, model)
                checked_draws(frame, np.random.default_rng(seed + 1))
                summary, recurrent, diag, null = old.permute_analysis(frame, rng)
                # Invariance of observed counts is checked against frozen results, not significance.
                if omitted == 'none':
                    prev = original_summary[(original_summary.species == species) & (original_summary.contrast == contrast) & (original_summary.threshold == threshold)]
                    z = summary.merge(prev[['metric','observed']], on='metric', suffixes=('', '_previous'))
                    assert len(z) == len(summary) and z.observed.eq(z.observed_previous).all()
                summary['enrichment_vs_mean'] = summary.observed / summary.null_mean.replace(0, np.nan)
                sums.append(summary.assign(seed=seed, **prefix))
                groupstats.append(diagnostics(frame).assign(**prefix))
                diags.append(diag.assign(**prefix))
                strata.append(frame.assign(model=model, omitted_study=omitted))
                sets.append(recurrent.assign(**prefix))
                for metric, values in null.items():archive[key + '__' + metric] = values
                r = summary[summary.metric.eq('N_ge2')].iloc[0]
                print(key, 'obs', r.observed, 'null', round(r.null_mean, 3), 'P', r.empirical_p,
                      'elapsed', round(time.monotonic()-start), flush=True)
    summary = pd.concat(sums, ignore_index=True)
    summary['q_eight_primary_settings'] = np.nan
    for _, f in summary[summary.omitted_study.eq('none') & summary.metric.eq('N_ge2')].groupby('model'):
        summary.loc[f.index,'q_eight_primary_settings'] = old.bh(f.empirical_p)
    summary['q_three_LOO'] = np.nan
    for _, f in summary[summary.omitted_study.ne('none') & summary.metric.eq('N_ge2')].groupby('model'):
        summary.loc[f.index,'q_three_LOO'] = old.bh(f.empirical_p)
    save(summary, 'recurrence_null_summary.tsv')
    save(pd.concat(sets, ignore_index=True), 'recurrent_gene_sets.tsv.gz')
    save(pd.concat(diags, ignore_index=True), 'unit_diagnostics.tsv')
    save(pd.concat(groupstats, ignore_index=True), 'stratum_diagnostics.tsv.gz')
    save(pd.concat(strata, ignore_index=True), 'gene_strata.tsv.gz')
    np.savez_compressed(O/'null_draws.npz', **archive)
    primary = summary[summary.omitted_study.eq('none') & summary.metric.eq('N_ge2')]
    prior = original_summary[original_summary.metric.eq('N_ge2')][['species','contrast','threshold','observed','null_mean','null_median','enrichment_vs_median','empirical_p']]
    comparison = primary.merge(prior, on=['species','contrast','threshold','observed'], suffixes=('','_frozen_original'), validate='many_to_one')
    save(comparison, 'original_vs_baseline.tsv')
    (O/'input_manifest.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in set(source_files)}, indent=2))
    (O/'run_config.json').write_text(json.dumps(dict(seed=SEED, permutations=N, models=MODELS, runtime_seconds=time.monotonic()-start,
        original_analysis=str(original_path), primary_model='baseline_q5_capacity', status='computed'), indent=2))


if __name__ == '__main__':
    run()
