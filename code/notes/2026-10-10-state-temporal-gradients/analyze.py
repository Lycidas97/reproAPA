"""Fixed-gene, actual-date curves and descriptive adjacent-state displacement."""
from pathlib import Path
import hashlib
import itertools
import json
import re

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
SETTING = 'supported50_p08'
SUBSETS = ['all_fixed', 'full_primary_tandem_UTR', 'P_S_recurrent_proximal_002',
           'P_S_recurrent_distal_002', 'P_M_recurrent_proximal_002', 'P_M_recurrent_distal_002']
ROUTES = ['GSE118258', '178325_0618', '178325_0809', 'GSE221739', 'CIP', 'WOT']
STAGE_ROUTES = ['GSE118258', 'CIP', 'WOT']
INPUTS = {}
CHECKS = {}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(name):
    path = TABLES / name
    INPUTS[str(path.relative_to(ROOT))] = sha(path)
    return pd.read_csv(path, sep='\t')


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.12g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def truth(values):
    return values.astype(str).str.lower().eq('true')


def fixed_library_genes(frame, libraries):
    """Missing one planned physical library makes a gene unavailable."""
    libraries = set(libraries)
    if not libraries:
        return set()
    frame = frame[frame.library_key.isin(libraries)]
    return set(frame.groupby('gene_id').library_key.nunique().loc[lambda x: x.eq(len(libraries))].index)


def daily_values(frame):
    """Library means within sampling unit, units within date; no cell pooling."""
    units = frame.groupby(['gene_id', 'time', 'sample_unit']).value.mean().reset_index()
    return units.groupby(['gene_id', 'time']).value.mean().reset_index()


def slope(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xc = x - x.mean()
    return float(xc @ y / (xc @ xc))


def date_permutation(x, y, seed, draws=10000, exact=False):
    x, y = np.asarray(x, float), np.asarray(y, float)
    xc = x - x.mean()
    denominator = xc @ xc
    observed = slope(x, y)
    if exact:
        simulated = np.array([np.asarray(p) @ xc / denominator for p in itertools.permutations(y)])
        p = float((np.abs(simulated) >= abs(observed) - 1e-14).mean())
    else:
        rng = np.random.default_rng(seed)
        shuffled = y[np.argsort(rng.random((draws, len(y))), axis=1)]
        simulated = shuffled @ xc / denominator
        p = (1 + int((np.abs(simulated) >= abs(observed) - 1e-14).sum())) / (draws + 1)
    return observed, simulated, p


def bh(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    out = np.empty(len(values))
    out[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return out


def numeric_day(value):
    match = re.fullmatch(r'D(\d+(?:\.\d+)?)', str(value))
    return float(match.group(1)) if match else np.nan


def summarize_values(values):
    values = np.asarray(values, float)
    return dict(mean=float(values.mean()) if len(values) else np.nan,
                median=float(np.median(values)) if len(values) else np.nan,
                proximal_002=int((values <= -.02).sum()), distal_002=int((values >= .02).sum()))


def make_subsets(ids, species, study, tandem, members):
    results = {'all_fixed': set(ids)}
    results['full_primary_tandem_UTR'] = set(ids) & set(tandem.loc[tandem.study_id.eq(study) & truth(tandem.tandem_utr_compatible), 'gene_id'])
    for contrast, direction in itertools.product(['P_S', 'P_M'], ['negative', 'positive']):
        name = contrast + '_recurrent_' + ('proximal' if direction == 'negative' else 'distal') + '_002'
        results[name] = set(ids) & set(members.loc[members.species.eq(species) & members.contrast.eq(contrast) & members.direction.eq(direction), 'gene_id'])
    return results


def curves(states, pairs, inventory, tandem, members):
    summaries, libs_out, dates_out, endpoints, universe, tests, inventory_out, sensitivity = [], [], [], [], [], [], [], []
    component_dates, component_summaries = [], []
    nulls = {}
    for route in ROUTES:
        meta = inventory[inventory.route.eq(route)]
        state = states[states.route.eq(route)]
        species = state.species.iloc[0]
        study = state.study_id.iloc[0]
        p_libs = set(meta.loc[meta.state.eq('P'), 'library_key'])
        m_libs = set(meta.loc[meta.state.eq('M'), 'library_key'])
        for metric in ['P_DWUI', 'P_M']:
            required = p_libs if metric == 'P_DWUI' else p_libs & m_libs
            frame = state[state.state.eq('P')].copy() if metric == 'P_DWUI' else pairs[pairs.route.eq(route) & pairs.contrast.eq('P_M')].rename(columns={'delta': 'value'}).copy()
            frame = frame[frame.library_key.isin(required)]
            assert not frame.duplicated(['gene_id', 'library_key']).any()
            fixed = fixed_library_genes(frame, required)
            meta_columns = ['library_key', 'sample_unit', 'time', 'replicate_status']
            library_meta = frame[meta_columns].drop_duplicates()
            if required:
                assert set(library_meta.library_key) == required, (route, metric, 'a planned library has no measurable genes')
            for row in library_meta.to_dict('records'):
                inventory_out.append(dict(species=species, study_id=study, route=route, metric=metric, **row))
            dates = sorted(frame.time.unique(), key=lambda t: (np.isnan(numeric_day(t)), numeric_day(t), str(t)))
            is_numeric = bool(dates) and all(np.isfinite(numeric_day(t)) for t in dates)
            status = 'measured' if len(dates) >= 2 and is_numeric else ('insufficient_P_dates' if required else 'no_same_library_P_M')
            subsets = make_subsets(fixed, species, study, tandem, members)
            for subset, ids in subsets.items():
                context = dict(species=species, study_id=study, route=route, metric=metric, subset=subset)
                z = frame[frame.gene_id.isin(ids)]
                universe.extend(context | dict(gene_id=gid) for gid in sorted(ids))
                record = context | dict(status=status if ids or not required else 'empty_fixed_subset', n_genes=len(ids),
                                        n_libraries=len(required), n_dates=len(dates), first_date=dates[0] if dates else '', last_date=dates[-1] if dates else '')
                if not ids:
                    summaries.append(record)
                    continue
                assert len(z) == len(ids) * len(required)
                span = tandem[tandem.study_id.eq(study)].set_index('gene_id').native_axis_span_nt
                assert pd.Index(sorted(ids)).isin(span.index).all()
                for key, part in z.groupby(meta_columns, dropna=False):
                    libs_out.append(context | dict(zip(meta_columns, key), n_genes=len(part), value=part.value.mean(),
                                                  retained_coordinate_nt=(part.value * part.gene_id.map(span)).mean() if subset == 'full_primary_tandem_UTR' else np.nan))
                daily = daily_values(z)
                assert daily.groupby('time').gene_id.nunique().eq(len(ids)).all()
                for time, part in daily.groupby('time'):
                    n_units = z[z.time.eq(time)].sample_unit.nunique()
                    dates_out.append(context | dict(time=time, day=numeric_day(time), n_genes=len(part), n_units=n_units, value=part.value.mean(),
                                                   retained_coordinate_nt=(part.value * part.gene_id.map(span)).mean() if subset == 'full_primary_tandem_UTR' else np.nan))
                wide = daily.pivot(index='gene_id', columns='time', values='value').reindex(columns=dates)
                assert wide.notna().all().all()
                if status != 'measured':
                    summaries.append(record)
                    continue
                if metric == 'P_M':
                    contributor_keys = ['gene_id', 'library_key', 'time', 'sample_unit']
                    native = state[state.state.isin(['P', 'M'])].pivot(index=contributor_keys, columns='state', values='value')
                    joined = z[contributor_keys + ['value']].merge(native.reset_index(), on=contributor_keys, validate='one_to_one')
                    assert joined[['P', 'M']].notna().all().all()
                    assert np.allclose(joined.value, joined.P - joined.M, rtol=0, atol=1e-12)
                    change = {}
                    for node in ['P', 'M']:
                        node_daily = daily_values(joined[contributor_keys + [node]].rename(columns={node: 'value'}))
                        means = node_daily.groupby('time').value.mean().reindex(dates)
                        change[node] = float(means.iloc[-1] - means.iloc[0])
                        for time, value in means.items():
                            component_dates.append(context | dict(state=node, time=time, day=numeric_day(time), n_genes=len(ids), value=value))
                    delta_PM = float(wide[dates[-1]].mean() - wide[dates[0]].mean())
                    error = delta_PM - change['P'] + change['M']
                    assert abs(error) < 1e-12
                    component_summaries.append(context | dict(n_genes=len(ids), n_dates=len(dates), n_libraries=len(required),
                        first_date=dates[0], last_date=dates[-1], P_last_minus_first=change['P'],
                        M_last_minus_first=change['M'], PM_last_minus_first=delta_PM, identity_error=error))
                diff = wide[dates[-1]] - wide[dates[0]]
                for gene, delta in diff.items():
                    endpoints.append(context | dict(gene_id=gene, first_date=dates[0], last_date=dates[-1], first_value=wide.at[gene, dates[0]],
                                                   last_value=wide.at[gene, dates[-1]], delta_DWUI=delta,
                                                   delta_retained_nt=delta * span[gene] if subset == 'full_primary_tandem_UTR' else np.nan))
                record.update({('last_minus_first_' + k): v for k, v in summarize_values(diff).items()})
                if subset == 'full_primary_tandem_UTR':
                    record['last_minus_first_mean_nt'] = float((diff * span.loc[diff.index]).mean())
                    record['last_minus_first_median_nt'] = float((diff * span.loc[diff.index]).median())
                if status == 'measured':
                    x = np.array([numeric_day(t) for t in dates])
                    y = wide.mean(0).to_numpy()
                    record.update(slope_DWUI_per_day=slope(x, y), pearson_r=float(np.corrcoef(x, y)[0, 1]),
                                  spearman_rho=float(np.corrcoef(pd.Series(x).rank(), pd.Series(y).rank())[0, 1]))
                    if len(dates) >= 3:
                        for omit in dates:
                            keep = [i for i, d in enumerate(dates) if d != omit]
                            sensitivity.append(context | dict(sensitivity='leave_one_date', omitted_date=omit,
                                n_genes=len(ids), n_dates=len(keep), slope_DWUI_per_day=slope(x[keep], y[keep]),
                                last_minus_first_mean=float(y[keep][-1] - y[keep][0])))
                    if route in ['WOT', 'GSE221739']:
                        balanced = [i for i, time in enumerate(dates) if z[z.time.eq(time)].sample_unit.nunique() >= 2]
                        if len(balanced) >= 2:
                            sensitivity.append(context | dict(sensitivity='at_least_two_units_per_date', omitted_date=';'.join(d for i, d in enumerate(dates) if i not in balanced),
                                n_genes=len(ids), n_dates=len(balanced), slope_DWUI_per_day=slope(x[balanced], y[balanced]),
                                last_minus_first_mean=float(y[balanced][-1] - y[balanced][0])))
                    if route == 'WOT' and subset == 'all_fixed':
                        assert library_meta.replicate_status.eq('biological_independent').all()
                        seed = 2026101021 + (0 if metric == 'P_DWUI' else 1)
                        observed, draws, p = date_permutation(x, y, seed)
                        tests.append(context | dict(n_genes=len(ids), n_dates=len(dates), n_libraries=len(required), statistic='date_mean_slope',
                                                    slope=observed, p_value=p, seed=seed, null_draws=len(draws),
                                                    status='exploratory_date_exchangeability_required'))
                        nulls[metric] = draws
                summaries.append(record)
    result = pd.DataFrame(tests)
    assert set(result.metric) == {'P_DWUI', 'P_M'} and len(result) == 2
    result['q_BH2'] = bh(result.p_value)
    save(result, 'date_level_tests.tsv')
    np.savez_compressed(HERE / 'date_slope_null_draws.npz', **nulls)
    for rows, filename in [(summaries, 'curve_summary.tsv'), (libs_out, 'curve_library_values.tsv'), (dates_out, 'curve_date_values.tsv'),
                           (endpoints, 'curve_gene_endpoints.tsv.gz'), (universe, 'curve_fixed_gene_membership.tsv.gz'), (inventory_out, 'curve_library_inventory.tsv')]:
        save(pd.DataFrame(rows), filename)
    save(pd.DataFrame(sensitivity), 'curve_date_sensitivity.tsv')
    save(pd.DataFrame(component_dates), 'PM_curve_component_dates.tsv')
    save(pd.DataFrame(component_summaries), 'PM_curve_component_summary.tsv')
    CHECKS['curve_composition_and_library_coverage_fixed'] = True
    CHECKS['PM_curves_same_gene_library_state_identity_verified'] = True
    return pd.DataFrame(summaries), result


def stage_profiles(states, tandem, members):
    daily = states.groupby(['species', 'study_id', 'route', 'gene_id', 'gene_symbol', 'state', 'time', 'sample_unit']).value.mean().reset_index()
    daily = daily.groupby(['species', 'study_id', 'route', 'gene_id', 'gene_symbol', 'state', 'time']).value.mean().reset_index()
    profiles, fixed_by_route = {}, {}
    for route in STAGE_ROUTES:
        q = daily[daily.route.eq(route) & daily.state.isin(list('SREP'))]
        required = q.groupby('state').time.nunique().reindex(list('SREP'))
        assert required.notna().all()
        complete = q.groupby(['gene_id', 'state']).time.nunique().unstack(fill_value=0).reindex(columns=list('SREP'), fill_value=0).eq(required, axis=1).all(1)
        ids = set(complete[complete].index)
        fixed_by_route[route] = ids
        z = q[q.gene_id.isin(ids)].groupby(['gene_id', 'state']).value.mean().unstack().reindex(columns=list('SREP'))
        assert z.notna().all().all()
        profiles[route] = z
    common_mouse = fixed_by_route['CIP'] & fixed_by_route['WOT']
    summary, genes_out, membership, coverage = [], [], [], []
    for route in STAGE_ROUTES:
        source = states[states.route.eq(route)]
        species, study = source.species.iloc[0], source.study_id.iloc[0]
        span = tandem[tandem.study_id.eq(study)].set_index('gene_id').native_axis_span_nt
        contexts = [('all_dates_within_route', fixed_by_route[route])]
        if species == 'mouse':
            contexts.append(('common_two_mouse_routes_all_dates', common_mouse))
        for background, ids in contexts:
            for subset, selected in make_subsets(ids, species, study, tandem, members).items():
                context = dict(species=species, study_id=study, route=route, background=background, subset=subset)
                z = profiles[route].loc[sorted(selected)]
                if z.empty:
                    summary.append(context | dict(n_genes=0, status='empty_fixed_subset'))
                    continue
                membership.extend(context | dict(gene_id=gene) for gene in z.index)
                steps = pd.DataFrame({'R_S': z.R - z.S, 'E_R': z.E - z.R, 'P_E': z.P - z.E}, index=z.index)
                lengths = steps.abs().sum(1)
                shares = steps.abs().div(lengths.replace(0, np.nan), axis=0)
                maximum = steps.abs().max(1)
                tied = steps.abs().sub(maximum, axis=0).abs().le(1e-14).sum(1).gt(1)
                dominant = steps.abs().idxmax(1).mask(tied, 'tie').mask(lengths.eq(0), 'zero_path')
                signed = steps.mean()
                record = context | dict(n_genes=len(z), status='descriptive_state_windows', P_S_mean=(z.P - z.S).mean(),
                                        largest_mean_proximal_step=signed.idxmin() if signed.min() < 0 else 'none_no_proximal_mean', largest_absolute_mean_step=signed.abs().idxmax(),
                                        path_length_mean=lengths.mean(), path_length_median=lengths.median())
                for step in steps:
                    record.update({step + '_' + k: v for k, v in summarize_values(steps[step]).items()})
                    record[step + '_median_absolute_path_share'] = shares[step].median()
                    record[step + '_dominant_genes'] = int(dominant.eq(step).sum())
                    if subset == 'full_primary_tandem_UTR':
                        record[step + '_mean_retained_nt'] = (steps[step] * span.loc[z.index]).mean()
                record['dominant_tie_genes'] = int(dominant.eq('tie').sum())
                record['zero_path_genes'] = int(dominant.eq('zero_path').sum())
                summary.append(record)
                for gene in z.index:
                    genes_out.append(context | dict(gene_id=gene, **z.loc[gene].to_dict(), **steps.loc[gene].to_dict(),
                                                    path_length=lengths[gene], dominant_step=dominant[gene],
                                                    **{step + '_absolute_path_share': shares.at[gene, step] for step in shares}))
                for state in 'SREP':
                    part = source[source.state.eq(state) & source.gene_id.isin(selected)]
                    coverage.append(context | dict(state=state, n_dates=part.time.nunique(), n_libraries=part.library_key.nunique(),
                                                  minimum_libraries_per_gene_date=int(part.groupby(['gene_id', 'time']).library_key.nunique().min()),
                                                  maximum_libraries_per_gene_date=int(part.groupby(['gene_id', 'time']).library_key.nunique().max())))
    for rows, filename in [(summary, 'stage_step_summary.tsv'), (genes_out, 'stage_gene_profiles.tsv.gz'),
                           (membership, 'stage_fixed_gene_membership.tsv.gz'), (coverage, 'stage_state_coverage.tsv')]:
        save(pd.DataFrame(rows), filename)
    CHECKS['stage_complete_all_state_dates'] = {route: len(ids) for route, ids in fixed_by_route.items()}
    CHECKS['common_mouse_stage_genes'] = len(common_mouse)
    return pd.DataFrame(summary)


def main():
    states = read('intermediate/library_state_values.tsv.gz')
    states = states[states.setting.eq(SETTING)]
    expected = read('intermediate/state_date_values.tsv.gz')
    expected = expected[expected.setting.eq(SETTING)]
    keys = ['species', 'study_id', 'route', 'setting', 'gene_id', 'gene_symbol', 'state', 'time']
    units = states.groupby(keys + ['sample_unit']).value.mean().reset_index()
    reproduced = units.groupby(keys).value.mean().reset_index()
    check = reproduced.merge(expected, on=keys, suffixes=('_new', '_old'), validate='one_to_one')
    assert len(check) == len(reproduced) == len(expected)
    error = float((check.value_new - check.value_old).abs().max())
    assert error < 1e-12, error
    CHECKS['primary_state_date_rows_reproduced'] = len(check)
    CHECKS['state_date_maximum_absolute_error'] = error
    pairs = read('intermediate/paired_library_effects.tsv.gz')
    pairs = pairs[pairs.setting.eq(SETTING)]
    assert pairs.n_joint_pas.ge(2).all()
    inventory = read('intermediate/state_library_inventory.tsv')
    inventory = inventory[inventory.setting.eq(SETTING)]
    observed = states[['study_id', 'route', 'library_key', 'state']].drop_duplicates()
    assert set(map(tuple, inventory[['study_id', 'route', 'library_key', 'state']].to_numpy())) == set(map(tuple, observed.to_numpy()))
    tandem = read('tandem/tandem_utr_eligibility.tsv.gz')
    assert not tandem.duplicated(['study_id', 'gene_id']).any()
    members = read('recurrence/recurrence_members.tsv.gz')
    members = members[members.metric.eq('WPP') & members.variant.eq('primary') & members.threshold.eq(.02) & truth(members.series_matching)]
    assert not members.duplicated(['species', 'contrast', 'gene_id']).any()
    CHECKS['fixed_recurrence_members'] = members.groupby(['species', 'contrast', 'direction']).size().to_dict().__str__()
    curves(states, pairs, inventory, tandem, members)
    stage_profiles(states, tandem, members)
    manifest = dict(inputs=INPUTS, protocol_sha256=sha(HERE / 'protocol.json'),
                    analysis_code_sha256=sha(HERE / 'analyze.py'),
                    upstream_code_sha256={name: sha(ROOT / name) for name in ['analysis_code/pas_rerun/downstream_intermediate.py',
                        'analysis_code/pas_rerun/downstream_runtime.py', 'analysis_code/pas_rerun/downstream_global.py']},
                    verification=CHECKS)
    (HERE / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(CHECKS, indent=2))


if __name__ == '__main__':
    main()
