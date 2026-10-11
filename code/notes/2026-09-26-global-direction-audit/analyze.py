"""Audit global WPP directions without changing frozen manuscript results."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROJECTS = Path('/Users/lycidas/work/科研/projects')
M = PROJECTS / '2026-09-02_mouse-time-aware-integration-pilot'
H = PROJECTS / '2026-09-09_human-time-aware-integration'
F1 = PROJECTS / '2026-09-10_reproapa-figure1'
F4 = PROJECTS / '2026-09-13_reproapa-figure4-real'
PROFILES = M / 'data/apa-majority-shared-v1/profiles'
STATE = M / 'outputs/figure3-apa-gex-state-analysis-2026-09-19'
RAW = M / 'analysis/figure3-apa-gex-state-analysis-2026-09-19/raw'
SOURCE_TIMES = ['D0', 'source', 'Stage_I_day_0']
KEYS = ['species', 'study_id', 'route', 'setting', 'gene_id', 'gene_symbol', 'contrast', 'design']
INPUTS = {}
CHECKS = {}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path, **kwargs):
    INPUTS[str(path)] = sha(path)
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.12g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def protected():
    paths = [p for d in ('main', 'supplement', 'Fig') for p in (ROOT / d).rglob('*') if p.is_file()]
    paths += [ROOT / p for p in ('main.tex', 'supplement.tex', 'output/pdf/main.pdf', 'output/pdf/supplement.pdf') if (ROOT / p).exists()]
    return {str(p.relative_to(ROOT)): sha(p) for p in paths}


def metrics(values):
    x = np.asarray(values, float)
    assert len(x) and np.isfinite(x).all()
    out = dict(n=len(x), mean=x.mean(), median=np.median(x), q10=np.quantile(x, .1),
               q25=np.quantile(x, .25), q75=np.quantile(x, .75), q90=np.quantile(x, .9),
               positive_contribution=np.maximum(x, 0).mean(), negative_contribution=np.minimum(x, 0).mean())
    for threshold in (.02, .05):
        suffix = str(threshold).replace('.', '')
        out['proximal_' + suffix] = int((x <= -threshold).sum())
        out['distal_' + suffix] = int((x >= threshold).sum())
        out['distal_excess_' + suffix] = ((x >= threshold).sum() - (x <= -threshold).sum()) / len(x)
    return out


def contrast_subset(x):
    return x[((x.contrast == 'P_S') & x.design.isin(['true_source', 'induced_D3_reference'])) |
             ((x.contrast == 'P_M') & (x.design == 'same_library'))].copy()


def collapse_pairs(pairs):
    """Equal source units/dates, then target units/dates; never pool cells."""
    units = pairs.groupby(KEYS + ['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'], dropna=False).APA_delta.mean()
    units = units.groupby(KEYS + ['time_a', 'sample_unit_a', 'time_b']).mean()
    units = units.groupby(KEYS + ['time_a', 'sample_unit_a']).mean()
    dates = units.groupby(KEYS + ['time_a']).mean().reset_index()
    routes = dates.groupby(KEYS).APA_delta.agg(['mean', 'size']).reset_index().rename(columns={'mean': 'delta', 'size': 'n_dates'})
    return routes, dates


def legacy():
    nodes, inventories, effects, reference_checks = [], [], [], []
    for path in sorted(PROFILES.glob('*/gene_node_positions.tsv.gz')):
        q = read(path, low_memory=False)
        meta = read(path.with_name('selected_metacells.tsv.gz'))
        assert meta.groupby('library_key').experimental_time.nunique().max() == 1
        for (sp, study, route, state, time), g in meta.groupby(['species', 'study_id', 'route', 'node', 'experimental_time']):
            inventories.append(dict(species=sp, study_id=study, route=route, state=state, time=time,
                n_mc=len(g), n_cells=int(g.n_cells.sum()), n_libraries=g.library_key.nunique(),
                libraries=';'.join(sorted(g.library_key.unique())), true_source=state == 'S' and time in SOURCE_TIMES))
        for layer in ['supported50', 'expanded_primary']:
            libs = read(path.with_name(layer + '_library_positions.tsv.gz'))
            got = libs.groupby(['route', 'node', 'gene_id']).position_index.mean().sort_index()
            expected = q[(q.layer == layer) & q.position_index.notna()].set_index(['route', 'node', 'gene_id']).position_index.sort_index()
            assert got.index.equals(expected.index) and np.allclose(got, expected, rtol=0, atol=1e-12)
            CHECKS.setdefault('library_node_records_reproduced', 0)
            CHECKS['library_node_records_reproduced'] += len(got)
            if layer == 'supported50':
                times = meta[['library_key', 'experimental_time']].drop_duplicates()
                libs = libs.merge(times, on='library_key', validate='many_to_one')
                for route_id, part in libs.groupby('route'):
                    if 'reference' in route_id or 'perturbation' in route_id:
                        continue
                    target = part[part.node.eq('P')].groupby('gene_id').position_index.mean()
                    source = part[part.node.eq('S')]
                    early = source[source.experimental_time.isin(['D3'] if route_id == 'GSE221739' else SOURCE_TIMES)]
                    if target.empty or early.empty:
                        continue
                    old_s = source.groupby('gene_id').position_index.mean()
                    early_s = early.groupby('gene_id').position_index.mean()
                    fixed = target.index.intersection(old_s.index).intersection(early_s.index)
                    for ref, s in [('all_S_libraries', old_s), ('earliest_eligible_S', early_s)]:
                        reference_checks.append(dict(species=part.species.iloc[0], study_id=part.study_id.iloc[0],
                            route=route_id, reference=ref, weighting='equal_libraries',
                            **metrics(target.loc[fixed] - s.loc[fixed])))
        nodes.append(q)
    data = pd.concat(nodes, ignore_index=True)
    for (sp, study, route, layer), q in data.groupby(['species', 'study_id', 'route', 'layer']):
        if 'reference' in route or 'perturbation' in route:
            continue
        wide = q.pivot(index='gene_id', columns='node', values='position_index')
        if not {'S', 'P'} <= set(wide) or wide.P.notna().sum() == 0:
            continue
        states = [s for s in 'SREPM' if s in wide and wide[s].notna().any()]
        for background, ids in [('all_state_intersection', wide[states].dropna().index),
                                ('endpoint_measurable', wide[['S', 'P']].dropna().index)]:
            out = pd.DataFrame({'gene_id': ids, 'delta': (wide.loc[ids, 'P'] - wide.loc[ids, 'S']).to_numpy()})
            effects.append(out.assign(species=sp, study_id=study, route=route, layer=layer, background=background))
    effects = pd.concat(effects, ignore_index=True)
    old = read(F1 / 'analysis/figure1d_m_trial_v19/gene_effects.tsv.gz')
    old = old[old.node.eq('P')]
    now = effects[effects.layer.eq('supported50') & effects.background.eq('all_state_intersection')]
    z = old.merge(now, on=['species', 'study_id', 'route', 'gene_id'], how='outer', suffixes=('_old', '_new'), indicator=True)
    assert z._merge.eq('both').all() and (z.delta_old - z.delta_new).abs().max() < 1e-12
    CHECKS['figure1_endpoint_records_reproduced'] = len(z)
    rows = [dict(zip(['species', 'study_id', 'route', 'layer', 'background'], keys), **metrics(g.delta))
            for keys, g in effects.groupby(['species', 'study_id', 'route', 'layer', 'background'])]
    save(pd.DataFrame(rows), 'legacy_summary.tsv')
    save(effects, 'legacy_gene_effects.tsv.gz')
    save(pd.DataFrame(inventories), 'state_source_inventory.tsv')
    save(pd.DataFrame(reference_checks), 'reference_only_fixed_gene_summary.tsv')
    return data


def paired_effects():
    cols = KEYS + ['APA_delta', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b',
                   'library_key_a', 'library_key_b', 'n_mc_a', 'n_mc_b', 'n_joint_PAS']
    pairs = contrast_subset(read(STATE / 'matched_library_pair_effects.tsv.gz', usecols=cols))
    pm = pairs[pairs.contrast.eq('P_M')]
    assert pm.library_key_a.equals(pm.library_key_b) and pm.time_a.equals(pm.time_b)
    assert pairs.n_joint_PAS.ge(2).all()
    ps = pairs[pairs.design.eq('true_source')]
    assert ps.time_b.isin(SOURCE_TIMES).all()
    route, dates = collapse_pairs(pairs)
    frozen = contrast_subset(read(STATE / 'matched_study_effects.tsv.gz'))
    q = route.merge(frozen, on=KEYS, how='outer', indicator=True, validate='one_to_one')
    assert q._merge.eq('both').all()
    error = (q.delta - q.APA_delta).abs().max()
    assert error < 1e-8
    CHECKS.update(qualified_pairs=len(pairs), route_effects_reproduced=len(q), max_route_error=float(error))
    # Independently recover primary pair effects from saved library/state APA values.
    states = []
    state_cols = ['species', 'study_id', 'route', 'setting', 'gene_id', 'library_key', 'time', 'state', 'APA']
    for path in sorted(RAW.glob('gse*/apa_gex_states.tsv.gz')):
        d = read(path, usecols=state_cols)
        states.append(d[d.setting.eq('supported50_p08') & d.state.isin(['S', 'P', 'M'])])
    state = pd.concat(states, ignore_index=True)
    primary = pairs[pairs.setting.eq('supported50_p08')].copy()
    for suffix in ['a', 'b']:
        primary['state_' + suffix] = 'P' if suffix == 'a' else primary.contrast.str[-1]
        rename = {k: k + '_' + suffix for k in ['library_key', 'time', 'state', 'APA']}
        on = [rename.get(k, k) for k in state_cols if k != 'APA']
        primary = primary.merge(state.rename(columns=rename), on=on, how='left', validate='many_to_one')
    assert primary[['APA_a', 'APA_b']].notna().all().all()
    error = (primary.APA_a - primary.APA_b - primary.APA_delta).abs().max()
    assert error < 1e-8
    CHECKS.update(primary_pairs_reproduced=len(primary), max_pair_error=float(error))
    coverage = pairs[pairs.setting.eq('supported50_p08') & pairs.n_mc_a.ge(5) & pairs.n_mc_b.ge(5)].copy()
    coverage['setting'] = 'supported50_p08_mc5'
    cr, cd = collapse_pairs(coverage)
    route = pd.concat([route, cr], ignore_index=True)
    dates = pd.concat([dates, cd], ignore_index=True)
    save(route, 'qualified_route_effects.tsv.gz')
    save(dates, 'qualified_date_effects.tsv.gz')
    # Fixed-gene leave-target-date diagnostic. Dates are not biological replicates.
    rows = []
    for keys, g in dates[dates.setting.eq('supported50_p08')].groupby(['species', 'study_id', 'route', 'contrast', 'design']):
        w = g.pivot(index='gene_id', columns='time_a', values='APA_delta').dropna()
        if len(w.columns) < 2 or w.empty:
            continue
        for omitted in ['none'] + list(w.columns):
            x = w if omitted == 'none' else w.drop(columns=omitted)
            rows.append(dict(zip(['species', 'study_id', 'route', 'contrast', 'design'], keys),
                omitted_date=omitted, available_dates=len(w.columns), **metrics(x.mean(axis=1))))
    save(pd.DataFrame(rows), 'leave_target_date_summary.tsv')
    # Both genuine source libraries in Guan0618 remain separate alternatives.
    rows = []
    for route_id, g in pairs[pairs.setting.eq('supported50_p08') & pairs.design.eq('true_source')].groupby('route'):
        all_libs = g.library_key_b.nunique()
        joint = g.groupby(['gene_id', 'library_key_a']).library_key_b.transform('nunique').eq(all_libs)
        pieces = []
        for lib, q in g.groupby('library_key_b'):
            r, _ = collapse_pairs(q)
            pieces.append(r.assign(source_library=lib, source_time=q.time_b.iloc[0], background='available'))
            r, _ = collapse_pairs(q[joint.loc[q.index]])
            pieces.append(r.assign(source_library=lib, source_time=q.time_b.iloc[0], background='matched_target_libraries_all_sources'))
        result = pd.concat(pieces, ignore_index=True)
        fixed = set.intersection(*(set(q.gene_id) for _, q in result[result.background.eq('available')].groupby('source_library')))
        result = pd.concat([result, result[result.background.eq('available') & result.gene_id.isin(fixed)].assign(background='fixed_genes_all_sources')], ignore_index=True)
        for (sp, study, con, lib, time, background), q in result.groupby(['species', 'study_id', 'contrast', 'source_library', 'source_time', 'background']):
            rows.append(dict(species=sp, study_id=study, route=route_id, contrast=con,
                source_library=lib, source_time=time, background=background, **metrics(q.delta)))
    save(pd.DataFrame(rows), 'single_source_library_summary.tsv')
    decomposition = []
    for key, g in route.groupby(['species', 'study_id', 'route', 'contrast', 'design']):
        a = g[g.setting.eq('supported50_p08')].set_index('gene_id').delta
        b = g[g.setting.eq('expanded_p08')].set_index('gene_id').delta
        common = a.index.intersection(b.index)
        for subset, values in [('supported_common', a.loc[common]), ('expanded_common', b.loc[common]),
                               ('supported_only', a.loc[a.index.difference(b.index)]),
                               ('expanded_only', b.loc[b.index.difference(a.index)])]:
            if len(values):
                decomposition.append(dict(zip(['species', 'study_id', 'route', 'contrast', 'design'], key), subset=subset, **metrics(values)))
    save(pd.DataFrame(decomposition), 'pas_layer_coverage_decomposition.tsv')
    return route


def orthologs():
    path = ROOT / 'analysis_code/release_2026-09-24/sources/MOUSE/sources/local/apa-cross-species-aligned-v1/ensembl-116-human-mouse-one-to-one-orthologs.tsv'
    x = read(path)
    x = x[x.documented_one_to_one & x.human_gene_id.notna()]
    assert x.human_gene_id.is_unique and x.mouse_gene_id.is_unique
    return x


def backgrounds(route, orth):
    groups = ['species', 'study_id', 'route', 'contrast', 'design']
    fixed = {k: set.intersection(*(set(q.gene_id) for _, q in g.groupby('setting')))
             for k, g in route.groupby(groups)}
    rows = []
    for key, g in route.groupby(groups + ['setting']):
        base = dict(zip(groups + ['setting'], key))
        rows.append(dict(base, background='available', **metrics(g.delta)))
        q = g[g.gene_id.isin(fixed[key[:-1]])]
        rows.append(dict(base, background='fixed_settings_within_route', **metrics(q.delta)))
    eligible = route[route.design.isin(['true_source', 'same_library'])]
    common_sets = []
    for (sp, con), g in eligible.groupby(['species', 'contrast']):
        ids = set.intersection(*(set(q.gene_id) for _, q in g.groupby(['route', 'setting'])))
        for gene in sorted(ids):
            common_sets.append(dict(species=sp, contrast=con, gene_id=gene))
        for key, q in g.groupby(groups + ['setting']):
            rows.append(dict(zip(groups + ['setting'], key), background='common_within_species_all_settings', **metrics(q[q.gene_id.isin(ids)].delta)))
    common = pd.DataFrame(common_sets)
    save(common, 'common_gene_backgrounds.tsv.gz')
    for con, c in common.groupby('contrast'):
        ids = {sp: set(q.gene_id) for sp, q in c.groupby('species')}
        mapping = orth[orth.human_gene_id.isin(ids.get('human', set())) & orth.mouse_gene_id.isin(ids.get('mouse', set()))]
        for sp in ['human', 'mouse']:
            for key, q in eligible[eligible.species.eq(sp) & eligible.contrast.eq(con)].groupby(groups + ['setting']):
                rows.append(dict(zip(groups + ['setting'], key), background='common_one_to_one_ortholog_all_settings', **metrics(q[q.gene_id.isin(mapping[sp + '_gene_id'])].delta)))
    summary = pd.DataFrame(rows)
    save(summary, 'global_direction_summary.tsv')
    # Descriptive stability, not significance or a post-hoc numerical threshold.
    rows = []
    for key, g in summary.groupby(groups):
        sign = lambda x: int(np.sign(x.iloc[0])) if (np.sign(x) == np.sign(x.iloc[0])).all() else 0
        fixed = g[g.background.ne('available')]
        rows.append(dict(zip(groups, key), mean_stable_sign=sign(g['mean']), fixed_background_mean_sign=sign(fixed['mean']), median_stable_sign=sign(g['median']),
            excess002_stable_sign=sign(g.distal_excess_002), excess005_stable_sign=sign(g.distal_excess_005),
            mean_min=g['mean'].min(), mean_max=g['mean'].max(), n_scenarios=len(g)))
    stability = pd.DataFrame(rows)
    sources = read(HERE / 'single_source_library_summary.tsv')
    for metric, column in [('mean', 'mean'), ('median', 'median'), ('excess002', 'distal_excess_002'), ('excess005', 'distal_excess_005')]:
        stability[metric + '_source_checked_sign'] = stability[metric + '_stable_sign']
        for i, row in stability[stability.design.eq('true_source')].iterrows():
            q = sources[sources.route.eq(row.route) & sources.background.ne('available')]
            if not np.sign(q[column]).eq(row[metric + '_stable_sign']).all():
                stability.loc[i, metric + '_source_checked_sign'] = 0
    save(stability, 'direction_stability.tsv')
    return summary, stability


def study_event(values, threshold):
    """Guan series are one study; a subthreshold/conflicting series blocks support."""
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return np.nan
    return 1 if (x >= threshold).all() else -1 if (x <= -threshold).all() else 0


def recurrence_mapping(route, stability):
    frozen = read(ROOT / 'analysis_code/recurrence_correction_2026-09-25/sources/MOUSE/outputs/figure3-baseline-wpp-null-2026-09-23/recurrent_gene_sets.tsv.gz')
    frozen = frozen[frozen.model.eq('original') & frozen.omitted_study.eq('none')]
    assert not frozen.duplicated(['species', 'contrast', 'threshold', 'gene_id']).any()
    primary = route[route.setting.eq('supported50_p08') & route.design.isin(['true_source', 'same_library'])]
    joined = []
    for threshold, f in frozen.groupby('threshold'):
        r = primary.merge(f[['species', 'contrast', 'gene_id', 'direction']], on=['species', 'contrast', 'gene_id'], validate='many_to_one')
        r['threshold'] = threshold
        r['recurrent_sign'] = r.direction.map({'positive': 1, 'negative': -1})
        ev = r.groupby(['species', 'contrast', 'study_id', 'gene_id']).delta.transform(lambda x: study_event(x, threshold))
        r['study_event'] = ev
        r['supports_frozen_event'] = ev.eq(r.recurrent_sign)
        r = r.merge(stability, on=['species', 'study_id', 'route', 'contrast', 'design'], validate='many_to_one')
        for rec in f.itertuples():
            z = r[r.species.eq(rec.species) & r.contrast.eq(rec.contrast) & r.gene_id.eq(rec.gene_id)]
            for study, q in z.groupby('study_id'):
                recorded = getattr(rec, study)
                actual = {1: '+', -1: '-', 0: '0'}[int(q.study_event.iloc[0])]
                assert recorded == actual, (study, rec.gene_id, recorded, actual)
        joined.append(r)
    joined = pd.concat(joined, ignore_index=True)
    save(joined, 'figure3_frozen_member_background_map.tsv.gz')
    outcomes = []
    for key, g in joined.groupby(['species', 'contrast', 'threshold', 'gene_id', 'gene_symbol']):
        support = g[g.supports_frozen_event]
        row = dict(zip(['species', 'contrast', 'threshold', 'gene_id', 'gene_symbol'], key), supporting_studies=support.study_id.nunique())
        for metric in ['mean', 'median', 'excess002', 'excess005']:
            col = metric + '_stable_sign'
            positive = set(support.loc[support[col].eq(1), 'study_id'])
            negative = set(support.loc[support[col].eq(-1), 'study_id'])
            row['cross_opposite_' + metric] = any(a != b for a in positive for b in negative)
        outcomes.append(row)
    outcomes = pd.DataFrame(outcomes)
    save(outcomes, 'figure3_member_outcomes.tsv.gz')
    save(outcomes.groupby(['species', 'contrast', 'threshold']).agg(
        frozen_genes=('gene_id', 'size'), cross_opposite_mean=('cross_opposite_mean', 'sum'),
        cross_opposite_median=('cross_opposite_median', 'sum'), cross_opposite_excess002=('cross_opposite_excess002', 'sum'),
        cross_opposite_excess005=('cross_opposite_excess005', 'sum')).reset_index(), 'figure3_mapping_summary.tsv')


def figure4_mapping(route, stability):
    frozen = read(F4 / 'data/figure4-ortholog-null-v24/ortholog_pair_features.tsv')
    members = []
    for con in ['PS', 'PM']:
        q = frozen[frozen['h_' + con] & frozen['m_' + con]].copy()
        assert len(q) == (18 if con == 'PS' else 10)
        q['contrast'] = 'P_S' if con == 'PS' else 'P_M'
        members.append(q[['human_gene_id', 'mouse_gene_id', 'human_symbol', 'contrast']])
    members = pd.concat(members, ignore_index=True)
    assert members.human_gene_id.nunique() == 26
    save(members, 'figure4_frozen_members.tsv')
    # Recover actual published support, including its S definition and Guan0618 only.
    codes = {'W': 'WOT', 'C': 'CIP', 'H1': 'GSE118258', 'H2': 'GSE150311', 'H3': '178325_0618', 'H4': 'GSE221739'}
    nodes = read(F4 / 'data/figure4-program-v9/all_mapped_nodes.tsv.gz')
    n = nodes[nodes.layer.eq('supported50') & nodes.measurement_pass & nodes.code.isin(codes)]
    wide = n.pivot(index=['species', 'study_id', 'route', 'gene_id', 'human_gene_id', 'code'], columns='node', values='position_index')
    ps = (wide.P - wide.S).dropna().rename('published_delta').reset_index().assign(contrast='P_S')
    pm = read(F4 / 'data/figure4-program-v9/paired_pm_study_effects.tsv.gz')
    pm = pm[pm.layer.eq('supported50') & pm.quality_pass & pm.code.isin(codes) &
            np.where(pm.species.eq('mouse'), pm.n_sample_units.ge(2), pm.n_sample_units.ge(1))].copy()
    pm = pm.rename(columns={'mean_delta': 'published_delta'}).assign(contrast='P_M')
    published = pd.concat([ps, pm[ps.columns]], ignore_index=True)
    published = members.merge(published, on=['human_gene_id', 'contrast'], validate='one_to_many')
    published['published_support'] = published.published_delta.le(-.02)
    for (gene, con), q in published.groupby(['human_gene_id', 'contrast']):
        sup = q[q.published_support]
        assert sup[sup.species.eq('mouse')].study_id.nunique() == 2
        assert sup[sup.species.eq('human')].study_id.nunique() >= 2
    save(published, 'figure4_published_support.tsv')
    primary = route[route.setting.eq('supported50_p08') & route.design.isin(['true_source', 'same_library'])]
    mapped = published.merge(primary, on=['species', 'study_id', 'route', 'gene_id', 'contrast'], how='left', validate='one_to_one')
    mapped['supports_frozen_proximal'] = mapped.published_support & mapped.delta.le(-.02)
    # Additional measured Guan series are a separate study-consensus sensitivity,
    # not silently substituted for the published Figure 4 route.
    consensus = primary.groupby(['species', 'study_id', 'gene_id', 'contrast']).delta.agg(lambda x: study_event(x, .02)).rename('all_series_event').reset_index()
    mapped = mapped.merge(consensus, on=['species', 'study_id', 'gene_id', 'contrast'], how='left', validate='many_to_one')
    mapped['all_series_support'] = mapped.supports_frozen_proximal & mapped.all_series_event.eq(-1)
    every_setting = route[route.design.isin(['true_source', 'same_library'])].groupby(
        ['species', 'study_id', 'gene_id', 'contrast', 'setting']).delta.agg(lambda x: study_event(x, .02)).reset_index()
    robust = every_setting.groupby(['species', 'study_id', 'gene_id', 'contrast']).delta.agg(
        all_settings_event=lambda x: -1 if len(x) == 4 and x.eq(-1).all() else 0).reset_index()
    mapped = mapped.merge(robust, on=['species', 'study_id', 'gene_id', 'contrast'], how='left', validate='many_to_one')
    mapped = mapped.merge(stability, on=['species', 'study_id', 'route', 'contrast', 'design'], how='left', validate='many_to_one')
    save(mapped, 'figure4_member_background_map.tsv')
    outcomes = []
    for key, g in mapped.groupby(['human_gene_id', 'human_symbol', 'contrast']):
        support = g[g.supports_frozen_proximal]
        row = dict(zip(['human_gene_id', 'human_symbol', 'contrast'], key),
            supporting_human_studies=support[support.species.eq('human')].study_id.nunique(),
            supporting_mouse_studies=support[support.species.eq('mouse')].study_id.nunique(),
            supporting_routes=';'.join(sorted(support.route.dropna().unique())))
        for metric in ['mean', 'median', 'excess002', 'excess005']:
            col = metric + '_stable_sign'
            pos = set(support.loc[support[col].eq(1), 'study_id'])
            neg = set(support.loc[support[col].eq(-1), 'study_id'])
            row['cross_opposite_' + metric] = any(a != b for a in pos for b in neg)
            row['positive_background_support_' + metric] = ';'.join(sorted(pos))
            row['negative_background_support_' + metric] = ';'.join(sorted(neg))
            col = metric + '_source_checked_sign'
            pos = set(support.loc[support[col].eq(1), 'study_id'])
            neg = set(support.loc[support[col].eq(-1), 'study_id'])
            row['cross_opposite_' + metric + '_source_checked'] = any(a != b for a in pos for b in neg)
        strict = support[support.all_series_support]
        pos = set(strict.loc[strict.mean_source_checked_sign.eq(1), 'study_id'])
        neg = set(strict.loc[strict.mean_source_checked_sign.eq(-1), 'study_id'])
        row['cross_opposite_mean_source_and_series_checked'] = any(a != b for a in pos for b in neg)
        strict = strict[strict.all_settings_event.eq(-1)]
        pos = set(strict.loc[strict.mean_source_checked_sign.eq(1), 'study_id'])
        neg = set(strict.loc[strict.mean_source_checked_sign.eq(-1), 'study_id'])
        row['cross_opposite_mean_all_four_settings'] = any(a != b for a in pos for b in neg)
        outcomes.append(row)
    outcomes = pd.DataFrame(outcomes)
    save(outcomes, 'figure4_member_outcomes.tsv')
    save(outcomes.groupby('contrast').agg(frozen_genes=('human_gene_id', 'size'),
        cross_opposite_mean=('cross_opposite_mean', 'sum'), cross_opposite_median=('cross_opposite_median', 'sum'),
        cross_opposite_excess002=('cross_opposite_excess002', 'sum'), cross_opposite_excess005=('cross_opposite_excess005', 'sum'),
        cross_opposite_mean_source_checked=('cross_opposite_mean_source_checked', 'sum'),
        cross_opposite_mean_source_and_series_checked=('cross_opposite_mean_source_and_series_checked', 'sum'),
        cross_opposite_mean_all_four_settings=('cross_opposite_mean_all_four_settings', 'sum')).reset_index(), 'figure4_mapping_summary.tsv')


def main():
    before = protected()
    (HERE / 'protected_sources_sha256.json').write_text(json.dumps(before, indent=2) + '\n')
    legacy()
    print('Original library means and Figure 1 reproduced.', flush=True)
    route = paired_effects()
    print('Qualified contrasts independently reconstructed.', flush=True)
    summary, stability = backgrounds(route, orthologs())
    recurrence_mapping(route, stability)
    figure4_mapping(route, stability)
    states = read(ROOT / 'notes/2026-09-22-integration-support-layout/source_data/states.tsv')
    save(states[states.family.eq('shared_tp')], 'endpoint_annotation_evidence.tsv')
    inventory = read(ROOT / 'table/S2-sequencing-libraries.tsv')
    save(inventory[inventory.study.isin(['GSE118258', 'GSE150311', 'GSE178325', 'GSE221739', 'GSE103221', 'GSE115943'])][
        ['study', 'species', 'sample_id', 'library_id', 'geo_sample_title', 'geo_source_name', 'geo_characteristics']], 'original_sample_metadata.tsv')
    assert protected() == before
    CHECKS['formal_sources_unchanged'] = True
    CHECKS['protected_files'] = len(before)
    (HERE / 'input_sha256.json').write_text(json.dumps(INPUTS, indent=2) + '\n')
    (HERE / 'verification.json').write_text(json.dumps(CHECKS, indent=2) + '\n')
    print(stability.to_string(index=False), flush=True)
    print(json.dumps(CHECKS, indent=2), flush=True)


if __name__ == '__main__':
    main()
