"""Connect source and coexisting-state contrasts on identical qualified contributors."""
from pathlib import Path
import hashlib
import json
import time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MOUSE = Path('/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot')
SOURCE = MOUSE/'outputs/figure3-apa-gex-state-analysis-2026-09-19'
RAW = MOUSE/'analysis/figure3-apa-gex-state-analysis-2026-09-19/raw'
RECURRENCE = MOUSE/'outputs/figure3-baseline-wpp-null-2026-09-23'
STUDIES = ['GSE118258', 'GSE103221', 'GSE115943']
SETTINGS = ['supported50_p08', 'supported50_p09', 'expanded_p08', 'supported50_p08_mc5']
KEY = ['species', 'study_id', 'route', 'setting', 'gene_id']
PAIR = KEY+['library_key_a', 'library_key_b']
SOURCES = ['D0', 'source', 'Stage_I_day_0']
CATEGORIES = ['dual_reference_aligned', 'reference_opposed', 'PM_event_PS_subthreshold',
              'shared_source_shift_PM_subthreshold', 'PS_event_PM_subthreshold_otherMS',
              'both_P_contrasts_subthreshold']
INPUTS = {}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path, **kwargs):
    INPUTS[str(path)] = sha(path)
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, name):
    frame.to_csv(HERE/name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def protected():
    return {str(p.relative_to(ROOT)): sha(p) for folder in ['main', 'supplement', 'Fig', 'table', 'output/pdf']
            for p in (ROOT/folder).rglob('*') if p.is_file()}


def aggregate(frame, values, keys=KEY):
    units = frame.groupby(keys+['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'], dropna=False)[values].mean()
    source_units = units.groupby(keys+['time_a', 'sample_unit_a', 'time_b'], dropna=False).mean()
    source_dates = source_units.groupby(keys+['time_a', 'sample_unit_a'], dropna=False).mean()
    target_units = source_dates.groupby(keys+['time_a'], dropna=False).mean()
    return target_units.groupby(keys, dropna=False).mean().reset_index()


def event(values, threshold):
    values = np.asarray(values, float)
    assert np.isfinite(values).all(), 'Missing measurements must not become subthreshold events'
    return np.where(values >= threshold, 1, np.where(values <= -threshold, -1, 0))


def classify(ps, ms, pm):
    ps, ms, pm = np.asarray(ps), np.asarray(ms), np.asarray(pm)
    assert all(np.isin(x, [-1, 0, 1]).all() for x in [ps, ms, pm])
    conditions = [(ps*pm == 1), (ps*pm == -1), (ps == 0) & (pm != 0),
                  (ps != 0) & (ms == ps) & (pm == 0),
                  (ps != 0) & (pm == 0) & (ms != ps), (ps == 0) & (pm == 0)]
    assert np.all(np.sum(conditions, axis=0) == 1)
    return np.select(conditions, CATEGORIES, default='invalid')


def make_triplets(pairs):
    ps = pairs[pairs.contrast.eq('P_S') & pairs.design.eq('true_source')].copy()
    ms = pairs[pairs.contrast.eq('M_S') & pairs.design.eq('true_source')].copy()
    pm = pairs[pairs.contrast.eq('P_M') & pairs.design.eq('same_library')].copy()
    assert ps.time_b.isin(SOURCES).all() and ms.time_b.isin(SOURCES).all()
    for col in ['library_key', 'time', 'sample_unit']:
        assert pm[col+'_a'].eq(pm[col+'_b']).all(), col
    assert not ps.duplicated(PAIR).any() and not ms.duplicated(PAIR).any()
    assert not pm.duplicated(KEY+['library_key_a']).any()
    t = ps.merge(ms, on=PAIR, validate='one_to_one', suffixes=('_PS', '_MS'))
    for col in ['time_a', 'time_b', 'sample_unit_a', 'sample_unit_b', 'replicate_status_a', 'replicate_status_b']:
        assert t[col+'_PS'].eq(t[col+'_MS']).all(), col
        t[col] = t[col+'_PS']
    for col in ['n_mc_b', 'n_cells_b', 'pas_counts_b']:
        assert t[col+'_PS'].eq(t[col+'_MS']).all(), col
    t = t.merge(pm.drop(columns='library_key_b'), on=KEY+['library_key_a'],
                validate='many_to_one', suffixes=('', '_PM'))
    for col in ['time_a', 'sample_unit_a']:
        assert t[col].eq(t[col+'_PM']).all(), col
    for col in ['n_mc', 'n_cells', 'pas_counts']:
        assert t[col+'_a_PS'].eq(t[col+'_a']).all(), col
        assert t[col+'_a_MS'].eq(t[col+'_b']).all(), col
        t['P_'+col] = t[col+'_a_PS']
        t['M_'+col] = t[col+'_a_MS']
        t['S_'+col] = t[col+'_b_PS']
    t['gene_symbol'] = t.gene_symbol_PS
    t['original_pair_PS'], t['original_pair_MS'], t['original_pair_PM'] = t.APA_delta_PS, t.APA_delta_MS, t.APA_delta
    t['joint_PAS_PS'], t['joint_PAS_MS'], t['joint_PAS_PM'] = t.n_joint_PAS_PS, t.n_joint_PAS_MS, t.n_joint_PAS
    keep = PAIR+['gene_symbol', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b',
                 'replicate_status_a', 'replicate_status_b']
    keep += [f'{st}_{col}' for st in ['P', 'M', 'S'] for col in ['n_mc', 'n_cells', 'pas_counts']]
    keep += ['original_pair_PS', 'original_pair_MS', 'original_pair_PM', 'joint_PAS_PS', 'joint_PAS_MS', 'joint_PAS_PM']
    return t[keep].copy()


def add_states(triplets):
    states = []
    cols = KEY+['library_key', 'time', 'sample_unit', 'state', 'APA', 'n_mc', 'n_cells', 'pas_counts']
    for study in STUDIES:
        matches = list(RAW.glob(study.lower()+'*/apa_gex_states.tsv.gz'))
        assert len(matches) == 1
        d = read(matches[0], usecols=cols)
        states.append(d[d.state.isin(['P', 'M', 'S'])])
    states = pd.concat(states, ignore_index=True)
    state_keys = KEY+['library_key', 'time', 'sample_unit', 'state']
    assert not states.duplicated(state_keys).any()
    t = triplets.copy()
    for st, side in [('P', 'a'), ('M', 'a'), ('S', 'b')]:
        part = states[states.state.eq(st)].drop(columns='state').rename(columns={
            'library_key': 'library_key_'+side, 'time': 'time_'+side, 'sample_unit': 'sample_unit_'+side,
            'APA': st+'_WPP', 'n_mc': st+'_raw_n_mc', 'n_cells': st+'_raw_n_cells',
            'pas_counts': st+'_raw_pas_counts'})
        t = t.merge(part, on=KEY+['library_key_'+side, 'time_'+side, 'sample_unit_'+side], how='left', validate='many_to_one')
        assert t[st+'_WPP'].notna().all()
        for col in ['n_mc', 'n_cells', 'pas_counts']:
            assert t[st+'_'+col].eq(t[st+'_raw_'+col]).all()
        t = t.drop(columns=[st+'_raw_'+col for col in ['n_mc', 'n_cells', 'pas_counts']])
    for contrast, a, b in [('PS', 'P', 'S'), ('MS', 'M', 'S'), ('PM', 'P', 'M')]:
        t['delta_'+contrast] = t[a+'_WPP']-t[b+'_WPP']
        assert np.max(np.abs(t['delta_'+contrast]-t['original_pair_'+contrast])) < 1e-8
    assert np.allclose(t.delta_PS, t.delta_PM+t.delta_MS, rtol=0, atol=1e-12)
    return t


def frozen_membership():
    f = read(RECURRENCE/'recurrent_gene_sets.tsv.gz')
    f = f[f.model.eq('baseline_q5_capacity') & f.omitted_study.eq('none') & f.threshold.eq(.02)].copy()
    f['frozen_direction'] = f.direction.map({'positive': 1, 'negative': -1})
    assert not f.duplicated(['species', 'contrast', 'gene_id']).any()
    assert f.groupby(['species', 'contrast']).size().to_dict() == {
        ('human', 'P_M'): 170, ('human', 'P_S'): 510, ('mouse', 'P_M'): 216, ('mouse', 'P_S'): 488}
    return f[['species', 'contrast', 'gene_id', 'gene_symbol', 'frozen_direction']].rename(columns={'contrast': 'cohort'})


def summaries(effect, cohorts):
    rows, configurations, tracking = [], [], []
    for setting in SETTINGS:
        e = effect[effect.setting.eq(setting)].copy()
        for threshold in [.02, .05]:
            for name in ['PS', 'MS', 'PM']:
                e['event_'+name] = event(e['delta_'+name], threshold)
            e['configuration'] = classify(e.event_PS, e.event_MS, e.event_PM)
            e['threshold'] = threshold
            configurations.append(e.copy())
            for cohort in ['all_jointly_measurable', 'P_S', 'P_M']:
                sub = e if cohort == 'all_jointly_measurable' else e.merge(
                    cohorts[cohorts.cohort.eq(cohort)][['species', 'gene_id']], on=['species', 'gene_id'], validate='many_to_one')
                for (species, study), g in sub.groupby(['species', 'study_id']):
                    for label in CATEGORIES:
                        rows.append(dict(setting=setting, threshold=threshold, cohort=cohort, species=species,
                                         study_id=study, configuration=label, n=int(g.configuration.eq(label).sum()),
                                         evaluable_genes=len(g)))
            # Explicitly retain original genes at studies lacking a genuine source or complete triple.
            study_frame = pd.DataFrame([('human', 'GSE118258'), ('human', 'GSE178325'),
                                        ('human', 'GSE150311'), ('human', 'GSE221739'),
                                        ('mouse', 'GSE103221'), ('mouse', 'GSE115943')], columns=['species', 'study_id'])
            grid = cohorts.merge(study_frame, on='species').merge(
                e.drop(columns='gene_symbol'), on=['species', 'gene_id', 'study_id'], how='left', validate='many_to_one')
            grid['setting'], grid['threshold'] = setting, threshold
            grid['evaluable'] = grid.delta_PM.notna()
            grid['missing_reason'] = np.where(grid.evaluable, 'none', np.where(
                grid.study_id.isin(['GSE150311', 'GSE221739']), 'no_qualified_true_source', np.where(
                    grid.study_id.eq('GSE178325'), 'no_qualified_same_library_PM', 'no_complete_qualified_triplet')))
            for name in ['PS', 'MS', 'PM']:
                grid[name+'_frozen_aligned'] = (grid['event_'+name]*grid.frozen_direction).eq(1).where(grid.evaluable)
            tracking.append(grid)
    configs, trace = pd.concat(configurations, ignore_index=True), pd.concat(tracking, ignore_index=True)
    save(configs, 'gene_configurations.tsv.gz')
    save(pd.DataFrame(rows), 'configuration_summary.tsv')
    save(configs.groupby(['setting', 'threshold', 'species', 'study_id', 'event_PS', 'event_MS', 'event_PM']).size().rename('n').reset_index(),
         'complete_sign_configurations.tsv')
    save(trace, 'frozen_member_trace.tsv.gz')
    counts = []
    for (setting, threshold, cohort, species), g in trace.groupby(['setting', 'threshold', 'cohort', 'species']):
        q = g.groupby('gene_id').agg(n_evaluable=('evaluable', 'sum'),
                n_PS_aligned=('PS_frozen_aligned', lambda x: x.fillna(False).astype(bool).sum()),
                n_PM_aligned=('PM_frozen_aligned', lambda x: x.fillna(False).astype(bool).sum()))
        g = g.assign(dual=g.PS_frozen_aligned.eq(True) & g.PM_frozen_aligned.eq(True))
        q['n_dual_aligned'] = g.groupby('gene_id').dual.sum()
        for minimum in [1, 2]:
            counts.append(dict(setting=setting, threshold=threshold, cohort=cohort, species=species,
                required_studies=minimum, frozen_genes=len(q), evaluable_genes=int(q.n_evaluable.ge(minimum).sum()),
                PS_frozen_aligned_genes=int(q.n_PS_aligned.ge(minimum).sum()),
                PM_frozen_aligned_genes=int(q.n_PM_aligned.ge(minimum).sum()),
                dual_frozen_aligned_genes=int(q.n_dual_aligned.ge(minimum).sum())))
    save(pd.DataFrame(counts), 'frozen_cohort_summary.tsv')
    pairs = configs[configs.study_id.eq('GSE103221')].merge(configs[configs.study_id.eq('GSE115943')],
        on=['species', 'gene_id', 'setting', 'threshold'], suffixes=('_CIP', '_WOT'), validate='one_to_one')
    pairs['same_configuration'] = pairs.configuration_CIP.eq(pairs.configuration_WOT)
    pairs['dual_same_direction'] = (pairs.configuration_CIP.eq('dual_reference_aligned') &
        pairs.configuration_WOT.eq('dual_reference_aligned') & pairs.event_PS_CIP.eq(pairs.event_PS_WOT))
    save(pairs, 'mouse_cross_study_configurations.tsv.gz')
    save(pairs.groupby(['setting', 'threshold']).agg(joint_genes=('gene_id', 'size'),
        same_configuration=('same_configuration', 'sum'), dual_same_direction=('dual_same_direction', 'sum')).reset_index(),
        'mouse_cross_study_summary.tsv')


def main():
    started = time.perf_counter()
    before = protected()
    pairs = read(SOURCE/'matched_library_pair_effects.tsv.gz')
    pairs = pairs[pairs.study_id.isin(STUDIES) & (
        (pairs.contrast.isin(['P_S', 'M_S']) & pairs.design.eq('true_source')) |
        (pairs.contrast.eq('P_M') & pairs.design.eq('same_library')))].copy()
    assert pairs.n_joint_PAS.ge(2).all()
    native = read(SOURCE/'matched_study_effects.tsv.gz')
    rebuilt = aggregate(pairs, ['APA_delta'], KEY+['contrast', 'design'])
    comparison = rebuilt.merge(native, on=KEY+['contrast', 'design'], validate='one_to_one', suffixes=('_rebuilt', '_native'))
    assert len(comparison) == len(rebuilt)
    native_error = np.max(np.abs(comparison.APA_delta_rebuilt-comparison.APA_delta_native))
    assert native_error < 1e-8
    t = add_states(make_triplets(pairs))
    high = t[t.setting.eq('supported50_p08') & t[['P_n_mc', 'M_n_mc', 'S_n_mc']].min(axis=1).ge(5)].copy()
    high['setting'] = 'supported50_p08_mc5'
    t = pd.concat([t, high], ignore_index=True)
    assert set(t.setting) == set(SETTINGS) and not t.duplicated(PAIR).any()
    save(t, 'matched_library_source_triplets.tsv.gz')
    values = ['P_WPP', 'M_WPP', 'S_WPP', 'delta_PS', 'delta_MS', 'delta_PM']
    effects = aggregate(t, values)
    assert np.allclose(effects.delta_PS, effects.delta_PM+effects.delta_MS, rtol=0, atol=1e-12)
    support = t.groupby(KEY).agg(n_triplets=('library_key_a', 'size'), n_target_libraries=('library_key_a', 'nunique'),
        n_source_libraries=('library_key_b', 'nunique'), n_target_dates=('time_a', 'nunique'),
        n_target_units=('sample_unit_a', 'nunique'), P_min_mc=('P_n_mc', 'min'),
        M_min_mc=('M_n_mc', 'min'), S_min_mc=('S_n_mc', 'min'), gene_symbol=('gene_symbol', 'first')).reset_index()
    effects = effects.merge(support, on=KEY, validate='one_to_one')
    for short, contrast in [('PS', 'P_S'), ('MS', 'M_S'), ('PM', 'P_M')]:
        old = native[native.contrast.eq(contrast) & native.design.eq('same_library' if short == 'PM' else 'true_source')]
        # MC5 removes contributors; its original comparator is the unfiltered primary setting.
        old = pd.concat([old, old[old.setting.eq('supported50_p08')].assign(setting='supported50_p08_mc5')])
        effects = effects.merge(old[KEY+['APA_delta']].rename(columns={'APA_delta': 'original_delta_'+short}),
                                on=KEY, how='left', validate='one_to_one')
    save(effects, 'matched_state_summaries.tsv')
    cohorts = frozen_membership()
    save(cohorts, 'frozen_cohorts.tsv')
    summaries(effects, cohorts)
    availability = native[native.study_id.isin(STUDIES) & (
        (native.contrast.isin(['P_S', 'M_S']) & native.design.eq('true_source')) |
        (native.contrast.eq('P_M') & native.design.eq('same_library')))].copy()
    availability = pd.concat([availability, availability[availability.setting.eq('supported50_p08')].assign(setting='supported50_p08_mc5')])
    coverage = availability.groupby(['setting', 'species', 'study_id', 'contrast']).gene_id.nunique().rename('original_genes').reset_index()
    matched = effects.groupby(['setting', 'species', 'study_id']).gene_id.nunique().rename('complete_triplet_genes').reset_index()
    coverage = coverage.merge(matched, on=['setting', 'species', 'study_id'], how='outer', validate='many_to_one')
    save(coverage, 'coverage_ledger.tsv')
    assert protected() == before
    manifest = dict(status='computed; manuscript and figures unchanged', inputs_sha256=INPUTS,
        protected_sha256=before, protocol_sha256=sha(HERE/'PROTOCOL.zh.md'), analysis_sha256=sha(Path(__file__)),
        reconstructed_native_records=len(rebuilt), maximum_native_error=float(native_error),
        maximum_triplet_identity_error=float(abs(effects.delta_PS-effects.delta_PM-effects.delta_MS).max()),
        settings=SETTINGS, thresholds=[.02, .05], source_qualified_studies=STUDIES,
        matched_triplets=len(t), matched_route_gene_records=len(effects), elapsed_seconds=time.perf_counter()-started)
    (HERE/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: v for k, v in manifest.items() if k not in ['inputs_sha256', 'protected_sha256']}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
