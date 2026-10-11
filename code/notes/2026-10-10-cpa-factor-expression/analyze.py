"""Descriptive CPA RNA expression using the frozen physical-library contract."""
from pathlib import Path
import ast
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
B = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
KEY = ['species', 'study_id', 'route', 'factor']
LIBKEY = KEY + ['library_key', 'state']
NUMERIC = ['count', 'total_gex', 'pooled_CPM', 'mean_log2CPM1', 'mean_CPM', 'detected_metacells', 'detected_fraction', 'n_metacells', 'n_cells']


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.12g', compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def legacy_library_rows():
    path = ROOT / 'notes/2026-10-02-rbp-expression/analyze.py'
    module = ast.parse(path.read_text())
    node = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'library_rows')
    namespace = {'np': np, 'pd': pd}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['library_rows']


def documented_source_flags(libraries, design):
    """Only design-verified starting libraries may label the primary source."""
    d = design[['study_id', 'route', 'library_key', 'state', 'time']].drop_duplicates()
    assert not d.duplicated(['study_id', 'route', 'library_key', 'state']).any()
    z = libraries.merge(d.rename(columns={'time': 'design_time'}), on=['study_id', 'route', 'library_key', 'state'], how='left', validate='many_to_one')
    z['documented_source'] = z.state.eq('S') & z.design_time.astype(str).isin(['D0', 'source', 'Stage_I_day_0'])
    return z


def summarize_states(libraries):
    out = []
    for key, raw in libraries.groupby(KEY):
        for state in ['S', 'R', 'E', 'P', 'M', 'N']:
            q = raw[raw.state.eq(state) & raw.qualified & raw.mapped]
            r = dict(zip(KEY, key), state=state, n_libraries=len(q), n_times=q.time.nunique(), n_metacells=int(q.n_metacells.sum()), n_cells=int(q.n_cells.sum()), status='available' if len(q) else 'absent_or_unmapped')
            for column in ['mean_log2CPM1', 'pooled_CPM', 'detected_fraction']:
                r[column] = q[column].mean()
            r.update(time_equal_log2CPM1=q.groupby('time').mean_log2CPM1.mean().mean(), time_equal_pooled_CPM=q.groupby('time').pooled_CPM.mean().mean(), library_keys=json.dumps(q.library_key.tolist()))
            out.append(r)
    return pd.DataFrame(out)


def comparison(key, name, target, reference, paired=False, required_target=0, required_reference=0):
    r = dict(zip(KEY, key), contrast=name, n_target_libraries=len(target), n_reference_libraries=len(reference), n_pairs=0, paired=paired, target_library_keys=json.dumps(target.library_key.tolist()), reference_library_keys=json.dumps(reference.library_key.tolist()), status='unavailable', delta_log2CPM1=np.nan, time_equal_delta=np.nan, target_log2CPM1=np.nan, reference_log2CPM1=np.nan, target_CPM=np.nan, reference_CPM=np.nan, delta_pooled_CPM=np.nan)
    pair_rows = []
    if target.empty or reference.empty or (required_target and len(target) != required_target) or (required_reference and len(reference) != required_reference):
        return r, pair_rows
    if paired:
        pairs = target.merge(reference, on='library_key', suffixes=('_target', '_reference'), validate='one_to_one')
        if pairs.empty:
            return r, pair_rows
        assert pairs.time_target.eq(pairs.time_reference).all()
        pairs['delta'] = pairs.mean_log2CPM1_target - pairs.mean_log2CPM1_reference
        for z in pairs.itertuples():
            pair_rows.append(dict(zip(KEY, key), contrast=name, library_key=z.library_key, time=z.time_target, target_log2CPM1=z.mean_log2CPM1_target, reference_log2CPM1=z.mean_log2CPM1_reference, delta_log2CPM1=z.delta, target_CPM=z.pooled_CPM_target, reference_CPM=z.pooled_CPM_reference))
        target = target[target.library_key.isin(pairs.library_key)]
        reference = reference[reference.library_key.isin(pairs.library_key)]
        r.update(n_pairs=len(pairs), n_target_libraries=len(pairs), n_reference_libraries=len(pairs), delta_log2CPM1=pairs.delta.mean(), time_equal_delta=pairs.groupby('time_target').delta.mean().mean(), status='descriptive_paired')
    else:
        r.update(delta_log2CPM1=target.mean_log2CPM1.mean()-reference.mean_log2CPM1.mean(), time_equal_delta=target.groupby('time').mean_log2CPM1.mean().mean()-reference.groupby('time').mean_log2CPM1.mean().mean(), status='descriptive_unpaired')
    r.update(target_library_keys=json.dumps(target.library_key.tolist()), reference_library_keys=json.dumps(reference.library_key.tolist()), target_log2CPM1=target.mean_log2CPM1.mean(), reference_log2CPM1=reference.mean_log2CPM1.mean(), target_CPM=target.pooled_CPM.mean(), reference_CPM=reference.pooled_CPM.mean(), delta_pooled_CPM=target.pooled_CPM.mean()-reference.pooled_CPM.mean())
    return r, pair_rows


def expression_contrasts(libraries):
    out, pair_rows = [], []
    for key, raw in libraries.groupby(KEY):
        q = raw[raw.qualified & raw.mapped]
        states = {state: q[q.state.eq(state)] for state in ['P', 'M', 'S', 'N']}
        source = q[q.documented_source]
        plans = [('P_S', states['P'], source, False), ('P_M', states['P'], states['M'], True), ('M_S', states['M'], source, False), ('P_S_observed', states['P'], states['S'], False), ('M_S_observed', states['M'], states['S'], False)]
        for name, target, reference, paired in plans:
            r, pairs = comparison(key, name, target, reference, paired=paired)
            out.append(r); pair_rows.extend(pairs)
        if key[1] == 'GSE150311':
            n6 = states['N'][states['N'].library_key.str.endswith('GSM4546546')]
            n8 = states['N'][states['N'].library_key.str.endswith('GSM4546548')]
            p7 = states['P'][states['P'].library_key.str.endswith('GSM4546547')]
            for name, target in [('hN_P_primed_condition', pd.concat([n6, n8])), ('hN6_P7', n6), ('hN8_P7', n8)]:
                r, pairs = comparison(key, name, target, p7, required_target=2 if name == 'hN_P_primed_condition' else 1, required_reference=1)
                out.append(r)
            r, pairs = comparison(key, 'hN_P_same_library', states['N'], states['P'], paired=True)
            r['planned_pairs'] = 2
            r['complete_contributors'] = r['n_pairs'] == 2
            r['comparison_scope'] = 'RSeT_N8_P8_only' if r['n_pairs'] == 1 else 'two_naive_conditions' if r['n_pairs'] == 2 else 'unavailable'
            out.append(r); pair_rows.extend(pairs)
    result = pd.DataFrame(out)
    result['sign_stable_time_weighting'] = np.isfinite(result.delta_log2CPM1) & np.isfinite(result.time_equal_delta) & np.sign(result.delta_log2CPM1).eq(np.sign(result.time_equal_delta))
    result['log_and_CPM_same_sign'] = np.isfinite(result.delta_log2CPM1) & np.isfinite(result.delta_pooled_CPM) & np.sign(result.delta_log2CPM1).eq(np.sign(result.delta_pooled_CPM))
    return result, pd.DataFrame(pair_rows)


def verify_frozen_rbp(libraries):
    old = pd.read_csv(B / 'rbp_expression/library_expression.tsv.gz', sep='\t')
    old['factor'] = old.factor.str.upper()
    old = old.drop_duplicates(LIBKEY)
    factors = set(libraries.factor)
    old = old[old.factor.isin(factors) & old.state.isin(['S', 'R', 'E', 'P', 'M', 'N'])]
    merged = libraries.merge(old, on=LIBKEY, suffixes=('_new', '_frozen'), validate='one_to_one')
    assert len(merged) == len(old), 'Every already-frozen factor/library in requested states must be recovered'
    rows = []
    for column in NUMERIC:
        a = merged[column+'_new'].to_numpy(float); b = merged[column+'_frozen'].to_numpy(float)
        assert np.allclose(a, b, rtol=1e-9, atol=1e-9, equal_nan=True), column
        rows.append(dict(check='frozen_RBP_'+column, n_rows=len(merged), max_abs_error=float(np.nanmax(np.abs(a-b))), passed=True))
    return rows


def shared_roster(protocol):
    path = ROOT / 'notes/2026-10-10-naive-mouse-correspondence/primary_shared_orthologs.tsv'
    shared = pd.read_csv(path, sep='\t')
    p = shared[shared.mouse_target.eq('P_M_proximal_recurrent')]
    d = shared[shared.mouse_target.eq('P_M_PS_same_contributor_dual_proximal')]
    assert len(p) == 79 and len(d) == 46
    rows = []
    for factor in protocol['primary_factors'] + protocol['secondary_shared_factors']:
        z = p[p.h_gene_symbol.str.upper().eq(factor)]
        rows.append(dict(factor=factor, in_79_PM_overlap=bool(len(z)), in_46_dual_reference_overlap=bool(d.h_gene_symbol.str.upper().eq(factor).any()), human_gene_id=';'.join(z.human_gene_id), mouse_gene_id=';'.join(z.mouse_gene_id), hN_DWUI=float(z.h_delta.iloc[0]) if len(z) else np.nan, human_full_primary_tandem_UTR_eligible=bool(z.h_tandem.iloc[0]) if len(z) else np.nan))
    return pd.DataFrame(rows)


def add_dwui(contrasts):
    global_summary = pd.read_csv(B / 'global/global_direction_summary.tsv', sep='\t')
    q = global_summary[global_summary.setting.eq('supported50_p08') & global_summary.background.eq('available') & (global_summary.contrast.eq('P_S') & global_summary.design.eq('true_source') | global_summary.contrast.eq('P_M') & global_summary.design.eq('same_library'))]
    dwui = q[['species', 'study_id', 'route', 'contrast', 'mean', 'n']].rename(columns={'mean': 'net_DWUI_mean', 'n': 'net_DWUI_n_genes'})
    assert not dwui.duplicated(['species', 'study_id', 'route', 'contrast']).any()
    result = contrasts.merge(dwui, on=['species', 'study_id', 'route', 'contrast'], how='left', validate='many_to_one')
    naive = pd.read_csv(B / 'naive_state/library_pair_effects.tsv.gz', sep='\t')
    z = naive[naive.setting.eq('supported50_p08') & naive.node_a.eq('N') & naive.node_b.eq('P') & naive.library_b.str.endswith('GSM4546547')]
    pivot = z.pivot(index='gene_id', columns='library_a', values='delta_DWUI').dropna()
    assert pivot.shape == (2575, 2)
    result.loc[result.contrast.eq('hN_P_primed_condition'), 'net_DWUI_mean'] = pivot.mean(1).mean()
    result.loc[result.contrast.eq('hN_P_primed_condition'), 'net_DWUI_n_genes'] = len(pivot)
    coefficients = {'NUDT21': 1, 'CPSF6': 1, 'PCF11': -1, 'FIP1L1': -1}
    result['transcript_abundance_direction_coefficient'] = result.factor.map(coefficients)
    result['predicted_DWUI_direction_from_RNA'] = np.sign(result.delta_log2CPM1 * result.transcript_abundance_direction_coefficient)
    valid = np.isfinite(result.net_DWUI_mean) & np.isfinite(result.predicted_DWUI_direction_from_RNA)
    result['RNA_direction_matches_net_DWUI'] = pd.array([pd.NA] * len(result), dtype='boolean')
    result.loc[valid, 'RNA_direction_matches_net_DWUI'] = np.sign(result.net_DWUI_mean[valid]).eq(result.predicted_DWUI_direction_from_RNA[valid])
    return result


def main():
    protocol = json.loads((HERE / 'protocol.json').read_text())
    archive = np.load(HERE / 'expression_inputs.npz', allow_pickle=False)
    libraries, designs = [], []
    library_rows = legacy_library_rows()
    for study, spec in protocol['studies'].items():
        if study == 'GSE150311':
            continue
        meta = pd.DataFrame(json.loads(str(archive[study+'_meta_json'])))
        libraries.append(library_rows(meta, archive[study+'_counts'], archive[study+'_total_gex'], archive[study+'_factors'], study, spec['species']))
        designs.append(pd.DataFrame(json.loads(str(archive[study+'_design_json']))))
    naive_path = ROOT / 'notes/2026-10-10-naive-gex-recovery/factor_library_expression.tsv.gz'
    naive = pd.read_csv(naive_path, sep='\t')
    assert set(naive.factor) == set(protocol['primary_factors'] + protocol['secondary_shared_factors'])
    libraries.append(naive)
    libraries = pd.concat(libraries, ignore_index=True)
    assert not libraries.duplicated(LIBKEY).any()
    libraries = documented_source_flags(libraries, pd.concat(designs, ignore_index=True))
    libraries['factor_role'] = np.where(libraries.factor.isin(protocol['primary_factors']), 'primary_CPA', 'secondary_shared')
    states = summarize_states(libraries)
    contrasts, pairs = expression_contrasts(libraries)
    contrasts = add_dwui(contrasts)
    roster = shared_roster(protocol)
    checks = verify_frozen_rbp(libraries)
    checks.extend([dict(check='no_Liu_documented_source', n_rows=int(libraries.study_id.eq('GSE150311').sum()), max_abs_error=0, passed=not libraries[libraries.study_id.eq('GSE150311')].documented_source.any()), dict(check='no_duplicate_physical_library_state_factor', n_rows=len(libraries), max_abs_error=0, passed=not libraries.duplicated(LIBKEY).any())])
    assert all(r['passed'] for r in checks)
    mapping = pd.DataFrame(json.loads(str(archive['mapping_json'])))
    mapping = pd.concat([mapping, naive[['species', 'study_id', 'factor', 'gene_id', 'gene_symbol']].drop_duplicates().rename(columns={'gene_id': 'gene_ids', 'gene_symbol': 'symbols'}).assign(status='unique_exact_casefold_symbol')], ignore_index=True)
    ids = mapping.set_index(['species', 'study_id', 'factor'])
    index = pd.MultiIndex.from_frame(libraries[['species', 'study_id', 'factor']])
    libraries['gene_id'] = ids.gene_ids.reindex(index).to_numpy()
    libraries['gene_symbol'] = ids.symbols.reindex(index).to_numpy()
    selected = contrasts[contrasts.contrast.isin(['P_S', 'P_M', 'hN_P_primed_condition']) & contrasts.status.ne('unavailable')].copy()
    selected['comparison'] = selected.route + '::' + selected.contrast
    matrix = selected.pivot(index='factor', columns='comparison', values='delta_log2CPM1').reset_index()
    save(matrix, 'primary_comparison_matrix.tsv')
    for frame, name in [(libraries, 'library_expression.tsv.gz'), (states, 'state_expression.tsv'), (contrasts, 'expression_contrasts.tsv'), (pairs, 'paired_library_effects.tsv'), (mapping, 'gene_mapping.tsv'), (roster, 'shared_roster_audit.tsv'), (pd.DataFrame(checks), 'verification.tsv')]:
        save(frame, name)
    (HERE / 'remote_input_sha256.json').write_text(str(archive['provenance_json'])+'\n')
    inputs = [HERE / name for name in ['protocol.json', 'expression_inputs.npz', 'analyze.py', 'collect.py', 'extract_remote.py', 'test_analysis.py', 'verify.py', 'REPORT.zh.md']] + [ROOT / 'notes/2026-10-02-rbp-expression/analyze.py', ROOT / 'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json', naive_path, B / 'rbp_expression/library_expression.tsv.gz', B / 'global/global_direction_summary.tsv', B / 'naive_state/library_pair_effects.tsv.gz', ROOT / 'notes/2026-10-10-naive-mouse-correspondence/primary_shared_orthologs.tsv']
    manifest = dict(input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs}, outputs={p.name: sha(p) for p in HERE.glob('*.tsv*')}, remote_read_only=True, mapped_status=mapping.status.value_counts().to_dict(), n_library_rows=len(libraries), n_state_rows=len(states), n_contrast_rows=len(contrasts), tests_skipped=0)
    (HERE / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k: v for k, v in manifest.items() if k not in ['input_sha256', 'outputs']}, indent=2))


if __name__ == '__main__':
    main()
