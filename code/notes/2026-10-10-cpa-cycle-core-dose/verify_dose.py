"""Independent source reconstruction of the frozen exploratory WOT dose test."""
from pathlib import Path
import hashlib
import itertools
import json

import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
PAIR_SOURCE = TABLES / 'intermediate/paired_library_effects.tsv.gz'
ROSTER = ROOT / 'notes/2026-10-10-rna-module-enrichment/target_members.tsv'
CYCLE_SOURCE = ROOT / 'notes/2026-10-10-proliferation-apa/inputs.npz'
CPA_SOURCE = ROOT / 'notes/2026-10-10-cpa-factor-expression/expression_inputs.npz'
CPA_PROTOCOL = ROOT / 'notes/2026-10-10-cpa-factor-expression/protocol.json'
PAS_SOURCE = ROOT / 'notes/2026-10-10-marker-reference-sensitivity/inputs.npz'
EXPECTED_PROTOCOL = 'af5c2c2486b0d6ef032c376b3c87d3142f5e672189ec9a6caeaf01a644756291'
checks = []


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify(name, condition, **details):
    checks.append(dict(check=name, passed=bool(condition), **details))
    if not condition:
        raise ValueError(name)


def read(path):
    return pd.read_csv(path, sep='\t')


def compare(name, a, b, tolerance=1e-10):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    good = a.shape == b.shape and np.allclose(a, b, atol=tolerance, rtol=tolerance, equal_nan=True)
    error = float(np.nanmax(np.abs(a - b))) if a.size and a.shape == b.shape else None
    verify(name, good, max_abs_error=error)


def coefficient(x, y, center=False):
    x, y = np.array(x, float), np.array(y, float)
    if center:
        x -= x.mean()
        y -= y.mean()
    xx, xy, yy = x @ x, x @ y, y @ y
    return float(xy / np.sqrt(xx * yy)), float(xy / xx)


def native_counts_recovery(core, libraries, native):
    """Rebuild all 648 pass/fail rows and the 232 accepted native shifts."""
    archive = np.load(PAS_SOURCE, allow_pickle=False)
    meta = pd.DataFrame(json.loads(str(archive['GSE115943_meta_json'])))
    features = pd.DataFrame(json.loads(str(archive['GSE115943_features_json'])))
    counts = sparse.csc_matrix((archive['GSE115943_pas_data'], archive['GSE115943_pas_indices'],
                                archive['GSE115943_pas_indptr']), shape=tuple(archive['GSE115943_pas_shape']))
    verify('raw_PAS_counts_nonnegative_integer', np.issubdtype(counts.dtype, np.integer) and (counts.data >= 0).all())
    rows = []
    for gene in core.mouse_gene_id:
        full = features[features.coordinate_eligible & features.primary_gene_id.eq(gene)]
        signed = full.cut_site.to_numpy() * np.where(full.strand.eq('+'), 1, -1)
        axis = (signed - signed.min()) / np.ptp(signed)
        compare(gene + '_raw_primary_axis_recovery', full.position, axis, tolerance=1e-13)
        support = full.matched_source_count_50nt.ge(1)
        verify(gene + '_supported_sites_coordinate_unique', len(full[support]) >= 2 and full[support].cut_site.nunique() == len(full[support]))
        data = counts[:, full[support].index.to_numpy()].toarray()
        position = axis[support]
        for library in libraries:
            states = []
            for state in ['P', 'M']:
                ix = np.flatnonzero(meta.library_key.eq(library) & meta.original_selected &
                                    meta.original_node.eq(state) & meta.original_purity.ge(.8))
                pre_qualified = len(ix) >= 3 and meta.n_cells.iloc[ix].sum() >= 30
                selected = data[ix].sum(1) >= 20
                ix = ix[selected]
                cc = data[ix]
                good = pre_qualified and len(ix) >= 3 and meta.n_cells.iloc[ix].sum() >= 30
                states.append(dict(qualified=good, cc=cc, mean=float(np.mean(cc @ position / cc.sum(1))) if len(cc) else np.nan))
            a, b = states
            joint = (a['cc'].sum(0) + b['cc'].sum(0) >= 10) & ((a['cc'] > 0).sum(0) + (b['cc'] > 0).sum(0) >= 2)
            passes = a['qualified'] and b['qualified'] and joint.sum() >= 2
            rows.append(dict(gene_id=gene, library_key=library, native_pass=passes,
                             delta=a['mean'] - b['mean'] if passes else np.nan))
    raw = pd.DataFrame(rows)
    verify('raw_648_gene_library_pairs_232_pass', len(raw) == 648 and raw.native_pass.sum() == 232)
    recovered = raw[raw.native_pass].merge(native[['gene_id', 'library_key', 'delta']], on=['gene_id', 'library_key'],
                                          suffixes=('_raw', '_native'), how='outer', validate='one_to_one', indicator=True)
    verify('raw_PAS_and_frozen_native_effect_inventory_agree', recovered._merge.eq('both').all())
    compare('raw_PAS_and_frozen_native_shifts_agree', recovered.delta_raw, recovered.delta_native, tolerance=1e-12)
    return raw


def reconstruct():
    verify('dose_protocol_matches_frozen_hash', digest(HERE / 'dose_protocol.json') == EXPECTED_PROTOCOL)
    roster = read(ROSTER)
    core = roster[roster.contrast.eq('P_M') & roster.foreground_002]
    verify('fixed27_unique_membership', len(core) == 27 and core.mouse_gene_id.is_unique)
    cycle = np.load(CYCLE_SOURCE, allow_pickle=False)
    meta = pd.DataFrame(json.loads(str(cycle['GSE115943_meta_json'])))
    design = pd.DataFrame(json.loads(str(cycle['GSE115943_design_json'])))
    raw = np.log1p(cycle['GSE115943_marker_counts'] / cycle['GSE115943_total_gex'][:, None] * 1e4)
    mc_score = (raw[:, :10].mean(1) + raw[:, 10:].mean(1)) / 2
    verify('18_cycle_markers_raw_recovery', raw.shape[1] == 18 and np.allclose(mc_score, meta.Cycle, atol=1e-13, rtol=1e-13))
    meta['independent_cycle'] = mc_score
    qualified = []
    for (library, state), z in meta.groupby(['library_key', 'node']):
        if len(z) >= 3 and z.n_cells.sum() >= 30:
            qualified.append(dict(library_key=library, state=state, value=z.independent_cycle.mean()))
    cs = pd.DataFrame(qualified).pivot(index='library_key', columns='state', values='value')
    both = cs.dropna(subset=['P', 'M']).index
    verify('24_same_library_cycle_pairs', len(both) == 24)
    libs = design[design.library_key.isin(both)][['library_key', 'time', 'sample_unit', 'replicate_status']].drop_duplicates().set_index('library_key')
    verify('24_independent_well_metadata_inventory', len(libs) == 24 and libs.sample_unit.nunique() == 24 and libs.replicate_status.eq('biological_independent').all())
    libs['cycle'] = cs.P - cs.M
    libs['day'] = libs.time.str.removeprefix('D').astype(float)
    libs['C_label'] = libs.sample_unit.str.extract('(C[12])$')[0]
    verify('frozen_C_labels_resolved', libs.C_label.notna().all())
    cpa = np.load(CPA_SOURCE, allow_pickle=False)
    cm = pd.DataFrame(json.loads(str(cpa['GSE115943_meta_json'])))
    factors = json.loads(CPA_PROTOCOL.read_text())['primary_factors']
    verify('exact_11_primary_factors', len(factors) == len(set(factors)) == 11)
    fcols = np.array([np.flatnonzero(cpa['GSE115943_factors'] == factor).item() for factor in factors])
    counts = cpa['GSE115943_counts'][:, fcols]
    verify('all11_raw_CPA_mapped', (counts >= 0).all() and np.issubdtype(counts.dtype, np.integer))
    values = np.log2(1 + counts / cpa['GSE115943_total_gex'][:, None] * 1e6)
    for library in both:
        states = []
        for state in ['P', 'M']:
            ix = np.flatnonzero(cm.library_key.eq(library) & cm.node.eq(state))
            verify(library + '_' + state + '_qualified_CPA_state', len(ix) >= 3 and cm.n_cells.iloc[ix].sum() >= 30)
            states.append(values[ix].mean(0))
        libs.loc[library, 'cpa'] = (states[0] - states[1]).mean()

    chunks = []
    for z in pd.read_csv(PAIR_SOURCE, sep='\t', chunksize=250000):
        q = z[z.study_id.eq('GSE115943') & z.route.eq('WOT') & z.setting.eq('supported50_p08') & z.contrast.eq('P_M') & z.library_key.isin(both) & z.gene_id.isin(core.mouse_gene_id)]
        if len(q):
            chunks.append(q)
    native = pd.concat(chunks, ignore_index=True)
    verify('232_native_core_gene_library_effects', len(native) == 232 and not native.duplicated(['library_key', 'gene_id']).any())
    verify('native_values_finite_and_correct_time', np.isfinite(native.delta).all() and native.time.eq(native.library_key.map(libs.time)).all())
    native_counts_recovery(core, both, native)
    verify('fixed27_all24_panel_unavailable', len(set.intersection(*(set(native[native.library_key.eq(k)].gene_id) for k in both))) == 0)
    late = libs[libs.day.ge(14)].index
    verify('fixed27_D14_D18_panel_unavailable', len(late) == 18 and len(set.intersection(*(set(native[native.library_key.eq(k)].gene_id) for k in late))) == 0)
    rows = []
    for day, z in libs.groupby('day', sort=True):
        if len(z) != 2:
            rows.append(dict(time=z.time.iloc[0], day=day, n_wells=len(z), n_common_genes=0))
            continue
        c1, c2 = [z[z.C_label.eq(label)].index.item() for label in ['C1', 'C2']]
        a = native[native.library_key.eq(c1)].set_index('gene_id').delta
        b = native[native.library_key.eq(c2)].set_index('gene_id').delta
        genes = sorted(set(a.index) & set(b.index))
        rows.append(dict(time=z.time.iloc[0], day=day, n_wells=2, n_common_genes=len(genes),
                         gene_ids=';'.join(genes), library_C1=c1, library_C2=c2,
                         x_cycle=libs.loc[c1, 'cycle'] - libs.loc[c2, 'cycle'],
                         x_cpa=libs.loc[c1, 'cpa'] - libs.loc[c2, 'cpa'],
                         y_mean=float(-a.reindex(genes).mean() + b.reindex(genes).mean()) if genes else np.nan,
                         y_median=float(-a.reindex(genes).median() + b.reindex(genes).median()) if genes else np.nan))
    dates = pd.DataFrame(rows)
    primary = dates[dates.n_wells.eq(2) & dates.n_common_genes.ge(5)]
    verify('seven_primary_blocks_and_fixed_coverage_sizes', primary.time.tolist() == ['D13', 'D13.5', 'D14', 'D14.5', 'D17', 'D17.5', 'D18'] and primary.n_common_genes.tolist() == [12, 5, 22, 14, 8, 5, 9])
    stats, assignments, loo, sensitivity = [], [], [], []
    for predictor in ['cycle', 'cpa']:
        x, y = primary['x_' + predictor].to_numpy(), primary.y_mean.to_numpy()
        r, beta = coefficient(x, y)
        tails = 0
        for index, signs in enumerate(itertools.product([-1, 1], repeat=len(primary))):
            perm_r, perm_beta = coefficient(x * np.array(signs), y)
            tail = abs(perm_r) >= abs(r) - 1e-12
            tails += tail
            assignments.append(dict(predictor=predictor, assignment_index=index, signs=','.join(map(str, signs)), r=perm_r, beta=perm_beta, in_tail=bool(tail)))
        stats.append(dict(predictor=predictor, r=r, beta=beta, exact_n=128, tail_n=int(tails), p=tails / 128))
        for index, time in enumerate(primary.time):
            omit = np.arange(len(primary)) != index
            rr, bb = coefficient(x[omit], y[omit])
            loo.append(dict(predictor=predictor, omitted_time=time, r=rr, beta=bb))
        plans = [('median', primary, 'y_median', False),
                 ('all_ge1', dates[dates.n_wells.eq(2) & dates.n_common_genes.ge(1)], 'y_mean', False),
                 ('ge10', dates[dates.n_wells.eq(2) & dates.n_common_genes.ge(10)], 'y_mean', False),
                 ('batch_centered', primary, 'y_mean', True)]
        for label, q, outcome, center in plans:
            rr, bb = coefficient(q['x_' + predictor].to_numpy(), q[outcome].to_numpy(), center=center)
            sensitivity.append(dict(predictor=predictor, sensitivity=label, n_blocks=len(q), r=rr, beta=bb))
    stats = pd.DataFrame(stats)
    ordering = np.argsort(stats.p.to_numpy(), kind='stable')
    ordered_q = np.minimum.accumulate((stats.p.to_numpy()[ordering] * 2 / np.arange(1, 3))[::-1])[::-1]
    q = np.empty(2); q[ordering] = np.minimum(ordered_q, 1)
    stats['q2'] = q
    return core, native, libs.reset_index(), dates, stats, pd.DataFrame(assignments), pd.DataFrame(loo), pd.DataFrame(sensitivity)


def compare_production(core, native, libraries, dates, stats, permutations, loo, sensitivity):
    genes = read(HERE / 'library_gene_effects.tsv')
    verify('production_gene_library_matrix_is_fixed24_by27', len(genes) == 648 and genes.gene_id.nunique() == 27 and genes.library_key.nunique() == 24 and not genes.duplicated(['library_key', 'gene_id']).any())
    verify('production_fixed27_roster_agrees', set(genes.gene_id) == set(core.mouse_gene_id))
    verify('production_native_unavailable_never_zero_imputed', genes.native_pass.sum() == 232 and genes[~genes.native_pass].delta_DWUI.isna().all() and genes[~genes.native_pass].proximal_strength.isna().all())
    q = native.merge(genes, on=['library_key', 'gene_id'], suffixes=('_original', '_production'), validate='one_to_one')
    compare('production_232_native_gene_library_shifts', q.delta, q.delta_DWUI, tolerance=1e-14)
    compare('production_proximal_strength_oriented_correctly', q.delta, -q.proximal_strength, tolerance=1e-14)
    pred = read(HERE / 'paired_library_predictors.tsv').set_index('library_key').reindex(libraries.library_key)
    verify('production_24_library_predictors_complete11', len(pred) == 24 and pred.n_CPA_factors.eq(11).all())
    compare('production_cycle_from_raw_counts', pred.delta_Cycle, libraries.cycle)
    compare('production_CPA_from_raw_counts', pred.delta_CPA, libraries.cpa)
    verify('production_predictor_time_well_labels', np.array_equal(pred.time, libraries.time) and np.array_equal(pred.C_label, libraries.C_label) and np.array_equal(pred.sample_unit, libraries.sample_unit))
    d = read(HERE / 'date_effects.tsv').set_index('time').reindex(dates.time)
    verify('production_all13_date_blocks_archived', len(d) == 13 and np.array_equal(d.n_wells, dates.n_wells) and np.array_equal(d.n_shared_genes, dates.n_common_genes))
    for original, expected in [('X_Cycle', 'x_cycle'), ('X_CPA', 'x_cpa'), ('Y_mean', 'y_mean'), ('Y_median', 'y_median')]:
        compare('production_' + original + '_date_differences', d[original], dates[expected])
    for index, row in dates.iterrows():
        expected_ids = row.gene_ids.split(';') if isinstance(row.get('gene_ids'), str) and row.gene_ids else []
        verify(row.time + '_production_same_date_gene_panel', json.loads(d.loc[row.time, 'gene_ids']) == expected_ids)
        if row.n_wells != 2:
            continue
        for label in ['C1', 'C2']:
            library = row['library_' + label]
            observed = native[native.library_key.eq(library)].set_index('gene_id').delta.reindex(expected_ids)
            for variant in ['mean', 'median']:
                target = -getattr(observed, variant)() if len(observed) else np.nan
                compare(row.time + '_' + label + '_' + variant + '_matched_gene_strength', [d.loc[row.time, 'strength_' + label + '_' + variant]], [target])
    verify('production_fixed_min5_eligibility', np.array_equal(d.eligible_primary_secondary, dates.n_wells.eq(2) & dates.n_common_genes.ge(5)))
    verify('production_fixed_min1_eligibility', np.array_equal(d.eligible_sensitivity_min1, dates.n_wells.eq(2) & dates.n_common_genes.ge(1)))
    verify('production_fixed_min10_eligibility_has_three_dates', np.array_equal(d.eligible_sensitivity_min10, dates.n_wells.eq(2) & dates.n_common_genes.ge(10)) and d.eligible_sensitivity_min10.sum() == 3)
    expected_common = {(row['library_' + label], gene) for _, row in dates[dates.n_wells.eq(2)].iterrows() for label in ['C1', 'C2'] for gene in row.gene_ids.split(';') if gene}
    actual_common = set(map(tuple, genes.loc[genes.in_same_date_common_panel, ['library_key', 'gene_id']].to_numpy()))
    verify('production_within_date_panel_memberships_on_full648', expected_common == actual_common)
    s = read(HERE / 'statistics.tsv')
    s['predictor'] = s.predictor.str.lower()
    s = s.set_index('predictor').reindex(stats.predictor)
    verify('production_two_primary_statistics_no_extra_predictors', len(s) == 2 and s.n_dates.eq(7).all() and s.n_wells.eq(14).all())
    for expected, original in [('r', 'r'), ('beta', 'beta'), ('exact_n', 'exact_n'), ('tail_n', 'tail_n'), ('p', 'p_exact_two_sided'), ('q2', 'q2_BH')]:
        compare('production_primary_' + expected, s[original], stats[expected])
    for predictor in stats.predictor:
        ll = loo[loo.predictor.eq(predictor)].r
        compare('production_' + predictor + '_leave_one_date_range',
                [s.loc[predictor, 'leave_one_date_r_min'], s.loc[predictor, 'leave_one_date_r_max']], [ll.min(), ll.max()])
        verify('production_' + predictor + '_leave_one_date_sign_counts', s.loc[predictor, 'leave_one_date_positive'] == (ll > 0).sum() and s.loc[predictor, 'leave_one_date_negative'] == (ll < 0).sum() and s.loc[predictor, 'leave_one_date_zero'] == (ll == 0).sum())
    assignments = read(HERE / 'exact_assignments.tsv')
    assignments['predictor'] = assignments.predictor.str.lower()
    verify('production_all256_exact_assignments_unique', len(assignments) == 256 and not assignments.duplicated(['predictor', 'signs']).any())
    a = permutations.merge(assignments, on=['predictor', 'signs'], suffixes=('_independent', '_production'), validate='one_to_one')
    verify('production_exact_assignment_inventory_agrees', len(a) == 256)
    compare('production_all256_exact_coefficients', a.r_independent, a.r_production)
    verify('production_all256_exact_tail_decisions', np.array_equal(a.in_tail, a['tail']))
    sx = read(HERE / 'sensitivities.tsv')
    sx['predictor'] = sx.predictor.str.lower()
    verify('production_all22_sensitivity_and_leaveout_rows', len(sx) == 22)
    mapping = {'median': 'median_min5', 'all_ge1': 'mean_min1', 'ge10': 'mean_min10', 'batch_centered': 'mean_min5_batch_centered'}
    for _, row in sensitivity.iterrows():
        original = sx[sx.predictor.eq(row.predictor) & sx.analysis.eq(mapping[row.sensitivity])]
        verify(row.predictor + '_' + row.sensitivity + '_one_row_and_block_count', len(original) == 1 and original.n_dates.iloc[0] == row.n_blocks)
        compare('production_' + row.predictor + '_' + row.sensitivity + '_coefficients', original[['r', 'beta']].to_numpy(), [[row.r, row.beta]])
    for _, row in loo.iterrows():
        original = sx[sx.predictor.eq(row.predictor) & sx.analysis.eq('mean_min5_leave_one_date_out') & sx.omitted_date.eq(row.omitted_time)]
        verify(row.predictor + '_' + row.omitted_time + '_leaveout_inventory', len(original) == 1 and original.n_dates.iloc[0] == 6)
        compare('production_' + row.predictor + '_' + row.omitted_time + '_leaveout_coefficients', original[['r', 'beta']].to_numpy(), [[row.r, row.beta]])
    verify('production_sensitivities_have_no_new_p_values', not any(name.startswith('p_') or name in ['p', 'q', 'q2'] for name in sx.columns))
    return ['library_gene_effects.tsv', 'paired_library_predictors.tsv', 'date_effects.tsv', 'statistics.tsv', 'sensitivities.tsv', 'exact_assignments.tsv']


def main():
    core, native, libraries, dates, stats, permutations, loo, sensitivity = reconstruct()
    production_files = compare_production(core, native, libraries, dates, stats, permutations, loo, sensitivity)
    for label, frame in [('library_predictors', libraries), ('date_effects', dates),
                         ('statistics', stats), ('exact_assignments', permutations),
                         ('leave_one_date_out', loo), ('sensitivities', sensitivity)]:
        frame.to_csv(HERE / ('independent_' + label + '.tsv'), sep='\t', index=False, float_format='%.17g', na_rep='NA')
    verify('orientation_changes_sign_as_intended', coefficient(dates[dates.n_common_genes.ge(5)].x_cycle, dates[dates.n_common_genes.ge(5)].y_mean)[0] == -coefficient(dates[dates.n_common_genes.ge(5)].x_cycle, -dates[dates.n_common_genes.ge(5)].y_mean)[0])
    verify('batch_centering_is_explicit_sensitivity', any(abs(sensitivity[sensitivity.predictor.eq(p) & sensitivity.sensitivity.eq('batch_centered')].r.item() - stats[stats.predictor.eq(p)].r.item()) > 1e-6 for p in ['cycle', 'cpa']))
    validation = dict(status='passed', n_checks=len(checks), tests_skipped=0, checks=checks,
                      n_native_core_pairs=len(native), n_primary_time_blocks=7, exact_assignments_per_predictor=128,
                      maximum_numeric_comparison_error=max(row['max_abs_error'] for row in checks if row.get('max_abs_error') is not None),
                      compared_production_sha256={name: digest(HERE / name) for name in production_files},
                      input_sha256={str(p): digest(p) for p in [PAIR_SOURCE, ROSTER, CYCLE_SOURCE, CPA_SOURCE, CPA_PROTOCOL, PAS_SOURCE, HERE / 'dose_protocol.json']})
    (HERE / 'dose_independent_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(stats.to_string(index=False))
    print(sensitivity.to_string(index=False))


if __name__ == '__main__':
    main()
