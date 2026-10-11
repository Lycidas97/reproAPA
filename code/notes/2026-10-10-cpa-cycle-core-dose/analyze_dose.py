"""Frozen, exploratory within-date WOT comparison; wells are units."""
from pathlib import Path
import hashlib
import itertools
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / 'dose_protocol.json'
EXPECTED_PROTOCOL_SHA = 'af5c2c2486b0d6ef032c376b3c87d3142f5e672189ec9a6caeaf01a644756291'
NATIVE = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/state_vectors/pair_effects.primary.tsv.gz'
ROSTER = ROOT / 'notes/2026-10-10-rna-module-enrichment/target_members.tsv'
CYCLE = ROOT / 'notes/2026-10-10-proliferation-apa/paired_library_score_differences.tsv'
CPA = ROOT / 'notes/2026-10-10-cpa-factor-expression/paired_library_effects.tsv'
CPA_PROTOCOL = ROOT / 'notes/2026-10-10-cpa-factor-expression/protocol.json'
INPUTS = [PROTOCOL, HERE / 'PROTOCOL-CLARIFICATION.zh.md', HERE / 'coverage_protocol.json', HERE / 'coverage_presence.tsv',
          HERE / 'coverage_preflight.json', HERE / 'within_date_coverage.json',
          HERE / 'raw_coverage_audit.tsv', HERE / 'table_audit_wot_units.tsv',
          NATIVE, ROSTER, CYCLE, CPA, CPA_PROTOCOL]
checks = []


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def check(name, passed, **details):
    checks.append(dict(check=name, passed=bool(passed), **details))
    if not passed:
        raise ValueError(name)


def read(path):
    return pd.read_csv(path, sep='\t')


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, float_format='%.17g', na_rep='NA')


def balanced_panel(a, b):
    return sorted(set(a).intersection(b))


def proximal(delta):
    return -np.asarray(delta, dtype=float)


def coefficient(x, y, center=False):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if center:
        x = x - x.mean()
        y = y - y.mean()
    xx, yy, xy = x @ x, y @ y, x @ y
    if xx <= 0 or yy <= 0:
        return dict(r=np.nan, beta=np.nan, status='constant_predictor_or_outcome')
    return dict(r=float(xy / np.sqrt(xx * yy)), beta=float(xy / xx), status='available')


def exact(x, y):
    observed = coefficient(x, y)['r']
    rows = []
    for assignment, signs in enumerate(itertools.product([-1, 1], repeat=len(x))):
        value = coefficient(np.asarray(x) * signs, y)['r']
        rows.append(dict(assignment=assignment, signs=','.join(map(str, signs)),
                         r=value, tail=abs(value) >= abs(observed) - 1e-12))
    tail = sum(row['tail'] for row in rows)
    return dict(exact_n=len(rows), tail_n=tail, p_exact_two_sided=tail / len(rows)), rows


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    result = np.empty(len(p))
    result[order] = q
    return result


def intent_checks():
    x = np.array([1., 2., 4.])
    y = np.array([2., 3., 2.])
    check('intent_positive_proximal_direction', proximal([-.08, .02]).tolist() == [.08, -.02]
          and coefficient(x, proximal(-x))['r'] > .999999999999)
    shared = balanced_panel(['shared', 'only_C1'], ['shared', 'only_C2'])
    a, b = {'shared': -.1, 'only_C1': -.9}, {'shared': -.1, 'only_C2': .9}
    check('intent_unbalanced_within_date_gene_excluded', shared == ['shared']
          and np.mean([a[g] for g in shared]) - np.mean([b[g] for g in shared]) == 0
          and np.mean(list(a.values())) - np.mean(list(b.values())) != 0)
    xx = np.ravel(np.column_stack([x / 2, -x / 2]))
    yy = np.ravel(np.column_stack([y / 2, -y / 2]))
    through = coefficient(x, y)['r']
    centered = coefficient(x, y, center=True)['r']
    check('intent_date_fixed_centering_is_not_oriented_difference_centering',
          np.isclose(through, np.corrcoef(xx, yy)[0, 1]) and not np.isclose(through, centered))
    result, assignments = exact(np.arange(1., 8.), np.arange(1., 8.))
    check('intent_exact_assignment_inventory', result['exact_n'] == len(assignments) == 128
          and result['tail_n'] == 2 and result['p_exact_two_sided'] == 2 / 128)


def load_native(core, units):
    columns = ['study_id', 'route', 'library_key_a', 'library_key_b', 'time_a', 'time_b',
               'sample_unit_a', 'sample_unit_b', 'contrast', 'design', 'gene_id',
               'gene_symbol', 'setting', 'delta_DWUI', 'n_mc_a', 'n_mc_b', 'n_cells_a',
               'n_cells_b', 'n_joint_PAS']
    pieces = []
    for q in pd.read_csv(NATIVE, sep='\t', usecols=columns, chunksize=10000):
        selected = q.study_id.eq('GSE115943') & q.route.eq('WOT') & q.contrast.eq('P_M') \
            & q.design.eq('same_library') & q.setting.eq('supported50_p08') \
            & q.gene_id.isin(core.mouse_gene_id) & q.library_key_a.isin(units.library_key)
        pieces.append(q[selected])
    native = pd.concat(pieces, ignore_index=True)
    check('native_gene_library_unique', not native.duplicated(['gene_id', 'library_key_a']).any())
    check('native_physical_library_and_time_agree', native.library_key_a.eq(native.library_key_b).all()
          and native.time_a.eq(native.time_b).all() and native.sample_unit_a.eq(native.sample_unit_b).all())
    check('native_finite_qualified_effects', np.isfinite(native.delta_DWUI).all()
          and native[['n_mc_a', 'n_mc_b']].ge(3).all().all()
          and native[['n_cells_a', 'n_cells_b']].ge(30).all().all()
          and native.n_joint_PAS.ge(2).all())
    presence = read(HERE / 'coverage_presence.tsv')
    key = ['gene_id', 'library_key_a']
    expected = set(map(tuple, presence[key].to_numpy()))
    observed = set(map(tuple, native[key].to_numpy()))
    check('native_outcome_inventory_equals_blind_presence', observed == expected,
          n_measured_gene_library=len(observed))
    native = native.rename(columns={'library_key_a': 'library_key', 'time_a': 'time',
                                    'sample_unit_a': 'native_sample_unit'})
    audit = read(HERE / 'raw_coverage_audit.tsv')
    check('all27_all24_raw_coverage_inventory', len(audit) == 27 * 24
          and not audit.duplicated(['gene_id', 'library_key']).any())
    q = audit.merge(native[['gene_id', 'library_key', 'delta_DWUI', 'n_mc_a', 'n_mc_b',
                           'n_cells_a', 'n_cells_b', 'n_joint_PAS', 'native_sample_unit']],
                    on=['gene_id', 'library_key'], how='left', validate='one_to_one')
    check('native_presence_equals_original_raw_gate', q.native_pass.eq(q.delta_DWUI.notna()).all())
    q = q.merge(units[['library_key', 'time', 'day', 'C_label', 'sample_unit']],
                on='library_key', validate='many_to_one')
    check('native_and_audited_sample_unit_agree', q[q.native_pass].native_sample_unit.eq(q[q.native_pass].sample_unit).all())
    q['proximal_strength'] = proximal(q.delta_DWUI)
    q['status'] = np.where(q.native_pass, 'measured_native_primary', q.native_failure_reason)
    return q


def load_predictors(units):
    cycle = read(CYCLE)
    cycle = cycle[cycle.study_id.eq('GSE115943') & cycle.route.eq('WOT')
                  & cycle.contrast.eq('P_M') & cycle.score.eq('Cycle')]
    check('cycle_24_unique_physical_pairs', len(cycle) == 24 and cycle.library_key.is_unique
          and set(cycle.library_key) == set(units.library_key))
    factors = json.loads(CPA_PROTOCOL.read_text())['primary_factors']
    check('exactly11_fixed_CPA_factors', len(factors) == len(set(factors)) == 11)
    cpa = read(CPA)
    cpa = cpa[cpa.study_id.eq('GSE115943') & cpa.route.eq('WOT')
              & cpa.contrast.eq('P_M') & cpa.factor.isin(factors)]
    check('all11_all24_CPA_pairs_no_selection', len(cpa) == 11 * 24
          and not cpa.duplicated(['library_key', 'factor']).any()
          and cpa.groupby('library_key').factor.nunique().eq(11).all()
          and set(cpa.library_key) == set(units.library_key)
          and np.isfinite(cpa.delta_log2CPM1).all())
    q = units[['library_key', 'time', 'day', 'C_label', 'sample_unit']].merge(
        cycle[['library_key', 'time', 'delta']], on=['library_key', 'time'], validate='one_to_one')
    q = q.rename(columns={'delta': 'delta_Cycle'})
    cpa_rows = cpa.groupby(['library_key', 'time']).agg(
        delta_CPA=('delta_log2CPM1', 'mean'), n_CPA_factors=('factor', 'nunique')).reset_index()
    q = q.merge(cpa_rows, on=['library_key', 'time'], validate='one_to_one')
    check('predictor24_time_and_unit_preserved', len(q) == 24 and np.isfinite(q[['delta_Cycle', 'delta_CPA']]).all().all())
    return q


def date_summary(effects, predictors):
    rows = []
    effects['in_same_date_common_panel'] = False
    for time, wells in predictors.groupby('time', sort=False):
        time_effects = effects[effects.time.eq(time)]
        row = dict(time=time, day=float(wells.day.iloc[0]), n_wells=len(wells),
                   library_C1='NA', library_C2='NA', n_shared_genes=0, gene_ids='[]',
                   strength_C1_mean=np.nan, strength_C2_mean=np.nan,
                   strength_C1_median=np.nan, strength_C2_median=np.nan,
                   Y_mean=np.nan, Y_median=np.nan, X_Cycle=np.nan, X_CPA=np.nan)
        if len(wells) != 2:
            row['status'] = 'excluded_singleton_date'
            rows.append(row)
            continue
        check(time + '_C1_C2_exactly_one_each', sorted(wells.C_label.tolist()) == ['C1', 'C2'])
        well = wells.set_index('C_label')
        one, two = well.loc['C1'], well.loc['C2']
        a = time_effects[time_effects.library_key.eq(one.library_key) & time_effects.native_pass]
        b = time_effects[time_effects.library_key.eq(two.library_key) & time_effects.native_pass]
        shared = balanced_panel(a.gene_id, b.gene_id)
        row.update(library_C1=one.library_key, library_C2=two.library_key,
                   n_shared_genes=len(shared), gene_ids=json.dumps(shared),
                   X_Cycle=float(one.delta_Cycle - two.delta_Cycle),
                   X_CPA=float(one.delta_CPA - two.delta_CPA))
        if shared:
            aa = a[a.gene_id.isin(shared)].sort_values('gene_id')
            bb = b[b.gene_id.isin(shared)].sort_values('gene_id')
            check(time + '_balanced_equal_gene_panels', aa.gene_id.tolist() == bb.gene_id.tolist() == shared)
            row.update(strength_C1_mean=float(aa.proximal_strength.mean()),
                       strength_C2_mean=float(bb.proximal_strength.mean()),
                       strength_C1_median=float(aa.proximal_strength.median()),
                       strength_C2_median=float(bb.proximal_strength.median()), status='balanced_observed_panel')
            row['Y_mean'] = row['strength_C1_mean'] - row['strength_C2_mean']
            row['Y_median'] = row['strength_C1_median'] - row['strength_C2_median']
            effects.loc[effects.time.eq(time) & effects.gene_id.isin(shared), 'in_same_date_common_panel'] = True
        else:
            row['status'] = 'excluded_no_common_measured_gene'
        rows.append(row)
    dates = pd.DataFrame(rows).sort_values('day').reset_index(drop=True)
    dates['eligible_primary_secondary'] = dates.n_shared_genes.ge(5) & dates.n_wells.eq(2)
    dates['eligible_sensitivity_min1'] = dates.n_shared_genes.ge(1) & dates.n_wells.eq(2)
    dates['eligible_sensitivity_min10'] = dates.n_shared_genes.ge(10) & dates.n_wells.eq(2)
    check('13dates7primary10min1three_min10_rule', len(dates) == 13
          and int(dates.eligible_primary_secondary.sum()) == 7
          and int(dates.eligible_sensitivity_min1.sum()) == 10
          and int(dates.eligible_sensitivity_min10.sum()) == 3)
    frozen = json.loads((HERE / 'within_date_coverage.json').read_text())
    for f in frozen:
        if f['n_wells'] == 2:
            d = dates[dates.time.eq(f['time'])].iloc[0]
            check(f['time'] + '_outcome_blind_panel_preserved', d.n_shared_genes == f['n_shared_genes']
                  and json.loads(d.gene_ids) == f['gene_ids'])
    return effects, dates


def analyze(dates):
    primary = dates[dates.eligible_primary_secondary]
    stats, exact_rows, sensitivity_rows = [], [], []
    for predictor in ['Cycle', 'CPA']:
        x = primary['X_' + predictor].to_numpy()
        y = primary.Y_mean.to_numpy()
        result = dict(analysis='secondary_exploratory_frozen_min5', predictor=predictor,
                      n_dates=len(primary), n_wells=2 * len(primary),
                      minimum_genes=int(primary.n_shared_genes.min()),
                      maximum_genes=int(primary.n_shared_genes.max()),
                      coefficient_definition='date_fixed_through_origin_C1_minus_C2', **coefficient(x, y))
        check(predictor + '_primary_nonconstant', result['status'] == 'available')
        exact_result, assignments = exact(x, y)
        check(predictor + '_exact128', exact_result['exact_n'] == 128)
        result.update(exact_result)
        leave = [coefficient(np.delete(x, i), np.delete(y, i))['r'] for i in range(len(x))]
        result.update(leave_one_date_r_min=float(min(leave)), leave_one_date_r_max=float(max(leave)),
                      leave_one_date_positive=int(np.sum(np.array(leave) > 0)),
                      leave_one_date_negative=int(np.sum(np.array(leave) < 0)),
                      leave_one_date_zero=int(np.sum(np.array(leave) == 0)))
        stats.append(result)
        exact_rows.extend(dict(predictor=predictor, **a) for a in assignments)
        variants = [('median_min5', primary, 'Y_median', False),
                    ('mean_min1', dates[dates.eligible_sensitivity_min1], 'Y_mean', False),
                    ('mean_min10', dates[dates.eligible_sensitivity_min10], 'Y_mean', False),
                    ('mean_min5_batch_centered', primary, 'Y_mean', True)]
        for name, q, outcome, center in variants:
            r = coefficient(q['X_' + predictor], q[outcome], center=center)
            sensitivity_rows.append(dict(analysis=name, predictor=predictor,
                n_dates=len(q), n_wells=2 * len(q), minimum_genes=int(q.n_shared_genes.min()),
                maximum_genes=int(q.n_shared_genes.max()), outcome=outcome,
                centered_oriented_differences=center, **r))
        for i, row in enumerate(primary.itertuples(index=False)):
            sensitivity_rows.append(dict(analysis='mean_min5_leave_one_date_out', predictor=predictor,
                omitted_date=row.time, n_dates=6, n_wells=12,
                outcome='Y_mean', centered_oriented_differences=False,
                **coefficient(np.delete(x, i), np.delete(y, i))))
    frame = pd.DataFrame(stats)
    frame['q2_BH'] = bh(frame.p_exact_two_sided)
    return frame, pd.DataFrame(sensitivity_rows), pd.DataFrame(exact_rows), dict(
        date_fixed_r=coefficient(primary.X_Cycle, primary.X_CPA)['r'],
        batch_centered_r=coefficient(primary.X_Cycle, primary.X_CPA, center=True)['r'])


def main():
    check('frozen_protocol_sha', sha(PROTOCOL) == EXPECTED_PROTOCOL_SHA)
    intent_checks()
    core = read(ROSTER)
    core = core[core.contrast.eq('P_M') & core.foreground_002]
    check('exactly27_fixed_unique_mouse_core', len(core) == 27 and core.mouse_gene_id.is_unique)
    units = read(HERE / 'table_audit_wot_units.tsv')
    check('all24_unique_audited_well_units', len(units) == units.library_key.nunique()
          == units.sample_unit.nunique() == 24)
    effects = load_native(core, units)
    predictors = load_predictors(units)
    effects, dates = date_summary(effects, predictors)
    stats, sensitivities, assignments, co_variation = analyze(dates)
    save(effects.sort_values(['day', 'C_label', 'gene_id']), 'library_gene_effects.tsv')
    save(predictors.sort_values(['day', 'C_label']), 'paired_library_predictors.tsv')
    save(dates, 'date_effects.tsv')
    save(stats, 'statistics.tsv')
    save(sensitivities, 'sensitivities.tsv')
    save(assignments, 'exact_assignments.tsv')
    exclusions = dates.loc[~dates.eligible_primary_secondary,
                           ['time', 'n_wells', 'n_shared_genes', 'status']].to_dict('records')
    validation = dict(status='passed', independent_validation='pending', checks=checks,
        n_checks=len(checks), tests_skipped=0, protocol_sha256=sha(PROTOCOL),
        primary_requested_fixed27_all24_status='unavailable_no_common_measured_gene',
        secondary_status='completed_exploratory_seven_dates_14_wells_varying_balanced_panels',
        native_gene_library_missing=int((~effects.native_pass).sum()),
        native_gene_library_measured=int(effects.native_pass.sum()),
        native_missing_reasons=effects.loc[~effects.native_pass, 'native_failure_reason'].value_counts().to_dict(),
        secondary_date_exclusions=exclusions, cycle_CPA_co_variation=co_variation,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in INPUTS})
    (HERE / 'analysis_validation.json').write_text(json.dumps(validation, indent=2, ensure_ascii=False) + '\n')
    produced = ['analyze_dose.py', 'library_gene_effects.tsv', 'paired_library_predictors.tsv',
                'date_effects.tsv', 'statistics.tsv', 'sensitivities.tsv',
                'exact_assignments.tsv', 'analysis_validation.json', 'ANALYSIS.zh.md']
    manifest = dict(analysis='WOT within-date exploratory core dose analysis',
        status='production_passed_independent_QA_separate', protocol_sha256=sha(PROTOCOL),
        software=dict(numpy=np.__version__, pandas=pd.__version__),
        input_sha256=validation['input_sha256'],
        output_sha256={name: sha(HERE / name) for name in produced},
        output_row_counts=dict(library_gene_effects=len(effects), paired_library_predictors=len(predictors),
                               date_effects=len(dates), statistics=len(stats),
                               sensitivities=len(sensitivities), exact_assignments=len(assignments)),
        no_manuscript_or_reference_changes=True)
    (HERE / 'production_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')
    print(stats.to_string(index=False))
    print(sensitivities[sensitivities.analysis.ne('mean_min5_leave_one_date_out')].to_string(index=False))
    print(json.dumps(dict(cycle_CPA_co_variation=co_variation,
                         native_measured=validation['native_gene_library_measured'],
                         native_missing=validation['native_gene_library_missing'],
                         n_checks=len(checks), tests_skipped=0), ensure_ascii=False))


if __name__ == '__main__':
    main()
