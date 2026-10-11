"""Expression/depth sensitivity on the frozen Figure 4c ortholog labels."""
import argparse
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, hypergeom, nchypergeom_fisher

from downstream_runtime import Context, Study, bh, seed
from downstream_cross import CODES, MOUSE, feature_frame, matrices

CONTROL = 'notes/2026-10-02-null-model-controls/'
BASES = {
    'Q5': ['effect_n', 'actual_effect_n', 'prox_capacity', 'ref_q5', 'pas_q3'],
    'Range': ['effect_n', 'actual_effect_n', 'prox_capacity', 'pas_q3', 'ref_min10', 'ref_max10'],
}


def paired_weights(frame):
    """Match Figure 4c: equal dates, then units within date, then libraries."""
    assert frame[['time', 'sample_unit']].notna().all().all()
    dates = frame.time.nunique()
    units = frame.groupby('time').sample_unit.transform('nunique')
    libraries = frame.groupby(['time', 'sample_unit']).time.transform('size')
    weights = 1. / (dates * units.to_numpy() * libraries.to_numpy())
    assert abs(weights.sum() - 1) < 1e-12
    return weights


def summarize_contributors(frame, contrast, design=None):
    """Keep the different observed-node and paired-library estimands intact."""
    columns = ['position', 'GEX_CPM', 'APA_coverage', 'n_mc']
    if contrast == 'PS':
        selected = frame[frame.node.isin(['P', 'S'])].copy()
        assert set(selected.node) == {'P', 'S'}
        counts = selected.groupby('node').node.transform('size')
        selected['weight'] = .5 / counts
        target = selected[selected.node.eq('P')][columns].mean(skipna=False)
        reference = selected[selected.node.eq('S')][columns].mean(skipna=False)
    else:
        assert contrast == 'PM'
        wide = frame.pivot(index='library_key', columns='node', values='position')
        eligible = wide.reindex(columns=['P', 'M']).dropna().index
        pairs = pd.DataFrame({'library_key': eligible}).merge(
            design, on='library_key', validate='one_to_one')
        assert len(pairs)
        pairs['pair_weight'] = paired_weights(pairs)
        selected = frame[frame.library_key.isin(eligible) & frame.node.isin(['P', 'M'])].copy()
        selected = selected.merge(pairs, on='library_key', validate='many_to_one')
        selected['weight'] = selected.pair_weight / 2
        target = selected[selected.node.eq('P')].set_index('library_key').loc[eligible]
        reference = selected[selected.node.eq('M')].set_index('library_key').loc[eligible]
        weights = pairs.set_index('library_key').pair_weight.loc[eligible].to_numpy()
        target = pd.Series(weights @ target[columns].to_numpy(), index=columns)
        reference = pd.Series(weights @ reference[columns].to_numpy(), index=columns)
    assert abs(selected.weight.sum() - 1) < 1e-12
    record = dict(
        effect=float(target.position - reference.position),
        mean_GEX=float((target.GEX_CPM + reference.GEX_CPM) / 2),
        APA_coverage=float((target.APA_coverage + reference.APA_coverage) / 2),
        target_GEX=float(target.GEX_CPM), reference_GEX=float(reference.GEX_CPM),
        target_coverage=float(target.APA_coverage), reference_coverage=float(reference.APA_coverage),
        min_endpoint_GEX=float(np.minimum(target.GEX_CPM, reference.GEX_CPM)),
        min_endpoint_coverage=float(min(target.APA_coverage, reference.APA_coverage)),
        min_library_coverage=float(selected.APA_coverage.min()),
        min_library_metacells=int(selected.n_mc.min()),
        contributor_states=len(selected),
    )
    return record, selected


def extract(ctx):
    source = ctx.analysis / 'cross_overlap'
    frozen = ctx.read(source / 'ortholog_matched_features.tsv.gz')
    assert len(frozen) == 16092 and frozen.human_gene_id.is_unique and frozen.mouse_gene_id.is_unique
    nodes = ctx.read(source / 'all_mapped_nodes.tsv.gz')
    paired = ctx.read(source / 'paired_pm_study_effects.tsv.gz')
    replay = feature_frame(frozen[['human_gene_id', 'mouse_gene_id']], nodes, paired, 'supported50')
    for column in replay.columns.intersection(frozen.columns):
        pd.testing.assert_series_equal(replay[column], frozen[column], check_dtype=False,
                                       check_exact=False, rtol=0, atol=1e-12)
    effects, _ = matrices(nodes, paired, 'supported50')
    libraries = ctx.read(source / 'all_mapped_library_nodes.tsv.gz', usecols=[
        'code', 'layer', 'gene_id', 'library_key', 'node', 'position_index',
        'qualified_metacells', 'qualified_cells', 'pas_counts', 'PAS_number'])
    libraries = libraries[libraries.layer.eq('supported50') & libraries.node.isin(['P', 'S', 'M'])]
    all_covariates, checks = [], []
    for code, (sid, route) in CODES.items():
        print('EXTRACT_START', code, sid, flush=True)
        study = Study(ctx, sid)
        study.load_gex()
        id_column = 'mouse_gene_id' if code in MOUSE else 'human_gene_id'
        orthologs = set(frozen[id_column])
        eligible = {co: set(effects[co][code].dropna().index) & orthologs for co in ['PS', 'PM']}
        needed = sorted(eligible['PS'] | eligible['PM'])
        lib = libraries[libraries.code.eq(code)].set_index('gene_id')
        design = study.design[study.design.route.eq(route)][['library_key', 'time', 'sample_unit']].drop_duplicates()
        assert design.library_key.is_unique
        covariates, contributors = [], []
        maximum_position_error = 0.
        for number, gene in enumerate(needed):
            axis = study.axis(gene, 'supported50_p08')
            data = study.gene_counts(axis)
            expected = lib.loc[[gene]].reset_index()
            states = []
            for row in expected.itertuples(index=False):
                state = study.state(data, route, row.library_key, row.node)
                assert state['profile_eligible']
                assert state['n_mc'] == row.qualified_metacells
                assert state['n_cells'] == row.qualified_cells
                assert int(state['site_counts'].sum()) == row.pas_counts
                assert len(axis) == row.PAS_number
                position = float(state['usage'] @ axis.WPP.to_numpy())
                error = abs(position - row.position_index)
                assert error < 1e-12, (code, gene, row.library_key, row.node, error)
                maximum_position_error = max(maximum_position_error, error)
                expression = study.expression_summary(gene, state['rows'])
                states.append(dict(
                    library_key=row.library_key, node=row.node, position=position,
                    n_mc=state['n_mc'], n_cells=state['n_cells'], pas_counts=row.pas_counts,
                    APA_coverage=row.pas_counts / state['n_mc'], **expression,
                    metacell_ids=json.dumps(study.meta.metacell_id.iloc[state['rows']].tolist()),
                ))
            frame = pd.DataFrame(states)
            for co in ['PS', 'PM']:
                if gene not in eligible[co]:
                    continue
                record, contribution = summarize_contributors(frame, co, design)
                assert abs(record['effect'] - effects[co].loc[gene, code]) < 1e-12
                context = dict(code=code, study_id=sid, route=route, gene_id=gene, contrast=co)
                covariates.append(context | record)
                contributors.append(contribution.assign(**context))
            if (number + 1) % 1000 == 0:
                print('EXTRACT_PROGRESS', code, number + 1, len(needed), flush=True)
        covariates = pd.DataFrame(covariates)
        ctx.table(covariates, f'{code}-study-covariates.tsv.gz')
        ctx.table(pd.concat(contributors, ignore_index=True), f'{code}-contributors.tsv.gz')
        all_covariates.append(covariates)
        checks.append(dict(code=code, genes=len(needed), PS=len(eligible['PS']), PM=len(eligible['PM']),
                           maximum_position_error=maximum_position_error,
                           missing_expression=int(covariates.mean_GEX.isna().sum())))
        print('EXTRACT_DONE', checks[-1], flush=True)
        del study, data, frame, lib, expected, contributors
        gc.collect()
    covariates = pd.concat(all_covariates, ignore_index=True)
    assert not covariates.duplicated(['code', 'gene_id', 'contrast']).any()
    ctx.table(covariates, 'study-covariates.tsv.gz')
    joined = frozen.copy()
    for side, codes, id_column in [
        ('h', ['H1', 'H2', 'H3', 'H4'], 'human_gene_id'),
        ('m', ['W', 'C'], 'mouse_gene_id'),
    ]:
        for co in ['PS', 'PM']:
            part = covariates[covariates.code.isin(codes) & covariates.contrast.eq(co)]
            grouped = part.groupby('gene_id')
            for column in ['mean_GEX', 'APA_coverage', 'min_endpoint_GEX',
                           'min_endpoint_coverage', 'min_library_coverage', 'min_library_metacells']:
                # Any unavailable contributing expression makes the species covariate NA.
                value = grouped[column].agg(lambda x: np.nan if x.isna().any() else x.median())
                joined[f'{side}_{co}_{column}'] = joined[id_column].map(value)
            joined[f'{side}_{co}_actual_effect_n'] = joined[id_column].map(grouped.size()).fillna(0).astype(int)
            joined[f'{side}_{co}_study_pattern'] = joined[id_column].map(
                grouped.code.agg(lambda x: ';'.join(sorted(x)))).fillna('unmeasured')
            joined[f'{side}_{co}_GEX_missing_studies'] = joined[id_column].map(
                grouped.mean_GEX.agg(lambda x: int(x.isna().sum()))).fillna(0).astype(int)
            measurable = joined[f'{side}_{co}_actual_effect_n'].gt(0)
            assert joined.loc[measurable, f'{side}_{co}_APA_coverage'].notna().all()
            assert joined.loc[joined[f'{side}_{co}'], f'{side}_{co}_mean_GEX'].notna().all()
    ctx.table(joined, 'ortholog-expression-coverage.tsv.gz')
    ctx.table(checks, 'extraction-checks.tsv')
    ctx.finish(dict(frozen_features_replayed=True, full_ortholog_background_preserved=True,
                    exact_gene_qualified_metacells=True, library_counts_and_positions_recovered=True,
                    observed_node_and_paired_effects_recovered=True, no_covariate_complete_case_deletion=True),
               aggregation='median of actual measurable studies; any contributing GEX NA makes aggregate NA')


def terciles(values):
    result = pd.Series(-1, index=values.index, dtype=int)
    valid = np.isfinite(values)
    if not valid.any():
        return result, []
    if values[valid].nunique() == 1:
        result.loc[valid] = 0
        return result, [float(values[valid].iloc[0])]
    codes, edges = pd.qcut(values[valid], 3, labels=False, retbins=True, duplicates='drop')
    result.loc[valid] = codes.fillna(0).astype(int)
    return result, edges.tolist()


def balance(values, labels, probabilities):
    valid = np.isfinite(values)
    observed = labels & valid
    weights = probabilities[valid]
    selected = values[observed]
    background = values[valid]
    if not len(selected) or not weights.sum():
        return dict(observed_mean=np.nan, null_mean=np.nan, standardized_difference=np.nan)
    mu = float(np.average(background, weights=weights))
    variance = float(np.average((background - mu) ** 2, weights=weights))
    scale = np.sqrt((selected.var() + variance) / 2)
    return dict(observed_mean=float(selected.mean()), null_mean=mu,
                standardized_difference=float((selected.mean() - mu) / scale) if scale else 0.,
                observed_missing=int((labels & ~valid).sum()), expected_missing=float(probabilities[~valid].sum()))


def infer(ctx, extraction):
    f = ctx.read(Path(extraction) / 'ortholog-expression-coverage.tsv.gz')
    original = ctx.read(ctx.analysis / 'cross_overlap/overlap_null_results.tsv')
    common = ctx.functions(CONTROL + 'common.py', ['groups', 'binary_draw', 'mobility', 'probability', 'rate'], {'beta': beta})
    numerical = ctx.functions(CONTROL + 'cross_species.py', ['independent', 'injected'],
                              dict(binary_draw=common.binary_draw, hypergeom=hypergeom,
                                   nchypergeom_fisher=nchypergeom_fisher))
    rows, strata, balances, powers, bin_rows, missing, replay = [], [], [], [], [], [], []
    archive = {}
    for co in ['PS', 'PM']:
        x, y = f['h_' + co].to_numpy(bool), f['m_' + co].to_numpy(bool)
        observed = int((x & y).sum())
        source_row = original[original.contrast.eq(co)].iloc[0]
        assert len(f) == 16092 and observed == source_row.observed
        for side, labels in [('h', x), ('m', y)]:
            assert not (labels & f[f'{side}_{co}_actual_effect_n'].lt(2).to_numpy()).any()
        for side, labels in [('h', x), ('m', y)]:
            for column in ['mean_GEX', 'APA_coverage']:
                values = f[f'{side}_{co}_{column}']
                codes, edges = terciles(values)
                f[f'{side}_{co}_{column}_q3'] = codes
                bin_rows.append(dict(contrast=co, side=side, covariate=column,
                                     edges=json.dumps(edges), counts=json.dumps(codes.value_counts().sort_index().to_dict())))
            measured = f[f'{side}_{co}_actual_effect_n'].gt(0).to_numpy()
            missing.append(dict(contrast=co, side=side, actual_measurable_genes=int(measured.sum()),
                                candidates=int(labels.sum()),
                                measured_missing_GEX=int((measured & f[f'{side}_{co}_mean_GEX'].isna()).sum()),
                                measured_missing_coverage=int((measured & f[f'{side}_{co}_APA_coverage'].isna()).sum()),
                                candidate_missing_GEX=int((labels & f[f'{side}_{co}_mean_GEX'].isna()).sum())))
        for base, fields in BASES.items():
            model = base + '_expression_coverage'
            groups = {}
            probabilities = {}
            for side, labels in [('h', x), ('m', y)]:
                old = common.groups(f, [f'{side}_{co}_{column}' for column in fields if column != 'actual_effect_n'])
                old_probabilities, _ = common.mobility(labels, old)
                probabilities[side + '_old'] = old_probabilities
                columns = fields + ['mean_GEX_q3', 'APA_coverage_q3']
                cells = common.groups(f, [f'{side}_{co}_{column}' for column in columns])
                groups[side] = cells
                prob, diagnostic = common.mobility(labels, cells)
                probabilities[side] = prob
                row_diagnostic = {side + '_' + key: value for key, value in diagnostic.items()}
                probabilities[side + '_diagnostic'] = row_diagnostic
                # Actual random draws verify every stratum margin, including fixed zeros/ones.
                sample = common.binary_draw(labels, cells, 17, np.random.default_rng(seed('margin_check', co, model, side)))
                for number, ids in enumerate(cells):
                    events = int(labels[ids].sum())
                    assert np.all(sample[:, ids].sum(1) == events)
                    strata.append(dict(contrast=co, model=model, side=side, stratum=number,
                                       genes=len(ids), events=events, mixed=0 < events < len(ids),
                                       keys=json.dumps([int(f.iloc[ids[0]][f'{side}_{co}_{column}']) for column in columns])))
                ineligible = f[f'{side}_{co}_actual_effect_n'].lt(2).to_numpy()
                assert not sample[:, ineligible].any()
                assert not prob[ineligible].any()
                for column in ['mean_GEX', 'APA_coverage', 'min_endpoint_GEX', 'min_endpoint_coverage', 'min_library_coverage']:
                    values = f[f'{side}_{co}_{column}'].to_numpy(float)
                    for scale in ['raw', 'log2plus1']:
                        transformed = values if scale == 'raw' else np.log2(1 + values)
                        for matching, weights in [('original', old_probabilities), ('expression_coverage', prob)]:
                            balances.append(dict(contrast=co, model=model, side=side, covariate=column,
                                                 scale=scale, matching=matching, **balance(transformed, labels, weights)))
            old_mean = float(probabilities['h_old'] @ probabilities['m_old'])
            source_model = 'independent_q5_reference_capacity' if base == 'Q5' else 'independent_range10_reference_capacity'
            source = original[original.contrast.eq(co) & original.model.eq(source_model)].iloc[0]
            assert abs(old_mean - source.analytic_null_mean) < 1e-10
            replay.append(dict(contrast=co, model=source_model, original_expected=source.analytic_null_mean,
                               replay_expected=old_mean, observed=observed))
            expected = float(probabilities['h'] @ probabilities['m'])
            draw_seed = seed('cross_expression_coverage', co, model)
            null = numerical.independent(x, y, groups['h'], groups['m'], 20000, np.random.default_rng(draw_seed))
            assert abs(null.mean() - expected) <= max(.01, 6 * null.std(ddof=1) / np.sqrt(len(null)))
            archive[co + '|' + model] = null
            fixed = int(((probabilities['h'] == 1) & (probabilities['m'] == 1)).sum())
            record = dict(contrast=co, model=model, observed=observed, orthologs=len(f),
                          human_events=int(x.sum()), mouse_events=int(y.sum()), analytic_expected=expected,
                          null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                          null_q025=float(np.quantile(null, .025)), null_q975=float(np.quantile(null, .975)),
                          null_unique=len(np.unique(null)), fixed_overlap=fixed,
                          enrichment=observed / expected if expected else np.nan,
                          p=common.probability(null, observed), seed=draw_seed, permutations=len(null),
                          **probabilities['h_diagnostic'], **probabilities['m_diagnostic'])
            rows.append(record)
            print('NULL_RESULT', record, flush=True)
            for odds in [1, 2, 4]:
                injection_seed = seed('cross_expression_coverage_injection', co, model, odds)
                simulated = numerical.injected(x, y, groups['h'], groups['m'], 3000, odds,
                                               np.random.default_rng(injection_seed))
                archive[co + '|' + model + '|OR' + str(odds)] = simulated
                for alpha in [.05, .025, .0125]:
                    powers.append(dict(contrast=co, model=model, odds_ratio=odds, alpha=alpha,
                                       seed=injection_seed, simulated_mean=float(simulated.mean()),
                                       **common.rate(common.probability(null, simulated) <= alpha)))
                print('CALIBRATION_POWER', co, model, odds, powers[-3], flush=True)
    summary = pd.DataFrame(rows)
    summary['q_two_contrasts'] = summary.groupby('model').p.transform(bh)
    summary['q_four_tests'] = bh(summary.p)
    ctx.table(summary, 'expression-coverage-null-results.tsv')
    ctx.table(f, 'ortholog-expression-coverage-binned.tsv.gz')
    ctx.table(strata, 'strata.tsv.gz')
    ctx.table(balances, 'balance.tsv')
    ctx.table(powers, 'calibration-power.tsv')
    ctx.table(bin_rows, 'bin-definitions.tsv')
    ctx.table(missing, 'coverage-and-missingness.tsv')
    ctx.table(replay, 'original-model-replay.tsv')
    np.savez_compressed(ctx.out / 'null-draws.npz', **archive)
    ctx.finish(dict(full_ortholog_background_preserved=True, original_labels_and_counts_preserved=True,
                    original_model_expectations_replayed=True, actual_permutation_margins_preserved=True,
                    no_events_outside_candidate_measurement_eligibility=True,
                    analytic_expected_agrees_with_monte_carlo=True, prespecified_four_test_family=True,
                    self_calibration_and_conditional_power=True),
               permutations=20000, injection_datasets=3000, odds_ratios=[1, 2, 4],
               missing_expression='explicit missing bin; no background deletion; actual zero stays zero')
    print(summary.to_string(index=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--stage', choices=['extract', 'infer'], required=True)
    parser.add_argument('--extraction')
    args = parser.parse_args()
    context = Context(args.config, 'cross_overlap_expression_coverage_' + args.stage, args.out)
    for name in ['cross_overlap_expression_coverage.py', 'downstream_runtime.py',
                 'downstream_cross.py', 'polyapipe_prepare.py', 'polyapipe_effects.py']:
        path = context.register(Path(__file__).parent / name)
        if name != 'cross_overlap_expression_coverage.py':
            expected = next(value for key, value in context.config['code_sha256'].items() if key.endswith('/' + name))
            assert context.inputs[str(path)] == expected, name
    if args.stage == 'extract':
        extract(context)
    else:
        assert args.extraction
        infer(context, args.extraction)
