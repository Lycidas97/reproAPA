"""Reconstruct all five null supports from effects, without the analysis helpers."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def adjusted(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    result = np.empty(len(values))
    result[order] = np.minimum(1, np.minimum.accumulate(
        (values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1])
    return result


def verify(root, config_path):
    config = json.loads(config_path.read_text())
    source = Path(config['analysis_root']) / 'cross_overlap'
    read = lambda path: pd.read_csv(path, sep='\t', float_precision='round_trip')
    frozen = read(source / 'ortholog_matched_features.tsv.gz')
    corrected = read(root / 'original-nulls/ortholog_matched_features.tsv.gz')
    extended = read(root / 'inference/ortholog-expression-coverage-binned.tsv.gz')
    pd.testing.assert_frame_equal(corrected[frozen.columns], frozen, check_dtype=False, rtol=0, atol=1e-12)
    pd.testing.assert_frame_equal(extended[frozen.columns], frozen, check_dtype=False, rtol=0, atol=1e-12)
    assert len(corrected) == 16092 and corrected.human_gene_id.is_unique and corrected.mouse_gene_id.is_unique
    nodes = read(source / 'all_mapped_nodes.tsv.gz')
    nodes = nodes[nodes.layer.eq('supported50') & nodes.measurement_pass]
    positions = nodes.pivot(index='gene_id', columns=['code', 'node'], values='position_index')
    paired = read(source / 'paired_pm_study_effects.tsv.gz')
    paired = paired[paired.layer.eq('supported50')]
    pairs = paired.pivot(index='gene_id', columns='code', values='mean_delta')
    for side, codes, gene_column in [('h', ['H1', 'H2', 'H3', 'H4'], 'human_gene_id'),
                                     ('m', ['W', 'C'], 'mouse_gene_id')]:
        ps = pd.DataFrame({code: positions.get((code, 'P'), pd.Series(dtype=float)) -
                                positions.get((code, 'S'), pd.Series(dtype=float)) for code in codes})
        for contrast, effects in [('PS', ps), ('PM', pairs.reindex(columns=codes))]:
            values = effects.reindex(corrected[gene_column]).reset_index(drop=True)
            counts = values.notna().sum(axis=1)
            labels = (values.le(-.02).all(axis=1) if side == 'm' else
                      values.le(-.02).sum(axis=1).ge(2) & values.ge(.02).sum(axis=1).eq(0))
            assert np.array_equal(counts, corrected[f'{side}_{contrast}_actual_effect_n'])
            assert np.array_equal(counts, extended[f'{side}_{contrast}_actual_effect_n'])
            assert np.array_equal(labels, frozen[f'{side}_{contrast}'])
    for stage in ['extraction', 'original-nulls', 'inference']:
        folder = root / stage
        result = json.loads((folder / 'result.json').read_text())
        assert result['status'] == 'complete' and all(result['checks'].values())
        for artifact in result['artifacts']:
            assert digest(folder / artifact['path']) == artifact['sha256']
        for path, fingerprint in json.loads((folder / 'audit.json').read_text())['inputs'].items():
            assert digest(path) == fingerprint, path
    for path, fingerprint in config['code_sha256'].items():
        if Path(path).parent == Path(__file__).parent:
            assert digest(path) == fingerprint
    original = read(root / 'original-nulls/overlap_null_results.tsv')
    extension = read(root / 'inference/expression-coverage-null-results.tsv')
    old_original = read(source / 'overlap_null_results.tsv')
    models = {
        'joint_measurement_PAS': None,
        'independent_q5_reference_capacity': ['effect_n', 'pas_q3', 'prox_capacity', 'ref_q5'],
        'independent_range10_reference_capacity': ['effect_n', 'pas_q3', 'prox_capacity', 'ref_min10', 'ref_max10'],
        'Q5_expression_coverage': ['effect_n', 'prox_capacity', 'ref_q5', 'pas_q3', 'mean_GEX_q3', 'APA_coverage_q3'],
        'Range_expression_coverage': ['effect_n', 'prox_capacity', 'pas_q3', 'ref_min10', 'ref_max10', 'mean_GEX_q3', 'APA_coverage_q3'],
    }
    draws = dict(np.load(root / 'original-nulls/overlap_null_draws.npz'))
    draws.update(dict(np.load(root / 'inference/null-draws.npz')))
    rows = []
    for contrast in ['PS', 'PM']:
        for model, fields in models.items():
            is_extension = model.endswith('expression_coverage')
            frame = extended if is_extension else corrected
            table = extension if is_extension else original
            row = table[table.contrast.eq(contrast) & table.model.eq(model)].iloc[0]
            probabilities, legacy = {}, {}
            for side in ['h', 'm']:
                label = f'{side}_{contrast}'
                old_columns = (['h_n', 'm_n', 'h_pas', 'm_pas'] if fields is None else
                               [label + '_' + field for field in fields])
                columns = old_columns + ([f'h_{contrast}_actual_effect_n', f'm_{contrast}_actual_effect_n']
                                         if fields is None else [label + '_actual_effect_n'])
                probabilities[side] = frame.groupby(columns, dropna=False)[label].transform('mean').to_numpy()
                legacy[side] = frame.groupby(old_columns, dropna=False)[label].transform('mean').to_numpy()
                ineligible = frame[label + '_actual_effect_n'].lt(2).to_numpy()
                assert not probabilities[side][ineligible].any(), (contrast, model, side)
                assert not frame.loc[ineligible, label].any()
            expected = float(np.sum(probabilities['h'] * probabilities['m']))
            reported = row.analytic_expected if is_extension else row.analytic_null_mean
            assert abs(expected - reported) < 1e-10
            old_expected = float(np.sum(legacy['h'] * legacy['m']))
            if not is_extension:
                old_row = old_original[old_original.contrast.eq(contrast) & old_original.model.eq(model)].iloc[0]
                assert abs(old_expected - old_row.analytic_null_mean) < 1e-10
            null = draws[contrast + '|' + model]
            observed = int((frame['h_' + contrast] & frame['m_' + contrast]).sum())
            assert observed == row.observed
            assert len(null) == (20000 if is_extension or fields is not None else 5000)
            assert abs(null.mean() - row.null_mean) < 1e-12
            p = (1 + np.count_nonzero(null >= observed)) / (len(null) + 1)
            assert abs(p - (row.p if is_extension else row.empirical_p)) < 1e-14
            assert abs(expected - null.mean()) <= max(.01, 6 * null.std(ddof=1) / np.sqrt(len(null)))
            record = dict(contrast=contrast, model=model, observed=observed,
                          human_events=int(frame['h_' + contrast].sum()), mouse_events=int(frame['m_' + contrast].sum()),
                          legacy_analytic_expected=old_expected, analytic_expected=expected,
                          null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                          null_q025=float(np.quantile(null, .025)), null_q975=float(np.quantile(null, .975)),
                          enrichment=observed / expected, p=p,
                          q_two_contrasts=row.q_two_contrasts if is_extension else row.q_within_model_two_contrasts,
                          q_four_tests=row.q_four_tests if is_extension else np.nan,
                          permutations=len(null))
            for side in ['h', 'm']:
                labels = frame[side + '_' + contrast].to_numpy(bool)
                ineligible = frame[f'{side}_{contrast}_actual_effect_n'].lt(2).to_numpy()
                prob = probabilities[side]
                record[side + '_legacy_expected_ineligible_events'] = float(legacy[side][ineligible].sum())
                record[side + '_expected_ineligible_events'] = float(prob[ineligible].sum())
                record[side + '_candidate_movable_fraction'] = float((labels & (prob > 0) & (prob < 1)).sum() / labels.sum())
            rows.append(record)
    for _, part in original.groupby('model'):
        assert np.allclose(part.q_within_model_two_contrasts, adjusted(part.empirical_p), rtol=0, atol=1e-14)
    for _, part in extension.groupby('model'):
        assert np.allclose(part.q_two_contrasts, adjusted(part.p), rtol=0, atol=1e-14)
    assert np.allclose(extension.q_four_tests, adjusted(extension.p), rtol=0, atol=1e-14)
    pd.DataFrame(rows).to_csv(root / 'all-five-null-comparison.tsv', sep='\t', index=False, float_format='%.17g')
    result = dict(status='PASS', null_tests=10, tests_skipped=0,
                  effect_counts_reconstructed_from_measured_effects=True,
                  frozen_features_and_candidate_members_unchanged=True,
                  all_five_nulls_have_zero_probability_outside_measurement_eligibility=True,
                  original_expectations_replayed=True, monte_carlo_and_tail_probabilities_verified=True,
                  separate_original_and_extended_BH_families_verified=True, artifact_and_input_hashes_verified=True)
    (root / 'all-five-verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    verify(args.root, args.config)
