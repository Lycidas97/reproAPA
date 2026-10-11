"""Independent checks of frozen cohorts, contributor weights and reported tests."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from polyapipe_prepare import digest


def adjusted(values):
    """BH from its definition, separately from the analysis helper."""
    values = np.asarray(values)
    order = np.argsort(values)
    answer = np.empty(len(values))
    for rank, index in enumerate(order, start=1):
        answer[index] = min(1., min(len(values) * values[j] / k
                                   for k, j in enumerate(order, start=1) if k >= rank))
    return answer


def verify(root, config):
    source = Path(json.loads(Path(config).read_text())['analysis_root'])
    read = lambda path: pd.read_csv(path, sep='\t', float_precision='round_trip')
    checks = []
    mobility = []
    for stage in ['extraction', 'inference']:
        folder = root / stage
        result = json.loads((folder / 'result.json').read_text())
        assert result['status'] == 'complete' and all(result['checks'].values())
        for artifact in result['artifacts']:
            assert digest(folder / artifact['path']) == artifact['sha256'], artifact['path']
        audit = json.loads((folder / 'audit.json').read_text())
        for path, fingerprint in audit['inputs'].items():
            assert digest(path) == fingerprint, path
        checks.append(dict(check=stage + '_output_and_input_hashes', status='PASS',
                           outputs=len(result['artifacts']), inputs=len(audit['inputs'])))
    frozen = read(source / 'cross_overlap/ortholog_matched_features.tsv.gz')
    frame = read(root / 'inference/ortholog-expression-coverage-binned.tsv.gz')
    pd.testing.assert_frame_equal(frame[frozen.columns], frozen, check_dtype=False)
    assert len(frame) == 16092
    checks.append(dict(check='all_frozen_columns_and_full_background_unchanged', status='PASS'))
    for code in ['W', 'C', 'H1', 'H2', 'H3', 'H4']:
        contributors = read(root / 'extraction' / (code + '-contributors.tsv.gz'))
        summaries = read(root / 'extraction' / (code + '-study-covariates.tsv.gz'))
        index = summaries.set_index(['contrast', 'gene_id'])
        for key, part in contributors.groupby(['contrast', 'gene_id']):
            assert abs(part.weight.sum() - 1) < 1e-12
            reference = 'S' if key[0] == 'PS' else 'M'
            assert abs(part[part.node.eq('P')].weight.sum() - .5) < 1e-12
            assert abs(part[part.node.eq(reference)].weight.sum() - .5) < 1e-12
            signed = np.where(part.node.eq('P'), 1., -1.)
            effect = float(2 * np.sum(part.weight * signed * part.position))
            assert abs(effect - index.loc[key, 'effect']) < 1e-12
            for column, contributor_column in [('mean_GEX', 'GEX_CPM'), ('APA_coverage', 'APA_coverage')]:
                values = part[contributor_column].to_numpy()
                mean = float(part.weight.to_numpy() @ values)
                saved = index.loc[key, column]
                assert (np.isnan(mean) and np.isnan(saved)) or abs(mean - saved) < 1e-10 * max(1, abs(mean))
        checks.append(dict(check=code + '_contributor_weight_and_covariate_reconstruction',
                           status='PASS', study_contrasts=len(summaries), contributor_states=len(contributors)))
    summary = read(root / 'inference/expression-coverage-null-results.tsv')
    diagnostics = read(root / 'inference/strata.tsv.gz')
    draws = np.load(root / 'inference/null-draws.npz')
    powers = read(root / 'inference/calibration-power.tsv')
    for row in summary.itertuples(index=False):
        base = (['effect_n', 'prox_capacity', 'ref_q5', 'pas_q3'] if row.model.startswith('Q5')
                else ['effect_n', 'prox_capacity', 'pas_q3', 'ref_min10', 'ref_max10'])
        fields = base + ['actual_effect_n', 'mean_GEX_q3', 'APA_coverage_q3']
        probabilities = {}
        for side in ['h', 'm']:
            labels = frame[side + '_' + row.contrast].astype(float)
            columns = [side + '_' + row.contrast + '_' + field for field in fields]
            probabilities[side] = frame.groupby(columns, dropna=False)[side + '_' + row.contrast].transform('mean').to_numpy()
            ineligible = frame[side + '_' + row.contrast + '_actual_effect_n'].lt(2)
            assert not labels[ineligible].any()
            assert not probabilities[side][ineligible].any()
            for covariate in ['mean_GEX_q3', 'APA_coverage_q3']:
                bins = frame[side + '_' + row.contrast + '_' + covariate]
                for code in bins.unique():
                    select = bins.eq(code).to_numpy()
                    assert abs(probabilities[side][select].sum() - labels[select].sum()) < 1e-10
        expected = float(np.sum(probabilities['h'] * probabilities['m']))
        assert abs(expected - row.analytic_expected) < 1e-10
        null = draws[row.contrast + '|' + row.model]
        assert len(null) == 20000
        assert abs(row.p - (1 + np.count_nonzero(null >= row.observed)) / 20001) < 1e-14
        assert abs(row.null_mean - null.mean()) < 1e-12
        assert abs(row.null_sd - null.std(ddof=1)) < 1e-12
        fixed = int(((probabilities['h'] == 1) & (probabilities['m'] == 1)).sum())
        assert fixed == row.fixed_overlap and null.min() >= fixed
        movement = dict(contrast=row.contrast, model=row.model, fixed_overlap=fixed,
                        observed=row.observed, null_unique=row.null_unique, null_sd=row.null_sd)
        for side in ['h', 'm']:
            labels = frame[side + '_' + row.contrast].to_numpy(bool)
            movable = (probabilities[side] > 0) & (probabilities[side] < 1)
            movement[side + '_candidate_movable_fraction'] = float((labels & movable).sum() / labels.sum())
            movement[side + '_expected_events_below_two_actual_studies'] = float(
                probabilities[side][frame[side + '_' + row.contrast + '_actual_effect_n'].lt(2)].sum())
        mobility.append(movement)
        for side in ['h', 'm']:
            table = diagnostics[diagnostics.contrast.eq(row.contrast) & diagnostics.model.eq(row.model) & diagnostics.side.eq(side)]
            mixed = table.events.gt(0) & table.events.lt(table.genes)
            assert abs(table.loc[mixed, 'genes'].sum() / 16092 - getattr(row, side + '_movable_fraction')) < 1e-12
        calibration = draws[row.contrast + '|' + row.model + '|OR1']
        mean_error = abs(calibration.mean() - null.mean())
        standard_error = np.sqrt(calibration.var(ddof=1) / len(calibration) + null.var(ddof=1) / len(null))
        assert mean_error <= max(.01, 6 * standard_error)
        for power in powers[powers.contrast.eq(row.contrast) & powers.model.eq(row.model)].itertuples(index=False):
            simulated = draws[row.contrast + '|' + row.model + '|OR' + str(power.odds_ratio)]
            p = np.array([(1 + np.count_nonzero(null >= value)) / 20001 for value in simulated])
            assert int((p <= power.alpha).sum()) == power.rejected
        checks.append(dict(check=row.contrast + '_' + row.model + '_independent_probability_and_balance', status='PASS'))
    assert np.allclose(summary.q_four_tests, adjusted(summary.p), rtol=0, atol=1e-14)
    for _, group in summary.groupby('model'):
        assert np.allclose(group.q_two_contrasts, adjusted(group.p), rtol=0, atol=1e-14)
    checks.append(dict(check='q2_and_q4_reconstructed_independently', status='PASS'))
    pd.DataFrame(mobility).to_csv(root / 'mobility-diagnostics.tsv', sep='\t', index=False, float_format='%.17g')
    result = dict(status='PASS', checks=checks, tests_skipped=0, manuscript_modified=False)
    (root / 'independent-verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--config', required=True)
    arguments = parser.parse_args()
    verify(arguments.root, arguments.config)
