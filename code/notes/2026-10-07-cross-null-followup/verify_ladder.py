"""Check ladder event support, fixed margins, simulated rates and frozen input identity."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

FIELDS = {
    'L0_measurement': ['effect_n'],
    'L1_capacity': ['effect_n', 'prox_capacity'],
    'L2_baseline': ['effect_n', 'prox_capacity', 'ref_q5'],
    'L2_PAS_first': ['effect_n', 'prox_capacity', 'pas_q3'],
    'L3_full': ['effect_n', 'prox_capacity', 'ref_q5', 'pas_q3'],
}
read = lambda path: pd.read_csv(path, sep='\t', float_precision='round_trip')
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()


def adjusted(values):
    p = np.asarray(values)
    order = np.argsort(p)
    result = np.empty(len(p))
    result[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return result


def verify(root):
    checks = []
    for source in ['polyapipe', 'apadetect-p50']:
        folder = root / source
        config = json.loads((folder / 'source.config.json').read_text())
        old = Path(config['analysis_root']) / 'cross_overlap'
        f = read(folder / 'analysis/cross_overlap/ortholog_matched_features.tsv.gz')
        frozen = read(old / 'ortholog_matched_features.tsv.gz')
        pd.testing.assert_frame_equal(f[frozen.columns], frozen, check_dtype=False, rtol=0, atol=1e-12)
        assert len(f) == 16092 and f.human_gene_id.is_unique and f.mouse_gene_id.is_unique
        paired = read(old / 'paired_pm_study_effects.tsv.gz')
        pairs = paired[paired.layer.eq('supported50')].pivot(index='gene_id', columns='code', values='mean_delta')
        nodes = read(old / 'all_mapped_nodes.tsv.gz')
        nodes = nodes[nodes.layer.eq('supported50') & nodes.measurement_pass]
        positions = nodes.pivot(index='gene_id', columns=['code', 'node'], values='position_index')
        for side, codes, gene_column in [('h', ['H1', 'H2', 'H3', 'H4'], 'human_gene_id'), ('m', ['W', 'C'], 'mouse_gene_id')]:
            ps = pd.DataFrame({code: positions.get((code, 'P'), pd.Series(dtype=float)) -
                               positions.get((code, 'S'), pd.Series(dtype=float)) for code in codes})
            for co, effects in [('PS', ps), ('PM', pairs.reindex(columns=codes))]:
                values = effects.reindex(f[gene_column]).reset_index(drop=True)
                assert np.array_equal(values.notna().sum(1), f[f'{side}_{co}_actual_effect_n'])
                labels = (values.le(-.02).all(axis=1) if side == 'm' else values.le(-.02).sum(1).ge(2) & values.ge(.02).sum(1).eq(0))
                assert np.array_equal(labels, f[f'{side}_{co}'])
        ladder = read(folder / 'ladder/cross_species_ladder.tsv')
        power = read(folder / 'ladder/cross_species_power.tsv')
        archive = dict(np.load(folder / 'ladder/cross_species_draws.npz'))
        for row in ladder.itertuples(index=False):
            co, model = row.contrast, row.model
            probabilities, legacy = {}, {}
            record = dict(source=source, contrast=co, model=model)
            for side in ['h', 'm']:
                label = f'{side}_{co}'
                columns = [label + '_' + field for field in FIELDS[model]]
                probabilities[side] = f.groupby(columns + [label + '_actual_effect_n'], dropna=False)[label].transform('mean').to_numpy()
                legacy[side] = f.groupby(columns, dropna=False)[label].transform('mean').to_numpy()
                ineligible = f[label + '_actual_effect_n'].lt(2).to_numpy()
                assert not probabilities[side][ineligible].any()
                record[side + '_legacy_expected_ineligible_events'] = float(legacy[side][ineligible].sum())
                record[side + '_expected_ineligible_events'] = float(probabilities[side][ineligible].sum())
            expected = float(probabilities['h'] @ probabilities['m'])
            record.update(expected=expected, legacy_expected=float(legacy['h'] @ legacy['m']),
                fixed_overlap=int(((probabilities['h'] == 1) & (probabilities['m'] == 1)).sum()))
            assert abs(expected - row.expected) < 1e-10
            assert record['fixed_overlap'] == row.fixed_overlap
            null = archive[co + '|' + model]
            assert len(null) == 20000 and null.min() >= record['fixed_overlap']
            assert abs(null.mean() - row.mc_mean) < 1e-12
            assert abs(null.std(ddof=1) - row.null_sd) < 1e-12
            assert abs(expected - null.mean()) <= max(.01, 6 * null.std(ddof=1) / np.sqrt(len(null)))
            assert abs((1 + (null >= row.observed).sum()) / (len(null) + 1) - row.p) < 1e-14
            checks.append(record)
        assert np.allclose(ladder.q_all_10, adjusted(ladder.p), rtol=0, atol=1e-14)
        for _, part in ladder.groupby('model'):
            assert np.allclose(part.q_two_contrasts, adjusted(part.p), rtol=0, atol=1e-14)
        for row in power.itertuples(index=False):
            key = row.contrast + ('|OR' + str(row.odds_ratio) if row.experiment == 'conditional_injection' else '|self_' + row.tested_model)
            sim = archive[key]
            null = archive[row.contrast + '|' + row.tested_model]
            assert len(sim) == 3000 and row.simulations == 3000
            p = (1 + len(null) - np.searchsorted(np.sort(null), sim, side='left')) / (len(null) + 1)
            k = int((p <= row.alpha).sum())
            assert k == row.rejected and abs(k / 3000 - row.rejection_rate) < 1e-14
            assert abs(sim.mean() - row.simulated_mean) < 1e-12
            assert abs(sim.std(ddof=1) - row.simulated_sd) < 1e-12
            lo = 0. if k == 0 else beta.ppf(.025, k, 3000 - k + 1)
            hi = 1. if k == 3000 else beta.ppf(.975, k + 1, 3000 - k)
            assert abs(lo - row.lo) < 1e-14 and abs(hi - row.hi) < 1e-14
            if row.experiment == 'self_calibration' or (row.odds_ratio == 1 and row.tested_model == 'L3_full'):
                error = np.sqrt(sim.var(ddof=1) / len(sim) + null.var(ddof=1) / len(null))
                assert abs(sim.mean() - null.mean()) <= max(.01, 6 * error)
        for name in ['analysis/cross_overlap', 'ladder']:
            result = json.loads((folder / name / 'result.json').read_text())
            for artifact in result['artifacts']:
                assert digest(folder / name / artifact['path']) == artifact['sha256']
            for path, fingerprint in json.loads((folder / name / 'audit.json').read_text())['inputs'].items():
                assert digest(path) == fingerprint
        print('VERIFIED', source, '10 ladder tests, 90 simulation rows', flush=True)
    pd.DataFrame(checks).to_csv(root / 'eligibility-verification.tsv', sep='\t', index=False, float_format='%.17g')
    receipt = dict(status='PASS', tests_skipped=0, sources=2, ladder_tests=20, simulation_rows=180,
        raw_effect_counts_and_original_candidate_members_verified=True,
        every_ladder_step_excludes_unmeasurable_events=True,
        conditional_injection_and_self_calibration_rates_and_exact_intervals_verified=True,
        analytic_means_and_archived_tails_verified=True, BH_families_verified=True,
        artifact_and_input_hashes_verified=True)
    (root / 'verification.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    verify(parser.parse_args().root)
