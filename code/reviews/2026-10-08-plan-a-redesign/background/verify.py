"""Read-only independent reconstruction of frozen A1 outputs and moments."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
ROUTES = {'human': ['GSE118258'], 'mouse': ['CIP', 'WOT']}
TAU = .02


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(name, **kwargs):
    return pd.read_csv(OUT / name, sep='\t', float_precision='round_trip', **kwargs)


def independent_q2(x):
    return np.minimum(np.ceil(rankdata(x, method='average') / len(x) * 2), 2).astype(int) - 1


def matrix(f, routes, model):
    p, m, s = [f[[r + '__' + state + '_DWUI' for r in routes]].to_numpy() for state in ['P', 'M', 'S']]
    if model == 'S_vector':
        synthetic_s = np.broadcast_to(s[None, :, :], (len(f), len(f), len(routes)))
    else:
        z = np.zeros_like(m)
        for i in range(len(f)):
            for j in range(len(routes)):
                if s[i, j] > m[i, j]:
                    z[i, j] = (s[i, j] - m[i, j]) / (1 - m[i, j])
                elif s[i, j] < m[i, j]:
                    z[i, j] = (s[i, j] - m[i, j]) / m[i, j]
        synthetic_s = m[:, None, :] + np.maximum(z, 0)[None, :, :] * (1 - m[:, None, :]) + np.minimum(z, 0)[None, :, :] * m[:, None, :]
    assert np.all((synthetic_s >= 0) & (synthetic_s <= 1))
    ps = p[:, None, :] - synthetic_s
    pm_qualified = f.candidate.to_numpy() & (((p - m) * f.direction.to_numpy()[:, None]) >= TAU).all(axis=1)
    np.testing.assert_array_equal(pm_qualified, f.pm_eligible)
    return (pm_qualified[:, None] & ((ps * f.direction.to_numpy()[:, None, None]) >= TAU).all(axis=2)).astype(float)


def main():
    freeze = json.loads((HERE / 'freeze-receipt.json').read_text())
    audit = json.loads((OUT / 'input-audit.json').read_text())
    assert audit['status'] == 'COMPLETE' and audit['inference_tests'] == 12
    assert freeze['frozen_at_utc'] < audit['started']
    assert sha(freeze['global_protocol']) == freeze['global_sha256'] == audit['frozen_global_sha']
    assert sha(freeze['module_protocol']) == freeze['module_sha256'] == audit['protocol_sha256']
    assert sha(HERE / 'analyze.py') == freeze['code_sha256'] == audit['code_sha256']
    assert sha(HERE / 'test_model.py') == freeze['tests_sha256']
    assert freeze['input_sha256'] == audit['input_sha256']
    for path, fingerprint in freeze['input_sha256'].items():
        assert sha(path) == fingerprint
    inputs = read('background-input.tsv.gz')
    tests = read('raw-tests.tsv', dtype={'seed': 'string'})
    layers = read('layer-null-diagnostics.tsv.gz')
    marginals = read('candidate-null-probabilities.tsv.gz')
    verification = []
    for (source, species), f in inputs.groupby(['pas_source', 'species']):
        f = f.reset_index(drop=True)
        routes = ROUTES[species]
        source_paths = [Path(p) for p in freeze['input_sha256'] if source in p]
        if source == 'polyapipe':
            source_paths = [Path(p) for p in freeze['input_sha256'] if 'contract-r2' in p]
        base = next(p for p in source_paths if p.name == 'dual_reference_candidates.tsv').parents[1]
        selected = pd.read_csv(base / 'candidate_screen/dual_reference_candidates.tsv', sep='\t')
        recurrence = pd.read_csv(base / 'recurrence/recurrence_members.tsv.gz', sep='\t')
        r = recurrence[recurrence.species.eq(species) & recurrence.metric.eq('WPP') & recurrence.threshold.eq(TAU) & recurrence.series_matching & recurrence.contrast.eq('P_M')]
        signs = r.set_index('gene_id').direction.map({'negative': -1, 'positive': 1})
        np.testing.assert_array_equal(f.direction, f.gene_id.map(signs).fillna(0).astype(int))
        assert set(f.loc[f.dual, 'gene_id']) == set(selected.loc[selected.species.eq(species), 'gene_id'])
        frozen_matched = pd.read_csv(base / 'matched_pms/matched_state_vectors.tsv.gz', sep='\t', float_precision='round_trip')
        frozen_matched = frozen_matched[frozen_matched.setting.eq('supported50_p08') & frozen_matched.species.eq(species)]
        available = frozen_matched[frozen_matched.route.isin(routes)].groupby('gene_id').route.nunique()
        assert set(f.gene_id) == set(available[available.eq(len(routes))].index)
        for route in routes:
            donor = frozen_matched[frozen_matched.route.eq(route)].set_index('gene_id').reindex(f.gene_id)
            for name in ['P_DWUI', 'M_DWUI', 'S_DWUI', 'DWUI_PM', 'DWUI_PS', 'DWUI_MS', 'PAS_number']:
                np.testing.assert_array_equal(f[route + '__' + name], donor[name])
        m = f[[r + '__M_DWUI' for r in routes]].to_numpy()
        pm = f[[r + '__DWUI_PM' for r in routes]].to_numpy()
        np.testing.assert_array_equal(f.M_median_q2, independent_q2(np.median(m, axis=1)))
        np.testing.assert_array_equal(f.abs_PM_median_q2, independent_q2(np.median(abs(pm), axis=1)))
        ranks = [rankdata(f[[r + '__' + x for r in routes]].mean(axis=1), method='average') / len(f) for x in ['mean_GEX', 'APA_coverage']]
        np.testing.assert_allclose(f.measurement_composite, np.mean(ranks, axis=0), atol=1e-15)
        np.testing.assert_array_equal(f.measurement_q2, independent_q2(np.mean(ranks, axis=0)))
        np.testing.assert_array_equal(f.PAS_median_q2, independent_q2(f[[r + '__PAS_number' for r in routes]].median(axis=1)))
        for model in ['bounded_MS', 'S_vector']:
            mean, variance = 0., 0.
            for key, q in f.groupby('stratum'):
                q = q.reset_index(drop=True)
                a = matrix(q, routes, model)
                np.testing.assert_array_equal(a.diagonal(), q.dual)
                n, total, squares = len(q), float(a.sum()), float((a * a).sum())
                mu = total / n
                cross = total ** 2 - np.square(a.sum(axis=0)).sum() - np.square(a.sum(axis=1)).sum() + squares
                second = squares / n + cross / (n * (n - 1)) if n > 1 else mu ** 2
                var = float(second - mu ** 2)
                mean += mu; variance += var
                record = layers[layers.pas_source.eq(source) & layers.species.eq(species) & layers.model.eq(model) & layers.stratum.eq(key)].iloc[0]
                assert abs(mu - record['mean']) < 1e-12
                assert abs(var - record.variance) < 1e-10
                margins = marginals[marginals.pas_source.eq(source) & marginals.species.eq(species) & marginals.model.eq(model)].set_index('gene_id')
                for index, gene in enumerate(q.gene_id):
                    if q.candidate.iloc[index]:
                        assert abs(a[index].mean() - margins.loc[gene, 'null_alignment_probability']) < 1e-12
            row = tests[tests.pas_source.eq(source) & tests.species.eq(species) & tests.test.eq(model)].iloc[0]
            assert abs(mean - row.analytic_null_mean) < 1e-10
            assert abs(variance - row.analytic_null_variance) < 1e-9
            archive = np.load(OUT / f'null-{source}-{species}-{model}.npz')
            draws = archive['draws']; exceed = int((draws >= row.observed).sum())
            assert len(draws) == 20000 and exceed == row.exceedances
            assert abs((1 + exceed) / 20001 - row.p_value) < 1e-15
            assert abs(draws.mean() - row.simulated_mean) < 1e-12
            assert str(archive['seed']) == row.seed
            expected_seed = int.from_bytes(hashlib.sha256(f'A1-redesign|{source}|{species}|{model}|20261008'.encode()).digest()[:8], 'big')
            assert str(expected_seed) == row.seed
            assert draws.min() >= row.attainable_min and draws.max() <= row.attainable_max
            assert abs(draws.mean() - mean) <= max(.02, 5 * np.sqrt(variance / len(draws)))
            verification.append(dict(pas_source=source, species=species, model=model, status='PASS',
                                     independent_mean=mean, independent_variance=variance, exceedances=exceed))
        row = tests[tests.pas_source.eq(source) & tests.species.eq(species) & tests.test.eq('direction_composition')].iloc[0]
        selected = f[f.candidate]
        a, b = int((selected.dual & selected.direction.eq(-1)).sum()), int((selected.dual & selected.direction.eq(1)).sum())
        c, d = int((~selected.dual & selected.direction.eq(-1)).sum()), int((~selected.dual & selected.direction.eq(1)).sum())
        assert [a, b, c, d] == [row.dual_proximal, row.dual_distal, row.complement_proximal, row.complement_distal]
        n, successes, draws = a + b + c + d, a + c, a + b
        mass = hypergeom.pmf(np.arange(max(0, draws - (n - successes)), min(draws, successes) + 1), n, successes, draws)
        observed_mass = hypergeom.pmf(a, n, successes, draws)
        assert abs(mass[mass <= observed_mass * (1 + 1e-12)].sum() - row.p_value) < 1e-12
    pd.DataFrame(verification).to_csv(OUT / 'independent-verification.tsv', sep='\t', index=False, float_format='%.17g')
    output_hashes = {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'verification.json'}
    report = dict(status='PASS', finished_at_utc=datetime.now(timezone.utc).isoformat(), tests_skipped=0,
                  fixed_original_candidate_and_dual_IDs=True, continuous_rank_bins_independently_replayed=True,
                  shared_P_and_PM_fixed=True, two_route_joint_donor=True, independent_permutation_second_moments=True,
                  saved_MC_p_and_seeds_replayed=True, independent_two_sided_Fisher=True, all_original_inputs_unchanged=True,
                  frozen_receipt_sha256=sha(HERE / 'freeze-receipt.json'), verification_code_sha256=sha(Path(__file__)),
                  outputs_sha256=output_hashes)
    (OUT / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
