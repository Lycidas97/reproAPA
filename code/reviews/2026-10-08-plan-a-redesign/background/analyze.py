"""A1 conditional joint reference randomization; inference requires a frozen protocol."""
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import linear_sum_assignment
from scipy.stats import beta, fisher_exact

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASES = {
    'polyapipe': ROOT / 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables',
    'apadetect-p50': ROOT / 'notes/2026-10-06-dual-main-plot-data/apadetect-p50/supplement/tables',
}
ROUTES = {'human': ['GSE118258'], 'mouse': ['CIP', 'WOT']}
TAU = .02
B = 20000
INPUTS = {}


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read(base, relative):
    path = base / relative
    INPUTS[str(path)] = digest(path)
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def save(frame, out, name):
    kwargs = {'compression': {'method': 'gzip', 'mtime': 0}} if name.endswith('.gz') else {}
    pd.DataFrame(frame).to_csv(out / name, sep='\t', index=False, na_rep='NA', float_format='%.17g', **kwargs)


def rank_bin(x):
    assert np.isfinite(x).all()
    return np.minimum(np.ceil(x.rank(method='average', pct=True) * 2).astype(int), 2) - 1


def bounded_residual(m, s):
    assert np.all((m >= 0) & (m <= 1) & (s >= 0) & (s <= 1))
    z = np.zeros_like(m, dtype=float)
    upper, lower = s > m, s < m
    z[upper] = (s[upper] - m[upper]) / (1 - m[upper])
    z[lower] = (s[lower] - m[lower]) / m[lower]
    assert np.all(np.isfinite(z) & (z >= -1) & (z <= 1))
    return z


def inverse_residual(m, z):
    s = m + np.where(z >= 0, (1 - m) * z, m * z)
    assert np.all(np.isfinite(s) & (s >= 0) & (s <= 1))
    return s


def build(source):
    base = BASES[source]
    matched = read(base, 'matched_pms/matched_state_vectors.tsv.gz')
    matched = matched[matched.setting.eq('supported50_p08')].copy()
    recurrent = read(base, 'recurrence/recurrence_members.tsv.gz')
    recurrent = recurrent[recurrent.metric.eq('WPP') & recurrent.threshold.eq(TAU) &
                          recurrent.series_matching & recurrent.contrast.eq('P_M')].copy()
    recurrent['sign'] = recurrent.direction.map({'negative': -1, 'positive': 1})
    dual = read(base, 'candidate_screen/dual_reference_candidates.tsv')
    ladder = read(base, 'candidate_screen/dual_reference_ladder.tsv')
    effects = read(base, 'state_vectors/route_effects.primary.tsv.gz')
    covariates = effects[effects.contrast.eq('P_M')][['species', 'route', 'gene_id', 'mean_GEX', 'APA_coverage']]
    matched = matched.merge(covariates, on=['species', 'route', 'gene_id'], validate='one_to_one')
    output = []
    for species, routes in ROUTES.items():
        q = matched[matched.species.eq(species) & matched.route.isin(routes)].copy()
        assert not q.duplicated(['gene_id', 'route']).any()
        available = q.groupby('gene_id').route.nunique()
        ids = available[available.eq(len(routes))].index
        q = q[q.gene_id.isin(ids)]
        columns = ['P_DWUI', 'M_DWUI', 'S_DWUI', 'DWUI_PM', 'DWUI_PS', 'DWUI_MS',
                   'mean_GEX', 'APA_coverage', 'PAS_number']
        f = pd.concat([q[q.route.eq(route)].set_index('gene_id')[columns].add_prefix(route + '__')
                       for route in routes], axis=1).sort_index().reset_index()
        assert np.isfinite(f.drop(columns='gene_id').to_numpy()).all()
        f['pas_source'], f['species'] = source, species
        f['gene_symbol'] = f.gene_id.map(q.drop_duplicates('gene_id').set_index('gene_id').gene_symbol)
        r = recurrent[recurrent.species.eq(species)].set_index('gene_id')
        assert r.index.is_unique
        f['candidate'] = f.gene_id.isin(r.index)
        f['direction'] = f.gene_id.map(r.sign).fillna(0).astype(int)
        p, m, s, pm, ps, ms = [f[[route + '__' + name for route in routes]].to_numpy(float)
                              for name in ['P_DWUI', 'M_DWUI', 'S_DWUI', 'DWUI_PM', 'DWUI_PS', 'DWUI_MS']]
        np.testing.assert_allclose(p - m, pm, atol=1e-12)
        np.testing.assert_allclose(p - s, ps, atol=1e-12)
        np.testing.assert_allclose(pm + ms, ps, atol=1e-12)
        z = bounded_residual(m, s)
        np.testing.assert_allclose(inverse_residual(m, z), s, atol=1e-12)
        sign = f.direction.to_numpy()[:, None]
        f['pm_eligible'] = f.candidate & (pm * sign >= TAU).all(axis=1)
        f['dual'] = f.pm_eligible & (ps * sign >= TAU).all(axis=1)
        assert set(f.loc[f.dual, 'gene_id']) == set(dual.loc[dual.species.eq(species), 'gene_id'])
        counts = ladder[ladder.species.eq(species)].iloc[0]
        assert int(f.candidate.sum()) == counts.matched_measurable
        assert int(f.dual.sum()) == counts.dual_reference
        f['M_median_q2'] = rank_bin(pd.Series(np.median(m, axis=1)))
        f['abs_PM_median_q2'] = rank_bin(pd.Series(np.median(abs(pm), axis=1)))
        measurement_ranks = [f[[route + '__' + name for route in routes]].mean(axis=1).rank(method='average', pct=True)
                             for name in ['mean_GEX', 'APA_coverage']]
        f['measurement_composite'] = pd.concat(measurement_ranks, axis=1).mean(axis=1)
        f['measurement_q2'] = rank_bin(f.measurement_composite)
        f['PAS_median_q2'] = rank_bin(f[[route + '__PAS_number' for route in routes]].median(axis=1))
        f['PM_sign_pattern'] = pd.DataFrame(np.sign(pm).astype(int)).astype(str).agg('|'.join, axis=1)
        boundary = np.where(m == 0, 'zero', np.where(m == 1, 'one', 'interior'))
        f['M_boundary_pattern'] = pd.DataFrame(boundary).agg('|'.join, axis=1)
        fields = ['M_median_q2', 'abs_PM_median_q2', 'measurement_q2', 'PAS_median_q2',
                  'PM_sign_pattern', 'M_boundary_pattern']
        f['stratum'] = f[fields].astype(str).agg('||'.join, axis=1)
        for j, route in enumerate(routes):
            f[route + '__bounded_S_minus_M'] = z[:, j]
        output.append(f)
    return pd.concat(output, ignore_index=True)


def structure(frame):
    summaries, layers = [], []
    for (source, species), f in frame.groupby(['pas_source', 'species'], sort=True):
        routes = ROUTES[species]
        for key, q in f.groupby('stratum', sort=True):
            z = q[[route + '__bounded_S_minus_M' for route in routes]].to_numpy()
            s = q[[route + '__S_DWUI' for route in routes]].to_numpy()
            layers.append(dict(pas_source=source, species=species, stratum=key, genes=len(q),
                               candidates=int(q.candidate.sum()), eligible_candidates=int(q.pm_eligible.sum()),
                               singleton=len(q) == 1, distinct_residual_vectors=len(np.unique(z, axis=0)),
                               distinct_S_vectors=len(np.unique(s, axis=0))))
        a = pd.DataFrame([x for x in layers if x['pas_source'] == source and x['species'] == species])
        movable = ~a.singleton
        summaries.append(dict(pas_source=source, species=species, background=len(f), candidates=int(f.candidate.sum()),
                              pm_eligible_candidates=int(f.pm_eligible.sum()), observed_dual=int(f.dual.sum()),
                              strata=len(a), singleton_genes=int(a.loc[a.singleton, 'genes'].sum()),
                              movable_candidates=int(a.loc[movable, 'candidates'].sum()),
                              movable_eligible_candidates=int(a.loc[movable, 'eligible_candidates'].sum()),
                              residual_variable_candidate_genes=int(a.loc[a.distinct_residual_vectors.gt(1), 'candidates'].sum()),
                              S_variable_candidate_genes=int(a.loc[a.distinct_S_vectors.gt(1), 'candidates'].sum())))
    return summaries, layers


def hit_matrix(f, species, model):
    routes = ROUTES[species]
    p = f[[route + '__P_DWUI' for route in routes]].to_numpy(float)
    m = f[[route + '__M_DWUI' for route in routes]].to_numpy(float)
    if model == 'bounded_MS':
        z = f[[route + '__bounded_S_minus_M' for route in routes]].to_numpy(float)
        s = inverse_residual(m[:, None, :], z[None, :, :])
    else:
        s = np.broadcast_to(f[[route + '__S_DWUI' for route in routes]].to_numpy(float)[None, :, :], (len(f), len(f), len(routes)))
    ps = p[:, None, :] - s
    pm = p - m
    np.testing.assert_allclose(ps, pm[:, None, :] + m[:, None, :] - s, atol=1e-12)
    hits = f.pm_eligible.to_numpy()[:, None] & (ps * f.direction.to_numpy()[:, None, None] >= TAU).all(axis=2)
    assert np.array_equal(np.diag(hits), f.dual.to_numpy())
    return hits.astype(np.int8)


def assignment_moments(hits):
    n = len(hits)
    mean = float(hits.sum() / n)
    centered = hits - hits.mean(axis=0)[None, :] - hits.mean(axis=1)[:, None] + hits.mean()
    variance = float((centered ** 2).sum() / (n - 1)) if n > 1 else 0.
    return mean, variance


def randomization(frame, source, species, model, out):
    seed = int.from_bytes(hashlib.sha256(f'A1-redesign|{source}|{species}|{model}|20261008'.encode()).digest()[:8], 'big')
    rng = np.random.default_rng(seed)
    draws = np.zeros(B, int)
    mean = variance = 0.
    support_min = support_max = 0
    rows, marginal = [], []
    for key, f in frame.groupby('stratum', sort=True):
        f = f.reset_index(drop=True)
        hits = hit_matrix(f, species, model)
        layer_mean, layer_variance = assignment_moments(hits)
        mean += layer_mean; variance += layer_variance
        i, j = linear_sum_assignment(hits)
        low = int(hits[i, j].sum())
        i, j = linear_sum_assignment(hits, maximize=True)
        high = int(hits[i, j].sum())
        support_min += low; support_max += high
        for start in range(0, B, 500):
            count = min(500, B - start)
            order = np.argsort(rng.random((count, len(f))), axis=1)
            draws[start:start + count] += hits[np.arange(len(f))[None, :], order].sum(axis=1)
        rows.append(dict(pas_source=source, species=species, model=model, stratum=key, genes=len(f),
                         candidates=int(f.candidate.sum()), observed=int(f.dual.sum()), mean=layer_mean,
                         variance=layer_variance, attainable_min=low, attainable_max=high))
        for index, gene in enumerate(f.gene_id):
            if f.candidate.iloc[index]:
                marginal.append(dict(pas_source=source, species=species, model=model, gene_id=gene,
                                     direction=int(f.direction.iloc[index]), pm_eligible=bool(f.pm_eligible.iloc[index]),
                                     observed=bool(f.dual.iloc[index]), null_alignment_probability=float(hits[index].mean())))
    assert np.max(draws) <= support_max and np.min(draws) >= support_min
    assert abs(draws.mean() - mean) <= max(.02, 5 * np.sqrt(variance / B))
    observed = int(frame.dual.sum())
    exceed = int((draws >= observed).sum())
    np.savez_compressed(out / f'null-{source}-{species}-{model}.npz', draws=draws, observed=observed, seed=np.array(str(seed)))
    result = dict(test_id=f'{source}|{species}|{model}', family='A1-conditional-reference-redesign',
                  pas_source=source, species=species, test=model, observed=observed, candidates=int(frame.candidate.sum()),
                  pm_eligible_candidates=int(frame.pm_eligible.sum()), background=len(frame),
                  p_value=(1 + exceed) / (B + 1), exceedances=exceed, permutations=B, seed=str(seed),
                  tail_probability_ci_low=float(beta.ppf(.025, exceed, B - exceed + 1)) if exceed else 0.,
                  tail_probability_ci_high=float(beta.ppf(.975, exceed + 1, B - exceed)) if exceed < B else 1.,
                  analytic_null_mean=mean, analytic_null_variance=variance, simulated_mean=float(draws.mean()),
                  null_p025=float(np.quantile(draws, .025)), null_p975=float(np.quantile(draws, .975)),
                  attainable_min=support_min, attainable_max=support_max,
                  interpretation='post-hoc conditional gene association; fixed P/M and candidate selection')
    return result, rows, marginal


def direction_test(frame, source, species):
    candidates = frame[frame.candidate]
    a = int((candidates.dual & candidates.direction.eq(-1)).sum())
    b = int((candidates.dual & candidates.direction.eq(1)).sum())
    c = int((~candidates.dual & candidates.direction.eq(-1)).sum())
    d = int((~candidates.dual & candidates.direction.eq(1)).sum())
    odds, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
    return dict(test_id=f'{source}|{species}|direction_composition', family='A1-direction-composition-redesign',
                pas_source=source, species=species, test='direction_composition', p_value=float(p), odds_ratio=float(odds),
                dual_proximal=a, dual_distal=b, complement_proximal=c, complement_distal=d,
                method='two-sided Fisher; fixed measurable PM complement')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=['structure', 'inference'], required=True)
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--frozen-protocol', type=Path)
    parser.add_argument('--frozen-sha')
    args = parser.parse_args()
    if args.phase == 'inference':
        assert args.frozen_protocol and args.frozen_sha, 'Root must provide the frozen global protocol and SHA before inference.'
        assert digest(args.frozen_protocol) == args.frozen_sha, 'Frozen global protocol changed.'
    args.out.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    frame = pd.concat([build(source) for source in BASES], ignore_index=True)
    summaries, layers = structure(frame)
    save(frame, args.out, 'background-input.tsv.gz')
    save(summaries, args.out, 'structure-summary.tsv')
    save(layers, args.out, 'strata-structure.tsv.gz')
    results, diagnostics, marginals = [], [], []
    if args.phase == 'inference':
        for (source, species), f in frame.groupby(['pas_source', 'species'], sort=True):
            for model in ['bounded_MS', 'S_vector']:
                result, rows, marginal = randomization(f, source, species, model, args.out)
                results.append(result); diagnostics.extend(rows); marginals.extend(marginal)
            results.append(direction_test(f, source, species))
        save(results, args.out, 'raw-tests.tsv')
        save(diagnostics, args.out, 'layer-null-diagnostics.tsv.gz')
        save(marginals, args.out, 'candidate-null-probabilities.tsv.gz')
    assert all(digest(Path(path)) == sha for path, sha in INPUTS.items())
    audit = dict(status='STRUCTURE_ONLY_NO_P' if args.phase == 'structure' else 'COMPLETE',
                 started=started, finished=datetime.now(timezone.utc).isoformat(), phase=args.phase,
                 code_sha256=digest(Path(__file__)), input_sha256=INPUTS,
                 protocol_sha256=digest(HERE / 'PROTOCOL.draft.zh.md'),
                 frozen_global_protocol=str(args.frozen_protocol) if args.frozen_protocol else None,
                 frozen_global_sha=args.frozen_sha, python=sys.version, numpy=np.__version__,
                 pandas=pd.__version__, scipy=scipy.__version__, original_inputs_unchanged=True,
                 inference_tests=len(results))
    (args.out / 'input-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(pd.DataFrame(summaries).to_string(index=False))
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == '__main__':
    main()
