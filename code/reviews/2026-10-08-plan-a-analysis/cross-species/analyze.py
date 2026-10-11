"""A3: fixed dual-reference sets, measurable ortholog background, conditional power."""
import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import beta, hypergeom, nchypergeom_fisher

SEED = 20261008
B = 20000
POWER_N = 3000
ODDS = [1, 1.5, 2, 3, 4, 6, 8, 12, 16, 32]
ALPHAS = [.05, .0125]
FAMILY = 'A3-cross-species'
ROUTES = {'human': ['GSE118258'], 'mouse': ['CIP', 'WOT']}
ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def random(key):
    seed = int.from_bytes(hashlib.sha256(f'{SEED}|{key}'.encode()).digest()[:8], 'little')
    return np.random.default_rng(seed), seed


def qbin(x, n):
    assert np.isfinite(x).all()
    return pd.qcut(x, n, labels=False, duplicates='drop').fillna(0).astype(int)


def groups(frame, fields):
    return [np.asarray(v, int) for v in frame.groupby(fields, sort=True, dropna=False).indices.values()]


def draw_labels(labels, strata, n, rng):
    result = np.zeros((n, len(labels)), dtype=np.int8)
    for ii in strata:
        values = labels[ii]
        if len(np.unique(values)) == 1:
            result[:, ii] = values[0]
        else:
            order = np.argsort(rng.random((n, len(ii))), axis=1)
            result[:, ii] = values[order]
    return result


def probability(reference, observed):
    ordered = np.sort(reference)
    return (1 + len(reference) - np.searchsorted(ordered, observed, side='left')) / (len(reference) + 1)


def rate(events):
    k, n = int(np.sum(events)), len(events)
    return dict(rejections=k, simulations=n, power=k/n,
                power_ci_low=float(beta.ppf(.025, k, n-k+1)) if k else 0.,
                power_ci_high=float(beta.ppf(.975, k+1, n-k)) if k < n else 1.)


def build(source, input_paths):
    base = ROOT / 'notes/2026-10-06-dual-main-plot-data' / source / 'supplement/tables'
    def read(rel):
        path = base / rel
        input_paths[str(path)] = digest(path)
        return pd.read_csv(path, sep='\t', float_precision='round_trip')
    matched = read('matched_pms/matched_state_vectors.tsv.gz')
    matched = matched[matched.setting.eq('supported50_p08')].copy()
    orth = read('cross_overlap/ortholog_matched_features.tsv.gz')
    candidates = read('candidate_screen/dual_reference_candidates.tsv')
    ladder = read('candidate_screen/dual_reference_ladder.tsv')
    recurrent = read('recurrence/recurrence_members.tsv.gz')
    recurrent = recurrent[recurrent.metric.eq('WPP') & recurrent.threshold.eq(.02) &
                          recurrent.series_matching & recurrent.contrast.eq('P_M')]
    assert orth.human_gene_id.is_unique and orth.mouse_gene_id.is_unique
    assert orth.homology_type.eq('ortholog_one2one').all()
    assert orth.human_partner_count.eq(1).all() and orth.mouse_partner_count.eq(1).all()
    records, eligible = {}, {}
    for species, routes in ROUTES.items():
        q = matched[matched.species.eq(species) & matched.route.isin(routes)]
        assert not q.duplicated(['gene_id', 'route']).any()
        assert np.isfinite(q[['P_DWUI', 'M_DWUI', 'S_DWUI']]).all().all()
        counts = q.groupby('gene_id').route.nunique()
        eligible[species] = set(counts[counts.eq(len(routes))].index)
        r = recurrent[recurrent.species.eq(species)]
        selected = candidates[candidates.species.eq(species)]
        assert selected.gene_id.is_unique and r.gene_id.is_unique
        assert set(selected.gene_id) <= set(r.gene_id) & eligible[species]
        expected = []
        for gene in sorted(set(r.gene_id) & eligible[species]):
            z = q[q.gene_id.eq(gene)]
            d = r.set_index('gene_id').loc[gene, 'direction']
            sign = -1 if d == 'negative' else 1
            if (z[['DWUI_PS', 'DWUI_PM']] * sign).ge(.02).all().all():
                expected.append(gene)
        assert set(expected) == set(selected.gene_id), (source, species)
        row = ladder[ladder.species.eq(species)].iloc[0]
        assert len(r) == row.recurrent_candidates
        assert len(set(r.gene_id) & eligible[species]) == row.matched_measurable
        assert len(selected) == row.dual_reference
        wide = q.pivot(index='gene_id', columns='route', values=['S_DWUI', 'M_DWUI', 'PAS_number'])
        records[species] = wide
    universe = orth[orth.human_gene_id.isin(eligible['human']) & orth.mouse_gene_id.isin(eligible['mouse'])].copy().reset_index(drop=True)
    for species, side in [('human', 'h'), ('mouse', 'm')]:
        ids = universe[species + '_gene_id']
        q = records[species].reindex(ids).reset_index(drop=True)
        s = q['S_DWUI'][ROUTES[species]]
        m = q['M_DWUI'][ROUTES[species]]
        p = q['PAS_number'][ROUTES[species]]
        for route in ROUTES[species]:
            universe[f'{side}_{route}_S_DWUI'] = s[route]
            universe[f'{side}_{route}_M_DWUI'] = m[route]
            universe[f'{side}_{route}_PAS_number'] = p[route]
        for name, ref in [('S', s), ('M', m)]:
            universe[f'{side}_{name}_prox_capacity'] = ref.ge(.02).sum(axis=1)
            universe[f'{side}_{name}_dist_capacity'] = ref.le(.98).sum(axis=1)
        universe[f'{side}_ref_median'] = pd.concat([s, m], axis=1).median(axis=1)
        universe[f'{side}_ref_q5'] = qbin(universe[f'{side}_ref_median'], 5)
        universe[f'{side}_PAS_median'] = p.median(axis=1)
        universe[f'{side}_pas_q3'] = qbin(universe[f'{side}_PAS_median'], 3)
        selected = candidates[candidates.species.eq(species)].set_index('gene_id')
        direction = selected.direction.map({'negative': -1, 'positive': 1})
        universe[f'{side}_label'] = ids.map(direction).fillna(0).astype(int)
        universe[f'{side}_symbol'] = ids.map(selected.gene_symbol)
        for name in ['S', 'M']:
            assert universe.loc[universe[f'{side}_label'].eq(-1), f'{side}_{name}_prox_capacity'].eq(len(ROUTES[species])).all()
            assert universe.loc[universe[f'{side}_label'].eq(1), f'{side}_{name}_dist_capacity'].eq(len(ROUTES[species])).all()
    universe.insert(0, 'pas_source', source)
    attrition = []
    for species, side in [('human', 'h'), ('mouse', 'm')]:
        selected = candidates[candidates.species.eq(species)]
        included = universe[f'{side}_label'].ne(0)
        attrition.append(dict(pas_source=source, species=species, all_recurrent=len(recurrent[recurrent.species.eq(species)]),
            original_dual=len(selected), dual_in_joint_ortholog=int(included.sum()),
            dual_outside_joint_ortholog=int(len(selected)-included.sum()),
            proximal_in_joint_ortholog=int(universe[f'{side}_label'].eq(-1).sum()),
            distal_in_joint_ortholog=int(universe[f'{side}_label'].eq(1).sum()),
            all_species_complete_triplet=len(eligible[species]), joint_ortholog=len(universe)))
    if source == 'polyapipe':
        contract = ROOT / 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables'
        for rel in ['matched_pms/matched_state_vectors.tsv.gz', 'candidate_screen/dual_reference_candidates.tsv',
                    'candidate_screen/dual_reference_ladder.tsv', 'cross_overlap/ortholog_matched_features.tsv.gz']:
            input_paths[str(contract / rel)] = digest(contract / rel)
            assert input_paths[str(contract / rel)] == input_paths[str(base / rel)], rel
    return universe, attrition


def stratum_info(frame, labels, side, conditional=False):
    cols = [f'{side}_{x}' for x in ['S_prox_capacity', 'M_prox_capacity', 'S_dist_capacity', 'M_dist_capacity', 'ref_q5', 'pas_q3']]
    base_groups = groups(frame, cols)
    use_groups = [ii[labels[ii] != 0] for ii in base_groups] if conditional else base_groups
    use_groups = [ii for ii in use_groups if len(ii)]
    selected_prob = np.zeros(len(labels), float)
    prox_prob = np.zeros(len(labels), float)
    rows = []
    for j, ii in enumerate(use_groups):
        n, prox, dist = len(ii), int((labels[ii] == -1).sum()), int((labels[ii] == 1).sum())
        selected_prob[ii] = (prox + dist)/n
        prox_prob[ii] = prox/n
        rows.append(dict(pas_source=frame.pas_source.iloc[0], side=side,
            mode='direction_conditional' if conditional else 'membership', stratum=j, n=n,
            proximal=prox, distal=dist, nonmember=n-prox-dist,
            stratum_signature='|'.join(str(v) for v in frame.loc[ii[0], cols]),
            membership_fixed=(prox+dist in [0,n]), direction_fixed=(prox in [0,n])))
    return use_groups, selected_prob, prox_prob, rows


def null_draws(x, y, hg, mg, n, rng, direction=False):
    out = []
    for start in range(0, n, 500):
        k = min(500, n-start)
        h = draw_labels(x, hg, k, rng)
        m = draw_labels(y, mg, k, rng)
        hits = (h != 0) & (m != 0)
        out.append((hits & (h == m)).sum(axis=1) if direction else hits.sum(axis=1))
    return np.concatenate(out)


def overlap_injection(x, y, hg, mg, n, odds, rng):
    out = []
    for start in range(0, n, 500):
        k = min(500, n-start)
        h = draw_labels(x, hg, k, rng) != 0
        overlap = np.zeros(k, int)
        for ii in mg:
            N, draws = len(ii), int((y[ii] != 0).sum())
            K = h[:, ii].sum(axis=1)
            for val in np.unique(K):
                use = K == val
                low, high = max(0, draws+int(val)-N), min(draws, int(val))
                if low == high:
                    overlap[use] += low
                else:
                    dist = hypergeom if odds == 1 else nchypergeom_fisher
                    args = (N, int(val), draws) if odds == 1 else (N, int(val), draws, odds)
                    overlap[use] += dist.rvs(*args, size=int(use.sum()), random_state=rng)
        out.append(overlap)
    return np.concatenate(out)


def tilted_direction_distribution(A, D, C, prox, odds):
    counts, logw = [], []
    for a in range(max(0, prox-D-C), min(A, prox)+1):
        for d in range(max(0, prox-a-C), min(D, prox-a)+1):
            c = prox-a-d
            if not 0 <= c <= C:
                continue
            counts.append(a + D-d)
            logw.append(math.log(math.comb(A,a)) + math.log(math.comb(D,d)) +
                        math.log(math.comb(C,c)) + (a-d)*math.log(odds))
    logw = np.asarray(logw)
    weights = np.exp(logw-logw.max())
    weights /= weights.sum()
    return np.asarray(counts), weights


def direction_injection(x, y, hg, mg, n, odds, rng):
    out, cache = [], {}
    for start in range(0, n, 500):
        k = min(500, n-start)
        h = draw_labels(x, hg, k, rng)
        same = np.zeros(k, int)
        for ii in mg:
            prox = int((y[ii] == -1).sum())
            A = (h[:,ii] == -1).sum(axis=1)
            D = (h[:,ii] == 1).sum(axis=1)
            pairs = np.column_stack([A,D])
            for a,d in np.unique(pairs, axis=0):
                key = int(a), int(d), len(ii)-int(a)-int(d), prox, odds
                if key not in cache:
                    cache[key] = tilted_direction_distribution(*key)
                values, weights = cache[key]
                use = (A == a) & (D == d)
                same[use] += rng.choice(values, size=int(use.sum()), p=weights)
        out.append(same)
    return np.concatenate(out)


def analyze(frame):
    source = frame.pas_source.iloc[0]
    x, y = frame.h_label.to_numpy(np.int8), frame.m_label.to_numpy(np.int8)
    hg, hp, _, hs = stratum_info(frame, x, 'h')
    mg, mp, _, ms = stratum_info(frame, y, 'm')
    hdg, _, hn, hd = stratum_info(frame, x, 'h', True)
    mdg, _, mn, md = stratum_info(frame, y, 'm', True)
    hits = (x != 0) & (y != 0)
    diagnostics, tests, powers, seeds, archives = hs+ms+hd+md, [], [], [], {}
    for test, human_groups, mouse_groups, expected in [
        ('overlap_enrichment', hg, mg, float(hp@mp)),
        ('direction_concordance', hdg, mdg, float((hn[hits]*mn[hits] + (1-hn[hits])*(1-mn[hits])).sum()))]:
        direction = test == 'direction_concordance'
        obs = int((hits & (x == y)).sum()) if direction else int(hits.sum())
        rng, sd = random(f'{source}|{test}|null')
        seeds.append(dict(pas_source=source, test=test, purpose='null', odds_ratio=1, seed=sd))
        null = null_draws(x, y, human_groups, mouse_groups, B, rng, direction)
        archives[f'{source}|{test}|null'] = null
        assert abs(null.mean()-expected) <= max(.02, 6*null.std(ddof=1)/np.sqrt(B))
        testrow = dict(family=FAMILY, pas_source=source, test=test, status='complete',
            universe=len(frame), human_members=int((x != 0).sum()), mouse_members=int((y != 0).sum()),
            overlap=int(hits.sum()), same_direction=int((hits & (x == y)).sum()),
            observed=obs, analytic_expected=expected, null_mean=float(null.mean()),
            null_sd=float(null.std(ddof=1)), null_q025=float(np.quantile(null,.025)),
            null_q975=float(np.quantile(null,.975)), enrichment=obs/expected if expected else np.nan,
            raw_p=float(probability(null, obs)) if hits.any() or not direction else np.nan,
            permutations=B, minimum_mc_p=1/(B+1), null_min=int(null.min()), null_max=int(null.max()),
            h_membership_movable=int(((x != 0) & (hp > 0) & (hp < 1)).sum()),
            m_membership_movable=int(((y != 0) & (mp > 0) & (mp < 1)).sum()),
            h_overlap_direction_movable=int((hits & (hn > 0) & (hn < 1)).sum()),
            m_overlap_direction_movable=int((hits & (mn > 0) & (mn < 1)).sum()))
        if direction and not hits.any():
            testrow['status'] = 'no_shared_members_unestimable'
        tests.append(testrow)
        print(json.dumps(testrow, ensure_ascii=False), flush=True)
        simulate = direction_injection if direction else overlap_injection
        for odds in ODDS:
            rng, sd = random(f'{source}|{test}|injection|{odds}')
            seeds.append(dict(pas_source=source, test=test, purpose='power', odds_ratio=odds, seed=sd))
            sim = simulate(x, y, human_groups, mouse_groups, POWER_N, odds, rng)
            archives[f'{source}|{test}|OR{odds}'] = sim
            p = probability(null, sim)
            for alpha in ALPHAS:
                powers.append(dict(pas_source=source, test=test, odds_ratio=odds, alpha=alpha,
                    fixed_universe=len(frame), fixed_human_members=int((x != 0).sum()),
                    fixed_mouse_members=int((y != 0).sum()), fixed_observed_overlap=int(hits.sum()),
                    simulated_mean=float(sim.mean()), simulated_sd=float(sim.std(ddof=1)),
                    expected_fraction=sim.mean()/int(hits.sum()) if direction and hits.any() else np.nan,
                    **rate(p <= alpha)))
    return tests, powers, diagnostics, seeds, archives


def save(frame, path):
    frame.to_csv(path, sep='\t', index=False, float_format='%.17g')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    inputs, universes, attrition = {}, {}, []
    for source in ['polyapipe', 'apadetect-p50']:
        f, rows = build(source, inputs)
        universes[source] = f
        attrition.extend(rows)
        save(f, args.out / (source + '-universe.tsv.gz'))
        hits = f[f.h_label.ne(0) & f.m_label.ne(0)].copy()
        hits['same_direction'] = hits.h_label.eq(hits.m_label)
        save(hits, args.out / (source + '-shared-members.tsv'))
    audit = dict(status='inputs_frozen_before_inference', inputs_sha256=inputs,
        code_sha256=digest(Path(__file__)), protocol_sha256=digest(HERE/'PROTOCOL.zh.md'),
        python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        platform=platform.platform(), master_seed=SEED, null_draws=B, power_draws=POWER_N,
        PAS_roles={'polyapipe':'current_manuscript_primary', 'apadetect-p50':'alternative_validation'})
    (args.out/'input-audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    common = set(universes['polyapipe'].human_gene_id) & set(universes['apadetect-p50'].human_gene_id)
    common_rows = []
    for source, f in universes.items():
        q = f[f.human_gene_id.isin(common)]
        common_rows.append(dict(pas_source=source, scope='cross_PAS_common_measurable_descriptive_only',
            universe=len(q), human_dual=int(q.h_label.ne(0).sum()), mouse_dual=int(q.m_label.ne(0).sum()),
            shared=int((q.h_label.ne(0)&q.m_label.ne(0)).sum()),
            same_direction=int((q.h_label.ne(0)&q.m_label.ne(0)&q.h_label.eq(q.m_label)).sum())))
    save(pd.DataFrame(common_rows), args.out/'cross-PAS-common-background.tsv')
    save(pd.DataFrame(attrition), args.out/'candidate-attrition.tsv')
    all_tests, all_powers, all_strata, all_seeds, archives = [], [], [], [], {}
    for f in universes.values():
        tests, powers, strata, seeds, draws = analyze(f)
        all_tests.extend(tests); all_powers.extend(powers); all_strata.extend(strata); all_seeds.extend(seeds)
        archives.update(draws)
    tests, powers = pd.DataFrame(all_tests), pd.DataFrame(all_powers)
    save(tests, args.out/'raw-pvalues.tsv')
    save(powers, args.out/'power-curves.tsv')
    save(pd.DataFrame(all_strata), args.out/'stratum-diagnostics.tsv')
    save(pd.DataFrame(all_seeds), args.out/'seeds.tsv')
    detectable = []
    for (source, test, alpha), q in powers.groupby(['pas_source','test','alpha']):
        passing = q[q.power.ge(.8)].sort_values('odds_ratio')
        r = passing.iloc[0] if len(passing) else q.sort_values('odds_ratio').iloc[-1]
        detectable.append(dict(pas_source=source, test=test, alpha=alpha,
            detectable_80_percent_in_grid=bool(len(passing)), first_OR_80=r.odds_ratio if len(passing) else np.nan,
            maximum_OR_tested=max(ODDS), reported_effect_OR=r.odds_ratio, expected_statistic=r.simulated_mean,
            expected_direction_fraction=r.expected_fraction, power_at_reported_effect=r.power,
            power_ci_low=r.power_ci_low, power_ci_high=r.power_ci_high))
    save(pd.DataFrame(detectable), args.out/'detectable-effects.tsv')
    np.savez_compressed(args.out/'simulation-draws.npz', **archives)
    outputs = {str(p.relative_to(args.out)): digest(p) for p in args.out.iterdir() if p.is_file()}
    result = dict(status='complete', tests=4, tests_skipped=0, family=FAMILY,
        checks=dict(one_to_one_ID_mapping=True, full_joint_triplet_background=True,
            manuscript_polyapipe_SHA_replayed=True, fixed_candidate_labels_replayed=True,
            reference_capacity_eligibility=True, independent_overlap_direction_tests=True,
            analytic_expectations_vs_MC=True, conditional_power_margin_preservation=True),
        outputs_sha256=outputs)
    (args.out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
