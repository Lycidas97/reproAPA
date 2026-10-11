"""Independent delivery verification from saved vectors, scores and count samples."""
from collections import defaultdict
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / 'notes/2026-10-10-marker-reference-sensitivity'
VALUE_COLS = ['delta_DWUI', 'delta_Cycle_exact', 'delta_Adjusted_Cycle_exact',
              'delta_Cycle_fullcohort', 'delta_Adjusted_Cycle_fullcohort']


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def independent_weights(q):
    # Count physical captures and all child units within each fixed parent.
    keys = list(q[['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b']].itertuples(index=False, name=None))
    captures = defaultdict(int)
    children = [defaultdict(set) for _ in range(4)]
    for key in keys:
        captures[key] += 1
        for level in range(4):
            children[level][key[:level]].add(key[:level + 1])
    w = np.array([1 / captures[key] / np.prod([len(children[level][key[:level]]) for level in range(4)]) for key in keys])
    assert np.isclose(w.sum(), 1, rtol=0, atol=1e-14)
    return w


def validate_pairs(pairs, cycles, z):
    maxima = dict(vector=0., exact_cycle=0., fullcohort_cycle=0., axis=0., raw_usage=0.)
    profiles, full_means, axis, qualified_means = {}, {}, {}, {}
    cycle_records = {key: (row.Cycle, row.Adjusted_Cycle, row.n_cells, row.node, row.library_key)
                     for key, row in cycles.iterrows()}
    for sid in pairs.study_id.unique():
        feature = pd.DataFrame(json.loads(str(z[sid + '_features_json'])))
        fullfeature = feature.copy()
        feature = feature[feature.coordinate_eligible & feature.matched_source_count_50nt.ge(1)]
        for gene, part in feature.groupby('primary_gene_id', sort=False):
            axis[(sid, gene)] = (part.atlas_pas_id.tolist(), part.position.to_numpy(float), part.index.to_numpy())
        selected = pd.DataFrame(json.loads(str(z[sid + '_selected_json'])))
        included_routes = set(pairs.loc[pairs.study_id.eq(sid), 'route'])
        selected = selected[selected.route.isin(included_routes) & selected.purity.ge(.8) & selected.node.isin(['S', 'M', 'P'])]
        assert selected.metacell_id.isin(cycles.index).all()
        for (route, lib, state), q in selected.groupby(['route', 'library_key', 'node']):
            full_means[(sid, route, lib, state)] = cycles.loc[q.metacell_id, ['Cycle', 'Adjusted_Cycle']].mean().to_numpy()
    for row in pairs.itertuples(index=False):
        pos = np.array(json.loads(row.position))
        ids = json.loads(row.PAS_ids)
        wanted_ids, wanted_pos, columns = axis[(row.study_id, row.gene_id)]
        assert ids == wanted_ids
        maxima['axis'] = max(maxima['axis'], float(np.max(np.abs(pos - wanted_pos))))
        ua, ub = np.array(json.loads(row.usage_a)), np.array(json.loads(row.usage_b))
        assert len(ua) == len(ub) == len(pos)
        assert min(ua.min(), ub.min()) >= 0 and max(ua.max(), ub.max()) <= 1 + 1e-14
        assert np.isclose(ua.sum(), 1, rtol=0, atol=1e-14) and np.isclose(ub.sum(), 1, rtol=0, atol=1e-14)
        sign = -1 if row.contrast == 'S_M' else 1
        maxima['vector'] = max(maxima['vector'], abs(row.delta_DWUI - sign * np.dot(ua - ub, pos)))
        ra, rb = json.loads(row.qualified_MC_ids_a), json.loads(row.qualified_MC_ids_b)
        assert not set(ra) & set(rb), 'Original state cohorts must be disjoint'
        means = {}
        for side, rr in [('a', ra), ('b', rb)]:
            state = getattr(row, 'target_state' if side == 'a' else 'reference_state')
            lib = getattr(row, 'library_key_' + side)
            key = (row.study_id, row.route, lib, state, row.gene_id)
            assert len(rr) == getattr(row, 'n_qualified_MC_' + side) >= 3
            if key not in qualified_means:
                values = [cycle_records[mc] for mc in rr]
                assert all(value[3] == state and value[4] == lib for value in values)
                cells = sum(value[2] for value in values)
                qualified_means[key] = (np.array([value[:2] for value in values]).mean(0), cells)
            mean, cells = qualified_means[key]
            assert int(cells) == getattr(row, 'n_qualified_cells_' + side) >= 30
            means[side] = mean
        am, bm = means['a'], means['b']
        calculated = sign * (am - bm)
        actual = np.array([row.delta_Cycle_exact, row.delta_Adjusted_Cycle_exact])
        maxima['exact_cycle'] = max(maxima['exact_cycle'], float(np.max(np.abs(actual - calculated))))
        ac = full_means[(row.study_id, row.route, row.library_key_a, row.target_state)]
        bc = full_means[(row.study_id, row.route, row.library_key_b, row.reference_state)]
        actual = np.array([row.delta_Cycle_fullcohort, row.delta_Adjusted_Cycle_fullcohort])
        maxima['fullcohort_cycle'] = max(maxima['fullcohort_cycle'], float(np.max(np.abs(actual - sign * (ac - bc)))))
        for side, rr in [('a', ra), ('b', rb)]:
            key = (row.study_id, row.route, getattr(row, 'library_key_' + side),
                   getattr(row, 'target_state' if side == 'a' else 'reference_state'), row.gene_id)
            value = (tuple(rr), tuple(np.asarray(json.loads(getattr(row, 'usage_' + side)))),
                     tuple(json.loads(getattr(row, 'pooled_molecules_' + side))),
                     tuple(json.loads(getattr(row, 'nonzero_MC_' + side))))
            if key in profiles:
                assert profiles[key] == value
            else:
                profiles[key] = value
        ma, mb = np.array(json.loads(row.pooled_molecules_a)), np.array(json.loads(row.pooled_molecules_b))
        na, nb = np.array(json.loads(row.nonzero_MC_a)), np.array(json.loads(row.nonzero_MC_b))
        assert int(((ma + mb >= 10) & (na + nb >= 2)).sum()) == row.joint_supported_PAS >= 2
    # Deterministic sample spans every study/state and thin/fuller cohorts.
    candidates = list(profiles)
    samples = []
    for sid, state in sorted(set((key[0], key[3]) for key in candidates)):
        part = [key for key in candidates if key[0] == sid and key[3] == state]
        part.sort(key=lambda key: hashlib.sha256(repr(key).encode()).digest())
        samples.extend(part[:8])
    for sid in sorted(set(key[0] for key in samples)):
        meta = pd.DataFrame(json.loads(str(z[sid + '_meta_json'])))
        selected = pd.DataFrame(json.loads(str(z[sid + '_selected_json'])))
        selected = selected[selected.purity.ge(.8)]
        d, ix, ptr = (z[sid + suffix] for suffix in ['_pas_data', '_pas_indices', '_pas_indptr'])
        for key in [key for key in samples if key[0] == sid]:
            _, route, lib, state, gene = key
            ids, position, cols = axis[(sid, gene)]
            q = selected[selected.route.eq(route) & selected.library_key.eq(lib) & selected.node.eq(state)]
            rr = q.matrix_row.to_numpy()
            counts = np.zeros((len(meta), len(cols)), dtype=d.dtype)
            for j, col in enumerate(cols):
                start, end = ptr[col:col + 2]
                counts[ix[start:end], j] = d[start:end]
            keep = rr[counts[rr].sum(1) >= 20]
            assert tuple(meta.metacell_id.iloc[keep]) == profiles[key][0]
            c = counts[keep]
            usage = (c / c.sum(1)[:, None]).mean(0)
            maxima['raw_usage'] = max(maxima['raw_usage'], float(np.max(np.abs(usage - np.array(profiles[key][1])))))
            assert np.array_equal(c.sum(0), np.array(profiles[key][2]))
            assert np.array_equal((c > 0).sum(0), np.array(profiles[key][3]))
    assert max(maxima.values()) < 1e-12, maxima
    return maxima, len(profiles), len(samples)


def validate_aggregates(pairs, effects, by):
    source = {}
    for key, q in pairs.groupby(by, sort=False):
        if not isinstance(key, tuple):
            key = (key,)
        q = q.sort_values(['library_key_a', 'library_key_b']).reset_index(drop=True)
        source[key] = q
    max_delta, max_weight = 0., 0.
    for row in effects.to_dict('records'):
        key = tuple(row[field] for field in by)
        q = source[key]
        w = independent_weights(q)
        assert json.loads(row['contributor_pairs']) == [list(key) for key in zip(q.library_key_a, q.library_key_b)]
        assert int(row['n_pairs']) == len(q)
        max_weight = max(max_weight, float(np.max(np.abs(w - np.array(json.loads(row['weights']))))))
        actual = np.array([row[field] for field in VALUE_COLS])
        expected = w @ q[VALUE_COLS].to_numpy()
        max_delta = max(max_delta, float(np.max(np.abs(actual - expected))))
    assert max(max_delta, max_weight) < 1e-12
    return dict(n_aggregates=len(effects), max_delta_error=max_delta, max_weight_error=max_weight)


def main():
    files = [HERE / 'SM_pair_effects.tsv.gz', HERE / 'observed_PS_pair_effects.tsv.gz',
             HERE / 'SM_route_gene_effects.tsv.gz', HERE / 'observed_PS_route_effects.tsv.gz',
             HERE / 'SM_route_summary.tsv', HERE / 'SM_candidate_pairs.tsv', HERE / 'fixed_gene_sets.tsv',
             OLD / 'inputs.npz', ROOT / 'notes/2026-10-10-proliferation-apa/metacell_cycle_scores.tsv.gz',
             HERE.parent / 'PROTOCOL.zh.md', HERE / 'SM_library_date_summary.tsv']
    hashes = {str(path): digest(path) for path in files}
    sm, ps, smr, psr, summaries, candidates, sets = [pd.read_csv(path, sep='\t', float_precision='round_trip') for path in files[:7]]
    cycles = pd.read_csv(files[8], sep='\t', float_precision='round_trip').set_index('metacell_id')
    stable = pd.read_csv(files[10], sep='\t', float_precision='round_trip')
    z = np.load(files[7], allow_pickle=False)
    maxima, profiles, raw_samples = validate_pairs(pd.concat([sm, ps], ignore_index=True), cycles, z)
    checks = []
    stability_checks, stability_max_error = 0, 0.
    pscheck = validate_aggregates(ps, psr, ['species', 'study_id', 'route', 'gene_id', 'design', 'source_type'])
    for row in summaries.to_dict('records'):
        q = sm[sm.study_id.eq(row['study_id']) & sm.route.eq(row['route']) & sm.design.eq(row['design'])]
        wanted = set(sets.loc[sets.species.eq(row['species']) & sets.gene_set.eq(row['gene_set']), 'gene_id'])
        q = q[q.gene_id.isin(wanted)]
        n_candidate = len(candidates[candidates.study_id.eq(row['study_id']) & candidates.route.eq(row['route']) & candidates.design.eq(row['design'])])
        assert n_candidate == row['n_candidate_pairs']
        if row['background'] == 'all_candidate_pairs_intersection':
            counts = q.groupby('gene_id').size()
            common = set(counts.index[counts.eq(n_candidate)]) if n_candidate else set()
            q = q[q.gene_id.isin(common)]
        current = smr[smr.study_id.eq(row['study_id']) & smr.route.eq(row['route']) & smr.design.eq(row['design']) &
                      smr.gene_set.eq(row['gene_set']) & smr.background.eq(row['background'])]
        assert set(current.gene_id) == set(q.gene_id) and len(current) == row['N']
        if len(q):
            checks.append(validate_aggregates(q, current, ['gene_id']))
            assert abs(current.delta_DWUI.mean() - row['mean_DWUI']) < 1e-12
            assert abs(current.delta_DWUI.median() - row['median_DWUI']) < 1e-12
            assert int(current.delta_DWUI.le(-.02).sum()) == row['proximal']
            assert int(current.delta_DWUI.ge(.02).sum()) == row['distal']
            for field in VALUE_COLS[1:]:
                assert abs(current[field].mean() - row['mean_' + field]) < 1e-12
                assert abs(current[field].median() - row['median_' + field]) < 1e-12
        else:
            assert row['status'] == 'not_estimable' and pd.isna(row['mean_DWUI'])
        st = stable[stable.study_id.eq(row['study_id']) & stable.route.eq(row['route']) & stable.design.eq(row['design']) &
                    stable.gene_set.eq(row['gene_set']) & stable.background.eq(row['background'])]
        for s in st.to_dict('records'):
            field = 'library_key_a' if s['level'] == 'target_library' else 'time_a'
            current_pairs = q[q[field].eq(s['target_label'])]
            values = []
            for _, gene_pairs in current_pairs.groupby('gene_id', sort=False):
                values.append(independent_weights(gene_pairs) @ gene_pairs[VALUE_COLS].to_numpy())
            assert len(values) == s['N']
            if values:
                values = np.array(values)
                assert int((values[:, 0] <= -.02).sum()) == s['proximal']
                assert int((values[:, 0] >= .02).sum()) == s['distal']
                expected = [values[:, 0].mean(), np.median(values[:, 0])]
                actual = [s['mean_DWUI'], s['median_DWUI']]
                for j, name in enumerate(VALUE_COLS[1:], 1):
                    expected.extend([values[:, j].mean(), np.median(values[:, j])])
                    actual.extend([s['mean_' + name], s['median_' + name]])
                stability_max_error = max(stability_max_error, float(np.max(np.abs(np.array(actual) - np.array(expected)))))
                assert stability_max_error < 1e-12
            else:
                assert s['status'] == 'not_estimable' and pd.isna(s['mean_DWUI'])
            stability_checks += 1
    assert len(summaries) == 7 * 4 * 2 * 2
    assert stability_checks == len(stable) == 756
    assert all(digest(path) == expected for path, expected in hashes.items())
    validation = dict(status='PASS', input_sha256=hashes, SM_pair_rows=len(sm), observed_PS_pair_rows=len(ps),
                      all_saved_vectors_recomputed=True, all_qualified_MC_score_recomputed=True,
                      vector_and_score_max_errors=maxima, distinct_gene_state_profiles=profiles, raw_count_profiles_sampled=raw_samples,
                      SM_route_summary_rows=len(summaries), SM_aggregate_rows_checked=sum(row['n_aggregates'] for row in checks),
                      SM_max_delta_error=max(row['max_delta_error'] for row in checks),
                      SM_max_weight_error=max(row['max_weight_error'] for row in checks), observed_PS_aggregate_check=pscheck,
                      library_date_summary_rows_verified=stability_checks, library_date_max_error=stability_max_error,
                      fixed_intersections_and_NE_verified=True, no_gene_MC_species_p=True, verify_script_sha256=digest(Path(__file__)))
    (HERE / 'delivery_validation.json').write_text(json.dumps(validation, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(validation | {'input_sha256': 'see_file'}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
