"""Independent contribution weights, gene rosters, and signed summaries."""
from collections import Counter, defaultdict
from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FOLDER = HERE / 'sm-core'


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def equal(a, b):
    assert np.allclose(a, b, atol=1e-11, rtol=0, equal_nan=True)


def weights(q):
    dates = set(q.time_a)
    targets, sources, units = defaultdict(set), defaultdict(set), defaultdict(set)
    captures = Counter()
    for row in q.itertuples():
        targets[row.time_a].add(row.sample_unit_a)
        sources[row.time_a, row.sample_unit_a].add(row.time_b)
        units[row.time_a, row.sample_unit_a, row.time_b].add(row.sample_unit_b)
        captures[row.time_a, row.sample_unit_a, row.time_b, row.sample_unit_b] += 1
    return np.array([1. / (len(dates) * len(targets[r.time_a]) * len(sources[r.time_a, r.sample_unit_a]) *
                          len(units[r.time_a, r.sample_unit_a, r.time_b]) *
                          captures[r.time_a, r.sample_unit_a, r.time_b, r.sample_unit_b]) for r in q.itertuples()])


def main():
    pair = read(FOLDER / 'SM_pair_effects.tsv.gz')
    gene = read(FOLDER / 'SM_route_gene_effects.tsv.gz')
    groups = {key: q.sort_values(['library_key_a', 'library_key_b']) for key, q in pair.groupby(['study_id', 'route', 'design', 'gene_id'])}
    fields = ['delta_DWUI', 'delta_Cycle_exact', 'delta_Adjusted_Cycle_exact', 'delta_Cycle_fullcohort', 'delta_Adjusted_Cycle_fullcohort']
    cache = {}
    for row in gene.itertuples():
        key = row.study_id, row.route, row.design, row.gene_id
        q = groups[key]
        if key not in cache:
            w = weights(q)
            equal(w.sum(), 1.)
            cache[key] = w, w @ q[fields].to_numpy()
        w, values = cache[key]
        assert json.loads(row.contributor_pairs) == [list(x) for x in zip(q.library_key_a, q.library_key_b)]
        assert row.n_pairs == len(q)
        equal(w, json.loads(row.weights))
        equal(values, [getattr(row, field) for field in fields])
    summary = read(FOLDER / 'SM_route_summary.tsv')
    for row in summary.itertuples():
        q = gene[gene.route.eq(row.route) & gene.design.eq(row.design) & gene.background.eq(row.background) & gene.gene_set.eq(row.gene_set)]
        assert len(q) == row.N
        if not len(q):
            assert row.status == 'not_estimable' and np.isnan(row.mean_DWUI)
            continue
        equal(q.delta_DWUI.mean(), row.mean_DWUI)
        equal(q.delta_DWUI.median(), row.median_DWUI)
        assert int(q.delta_DWUI.le(-.02).sum()) == row.proximal
        assert int(q.delta_DWUI.ge(.02).sum()) == row.distal
        assert row.proximal + row.distal + row.below_threshold == row.N
        for field in fields[1:]:
            equal(q[field].mean(), getattr(row, 'mean_' + field))
        if row.background == 'all_candidate_pairs_intersection':
            assert q.n_pairs.eq(row.n_candidate_pairs).all()
    # Observed P-S is reconstructed separately and retains its observed label.
    observed = read(FOLDER / 'observed_PS_route_effects.tsv.gz')
    op = read(FOLDER / 'observed_PS_pair_effects.tsv.gz')
    og = {key: q.sort_values(['library_key_a', 'library_key_b']) for key, q in op.groupby(['study_id', 'route', 'gene_id'])}
    for row in observed.itertuples():
        q = og[row.study_id, row.route, row.gene_id]
        equal(weights(q) @ q[fields].to_numpy(), [getattr(row, field) for field in fields])
        assert row.source_type == 'observed_S'
    result = dict(status='PASS', SM_gene_route_rows=len(gene), SM_summary_rows=len(summary),
                  observed_PS_gene_route_rows=len(observed), independent_weights=True,
                  threshold_sign_and_unavailable_checked=True, tests_skipped=0)
    (HERE / 'sm_aggregate_independent_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
