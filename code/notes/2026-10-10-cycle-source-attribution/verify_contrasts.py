"""Independent contributor arithmetic and small-sample model checks."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FOLDER = HERE / 'contrast-dose'
ROOT = HERE.parents[1]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'


def read(path, **kwargs):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kwargs)


def equal(a, b, tolerance=1e-11):
    assert np.allclose(a, b, atol=tolerance, rtol=0, equal_nan=True)


def main():
    scores = read(HERE / 'independent_mc_scores.tsv.gz').set_index('metacell_id')
    pairs = read(FOLDER / 'matched_pair_cycles.tsv.gz')
    key = ['species', 'study_id', 'route', 'contrast', 'design', 'gene_id']
    frozen = read(TABLES / 'state_vectors/pair_effects.primary.tsv.gz', usecols=key + ['setting', 'library_key_a', 'library_key_b', 'metacell_ids_a', 'metacell_ids_b', 'delta_DWUI', 'weight'])
    frozen = frozen[frozen.setting.eq('supported50_p08')].drop(columns='setting')
    matched = pairs.merge(frozen, on=key + ['library_key_a', 'library_key_b'], suffixes=('_new', '_frozen'), validate='one_to_one', how='outer', indicator=True)
    assert matched._merge.eq('both').all()
    for col in ['metacell_ids_a', 'metacell_ids_b']:
        assert matched[col + '_new'].eq(matched[col + '_frozen']).all()
    equal(matched.delta_DWUI_new, matched.delta_DWUI_frozen)
    equal(matched.weight_new, matched.weight_frozen)
    full = scores.groupby(['study_id', 'route', 'library_key', 'node'])[['Cycle', 'Adjusted_Cycle']].mean()
    cache = {}
    for side in 'ab':
        expected = []
        cohort = []
        for row in pairs.itertuples():
            ids = getattr(row, 'metacell_ids_' + side)
            if ids not in cache:
                members = scores.loc[json.loads(ids)]
                cache[ids] = members[['Cycle', 'Adjusted_Cycle']].mean().to_numpy()
            expected.append(cache[ids])
            node = 'P' if side == 'a' else row.contrast[-1]
            cohort.append(full.loc[(row.study_id, row.route, getattr(row, 'library_key_' + side), node)].to_numpy())
        expected, cohort = np.array(expected), np.array(cohort)
        for i, name in enumerate(['Cycle', 'Adjusted_Cycle']):
            equal(expected[:, i], pairs['exact_' + name + '_' + side])
            equal(cohort[:, i], pairs['cohort_' + name + '_' + side])
    variants = ['exact_Cycle', 'exact_Adjusted_Cycle', 'cohort_Cycle', 'cohort_Adjusted_Cycle']
    for variant in variants:
        equal(pairs[variant + '_a'] - pairs[variant + '_b'], pairs[variant + '_delta'])
    for keys, q in pairs.groupby(key, sort=False):
        # Count children at each level rather than borrowing the analysis weights function.
        target_dates = q.time_a.nunique()
        targets = q.groupby('time_a').sample_unit_a.nunique()
        sources = q.groupby(['time_a', 'sample_unit_a']).time_b.nunique()
        units = q.groupby(['time_a', 'sample_unit_a', 'time_b']).sample_unit_b.nunique()
        captures = q.groupby(['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b']).size()
        w = np.array([1. / (target_dates * targets[row.time_a] * sources[row.time_a, row.sample_unit_a] *
                            units[row.time_a, row.sample_unit_a, row.time_b] *
                            captures[row.time_a, row.sample_unit_a, row.time_b, row.sample_unit_b]) for row in q.itertuples()])
        equal(w, q.weight)
        equal(w.sum(), 1.)
    genes = read(FOLDER / 'route_gene_cycles.tsv.gz').set_index(key)
    for keys, q in pairs.groupby(key, sort=False):
        row = genes.loc[keys]
        equal(np.dot(q.weight, q.delta_DWUI), row.delta_DWUI)
        for variant in variants:
            equal(np.dot(q.weight, q[variant + '_delta']), row[variant + '_delta'])
    summary = read(FOLDER / 'route_contrast_cycle_summary.tsv')
    for row in summary.itertuples():
        q = genes.reset_index()
        q = q[q.route.eq(row.route) & q.contrast.eq(row.contrast)]
        if row.background == 'common_ortholog':
            orthologs = read(FOLDER / ('common_' + row.contrast + '_orthologs.tsv'))
            q = q[q.gene_id.isin(orthologs[row.species + '_gene_id'])]
        assert len(q) == row.genes
        equal(q.delta_DWUI.mean(), row.mean_DWUI)
        equal(q[row.cycle_variant + '_delta'].mean(), row.mean_Cycle_delta)
    models = read(FOLDER / 'source_PS_descriptive_models.tsv')
    for row in models.itertuples():
        q = summary[summary.contrast.eq('P_S') & summary.background.eq(row.background) & summary.cycle_variant.eq(row.cycle_variant)].sort_values('route')
        x, y = q.mean_Cycle_delta.to_numpy(), q.mean_DWUI.to_numpy()
        w = np.where(q.study_id.eq('GSE178325'), .5, 1.)
        mouse = q.species.eq('mouse').to_numpy(float)
        X = np.column_stack([np.ones(len(q)), x])
        if row.model == 'cycle_plus_mouse':
            X = np.column_stack([X, mouse])
        beta = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * y))
        equal(beta[:2], [row.intercept, row.cycle_beta], tolerance=1e-9)
        if row.model == 'cycle_plus_mouse':
            equal(beta[2], row.mouse_beta, tolerance=1e-9)
            # Frisch-Waugh-Lovell independently isolates the numerical mouse coefficient.
            C = X[:, :2]
            rx = mouse - C @ np.linalg.solve(C.T @ (w[:, None] * C), C.T @ (w * mouse))
            ry = y - C @ np.linalg.solve(C.T @ (w[:, None] * C), C.T @ (w * y))
            equal(np.dot(w * rx, ry) / np.dot(w * rx, rx), row.mouse_beta, tolerance=1e-9)
        err = y - X @ beta
        equal(1 - np.dot(w, err ** 2) / np.dot(w, (y - np.average(y, weights=w)) ** 2), row.weighted_R2)
    triangles = read(FOLDER / 'matched_triangle_gene_cycles.tsv.gz')
    for variant in variants:
        equal(triangles[variant + '_PS'], triangles[variant + '_MS'] + triangles[variant + '_PM'])
    equal(triangles.DWUI_PS, triangles.DWUI_MS + triangles.DWUI_PM)
    assert models.p.isna().all()
    result = dict(status='PASS', independent_raw_cycle=True, pairs_checked=len(pairs), gene_route_rows_checked=len(genes),
                  route_summary_rows_checked=len(summary), numerical_models_checked=len(models),
                  triangle_gene_routes_checked=len(triangles), no_species_inference_p=True, tests_skipped=0)
    (HERE / 'contrast_independent_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
