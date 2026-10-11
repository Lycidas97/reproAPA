"""Independent date contrasts, fixed panels, and exact conditional statistics."""
from pathlib import Path
import itertools
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FOLDER = HERE / 'wot-dose'
ROOT = HERE.parents[1]


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def equal(a, b):
    assert np.allclose(a, b, atol=1e-11, rtol=0, equal_nan=True)


def coefficient(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return np.dot(x, y) / np.sqrt(np.dot(x, x) * np.dot(y, y)), np.dot(x, y) / np.dot(x, x)


def main():
    scores = read(HERE / 'independent_mc_scores.tsv.gz').set_index('metacell_id')
    genes = read(FOLDER / 'native_gene_effects_exact_MC_scores.tsv.gz')
    cache = {}
    for side in 'ab':
        expected = []
        for ids in genes['metacell_ids_' + side]:
            if ids not in cache:
                cache[ids] = scores.loc[json.loads(ids), ['Cycle', 'Adjusted_Cycle']].mean().to_numpy()
            expected.append(cache[ids])
        expected = np.array(expected)
        equal(expected[:, 0], genes['raw_Cycle_' + side])
        equal(expected[:, 1], genes['Adjusted_Cycle_' + side])
    equal(genes.raw_Cycle_a - genes.raw_Cycle_b, genes.delta_Cycle_exact_MC)
    libs = read(FOLDER / 'fixed_panel_library_summary.tsv')
    for row in libs.itertuples():
        panel = json.loads(row.gene_ids)
        q = genes[genes.library_key.eq(row.library_key) & genes.gene_id.isin(panel)]
        assert len(q) == row.n_genes == len(panel)
        equal(-q.delta_DWUI.mean(), row.proximal_mean)
        equal(-q.delta_DWUI.median(), row.proximal_median)
        equal(q.delta_Cycle_exact_MC.mean(), row.delta_Cycle_exact_MC)
        equal(q.delta_Adjusted_Cycle_exact_MC.mean(), row.delta_Adjusted_Cycle_exact_MC)
    for panel, q in libs.groupby('panel'):
        assert len(q) == 24 and q.gene_ids.nunique() == 1
    dates = read(FOLDER / 'fixed_panel_date_dose.tsv')
    for row in dates.itertuples():
        q = libs[libs.panel.eq(row.panel) & libs.time.eq(row.time)]
        assert len(q) == row.n_wells
        if len(q) == 1:
            assert not row.eligible_date and np.isnan(row.Y_mean)
            continue
        a, b = (q[q.C_label.eq(label)].iloc[0] for label in ['C1', 'C2'])
        equal(a.proximal_mean - b.proximal_mean, row.Y_mean)
        equal(a.proximal_median - b.proximal_median, row.Y_median)
        for field in ['delta_Cycle_exact_MC', 'delta_Adjusted_Cycle_exact_MC', 'delta_raw_Cycle_cohort', 'delta_Adjusted_Cycle_cohort']:
            equal(a[field] - b[field], getattr(row, 'X_' + field))
    stats = read(FOLDER / 'fixed_panel_descriptive_statistics.tsv')
    for row in stats.itertuples():
        q = dates[dates.panel.eq(row.panel) & dates.eligible_date]
        assert len(q) == row.n_dates == 11
        x, y = q['X_' + row.predictor].to_numpy(), q[row.outcome].to_numpy()
        if row.batch_constant_centered:
            x, y = x - x.mean(), y - y.mean()
        r, beta = coefficient(x, y)
        equal(r, row.r)
        equal(beta, row.beta)
    core = read(FOLDER / 'core_date_replay.tsv')
    core = core[core.eligible_primary_secondary]
    assert len(core) == 7 and core.n_common_genes.tolist() == [12, 5, 22, 14, 8, 5, 9]
    formal = read(FOLDER / 'core_formal_statistics_replay.tsv')
    assignments = read(FOLDER / 'core_exact_assignments_replay.tsv')
    recomputed_p = []
    for row in formal.itertuples():
        field = 'X_Cycle' if row.predictor == 'Cycle' else 'X_CPA'
        # Accept the old predictor's literal capitalization as a saved label only.
        if field not in core:
            raise AssertionError(field)
        x, y = core[field].to_numpy(), core.Y_mean.to_numpy()
        r, beta = coefficient(x, y)
        equal(r, row.r)
        equal(beta, row.beta)
        observed = assignments[assignments.predictor.eq(row.predictor)].sort_values('assignment')
        assert len(observed) == 128
        null = []
        for signs in itertools.product([-1, 1], repeat=7):
            null.append(coefficient(x * np.array(signs), y)[0])
        equal(null, observed.r)
        p = np.mean(np.abs(null) >= abs(r) - 1e-12)
        equal(p, row.p_exact_two_sided)
        recomputed_p.append(p)
    p = np.array(recomputed_p)
    order = np.argsort(p)
    q = np.empty(2)
    q[order] = np.minimum(1, np.minimum.accumulate((2 * p[order] / np.arange(1, 3))[::-1])[::-1])
    equal(q, formal.q2_BH)
    # Reconstruct the core date predictors from independent raw cycle and CPA counts.
    independent = scores[scores.study_id.eq('GSE115943')].reset_index()
    independent['short_library'] = independent.library_key.str.extract('(GSM[0-9]+)')
    group = independent.groupby(['short_library', 'node'])[['Cycle', 'Adjusted_Cycle']].mean()
    saved = read(FOLDER / 'raw_reconstructed_library_predictors.tsv')
    for row in saved.itertuples():
        equal(group.loc[(row.library_key, 'P'), 'Cycle'] - group.loc[(row.library_key, 'M'), 'Cycle'], row.delta_raw_Cycle_cohort)
    archive = np.load(ROOT / 'notes/2026-10-10-cpa-factor-expression/expression_inputs.npz', allow_pickle=False)
    meta = pd.DataFrame(json.loads(str(archive['GSE115943_meta_json'])))
    factors = json.loads((ROOT / 'notes/2026-10-10-cpa-factor-expression/protocol.json').read_text())['primary_factors']
    cols = [np.flatnonzero(archive['GSE115943_factors'] == factor).item() for factor in factors]
    values = np.log2(1 + 1000000. * archive['GSE115943_counts'][:, cols] / archive['GSE115943_total_gex'][:, None])
    for row in saved.itertuples():
        means = []
        for state in ['P', 'M']:
            ids = np.flatnonzero(meta.library_key.eq(row.library_key) & meta.node.eq(state) & meta.purity.ge(.8))
            assert len(ids) >= 3 and meta.n_cells.iloc[ids].sum() >= 30
            means.append(values[ids].mean(0))
        equal((means[0] - means[1]).mean(), row.delta_CPA_cohort)
    result = dict(status='PASS', gene_library_cycles_checked=len(genes), fixed_panel_library_rows=len(libs),
                  fixed_panel_statistic_rows=len(stats), all24_listed_only22_date_informative=True,
                  exact_assignments_checked=len(assignments), raw_cycle_and_CPA_reconstructed=True,
                  original_family2_replayed=True, new_p_values=0, tests_skipped=0)
    (HERE / 'wot_independent_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
