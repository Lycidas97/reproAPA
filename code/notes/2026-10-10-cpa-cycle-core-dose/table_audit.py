"""Independently recover existing cycle/CPA summaries; no APA dose analysis."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CYCLE = ROOT / 'notes/2026-10-10-proliferation-apa'
CPA = ROOT / 'notes/2026-10-10-cpa-factor-expression'
NAIVE = ROOT / 'notes/2026-10-10-naive-gex-recovery'
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
OLD_PROJECT = Path('/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot')
INPUTS = [CYCLE / 'inputs.npz', CYCLE / 'paired_library_score_differences.tsv',
          CPA / 'expression_inputs.npz', CPA / 'protocol.json',
          CPA / 'paired_library_effects.tsv', CPA / 'expression_contrasts.tsv',
          NAIVE / 'raw/metacell_by_gene_counts.npz', NAIVE / 'raw/metacells.tsv.gz',
          NAIVE / 'raw/genes.tsv.gz', TABLES / 'naive_state/frozen_selected_metadata.tsv.gz',
          OLD_PROJECT / 'analysis/tandem-apa-v1/library_replicate_design.tsv',
          OLD_PROJECT / 'outputs/figure3-apa-gex-state-analysis-2026-09-19/sampling_design.tsv']
KEY = ['species', 'study_id', 'route']
checks = []


def check(name, condition, **detail):
    checks.append(dict(check=name, passed=bool(condition), **detail))
    if not condition:
        raise ValueError(name)


def tab(path):
    return pd.read_csv(path, sep='\t')


def save(frame, filename):
    frame.to_csv(HERE / filename, sep='\t', index=False, float_format='%.17g', na_rep='NA')


def state_means(meta, values, names):
    pieces = []
    for (route, library, state), z in meta.groupby(['route', 'library_key', 'node']):
        ix = z.index.to_numpy()
        if len(z) < 3 or z.n_cells.sum() < 30:
            continue
        q = pd.DataFrame({'factor': names, 'value': values[ix].mean(0)})
        for key, value in dict(species=z.species.iloc[0], study_id=z.study_id.iloc[0],
                               route=route, library_key=library, state=state,
                               time=str(z.experimental_time.iloc[0])).items():
            q[key] = value
        pieces.append(q)
    return pd.concat(pieces, ignore_index=True)


def pair_effects(states):
    keys = KEY + ['factor', 'library_key', 'time']
    p = states[states.state.eq('P')]
    m = states[states.state.eq('M')]
    q = p.merge(m, on=keys, suffixes=('_P', '_M'), validate='one_to_one')
    q['delta'] = q.value_P - q.value_M
    return q[keys + ['delta']]


def main():
    factors = json.loads((CPA / 'protocol.json').read_text())['primary_factors']
    ca = np.load(CYCLE / 'inputs.npz', allow_pickle=False)
    cycle_states = []
    for key in ca.files:
        if not key.endswith('_marker_counts'):
            continue
        sid = key[:-14]
        meta = pd.DataFrame(json.loads(str(ca[sid + '_meta_json'])))
        meta['study_id'] = sid
        counts = ca[key]
        check(sid + '_cycle_marker_integer_counts', np.issubdtype(counts.dtype, np.integer) and (counts >= 0).all())
        raw = np.log1p(counts / ca[sid + '_total_gex'][:, None] * 1e4)
        check(sid + '_18_fixed_markers', raw.shape[1] == 18)
        cycle = (raw[:, :10].mean(1) + raw[:, 10:].mean(1)) / 2
        check(sid + '_cycle_raw_recovery', np.allclose(cycle, meta.Cycle, atol=1e-13, rtol=1e-13),
              max_abs_error=float(np.max(np.abs(cycle - meta.Cycle))))
        cycle_states.append(state_means(meta, cycle[:, None], ['Cycle']))
    cycle_pairs = pair_effects(pd.concat(cycle_states, ignore_index=True))
    frozen_cycle = tab(CYCLE / 'paired_library_score_differences.tsv')
    frozen_cycle = frozen_cycle[frozen_cycle.contrast.eq('P_M') & frozen_cycle.score.eq('Cycle')]
    joined = cycle_pairs.merge(frozen_cycle, on=KEY + ['library_key', 'time'], suffixes=('_new', '_old'), validate='one_to_one')
    check('cycle_pair_inventory_recovered', len(joined) == len(cycle_pairs) == len(frozen_cycle))
    check('cycle_pair_effects_recovered', np.allclose(joined.delta_new, joined.delta_old, atol=1e-13, rtol=1e-13))

    cp = np.load(CPA / 'expression_inputs.npz', allow_pickle=False)
    expression_states = []
    for sid in ['GSE103221', 'GSE115943', 'GSE118258', 'GSE178325', 'GSE221739']:
        meta = pd.DataFrame(json.loads(str(cp[sid + '_meta_json'])))
        names = cp[sid + '_factors'].astype(str)
        columns = np.array([np.flatnonzero(names == factor).item() for factor in factors])
        counts = cp[sid + '_counts'][:, columns]
        check(sid + '_CPA_all_11_mapped', (counts >= 0).all() and np.issubdtype(counts.dtype, np.integer))
        expression_states.append(state_means(meta, np.log2(1 + counts / cp[sid + '_total_gex'][:, None] * 1e6), factors))
    meta = tab(TABLES / 'naive_state/frozen_selected_metadata.tsv.gz')
    meta = meta[meta.purity.ge(.8)].reset_index(drop=True)
    gmeta = tab(NAIVE / 'raw/metacells.tsv.gz')
    genes = tab(NAIVE / 'raw/genes.tsv.gz')
    rows = pd.Index(gmeta.metacell_id).get_indexer(meta.metacell_id)
    check('Liu_raw_MC_ID_unique_and_joined', gmeta.metacell_id.is_unique and (rows >= 0).all())
    check('Liu_raw_MC_library_and_cells', np.array_equal(gmeta.library_key.iloc[rows], meta.library_key) and np.array_equal(gmeta.n_cells.iloc[rows], meta.n_cells))
    cols = np.array([np.flatnonzero(genes.gene_symbol.str.upper().eq(factor)).item() for factor in factors])
    mat = sparse.load_npz(NAIVE / 'raw/metacell_by_gene_counts.npz').tocsr()
    check('Liu_raw_total_UMI_conservation', np.array_equal(np.asarray(mat.sum(1)).ravel(), gmeta.total_gex))
    counts = mat[rows][:, cols].toarray()
    total = gmeta.total_gex.to_numpy()[rows]
    expression_states.append(state_means(meta, np.log2(1 + counts / total[:, None] * 1e6), factors))
    expression_pairs = pair_effects(pd.concat(expression_states, ignore_index=True))
    frozen_cpa = tab(CPA / 'paired_library_effects.tsv')
    frozen_cpa = frozen_cpa[frozen_cpa.contrast.eq('P_M') & frozen_cpa.factor.isin(factors)]
    joined = expression_pairs.merge(frozen_cpa, on=KEY + ['factor', 'library_key', 'time'], validate='one_to_one')
    check('CPA_pair_inventory_recovered', len(joined) == len(expression_pairs) == len(frozen_cpa))
    error = np.max(np.abs(joined.delta - joined.delta_log2CPM1))
    check('CPA_pair_raw_recovery', np.allclose(joined.delta, joined.delta_log2CPM1, atol=1e-10, rtol=1e-10), max_abs_error=float(error))
    factor_rows = []
    for key, q in expression_pairs.groupby(KEY + ['factor']):
        factor_rows.append(dict(zip(KEY + ['factor'], key), n_pairs=len(q), delta=q.delta.mean(),
                                time_equal_delta=q.groupby('time').delta.mean().mean()))
    factor_summary = pd.DataFrame(factor_rows)
    route_rows = []
    for key, q in factor_summary.groupby(KEY):
        c = cycle_pairs[(cycle_pairs[KEY] == pd.Series(key, index=KEY)).all(axis=1)]
        check(str(key) + '_11_factors_same_pair_count', len(q) == 11 and q.n_pairs.nunique() == 1 and q.n_pairs.iloc[0] == len(c))
        route_rows.append(dict(zip(KEY, key), n_pairs=len(c), delta_cycle=c.delta.mean(),
                              time_equal_delta_cycle=c.groupby('time').delta.mean().mean(),
                              n_CPA_factors=len(q), n_CPA_positive=int(q.delta.gt(0).sum()),
                              delta_CPA_mean=q.delta.mean(), time_equal_delta_CPA_mean=q.time_equal_delta.mean(),
                              n_CPA_time_equal_positive=int(q.time_equal_delta.gt(0).sum())))
    route_summary = pd.DataFrame(route_rows)
    wot = frozen_cycle[frozen_cycle.study_id.eq('GSE115943')].copy()
    old = tab(OLD_PROJECT / 'analysis/tandem-apa-v1/library_replicate_design.tsv')
    wot = wot.merge(old, on='library_key', suffixes=('', '_original_design'), validate='one_to_one')
    check('WOT_24_unique_wells_and_libraries', len(wot) == wot.library_key.nunique() == wot.sample_unit.nunique() == 24)
    check('WOT_frozen_and_original_well_agreement', wot.sample_unit.eq(wot.biological_replicate_id).all())
    check('WOT_13_times_11_doubletons_2_singletons', wot.time.nunique() == 13 and (wot.groupby('time').size() == 2).sum() == 11 and (wot.groupby('time').size() == 1).sum() == 2)
    wot['day'] = wot.time.str.removeprefix('D').astype(float)
    wot['C_label'] = wot.reported_batch
    save(route_summary, 'table_audit_route_summary.tsv')
    save(factor_summary, 'table_audit_factor_summary.tsv')
    save(wot[['library_key', 'time', 'day', 'sample_unit', 'replicate_status', 'C_label',
              'condition', 'reported_donor_animal_culture', 'pooling_status', 'capture_relationship',
              'evidence_level', 'confidence', 'source_locator', 'decision_reason', 'inference_scope']], 'table_audit_wot_units.tsv')
    validation = dict(status='passed', checks=checks, n_checks=len(checks), tests_skipped=0,
                      reconstructed_CPA_pairs=len(expression_pairs), reconstructed_cycle_pairs=len(cycle_pairs),
                      input_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUTS},
                      analysis_scope='Existing cycle/CPA summaries only; no shared-27 APA associations calculated')
    (HERE / 'table_audit_validation.json').write_text(json.dumps(validation, indent=2, ensure_ascii=False) + '\n')
    print(route_summary.to_string(index=False))
    print(json.dumps({k: v for k, v in validation.items() if k not in ['checks', 'input_sha256']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
