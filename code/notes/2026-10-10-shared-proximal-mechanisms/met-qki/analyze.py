"""Frozen QKI/EMT target transfer and descriptive QKI expression audit."""
from pathlib import Path
import datetime
import hashlib
import json
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
from preflight import read_sheet
sys.path.insert(0, str(ROOT / 'notes/2026-10-10-rna-module-enrichment'))
from exact_null import distribution, upper_tail, bh

COARSE = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
FULL = COARSE + ['mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']

def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.12g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)

def inputs():
    roster = pd.read_csv(ROOT / 'notes/2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz', sep='\t')
    roster = roster[roster.contrast.eq('P_M')].copy()
    features = pd.read_csv(ROOT / 'notes/2026-10-10-recurrence-utr-ortholog/ortholog_assay_features.tsv.gz', sep='\t')
    features = features[features.contrast.eq('P_M') & np.isclose(features.threshold, .02)]
    for species in ['human', 'mouse']:
        part = features[features.species.eq(species)][['gene_id'] + FULL]
        part = part.rename(columns={c: f'{species}_{c}' for c in part.columns})
        roster = roster.merge(part, on=f'{species}_gene_id', how='left', validate='one_to_one')
    assert len(roster) == 751 and roster.foreground_002.sum() == 27
    assert not roster[[f'{s}_{c}' for s in ['human', 'mouse'] for c in FULL]].isna().any().any()
    events = pd.read_csv(HERE.parent / 'event_inventory.tsv.gz', sep='\t')
    for source, value, field in [('driver_class', 'strict_tandem_UTR', 'driver_strict'),
                                 ('full_primary_axis_tandem', True, 'axis_strict')]:
        labels = events.groupby(['species', 'gene_id'])[source].apply(lambda x: bool(x.eq(value).all()))
        roster[field] = [bool(labels.get(('human', r.human_gene_id), False) and labels.get(('mouse', r.mouse_gene_id), False))
                         for r in roster.itertuples()]
    return roster

def expression():
    old = ROOT / 'notes/2026-10-02-rbp-expression'
    archive = np.load(old / 'expression_inputs.npz', allow_pickle=False)
    pairs, summaries = [], []
    for study in ['GSE118258', 'GSE150311', 'GSE103221', 'GSE115943']:
        meta = pd.DataFrame(json.loads(str(archive[study + '_meta_json'])))
        factors = list(archive[study + '_factors'])
        factor = 'QKI' if study in ['GSE118258', 'GSE150311'] else 'Qki'
        values = archive[study + '_counts'][:, factors.index(factor)]
        totals = archive[study + '_total_gex']
        assert (values >= 0).all() and (totals > 0).all()
        qualified = {}
        for (library, state), rows in meta.groupby(['library_key', 'node']).groups.items():
            indices = np.array(list(rows))
            if state not in ['P', 'M'] or len(indices) < 3 or meta.loc[indices, 'n_cells'].sum() < 30:
                continue
            qualified[library, state] = dict(mean=float(np.log2(1 + values[indices] / totals[indices] * 1e6).mean()),
                time=str(meta.loc[indices, 'experimental_time'].iloc[0]), n_metacells=len(indices),
                n_cells=int(meta.loc[indices, 'n_cells'].sum()))
        for library in sorted(set(k[0] for k in qualified)):
            if (library, 'P') not in qualified or (library, 'M') not in qualified:
                continue
            p, m = qualified[library, 'P'], qualified[library, 'M']
            assert p['time'] == m['time']
            pairs.append(dict(study_id=study, route=meta.route.iloc[0], factor=factor, library_key=library,
                time=p['time'], P_log2CPM1=p['mean'], M_log2CPM1=m['mean'],
                delta_log2CPM1=p['mean']-m['mean'], n_metacells_P=p['n_metacells'],
                n_metacells_M=m['n_metacells'], n_cells_P=p['n_cells'], n_cells_M=m['n_cells']))
        group = pd.DataFrame([r for r in pairs if r['study_id'] == study])
        summaries.append(dict(study_id=study, route=meta.route.iloc[0], factor=factor,
            n_pairs=len(group), delta_log2CPM1=group.delta_log2CPM1.mean(),
            time_equal_delta_log2CPM1=group.groupby('time').delta_log2CPM1.mean().mean(),
            sign='decrease' if group.delta_log2CPM1.mean() < 0 else 'increase'))
    save(pd.DataFrame(pairs), 'qki_expression_pairs.tsv')
    save(pd.DataFrame(summaries), 'qki_expression_summary.tsv')
    return summaries

def main():
    freeze = json.loads((HERE / 'protocol_freeze.json').read_text())
    for name, digest in freeze['files'].items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
    roster = inputs()
    results, descriptions, global_rows, members = [], [], [], []
    panels = {'all': np.ones(len(roster), dtype=bool), 'driver_strict': roster.driver_strict,
              'axis_strict': roster.axis_strict}
    for filename, label in [('KRNB_A_2294222_SM0790.xls', 'QKI_KD'), ('KRNB_A_2294222_SM0800.xls', 'EMT')]:
        frame = read_sheet(filename)
        for c in ['effect', 'fdr_zero']:
            frame[c] = pd.to_numeric(frame[c], errors='coerce')
        duplicated = frame.name.duplicated(keep=False)
        finite = np.isfinite(frame.effect) & np.isfinite(frame.fdr_zero)
        valid = frame[~duplicated & finite]
        assert valid.name.is_unique
        save(frame, f'external_{label}.tsv.gz')
        for direction, predicate in [('proximal', valid.effect.lt(0)), ('distal', valid.effect.gt(0))]:
            global_rows.append(dict(comparison=label, direction=direction, evaluated_genes=len(valid),
                                    n_response=int((predicate & valid.fdr_zero.lt(.05)).sum())))
        selected = valid[['name', 'gene', 'effect', 'fdr_zero', 'AveExpr', 'total_reads', 'Length']].rename(
            columns={'name': 'human_gene_id', 'gene': 'external_symbol'})
        full = roster.merge(selected, on='human_gene_id', how='left', validate='one_to_one')
        full['external_evaluable'] = full.effect.notna() & full.fdr_zero.notna()
        duplicate_ids = set(frame.loc[duplicated, 'name'])
        external_ids = set(frame.name)
        full['missing_reason'] = ['' if r.external_evaluable else 'duplicate_external_id' if r.human_gene_id in duplicate_ids else
                                  'nonfinite_external_APA_test' if r.human_gene_id in external_ids else
                                  'absent_external_APA_test' for r in full.itertuples()]
        full['comparison'] = label
        members.append(full)
        for panel, flags in panels.items():
            eligible = full[flags & full.external_evaluable].copy()
            fg = eligible.foreground_002.to_numpy(bool)
            for strict in [False, True]:
                selected_fg = eligible.foreground_005.to_numpy(bool) if strict else fg
                for fdr in [.01, .05, .10]:
                    for name, part in [('background', eligible), ('core', eligible[selected_fg])]:
                        descriptions.append(dict(comparison=label, panel=panel, discovery_threshold=.05 if strict else .02,
                            group=name, external_fdr=fdr, n_evaluable=len(part),
                            median_effect=part.effect.median(), n_negative=int(part.effect.lt(0).sum()),
                            n_positive=int(part.effect.gt(0).sum()), n_zero=int(part.effect.eq(0).sum()),
                            n_proximal_FDR=int((part.effect.lt(0) & part.fdr_zero.lt(fdr)).sum()),
                            n_distal_FDR=int((part.effect.gt(0) & part.fdr_zero.lt(fdr)).sum())))
            for direction in ['proximal', 'distal']:
                response = eligible.fdr_zero.lt(.05) & (eligible.effect.lt(0) if direction == 'proximal' else eligible.effect.gt(0))
                eligible['response'] = response
                observed = int((fg & response.to_numpy(bool)).sum())
                for model in ['unmatched', 'COARSE', 'FULL']:
                    params, variable_fg = [], 0
                    if model == 'unmatched':
                        groups = [(None, eligible)]
                    else:
                        columns = [f'{s}_{c}' for s in ['human', 'mouse'] for c in (COARSE if model == 'COARSE' else FULL)]
                        groups = list(eligible.groupby(columns, dropna=False, sort=False))
                    for _, part in groups:
                        N, K, n = len(part), int(part.response.sum()), int(part.foreground_002.sum())
                        params.append((N, K, n))
                        if min(n, K) > max(0, n + K - N):
                            variable_fg += n
                    null = distribution(tuple(params))
                    x = np.arange(null[0], null[0] + len(null[1]))
                    expected = float(x @ null[1])
                    variance = float(((x-expected)**2) @ null[1])
                    p = upper_tail(null, observed) if variance > 1e-14 else np.nan
                    results.append(dict(comparison=label, panel=panel, direction=direction, model=model,
                        N=len(eligible), n_core=int(fg.sum()), n_external_response=int(response.sum()),
                        observed=observed, expected=expected, fold=observed/expected if expected else np.nan,
                        null_variance=variance, n_strata=len(params), variable_foreground=variable_fg,
                        p_value=p, status='estimable' if np.isfinite(p) else 'zero_null_variance'))
    global_table = pd.DataFrame(global_rows)
    expected_counts = {('QKI_KD', 'proximal'): 94, ('QKI_KD', 'distal'): 25, ('EMT', 'proximal'): 27, ('EMT', 'distal'): 52}
    for row in global_table.itertuples():
        assert row.n_response == expected_counts[(row.comparison, row.direction)], (row.comparison, row.direction, row.n_response)
    result = pd.DataFrame(results)
    assert len(result) == 36
    result['q_BH36'] = bh(result.p_value.fillna(1).to_numpy())
    save(result, 'overlap_tests.tsv')
    save(global_table, 'author_response_count_check.tsv')
    save(pd.DataFrame(descriptions), 'signed_response_descriptive.tsv')
    save(pd.concat(members, ignore_index=True), 'ortholog_external_roster.tsv.gz')
    save(pd.concat(members, ignore_index=True).query('foreground_002'), 'all_core_external_members.tsv')
    expr = expression()
    audit = dict(status='production_completed', timestamp_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 protocol_sha256=freeze['files']['PROTOCOL.zh.md'], tests_family=36,
                 n_estimable=int(result.p_value.notna().sum()), n_zero_variance=int(result.p_value.isna().sum()),
                 n_q_below_005=int(result.q_BH36.lt(.05).sum()), author_counts_reproduced=True,
                 expression_summary=expr)
    (HERE / 'analysis_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(result.query("panel == 'all'").to_string(index=False))
    print(json.dumps(audit, indent=2))

if __name__ == '__main__':
    main()
