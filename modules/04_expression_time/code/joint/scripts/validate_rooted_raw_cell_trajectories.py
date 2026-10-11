#!/usr/bin/env python3
"""Independent closure checks and direct reaggregation of raw-cell couplings."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bundle', type=Path, required=True)
    args = parser.parse_args()
    output = args.bundle
    cfg = json.loads((output/'config.json').read_text())
    raw = pd.read_csv(Path(cfg['raw_cell_result'])/'raw_cell_nodes.tsv.gz', sep='\t',
                      usecols=['cell_id', 'dataset', 'local_state_id', 'day', 'UMAP1', 'UMAP2'], low_memory=False)
    assert raw.cell_id.is_unique
    raw = raw.set_index('cell_id')
    flow = pd.read_csv(output/'all_frozen_mouse_raw_cell_flows.tsv.gz', sep='\t')
    edges = pd.read_csv(output/'rendered_ot_edges.tsv', sep='\t')
    candidates = pd.read_csv(output/'all_mouse_timed_candidates.tsv', sep='\t')
    nodes = pd.read_csv(output/'mouse_timed_nodes.tsv', sep='\t').set_index('node')
    roots = pd.read_csv(output/'observed_mouse_sources.tsv', sep='\t').set_index('route')
    endpoints = pd.read_csv(output/'terminal_selection.tsv', sep='\t')
    paths = pd.read_csv(output/'terminal_paths.tsv', sep='\t')
    checked = []
    for (route, pair), expected in flow.loc[flow.route.isin(cfg['mouse_routes'])].groupby(['route', 'pair']):
        path = Path(cfg['mouse_coupling_root'])/route/f'pair_{pair}'/'cells'/'coupling_eps_0.05.npz'
        record = json.loads(path.with_name('summary.json').read_text())['results'][0]
        assert record['unit'] == 'cell' and record['converged']
        assert record['day_from'] < record['day_to']
        assert set(expected.day_from) == {record['day_from']} and set(expected.day_to) == {record['day_to']}
        with np.load(path) as data:
            si, ti = raw.loc[data['source']], raw.loc[data['target']]
            assert si.day.eq(record['day_from']).all() and ti.day.eq(record['day_to']).all()
            assert si.dataset.eq(route).all() and ti.dataset.eq(route).all()
            s = si.dataset+'|'+si.local_state_id
            t = ti.dataset+'|'+ti.local_state_id
            sc, sn = pd.factorize(s, sort=True)
            tc, tn = pd.factorize(t, sort=True)
            a = sparse.csr_matrix((np.ones(len(si)), (sc, np.arange(len(si)))), shape=(len(sn), len(si)))
            b = sparse.csr_matrix((np.ones(len(ti)), (np.arange(len(ti)), tc)), shape=(len(ti), len(tn)))
            matrix = np.asarray((a @ data['transport']) @ b)
            matrix /= matrix.sum()
        observed = pd.DataFrame(matrix, index=sn, columns=tn)
        residual = observed.copy()
        for r in expected.itertuples():
            residual.loc[r.source, r.target] -= r.mass
        error = float(np.abs(residual.to_numpy()).max())
        assert error < 1e-8, (route, pair, error)
        checked.append({'route': route, 'pair': int(pair), 'raw_source_cells': len(si),
                        'raw_target_cells': len(ti), 'max_mass_error': error})
        print(f'Raw-cell OT verified: {route} pair {pair}, error {error:.2g}', flush=True)
    by_id = flow.set_index('raw_flow_id')
    for e in candidates.itertuples():
        original = by_id.loc[e.raw_flow_ids.split(';')]
        assert set(original.route) == {e.route} and set(original.pair) == {e.pair}
        assert abs(original.mass.sum()-e.mass) < 1e-10
    cellkeys = raw.dataset.fillna('')+'|'+raw.local_state_id.fillna('')
    verified_nodes = set()
    for route, es in edges.loc[edges.species.eq('mouse')].groupby('route'):
        root = roots.loc[route, 'node']
        assert set(es.source)-set(es.target) == {root}
        assert es.target.is_unique and (es.day_from < es.day_to).all()
        reachable = {root}
        for e in es.sort_values('day_to').itertuples():
            assert e.source in reachable
            reachable.add(e.target)
            for end, day in [('source', e.day_from), ('target', e.day_to)]:
                node_id = getattr(e, end)
                node = nodes.loc[node_id]
                if node_id not in verified_nodes:
                    members = node.original_states.split(';')
                    cs = raw.loc[cellkeys.isin(members) & raw.day.eq(day)]
                    assert len(cs) == node.n_cells
                    assert np.allclose(cs[['UMAP1', 'UMAP2']].median(), [node.x, node.y], atol=1e-7)
                    verified_nodes.add(node_id)
                assert np.allclose([node.x, node.y],
                                   [getattr(e, end+'_x'), getattr(e, end+'_y')], atol=1e-7)
        tips = set(es.target)-set(es.source)
        assert tips == set(endpoints.loc[endpoints.route.eq(route) & endpoints.selected_terminal, 'node'])
        for _, path in paths.loc[paths.route.eq(route)].groupby('terminal'):
            sequence = es.set_index('edge_id').loc[path.sort_values('step').edge_id]
            assert sequence.iloc[0].source == root
            assert list(sequence.target)[:-1] == list(sequence.source)[1:]
    previous = pd.read_csv(Path(cfg['raw_cell_result'])/'state_and_trajectory_display_v22/rendered_ot_edges.tsv', sep='\t')
    human_cols = ['route', 'source', 'target', 'source_x', 'source_y', 'target_x', 'target_y', 'source_cells', 'target_cells']
    pd.testing.assert_frame_equal(edges.loc[edges.species.eq('human'), human_cols].reset_index(drop=True),
                                  previous.loc[previous.species.eq('human'), human_cols].reset_index(drop=True),
                                  check_dtype=False)
    result = {'status': 'PASS', 'raw_cell_ot_pairs_reaggregated': len(checked),
              'max_raw_flow_mass_error': max(r['max_mass_error'] for r in checked),
              'all_selected_paths_start_at_earliest_source': True,
              'all_selected_mouse_edges_forward': True, 'same_time_endpoint_medians_verified': True,
              'human_paths_unchanged': True, 'all_candidate_flows_trace_to_raw_couplings': True,
              'metacells_used': False,
              'rendered_edges_sha256': hashlib.sha256((output/'rendered_ot_edges.tsv').read_bytes()).hexdigest()}
    pd.DataFrame(checked).to_csv(output/'raw_coupling_validation.tsv', sep='\t', index=False)
    (output/'independent_validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
