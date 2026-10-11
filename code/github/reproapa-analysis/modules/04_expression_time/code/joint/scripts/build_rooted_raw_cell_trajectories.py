#!/usr/bin/env python3
"""Extract source-connected display paths from complete frozen raw-cell OT flows.

These are representative state-time paths, not fitted lineage trees. Only the
earliest same-route S cells are pooled for the source glyph; every later node
retains its original local state and exact experimental day.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np
import pandas as pd


def timed_node(state, day):
    return f'{state}|day={day:g}'


def best_rooted_paths(edges, root, terminals, min_mass, min_target_fraction):
    """Maximise summed log backward-OT support on a strictly forward-time DAG.

    Each timed node keeps one globally best predecessor. The union of terminal
    paths therefore shares its trunk, cannot acquire intermediate roots, and
    cannot merge again after splitting. This is a display selection rule only.
    """
    if not (edges.day_from < edges.day_to).all():
        raise ValueError('All candidate OT edges must be strictly forward in time')
    if edges.duplicated(['source', 'target']).any():
        raise ValueError('Duplicate timed edge')
    scores, parents = {root: 0.}, {}
    eligible = edges.loc[edges.mass.ge(min_mass) & edges.target_fraction.ge(min_target_fraction)]
    for e in eligible.sort_values(['day_to', 'target', 'source']).itertuples():
        if e.source not in scores:
            continue
        if not 0 < e.target_fraction <= 1 + 1e-8:
            raise ValueError('Invalid backward OT fraction')
        score = scores[e.source] + math.log(min(e.target_fraction, 1.))
        if score > scores.get(e.target, -np.inf):
            scores[e.target], parents[e.target] = score, e.Index
    selected, paths = set(), []
    for terminal in sorted(terminals):
        if terminal not in scores:
            raise ValueError(f'Terminal cannot be reached from the observed source: {terminal}')
        current, path = terminal, []
        while current != root:
            index = parents[current]
            selected.add(index)
            path.append(index)
            current = edges.loc[index, 'source']
        paths.append({'terminal': terminal, 'edge_indices': path[::-1],
                      'log_backward_support': scores[terminal]})
    chosen = edges.loc[sorted(selected)].copy()
    if set(chosen.source) - set(chosen.target) != {root}:
        raise ValueError('An internal segment was incorrectly treated as a root')
    return chosen, paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, type=Path)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    base = Path(cfg['raw_cell_result'])
    output = base / cfg['output_subdirectory']
    output.mkdir(exist_ok=False)
    raw = pd.read_csv(base/'raw_cell_nodes.tsv.gz', sep='\t', low_memory=False,
                      usecols=['cell_id', 'species_key', 'dataset', 'day', 'local_state_id', 'UMAP1', 'UMAP2'])
    labels = pd.read_csv(base/'state_and_trajectory_display_v13/raw_cell_stateflow_groups.tsv.gz',
                         sep='\t', usecols=['cell_id', 'state_group'])
    cells = raw.merge(labels, on='cell_id', validate='one_to_one')
    assert len(cells) == 551201 and cells.cell_id.is_unique
    flows = pd.read_csv(cfg['mouse_all_flows'], sep='\t')
    flows['raw_flow_id'] = [f'MF{i+1:04}' for i in range(len(flows))]
    assert (flows.day_from < flows.day_to).all()
    assert np.allclose(flows.groupby(['route', 'pair']).mass.sum(), 1., atol=1e-8)
    summaries, node_tables, candidate_tables, selected_tables, path_rows, root_rows, endpoint_rows = [], [], [], [], [], [], []
    for route in cfg['mouse_routes']:
        current = cells.loc[cells.dataset.eq(route)].copy()
        route_flows = flows.loc[flows.route.eq(route)].copy()
        first, last = route_flows.day_from.min(), route_flows.day_to.max()
        pairs = route_flows[['pair', 'day_from', 'day_to']].drop_duplicates().sort_values('day_from')
        if not np.array_equal(pairs.day_to.to_numpy()[:-1], pairs.day_from.to_numpy()[1:]):
            raise ValueError(f'Missing consecutive frozen time pair in {route}')
        assert first == current.day.min() and last == current.day.max()
        source_mask = current.day.eq(first) & current.state_group.eq(cfg['source_state_group'])
        source_cells = current.loc[source_mask]
        if len(source_cells) < cfg['minimum_source_cells']:
            raise ValueError(f'No adequate annotated early source in {route}')
        root_states = set(route+'|'+source_cells.local_state_id)
        root = timed_node(route+'|SOURCE_S', first)
        current['state'] = route+'|'+current.local_state_id
        current['node'] = [timed_node(s, d) for s, d in zip(current.state, current.day)]
        current.loc[source_mask, 'node'] = root
        nodes = current.groupby('node', sort=True).agg(
            day=('day', 'first'), state_group=('state_group', 'first'), n_cells=('cell_id', 'size'),
            x=('UMAP1', 'median'), y=('UMAP2', 'median'),
            original_states=('state', lambda s: ';'.join(sorted(set(s)))))
        nodes['route'] = route
        node_tables.append(nodes.reset_index())
        root_rows.append({'route': route, 'node': root, 'day': first,
                          'original_states': ';'.join(sorted(root_states)), 'n_source_cells': len(source_cells),
                          'n_earliest_cells': int(current.day.eq(first).sum()),
                          'source_cell_fraction': len(source_cells)/current.day.eq(first).sum(),
                          'x': nodes.loc[root, 'x'], 'y': nodes.loc[root, 'y']})
        route_flows['original_source'] = route_flows.source
        route_flows['original_target'] = route_flows.target
        route_flows['source'] = [root if d == first and s in root_states else timed_node(s, d)
                                for s, d in zip(route_flows.source, route_flows.day_from)]
        route_flows['target'] = [timed_node(s, d) for s, d in zip(route_flows.target, route_flows.day_to)]
        # Summation is only over earliest S substates with the same target/time.
        edges = route_flows.groupby(['route', 'pair', 'day_from', 'day_to', 'source', 'target'],
                                   as_index=False, sort=True).agg(
            mass=('mass', 'sum'), raw_flow_ids=('raw_flow_id', lambda s: ';'.join(s)))
        edges['source_fraction'] = edges.mass / edges.groupby(['pair', 'source']).mass.transform('sum')
        edges['target_fraction'] = edges.mass / edges.groupby(['pair', 'target']).mass.transform('sum')
        edges['source_cells'] = edges.source.map(nodes.n_cells)
        edges['target_cells'] = edges.target.map(nodes.n_cells)
        if edges[['source_cells', 'target_cells']].isna().any().any():
            raise ValueError('A frozen OT state-time endpoint has no original cells')
        edges['edge_id'] = [f'{route}::E{i+1:04}' for i in range(len(edges))]
        edges['eligible'] = (edges.mass.ge(cfg['minimum_joint_mass']) &
                             edges.target_fraction.ge(cfg['minimum_target_fraction']) &
                             edges.source_cells.gt(0) & edges.target_cells.gt(0))
        last_nodes = nodes.loc[nodes.day.eq(last)].copy()
        last_nodes['terminal_cell_fraction'] = last_nodes.n_cells / last_nodes.n_cells.sum()
        last_nodes['eligible_terminal'] = (last_nodes.n_cells.ge(cfg['minimum_terminal_cells']) &
                                           last_nodes.terminal_cell_fraction.ge(cfg['minimum_terminal_cell_fraction']))
        reachable = {root}
        for edge in edges.loc[edges.eligible].sort_values('day_to').itertuples():
            if edge.source in reachable:
                reachable.add(edge.target)
        last_nodes['source_reachable'] = last_nodes.index.isin(reachable)
        last_nodes['selected_terminal'] = last_nodes.eligible_terminal & last_nodes.source_reachable
        last_nodes['exclusion_reason'] = np.where(last_nodes.selected_terminal, '',
            np.where(~last_nodes.eligible_terminal, 'below final-time cell count/fraction gate',
                     'no complete source path under mass/support gates'))
        terminals = last_nodes.index[last_nodes.selected_terminal]
        if len(terminals) == 0:
            raise ValueError(f'No adequately supported source-to-terminal path in {route}')
        endpoint_rows.append(last_nodes.reset_index())
        chosen, paths = best_rooted_paths(edges.loc[edges.eligible], root, terminals,
                                          cfg['minimum_joint_mass'], cfg['minimum_target_fraction'])
        edges['selected'] = edges.index.isin(chosen.index)
        candidate_tables.append(edges)
        for endpoint in ['source', 'target']:
            for axis in ['x', 'y']:
                chosen[f'{endpoint}_{axis}'] = chosen[endpoint].map(nodes[axis])
                chosen[f'display_{endpoint}_{axis}'] = chosen[f'{endpoint}_{axis}']
        chosen['species'] = 'mouse'
        chosen['displayed'] = True
        chosen['frozen_edge_id'] = chosen.edge_id
        chosen['evidence'] = 'complete frozen raw-cell OT, source-connected representative time path'
        selected_tables.append(chosen)
        for path in paths:
            for step, index in enumerate(path['edge_indices']):
                path_rows.append({'route': route, 'terminal': path['terminal'], 'step': step,
                                  'edge_id': edges.loc[index, 'edge_id'],
                                  'log_backward_support': path['log_backward_support']})
        summaries.append({'route': route, 'observed_start': first, 'observed_end': last,
                          'time_pairs': len(pairs), 'raw_flow_rows': len(route_flows),
                          'eligible_edges': int(edges.eligible.sum()), 'selected_edges': len(chosen),
                          'roots': 1, 'terminals': len(terminals),
                          'terminal_candidates_excluded_by_connectivity': int((last_nodes.eligible_terminal & ~last_nodes.source_reachable).sum()),
                          'minimum_selected_mass': float(chosen.mass.min()),
                          'minimum_selected_target_fraction': float(chosen.target_fraction.min()),
                          'terminal_cell_coverage': float(last_nodes.loc[terminals, 'terminal_cell_fraction'].sum())})
        print(summaries[-1], flush=True)
    previous = pd.read_csv(base/'state_and_trajectory_display_v22/rendered_ot_edges.tsv', sep='\t')
    human = previous.loc[previous.species.eq('human')].copy()
    selected = pd.concat([human, *selected_tables], ignore_index=True)
    selected.to_csv(output/'rendered_ot_edges.tsv', sep='\t', index=False)
    pd.concat(candidate_tables).to_csv(output/'all_mouse_timed_candidates.tsv', sep='\t', index=False)
    pd.concat(node_tables).to_csv(output/'mouse_timed_nodes.tsv', sep='\t', index=False)
    pd.concat(endpoint_rows).to_csv(output/'terminal_selection.tsv', sep='\t', index=False)
    pd.DataFrame(root_rows).to_csv(output/'observed_mouse_sources.tsv', sep='\t', index=False)
    pd.DataFrame(path_rows).to_csv(output/'terminal_paths.tsv', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(output/'route_summary.tsv', sep='\t', index=False)
    flows.to_csv(output/'all_frozen_mouse_raw_cell_flows.tsv.gz', sep='\t', index=False)
    shutil.copyfile(args.config, output/'config.json')
    summary = {'raw_cells': len(cells), 'raw_cell_coordinates_changed': False, 'metacells_used': False,
               'human_visible_edges_unchanged': len(human), 'displayed_edges': len(selected),
               'displayed_routes': selected.route.nunique(), 'mouse_routes': summaries,
               'mouse_source_count': len(root_rows), 'mouse_candidate_flow_rows': len(flows),
               'source_pool_rule': 'earliest observed same-route cells already assigned S only',
               'node_identity': 'route + original local state + exact day; earliest S pooled',
               'path_rule': 'maximum sum log target_fraction from source to each eligible final-time state; union of paths',
               'terminal_rule': 'at least 50 raw cells and at least 1% of final observed time cells',
               'not_lineage_probability': True,
               'mouse_flow_source': cfg['mouse_all_flows'],
               'mouse_flow_sha256': hashlib.sha256(Path(cfg['mouse_all_flows']).read_bytes()).hexdigest()}
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
