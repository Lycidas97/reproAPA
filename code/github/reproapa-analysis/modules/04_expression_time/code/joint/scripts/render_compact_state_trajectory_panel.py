#!/usr/bin/env python3
"""Re-layout the accepted v13 raw-cell state/OT display at publication size."""

import argparse
import json
import shutil
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, PathPatch
from matplotlib.path import Path as PlotPath
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from render_state_and_trajectory_display import (
    configure_publication_font, maximal_trajectory_chains, save_all, setup_axis,
)
from smooth_trajectory_display import smooth_chain


def prepare_display_edges(arrows, mouse_times, ax, dpi):
    """Hide tiny isolated components and preserve the time identity of every endpoint."""
    display = arrows.copy()
    display['frozen_edge_id'] = [f'OT{i+1:03}' for i in range(len(display))]
    start_xy = display[['source_x', 'source_y']].to_numpy()
    end_xy = display[['target_x', 'target_y']].to_numpy()
    display['original_length_mm'] = np.linalg.norm(
        ax.transData.transform(end_xy)-ax.transData.transform(start_xy), axis=1
    ) * 25.4 / dpi
    # Apply the visual cleanup to the existing topology, before splitting timed nodes.
    # Short pieces of a longer trajectory must not be removed by this rule.
    display['displayed'] = True
    for _, route_edges in display.groupby('route'):
        sources = route_edges.source.value_counts()
        targets = route_edges.target.value_counts()
        for edge in route_edges.itertuples():
            isolated = (edge.source not in targets and edge.target not in sources
                        and sources[edge.source] == 1 and targets[edge.target] == 1)
            if isolated and edge.original_length_mm < 2.0:
                display.loc[edge.Index, 'displayed'] = False

    keys = ['route', 'source', 'target']
    mouse = display.species.eq('mouse')
    timing = display.loc[mouse, keys].merge(
        mouse_times[keys+['day_from', 'day_to', 'source_fraction', 'target_fraction']],
        on=keys, how='left', sort=False, validate='one_to_one',
    )
    if timing.day_from.isna().any() or not (timing.day_from < timing.day_to).all():
        raise ValueError('Every mouse edge must resolve to its frozen forward time pair')
    display['source_state'] = display.source
    display['target_state'] = display.target
    for column in ['day_from', 'day_to', 'source_fraction', 'target_fraction']:
        display.loc[mouse, column] = timing[column].to_numpy()
    display.loc[mouse, 'source'] = [f'{state}|day={day:g}' for state, day in
                                  zip(timing.source, timing.day_from)]
    display.loc[mouse, 'target'] = [f'{state}|day={day:g}' for state, day in
                                  zip(timing.target, timing.day_to)]
    for endpoint in ['source', 'target']:
        for axis in ['x', 'y']:
            display[f'display_{endpoint}_{axis}'] = display[f'{endpoint}_{axis}']

    endpoints = pd.concat([
        display.loc[mouse, ['route', 'source_state', 'day_from', 'source_x', 'source_y']].rename(
            columns={'source_state': 'state', 'day_from': 'day', 'source_x': 'x', 'source_y': 'y'}),
        display.loc[mouse, ['route', 'target_state', 'day_to', 'target_x', 'target_y']].rename(
            columns={'target_state': 'state', 'day_to': 'day', 'target_x': 'x', 'target_y': 'y'}),
    ], ignore_index=True).drop_duplicates()
    if endpoints.duplicated(['route', 'state', 'day']).any():
        raise ValueError('One timed OT node has inconsistent raw-cell endpoint coordinates')
    counts = endpoints.groupby(['route', 'state']).day.nunique().rename('distinct_times')
    audit = endpoints.merge(counts, on=['route', 'state'])
    return display, audit.loc[audit.distinct_times > 1].sort_values(['route', 'state', 'day'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', default='state_and_trajectory_display_v22')
    parser.add_argument('--rooted-trajectories', type=Path)
    parser.add_argument('--smooth-trajectories', action='store_true')
    args = parser.parse_args()
    if args.smooth_trajectories and not args.rooted_trajectories:
        parser.error('Smoothing requires the accepted rooted trajectory input')
    cfg = json.loads(args.config.read_text())
    font = configure_publication_font()
    root = Path(cfg['raw_cell_result'])
    source = root / 'state_and_trajectory_display_v13'
    output = root / args.output
    output.mkdir(exist_ok=False)
    groups = cfg['state_groups']
    display_labels = {g['id']: g['display_label'] for g in groups}
    display_labels.update({'shared_s': 'S', 'shared_r': 'R', 'shared_e': 'E',
                           'shared_tp': 'P', 'shared_m': 'M', 'shared_k': 'K'})
    nodes = pd.read_csv(source / 'raw_cell_stateflow_groups.tsv.gz', sep='\t', low_memory=False)
    arrows = pd.read_csv(source / 'direct_raw_cell_ot_arrows_display_coordinates.tsv', sep='\t')
    assert len(nodes) == 551201 and nodes.cell_id.is_unique
    assert len(arrows) == 61 and arrows.route.nunique() == 12
    assert set(nodes.state_group) == {g['id'] for g in groups}
    xy = nodes[['UMAP1', 'UMAP2']].to_numpy()
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    pad = .035 * (hi - lo)
    limits = (lo[0]-pad[0], hi[0]+pad[0], lo[1]-pad[1], hi[1]+pad[1])
    fig = plt.figure(figsize=(170/25.4, 112/25.4), dpi=300)
    ax = fig.add_axes([.0, .025, .683, .95])
    setup_axis(ax, limits)
    fig.canvas.draw()
    mouse_times = pd.read_csv(cfg['mouse_temporal_display_edges'], sep='\t')
    display_edges, time_audit = prepare_display_edges(arrows, mouse_times, ax, fig.dpi)
    visible_edges = display_edges.loc[display_edges.displayed]
    suppressed = display_edges.loc[~display_edges.displayed].copy()
    suppressed['display_exclusion_reason'] = 'isolated one-edge component shorter than 2 mm'
    rooted_summary = None
    if args.rooted_trajectories:
        rooted_summary = json.loads((args.rooted_trajectories/'summary.json').read_text())
        visible_edges = pd.read_csv(args.rooted_trajectories/'rendered_ot_edges.tsv', sep='\t')
        if not visible_edges.displayed.all():
            raise ValueError('Rooted input contains a hidden edge')
    colors = nodes.state_group.map({g['id']: g['color'] for g in groups}).to_numpy()
    order = np.random.default_rng(cfg['seed']).permutation(len(nodes))
    # Interleave all groups and species to prevent the last category covering the others.
    for chunk in np.array_split(order, 12):
        for species, marker in [('human', 'o'), ('mouse', 's')]:
            idx = chunk[nodes.species_key.to_numpy()[chunk] == species]
            ax.scatter(*xy[idx].T, c=colors[idx], marker=marker, s=.43, alpha=.46,
                       linewidths=0, rasterized=True, zorder=1)

    routes = cfg['trajectory_routes']
    trajectory_colors = {item['route']: item['color'] for item in routes}
    if (set(trajectory_colors) != set(visible_edges.route)
            or len(trajectory_colors) != len(routes)
            or len(set(trajectory_colors.values())) != len(routes)):
        raise ValueError('Visible routes must have exactly one distinct colour and legend entry each')
    for item in routes:
        if set(visible_edges.loc[visible_edges.route.eq(item['route']), 'species']) != {item['species']}:
            raise ValueError(f"Trajectory legend species mismatch for {item['route']}")
    chain_count = terminal_count = 0
    trajectory_markers = []
    smoothing_audit, curve_vertices = [], []
    for route, edges in visible_edges.groupby('route', sort=True):
        chains, starts, terminals = maximal_trajectory_chains(edges)
        if rooted_summary and len(starts) != 1:
            raise ValueError(f'Unexpected internal display root in {route}')
        coords = {}
        color = trajectory_colors[route]
        for edge in edges.itertuples():
            coords[edge.source] = (edge.display_source_x, edge.display_source_y)
            coords[edge.target] = (edge.display_target_x, edge.display_target_y)
        for chain in chains:
            chain_count += 1
            route_xy = np.asarray([coords[node] for node in chain])
            if args.smooth_trajectories:
                path_mm, audit = smooth_chain(ax.transData.transform(route_xy)*25.4/fig.dpi)
                route_xy = ax.transData.inverted().transform(path_mm.vertices*fig.dpi/25.4)
                # Preserve pinned coordinates exactly, including floating-point roundoff.
                route_xy[0], route_xy[-1] = coords[chain[0]], coords[chain[-1]]
                smooth_path = PlotPath(route_xy, path_mm.codes)
                ax.add_patch(PathPatch(smooth_path, facecolor='none', edgecolor=color,
                    lw=1.05, zorder=3, capstyle='round', joinstyle='round',
                    path_effects=[pe.Stroke(linewidth=1.45, foreground='white'), pe.Normal()]))
                smoothing_audit.append({'route': route, 'chain': chain_count,
                    'start_node': chain[0], 'end_node': chain[-1], **audit})
                curve_vertices.extend({'route': route, 'chain': chain_count, 'vertex': i,
                    'path_code': int(code), 'x': float(point[0]), 'y': float(point[1])}
                    for i, (point, code) in enumerate(zip(route_xy, path_mm.codes)))
            else:
                ax.plot(*route_xy.T, color=color, lw=1.05, alpha=1, zorder=3,
                        solid_capstyle='round', solid_joinstyle='round',
                        path_effects=[pe.Stroke(linewidth=1.45, foreground='white'), pe.Normal()])
            if chain[-1] in terminals:
                terminal_count += 1
                # Outline only the triangular head, not the final trajectory segment.
                tail, tip = ax.transData.transform(route_xy[-2:])
                direction = (tip-tail) / np.linalg.norm(tip-tail)
                normal = np.array([-direction[1], direction[0]])
                unit = 7.5 * fig.dpi / 72
                base = tip - .4 * unit * direction
                vertices = np.array([tip, base + .2*unit*normal, base - .2*unit*normal])
                ax.add_patch(Polygon(ax.transData.inverted().transform(vertices), closed=True,
                    facecolor=color, edgecolor='white', linewidth=.4, joinstyle='round', zorder=4))
                trajectory_markers.append((route, 'terminal', Bbox.from_extents(
                    vertices[:, 0].min()-2, vertices[:, 1].min()-2,
                    vertices[:, 0].max()+2, vertices[:, 1].max()+2)))
        for start in starts:
            ax.scatter(*coords[start], s=7.5, facecolor='white', edgecolor=color, lw=.65, zorder=5)
            px, py = ax.transData.transform(coords[start])
            radius = (np.sqrt(7.5)+.65)*fig.dpi/144
            trajectory_markers.append((route, 'start', Bbox.from_extents(
                px-radius, py-radius, px+radius, py+radius)))

    # Density weighted by local state purity places labels within supported state regions.
    bins = 110
    ranges = [[limits[0], limits[1]], [limits[2], limits[3]]]
    total, xedges, yedges = np.histogram2d(*xy.T, bins=bins, range=ranges)
    total = gaussian_filter(total, 1.2)
    annotations = []
    anchors = []
    offsets = {'human_u': (0, 8), 'shared_k': (10, -6)}
    if rooted_summary:
        offsets.update({'shared_s': (0, 8), 'shared_tp': (-9, 0), 'shared_m': (7, 3), 'human_u': (6, 8)})
    if args.smooth_trajectories:
        offsets['shared_m'] = (10, 3)
    for group in groups:
        selected = nodes.state_group.eq(group['id']).to_numpy()
        density = gaussian_filter(np.histogram2d(*xy[selected].T, bins=bins, range=ranges)[0], 1.2)
        score = density ** 2 / np.maximum(total, 1e-6)
        i, j = np.unravel_index(np.argmax(score), score.shape)
        anchor = ((xedges[i]+xedges[i+1])/2, (yedges[j]+yedges[j+1])/2)
        label = display_labels[group['id']]
        offset = offsets.get(group['id'], (0, 0))
        text = ax.annotate(label, xy=anchor, xytext=offset, textcoords='offset points',
            ha='center', va='center', fontsize=7, weight='bold', color='#17212B', zorder=6,
            path_effects=[pe.withStroke(linewidth=2.5, foreground='white', alpha=.98)])
        annotations.append(text)
        anchors.append({'state_group': group['id'], 'label': label,
                        'anchor_x': anchor[0], 'anchor_y': anchor[1],
                        'offset_x_pt': offset[0], 'offset_y_pt': offset[1]})

    key = fig.add_axes([.714, .04, .25, .92])
    key.set(xlim=(0, 42.5), ylim=(103.04, 0))
    key.set_axis_off()
    key.text(0, 0, 'State groups', fontsize=8, va='top', color='#17212B')
    legend_sections = [
        ('Shared / broad-process groups', [
            ('shared_s', 'Source-like / fibroblast'),
            ('shared_r', 'Source exit / remodeling'),
            ('shared_e', 'Epithelial / mixed transition'),
            ('shared_tp', 'Pluripotency acquisition'),
            ('shared_m', 'Mesenchymal / ECM retention'),
            ('shared_k', 'Keratin-related candidate')]),
        ('Mouse groups', [
            ('mouse_b', 'Trophoblast-like'),
            ('mouse_n', 'Neural / neurokeratin'),
            ('mouse_i', 'IFN-high transition'),
            ('mouse_u', 'Unresolved / cycle-dominant'),
            ('mouse_y', 'Myeloid-local')]),
        ('Human groups', [
            ('human_x', 'Endoderm-secretory'),
            ('human_n', 'Naive-related candidate'),
            ('human_u', 'Mixed / unresolved')]),
    ]
    legend_labels = {group_id: f'{display_labels[group_id]}  {description}'
                     for _, items in legend_sections for group_id, description in items}
    palette = {g['id']: g['color'] for g in groups}
    y = 6
    for title, items in legend_sections:
        key.text(0, y, title, fontsize=6.5, color='#68717A', va='center')
        y += 4.0
        for group_id, _ in items:
            key.plot(1.1, y, marker='o', ms=3.8, color=palette[group_id], markeredgewidth=0)
            key.text(3.4, y, legend_labels[group_id], fontsize=6.5, va='center', color='#27313A')
            y += 3.1
        y += 1.7
    key.plot([0, 41], [y-1, y-1], color='#DDE1E4', lw=.5)
    key.text(0, y+3, 'Trajectories · dataset / route', fontsize=6.5, va='center', color='#68717A')
    for i, item in enumerate(routes):
        row_y = y + 8 + i*3.3
        key.annotate('', (4.2, row_y), (0, row_y),
                     arrowprops={'arrowstyle': '-|>', 'lw': 1.05, 'color': item['color'],
                                 'mutation_scale': 6, 'shrinkA': 0, 'shrinkB': 0})
        key.text(5.4, row_y, item['label'], fontsize=6.5, va='center', color='#27313A')
    key.plot(1.2, y+32, marker='o', ms=3.5, mfc='white', mec='#485560', mew=.6)
    key.text(3.4, y+32, 'Observed start' if rooted_summary else 'Segment start',
             fontsize=6.5, va='center', color='#27313A')
    key.text(0, y+36.2, 'Cells', fontsize=6.5, va='center', color='#68717A')
    for x, marker, label in [(10, 'o', 'Human'), (27, 's', 'Mouse')]:
        key.plot(x, y+36.2, marker=marker, ms=2.6, color='#68717A', markeredgewidth=0)
        key.text(x+2.4, y+36.2, label, fontsize=6.5, va='center', color='#27313A')

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [text.get_window_extent(renderer).expanded(1.06, 1.15) for text in annotations]
    collisions = [(anchors[i]['state_group'], anchors[j]['state_group'])
                  for i in range(len(boxes)) for j in range(i+1, len(boxes)) if boxes[i].overlaps(boxes[j])]
    marker_collisions = [(anchors[i]['state_group'], route, kind)
                         for i, box in enumerate(boxes) for route, kind, marker in trajectory_markers
                         if box.padded(1.25*fig.dpi/72).overlaps(marker)]
    if rooted_summary and marker_collisions:
        raise ValueError(f'State label obscures a trajectory start/head: {marker_collisions}')
    key_bounds = [text.get_window_extent(renderer) for text in key.texts]
    legend_collisions = [(i, j) for i in range(len(key_bounds)) for j in range(i+1, len(key_bounds))
                         if key_bounds[i].overlaps(key_bounds[j])]
    if legend_collisions:
        raise ValueError(f'Legend text overlaps: {legend_collisions}')
    key_right = max(box.x1 for box in key_bounds)
    if key_right > fig.bbox.x1 - 6:
        raise ValueError('Legend text exceeds the publication margin')
    review = save_all(fig, output, 'human_mouse_raw_cell_states_with_direct_ot_trajectories')
    review.update({'status': 'PASS' if not collisions else 'LABEL_REVIEW',
        'source_display': str(source), 'raw_cells': len(nodes), 'state_groups': 14,
        'raw_cell_coordinates_and_assignments_changed': False, 'metacells_used': False,
        'direct_ot_edges': len(arrows), 'routes': arrows.route.nunique(),
        'displayed_ot_edges': len(visible_edges), 'displayed_routes': visible_edges.route.nunique(),
        'trajectory_colour_unit': 'dataset / frozen study route',
        'trajectory_colours': trajectory_colors, 'trajectory_legend_entries': len(routes),
        'trajectory_geo_datasets': len({item['dataset'] for item in routes}),
        'legend_text_collisions': legend_collisions,
        'suppressed_isolated_short_edges': len(suppressed), 'short_edge_threshold_mm': 2.0,
        'mouse_node_identity': 'route + original state + experimental day',
        'open_circle_meaning': 'start of a displayed OT segment, not necessarily start of reprogramming',
        'ot_endpoint_coordinates': 'original per-time raw-cell medians, without cross-time averaging',
        'previously_collapsed_multitime_mouse_states': int(time_audit[['route', 'state']].drop_duplicates().shape[0]),
        'maximal_chains': chain_count, 'terminal_arrowheads': terminal_count,
        'panel_mm': [170, 112], 'font': 'Nimbus Sans', 'font_file': str(font),
        'minimum_text_pt': 6.5, 'cell_alpha': .46, 'state_label_collisions': collisions,
        'trajectory_white_outline': True, 'trajectory_linewidth_pt': 1.05,
        'trajectory_understroke_linewidth_pt': 1.45, 'trajectory_visible_outline_per_side_pt': .2,
        'arrowhead_white_outline': True, 'arrowhead_outline_linewidth_pt': .4,
        'display_labels': display_labels,
        'layout': 'left UMAP, compact right legend, unboxed state labels',
        'label_anchor_rule': 'smoothed state density multiplied by local state purity'})
    if rooted_summary:
        for field in ['direct_ot_edges', 'routes', 'suppressed_isolated_short_edges',
                      'short_edge_threshold_mm', 'previously_collapsed_multitime_mouse_states']:
            review['legacy_'+field] = review.pop(field)
        review.update({'rooted_trajectory_summary': rooted_summary,
            'direct_ot_edges': len(visible_edges), 'routes': visible_edges.route.nunique(),
            'open_circle_meaning': 'earliest observed route source; human original observed time root retained',
            'mouse_node_identity': rooted_summary['node_identity'],
            'observed_route_starts': len(routes),
            'state_label_trajectory_marker_collisions': marker_collisions,
            'trajectory_input': str(args.rooted_trajectories),
            'old_display_edges_are_historical': True})
    pd.DataFrame(anchors).to_csv(output/'state_label_positions.tsv', sep='\t', index=False)
    if args.smooth_trajectories:
        pd.DataFrame(smoothing_audit).to_csv(output/'smoothing_chain_audit.tsv', sep='\t', index=False)
        pd.DataFrame(curve_vertices).to_csv(output/'smoothed_trajectory_path_vertices.tsv', sep='\t', index=False)
        review.update({'trajectory_display_smoothing': 'bounded simplification plus quadratic corner rounding',
            'smoothing_locked_nodes': 'all roots, forks and terminals',
            'smoothing_changes_ot_edges': False,
            'smoothing_max_deviation_mm': max(a['max_sampled_deviation_mm'] for a in smoothing_audit),
            'smoothing_deviation_upper_bound_mm': max(a['deviation_upper_bound_mm'] for a in smoothing_audit),
            'smoothing_original_waypoints': sum(a['original_vertices'] for a in smoothing_audit),
            'smoothing_retained_waypoints': sum(a['retained_waypoints'] for a in smoothing_audit),
            'intermediate_ot_node_coordinates_are_data_not_curve_interpolation_constraints': True})
    pd.DataFrame([{**item, 'displayed_ot_edges': int(visible_edges.route.eq(item['route']).sum())}
                  for item in routes]).to_csv(output/'trajectory_legend.tsv', sep='\t', index=False)
    (visible_edges if rooted_summary else display_edges).to_csv(
        output/'direct_raw_cell_ot_arrows_display_coordinates.tsv', sep='\t', index=False)
    visible_edges.to_csv(output/'rendered_ot_edges.tsv', sep='\t', index=False)
    suppressed.to_csv(output/'suppressed_isolated_short_arrows.tsv', sep='\t', index=False)
    time_audit.to_csv(output/'time_resolved_node_audit.tsv', sep='\t', index=False)
    mouse_times.to_csv(output/'frozen_mouse_temporal_display_edges.tsv', sep='\t', index=False)
    for name in ['stateflow_cross_species_group_contract.tsv', 'state_group_summary.tsv']:
        table = pd.read_csv(source/name, sep='\t')
        table['display_label'] = table.state_group.map(display_labels)
        table['legend_label'] = table.state_group.map(legend_labels)
        table.to_csv(output/name, sep='\t', index=False)
    if rooted_summary:
        visible_edges.to_csv(output/'direct_raw_cell_ot_arrows.tsv', sep='\t', index=False)
        shutil.copyfile(source/'direct_raw_cell_ot_arrows.tsv', output/'legacy_direct_raw_cell_ot_arrows.tsv')
        display_edges.to_csv(output/'legacy_time_resolved_display_edges.tsv', sep='\t', index=False)
    else:
        shutil.copyfile(source/'direct_raw_cell_ot_arrows.tsv', output/'direct_raw_cell_ot_arrows.tsv')
    shutil.copyfile(source/'direct_raw_cell_ot_arrows_display_coordinates.tsv',
                    output/'legacy_state_collapsed_display_coordinates.tsv')
    (output/'summary.json').write_text(json.dumps(review, indent=2)+'\n')
    print(json.dumps(review, indent=2))


if __name__ == '__main__':
    main()
