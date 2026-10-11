"""Matched layouts and observed-time summaries; arrows are OT relations, not lineages."""
import argparse
from pathlib import Path
import json
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import yaml

PALETTE = {'GSE118258': '#188e8a', 'GSE150311': '#dc8641', 'GSE178325': '#6468bd', 'GSE221739': '#bc4c78'}


def clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def short_time(label):
    return label.replace('Stage_', '').replace('_day_', ':').replace('source', 'start').replace('endpoint', 'end')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    args = p.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    root = Path(cfg['output_root'])
    out = root / 'figures'
    out.mkdir(exist_ok=True)
    a = ad.read_h5ad(root / 'layout/human_time_constrained.h5ad')
    obs = a.obs
    xy = a.obsm['X_umap']
    order = np.random.default_rng(cfg['seed']).permutation(a.n_obs)
    study = obs.study_id.astype(str).to_numpy()
    plan = pd.read_csv(Path(cfg['input_root']) / 'prepared/temporal_edge_plan.tsv', sep='\t').fillna('')
    colors = np.array([PALETTE[x] for x in study])
    fig, axes = plt.subplots(2, 3, figsize=(18, 11), layout='constrained')
    for ax, key, title in zip(axes[0], ['X_umap_expression', 'X_umap_control', 'X_umap'],
            ['Previous expression layout', 'Matched control: expression + cross-study MNN', 'Time constraints: expression + OT + MNN']):
        coords = a.obsm[key]
        ax.scatter(*coords[order].T, c=colors[order], s=.35, linewidths=0, rasterized=True)
        ax.set_title(title, fontsize=11)
        clean(ax)
    for ax, name in zip(axes[1, :2], ['GSE118258', 'GSE221739']):
        mask = study == name
        ax.scatter(*xy[order].T, c='#e6e6e6', s=.15, linewidths=0, rasterized=True)
        ix = order[mask[order]]
        im = ax.scatter(*xy[ix].T, c=pd.to_numeric(obs.iloc[ix].actual_day), cmap='viridis', s=.5, linewidths=0, rasterized=True)
        fig.colorbar(im, ax=ax, shrink=.65, label='Observed day within this study')
        ax.set_title(name + ('; D16 sorted strata retained' if name == 'GSE118258' else '; D3-D15 only'), fontsize=11)
        clean(ax)
    ax = axes[1, 2]
    ax.scatter(*xy[order].T, c='#e6e6e6', s=.15, linewidths=0, rasterized=True)
    mask = study == 'GSE178325'
    stage = pd.to_numeric(obs.stage_order.astype(str), errors='coerce').to_numpy()
    ix = order[mask[order] & np.isfinite(stage[order])]
    im = ax.scatter(*xy[ix].T, c=stage[ix], cmap='viridis', vmin=1, vmax=4, s=.5, linewidths=0, rasterized=True)
    fig.colorbar(im, ax=ax, ticks=[1, 2, 3, 4], shrink=.65, label='Stage order, not elapsed days')
    ax.set_title('GSE178325: ordered chemical stages', fontsize=11)
    clean(ax)
    fig.legend(handles=[Line2D([], [], marker='o', ls='', color=v, label=k) for k, v in PALETTE.items()],
               loc='outside upper center', ncol=4, frameon=False)
    fig.supxlabel('215,114 primary cells; all coordinates free during fitting. GSE150311 contributes expression only.', fontsize=10)
    fig.savefig(out / 'time_constraint_overview.png', dpi=180)
    plt.close(fig)

    # Every node below summarizes a registered observed time and capture condition.
    node_rows, edge_rows = [], []
    for i, row in plan.iterrows():
        sides = []
        for side, time_key in [('source', 'time_from'), ('target', 'time_to')]:
            mask = obs.library_key.isin(json.loads(row[side + '_libraries'])).to_numpy()
            groups = []
            for condition in sorted(obs.loc[mask, 'condition'].astype(str).unique()):
                ix = np.flatnonzero(mask & obs.condition.astype(str).eq(condition).to_numpy())
                node = dict(route_id=row.route_id, time=row[time_key], condition=condition, cells=len(ix),
                            UMAP1=float(np.median(xy[ix, 0])), UMAP2=float(np.median(xy[ix, 1])))
                node_rows.append(node); groups.append(node)
            sides.append(groups)
        for s in sides[0]:
            for t in sides[1]:
                edge_rows.append(dict(pair=i, route_id=row.route_id, source_time=s['time'], target_time=t['time'],
                    source_condition=s['condition'], target_condition=t['condition'],
                    x0=s['UMAP1'], y0=s['UMAP2'], x1=t['UMAP1'], y1=t['UMAP2']))
    nodes = pd.DataFrame(node_rows).drop_duplicates(['route_id', 'time', 'condition'])
    edges = pd.DataFrame(edge_rows)
    nodes.to_csv(root / 'layout/observed_time_nodes.tsv', sep='\t', index=False)
    edges.to_csv(root / 'layout/observed_time_display_edges.tsv', sep='\t', index=False)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), layout='constrained')
    panels = [(['GSE118258'], 'GSE118258: actual D0-D16'), (['GSE221739'], 'GSE221739: actual D3-D15'),
              (['gse178325_0618'], 'Chemical 0618: within-route stage order'),
              (['gse178325_0605'], 'Chemical 0605: start -> I -> II'),
              (list(plan.loc[plan.route_id.str.startswith('gse178325_0330'), 'route_id']), 'Chemical 0330: shared start, four conditions')]
    for ax, (routes, title) in zip(axes.flat, panels):
        relevant = plan[plan.route_id.isin(routes)]
        libraries = set()
        for values in relevant[['source_libraries', 'target_libraries']].to_numpy().ravel():
            libraries.update(json.loads(values))
        ix = np.flatnonzero(obs.library_key.isin(libraries))
        ax.scatter(*xy[ix].T, c='#c9c9d8', s=.4, linewidths=0, alpha=.5, rasterized=True)
        selected = nodes[nodes.route_id.isin(routes)].drop_duplicates(['time', 'condition', 'UMAP1', 'UMAP2'])
        for _, edge in edges[edges.route_id.isin(routes)].iterrows():
            ax.annotate('', xy=(edge.x1, edge.y1), xytext=(edge.x0, edge.y0),
                        arrowprops=dict(arrowstyle='->', color='#343850', lw=1.1, alpha=.8))
        for i, (_, node) in enumerate(selected.iterrows()):
            ax.scatter(node.UMAP1, node.UMAP2, s=24, color='#0c807c', edgecolor='white', lw=.5, zorder=4)
            label = short_time(node.time)
            if title.startswith('Chemical 0330') and node.time != 'source':
                label = node.route_id.split('0330_')[-1]
            if 'D16' in node.time:
                label += '+' if 'positive' in node.condition.lower() else '-'
            offset = [(8, 10), (8, -14), (-8, 10), (-8, -14)][i % 4]
            offset = {('GSE118258', 'D0'): (10, -14), ('GSE118258', 'D2'): (10, 10),
                      ('gse178325_0618', 'Stage_III_day_8'): (-12, 12),
                      ('gse178325_0618', 'Stage_III_day_12'): (12, -14)}.get((node.route_id, node.time), offset)
            ax.annotate(label, (node.UMAP1, node.UMAP2), xytext=offset,
                        textcoords='offset points', fontsize=8, color='#222222', zorder=5,
                        ha='left' if offset[0] > 0 else 'right',
                        bbox=dict(facecolor='white', alpha=.7, edgecolor='none', pad=.2),
                        arrowprops=dict(arrowstyle='-', color='#777777', lw=.5))
        lower, upper = np.quantile(xy[ix], [.005, .995], axis=0)
        lower = np.minimum(lower, selected[['UMAP1', 'UMAP2']].min().to_numpy())
        upper = np.maximum(upper, selected[['UMAP1', 'UMAP2']].max().to_numpy())
        padding = .1 * np.maximum(upper - lower, .1)
        ax.set_xlim(lower[0] - padding[0], upper[0] + padding[0])
        ax.set_ylim(lower[1] - padding[1], upper[1] + padding[1])
        ax.set_title(title, fontsize=11)
        ax.margins(.10)
        clean(ax)
    ax = axes.flat[5]
    ax.axis('off')
    ax.text(.03, .95, 'How to read these panels\n\n'
        'Dots: medians of observed time / condition groups.\n'
        'Arrows: registered forward OT comparisons.\n'
        'Full cell-to-cell couplings retain distributed paths.\n\n'
        'Arrow segments are not observed intermediate cells\n'
        'and do not establish a unique lineage.\n\n'
        'Chemical I:0 / I:0.5 / II:4 labels are\n'
        'within-stage days, not a shared global clock.\n\n'
        'D16+ / D16- are sorted captures; arrow mass\n'
        'does not estimate reprogramming success.\n\n'
        'GSE150311 has no temporal OT edges.', va='top', fontsize=11, linespacing=1.5)
    fig.suptitle('Observed-time relationships on the full-cell time-constrained layout', fontsize=15)
    fig.supxlabel('Route panels zoom to central 99% coordinate ranges and retain all observed-time nodes; the overview shows all cells.', fontsize=10)
    fig.savefig(out / 'observed_time_relations.png', dpi=190)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(15, 12), layout='constrained')
    for ax, study_name in zip(axes.flat, PALETTE):
        ax.scatter(*xy[order].T, c='#ececec', s=.12, linewidths=0, rasterized=True)
        ix = order[study[order] == study_name]
        labels = pd.Categorical(obs.iloc[ix].cluster_id).remove_unused_categories()
        assert len(labels.categories) == obs.iloc[ix].original_state_id.nunique()
        ax.scatter(*xy[ix].T, c=labels.codes, cmap='gist_ncar', s=.5, linewidths=0, rasterized=True)
        for cluster in labels.categories:
            rows = ix[obs.iloc[ix].cluster_id.astype(str).eq(str(cluster)).to_numpy()]
            center = np.median(xy[rows], axis=0)
            ax.text(*center, str(cluster), fontsize=7, ha='center',
                    bbox=dict(boxstyle='round,pad=.1', facecolor='white', alpha=.7, edgecolor='none'))
        ax.set_title(f'{study_name}: {len(labels.categories)} unchanged original clusters')
        clean(ax)
    fig.savefig(out / 'original_states_time_layout.png', dpi=180)
    plt.close(fig)
    print('Three figure panels and observed-time display tables saved', flush=True)


if __name__ == '__main__':
    main()
