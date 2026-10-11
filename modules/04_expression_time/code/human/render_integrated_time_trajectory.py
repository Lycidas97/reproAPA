"""One mouse-style human overview using frozen cells and registered time relations."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path('results_temporal_v1')
OUT = ROOT / 'trajectory_display'
PALETTE = {'GSE118258': '#188e8a', 'GSE150311': '#dc8641', 'GSE178325': '#6468bd', 'GSE221739': '#bc4c78'}
FONT_FILES = ['.figure-fonts/Arial.ttf',
              '/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf']
for path in FONT_FILES:
    fontManager.addfont(path)
FONT = FontProperties(family=['Arial', 'Droid Sans Fallback'])
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'path'})


def canonical(route, time):
    return 'gse178325_0330_shared' if route.startswith('gse178325_0330_') and time == 'source' else route


def key(route, time, condition):
    return canonical(route, time) + '|' + time + '|' + condition


def main():
    OUT.mkdir(exist_ok=True)
    coordinates_path = ROOT / 'layout/time_coordinates.npy'
    coordinate_hash = hashlib.sha256(coordinates_path.read_bytes()).hexdigest()
    assert coordinate_hash == 'ea6ad7cfe00c76e60b1a709fbc110e688b9126addeb0dc5dbfab633ee4dfbc93'
    obs = pd.read_csv(ROOT / 'layout/time_coordinates.tsv.gz', sep='\t', dtype=str).fillna('')
    xy = np.load(coordinates_path)
    assert len(obs) == 215114 and obs.cell_id.is_unique and obs.original_state_id.nunique() == 62
    assert np.allclose(obs[['UMAP1', 'UMAP2']].to_numpy(float), xy, atol=1e-6)
    nodes = pd.read_csv(ROOT / 'layout/observed_time_nodes.tsv', sep='\t')
    edges = pd.read_csv(ROOT / 'layout/observed_time_display_edges.tsv', sep='\t')
    plan = pd.read_csv('results_v1/prepared/temporal_edge_plan.tsv', sep='\t')
    assert len(plan) == edges.pair.nunique() == 34
    nodes['node_id'] = [key(r.route_id, r.time, r.condition) for r in nodes.itertuples()]
    # The four 0330 comparisons share the exact same source library/cells.
    for _, g in nodes.groupby('node_id'):
        assert g[['cells', 'UMAP1', 'UMAP2']].drop_duplicates().shape[0] == 1
    nodes = nodes.drop_duplicates('node_id').set_index('node_id')
    nodes['study_id'] = np.where(nodes.route_id.str.startswith('gse178325'), 'GSE178325', nodes.route_id)
    edges['source'] = [key(r.route_id, r.source_time, r.source_condition) for r in edges.itertuples()]
    edges['target'] = [key(r.route_id, r.target_time, r.target_condition) for r in edges.itertuples()]
    checks = []
    for r in edges.itertuples():
        contract = plan.iloc[r.pair]
        assert contract.route_id == r.route_id and contract.time_from == r.source_time and contract.time_to == r.target_time
        for side, time, condition in [('source', r.source_time, r.source_condition), ('target', r.target_time, r.target_condition)]:
            ix = np.flatnonzero(obs.library_key.isin(json.loads(contract[side + '_libraries'])) & obs.condition.eq(condition))
            center = np.median(xy[ix], axis=0)
            n = nodes.loc[key(r.route_id, time, condition)]
            assert len(ix) == n.cells and np.allclose(center, [n.UMAP1, n.UMAP2], atol=1e-6)
        assert np.allclose([r.x0, r.y0, r.x1, r.y1], [nodes.loc[r.source].UMAP1, nodes.loc[r.source].UMAP2, nodes.loc[r.target].UMAP1, nodes.loc[r.target].UMAP2])
        checks.append(dict(pair=r.pair, route=r.route_id, source=r.source, target=r.target, contract_and_centers_match=True))
    nodes.to_csv(OUT / 'display_nodes.tsv', sep='\t')
    edges.to_csv(OUT / 'display_temporal_arrows.tsv', sep='\t', index=False)
    pd.DataFrame(checks).to_csv(OUT / 'display_checks.tsv', sep='\t', index=False)

    fig, ax = plt.subplots(figsize=(11, 12))
    fig.subplots_adjust(left=.045, right=.955, bottom=.14, top=.875)
    order = np.random.default_rng(20260909).permutation(len(obs))
    colors = obs.study_id.map(PALETTE).to_numpy()
    cloud = ax.scatter(*xy[order].T, c=colors[order], s=.55, alpha=.24, linewidths=0, rasterized=True)
    assert cloud.get_offsets().shape[0] == 215114
    for r in edges.itertuples():
        a, b = nodes.loc[r.source], nodes.loc[r.target]
        style, width, bend = '-', 1.9, .035
        if '0605' in r.route_id:
            style, width, bend = '--', 1.7, .06
        elif '0330' in r.route_id:
            style, width, bend = ':', 1.4, .03
        ax.annotate('', xy=(b.UMAP1, b.UMAP2), xytext=(a.UMAP1, a.UMAP2),
                    arrowprops=dict(arrowstyle='-|>', color=PALETTE[a.study_id], lw=width, linestyle=style,
                                    mutation_scale=15, shrinkA=3, shrinkB=4, connectionstyle=f'arc3,rad={bend}'), zorder=5)
    for _, r in nodes.iterrows():
        ax.scatter(r.UMAP1, r.UMAP2, s=17, facecolor='white', edgecolor=PALETTE[r.study_id], linewidth=.9, zorder=6)

    # Labels are observed sample/condition landmarks; intermediate nodes and arrows remain visible.
    landmarks = [
        ('GSE118258', 'D0', '', 'OSKM D0', (18, 12)),
        ('GSE118258', 'D8', '', 'OSKM D8', (16, -20)),
        ('GSE118258', 'D16', 'TRA-1-60_negative', 'D16 TRA-1-60−', (22, 8)),
        ('GSE118258', 'D16', 'TRA-1-60_positive', 'D16 TRA-1-60+', (12, 13)),
        ('GSE221739', 'D3', '', 'mRNA D3', (12, -27)),
        ('GSE221739', 'D9', '', 'mRNA D9', (-18, -17)),
        ('GSE221739', 'D15', '', 'mRNA D15', (18, 12)),
        ('gse178325_0618', 'source', '', '0618 起始', (15, -22)),
        ('gse178325_0618', 'Stage_I_day_8', '', '0618 阶段 I', (12, -26)),
        ('gse178325_0618', 'Stage_II_day_20', '', '阶段 II', (-22, -20)),
        ('gse178325_0618', 'Stage_III_day_12', '', '阶段 III', (-22, 0)),
        ('gse178325_0618', 'Stage_IV_day_4', '', '阶段 IV', (-26, 2)),
        ('gse178325_0618', 'endpoint', '', '0618 终点', (-25, 18)),
        ('gse178325_0605', 'source', '', '0605 起始', (15, -22)),
        ('gse178325_0605', 'Stage_II', '', '0605 阶段 II', (-22, -25)),
        ('gse178325_0330_5azac', 'source', '', '0330 共享起始', (5, -32)),
        ('gse178325_0330_unperturbed', 'Stage_II', '', '0330 阶段 II\n四种处理分支', (10, 22)),
    ]
    label_artists, label_rows = [], []
    for route, time, condition, label, offset in landmarks:
        selected = nodes[(nodes.index.str.startswith(canonical(route, time) + '|' + time + '|'))]
        if condition:
            selected = selected[selected.condition.eq(condition)]
        assert len(selected) == 1, (route, time, condition)
        node_id, r = next(selected.iterrows())
        text = ax.annotate(label, (r.UMAP1, r.UMAP2), xytext=offset, textcoords='offset points',
                ha='left' if offset[0] >= 0 else 'right', va='center', fontproperties=FONT, fontsize=10,
                color=PALETTE[r.study_id], linespacing=1.2,
                bbox=dict(fc='white', ec='none', alpha=.88, pad=2),
                arrowprops=dict(arrowstyle='-', lw=.6, color=PALETTE[r.study_id]), zorder=8)
        label_artists.append(text)
        label_rows.append(dict(node_id=node_id, label=label, offset_x=offset[0], offset_y=offset[1]))
    pd.DataFrame(label_rows).to_csv(OUT / 'display_labels.tsv', sep='\t', index=False)
    low, high = xy.min(axis=0), xy.max(axis=0)
    pad = .065 * (high - low)
    ax.set_xlim(low[0] - pad[0], high[0] + pad[0]); ax.set_ylim(low[1] - pad[1], high[1] + pad[1])
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    names = {'GSE118258': 'GSE118258 · OSKM', 'GSE221739': 'GSE221739 · mRNA',
             'GSE178325': 'GSE178325 · 化学重编程', 'GSE150311': 'GSE150311 · 表达整合'}
    handles = [Line2D([], [], color=PALETTE[s], lw=2.5, label=names[s]) for s in ['GSE118258', 'GSE221739', 'GSE178325', 'GSE150311']]
    legend = fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .945), ncol=2,
                        frameon=False, prop=FontProperties(family=FONT.get_family(), size=11), columnspacing=2.5)
    title = fig.suptitle('人类重编程：整合数据集中的主要时间轨迹', fontproperties=FONT, fontsize=18, y=.985)
    foot = fig.text(.5, .074, '215,114 个单细胞 · 4 个研究 · 34 对时间 / 阶段关系\n'
        '节点：观测时间与条件中心；箭头：研究内前向关系。化学路线：0618 实线 / 0605 虚线 / 0330 点线。\n'
        '阶段 I–IV 不等于统一采样日；跨空白连线不代表观测到连续中间细胞。',
        ha='center', va='center', fontproperties=FONT, fontsize=10, linespacing=1.6)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_boxes = [a.get_bbox_patch().get_window_extent(renderer) for a in label_artists]
    overlaps = []
    for i, a in enumerate(text_boxes):
        for j in range(i+1, len(text_boxes)):
            if a.overlaps(text_boxes[j]):
                overlaps.append([label_rows[i]['label'], label_rows[j]['label']])
    for name, artist in [('title', title), ('legend', legend), ('foot', foot)]:
        bbox = artist.get_window_extent(renderer)
        assert bbox.x0 > 0 and bbox.x1 < fig.bbox.width and bbox.y0 > 0 and bbox.y1 < fig.bbox.height, name
    fig.savefig(OUT / 'integrated_time_trajectory.png', dpi=300, facecolor='white')
    fig.savefig(OUT / 'integrated_time_trajectory.pdf', dpi=300, facecolor='white')
    fig.savefig(OUT / 'integrated_time_trajectory.svg', dpi=300, facecolor='white')
    plt.close(fig)
    pixels = np.asarray(Image.open(OUT / 'integrated_time_trajectory.png').convert('RGB'))
    mask = np.min(pixels, axis=2) < 250
    yy, xx = np.where(mask)
    margins = dict(left_px=int(xx.min()), right_px=int(pixels.shape[1]-1-xx.max()),
                   top_px=int(yy.min()), bottom_px=int(pixels.shape[0]-1-yy.max()))
    assert hashlib.sha256(coordinates_path.read_bytes()).hexdigest() == coordinate_hash
    result = dict(status='COMPUTED_PENDING_VISUAL_REVIEW', cells=len(obs), studies=obs.study_id.nunique(), original_states=62,
        observed_time_condition_nodes=len(nodes), displayed_arrows=len(edges), temporal_pairs=34, labels=len(label_rows),
        all_cells_displayed=True, coordinates_changed=False, coordinate_sha256=coordinate_hash,
        coordinate_limits_include_all_cells=True, duplicate_0330_sources_merged=3,
        nodes='median of registered time and condition; shared 0330 source displayed once',
        arrows='all registered forward time/condition comparisons; constant-width summary, not inferred lineage probabilities',
        style_reference='mouse render_final_time_trajectory.py: full point cloud, arrows, sparse landmark labels',
        label_overlaps=overlaps, png_pixels=list(pixels.shape[:2][::-1]), rendered_margins=margins,
        fonts=FONT_FILES, font_substitution='Arial TrueType replaces Nimbus Sans L CFF because the latter triggered a PDF font-type mismatch; Droid Sans Fallback covers Chinese',
        exports=['integrated_time_trajectory.png','integrated_time_trajectory.pdf','integrated_time_trajectory.svg'])
    (OUT / 'summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
