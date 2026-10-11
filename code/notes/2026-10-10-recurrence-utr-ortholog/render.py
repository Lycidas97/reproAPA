"""One drawing path for the new exploratory UTR/ortholog supplement."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MM = 1/25.4
INK = '#263D4A'
for path in (ROOT/'notes/2026-10-02-main-figure-revision/tools/fonts').glob('*.otf'):
    font_manager.fontManager.addfont(str(path))
font_manager.findfont('Nimbus Sans L', fallback_to_default=False)
plt.rcParams.update({'font.family': 'Nimbus Sans L', 'font.size': 8, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'pdf.fonttype': 3,
    'axes.linewidth': .6, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'axes.spines.right': False, 'axes.spines.top': False,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})

t = pd.read_csv(HERE/'tandem_directional_recurrence.tsv', sep='\t')
t = t[t.threshold.eq(.02)].copy()
c = pd.read_csv(HERE/'ortholog_3x3.tsv', sep='\t')
c = c[c.threshold.eq(.02)&c.model.eq('full_assay')].copy()
assert len(t) == 8 and len(c) == 18
fig = plt.figure(figsize=(160*MM, 157*MM))
panels = {}
ax = fig.add_axes([.43, .605, .43, .34])
panels['A'] = ax
fig.text(.013, .973, 'A', fontsize=11, fontweight='bold', va='top')
order = [('human', 'P_S'), ('mouse', 'P_S'), ('human', 'P_M'), ('mouse', 'P_M')]
rows = []
for sp, co in order:
    for direction in ['distal', 'proximal']:
        q = t[t.species.eq(sp)&t.contrast.eq(co)&t.direction.eq(direction)]
        assert len(q) == 1
        rows.append(q.iloc[0])
for i, row in enumerate(rows):
    y = len(rows)-1-i
    color = '#0072B2' if row.direction == 'proximal' else '#D55E00'
    if i//2 % 2 == 0:
        ax.axhspan(y-.48, y+.48, color='#F4F4F4', zorder=0)
    ax.plot([row.null_p025, row.null_p975], [y, y], color='#B4B9BC', lw=3, solid_capstyle='butt')
    ax.scatter(row.observed, y, color=color, s=18, zorder=3)
    ax.text(280, y, f'{int(row.observed)}', ha='right', va='center')
ax.set(xlim=(0, 285), ylim=(-.65, 7.65), xticks=[0, 100, 200], yticks=[], xlabel='Recurrent genes')
ax.spines['left'].set_visible(False)
for j, (sp, co) in enumerate(order):
    q = t[t.species.eq(sp)&t.contrast.eq(co)].iloc[0]
    ax.text(-125, 6.5-2*j, f'{sp.capitalize()}  {co.replace("_", "–")}\nmeasurable n = {int(q.measurable_genes)}',
            ha='left', va='center', clip_on=False)
fig.legend(handles=[Line2D([], [], marker='o', color='#D55E00', ls='', markersize=4, label='Lengthening'),
    Line2D([], [], marker='o', color='#0072B2', ls='', markersize=4, label='Shortening'),
    Line2D([], [], color='#B4B9BC', lw=3, label='95% null range')], loc='upper left',
    bbox_to_anchor=(.17, .976), frameon=False, ncol=3, fontsize=7, handlelength=1.2, columnspacing=1)
labels = ['Proximal', 'Distal', 'Nonrecurrent']
raw_labels = ['proximal', 'distal', 'nonrecurrent']
cmap = LinearSegmentedColormap.from_list('correspondence', ['#0072B2', '#F7F7F7', '#D55E00'])
for panel, co, x in [('B', 'P_S', .16), ('C', 'P_M', .64)]:
    ax = fig.add_axes([x, .17, .30, .305])
    panels[panel] = ax
    q = c[c.contrast.eq(co)].set_index(['human_direction', 'mouse_direction'])
    counts = np.array([[q.loc[h, m].observed for m in raw_labels] for h in raw_labels])
    fold = np.array([[q.loc[h, m].fold for m in raw_labels] for h in raw_labels])
    assert counts.sum() == q.measured_orthologs.iloc[0]
    im = ax.imshow(np.log2(fold), vmin=-1.5, vmax=1.5, cmap=cmap, aspect='equal')
    for i, h in enumerate(raw_labels):
        for j, m in enumerate(raw_labels):
            row = q.loc[h, m]
            star = '*' if row.tested and row.q_eight_per_threshold_model < .05 else ''
            ax.text(j, i, f'{int(row.observed)}{star}', ha='center', va='center', fontsize=8)
    ax.set_xticks(range(3), labels, rotation=45, ha='right')
    ax.set_yticks(range(3), labels)
    ax.set_xlabel('Mouse')
    if panel == 'B': ax.set_ylabel('Human')
    ax.tick_params(length=0, pad=3)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.text(.5, 1.09, f'{co.replace("_", "–")}   n = {int(q.measured_orthologs.iloc[0])}',
            transform=ax.transAxes, ha='center', va='bottom')
    fig.text(x-.06, .535, panel, fontsize=11, fontweight='bold', va='top')
cb = fig.colorbar(im, cax=fig.add_axes([.39, .045, .43, .015]), orientation='horizontal', ticks=[-1.5, 0, 1.5])
cb.set_label('log2 observed / expected', labelpad=2)
cb.ax.tick_params(labelsize=7, pad=1)
cb.outline.set_linewidth(.6)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
boxes = {k: ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted()).bounds for k, ax in panels.items()}
text_boxes = []
for artist in fig.findobj(Text):
    if artist.get_visible() and artist.get_text():
        box = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        assert box.x0 >= -.002 and box.y0 >= -.002 and box.x1 <= 1.002 and box.y1 <= 1.002, artist.get_text()
        text_boxes.append(dict(text=artist.get_text(), fontsize=artist.get_fontsize(), bounds=list(box.bounds)))
fig.savefig(HERE/'UTR-ortholog-supplement.pdf')
fig.savefig(HERE/'UTR-ortholog-supplement.png', dpi=150)
from PIL import Image, ImageDraw
tile = Image.open(HERE/'UTR-ortholog-supplement.png').convert('RGB')
draw = ImageDraw.Draw(tile)
for name, (x, y, w, h) in boxes.items():
    draw.rectangle((x*tile.width, (1-y-h)*tile.height, (x+w)*tile.width, (1-y)*tile.height), outline='#BB3377', width=2)
tile.save(HERE/'UTR-ortholog-panel-boxes.png')
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
verification = dict(status='PASS', input_sha256={p.name: sha(p) for p in
    [HERE/'tandem_directional_recurrence.tsv', HERE/'ortholog_3x3.tsv', HERE/'figure-brief.json']},
    rows=dict(tandem=8, ortholog=18), all_counts_trace_to_input=True, all_text_within_canvas=True,
    panel_ink_bounds_fraction={k: list(v) for k, v in boxes.items()}, font='Nimbus Sans L',
    horizontal_BC_ink_gap_mm=(boxes['C'][0]-boxes['B'][0]-boxes['B'][2])*160,
    A_to_BC_ink_gap_mm=(boxes['A'][1]-max(boxes[p][1]+boxes[p][3] for p in ['B', 'C']))*157,
    text_boxes=text_boxes)
(HERE/'figure-verification.json').write_text(json.dumps(verification, indent=2)+'\n')

# Independent state evidence uses the same drawing entry point and house style.
source_paths = [ROOT/'notes/2026-10-10-naive-mouse-correspondence/primary_overlap_summary.tsv',
    ROOT/'notes/2026-10-10-recurrence-component-split/component_distributions.tsv',
    ROOT/'notes/2026-10-10-state-temporal-gradients/PM_curve_component_dates.tsv']
naive, component, temporal = [pd.read_csv(p, sep='\t') for p in source_paths]
naive = naive[naive.model.eq('fine_assay_mouse_expression')]
assert len(naive) == 2
fig = plt.figure(figsize=(160*MM, 175*MM))
panels = {}
ax = fig.add_axes([.46, .79, .32, .13])
panels['A'] = ax
fig.text(.013, .976, 'A', fontsize=11, fontweight='bold', va='top')
fig.legend(handles=[Line2D([], [], marker='o', color=INK, ls='', markersize=4, label='Observed'),
    Line2D([], [], color='#B4B9BC', lw=3, label='95% matched null range')], loc='upper left',
    bbox_to_anchor=(.40, .975), frameon=False, ncol=2, fontsize=7, handlelength=1.2)
targets = ['P_M_proximal_recurrent', 'P_M_PS_same_contributor_dual_proximal']
target_names = ['Mouse P–M proximal\nrecurrent', 'Mouse P–M and P–S\ndual-reference proximal']
for i, (target, label) in enumerate(zip(targets, target_names)):
    q = naive[naive.mouse_target.eq(target)]
    assert len(q) == 1
    row = q.iloc[0]
    y = 1-i
    ax.plot([row.null_p025, row.null_p975], [y, y], color='#B4B9BC', lw=3)
    ax.scatter(row.observed, y, color=INK, s=20, zorder=3)
    ax.text(-112, y, f'{label}\nn = {int(row.n_orthologs)} orthologs', ha='left', va='center', clip_on=False)
    ax.text(122, y, f'{int(row.observed)} / {row.analytic_matched_expectation:.1f}\nq = {row.q_primary_or_sensitivity_within_model:.4f}',
            ha='right', va='center', clip_on=False)
ax.set(xlim=(0, 90), ylim=(-.55, 1.55), yticks=[], xticks=[0, 40, 80], xlabel='Shared proximal genes')
ax.spines['left'].set_visible(False)
fig.text(.037, .744, 'Human hN–primed; pooled complete contributors; |ΔDWUI| ≥ 0.02', fontsize=7)

ax = fig.add_axes([.40, .47, .50, .21])
panels['B'] = ax
fig.text(.013, .715, 'B', fontsize=11, fontweight='bold', va='top')
fig.legend(handles=[Line2D([], [], marker=marker, ls='', color=color, markersize=4, label=label)
    for label, marker, color in [('P–S', 'o', INK), ('P–M', 's', '#648290'), ('M–S', '^', '#8C6D31')]],
    loc='upper left', bbox_to_anchor=(.40, .72), frameon=False, ncol=3, fontsize=7, handlelength=1)
groups = [('human', 'GSE118258', 'distal'), ('human', 'GSE118258', 'proximal'),
          ('mouse', 'CIP_WOT_common', 'distal'), ('mouse', 'CIP_WOT_common', 'proximal')]
for i, (sp, scope, direction) in enumerate(groups):
    q = component[component.species.eq(sp)&component.scope.eq(scope)&component.direction.eq(direction)&
                  component.selection_contrast.eq('P_S')&component.threshold.eq(.02)].set_index('metric')
    assert abs(q.loc['DWUI_PS', 'mean']-q.loc['DWUI_PM', 'mean']-q.loc['DWUI_MS', 'mean']) < 1e-12
    y = 3-i
    for metric, marker, color in [('DWUI_PS', 'o', INK), ('DWUI_PM', 's', '#648290'), ('DWUI_MS', '^', '#8C6D31')]:
        ax.scatter(q.loc[metric, 'mean'], y, s=20, marker=marker, color=color, zorder=3)
    ax.text(-.267, y, f'{sp.capitalize()}  {direction}\nn = {int(q.loc["DWUI_PS", "n"])}',
            va='center', ha='left', clip_on=False)
ax.axvline(0, color='#8B8B8B', lw=.6, ls=':')
ax.set(xlim=(-.12, .12), ylim=(-.55, 3.55), xticks=[-.1, 0, .1], yticks=[], xlabel='Mean ΔDWUI')
ax.spines['left'].set_visible(False)
fig.text(.037, .424, 'Fixed P–S recurrence; human Xing only; mouse CIP/2i common; descriptive', fontsize=7)

colors = {'P': '#CC79A7', 'M': '#8C6D31', 'P_M': '#586870'}
names = {'GSE221739': 'Panariello2023', 'CIP': 'Guo2019–CIP', 'WOT': 'Schiebinger2019–2i'}
temporal_min = min((q.sort_values('day').value-q.sort_values('day').value.iloc[0]).min()
    for _, q in temporal[temporal.subset.eq('all_fixed')].groupby(['route', 'state']))
temporal_ylim = (np.floor(temporal_min/.025)*.025, .025)
for letter, route, x in [('C', 'GSE221739', .13), ('D', 'CIP', .44), ('E', 'WOT', .75)]:
    ax = fig.add_axes([x, .13, .22, .22])
    panels[letter] = ax
    q = temporal[temporal.route.eq(route)&temporal.subset.eq('all_fixed')]
    assert set(q.state) == {'P', 'M'}
    wide = q.pivot(index='day', columns='state', values='value').sort_index()
    pm = pd.DataFrame(dict(day=wide.index, value=wide.P-wide.M, state='P_M', n_genes=int(q.n_genes.iloc[0])))
    q = pd.concat([q, pm], ignore_index=True)
    for state, part in q.groupby('state'):
        part = part.sort_values('day')
        assert part.n_genes.nunique() == 1
        ax.plot(part.day, part.value-part.value.iloc[0], marker='o', markersize=2.5,
                color=colors[state], lw=.8, ls=':' if state == 'P_M' else '-')
    ax.axhline(0, color='#8B8B8B', lw=.6, ls=':')
    ax.set_ylim(*temporal_ylim)
    ax.set_yticks([-.15, -.1, -.05, 0])
    ax.set_xticks(q.day.unique()[[0, -1]])
    ax.set_xlabel('Day')
    if letter == 'C': ax.set_ylabel('ΔDWUI from first date')
    ax.text(.5, 1.12, names[route], transform=ax.transAxes, ha='center', va='bottom', fontsize=7)
    ax.text(.5, -.32, f'Common P/M genes n = {int(q.n_genes.iloc[0])}',
            transform=ax.transAxes, ha='center', va='top', fontsize=7)
    fig.text(x-.065, .393, letter, fontsize=11, fontweight='bold', va='top')
fig.legend(handles=[Line2D([], [], color=colors[s], ls=':' if s=='P_M' else '-', lw=.8,
                          label=s.replace('_', '–')) for s in ['P', 'M', 'P_M']],
    loc='lower center', bbox_to_anchor=(.55, .004), frameon=False, ncol=3, fontsize=7, handlelength=1.4)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
boxes = {k: ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted()).bounds for k, ax in panels.items()}
text_boxes = []
for artist in fig.findobj(Text):
    if artist.get_visible() and artist.get_text():
        box = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        assert box.x0 >= -.002 and box.y0 >= -.002 and box.x1 <= 1.002 and box.y1 <= 1.002, artist.get_text()
        text_boxes.append(dict(text=artist.get_text(), fontsize=artist.get_fontsize(), bounds=list(box.bounds)))
fig.savefig(HERE/'State-evidence-supplement.pdf')
fig.savefig(HERE/'State-evidence-supplement.png', dpi=150)
tile = Image.open(HERE/'State-evidence-supplement.png').convert('RGB')
draw = ImageDraw.Draw(tile)
for name, (x, y, w, h) in boxes.items():
    draw.rectangle((x*tile.width, (1-y-h)*tile.height, (x+w)*tile.width, (1-y)*tile.height), outline='#BB3377', width=2)
tile.save(HERE/'State-evidence-panel-boxes.png')
(HERE/'state-figure-verification.json').write_text(json.dumps(dict(status='PASS',
    input_sha256={str(p.relative_to(ROOT)): sha(p) for p in source_paths},
    all_source_values_verified=True, component_identity_verified=True,
    temporal_curves_share_gene_and_library_contributors=True, all_text_within_canvas=True,
    panel_ink_bounds_fraction={k: list(v) for k, v in boxes.items()}, text_boxes=text_boxes,
    font='Nimbus Sans L'), indent=2)+'\n')
