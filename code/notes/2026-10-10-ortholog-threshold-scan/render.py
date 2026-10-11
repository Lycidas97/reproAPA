"""Single drawing path for the complete planned directional threshold scan."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for path in (ROOT / 'notes/2026-10-02-main-figure-revision/tools/fonts').glob('*.otf'):
    font_manager.fontManager.addfont(str(path))
font = font_manager.findfont('Nimbus Sans L', fallback_to_default=False)
INK = '#263D4A'
plt.rcParams.update({'font.family': 'Nimbus Sans L', 'font.size': 8, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'pdf.fonttype': 3,
    'axes.linewidth': .6, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'axes.spines.right': False, 'axes.spines.top': False,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})
data = pd.read_csv(HERE / 'directional_tests.tsv', sep='\t', float_precision='round_trip')
assert len(data) == 144 and data.tested.all()
assert set(data.threshold) == {.01, .015, .02, .025, .03, .035, .04, .045, .05}
order = [('proximal', 'proximal'), ('distal', 'distal'),
         ('distal', 'proximal'), ('proximal', 'distal')]
models = [('full_assay', '#0072B2', '-', 'o'),
          ('measurement_PAS_capacity', '#586870', '--', 's')]
fold_max = float(data.fold.max())
top = max(2., np.ceil(fold_max*1.27*2)/2)
fig = plt.figure(figsize=(160/25.4, 210/25.4))
fig.text(.30875, .982, 'P–S   n = 1656', ha='center', va='top')
fig.text(.77125, .982, 'P–M   n = 751', ha='center', va='top')
fig.legend(handles=[Line2D([], [], color=color, ls=ls, marker=marker, markersize=3,
    label='Full assay' if model == 'full_assay' else 'Measurement / PAS / capacity')
    for model, color, ls, marker in models], loc='upper center', bbox_to_anchor=(.5, .960),
    frameon=False, ncol=2, fontsize=7, columnspacing=1, handlelength=2.2)
fig.text(.5, .920, 'Row direction   human / mouse        n   observed pairs        *   q72 < 0.05',
         ha='center', va='bottom', fontsize=7)
verification, axes = [], []
for i, (human, mouse) in enumerate(order):
    for j, contrast in enumerate(['P_S', 'P_M']):
        x, y = [.10, .5625][j], [.735, .525, .315, .105][i]
        ax = fig.add_axes([x, y, .4175, .165])
        axes.append(ax)
        identity = f'{human.capitalize()} / {mouse}'
        ax.text(.5, 1.022, identity, transform=ax.transAxes, ha='center', va='bottom', fontsize=8)
        fig.text(x-.046, y+.173, chr(ord('A')+2*i+j), fontweight='bold', fontsize=11,
                 ha='left', va='top')
        part = data[data.contrast.eq(contrast) & data.human_direction.eq(human) &
                    data.mouse_direction.eq(mouse)]
        counts = None
        for model, color, linestyle, marker in models:
            q = part[part.model.eq(model)].sort_values('threshold')
            assert len(q) == 9 and q.measured_orthologs.eq(1656 if contrast == 'P_S' else 751).all()
            if counts is None:
                counts = q[['threshold', 'observed']].reset_index(drop=True)
            else:
                assert counts.equals(q[['threshold', 'observed']].reset_index(drop=True))
            ax.plot(q.threshold, q.fold, color=color, ls=linestyle, lw=.8,
                    marker=marker, markersize=3.2, markeredgewidth=.6, zorder=2)
            for row in q.itertuples():
                if np.isfinite(row.fold) and row.q_scan_72_per_model < .05:
                    offset = .022*top if model == 'full_assay' else -.055*top
                    ax.text(row.threshold, row.fold+offset, '*', fontsize=8, color=color,
                            ha='center', va='bottom', zorder=3)
                elif not np.isfinite(row.fold):
                    ax.text(row.threshold, .025*top, 'NA', color='#586870', fontsize=7, ha='center')
                verification.append(dict(panel=chr(ord('A')+2*i+j), contrast=contrast,
                    human_direction=human, mouse_direction=mouse, model=model,
                    threshold=float(row.threshold), observed=int(row.observed),
                    analytic_expected=float(row.analytic_expected), fold=float(row.fold),
                    q72=float(row.q_scan_72_per_model), star=bool(row.q_scan_72_per_model < .05)))
        for row in counts.itertuples():
            ax.text(row.threshold, .94*top, str(int(row.observed)), ha='center', va='center', fontsize=7)
        ax.axhline(1, color='#B4B9BC', ls=':', lw=.7, zorder=0)
        ax.set(xlim=(.0085, .0515), ylim=(0, top), xticks=[.01, .03, .05],
               yticks=np.linspace(0, np.floor(top), 4))
        if j == 0:
            ax.set_ylabel('Observed / expected')
        if i == 3:
            ax.set_xlabel('|ΔDWUI| threshold')
        else:
            ax.tick_params(labelbottom=False)
        ax.tick_params(width=.6)
fig.text(.5, .027, 'Fixed measured backgrounds and class margins   100000 permutations per context',
         fontsize=7, ha='center', va='center')
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
canvas = fig.bbox
outside = []
for item in fig.findobj(match=Text):
    if item.get_visible() and item.get_text():
        box = item.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > canvas.x1+1 or box.y1 > canvas.y1+1:
            outside.append(item.get_text())
assert not outside, outside
fig.savefig(HERE / 'Threshold-scan.pdf')
fig.savefig(HERE / 'Threshold-scan.png', dpi=150)
fig.savefig(HERE / 'Threshold-scan-2x.png', dpi=300)
# Bounding boxes are retained for visual QA; this does not infer scientific validity.
boxes = [dict(panel=chr(ord('A')+j), axes_fraction=list(ax.get_position().bounds))
         for j, ax in enumerate(axes)]
ink_boxes = [ax.get_tightbbox(renderer) for ax in axes]
scale_mm = 25.4 / fig.dpi
gutters = [dict(left=chr(ord('A')+2*i), right=chr(ord('B')+2*i),
    visible_gutter_mm=float((ink_boxes[2*i+1].x0-ink_boxes[2*i].x1)*scale_mm)) for i in range(4)]
vertical = [dict(above=chr(ord('A')+j), below=chr(ord('A')+j+2),
    visible_gutter_mm=float((ink_boxes[j].y0-ink_boxes[j+2].y1)*scale_mm)) for j in range(6)]
assert all(row['visible_gutter_mm'] > 0 for row in gutters + vertical)
(HERE / 'figure-verification.json').write_text(json.dumps(dict(status='PASS',
    source_sha256=hashlib.sha256((HERE / 'directional_tests.tsv').read_bytes()).hexdigest(),
    source_rows=144, plotted_rows=len(verification), source_values=verification,
    font_path=font, no_canvas_text_clipping=True, panel_boxes=boxes,
    horizontal_ink_gutters=gutters, vertical_ink_gutters=vertical,
    width_mm=160, height_mm=210, matplotlib=matplotlib.__version__), indent=2) + '\n')
for box in ink_boxes:
    bounds = box.transformed(fig.transFigure.inverted()).bounds
    fig.add_artist(Rectangle(bounds[:2], bounds[2], bounds[3], transform=fig.transFigure,
                             fill=False, edgecolor='#D55E00', linewidth=.6))
fig.savefig(HERE / 'Threshold-scan-panel-boxes.png', dpi=150)
