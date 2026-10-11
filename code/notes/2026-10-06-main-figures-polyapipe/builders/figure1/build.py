"""Figure 1 · records adapter of the v5.1 production builder (layout unchanged; data from one result bundle).

Usage: build.py --records <bundle>/result-records.json --out <new dir>
Panel notes below describe the production figure drawn from the archived v2 data; their numbers are history, not checks.

Base: locked master Fig/Figure.1.pdf. Contract rows 1A–1D of main-panel-replacement-v5.tsv.
  A, B  kept as in the master (frozen expression atlas and state relations).
  C  current gene-centred DWUI (figure_scores/metacell_scores centered_DWUI) on the frozen cell embedding;
     state letters re-placed at the state-group medians of the frozen coordinates.
  D  top: mean ΔDWUI of R/E/P/M relative to observed S per route (supported50_p08, depth 20, all measured genes);
     bottom: source-qualified P–S means — open circle, available genes (supported50_p08); tick, 932 shared
     one-to-one orthologs (same setting); grey bar, range of the mean over the 16 settings in
     route_direction_audit (4 PAS settings x 4 gene backgrounds). No method tag; "Setting range" legend.
"""
import csv, gzip, json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools-local'))
from v5draw import *  # noqa
from figrec import Run, MASTERS
HERE = Path(__file__).resolve().parent
R = Run('F1')
import numpy as np
import pandas as pd
from matplotlib import patheffects as pe
from matplotlib.colors import LinearSegmentedColormap

MASTER = MASTERS/'Figure.1.pdf'      # locked master copy (Fig/ now holds promoted production files)
assert sha(MASTER) == 'ad248e74a0b3d782e6b4a03d44136bdb2f4f37283563d375d300664a25f99f63'
OUT = R.out/'Figure.1.draft.pdf'
def rd(panel, table, **kw):
    return R.df(panel, table, **kw)


pl = pdfplumber.open(MASTER).pages[0]; W, H = pl.width, pl.height
pdf = pikepdf.open(MASTER)
# 1A atlas block: the master leaves ~25 pt empty under the embedding. The block (raster, route paths and
# start markers) is placed again at 1.10x and centred in the space above the route legend; the 14 state
# labels are removed from the copied block and redrawn at their scaled centres at the master size (8 pt bold).
UB = (14.5, 5.3, 282.1, 233.8)                  # atlas image box, top-origin
ALAB = [('S', 202.3, 130.7, 207.7, 138.7), ('E', 75.3, 146.5, 80.6, 154.5), ('R', 116.7, 168.7, 122.4, 176.7),
        ('P', 118.4, 57.8, 123.7, 65.8), ('mB', 26.0, 88.9, 38.9, 96.9), ('K', 168.2, 125.2, 174.0, 133.2),
        ('M', 164.7, 183.8, 171.4, 191.8), ('mI', 208.8, 66.7, 218.2, 74.7), ('mN', 135.0, 82.2, 147.9, 90.2),
        ('hX', 92.0, 102.2, 102.2, 110.2), ('mY', 252.2, 100.0, 264.6, 108.0), ('mU', 184.3, 219.8, 197.2, 227.8),
        ('hN', 96.8, 37.9, 107.5, 45.9), ('hU', 133.8, 39.8, 144.5, 47.8)]
SA = 1.04
CONT = (15.0, 7.0, 262.0, 228.0)               # visible content box of the atlas (cells and labels), top-origin
TGT = (8.0, 4.5, 303.0, 246.0)                 # space between the panel edge, the family legend and the route legend
ax0 = TGT[0]+((TGT[2]-TGT[0])-(CONT[2]-CONT[0])*SA)/2-(CONT[0]-UB[0])*SA
ay0 = TGT[1]+((TGT[3]-TGT[1])-(CONT[3]-CONT[1])*SA)/2-(CONT[1]-UB[1])*SA
# round 4: rescale the r10 placement about its measured visible content box (28, 8)-(294, 246) so that the top gap
# (~9 pt) matches the left gap (~12 pt) and the family legend (x 307) gets a ~32 pt gap; bottom stays above the route legend (y 252)
CN, CNEW, R_ = (28.0, 8.0, 294.0, 246.0), (16.5, 15.5), 0.96
ax0 = CNEW[0]+(ax0-CN[0])*R_; ay0 = CNEW[1]+(ay0-CN[1])*R_; SA = SA*R_
AX = lambda x: ax0+(x-UB[0])*SA; AY = lambda y: ay0+(y-UB[1])*SA
CBOX = (CNEW[0], CNEW[1], CNEW[0]+(CN[2]-CN[0])*R_, CNEW[1]+(CN[3]-CN[1])*R_)
assert CBOX[2] < 307.0-28.0 and CBOX[3] < 252.0-6.0, CBOX
blk = pikepdf.open(MASTER)
nlab = len(clear(blk, H, [(x0-0.6, y0-0.6, x1+0.6, y1+0.6) for _, x0, y0, x1, y1 in ALAB]))
assert nlab >= len(ALAB), nlab
# keep only the atlas block in the copy (everything else would ride along hidden behind the clip)
nout = len(clear(blk, H, [(-20, -20, UB[0]-0.05, H+20), (UB[2]+0.05, -20, W+20, H+20), (-20, -20, W+20, UB[1]-0.05),
                           (-20, UB[3]+0.05, W+20, H+20)]))
prune_xobjects(blk)
bp = blk.pages[0]
box = [UB[0], H-UB[3], UB[2], H-UB[1]]
bp.MediaBox = box; bp.CropBox = box
CLEAR = {
    'C': [(29.0, 481.8, 108.8, 570.2), (129.7, 483.9, 216.4, 570.5), (115.0, 585.0, 210.0, 595.5)],
    'D': [(227.5, 331.6, 453.5, 599.8), (235.6, 320.5, 453.5, 331.8)],          # everything except the panel letter
}
removed = {k: len(clear(pdf, H, v)) for k, v in CLEAR.items()}
removed['A_block'] = len(clear(pdf, H, [(UB[0]-0.5, UB[1]-0.5, UB[2]+0.5, UB[3]+0.5)]))
pdf.pages[0].add_overlay(bp, pikepdf.Rectangle(AX(UB[0]), H-AY(UB[3]), AX(UB[2]), H-AY(UB[1])))
# round 5: the human route legend moves 12 pt right (level with the embedding's left edge, closer to the mouse legend)
RL, RDX = (0.0, 250.0, 222.0, 312.0), 12.0
blk2 = pikepdf.open(MASTER)
clear(blk2, H, [(-20, -20, RL[0]-0.05, H+20), (RL[2]+0.05, -20, W+20, H+20), (-20, -20, W+20, RL[1]-0.05), (-20, RL[3]+0.05, W+20, H+20)])
prune_xobjects(blk2)
bp2 = blk2.pages[0]; bx2 = [RL[0], H-RL[3], RL[2], H-RL[1]]; bp2.MediaBox = bx2; bp2.CropBox = bx2
removed['A_route_legend'] = len(clear(pdf, H, [RL]))
pdf.pages[0].add_overlay(bp2, pikepdf.Rectangle(RL[0]+RDX, H-RL[3], RL[2]+RDX, H-RL[1]))
removed['A_legend_title'] = len(clear(pdf, H, [(306.0, 22.5, 450.5, 31.5)]))      # 'h = human; m = mouse' goes to the caption
t = Top(W, H); ver = {'master_sha256': sha(MASTER), 'removed_objects': removed}
t.text_top(307.1, 23.0, 'Display families', 8)
for lab, x0, y0, x1, y1 in ALAB:
    t.s.ax.text(AX((x0+x1)/2), H-AY((y0+y1)/2), lab, fontproperties=FONTS['Bol'], fontsize=8, color=INK, ha='center', va='center',
                zorder=6, path_effects=[pe.withStroke(linewidth=1.6, foreground='white')])
ver['A'] = dict(scale=SA, visible_content_box=[round(v, 1) for v in CBOX], image_box_top_origin=[round(v, 2) for v in (AX(UB[0]), AY(UB[1]), AX(UB[2]), AY(UB[3]))],
                content_box_top_origin=[round(v, 2) for v in (AX(CONT[0]), AY(CONT[1]), AX(CONT[2]), AY(CONT[3]))], labels_removed_from_copy=nlab)
# the 1C colour bar in the master is kept; its map, sampled at 600 dpi, is #377eb8 -> #f2f1f0 -> #cc8026
CMAP = LinearSegmentedColormap.from_list('d1c', ['#377eb8', '#f2f1f0', '#cc8026'])
PROX, DIST = '#377eb8', '#cc8026'   # master 1D sign colours

# ------------------------------------------------------------------ C
cells = rd('F1C', 'frozen/frozen_cell_embedding.tsv.gz')
ms = rd('F1C', 'figure_scores/metacell_scores.tsv.gz')
alias = {'WOT_GSE115943': 'GSE115943'}
cells['_s'] = cells.study_id.replace(alias); ms['_s'] = ms.study_id.replace(alias)
j_ = cells.merge(ms[['metacell_id', '_s', 'centered_DWUI']], on=['metacell_id', '_s'], how='left', validate='many_to_one', indicator=True)
assert len(j_) == len(cells) and j_.cell_id.tolist() == cells.cell_id.tolist()
BOX = {'human': (29.6, 108.2, 482.3, 569.6), 'mouse': (130.2, 215.9, 484.4, 570.0)}   # master embedding boxes; master titles kept
LET = dict(shared_s='S', shared_r='R', shared_e='E', shared_tp='P', shared_m='M')
LIM = 0.03; Cd = {}
for sp in ['human', 'mouse']:
    q = j_[j_.species_key == sp]
    x0, x1, y0, y1 = BOX[sp]
    ux0, ux1 = q.UMAP1.min(), q.UMAP1.max(); uy0, uy1 = q.UMAP2.min(), q.UMAP2.max()
    kx = (x1-x0)/(ux1-ux0); ky = (y1-y0)/(uy1-uy0)
    assert abs(kx/ky-1) < 0.02, (sp, kx, ky)
    k_ = min(kx, ky)
    PX = lambda u: x0+(u-ux0)*k_+((x1-x0)-(ux1-ux0)*k_)/2
    PY = lambda u: H-(y1-(u-uy0)*k_-((y1-y0)-(uy1-uy0)*k_)/2)
    val = q.centered_DWUI.to_numpy(float); ok = np.isfinite(val)
    t.ax.scatter(PX(q.UMAP1.to_numpy()), PY(q.UMAP2.to_numpy()), s=0.04, c='#DADFE3', linewidths=0, rasterized=True, zorder=2)
    o = np.argsort(np.abs(val[ok]))
    t.ax.scatter(PX(q.UMAP1.to_numpy()[ok][o]), PY(q.UMAP2.to_numpy()[ok][o]), s=0.04, c=val[ok][o], cmap=CMAP,
                 vmin=-LIM, vmax=LIM, linewidths=0, rasterized=True, zorder=3)
    let = {}
    for grp, L in LET.items():
        m_ = q.state_group == grp
        lx, ly = PX(float(q.UMAP1[m_].median())), PY(float(q.UMAP2[m_].median()))
        t.s.ax.text(lx, ly, L, fontproperties=FONTS['Bol'], fontsize=8, color=INK, ha='center', va='center', zorder=6,
                    path_effects=[pe.withStroke(linewidth=1.6, foreground='white')])
        let[L] = [round(lx, 1), round(H-ly, 1)]
    if ok.sum() and (np.abs(val[ok]) > LIM).mean() > 0.05:
        R.flag('F1C', 'colour_saturation', f'{sp}: {100*(np.abs(val[ok]) > LIM).mean():.1f}% of scored cells beyond the locked ±{LIM} colour bar')
    if not ok.any():
        R.unavailable('F1C', f'{sp}: no metacell score joined')
    Cd[sp] = dict(cells=int(len(q)), cells_with_score=int(ok.sum()), clipped=int((np.abs(val[ok]) > LIM).sum()),
                  score_range=[float(np.nanmin(val)), float(np.nanmax(val))], letters_top_origin=let)
for xc, lab in zip((128.8, 163.1, 197.4), ['−0.03', '0', '0.03']):
    t.text_top(xc, 585.9, lab, 8, ha='center')
ver['C'] = Cd

# ------------------------------------------------------------------ D top: R/E/P/M relative to observed S
srs = rd('F1D', 'geometry/state-reference-summary.tsv')
srs = srs[(srs.background == 'available') & srs.state.isin(list('REPM'))].copy()
assert not srs.duplicated(['species', 'route', 'state']).any()
srs['n_'] = np.where(srs.status == 'measured', srs.n.fillna(0), 0).astype(int)
mean = srs.pivot(index=['species', 'route'], columns='state', values='mean').reindex(columns=list('REPM'))
cnt = srs.pivot(index=['species', 'route'], columns='state', values='n_').reindex(columns=list('REPM')).fillna(0).astype(int)
# layout: one row pitch (11.5 pt) in both halves; species groups marked by a rotated label and bracket, with a gap between groups
P_ = 11.5; LX = 239.0; GAP = 4.0
def bracket(y0_, y1_, lab_):
    t.line([235.2, 235.2], [y0_, y1_], lw=0.6, color=MUTED)
    t.vtext(230.6, (y0_+y1_)/2, lab_, 8)
t.text_top(338.0, 322.2, 'Mean ΔDWUI relative to observed S', 8, ha='center')
CX = dict(R=344.5, E=373.0, P=401.5, M=434.0); KX = 12.55/0.03
# legend: bar length = ±0.03
lgx = 356.0
t.text_top(lgx-14.6, 334.6, '−0.03', 7.2, ha='right'); t.text_top(lgx-14.6-textwidth('−0.03', 7.2)-8.0, 334.6, 'Proximal', 7.2, ha='right')
t.rect(lgx-12.55, 335.2, 12.55, 5.5, fc=PROX); t.rect(lgx, 335.2, 12.55, 5.5, fc=DIST)
t.text_top(lgx+14.6, 334.6, '+0.03', 7.2); t.text_top(lgx+14.6+textwidth('+0.03', 7.2)+8.0, 334.6, 'Distal', 7.2)
for k, c_ in [('R', '#0072B2'), ('E', '#009E73'), ('P', '#CC79A7'), ('M', '#8C6D31')]:
    t.text_top(CX[k], 345.4, k, 8, ha='center')
    t.line([CX[k]-12.5, CX[k]+12.5], [355.0, 355.0], lw=1.6, color=c_)
UP = [('human', 'GSE118258', 'Xing2020', 0), ('human', 'GSE150311', 'Liu2020*', 0), (None, None, 'Guan2022', 0),
      ('human', '178325_0330', '0330 Embryonic fibroblasts', 1), ('human', '178325_0605', '0605 Adult skin fibroblasts', 1),
      ('human', '178325_0618', '0618 ASC time course', 1), ('human', '178325_0809', '0809 ASC start/end', 1),
      ('human', 'GSE221739', 'Panariello2023', 0), ('mouse', 'OSK', 'Guo2019–OSK', 0), ('mouse', 'CIP', 'Guo2019–CIP', 0),
      ('mouse', 'WOT', 'Schiebinger2019–2i', 0)]
yc = 362.5; ROWS = []; yh = {}
for sp, r, lab, ind in UP:
    if sp == 'mouse' and 'mouse' not in yh:
        yc += GAP
    t.text_top(LX+2.7*ind, yc-3.6, lab, 7.2)
    assert LX+2.7*ind+textwidth(lab, 7.2) < CX['R']-12.5-2.0, lab
    if sp:
        ROWS.append((sp, r, yc)); yh.setdefault(sp, []).append(yc)
    else:
        yh.setdefault('human', []).append(yc)
    yc += P_
for sp, lab in [('human', 'Human'), ('mouse', 'Mouse')]:
    bracket(min(yh[sp])-4.5, max(yh[sp])+4.5, lab)
ytop_end = max(yh['mouse'])+4.5
t.line([(CX['P']+CX['M'])/2]*2, [352.0, ytop_end+1.0], lw=0.5, color='#c5ced4')
Dtop = {}
for sp, r, y in ROWS:
    row = {}
    if (sp, r) not in cnt.index:
        R.unavailable('F1D', f'top row {sp}|{r}: route absent from state-reference summary (row left empty)')
        Dtop[f'{sp}|{r}'] = row; continue
    for k in 'REPM':
        if cnt.loc[(sp, r), k] == 0:
            continue
        v = float(mean.loc[(sp, r), k])
        if not R.axis('F1D', f'top {sp}|{r}|{k} mean ΔDWUI vs S', [v], -0.03, 0.03):
            v = max(-0.03, min(0.03, v))
        t.rect(CX[k]-12.5, y-3.7, 25.1, 7.4, fc='#eef1f3')
        t.rect(min(CX[k], CX[k]+v*KX), y-2.75, abs(v)*KX, 5.5, fc=DIST if v > 0 else PROX)
        t.line([CX[k]]*2, [y-4.3, y+4.3], lw=0.3, color='#8c9da7')
        row[k] = [round(v, 5), int(cnt.loc[(sp, r), k])]
    Dtop[f'{sp}|{r}'] = row
# every route with any S-relative node is drawn, except Jaber2022 (M only; not in the master grid either)
drawn = {(sp, r) for sp, r, _ in ROWS}
extra = sorted(k for k in cnt.index if cnt.loc[k].sum() > 0 and k not in drawn)
for k in extra:
    if cnt.loc[k][list('REP')].sum() > 0:
        R.flag('F1D', 'row_capacity', f'top: {k[0]}|{k[1]} has R/E/P estimates but no row in the locked grid')
ver['D_top'] = dict(rows=Dtop, not_drawn_M_only=[r for _, r in extra], row_pitch=P_)

# ------------------------------------------------------------------ D bottom: source-qualified P–S
ysep = ytop_end+5.5
t.line([228.0, 448.0], [ysep, ysep], lw=0.5, color='#c5ced4')
gds = rd('F1D', 'global/global_direction_summary.tsv')
aud = rd('F1D', 'axis_extensions/route_direction_audit.tsv')
n_orth = set()
yt = ysep+3.6
y0b = yt+16.5
PS = [('human', 'GSE118258', 'Xing2020', y0b), ('human', '178325_0618', 'Guan2022–0618', y0b+P_), ('human', '178325_0809', 'Guan2022–0809', y0b+2*P_),
      ('mouse', 'CIP', 'Guo2019–CIP', y0b+3*P_+GAP), ('mouse', 'WOT', 'Schiebinger2019–2i', y0b+4*P_+GAP)]
FA0, FA1 = 306.0, 388.0; KF = (FA1-FA0)/0.08; Z = FA0+0.05*KF; FX = lambda v: Z+v*KF
EA0, EA1, E0, E1 = 395.0, 448.0, -0.35, 0.25; EX = lambda v: EA0+(v-E0)/(E1-E0)*(EA1-EA0)
t.text_top(LX, yt, 'Documented-source P–S mean ΔDWUI', 8)
eh = 'Distal excess, %'
assert textwidth(eh, 7.2) <= EA1-EA0+1.0 and LX+textwidth('Documented-source P–S mean ΔDWUI', 8) < EA0-3
t.text_top((EA0+EA1)/2, yt+0.6, eh, 7.2, ha='center')
for sp, lab in [('human', 'Human'), ('mouse', 'Mouse')]:
    ys_ = [y for s_, _, _, y in PS if s_ == sp]; bracket(min(ys_)-4.5, max(ys_)+4.5, lab)
ylast = PS[-1][3]
t.line([EX(0)]*2, [y0b-6.0, ylast+5.5], lw=0.4, color='#8c9da7', ls=(0, (1, 1.6)))
t.line([Z, Z], [y0b-6.0, ylast+5.5], lw=0.4, color='#8c9da7', ls=(0, (1, 1.6)))
Dbot = {}
clipF = lambda v: min(max(v, -0.05+1e-9), 0.03-1e-9)
clipE = lambda v: min(max(v, E0), E1)
for sp, r, lab, y in PS:
    t.text_top(LX, y-3.6, lab, 7.2)
    assert LX+textwidth(lab, 7.2) < FA0-2.0, lab
    g_ = gds[(gds.route == r) & (gds.contrast == 'P_S') & (gds.design == 'true_source') & (gds.setting == 'supported50_p08')]
    av = g_[(g_.background == 'available') & (g_.status == 'measured')]; orth = g_[(g_.background == 'common_one_to_one_ortholog_all_settings') & (g_.status == 'measured')]
    a_ = aud[(aud.route == r) & (aud.contrast == 'P_S') & (aud.design == 'true_source')]
    if len(av) != 1 or len(orth) != 1 or len(a_) != 1:
        R.unavailable('F1D', f'bottom row {r}: available/ortholog/setting-range rows {len(av)}/{len(orth)}/{len(a_)}; row left empty')
        continue
    av, orth, a_ = av.iloc[0], orth.iloc[0], a_.iloc[0]
    n_orth.add(int(orth.n))
    if int(a_.scenarios) != 16 or int(a_.empty_scenarios) != 0:
        R.flag('F1D', 'caption_value', f'{r}: setting range covers {int(a_.scenarios)} settings ({int(a_.empty_scenarios)} empty); caption says 16')
    lo, hi = float(a_.mean_min), float(a_.mean_max); lo_raw, hi_raw = lo, hi
    for v in (float(av['mean']), float(orth['mean'])):
        if not lo-1e-12 <= v <= hi+1e-12:
            R.flag('F1D', 'consistency', f'{r}: primary mean {v:.4g} outside its own setting range {lo:.4g}..{hi:.4g}')
    if not R.axis('F1D', f'bottom {r} P–S mean / setting range', [lo, hi, float(av['mean']), float(orth['mean'])], -0.05, 0.03):
        lo, hi = clipF(lo), clipF(hi)
    t.line([FX(lo), FX(hi)], [y, y], lw=2.1, color='#8c9da7')
    for v in (lo, hi):
        t.line([FX(v)]*2, [y-2.3, y+2.3], lw=0.6, color='#8c9da7')
    va, vo = clipF(float(av['mean'])), clipF(float(orth['mean']))
    t.line([FX(vo), FX(va)], [y, y], lw=0.35, color='#8c9da7')
    t.line([FX(vo)]*2, [y-2.8, y+2.8], lw=1.3, color=INK)
    t.marker(FX(va), y, 'o', 3.6, fc='white', ec=DIST if va > 0 else PROX, mew=0.65, zorder=5)
    ex_ = float(av['distal_excess_002'])
    exr = gds[(gds.route == r) & (gds.contrast == 'P_S') & (gds.design == 'true_source') & (gds.status == 'measured')].distal_excess_002.astype(float)
    if len(exr) != 16:
        R.flag('F1D', 'caption_value', f'{r}: distal-excess range over {len(exr)} settings; caption says 16')
    if not R.axis('F1D', f'bottom {r} distal excess', list(exr)+[ex_], E0, E1):
        exr = exr.clip(E0, E1); ex_ = clipE(ex_)
    t.line([EX(exr.min()), EX(exr.max())], [y, y], lw=2.1, color='#c5ced4')
    t.rect(min(EX(0), EX(ex_)), y-2.4, abs(EX(ex_)-EX(0)), 4.8, fc=DIST if ex_ > 0 else PROX)
    Dbot[r] = dict(available_mean=va, available_n=int(av.n), ortholog932_mean=vo, setting_range=[lo_raw, hi_raw],
                   distal_excess_002=ex_, distal_excess_range16=[float(exr.min()), float(exr.max())],
                   proximal_002=int(av.proximal_002), distal_002=int(av.distal_002))
PS = [(sp, r, y) for sp, r, _, y in PS]
hum = [Dbot[r]['setting_range'] for _, r, _ in PS[:3] if r in Dbot]; mou = [Dbot[r]['setting_range'] for _, r, _ in PS[3:] if r in Dbot]
R.claim('F1D', 'every human route lies above every mouse route in all settings',
        bool(hum and mou) and min(v[0] for v in hum) > max(v[1] for v in mou))
if len(n_orth) != 1:
    R.flag('F1D', 'consistency', f'common-ortholog background size differs between routes: {sorted(n_orth)}')
n_orth_lab = f'{min(n_orth)} orthologs' if n_orth else 'Common orthologs'
ya = ylast+8.5
t.line([FX(-0.05), FX(0.03)], [ya]*2, lw=0.4, color='#8c9da7')
t.line([EA0, EA1], [ya]*2, lw=0.4, color='#8c9da7')
for v, lab in [(-0.3, '−30'), (0.0, '0'), (0.2, '+20')]:  # polyApipe redraw 2026-10-06: axis rescaled within the same span
    t.line([EX(v)]*2, [ya, ya+2.4], lw=0.4, color='#8c9da7'); t.text_top(EX(v), ya+4.1, lab, 7.2, ha='center')
for v, lab in [(-0.04, '−0.04'), (-0.02, '−0.02'), (0.0, '0'), (0.02, '+0.02')]:
    t.line([FX(v)]*2, [ya, ya+2.4], lw=0.4, color='#8c9da7')
    t.text_top(FX(v), ya+4.1, lab, 7.2, ha='center')
yl = ya+15.0; x_ = LX
t.marker(x_+2.0, yl+3.6, 'o', 3.6, fc='white', ec=INK, mew=0.65); x_ += 6.5; t.text_top(x_, yl, 'Available genes', 7.2); x_ += textwidth('Available genes', 7.2)+10.0
t.line([x_+1.0]*2, [yl+0.8, yl+6.4], lw=1.3, color=INK); x_ += 4.5; t.text_top(x_, yl, n_orth_lab, 7.2); x_ += textwidth(n_orth_lab, 7.2)+10.0
t.line([x_, x_+9.0], [yl+3.6]*2, lw=2.1, color='#8c9da7'); x_ += 12.0; t.text_top(x_, yl, 'Setting range', 7.2); x_ += textwidth('Setting range', 7.2)
assert x_ < 450.0 and yl+7.5 < 598.0, (x_, yl)
ver['D_bottom'] = Dbot

t.save(R.out/'overlay.pdf'); overlay(pdf, R.out/'overlay.pdf', 0, 0)
prune_xobjects(pdf); pdf.save(OUT)
render_png(OUT, R.out/'Figure.1.draft.png', 150)
ver['output_sha256'] = sha(OUT)
R.finish(1, OUT, json.loads(json.dumps(ver, default=str)))
