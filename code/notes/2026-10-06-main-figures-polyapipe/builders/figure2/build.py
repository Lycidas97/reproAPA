"""Figure 2 · records adapter of the v5.1 production builder (layout unchanged; data from one result bundle).

Usage: build.py --records <bundle>/result-records.json --out <new dir>
Panel notes below describe the production figure drawn from the archived v2 data; their numbers are history, not checks.

Base: locked master Fig/Figure.2.pdf. Contract rows 2A–2F of main-panel-replacement-v5.tsv.
  A  current k = 3 classes (human 807 = 264/147/396; mouse 485 = 241/148/96): route-mean S-relative curves and
     gene x (state, route) heatmaps, supported50_p08 depth 20. Master grammar kept (same routes, columns, boxes).
     Routes without any S-relative estimate keep their legend entry, set in grey.
  B  L/D/T schematic now drawn from the C2 k = 3 shape centre (early proximal), plus P(T >= t) per route x class.
  C  mean step ΔDWUI cells, mean stage contributions to L (contribution_*_mean), per-route median D/L with IQR.
  D  class-centred DWUI on the frozen embedding (current k3 metacell scores), state letters at state medians,
     studies with Moran q < .05 (k 20, state + coverage + library) beside each class: 17/21.
  E  leave-route-out: 3 x 3 held-out centroid r (RMSE beneath) + assignment agreement / ARI vs full and
     remaining-route fits. Old APA–GEX correlations moved to the supplement.
  F  class proportions of the candidate genes per species, with gene counts.
"""
import csv, gzip, json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools-local'))
from v5draw import *  # noqa
from figrec import Run, MASTERS
HERE = Path(__file__).resolve().parent
R = Run('F2')
import numpy as np
import pandas as pd
from matplotlib import patheffects as pe
from matplotlib.patches import Polygon

MASTER = MASTERS/'Figure.2.pdf'      # locked master copy (Fig/ now holds promoted production files)
assert sha(MASTER) == '9fd2854a37baeb71c5d141cda8fa2f8e3568361b3dbce9f91f6d0fa256db9a70'
OUT = R.out/'Figure.2.draft.pdf'
def rd(panel, table, **kw):
    return R.df(panel, table, **kw)

pl = pdfplumber.open(MASTER).pages[0]; W, H = pl.width, pl.height
pdf = pikepdf.open(MASTER)
HX = 223.1   # human offset of the A curve axes and heatmaps relative to mouse
GREYED = [(162.0, 171.0, 246.5, 180.0), (312.4, 182.2, 402.6, 191.0), (305.2, 193.4, 365.7, 202.2),
          (65.2, 204.5, 125.7, 213.3), (145.0, 204.5, 223.7, 213.3)]
CLEAR = {
    'A': [(21.5, 27.5, 81.0, 54.0), (21.5, 60.6, 81.0, 87.0), (21.5, 93.6, 81.0, 120.0),
          (21.5+HX, 27.5, 81.0+HX, 54.0), (21.5+HX, 60.6, 81.0+HX, 87.0), (21.5+HX, 93.6, 81.0+HX, 120.0),
          (0, 28, 21.4, 121), (224, 28, 244.4, 121),
          (88.5, 27.5, 183.0, 121.0), (311.7, 27.5, 406.2, 121.0),
          (186.0, 42.5, 202.0, 51.5), (186.0, 75.5, 202.0, 84.5), (186.0, 108.5, 202.0, 117.5),
          (409.5, 42.5, 425.0, 51.5), (409.5, 75.5, 425.0, 84.5), (409.5, 108.5, 425.0, 117.5),
          (428.0, 140.0, 446.0, 149.5), (428.0, 170.5, 446.0, 179.5), (428.0, 210.5, 450.0, 219.5), (420.0, 140.5, 426.0, 221.0), (407.5, 152.0, 417.0, 198.5),
          (0.0, 136.5, 405.0, 216.0), (186.0, 8.5, 218.5, 30.5), (409.0, 10.5, 441.0, 30.5)],
    'B': [(14.0, 233.0, 232.0, 372.0), (0.0, 250.0, 14.0, 372.0)],
    'C': [(265.0, 245.5, 452.0, 320.2), (425.0, 234.5, 442.0, 244.0), (413.0, 320.0, 452.0, 329.5)],
    'D': [(4.5, 392.5, 69.0, 456.7), (79.0, 392.5, 143.3, 456.7), (153.5, 392.5, 217.8, 456.7),
          (241.8, 391.0, 300.9, 456.3), (316.3, 391.0, 375.4, 456.3), (390.7, 391.0, 449.9, 456.3),
          (0.0, 459.0, 453.5, 491.0), (228.0, 387.5, 400.0, 398.5)],
    'E': [(0, 515.0, 330.0, 588.5), (85.0, 500.0, 330.0, 515.0), (85.0, 588.5, 330.0, 612.0)],
    'F': [(336.5, 505.0, 454.0, 612.0)],
}
removed = {k: len(clear(pdf, H, v)) for k, v in CLEAR.items()}
t = Top(W, H); ver = {'master_sha256': sha(MASTER), 'removed_objects': removed}
SC = dict(S='#E69F00', R='#0072B2', E='#009E73', P='#CC79A7', M='#8C6D31')
ROUTE = {  # master legend colours / markers
    ('mouse', 'WOT'): ('#146b4a', 'o'), ('mouse', 'CIP'): ('#cc4c02', 's'), ('mouse', 'OSK'): ('#255fbc', '^'),
    ('mouse', '171127_GETM'): ('#766a8a', 'D'), ('mouse', '171127_OSKM'): ('#b3a4c4', 'v'), ('mouse', '171127_late'): ('#8c8c8c', 'P'),
    ('human', 'GSE118258'): ('#6f2c91', 'o'), ('human', 'GSE150311'): ('#5e7c8a', 's'), ('human', 'GSE221739'): ('#30363d', '^'),
    ('human', '178325_0330'): ('#73501a', 'D'), ('human', '178325_0330_perturbation'): ('#a58b61', 'v'),
    ('human', '178325_0605'): ('#a66b00', 'P'), ('human', '178325_0618'): ('#a22753', 'X'), ('human', '178325_0809'): ('#bb702f', '>'),
    ('human', '178325_1117'): ('#92816c', '<'), ('human', '178325_1230'): ('#b0a799', 'h'), ('human', '178325_reference'): ('#606060', '*')}
FIT = {('mouse', 'WOT'), ('mouse', 'CIP'), ('human', 'GSE118258')}
CL = ['C1', 'C2', 'C3']

# ------------------------------------------------------------------ data shared by A
comp = rd('F2A', 'panel_statistics/profile-species-composition.tsv')
cnt = comp.set_index(['species', 'cluster']).genes.astype(int)
TOT = comp.groupby('species').species_total.first().astype(int).to_dict()
hm = rd('F2A', 'panel_statistics/profile-centred-heatmap.tsv.gz')
hm = hm[hm.node.isin(list('REPM'))]
cm = rd('F2A', 'panel_statistics/profile-centroid-means.tsv')
cm = cm[cm.node.isin(list('SREPM'))]
means = cm.pivot_table(index=['species', 'route', 'cluster'], columns='node', values='mean_delta', aggfunc='first').reindex(columns=list('SREPM'))
ngen = cm.pivot_table(index=['species', 'route', 'cluster'], columns='node', values='n', aggfunc='first').reindex(columns=list('SREPM')).fillna(0).astype(int)
assert (means['S'].dropna().abs() < 1e-12).all(), 'S is the reference node (ΔDWUI 0)'
estimable = {(s_, r_) for (s_, r_, c_), v_ in ngen[list('REPM')].sum(axis=1).items() if v_ > 0}
none = sorted(set(ROUTE)-estimable)
for k_ in sorted(estimable-set(ROUTE)):
    R.flag('F2A', 'row_capacity', f'route {k_[0]}|{k_[1]} has S-relative estimates but no colour/marker or legend slot')
for sp_, c_ in cnt.index:
    if int(cnt[(sp_, c_)]) != hm[(hm.species == sp_) & (hm.cluster == c_)].gene_id.nunique():
        R.flag('F2A', 'consistency', f'{sp_}|{c_}: heatmap genes differ from class composition')

# ------------------------------------------------------------------ A curves
NX = dict(S=23.7, R=35.8, E=48.0, P=60.1, M=77.8)
ROWTOP = [28.2, 61.3, 94.3]; RH = 26.1
YR = {'C1': (-0.15, 0.05, [-0.1, 0.0]), 'C2': (-0.11, 0.06, [-0.1, 0.0]), 'C3': (-0.03, 0.13, [0.0, 0.1])}
Acur = {}
for sp, ox in [('mouse', 0.0), ('human', HX)]:
    for ci, c in enumerate(CL):
        lo, hi, ticks = YR[c]; top = ROWTOP[ci]
        Y = lambda v, lo=lo, hi=hi, top=top: top+(hi-v)/(hi-lo)*RH
        t.line([21.9+ox, 80.2+ox], [Y(0)]*2, lw=0.5, color='#bec6cb')
        for v in ticks:
            t.line([21.9+ox, 23.9+ox], [Y(v)]*2, lw=0.8, color=INK)
            t.text_top(20.3+ox, Y(v)-2.9, num(v, 1), 8, ha='right')
        routes = [r for (s_, r) in ROUTE if s_ == sp and (s_, r) in estimable and (sp, r, c) in ngen.index]
        routes = sorted(routes, key=lambda r: (sp, r) in FIT)   # fitted routes drawn last
        vals_ = [float(means.loc[(sp, r, c)][k]) for r in routes for k in 'SREPM' if ngen.loc[(sp, r, c)][k] > 0]
        R.axis('F2A', f'{sp} {c} route-mean ΔDWUI vs S', vals_, lo, hi)
        for r in routes:
            col, mk = ROUTE[(sp, r)]; mv = means.loc[(sp, r, c)].copy(); nv = ngen.loc[(sp, r, c)]
            mv[list('SREPM')] = mv[list('SREPM')].clip(lo+1e-9, hi-1e-9)
            pts = [(k, float(mv[k])) for k in 'SREP' if nv[k] > 0]
            for (k0, v0), (k1, v1) in zip(pts, pts[1:]):
                adj = 'SREP'.index(k1)-'SREP'.index(k0) == 1
                t.line([NX[k0]+ox, NX[k1]+ox], [Y(v0), Y(v1)], lw=0.8, color=col, ls='-' if adj else (0, (2.2, 1.4)))
            for k, v in pts + ([('M', float(mv['M']))] if nv['M'] > 0 else []):
                t.marker(NX[k]+ox, Y(v), mk, 2.6 if mk not in 'P*X' else 3.0, fc=col, ec=col, mew=0.45, zorder=4)
            Acur[f'{sp}|{r}|{c}'] = {k: (round(float(mv[k]), 5), int(nv[k])) for k in 'SREPM' if nv[k] > 0}
ver['A_curves'] = Acur

# ------------------------------------------------------------------ A heatmaps
from matplotlib.colors import LinearSegmentedColormap
CMAP = LinearSegmentedColormap.from_list('d', DIV)
COLS = {'mouse': [('R', ['WOT', 'CIP', 'OSK']), ('E', ['WOT', 'CIP', 'OSK']), ('P', ['WOT', 'CIP']),
                  ('M', ['WOT', 'CIP', 'OSK', '171127_GETM', '171127_OSKM'])],
        'human': [('R', ['GSE118258', '178325_0330', '178325_0605', '178325_0618']), ('E', ['GSE118258', 'GSE150311']),
                  ('P', ['GSE118258', 'GSE150311', 'GSE221739', '178325_0618', '178325_0809']),
                  ('M', ['GSE118258', 'GSE150311', 'GSE221739', '178325_0618'])]}
HB = {'mouse': (89.0, 182.3), 'human': (312.2, 405.5)}
LIM = 0.10
Aheat = {}
for sp in ['mouse', 'human']:
    cols = [(k, r) for k, rs in COLS[sp] for r in rs]
    # every column with any data, and no data column left out
    have = {(k, r) for (s_, r, c_), row_ in ngen.iterrows() if s_ == sp for k in 'REPM' if row_[k] > 0}
    for k_, r_ in sorted(have-set(cols)):
        R.flag('F2A', 'row_capacity', f'{sp} heatmap: {k_}:{r_} has data but no column in the locked grid')
    for k_, r_ in sorted(set(cols)-have):
        R.unavailable('F2A', f'{sp} heatmap column {k_}:{r_} has no data (drawn grey)')
    x0, x1 = HB[sp]; cw = (x1-x0)/len(cols)
    for ci, c in enumerate(CL):
        h_ = hm[(hm.species == sp) & (hm.cluster == c)]
        order_ = h_.drop_duplicates('gene_id').sort_values('class_order', kind='stable')   # one order across all route columns
        genes = order_.gene_id.tolist()
        M = np.full((len(genes), len(cols)), np.nan); gi = {g: i for i, g in enumerate(genes)}
        ci_ = {kr: j for j, kr in enumerate(cols)}
        for g, k, r, v in zip(h_.gene_id, h_.node, h_.route, h_.delta_from_S):
            if (k, r) in ci_ and v == v:
                M[gi[g], ci_[(k, r)]] = v
        rgba = CMAP((np.clip(np.nan_to_num(M), -LIM, LIM)+LIM)/(2*LIM))
        rgba[~np.isfinite(M)] = (0xDE/255, 0xDE/255, 0xDE/255, 1)
        img = np.repeat(rgba, 24, axis=1)
        top = ROWTOP[ci]
        t.ax.imshow(img, extent=(x0, x1, H-(top+RH), H-top), interpolation='none', aspect='auto', zorder=2)
        b = 0
        for k, rs in COLS[sp][:-1]:
            b += len(rs)
            t.line([x0+b*cw]*2, [top, top+RH], lw=1.5 if k == 'P' else 0.7, color='white', zorder=3)
        Aheat[f'{sp}|{c}'] = dict(genes=len(genes), columns=[f'{k}:{r}' for k, r in cols],
                                  measured_cells=int(np.isfinite(M).sum()), clipped_cells=int((np.abs(M) > LIM).sum()))
for (x, sp) in [(186.7, 'mouse'), (409.9, 'human')]:
    for ci, c in enumerate(CL):
        t.text_top(x, [43.0, 76.0, 109.0][ci], str(int(cnt[(sp, c)])), 8)
# colour bar rebuilt: the master frame and lower extend triangle sat 8–9 pt below the gradient, and the
# lowest label was centred on the triangle tip instead of the gradient end
cb0, cb1, cx0, cx1 = 141.5, 184.5, 420.8, 424.5          # same height as the route legend beside it
grad = CMAP(np.linspace(1, 0, 256))[:, None, :]
t.ax.imshow(np.repeat(grad, 4, axis=1), extent=(cx0, cx1, H-cb1, H-cb0), interpolation='bilinear', aspect='auto', zorder=2)
t.ax.add_patch(Polygon([(cx0, H-cb0), (cx1, H-cb0), ((cx0+cx1)/2, H-(cb0-3.0))], closed=True, fc=DIV[2], ec='none'))
t.ax.add_patch(Polygon([(cx0, H-cb1), (cx1, H-cb1), ((cx0+cx1)/2, H-(cb1+3.0))], closed=True, fc=DIV[0], ec='none'))
t.ax.add_patch(Polygon([(cx0, H-cb0), ((cx0+cx1)/2, H-(cb0-3.0)), (cx1, H-cb0), (cx1, H-cb1), ((cx0+cx1)/2, H-(cb1+3.0)), (cx0, H-cb1)],
                       closed=True, fc='none', ec='#87959c', lw=0.32))
for yv in (cb0, (cb0+cb1)/2, cb1):
    t.line([cx0, cx0+2.0], [yv, yv], lw=0.51, color=INK)
t.text_top(428.8, cb0-4.0, '0.10', 8); t.text_top(428.8, cb1-4.0, '−0.10', 8); t.text_top(428.8, (cb0+cb1)/2-4.0, '0.00', 8)
t.vtext(412.2, (cb0+cb1)/2, 'ΔDWUI vs S', 8)
ver['A_heatmaps'] = dict(limit=LIM, cells=Aheat)
# route legend: only routes with an S-relative estimate (the five without one are not drawn anywhere in A)
NAME = {('mouse', 'WOT'): 'Schiebinger2019–2i', ('mouse', 'CIP'): 'Guo2019–CIP', ('mouse', 'OSK'): 'Guo2019–OSK',
        ('mouse', '171127_GETM'): 'Jaber2022–GETM D3', ('mouse', '171127_OSKM'): 'Jaber2022–OSKM D3',
        ('human', 'GSE118258'): 'Xing2020', ('human', 'GSE150311'): 'Liu2020', ('human', 'GSE221739'): 'Panariello2023',
        ('human', '178325_0330'): 'Guan2022–0330', ('human', '178325_0605'): 'Guan2022–0605', ('human', '178325_0618'): 'Guan2022–0618',
        ('human', '178325_0809'): 'Guan2022–0809'}
for k_ in sorted(set(NAME)-estimable):
    R.unavailable('F2A', f'legend: {NAME[k_]} has no S-relative estimate (legend slot left empty)')
LEGR = [('Mouse', 141.6, [('mouse', 'WOT'), ('mouse', 'CIP'), ('mouse', 'OSK'), ('mouse', '171127_GETM')]),
        (None, 152.6, [('mouse', '171127_OSKM')]),
        ('Human', 165.6, [('human', 'GSE118258'), ('human', 'GSE150311'), ('human', 'GSE221739'), ('human', '178325_0330')]),
        (None, 176.6, [('human', '178325_0605'), ('human', '178325_0618'), ('human', '178325_0809')])]
COLX = [60.0, 144.0, 228.0, 312.0]
for lab_, yt_, items in LEGR:
    if lab_:
        t.text_top(4.3, yt_-0.6, lab_, 9)
    for x_, key in zip(COLX, items):
        if key not in estimable:
            continue
        col, mk = ROUTE[key]
        t.marker(x_+3.0, yt_+3.9, mk, 3.6 if mk not in 'P*X' else 4.0, fc=col, ec=col, mew=0.45)
        t.text_top(x_+8.0, yt_, NAME[key], 8)
        assert x_+8.0+textwidth(NAME[key], 8) < min([c_ for c_ in COLX[:len(items)] if c_ > x_] + [405.0])-3.0, NAME[key]
ver['A_routes_without_S_relative_estimate_not_drawn'] = [r for _, r in none]

# ------------------------------------------------------------------ B
# schematic = the k = 3 C2 shape centre, drawn level with the formula block beside it; survival panels redrawn taller
sc = rd('F2B', 'profiles/shape_centers.tsv')
c2 = sc[(sc.k == 3) & (sc.cluster == 'C2')].iloc[0]
v = {k: float(c2[k]) for k in 'SREP'}
FT = [242.4, 252.4, 262.5, 274.7]
for s_, top in zip(['L = |R−S| + |E−R| + |P−E|', 'D = |P−S|', 'T = 1 − D/L'], FT):
    t.text_top(122.9, top, s_, 7.5)
SXs = dict(S=31.7, R=55.0, E=78.3, P=101.7)
span = max(v.values())-min(v.values())
YS = lambda q: (FT[0]+3.8)+(max(v.values())-q)/span*(FT[2]+3.8-(FT[0]+3.8)+2.0)
t.line([29.5, 112.0], [YS(v['S'])]*2, lw=0.5, color='#87969E', ls=(0, (2, 1.5)))
t.line([SXs[k] for k in 'SREP'], [YS(v[k]) for k in 'SREP'], lw=1.0, color=INK)
for k in 'SREP':
    t.marker(SXs[k], YS(v[k]), 'o', 4.0, fc=SC[k], ec='white', mew=0.5, zorder=4)
    t.text_top(SXs[k], FT[3], k, 8, ha='center', color=SC[k])
ya, yb = YS(v['S'])+0.8, YS(v['P'])-0.8
t.line([109.1, 109.1], [ya, yb], lw=0.7, color=INK)
from matplotlib.patches import Polygon
for yy, d in [(ya, 1), (yb, -1)]:
    t.ax.add_patch(Polygon([(107.9, H-(yy+d*1.8)), (110.3, H-(yy+d*1.8)), (109.1, H-yy)], closed=True, fc=INK, ec='none'))
t.text_top(112.4, (ya+yb)/2-4.0, 'D', 8)
mx, my = (SXs['S']+SXs['R'])/2, (YS(v['S'])+YS(v['R']))/2
t.text_top(mx+3.0, my-8.5, 'L', 8)
t.vtext(21.0, 264.5, 'DWUI', 8)                 # below the letter slot, same x as P(T >= t)
surv = rd('F2B', 'panel_statistics/transient-survival.tsv')
Bd = {}
RT = [('mouse', 'WOT'), ('mouse', 'CIP'), ('human', 'GSE118258')]
SY0, SY1 = 299.0, 348.0                     # 100 % and 0 %
for j, c in enumerate(CL):
    dx = [0.0, 63.6, 127.1][j]
    X = lambda q, dx=dx: 46.4+dx+47.6*q; Yp = lambda p: SY1-(SY1-SY0)*p
    t.text_top(X(0.5), 287.6, c, 8, ha='center')
    t.line([X(0)-1.2]*2, [SY0, SY1+1.0], lw=0.5, color=AXC); t.line([X(0)-1.2, X(1)+1.2], [SY1+1.0]*2, lw=0.5, color=AXC)
    t.line([X(0), X(1)], [Yp(0.5)]*2, lw=0.5, color='#d5dde1')
    t.line([X(0.5)]*2, [SY0, SY1], lw=0.6, color='#aab5bc', ls=(0, (1, 1.5)))
    for q_, s_ in [(0, '0'), (0.5, '0.5'), (1, '1')]:
        t.line([X(q_)]*2, [SY1+1.0, SY1+3.0], lw=0.5, color=AXC); t.text_top(X(q_), SY1+4.2, s_, 8, ha='center')
    for p_, s_ in [(1, '100'), (0.5, '50'), (0, '0')]:
        t.line([X(0)-3.2, X(0)-1.2], [Yp(p_)]*2, lw=0.5, color=AXC)
        if j == 0: t.text_top(X(0)-4.4, Yp(p_)-4.0, s_, 8, ha='right')
    for sp, r in [('human', 'GSE118258'), ('mouse', 'WOT'), ('mouse', 'CIP')]:
        s_ = surv[(surv.species == sp) & (surv.route == r) & (surv.cluster == c)].sort_values('threshold')
        if s_.empty:
            R.unavailable('F2B', f'{sp}|{r}|{c}: no survival curve'); continue
        nT = int(s_.n.iloc[0])
        if nT != int(cnt.get((sp, c), -1)):
            R.flag('F2B', 'consistency', f'{sp}|{r}|{c}: survival n {nT} != class size {int(cnt.get((sp, c), -1))}')
        grid = s_.threshold.to_numpy(float); sv = s_.survival.to_numpy(float)
        col, mk = ROUTE[(sp, r)]
        t.s.ax.step([X(q) for q in grid], [H-Yp(p) for p in sv], where='post', lw=0.9, color=col, zorder=3)
        i5 = np.where(np.isclose(grid, 0.5))[0]; assert len(i5) == 1, 'survival grid lacks t = 0.5'
        p5 = float(sv[i5[0]])
        t.marker(X(0.5), Yp(p5), mk, 3.4, fc=col, ec='white', mew=0.3, zorder=5)
        Bd[f'{sp}|{r}|{c}'] = dict(n=nT, P_T_ge_0p5=round(p5, 4))      # median T is a 2C quantity, not a 2B mark
t.vtext(21.0, (SY0+SY1)/2, 'P(T ≥ t) (%)', 8)
t.text_top(46.4+63.6+47.6/2, SY1+14.0, 'Threshold t', 8, ha='center')
assert SY1+14.0+7.5 < 372.0
ver['B'] = dict(C2_shape_centre=v, survival=Bd)

# ------------------------------------------------------------------ C
steps = rd('F2C', 'panel_statistics/signed-stage-means.tsv').rename(columns={'signed_SR': 'sSR', 'signed_RE': 'sRE', 'signed_EP': 'sEP'})
steps = steps.set_index(['species', 'route', 'cluster'])
R.path('F2C', 'panel_statistics/pooled-class-geometry.tsv')      # registered for 2C; not drawn in this layout
YC = {'C1': [250.6, 257.5, 264.4], 'C2': [275.8, 282.7, 289.6], 'C3': [301.0, 307.9, 314.8]}
CLIM = 0.08; DLX = lambda q: 416.6+29.5*q
stq = rd('F2C', 'panel_statistics/stage-contribution-means-and-quantiles.tsv')
stg = lambda sp, r, c, m_: stq[(stq.species == sp) & (stq.route == r) & (stq.cluster == c) & (stq.measure == m_)]
t.line([DLX(0.5)]*2, [245.0, 319.6], lw=0.55, color='#bcc6cd', ls=(0, (1, 1.5)))
Cd = {}; clipped = []
for c in CL:
    for (sp, r), y in zip(RT, YC[c]):
        col, mk = ROUTE[(sp, r)]
        t.marker(271.6, y, mk, 3.6, fc=col, ec=col, mew=0.3)
        if (sp, r, c) not in steps.index or not len(stg(sp, r, c, 'contribution_SR')):
            R.unavailable('F2C', f'{sp}|{r}|{c}: no stage summary (row left empty)'); continue
        st = steps.loc[(sp, r, c)]
        for k, key in enumerate(['sSR', 'sRE', 'sEP']):
            val = float(st[key]); x0 = 286.3+21.6*k
            if abs(val) > CLIM:
                clipped.append((sp, r, c, key, round(val, 4)))
            t.rect(x0, y-3.15, 21.6, 6.3, fc=diverge(val, CLIM), ec='white', lw=0.6)
            ac = 'white' if abs(val)/CLIM > 0.6 else INK; xc = x0+10.8; d = 1 if val > 0 else -1
            t.line([xc, xc], [y+d*1.5, y-d*0.6], lw=0.65, color=ac)
            t.ax.add_patch(Polygon([(xc-1.0, H-(y-d*0.4)), (xc+1.0, H-(y-d*0.4)), (xc, H-(y-d*1.9))], closed=True, fc=ac, ec='none'))
        fr = [float(stg(sp, r, c, m_).iloc[0]['mean']) for m_ in ('contribution_SR', 'contribution_RE', 'contribution_EP')]
        if abs(sum(fr)-1) > 1e-6:
            R.flag('F2C', 'consistency', f'{sp}|{r}|{c}: mean stage contributions sum to {sum(fr):.4f}')
        xx = 359.9
        for f_, cc in zip(fr, [SC['R'], SC['E'], SC['P']]):
            t.rect(xx, y-2.95, 50.1*f_, 5.9, fc=cc, ec='white', lw=0.3); xx += 50.1*f_
        dlr = stg(sp, r, c, 'endpoint_retention_D_over_L')
        tq = stg(sp, r, c, 'transient_fraction_T')
        if len(dlr):
            med, q1, q3 = (float(dlr.iloc[0][k]) for k in ['q50', 'q25', 'q75'])
        elif len(tq):                                   # polyApipe bundle: stage table carries T but not D/L; T = 1 - D/L exactly
            med, q1, q3 = 1-float(tq.iloc[0]['q50']), 1-float(tq.iloc[0]['q75']), 1-float(tq.iloc[0]['q25'])
            if not getattr(R, '_dl_from_T', False):
                R.derive('F2C', 'D/L median and quartiles taken as 1 - T quartiles (T = 1 - D/L); the stage table of this bundle lacks the endpoint_retention_D_over_L rows')
                R._dl_from_T = True
        else:
            R.unavailable('F2C', f'{sp}|{r}|{c}: endpoint_retention_D_over_L quantiles not exported'); med = q1 = q3 = float('nan')
        t.line([DLX(0), DLX(1)], [y, y], lw=0.55, color='#e9edf0')
        t.line([DLX(q1), DLX(q3)], [y, y], lw=1.1, color=col)
        t.marker(DLX(med), y, mk, 2.8, fc=col, ec=col, mew=0.3, zorder=5)
        Cd[f'{sp}|{r}|{c}'] = dict(step_mean=[round(float(st[k]), 4) for k in ['sSR', 'sRE', 'sEP']],
                                   contribution_mean=[round(q, 4) for q in fr], D_over_L_median=med, D_over_L_IQR=[q1, q3], n=int(st['n']))
t.text_top(449.4, 235.2, 'Median D/L', 8, ha='right')
for q_ in (0, 1): t.text_top(DLX(q_), 320.8, str(q_), 8, ha='center')
if any(v_[4] < -CLIM for v_ in clipped):
    t.ax.add_patch(Polygon([(377.6, H-344.9), (377.6, H-348.0), (374.6, H-346.45)], closed=True, fc=DIV[0], ec='none'))
if any(v_[4] > CLIM for v_ in clipped):
    t.ax.add_patch(Polygon([(436.6, H-344.9), (436.6, H-348.0), (439.6, H-346.45)], closed=True, fc=DIV[2], ec='none'))
ver['C'] = dict(rows=Cd, colour_limit=CLIM, clipped=clipped)

# ------------------------------------------------------------------ D
cells = rd('F2D', 'frozen/frozen_cell_embedding.tsv.gz')
scores = rd('F2D', 'moran/k3_group_wpp_metacells.tsv.gz')
alias = {'WOT_GSE115943': 'GSE115943'}
cells['_s'] = cells.study_id.replace(alias); scores['_s'] = scores.study_id.replace(alias)
j_ = cells.merge(scores.drop(columns=['study_id', 'species']), on=['metacell_id', '_s'], how='left', validate='many_to_one', indicator=True)
assert len(j_) == len(cells)
if not j_._merge.eq('both').all():
    R.flag('F2D', 'unavailable', f'{int((~j_._merge.eq("both")).sum())} cells without a k3 metacell score (drawn grey)')
mor = rd('F2D', 'moran/moran_results.tsv')
mb = mor[(mor.k == 20) & (mor.test == 'state_coverage_library_state') & (mor.status == 'tested')]
BOX = {'mouse': [(5.1, 68.3), (79.6, 142.7), (154.1, 217.2)], 'human': [(242.4, 300.3), (316.9, 374.8), (391.3, 449.3)]}
YB = {'mouse': (399.6, 462.6), 'human': (398.1, 462.2)}
LBL = {'mouse': [13.3, 78.8, 153.3], 'human': [232.0, 306.6, 381.1]}
STUDY_ORDER = {'mouse': ['GSE103221', 'GSE171127', 'GSE115943'], 'human': ['GSE118258', 'GSE150311', 'GSE178325', 'GSE221739']}
STUDY_COL = {'GSE103221': '#cc4c02', 'GSE171127': '#766a8a', 'GSE115943': '#146b4a', 'GSE118258': '#6f2c91',
             'GSE150311': '#5e7c8a', 'GSE178325': '#a22753', 'GSE221739': '#30363d'}
LET = dict(shared_s='S', shared_r='R', shared_e='E', shared_tp='P', shared_m='M')
Dd = {}
for sp in ['mouse', 'human']:
    q = j_[j_.species_key == sp]
    ux0, ux1 = q.UMAP1.min(), q.UMAP1.max(); uy0, uy1 = q.UMAP2.min(), q.UMAP2.max()
    y0, y1 = YB[sp]
    for ci, c in enumerate(CL):
        x0, x1 = BOX[sp][ci]
        kx = (x1-x0)/(ux1-ux0); ky = (y1-y0)/(uy1-uy0)
        assert abs(kx/ky-1) < 0.02, (sp, kx, ky)
        k_ = min(kx, ky)
        PX = lambda u: x0+(u-ux0)*k_ + ((x1-x0)-(ux1-ux0)*k_)/2
        PY = lambda u: H-(y1-(u-uy0)*k_ - ((y1-y0)-(uy1-uy0)*k_)/2)
        val = q[c+'_centered'].to_numpy(float)
        t.ax.scatter(PX(q.UMAP1.to_numpy()), PY(q.UMAP2.to_numpy()), s=0.035, c='#DADFE3', linewidths=0, rasterized=True, zorder=2)
        ok = np.isfinite(val); o = np.argsort(np.abs(val[ok]))
        xs = PX(q.UMAP1.to_numpy()[ok][o]); ys = PY(q.UMAP2.to_numpy()[ok][o])
        t.ax.scatter(xs, ys, s=0.035, c=val[ok][o], cmap=CMAP, vmin=-0.10, vmax=0.10, linewidths=0, rasterized=True, zorder=3)
        for grp, L in LET.items():
            m_ = q.state_group == grp
            lx, ly = PX(float(q.UMAP1[m_].median())), PY(float(q.UMAP2[m_].median()))
            t.s.ax.text(lx, ly, L, fontproperties=FONTS['Bol'], fontsize=8, color=INK, ha='center', va='center', zorder=6,
                        path_effects=[pe.withStroke(linewidth=1.6, foreground='white')])
        studies = mb[(mb.species == sp) & (mb.cluster == c)]
        k_sig = int((studies.q_BH < 0.05).sum()); n_st = len(studies)
        t.text_top(LBL[sp][ci], 388.6, c, 8)
        xsq = LBL[sp][ci]+textwidth(c, 8)+3.5
        for sid in STUDY_ORDER[sp]:
            rq = studies[studies.study_id == sid]; col_ = STUDY_COL[sid]
            if len(rq) != 1:
                R.unavailable('F2D', f'{sp}|{c}: no Moran test for {sid} (grey square)')
                t.rect(xsq, 390.1, 4.6, 4.6, fc=MISSING, ec=MISSING, lw=0.8, zorder=7); xsq += 6.2; continue
            sig = float(rq.q_BH.iloc[0]) < 0.05
            t.rect(xsq, 390.1, 4.6, 4.6, fc=col_ if sig else 'white', ec=col_, lw=0.8, zorder=7); xsq += 6.2
        for sid in sorted(set(studies.study_id)-set(STUDY_ORDER[sp])):
            R.flag('F2D', 'row_capacity', f'{sp}|{c}: tested study {sid} has no square in the layout')
        if ok.sum() and (np.abs(val[ok]) > 0.10).mean() > 0.05:
            R.flag('F2D', 'colour_saturation', f'{sp}|{c}: {100*(np.abs(val[ok]) > 0.10).mean():.1f}% of scored cells beyond ±0.10')
        Dd[f'{sp}|{c}'] = dict(cells=int(len(q)), cells_with_score=int(ok.sum()), clipped=int((np.abs(val[ok]) > 0.10).sum()),
                               moran_q_lt_0p05=f'{k_sig}/{n_st}', studies=studies.study_id.tolist())
lab_ = 'Gene-centered DWUI, within study'; yr = 470.0
t.text_top(6.0, yr, lab_, 8); x = 6.0+textwidth(lab_, 8)+6.0
t.text_top(x, yr, '−0.10', 8); x += textwidth('−0.10', 8)+5.5
cbx0, cbx1 = x, x+64.0
grad = CMAP(np.linspace(0, 1, 256))[None, :, :]
t.ax.imshow(np.repeat(grad, 3, axis=0), extent=(cbx0, cbx1, H-(yr+6.3), H-(yr+1.3)), interpolation='bilinear', aspect='auto', zorder=2)
t.ax.add_patch(Polygon([(cbx0, H-(yr+1.3)), (cbx0, H-(yr+6.3)), (cbx0-3.6, H-(yr+3.8))], closed=True, fc=DIV[0], ec='none'))
t.ax.add_patch(Polygon([(cbx1, H-(yr+1.3)), (cbx1, H-(yr+6.3)), (cbx1+3.6, H-(yr+3.8))], closed=True, fc=DIV[2], ec='none'))
t.line([(cbx0+cbx1)/2]*2, [yr+0.6, yr+7.0], lw=0.5, color=INK)
x = cbx1+5.5; t.text_top(x, yr, '0.10', 8); x += textwidth('0.10', 8)+18.0
t.text_top(x, yr, 'Local structure per study', 8); x += textwidth('Local structure per study', 8)+6.0
t.rect(x, yr+1.4, 4.6, 4.6, fc=INK, ec=INK, lw=0.8); t.text_top(x+7.0, yr, 'q < .05', 8); x += 7.0+textwidth('q < .05', 8)+6.0
t.rect(x, yr+1.4, 4.6, 4.6, fc='white', ec=INK, lw=0.8); t.text_top(x+7.0, yr, 'q ≥ .05', 8); x += 7.0+textwidth('q ≥ .05', 8)
assert x < 451.0, x
SNAME = {'GSE103221': 'Guo2019', 'GSE171127': 'Jaber2022', 'GSE115943': 'Schiebinger2019', 'GSE118258': 'Xing2020',
         'GSE150311': 'Liu2020', 'GSE178325': 'Guan2022', 'GSE221739': 'Panariello2023'}
yl = 481.4; x = 6.0
for sp in ['mouse', 'human']:
    for sid in STUDY_ORDER[sp]:
        t.rect(x, yl+1.4, 4.6, 4.6, fc=STUDY_COL[sid], ec=STUDY_COL[sid], lw=0.8); t.text_top(x+7.0, yl, SNAME[sid], 8)
        x += 7.0+textwidth(SNAME[sid], 8)+9.0
    x += 9.0
assert x < 460.0, x
ver['D'] = Dd

# ------------------------------------------------------------------ E
# rows = route as a colour block (legend below); two 3 x 3 matrices of class-mean profile r (RMSE in verification.json):
# held-out route class means vs centroids refitted on the remaining routes, and each route's class means vs the full-fit centroids (no refit); right: assignment agreement and ARI
lc = rd('F2E', 'profiles/leave_route_centroids.tsv'); lp = rd('F2E', 'profiles/leave_route_projection.tsv')
ra = rd('F2E', 'profiles/route_shape_support.tsv')
ER = [('mouse', 'WOT', 'Schiebinger2019–2i'), ('mouse', 'CIP', 'Guo2019–CIP'), ('human', 'GSE118258', 'Xing2020')]
from matplotlib.colors import LinearSegmentedColormap as _LSC
SEQ = _LSC.from_list('seq', ['#EEF3F7', '#9CC3DD', '#2D6A96', '#163A57']); R0, R1 = 0.75, 1.0
seqc = lambda r_: SEQ((min(max(r_, R0), R1)-R0)/(R1-R0))
cw, rh, ty0 = 27.0, 15.0, 517.0; MX = [26.0, 110.0]
ty1 = ty0+3*rh
for m_, head in zip(MX, ['Route held out', 'Route vs full fit']):
    t.text_top(m_+1.5*cw, 496.6, head, 8, ha='center')
    for k, c in enumerate(CL):
        t.text_top(m_+cw*(k+0.5), 506.7, c, 8, ha='center')
Ed = {}
for i, (sp, r, name) in enumerate(ER):
    yc = ty0+rh*(i+0.5); col, mk = ROUTE[(sp, r)]
    t.rect(14.0, yc-3.0, 6.0, 6.0, fc=col)
    for m_, tab, key in [(MX[0], lc, 'heldout'), (MX[1], ra, 'alone')]:
        for k, c in enumerate(CL):
            sel = (tab.heldout_route == r) if key == 'heldout' else (tab.route == r)
            row = tab[sel & (tab.cluster == c)]
            if len(row) != 1:
                R.unavailable('F2E', f'{r}|{c}|{key}: no centroid row (cell left empty)'); continue
            row = row.iloc[0]
            rr, rm = float(row.centroid_correlation), float(row.centroid_RMSE)
            if key == 'heldout':
                R.claim('F2E', f'{r} {c} shape recovered when the route is held out', bool(row.shape_recovered) and row.get('nearest_LOO_cluster', c) == c)
            else:
                R.claim('F2E', f'{r} {c} shape supported against the full fit', bool(row.shape_supported) and row.get('nearest_full_cluster', c) == c)
            R.axis('F2E', f'{r}|{c}|{key} centroid r (colour)', [rr], R0, R1)
            t.rect(m_+cw*k, ty0+rh*i, cw, rh, fc=seqc(rr), ec='white', lw=0.8)
            tc = 'white' if (rr-R0)/(R1-R0) > 0.45 else INK
            t.text_top(m_+cw*(k+0.5), yc-4.0, f'{rr:.3f}', 8, ha='center', color=tc)     # r only; RMSE in verification.json
            Ed[f'{r}|{c}|{key}'] = dict(r=rr, RMSE=rm, n=int(row.n))
# colour scale, inline
t.text_top(MX[0], ty1+4.6, 'r', 8); x = MX[0]+textwidth('r', 8)+5.0
t.text_top(x, ty1+4.6, '0.75', 7.5); x += textwidth('0.75', 7.5)+3.0
for i_ in range(60):
    v_ = R0+(R1-R0)*(i_+0.5)/60
    t.rect(x+i_*70.0/60, ty1+5.6, 70.0/60+0.1, 4.0, fc=seqc(v_))
x += 73.0; t.text_top(x, ty1+4.6, '1', 7.5)
# assignment agreement vs full / remaining-route fit, ARI columns
AX0, AX1, LO, HI = 205.0, 270.0, 0.6, 1.0
AXx = lambda q: AX0+(q-LO)/(HI-LO)*(AX1-AX0)
t.text_top((AX0+AX1)/2, 496.6, 'Assignment agreement', 8, ha='center')
t.text_top(304.5, 496.6, 'ARI', 8, ha='center')
t.marker(290.0, 510.6, 'o', 3.6, fc=INK, ec=INK, mew=0.6); t.marker(310.5, 510.6, 'o', 3.6, fc='white', ec=INK, mew=0.7)
for i, (sp, r, name) in enumerate(ER):
    yc = ty0+rh*(i+0.5); col, mk = ROUTE[(sp, r)]
    t.rect(AX0-4, yc-6.0, AX1-AX0+8, 12.0, fc=BAND)
    assert AX0-4 > MX[1]+3*cw+4
    row = lp[lp.heldout_route == r]
    if len(row) != 1:
        R.unavailable('F2E', f'{r}: no projection row'); continue
    row = row.iloc[0]
    af, arf = float(row.assignment_agreement_vs_full), float(row.ARI_vs_full)
    t.marker(AXx(af), yc, 'o', 4.0, fc=col, ec=col, mew=0.6, zorder=5)
    t.text_top(297.0, yc-4.0, f'{arf:.2f}', 8, ha='right')
    ar = row.assignment_agreement_vs_remaining_route
    if pd.notna(ar):
        t.marker(AXx(float(ar)), yc, 'o', 4.0, fc='white', ec=col, mew=0.8, zorder=5)
        t.text_top(317.5, yc-4.0, f'{float(row.ARI_vs_remaining_route):.2f}', 8, ha='right')
    else:
        t.text_top(317.5, yc-4.0, 'n/a', 8, ha='right', color=MUTED)
    R.axis('F2E', f'{r} assignment agreement', [af]+([float(ar)] if pd.notna(ar) else []), LO, HI)
    Ed[f'{r}|projection'] = dict(n=int(row.n_heldout), agreement_full=af, ARI_full=arf,
                                 agreement_remaining=None if pd.isna(ar) else float(ar),
                                 ARI_remaining=None if pd.isna(ar) else float(row.ARI_vs_remaining_route),
                                 training=row.get('training_species'))
t.line([AX0, AX1], [ty1+1.0]*2, lw=0.5, color=AXC)
for q in (0.6, 0.8, 1.0):
    t.line([AXx(q)]*2, [ty1+1.0, ty1+3.0], lw=0.5, color=AXC)
    t.text_top(AXx(q), ty1+4.0, f'{q:.1f}' if q < 1 else '1', 8, ha='center')
# legend row: route colour blocks and fit markers
yl = ty1+17.0; x = 6.0
for sp, r, name in ER:
    col, mk = ROUTE[(sp, r)]
    t.rect(x, yl+1.4, 4.6, 4.6, fc=col); t.text_top(x+7.0, yl, name, 8); x += 7.0+textwidth(name, 8)+9.0
assert x < 318.0, x
x = 6.0; yl += 10.5
t.marker(x+2.0, yl+3.7, 'o', 3.6, fc=INK, ec=INK, mew=0.6); t.text_top(x+6.5, yl, 'full fit', 8); x += 6.5+textwidth('full fit', 8)+8.0
t.marker(x+2.0, yl+3.7, 'o', 3.6, fc='white', ec=INK, mew=0.7); t.text_top(x+6.5, yl, 'remaining route', 8); x += 6.5+textwidth('remaining route', 8)
assert x < 318.0, x
ver['E'] = Ed

# ------------------------------------------------------------------ F
# left: class shares per species; right: class correspondence of the 82 one-to-one orthologs selected in both species
bsh = rd('F2F', 'shape_ortholog/shape_contingency_3x3.tsv')
cont = bsh[bsh.model == 'selected_overlap_margins']
assert len(cont) == 9, 'contingency must be 3 x 3'
N_ORTH = int(cont.observed.sum())
R.path('F2F', 'shape_ortholog/shape_tests.tsv')      # tests belong to the caption; registered here
compF = rd('F2F', 'panel_statistics/profile-species-composition.tsv').set_index(['species', 'cluster'])
TOTF = compF.groupby(level='species').species_total.first().astype(int).to_dict()
MC, HC = '#146b4a', '#6f2c91'
_x = 350.5                                       # clear of the panel-letter slot
_lm, _lh = f'Mouse {TOTF["mouse"]}', f'Human {TOTF["human"]}'
t.marker(_x+2.0, 510.6, 'o', 4.2, fc=MC, ec=MC, mew=0.5); t.text_top(_x+6.5, 506.6, _lm, 8); _x += 6.5+textwidth(_lm, 8)+6.0
t.marker(_x+2.0, 510.6, 'D', 4.0, fc=HC, ec=HC, mew=0.5); t.text_top(_x+6.5, 506.6, _lh, 8); _x += 6.5+textwidth(_lh, 8)
assert _x < 451.5, _x
FX = lambda p: 352.0+p*48.0/70  # polyApipe redraw 2026-10-06: axis rescaled within the same span
YF = [541.0, 559.0, 577.0]
Fd = {}
for i, c in enumerate(CL):
    y = YF[i]
    t.rect(350.5, y-8.0, 52.0, 16.0, fc=BAND)
    t.text_top(348.5, y-4.0, c, 8, ha='right')
    pm = 100*float(compF.loc[('mouse', c)].fraction); ph = 100*float(compF.loc[('human', c)].fraction)
    if not R.axis('F2F', f'{c} class share (%)', [pm, ph], 0, 70):
        pm, ph = min(pm, 60.0), min(ph, 60.0)
    t.line([FX(pm), FX(ph)], [y, y], lw=1.6, color='#C9CFD3')
    t.marker(FX(pm), y, 'o', 4.2, fc=MC, ec=MC, mew=0.5, zorder=5)
    t.marker(FX(ph), y, 'D', 4.0, fc=HC, ec=HC, mew=0.5, zorder=5)
    Fd[c] = dict(mouse=[int(compF.loc[('mouse', c)].genes), TOTF['mouse'], round(pm, 2)], human=[int(compF.loc[('human', c)].genes), TOTF['human'], round(ph, 2)])
t.line([FX(0), FX(70)], [587.0]*2, lw=0.5, color=AXC)
for p in (0, 30, 60):
    t.line([FX(p)]*2, [587.0, 589.0], lw=0.5, color=AXC); t.text_top(FX(p), 589.6, str(p), 8, ha='center')
t.text_top(FX(30), 600.6, 'Genes in class (%)', 8, ha='center')
# 3 x 3 correspondence, rows human class (shared row labels), columns mouse class; shading observed / expected
mx0, mcw = 409.0, 14.3
t.text_top(mx0+1.5*mcw, 515.6, 'Mouse class', 7.5, ha='center', color=MUTED)
for k, c in enumerate(CL):
    t.text_top(mx0+(k+0.5)*mcw, 524.2, c, 7.5, ha='center')
t.vtext(404.6, YF[1], 'Human', 7.5, color=MUTED)
Fm = {}
for i, hc in enumerate(CL):
    for k, mc in enumerate(CL):
        r_ = cont[(cont.human_cluster == hc) & (cont.mouse_cluster == mc)].iloc[0]
        o, e = int(r_.observed), float(r_.expected); lr2 = float(np.log2(o/e)) if o > 0 else -1.0
        x, y = mx0+k*mcw, YF[i]-8.0
        t.rect(x, y, mcw-0.8, 16.0, fc=diverge(lr2, 1.0), ec=INK if hc == mc else 'none', lw=0.8 if hc == mc else 0)
        t.text_top(x+(mcw-0.8)/2, y+4.2, str(o), 8, ha='center', color='white' if abs(lr2) > 0.7 else INK)
        Fm[f'{hc}|{mc}'] = dict(observed=o, expected=round(e, 3), log2_obs_over_exp=round(lr2, 3))
lab = f'{N_ORTH} orthologs'
if textwidth(lab, 7.5) >= 3*mcw+1.0:              # polyApipe redraw: 3-digit counts no longer fit under the 3x3 block
    lab = f'n = {N_ORTH}'
    R.flag('F2F', 'label_changed', f'"{N_ORTH} orthologs" wider than the 3x3 block; printed as "{lab}" (caption must name the orthologs)')
assert textwidth(lab, 7.5) < 3*mcw+1.0
t.text_top(mx0+1.5*mcw-0.4, 589.6, lab, 7.5, ha='center', color=MUTED)
ver['F'] = dict(shares=Fd, ortholog_correspondence=Fm,
                diagonal=int(sum(Fm[f'{c}|{c}']['observed'] for c in CL)), diagonal_expected=round(sum(Fm[f'{c}|{c}']['expected'] for c in CL), 3))

t.save(R.out/'overlay.pdf'); overlay(pdf, R.out/'overlay.pdf', 0, 0); pdf.save(OUT)
render_png(OUT, R.out/'Figure.2.draft.png', 150)
ver['output_sha256'] = sha(OUT)
R.finish(2, OUT, json.loads(json.dumps(ver, default=str)))
