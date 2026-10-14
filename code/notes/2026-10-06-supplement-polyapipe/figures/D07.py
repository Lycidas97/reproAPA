"""Supplementary Figure D07 (printed number from supdraw.SNUM) — class identity, stage structure and endpoint retention (R2.3–R2.5), polyApipe v2 primary.

A  ortholog pairs selected in both species: same-class null (two models) and the nine human×mouse class cells
B  class composition of the full candidate sets and of the 139 paired orthologs
C  per class × fitting route: median/IQR share of path length L per step; per-gene endpoint retention D/L with median/IQR
D  sensitivity of the route median T to a minimum path length L (descriptive)
E  mouse dominant-stage agreement between Guo2019–CIP and Schiebinger2019–2i against the marginal expectation
Layout follows STYLE.zh.md (legends at the panel bottom, aligned across a row).
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

D = Data(); OUT = outdir('D07')
X0, R = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 110.5, 146.8, 95.5
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
WL = 292.0; XB = X0+WL+GAP; WB = R-XB
F = Fig(H, figname('D07'))
CL = ['C1', 'C2', 'C3']
CLS = {'C1': None, 'C2': 'dash', 'C3': 'dot'}
ROUTES = FIT_ROUTES
LEG1 = Y1+102.5                      # legend line of row 1 (A and B)

# ================================================================ A  ortholog class correspondence
x, y, w, h = F.panel('A', X0, Y1, WL, H1)
pmf = D.df('D07B', 'supplement_statistics/shape-null-PMF.tsv')
tests = D.df('D07A', 'shape_ortholog/shape_tests.tsv')
sx = Lin(18, 62, x+32, x+150); sy = Lin(0, 0.10, y+76, y+6)
F.lineplot(sx, sy, 'null PMF', standalone=False)
F.yaxis(sy, x+30, [0, 0.05, 0.10], ['0', '.05', '.10'], 'Null probability', title_y=y+44.5)
F.xaxis(sx, y+78, [20, 30, 40, 50, 60], title='Same-class pairs')
Av = {}
for model, col, ls in [('selected_overlap_margins', INK, '-'), ('measurement_q2', MUTED, 'dash')]:
    m = pmf[pmf.model == model].sort_values('same_class_genes')
    F.line([sx(v) for v in m.same_class_genes], [sy(v) for v in m.probability], lw='data', color=col, ls=ls)
    dg = tests[(tests.model == model) & (tests.test == 'diagonal')].iloc[0]
    Av[model] = dict(observed=int(dg.observed), expected=float(dg.expected), q20=float(dg.q_exploratory_20), pmf_mass=float(m.probability.sum()))
obs = Av['selected_overlap_margins']['observed']
F.line([sx(obs)]*2, [sy(0), sy(0.10)], lw='data', color=INK)
F.text(sx(obs)+1.8, sy(0.10)-0.3, f'{obs} observed', 7.0)
F.text(sx(obs)+1.8, sy(0.10)+7.4, f'q = {qfmt(Av["selected_overlap_margins"]["q20"])}', 7.0, color=MUTED)
# nine class cells + diagonal (composition null): observed / null mean with the central 95% null range
cells = ['diagonal']+[f'{a}_{b}' for a in CL for b in CL]
fx = Lin(math.log(0.25), math.log(2.5), x+186, x+276)
F.text(x+157, y+0.4, 'Human–mouse class', 7.0, color=MUTED)
F.text(x+w, y+0.4, 'q', 7.0, color=MUTED, ha='right')
rows_y = [y+13.5]+[y+22.8+6.6*i for i in range(9)]
F.line([fx(0.0)]*2, [rows_y[0]-3.6, rows_y[-1]+3.6], lw='ref', color=AXC, ls='ref')
cm = tests[tests.model == 'selected_overlap_margins'].set_index('test')
Av['cells'] = {}
for c, yy in zip(cells, rows_y):
    r_ = cm.loc[c]; lo, hi = r_.null_p025/r_.null_mean, r_.null_p975/r_.null_mean; ratio = r_.observed/r_.null_mean
    F.text(x+157, yy-3.2, 'Diagonal' if c == 'diagonal' else c.replace('_', '–'), 6.8)
    F.line([fx.clip(math.log(lo)), fx.clip(math.log(hi))], [yy, yy], lw='range', color=NULLBAR)
    sig = r_.q_exploratory_20 < 0.05
    F.marker(fx.clip(math.log(ratio)), yy, 'D' if c == 'diagonal' else 'o', 'small', fc=INK if sig else 'white', ec=INK)
    F.text(x+w, yy-3.2, qfmt(r_.q_exploratory_20), 6.8, ha='right')
    Av['cells'][c] = dict(observed=int(r_.observed), null_mean=float(r_.null_mean), ratio=float(ratio), null_ratio_95=[float(lo), float(hi)],
                          q20=float(r_.q_exploratory_20))
F.xaxis(fx, y+78, [math.log(v) for v in (0.25, 0.5, 1, 2)], ['.25', '.5', '1', '2'], 'Observed / null mean')
F.key_row(None, LEG1, [('line', (INK, None, 'data'), 'Composition null'), ('line', (MUTED, 'dash', 'data'), 'Measurement null')],
          gap=7.0, center=(sx.p0+sx.p1)/2)
F.key_row(None, LEG1, [('line', (NULLBAR, None, 'range'), '95% null range'), ('marker', ('o', 'white', INK), 'q ≥ .05')],
          gap=7.0, center=(fx.p0+fx.p1)/2)
F.values['A'] = Av

# ================================================================ B  class composition (form of main Figure 2f, left)
x, y, w, h = F.panel('B', XB, Y1, WB, H1)
mem = D.df('D07C', 'profiles/all_resolution_memberships.tsv.gz')
ct = D.df('D07C', 'shape_ortholog/shape_contingency_3x3.tsv'); ct = ct[ct.model == 'selected_overlap_margins']
kcol = [c for c in mem.columns if c in ('k', 'resolution')][0]
m3 = mem[mem[kcol] == 3]
bx_ = Lin(0, 70, x+20.0, x+w-5.0)
F.text((bx_.p0+bx_.p1)/2, y+3.0, 'Class composition', 7.5, ha='center')
Bv, N_ = {}, {}
for sp in ('human', 'mouse'):
    for kind in ('all', 'pairs'):
        if kind == 'all':
            cnt = m3[m3.species == sp].groupby('cluster').gene_id.nunique().reindex(CL).fillna(0).astype(int)
        else:
            cnt = ct.groupby('human_cluster' if sp == 'human' else 'mouse_cluster')['observed'].sum().reindex(CL).fillna(0).astype(int)
        Bv[f'{sp}|{kind}'] = {c: int(cnt[c]) for c in CL}; N_[(sp, kind)] = int(cnt.sum())
MK = {'mouse': 'o', 'human': 'D'}
for i, c in enumerate(CL):
    yy = y+24.0+i*15.0
    F.rect(bx_.p0-2.0, yy-6.8, bx_.p1-bx_.p0+4.0, 13.6, fc=BAND)
    F.text(bx_.p0-4.0, yy-3.4, c, 7.5, ha='right')
    for kind, dy, lw_ in (('all', -2.6, 'bold'), ('pairs', 2.6, 'thin')):
        pc = {sp: 100*Bv[f'{sp}|{kind}'][c]/N_[(sp, kind)] for sp in ('human', 'mouse')}
        F.line([bx_(pc['mouse']), bx_(pc['human'])], [yy+dy]*2, lw=lw_, color=NULLBAR)
        for sp in ('mouse', 'human'):
            F.marker(bx_(pc[sp]), yy+dy, MK[sp], 'mark' if kind == 'all' else 'small',
                     fc=SPECIES[sp] if kind == 'all' else 'white', ec=SPECIES[sp])
F.xaxis(bx_, y+24.0+2*15.0+8.0, [0, 30, 60], title='Genes in class (%)')
F.key_row(None, LEG1-10.0, [('marker', ('o', SPECIES['mouse']), f"Mouse {N_[('mouse', 'all')]}"),
                            ('marker', ('D', SPECIES['human']), f"Human {N_[('human', 'all')]}")], gap=6.0, center=x+w/2)
F.key_row(None, LEG1, [('marker', ('o', 'white', SPECIES['mouse']), ''), ('marker', ('D', 'white', SPECIES['human']), f"{N_[('human', 'pairs')]} ortholog pairs")],
          gap=1.0, center=x+w/2)
F.values['B'] = Bv
F.notes.append('B: form of main Figure 2f (left): filled = all genes of each species, open = the ortholog pairs selected in both species')

# ================================================================ C  stage structure and endpoint retention
x, y, w, h = F.panel('C', X0, Y2, R-X0, H2)
rs = D.df('D07D', 'profiles/architecture_route_summary.tsv').set_index(['route', 'cluster'])
geo = D.df('D07D', 'profiles/trajectory_geometry.tsv.gz')
STEPS = [('SR', 'S→R', SC['R']), ('RE', 'R→E', SC['E']), ('EP', 'E→P', SC['P'])]
cx = [Lin(0, 1, x+50+70*i, x+114+70*i) for i in range(3)]
dx = Lin(0, 1, x+262, x+w-4.6)
for (k, lab, col), s_ in zip(STEPS, cx):
    F.text((s_.p0+s_.p1)/2, y+3.0, lab, 7.5, ha='center', color=col)
F.text(x+44, y+3.0, 'n', 7.0, color=MUTED, ha='right')
# 2026-10-12: 'Genes, median and IQR' header moved to the caption
pitch, ggap = 10.22, 4.0
rows = []; yy = y+21.2          # first band starts below the letter clear zone
for c in CL:
    for r in ROUTES:
        rows.append((c, r, yy)); yy += pitch
    yy += ggap
ylast = rows[-1][2]
for ci, c in enumerate(CL):
    ys_ = [r[2] for r in rows if r[0] == c]
    if ci % 2 == 0:
        F.rect(x, ys_[0]-5.3, w, ys_[-1]-ys_[0]+10.6, fc=BAND)
    F.text(x+1.0, ys_[0]-3.8, c, 8.0, font='Bol')
for s_ in cx:
    F.line([s_(1/3)]*2, [rows[0][2]-5.3, ylast+5.3], lw='ref', color=AXC, ls='ref')
F.line([dx(0.5)]*2, [rows[0][2]-5.3, ylast+5.3], lw='ref', color=AXC, ls='ref')
rng = np.random.default_rng(20261006)
Cv = {}
for c, r, yy in rows:
    nm, col, mk = STUDY[r]
    F.marker(x+22, yy, mk, 'mark', fc=col)
    q = rs.loc[(r, c)]
    F.text(x+44, yy-3.4, f'{int(q.n)}', 7.0, ha='right', color=MUTED)
    for (k, lab, sc_), s_ in zip(STEPS, cx):
        F.line([s_(q[f'contribution_{k}_q25']), s_(q[f'contribution_{k}_q75'])], [yy, yy], lw='bold', color=sc_)
        F.marker(s_(q[f'contribution_{k}_median']), yy, 'o', 'small', fc='white', ec=sc_)
    g = geo[(geo.route == r) & (geo.cluster == c)]
    dl = 1.0-g['T'].to_numpy()
    F.points([dx(v) for v in dl], yy+rng.uniform(-3.5, 3.5, len(dl)), color=col, alpha='dense')
    lo, me, hi = q.endpoint_retention_D_over_L_q25, q.endpoint_retention_D_over_L_median, q.endpoint_retention_D_over_L_q75
    F.line([dx(lo), dx(hi)], [yy, yy], lw='bold', color=INK)
    F.line([dx(me)]*2, [yy-3.6, yy+3.6], lw='bold', color=INK)
    Cv[f'{r}|{c}'] = dict(n=int(q.n), genes_drawn=int(len(g)), contribution={k: [float(q[f'contribution_{k}_q25']), float(q[f'contribution_{k}_median']),
                          float(q[f'contribution_{k}_q75'])] for k, _, _ in STEPS}, D_over_L=[float(lo), float(me), float(hi)])
    if len(g) != int(q.n): F.notes.append(f'C {r}|{c}: {len(g)} geometry rows vs n={int(q.n)} in the route summary')
ax_y = ylast+5.8
for s_ in cx:
    F.xaxis(s_, ax_y, [0, 0.5, 1], ['0', '.5', ''])
F.xaxis(dx, ax_y, [0, 0.25, 0.5, 0.75, 1], ['0', '.25', '.5', '.75', '1'])
F.text((cx[0].p0+cx[2].p1)/2, ax_y+12.0, 'Within-gene share of path length L (median, IQR)', 8.0, ha='center')
F.text((dx.p0+dx.p1)/2, ax_y+12.0, 'Endpoint retention D/L = 1 − T', 8.0, ha='center')
F.key_row(None, ax_y+23.0, [('marker', (STUDY[r][2], STUDY[r][1]), STUDY[r][0]) for r in ROUTES], gap=8.0, center=x+w/2)
F.values['C'] = Cv

# ================================================================ D  amplitude sensitivity of T
LEG3 = Y3+H3-7.4
x, y, w, h = F.panel('D', X0, Y3, WL, H3)
amp = D.df('D07E', 'profiles/trajectory_amplitude_sensitivity.tsv')
amp = amp[amp.measure == 'path_length_L']
Dv = {}
sy = Lin(0, 1, y+62, y+13)
pw = (w-24-2*11-4.5)/3
for i, r in enumerate(ROUTES):
    s_ = Lin(0, 0.2, x+24+(pw+11)*i, x+24+(pw+11)*i+pw)
    nm, col, mk = STUDY[r]
    F.marker(s_.p0+2.8, y+4.2, mk, 'mark', fc=col); F.text(s_.p0+7.6, y+0.4, nm, 7.0)
    F.lineplot(s_, sy, r)
    F.line([s_.p0, s_.p0], [sy.p0, sy.p1], lw='axis', color=AXC)
    if i == 0:
        F.yaxis(sy, s_.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Median T', line=False)
    Dv[r] = {}
    for c in CL:
        m = amp[(amp.route == r) & (amp.cluster == c)].sort_values('cutoff')
        F.line([s_(v) for v in m.cutoff], [sy(v) for v in m.median_T], lw='data', color=INK, ls=CLS[c])
        Dv[r][c] = dict(cutoff=[round(v, 2) for v in m.cutoff], median_T=[float(v) for v in m.median_T], retained_n=[int(v) for v in m.retained_n])
    F.xaxis(s_, y+64, [0, 0.1, 0.2], ['0', '.1', '.2'])
F.text(x+24+(3*pw+22)/2, y+76.0, 'Minimum path length L', 8.0, ha='center')
F.key_row(None, LEG3, [('line', (INK, CLS[c], 'data'), c) for c in CL], gap=8.0, center=x+24+(3*pw+22)/2)
F.values['D'] = Dv

# ================================================================ E  mouse dominant-stage agreement
x, y, w, h = F.panel('E', XB, Y3, WB, H3)
dm = D.df('D07F', 'panel_statistics/dominant-stage-marginal-expectation.tsv').set_index('cluster')
ex = Lin(0.2, 0.8, x+36, x+w-4)
F.text(x+30, y+3.0, 'n', 7.0, color=MUTED, ha='right')
Ev = {}
for i, c in enumerate(CL):
    yy = y+22.6+14.7*i; r_ = dm.loc[c]
    F.rect(x, yy-6.5, w, 13.0, fc=BAND)
    F.text(x+1.0, yy-3.8, c, 8.0, font='Bol'); F.text(x+30, yy-3.4, f'{int(r_.n_genes)}', 7.0, color=MUTED, ha='right')
    F.line([ex(r_.marginal_expected), ex(r_.observed)], [yy, yy], lw='data', color=NULLBAR)
    F.marker(ex(r_.marginal_expected), yy, 'D', 'mark', fc='white', ec=INK)
    F.marker(ex(r_.observed), yy, 'o', 'mark', fc=INK)
    Ev[c] = dict(n=int(r_.n_genes), observed=float(r_.observed), marginal_expected=float(r_.marginal_expected))
F.xaxis(ex, y+64, [0.2, 0.4, 0.6, 0.8], ['.2', '.4', '.6', '.8'], 'Same dominant step')
F.key_row(None, LEG3, [('marker', ('o', INK), 'Observed'), ('marker', ('D', 'white', INK), 'Marginal')], gap=6.0, center=(ex.p0+ex.p1)/2)
F.values['E'] = Ev

F.save(OUT, D)
