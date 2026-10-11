"""Supplementary Figure D06 (printed number from supdraw.SNUM) — resolution, shape robustness and route hold-out (R2.1–R2.2), polyApipe v2 primary.

A  silhouette by k and the normalized class centres for k = 2–6; up-then-down (reversal) centres drawn in ink
B  clustering stability (median ARI over replicates) across nine perturbations and k = 2–6
C  observed route median T per class against the fixed-class null and 500 full-background reselection draws
D  class transitions when a route is held out (full-fit class -> projected class), with agreement and ARI
E  held-out route mean minus aligned refit centre per state, with the centroid correlation
Layout follows STYLE.zh.md.
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

D = Data(); OUT = outdir('D06')
import os
if os.environ.get('SUP_OUT'): OUT = Path(os.environ['SUP_OUT'])
X0, R = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 95.0, 138.5, 111.5
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname('D06'))
CL = ['C1', 'C2', 'C3']; ST = ['S', 'R', 'E', 'P']
ROUTES = FIT_ROUTES
REV = lambda c: max(c['R'], c['E'])-c['S'] >= 0.15 and c['P'] <= max(c['R'], c['E'])-0.5


colorbar = F.cbar

# ================================================================ A  resolution and class shapes
x, y, w, h = F.panel('A', X0, Y1, R-X0, H1)
LEGA = y+87.0
rc = D.df('D06A', 'profiles/resolution_comparison.tsv').sort_values('k')
scn = D.df('D06A', 'profiles/shape_centers.tsv')
ks = Lin(2, 6, x+30, x+100); ss = Lin(0.3, 0.6, y+64, y+8)
F.yaxis(ss, x+26, [0.3, 0.4, 0.5, 0.6], ['.3', '.4', '.5', '.6'], 'Silhouette', title_y=y+40)
F.lineplot(ks, ss, 'silhouette', standalone=False)
F.line([ks(k) for k in rc.k], [ss(v) for v in rc.species_balanced_silhouette], lw='thin', color=AXC)
Av = dict(silhouette={}, centres={}, reversal_rule='max(R,E) − S ≥ 0.15 and P ≤ max(R,E) − 0.5 (normalized centre)')
for k, v in zip(rc.k, rc.species_balanced_silhouette):
    F.marker(ks(k), ss(v), 'o', 'mark', fc=INK if k == 3 else 'white', ec=INK)
    Av['silhouette'][int(k)] = float(v)
F.xaxis(ks, y+66, [2, 3, 4, 5, 6], title='k')
cy = Lin(-0.85, 0.85, y+64, y+10)
c0 = x+130; pw = (x+w-c0-4*7.0)/5
for i, k in enumerate([2, 3, 4, 5, 6]):
    cx = Lin(0, 3, c0+(pw+7.0)*i+3, c0+(pw+7.0)*i+pw-3)
    F.text((cx.p0+cx.p1)/2, y+0.4, f'k = {k}', 7.5, ha='center', font='Bol' if k == 3 else 'Reg')
    F.lineplot(Lin(0, 1, cx.p0-3, cx.p1+3), cy, f'centres k={k}', standalone=False)
    F.line([cx.p0-3, cx.p1+3], [cy(0)]*2, lw='ref', color=AXC, ls='ref')
    if i == 0:
        F.yaxis(cy, cx.p0-4, [-0.5, 0, 0.5], ['−.5', '0', '.5'], 'Center')
    sub = scn[scn.k == k]; Av['centres'][k] = {}; LAB = []
    for _, r_ in sorted(sub.iterrows(), key=lambda t: REV(t[1])):
        rev = REV(r_)
        F.line([cx(j) for j in range(4)], [cy(r_[s]) for s in ST], lw='bold' if rev else 0.8, color=INK if rev else NULLW)
        Av['centres'][k][r_.cluster] = dict(**{s: float(r_[s]) for s in ST}, reversal=bool(rev))
        if k == 3: LAB.append([cy(r_['P'])-3.4, r_.cluster, rev])
    LAB.sort()
    for m_ in range(1, len(LAB)):
        LAB[m_][0] = max(LAB[m_][0], LAB[m_-1][0]+7.6)
    for j, s in enumerate(ST):
        F.line([cx(j)]*2, [y+66, y+68], lw='axis', color=AXC)
        F.text(cx(j), y+69.0, s, 7.5, ha='center', color=SC[s])
F.key_row(None, LEGA, [('marker', ('o', INK), 'k used'), ('marker', ('o', 'white', INK), 'other k')], gap=7.0, center=(ks.p0+ks.p1)/2)
F.key_row(None, LEGA, [('line', (INK, None, 'bold'), 'Reversing center'), ('line', (NULLW, None, 'thin'), 'Other center')], gap=8.0,
          center=(c0+x+w)/2)
F.notes.append('A: reversal highlight is a display rule on the statistics-side centres: '+Av['reversal_rule'])
F.values['A'] = Av

# ================================================================ B  clustering stability across perturbations
WB = 200.0
LEG2 = Y2+H2-7.4
x, y, w, h = F.panel('B', X0, Y2, WB, H2)
cr = D.df('D06A', 'profile_controls/clustering_robustness.tsv')
KINDS = [('seed', 'Random seed'), ('80pct_gene_subsample', '80% of genes'), ('library_reweight', 'Library reweight'),
         ('leave_route_CIP', 'Leave out CIP'), ('leave_route_WOT', 'Leave out 2i'), ('S_origin', 'S at origin'),
         ('expanded_primary_p08_depth20', 'Expanded layer'), ('supported50_p08_depth40', '40 PAS UMIs'),
         ('supported50_p08_depth80', '80 PAS UMIs')]
kx0 = x+1.0+textwidth('Library reweight', 7.0)+3.0; cw = (x+w-kx0)/5; chh = 11.2; B0 = y+15.5          # first row below the letter clear zone
for i, k in enumerate(range(2, 7)):
    F.text(kx0+cw*i+cw/2, y+3.0, f'k = {k}', 7.5, ha='center', font='Bol' if k == 3 else 'Reg')
Bv = {}
for j, (kind, lab) in enumerate(KINDS):
    yy = B0+chh*j
    sub = cr[cr.kind == kind]
    meth = 'centered' if (sub.method == 'centered').any() else 'S_zero'
    sub = sub[sub.method == meth]
    F.text(kx0-2.0, yy+2.3, lab, 7.0, ha='right')
    Bv[kind] = dict(method=meth, replicates=int(sub.replicate.nunique()), median_ARI={})
    for i, k in enumerate(range(2, 7)):
        v = float(sub[sub.k == k].ARI.median())
        F.rect(kx0+cw*i, yy, cw-1.0, chh-1.0, fc=seq(v, 0.3, 1.0))
        F.text(kx0+cw*i+(cw-1.0)/2, yy+2.3, f'{v:.2f}'.replace('0.', '.'), 7.0, ha='center', color=on(seq(v, 0.3, 1.0)))
        Bv[kind]['median_ARI'][k] = v
F.rect(kx0+cw-0.5, B0-0.5, cw, chh*len(KINDS), fc='none', ec=INK, lw='thin')
ybar = B0+chh*len(KINDS)+5.0
colorbar(kx0+cw*1.0, ybar, cw*3-1.0, '.3', '1', lambda f: seq(0.3+0.7*f, 0.3, 1.0))
F.text(kx0+cw*2.5-0.5, LEG2, 'Median ARI against the full fit', 7.5, ha='center')
F.notes.append('B: cells are medians over replicates of clustering_robustness (centered method where run, S_zero otherwise); display summary only')
F.values['B'] = Bv

# ================================================================ C  observed T against null and reselection draws
xc = X0+WB+GAP; wc = R-xc
x, y, w, h = F.panel('C', xc, Y2, wc, H2)
fc_ = D.df('D06B', 'profile_controls/trajectory_fixed_class_T.tsv').set_index(['route', 'cluster'])
dr = D.df('D06B', 'profile_controls/trajectory_reselected_draws.tsv.gz')
tx = Lin(0, 0.8, x+40, x+w-3)
F.text((tx.p0+tx.p1)/2, y+3.0, 'Route median T', 7.5, ha='center', color=MUTED)
rows = []; yy = y+21.1          # first band starts below the letter clear zone
for c in CL:
    for r in ROUTES:
        rows.append((c, r, yy)); yy += 9.22
    yy += 4.0
rng = np.random.default_rng(6)
for ci, c in enumerate(CL):
    ys_ = [r[2] for r in rows if r[0] == c]
    if ci % 2 == 0: F.rect(x, ys_[0]-5.0, w, ys_[-1]-ys_[0]+10.0, fc=BAND)
    F.text(x+1.0, ys_[0]-3.8, c, 8.0, font='Bol')
Cv = {}
for c, r, yy in rows:
    nm, col, mk = STUDY[r]
    F.marker(x+26, yy, mk, 'mark', fc=col)
    f_ = fc_.loc[(r, c)]
    F.line([tx(f_.null_lo), max(tx(f_.null_hi), tx(f_.null_lo)+1.2)], [yy, yy], lw='range', color=NULLBAR)
    d_ = dr[(dr.route == r) & (dr.cluster == c)].median_T.to_numpy()
    F.points([tx(v) for v in d_], yy+rng.uniform(-3.0, 3.0, len(d_)), color=NULLBAR, alpha='secondary')
    F.marker(tx(f_.observed_median_T), yy, 'o', 'mark', fc=INK)
    Cv[f'{r}|{c}'] = dict(observed=float(f_.observed_median_T), null=[float(f_.null_lo), float(f_.null_hi)],
                          draws=int(len(d_)), draws_max=float(d_.max()), draws_at_zero=int((d_ == 0).sum()))
F.xaxis(tx, rows[-1][2]+6.0, [0, 0.2, 0.4, 0.6, 0.8], ['0', '.2', '.4', '.6', '.8'], 'T = 1 − D/L')
F.key_row(None, LEG2, [('marker', ('o', INK), 'Observed'), ('line', (NULLBAR, None, 'range'), 'Fixed-class null'),
                        ('marker', ('o', NULLBAR), 'Reselection draws')] +
          [('marker', (STUDY[r][2], STUDY[r][1]), STUDY[r][0].split('–')[0].replace('2019', '').replace('2020', '')) for r in []], gap=7.0, center=(tx.p0+tx.p1)/2)
F.values['C'] = Cv

# ================================================================ D  leave-route class transitions
WD = 246.0
LEG3 = Y3+H3-7.4
x, y, w, h = F.panel('D', X0, Y3, WD, H3)
fl = D.df('D06D', 'supplement_statistics/LOO-assignment-flows.tsv')
lp = D.df('D06C', 'profiles/leave_route_projection.tsv').set_index('heldout_route')
x_m = x+40.0; gm = 7.0; cs = (x+w-x_m-2*gm+0.8)/9; bw = 3*cs+gm
F.text(x+1.0, y+12+3*cs+9.0, 'Full fit', 7.0); F.text(x+1.0, y+12+3*cs+18.0, 'Remaining', 7.0, color=MUTED)
Dv = {}
for i, r in enumerate(ROUTES):
    nm, col, mk = STUDY[r]
    mx = x_m+bw*i; my = y+12
    F.marker(mx+2.8, y+4.2, mk, 'mark', fc=col); F.text(mx+7.6, y+0.4, nm, 7.0)
    sub = fl[fl.heldout_route == r]
    tot = {c: int(sub[sub.full_cluster == c].genes.sum()) for c in CL}
    for a_, ca in enumerate(CL):
        if i == 0: F.text(mx-1.8, my+cs*a_+cs/2-3.2, ca, 7.0, ha='right', color=MUTED)
        for b_, cb in enumerate(CL):
            g = sub[(sub.full_cluster == ca) & (sub.projected_cluster == cb)].genes
            n = int(g.iloc[0]) if len(g) else 0
            fr = n/tot[ca] if tot[ca] else 0
            F.rect(mx+cs*b_, my+cs*a_, cs-0.8, cs-0.8, fc=seq(fr, 0, 1), ec=INK if a_ == b_ else 'none', lw='thin' if a_ == b_ else 0)
            F.text(mx+cs*b_+(cs-0.8)/2, my+cs*a_+cs/2-3.6, str(n), 7.0, ha='center', color=on(seq(fr, 0, 1)))
    for b_, cb in enumerate(CL):
        F.text(mx+cs*b_+(cs-0.8)/2, my+3*cs+0.6, cb, 7.0, ha='center', color=MUTED)
    q = lp.loc[r]
    F.text(mx, my+3*cs+9.0, f'{100*q.assignment_agreement_vs_full:.1f}%, ARI {q.ARI_vs_full:.2f}'.replace(' 0.', ' .'), 7.0)
    rem = q.assignment_agreement_vs_remaining_route
    F.text(mx, my+3*cs+18.0, ('n/a' if math.isnan(rem) else f'{100*rem:.1f}%, ARI {q.ARI_vs_remaining_route:.2f}'.replace(' 0.', ' .')),
           7.0, color=MUTED)
    Dv[r] = dict(flows={f'{a}->{b}': int(v) for a, b, v in zip(sub.full_cluster, sub.projected_cluster, sub.genes)},
                 agreement_full=float(q.assignment_agreement_vs_full), ARI_full=float(q.ARI_vs_full),
                 agreement_remaining=None if math.isnan(rem) else float(rem))
F.cbar_key((x_m+x_m+3*bw-gm)/2, LEG3, 'Row fraction', '0', '1', lambda f: seq(f, 0, 1))
F.values['D'] = Dv

# ================================================================ E  held-out minus aligned centre
xe = X0+WD+GAP; we = R-xe
x, y, w, h = F.panel('E', xe, Y3, we, H3)
cd = D.df('D06E', 'supplement_statistics/LOO-center-display.tsv')
lc = D.df('D06E', 'profiles/leave_route_centroids.tsv').set_index(['heldout_route', 'cluster'])
ex0 = x+23; ecw = (w-23-17)/4
for j, s in enumerate(ST):
    F.text(ex0+ecw*j+ecw/2, y+3.0, s, 7.5, ha='center', color=SC[s])
F.text(x+w, y+3.0, 'r', 7.5, ha='right', color=MUTED)
Ev = {}; yy = y+15.5; rh = 9.0          # first row below the letter clear zone
for ci, c in enumerate(CL):
    F.text(x+1.0, yy+0.6, c, 8.0, font='Bol')
    for r in ROUTES:
        nm, col, mk = STUDY[r]
        F.marker(x+16, yy+4.4, mk, 'mark', fc=col)
        for j, s in enumerate(ST):
            v = float(cd[(cd.heldout_route == r) & (cd.cluster == c) & (cd.state == s)].difference.iloc[0])
            F.rect(ex0+ecw*j, yy, ecw-0.8, rh-0.8, fc=diverge(v, 0.2))
            F.text(ex0+ecw*j+(ecw-0.8)/2, yy+1.4, num(v, 2, plus=True).replace('0.', '.'), 6.8, ha='center',
                   color=on(diverge(v, 0.2)))
            Ev.setdefault(f'{r}|{c}', {})[s] = v
        rr = float(lc.loc[(r, c)].centroid_correlation)
        F.text(x+w, yy+1.4, f'{rr:.3f}'.replace('0.', '.'), 7.0, ha='right')
        Ev[f'{r}|{c}']['r'] = rr
        yy += rh
    yy += 2.5
F.cbar_key(ex0+ecw*2, LEG3, 'Held-out − refit', '−.2', '+.2', lambda f: diverge(-0.2+0.4*f, 0.2), wbar=60.0)
F.values['E'] = Ev

F.save(OUT, D)
