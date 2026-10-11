"""Supplementary Figure D08 (printed number from supdraw.SNUM) — local APA structure, neighbourhood and panel controls, and
APA–GEX profile agreement (R2), polyApipe v2 primary.

a  excess Moran's I for 7 studies x 3 classes under six null configurations (k = 20; same-node and fixed-panel tests at
   threshold .8); cells with q >= .05 left white
b  neighbourhood size k (10–30) for the primary conditional configuration (coverage + library + state)
c  per-gene APA–GEX profile correlation (Pearson) and adjacent-transition same-direction fraction, 3 fitting routes x 3 classes
d  fixed-gene-panel size against complete metacells over thresholds .7/.8/.9
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D08'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 232.0, 132.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
CL = ['C1', 'C2', 'C3']
STUDIES = ['GSE118258', 'GSE150311', 'GSE178325', 'GSE221739', 'GSE103221', 'GSE115943', 'GSE171127']
TESTS = [('total_unrestricted', 'Total', ''), ('state_residual_within_state', 'State', 'resid.'),
         ('state_residual_within_library_state', 'Lib. ×', 'state'), ('state_coverage_library_state', 'Cov., lib.,', 'state'),
         ('original_same_nodes', 'Same', 'nodes'), ('fixed_panel', 'Fixed', 'panel')]
PRIMARY = 'state_coverage_library_state'
m = D.df('D08A', 'moran/moran_results.tsv')

# ================================================================ a  six null configurations
WA = 300.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
labw = max(textwidth(DATASET[s][0], 7.0) for s in STUDIES)
xl = x+0.4; xcl = xl+7.6+labw+4.0
nw = 21.0
mx0 = xcl+textwidth('C3', 7.0)+4.0; cw = (x+w-nw-mx0)/len(TESTS)
rh = 9.0; r0 = y+21.0
for k, (t, l1, l2) in enumerate(TESTS):
    bold = 'Bol' if t == PRIMARY else 'Reg'
    F.text(mx0+cw*(k+0.5), y+3.0, l1, 7.0, ha='center', font=bold); F.text(mx0+cw*(k+0.5), y+10.6, l2, 7.0, ha='center', font=bold)
F.text(x+w, y+3.0, 'Meta-', 7.0, ha='right', color=MUTED); F.text(x+w, y+10.6, 'cells', 7.0, ha='right', color=MUTED)
SMAX = 0.6
Av = {}; nsig = {t: 0 for t, *_ in TESTS}
for i, (s, c) in enumerate([(s, c) for s in STUDIES for c in CL]):
    yy = r0+rh*i+(3.0 if i >= 12 else 0.0)
    if c == 'C1':
        F.rect(xl, yy+rh*1.5-2.8, 5.6, 5.6, fc=DATASET[s][1]); F.text(xl+7.6, yy+rh*1.5-3.4, DATASET[s][0], 7.0)
    F.text(xcl, yy+rh/2-3.4, c, 7.0, color=MUTED)
    Av[f'{s}|{c}'] = {}
    for k, (t, *_) in enumerate(TESTS):
        q = m[(m.study_id == s) & (m.cluster == c) & (m.test == t) & (m.k == 20)]
        if t in ('original_same_nodes', 'fixed_panel'): q = q[(q.threshold-0.8).abs() < 1e-6]
        assert len(q) == 1, (s, c, t); q = q.iloc[0]
        sig = q.q_BH < 0.05; nsig[t] += int(sig)
        fc = seq(min(q.excess_i, SMAX), 0, SMAX) if sig else 'white'
        F.rect(mx0+cw*k+0.4, yy+0.4, cw-0.8, rh-0.8, fc=fc, ec='none' if sig else AXC, lw=0 if sig else 'ref')
        F.text(mx0+cw*(k+0.5), yy+rh/2-3.2, f'{q.excess_i:.2f}'.replace('0.', '.').replace('-.', '−.'), 6.8, ha='center',
               color=on(fc) if sig else MUTED)
        Av[f'{s}|{c}'][t] = dict(excess_I=float(q.excess_i), q_BH=float(q.q_BH), z=float(q.permutation_z), n_metacells=float(q.n_metacells))
    nm = m[(m.study_id == s) & (m.cluster == c) & (m.test == PRIMARY) & (m.k == 20)].n_metacells.iloc[0]
    F.text(x+w, yy+rh/2-3.4, f'{int(nm):,}', 7.0, ha='right', color=MUTED)
yb = r0+rh*21+3.0
for k, (t, *_) in enumerate(TESTS):
    F.text(mx0+cw*(k+0.5), yb+1.6, f'{nsig[t]}/21', 7.0, ha='center', font='Bol' if t == PRIMARY else 'Reg')
F.text(mx0-3.0, yb+1.6, 'q < .05', 7.0, ha='right', color=MUTED)
F.cbar_key((mx0+x+w-nw)/2, yb+12.0, 'Excess Moran’s I', '0', '≥ .6', lambda f_: seq(f_*SMAX, 0, SMAX), wbar=60)
F.values['a'] = dict(rows=Av, q_lt_05=nsig)

# ================================================================ b  neighbourhood size
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
kx = Lin(10, 30, x+26, x+w-4.0)
Bv = {}
ph = 55.4; tops = [y+3.0+(ph+17.6)*i for i in range(3)]
for i, c in enumerate(CL):
    sy = Lin(0, 0.35, tops[i]+ph, tops[i])
    F.lineplot(kx, sy, c)
    F.text(kx.p1, tops[i]+0.6, c, 7.5, ha='right', font='Bol')
    F.yaxis(sy, kx.p0, [0, 0.1, 0.2, 0.3], ['0', '.1', '.2', '.3'])
    F.line([kx.p0, kx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    for s in STUDIES:
        q = m[(m.study_id == s) & (m.cluster == c) & (m.test == PRIMARY)].sort_values('k')
        F.line([kx(v) for v in q.k], [sy.clip(v) for v in q.excess_i], lw='thin', color=DATASET[s][1])
        Bv[f'{s}|{c}'] = {int(k_): float(v) for k_, v in zip(q.k, q.excess_i)}
    ab = F.xaxis(kx, sy.p0+1.0, [10, 15, 20, 30], line=True)
F.vtext(x+10.0, (tops[0]+tops[2]+ph)/2, 'Excess Moran’s I (coverage, library, state)', 8.0)
F.text((kx.p0+kx.p1)/2, ab+1.6, 'Neighbors k', 8.0, ha='center')
F.text((kx.p0+kx.p1)/2, y+h-7.4, 'Line colors: studies as in a', 7.5, ha='center', color=MUTED)
F.values['b'] = Bv

# ================================================================ c  APA–GEX
WC = 290.0
x, y, w, h = F.panel('c', X0, Y2, WC, H2)
gc = D.df('D08C', 'profile_gex/APA_GEX_gene_correlations.tsv.gz')
ad = D.df('D08C', 'profile_gex/adjacent_direction_summary.tsv')
cs = D.path('D08C', 'profile_gex/APA_GEX_correlation_summary.tsv')
import pandas as pd
cs = pd.read_csv(cs, sep='\t'); cs = cs[(cs.grouping == 'species;route;cluster') & (cs.method == 'Pearson')].set_index(['route', 'cluster'])
RTS = FIT_ROUTES
rx = Lin(-1, 1, x+30, x+170)
fx = Lin(0.3, 0.7, x+186, x+w-4)
TR = [('S_R', 'S→R', SC['R'], 'o'), ('R_E', 'R→E', SC['E'], 's'), ('E_P', 'E→P', SC['P'], 'D')]
F.text((rx.p0+rx.p1)/2, y+3.0, 'Per-gene APA–GEX profile r', 7.5, ha='center')
F.text((fx.p0+fx.p1)/2, y+3.0, 'Same-direction fraction', 7.5, ha='center')
rh = 9.6; r0 = y+21.0
rng = np.random.default_rng(8)
Cv = {}
for i, (r, c) in enumerate([(r, c) for r in RTS for c in CL]):
    yy = r0+rh*i+2.0*(i//3)
    if i % 2 == 0: F.rect(x, yy-rh/2, w, rh, fc=BAND)
    nm_, col, mk = STUDY[r]
    F.marker(x+4.0, yy, mk, 'mark', fc=col); F.text(x+9.5, yy-3.4, c, 7.0)
    v = gc[(gc.route == r) & (gc.cluster == c)].APA_GEX_r.dropna().to_numpy()
    F.points([rx(t) for t in v], yy+rng.uniform(-2.8, 2.8, len(v)), color=col, alpha='dense')
    s_ = cs.loc[(r, c)]
    F.line([rx(s_.q25_r), rx(s_.q75_r)], [yy, yy], lw='bold', color=INK)
    F.line([rx(s_.median_r)]*2, [yy-3.4, yy+3.4], lw='bold', color=INK)
    Cv[f'{r}|{c}'] = dict(n=int(s_.n_estimable), median=float(s_.median_r), q25=float(s_.q25_r), q75=float(s_.q75_r))
    for tk, tl, tc, tm in TR:
        a_ = ad[(ad.route == r) & (ad.cluster == c) & (ad.transition == tk)].iloc[0]
        fr = a_.same_direction/(a_.genes-a_.zero_changes)
        F.marker(fx.clip(fr), yy, tm, 'small', fc=tc)
        Cv[f'{r}|{c}'][tk] = dict(same=int(a_.same_direction), genes=int(a_.genes), zero=int(a_.zero_changes), fraction=float(fr))
yb = r0+rh*8+4.0+rh/2
for s_, v in ((rx, 0.0), (fx, 0.5)):
    F.line([s_(v)]*2, [r0-rh/2, yb], lw='ref', color=AXC, ls='ref')
F.xaxis(rx, yb+0.6, [-1, -0.5, 0, 0.5, 1], ['−1', '−.5', '0', '.5', '1'])
F.xaxis(fx, yb+0.6, [0.3, 0.4, 0.5, 0.6, 0.7], ['.3', '.4', '.5', '.6', '.7'])
F.key_row(None, y+h-7.4, [('line', (INK, None, 'bold'), 'Median, IQR'), ('marker', ('o', AXC), 'One gene')], gap=7.0,
          center=(rx.p0+rx.p1)/2)
F.key_row(None, y+h-7.4, [('marker', (tm, tc), tl) for _, tl, tc, tm in TR], gap=6.0, center=(fx.p0+fx.p1)/2)
F.values['c'] = Cv
F.notes.append('c: points = all genes with an estimable correlation (Pearson over the observed states); median/IQR from the '
               'statistics summary; same-direction fraction excludes genes with a zero change')

# ================================================================ d  fixed-panel feasibility
x, y, w, h = F.panel('d', X0+WC+GAP, Y2, RR-(X0+WC+GAP), H2)
fp = D.df('D08D', 'moran/fixed_panel_feasibility.tsv')
gx = Lin(0, 330, x+28, x+w-1.0)
F.text((gx.p0+gx.p1)/2, y+3.0, 'Fixed gene panel', 7.5, ha='center')
ly = Lin(math.log10(70), math.log10(3000), y+98, y+19)
F.lineplot(gx, ly, 'fixed panel')
F.yaxis(ly, gx.p0, [math.log10(v) for v in (100, 300, 1000, 3000)], ['100', '300', '1k', '3k'], 'Complete metacells', title_y=(ly.p0+ly.p1)/2+3)
Dv = {}
for s in STUDIES:
    for c in CL:
        q = fp[(fp.study_id == s) & (fp.cluster == c)].sort_values('threshold')
        xs_ = [gx(v) for v in q.selected_genes]; ys_ = [ly.clip(math.log10(v)) for v in q.complete_metacells]
        F.line(xs_, ys_, lw='thin', color=DATASET[s][1], alpha='secondary')
        for (xx, yy_), th in zip(zip(xs_, ys_), q.threshold):
            F.marker(xx, yy_, 'o', 'small', fc=DATASET[s][1] if abs(th-0.8) < 1e-6 else 'white', ec=DATASET[s][1])
        Dv[f'{s}|{c}'] = [dict(threshold=round(float(t_), 2), genes=int(g_), metacells=int(n_)) for t_, g_, n_ in zip(q.threshold, q.selected_genes, q.complete_metacells)]
F.xaxis(gx, ly.p0+1.0, [0, 100, 200, 300], title='Fixed-panel genes')
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'Threshold .8'), ('marker', ('o', 'white', INK), '.7 / .9')], gap=6.0, center=(gx.p0+gx.p1)/2)
F.values['d'] = Dv
F.save(OUT, D)
