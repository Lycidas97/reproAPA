"""Supplementary Figure S18 (v11): cell-cycle and CPA-factor expression context of P, and the upstream UGUA test.
a  cell-cycle score difference (mean of 10 S and 8 G2/M markers, ln(1+CP10K)): P–M in the same physical libraries
   (per-library points, route mean) and P–S against the documented source (route value, unpaired libraries)
b  P–S cycle-score difference against the global APA direction (net mean dDWUI, same contrast) per documented-source series
c  CPA-factor expression, difference in mean log2(CPM+1): P–M (same libraries) and P–S (documented source), 11 core factors
d  distal minus proximal UGUA count in [-100, -40) for proximal responders against matched nulls (8 primary tests)
e  within-state partial correlation of cycle score with mean DWUI per route and state (library bootstrap 95% range)
Data: notes/2026-10-10-proliferation-apa, cpa-factor-expression, cfim-sequence-signature (read-only)."""
import numpy as np
from common import *                                     # noqa
from common import supdraw, Reads, name, OUT

R = Reads()
X0, XR = 4.5, PAGE_W-4.5
ROUTES = ['GSE118258', 'GSE150311', '178325_0618', '178325_0809', 'GSE221739', 'CIP', 'WOT']
DOC = ['GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT']
PAIRED = ['GSE118258', 'GSE150311', 'GSE221739', 'CIP', 'WOT']
FACT = ['NUDT21', 'CPSF6', 'CSTF2', 'CSTF2T', 'FIP1L1', 'PCF11', 'CPSF1', 'CPSF2', 'CPSF3', 'CPSF4', 'PABPN1']
sc_ = R.df('notes/2026-10-10-proliferation-apa/state_score_contrasts.tsv')
pl = R.df('notes/2026-10-10-proliferation-apa/paired_library_score_differences.tsv')
ex = R.df('notes/2026-10-10-cpa-factor-expression/expression_contrasts.tsv')
ra = R.df('notes/2026-10-10-proliferation-apa/route_associations.tsv')
sq = R.df('notes/2026-10-10-cfim-sequence-signature/sequence_signature_tests.tsv')
H1, H2, H3 = 107.3, 154.5, 108.1
H = 4.5+H1+GAP+H2+GAP+H3+4.5
F = supdraw.Fig(H, name('cycle'))
y1 = 4.5; y2 = y1+H1+GAP; y3 = y2+H2+GAP

def route_label(x, yy, r, size=7.0):
    F.rect(x, yy-2.8, 5.6, 5.6, fc=STUDY[r][1])
    return F.text(x+7.6, yy-3.3, STUDY[r][0], size)

def cval(route, contrast):
    r = sc_[(sc_.route == route) & (sc_.contrast == contrast) & (sc_.score == 'Cycle')].iloc[0]
    return None if r.status.startswith('absent') or np.isnan(r.delta) else float(r.delta)

# ---------------------------------------------------------------- a
WA = 250.0
x, y, w, h = F.panel('a', X0, y1, WA, H1)
lab_w = 78.0
s1 = Lin(-0.1, 0.6, x+lab_w+2.0, x+lab_w+2.0+74.0)
s2 = Lin(-0.1, 0.6, s1.p1+16.0, x+WA-5.0)
F.text(x+16.0, y+0.6, 'Cell-cycle score difference', 7.5)
F.text((s1.p0+s1.p1)/2, y+11.0, 'P–M, same libraries', 7.0, ha='center', color=MUTED)
F.text((s2.p0+s2.p1)/2, y+11.0, 'P–S, documented source', 7.0, ha='center', color=MUTED)
top, pitch = y+24.0, 9.6
for k, r in enumerate(ROUTES):
    yy = top+k*pitch
    if k % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
    route_label(x+0.2, yy, r)
    lib = pl[(pl.route == r) & (pl.contrast == 'P_M') & (pl.score == 'Cycle')]
    v = cval(r, 'P_M')
    if len(lib):
        F.points([s1(u) for u in lib.delta], [yy]*len(lib), ms='small', color=SC['P'], alpha='secondary')
    if v is not None:
        F.marker(s1(v), yy, '|', 'mark', fc=INK)
    else:
        F.text((s1.p0+s1.p1)/2, yy-3.3, 'no paired libraries', 7.0, ha='center', color=MUTED)
    vd = cval(r, 'P_S_documented_source')
    if vd is not None: F.marker(s2(vd), yy, 'o', 'mark', fc=INK)
    else: F.text((s2.p0+s2.p1)/2, yy-3.3, 'no documented S', 7.0, ha='center', color=MUTED)
    F.values[f'a {r}'] = dict(P_M=v, P_M_libraries=[float(u) for u in lib.delta], P_S_documented=vd)
yb = top+len(ROUTES)*pitch-pitch/2+1.0
for s in (s1, s2):
    F.line([s(0)]*2, [top-pitch/2, yb], lw='ref', color=AXC, ls='ref')
    ya = F.xaxis(s, yb, [0, 0.2, 0.4, 0.6], ['0', '.2', '.4', '.6'])
F.key_row(0, ya+2.5, [('marker', ('o', SC['P']), 'Library'), ('marker', ('|', INK), 'Route mean'), ('marker', ('o', INK), 'Route')],
          center=x+lab_w+(WA-lab_w)/2)
a_end = ya+2.5+7.6

# ---------------------------------------------------------------- b
xb0 = X0+WA+GAP
x, y, w, h = F.panel('b', xb0, y1, XR-xb0, H1)
sx = Lin(-0.1, 0.5, x+27.5, XR-1.0); sy = Lin(-0.04, 0.03, y+H1-19.0, y+12.0)
F.lineplot(sx, sy, label='cycle vs direction', standalone=False)
F.text(x+16.0, y+0.6, 'P–S, documented source', 7.5)
F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
pts = {}
for r in DOC:
    cy = cval(r, 'P_S_documented_source')
    nd = ex[(ex.route == r) & (ex.contrast == 'P_S') & (ex.factor == 'NUDT21')].iloc[0].net_DWUI_mean
    pts[r] = (cy, float(nd))
    F.marker(sx(cy), sy(nd), STUDY[r][2], 'mark', fc=STUDY[r][1])
    F.values[f'b {r}'] = dict(cycle_delta=cy, net_dDWUI=float(nd))
for r, (dx, dy, ha) in {'GSE118258': (4.0, -3.4, 'left'), '178325_0618': (4.0, -8.0, 'left'), '178325_0809': (4.0, 0.0, 'left'),
                        'CIP': (4.0, -3.4, 'left'), 'WOT': (-4.0, -3.4, 'right')}.items():
    cy, nd = pts[r]
    F.text(sx(cy)+dx, sy(nd)+dy, STUDY[r][0], 7.0, ha=ha)
F.yaxis(sy, sx.p0, [-0.04, -0.02, 0, 0.02], ['−.04', '−.02', '0', '.02'], title='Global direction')
F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref')
F.xaxis(sx, sy.p0, [0, 0.2, 0.4], ['0', '.2', '.4'], title='Cell-cycle score difference')

# ---------------------------------------------------------------- c
WC = 262.0
x, y, w, h = F.panel('c', X0, y2, WC, H2)
F.text(x+16.0, y+0.6, 'CPA-factor expression, Δ mean log2(CPM + 1)', 7.5)
cx0, cw, ch = x+40.0, 21.65, 8.6
cols = [('P_M', r) for r in PAIRED]+[('P_S', r) for r in DOC]
CX = lambda j: cx0+j*cw+(5.0 if j >= 5 else 0.0)
F.text(CX(0)+2.5*cw, y+11.4, 'P–M, same libraries', 7.0, ha='center', color=MUTED)
F.text(CX(5)+2.5*cw, y+11.4, 'P–S, documented source', 7.0, ha='center', color=MUTED)
ht = y+21.0
for j, (con, r) in enumerate(cols):
    F.rect(CX(j)+0.5, ht, cw-1.5, 4.0, fc=STUDY[r][1])
top = ht+6.0
LIM = 1.5
for i, fct in enumerate(FACT):
    yy = top+i*ch
    F.text(x+0.2, yy+1.0, fct, 7.0)
    for j, (con, r) in enumerate(cols):
        row = ex[(ex.route == r) & (ex.contrast == con) & (ex.factor == fct)].iloc[0]
        v = row.delta_log2CPM1
        if row.status == 'unavailable' or np.isnan(v):
            F.rect(CX(j)+0.5, yy, cw-1.5, ch-1.0, fc=MISSING); v = None
        else:
            F.rect(CX(j)+0.5, yy, cw-1.5, ch-1.0, fc=supdraw.diverge(v, LIM))
        F.values[f'c {con}|{r}|{fct}'] = None if v is None else float(v)
yk = top+len(FACT)*ch+3.5
F.cbar_key(cx0+(CX(9)+cw-cx0)/2, yk, 'Δ', '−1.5', '+1.5', lambda u: supdraw.diverge(-LIM+2*LIM*u, LIM), wbar=60.0)
keys = [('block', STUDY[r][1], STUDY[r][0]) for r in ['GSE118258', 'GSE150311', '178325_0618', '178325_0809', 'GSE221739', 'CIP', 'WOT']]
F.key_row(0, yk+10.5, keys[:4], center=x+WC/2)
F.key_row(0, yk+21.0, keys[4:], center=x+WC/2)
c_end = yk+21.0+7.6

# ---------------------------------------------------------------- d  UGUA
xd0 = X0+WC+GAP
x, y, w, h = F.panel('d', xd0, y2, XR-xd0, H2)
F.text(x+16.0, y+0.6, 'Upstream UGUA, distal − proximal', 7.5)
p = sq[(sq.primary == True) & (sq.motif == 'TGTA') & (sq.window == 'upstream40_100')]
TESTS = [('paired_strict', 'human', 'proximal_vs_nonrecurrent'), ('paired_strict', 'human', 'proximal_vs_distal'),
         ('paired_strict', 'mouse', 'proximal_vs_nonrecurrent'), ('paired_strict', 'mouse', 'proximal_vs_distal')]
BG = [('paired_strict', 'Common orthologs'), ('species_strict', 'Each species')]
bgs = sorted(p.background.unique())
assert len(p) == 8, (len(p), bgs)
other = [b for b in bgs if b != 'paired_strict'][0]
sd = Lin(-1.0, 1.0, x+70.0, XR-20.0)
F.text(XR-0.2, y+11.0, 'q', 7.0, ha='right', color=MUTED)
top, pitch = y+28.0, 11.2
k = 0
for bi, (bg, bl) in enumerate(((('paired_strict'), 'Common orthologs'), (other, 'Each species'))):
    yh = top+bi*(4*pitch+12.0)-12.0
    F.text(x+0.2, yh, bl+f', {int(p[p.background == bg].n_genes.iloc[0]) if bg == "paired_strict" else ""}'.rstrip(', '), 7.0, color=MUTED)
    for sp in ('human', 'mouse'):
        for comp, cl in (('proximal_vs_nonrecurrent', 'vs not recurrent'), ('proximal_vs_distal', 'vs distal')):
            r = p[(p.background == bg) & (p.species == sp) & (p.comparison == comp)].iloc[0]
            yy = top+bi*(4*pitch+12.0)+k % 4*pitch
            k += 1
            if (k-1) % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
            F.text(x+0.2, yy-3.3, ('Hs ' if sp == 'human' else 'Mm ')+cl, 7.0)
            F.line([sd(r.null_p025), sd(r.null_p975)], [yy, yy], lw='range', color=NULLBAR)
            F.line([sd(r.analytic_expected)]*2, [yy-2.6, yy+2.6], lw='thin', color=MUTED)
            F.marker(sd(r.observed), yy, 'o', 'mark', fc=INK if r.q_family < 0.05 else 'white', ec=INK)
            F.text(XR-0.2, yy-3.3, f'{r.q_family:.2f}'.lstrip('0'), 7.0, ha='right')
            F.values[f'd {bg}|{sp}|{comp}'] = dict(observed=float(r.observed), expected=float(r.analytic_expected),
                                                  null=[float(r.null_p025), float(r.null_p975)], p=float(r.p_enrichment), q8=float(r.q_family),
                                                  n_proximal=int(r.n_proximal), n_control=int(r.n_control))
yb = top+(4*pitch+12.0)+3*pitch+pitch/2+1.0
F.line([sd(0)]*2, [top-pitch/2, yb], lw='ref', color=AXC, ls='ref')
yy = F.xaxis(sd, yb, [-1, 0, 1], ['−1', '0', '1'], title='UGUA count difference')
F.key_row(0, yy+2.5, [('marker', ('o', 'white', INK), 'q ≥ .05'), ('line', (NULLBAR, None, 'range'), 'Null 95%'), ('line', (MUTED, None, 'thin'), 'Expected')],
          center=(x+XR)/2)
d_end = yy+2.5+7.6

# ---------------------------------------------------------------- e  within-state association
x, y, w, h = F.panel('e', X0, y3, XR-X0, H3)
F.text(x+16.0, y+0.6, 'Within-state partial correlation of cell-cycle score with mean DWUI', 7.5)
lab_w = 82.0
fw = (XR-x-lab_w-2*10.0)/3
top, pitch = y+24.0, 8.4
for k, r in enumerate(ROUTES):
    yy = top+k*pitch
    if k % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
    route_label(x+0.2, yy, r)
ye = 0
for si, st in enumerate(('S', 'M', 'P')):
    s = Lin(-0.6, 0.6, x+lab_w+si*(fw+10.0), x+lab_w+si*(fw+10.0)+fw)
    F.text((s.p0+s.p1)/2, y+12.0, st, 7.5, ha='center', color=SC[st])
    for k, r in enumerate(ROUTES):
        yy = top+k*pitch
        a = ra[(ra.route == r) & (ra.score == 'Cycle') & (ra.APA_variant == 'mean_DWUI') & (ra.scope == st)].iloc[0]
        if a.status == 'unavailable':
            F.text(s(0), yy-3.3, '–', 7.0, ha='center', color=MUTED); continue
        if not np.isnan(a.library_boot_p025):
            F.line([s.clip(a.library_boot_p025), s.clip(a.library_boot_p975)], [yy, yy], lw='bold', color=NULLBAR)
        F.marker(s(a.partial_correlation), yy, 'o', 'small', fc=SC[st], ec=SC[st])
        F.values[f'e {r}|{st}'] = dict(r=float(a.partial_correlation), n_libraries=int(a.n_libraries),
                                       boot=None if np.isnan(a.library_boot_p025) else [float(a.library_boot_p025), float(a.library_boot_p975)])
    yb = top+len(ROUTES)*pitch-pitch/2+1.0
    F.line([s(0)]*2, [top-pitch/2, yb], lw='ref', color=AXC, ls='ref')
    ye = max(ye, F.xaxis(s, yb, [-0.5, 0, 0.5], ['−.5', '0', '.5'], title='Partial r' if si == 1 else None))
F.key_row(0, ye+2.5, [('line', (NULLBAR, None, 'bold'), 'Library bootstrap 95%')], center=(x+lab_w+XR)/2)
print('ends', round(a_end-y1, 1), round(c_end-y2, 1), round(d_end-y2, 1), round(ye+2.5+7.6-y3, 1))
F.save(OUT, R)
