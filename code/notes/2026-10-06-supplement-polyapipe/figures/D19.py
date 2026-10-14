"""Supplementary Figure D19 (printed number from supdraw.SNUM) — cis features of the cross-species validation model (R4),
polyApipe v2 primary.

a  definitions of the five cis features and the three model groups they enter (PAS, Architecture, Sequence); schematic,
   no measured coordinates
b  model cohorts: genes with the validation-species effect in the same direction (Y = 1) and all others (Y = 0)
c  weighted IQR and median of each feature in the descriptive matched cohorts (Shared, Recurrent-only, Matched background)
d  adjusted logistic coefficients of the five features in the All-cis model (primary setting: supported layer, τ = .02)
e  pairwise Spearman correlations between the five features in each model cohort
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D19'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 98.0, 134.0, 125.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname(FID))
CO = [('human', 'P-S'), ('human', 'P-M'), ('mouse', 'P-S'), ('mouse', 'P-M')]
CL = {c: f'{c[0].capitalize()} {c[1].replace("-", "–")}' for c in CO}
FEAT = ['driver_conserved', 'frac_conserved_discovery_pas', 'A_g', 'phastcons_driver_mean', 'motif_rank_driver_mean']
FL = {'driver_conserved': 'Driver', 'frac_conserved_discovery_pas': 'PAS fraction', 'A_g': 'A_g',
      'phastcons_driver_mean': 'phastCons', 'motif_rank_driver_mean': 'Motif rank'}
FS = {'driver_conserved': 'Driver', 'frac_conserved_discovery_pas': 'PAS frac.', 'A_g': 'A_g',
      'phastcons_driver_mean': 'phastCons', 'motif_rank_driver_mean': 'Motif'}
PRIMARY = dict(layer='supported50', tau=0.02)
WIN = SHADE

# ================================================================ a  feature definitions (schematic)
WA = 150.0   # 2026-10-12: feature definitions and worked-example values moved to the caption; schematic only
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
ax0, ax1 = x+30.0, x+118.0
YH, YM = y+30.0, y+56.0
F.text(x+0.4, YH-3.4, 'Human', 7.0); F.text(x+0.4, YM-3.4, 'Mouse', 7.0)
for yy in (YH, YM):
    F.line([ax0, ax1], [yy, yy], lw='thin', color=AXC)
hs = [0.06, 0.30, 0.52, 0.78, 0.95]; ms_ = [0.10, 0.36, 0.60, 0.86]
matched = [(0, 0), (2, 2)]
gx = lambda f: ax0+f*(ax1-ax0)
F.rect(gx(hs[0])-4.0, YH-6.0, 6.0, 12.0, fc=WIN)            # sequence window around a driver (−100 to +50 nt)
for i, f in enumerate(hs):
    m = any(i == a for a, _ in matched)
    F.line([gx(f)]*2, [YH-4.0, YH+4.0], lw='bold', color=INK if m else AXC)
for j, f in enumerate(ms_):
    m = any(j == b for _, b in matched)
    F.line([gx(f)]*2, [YM-4.0, YM+4.0], lw='bold', color=INK if m else AXC)
for a, b in matched:
    F.line([gx(hs[a]), gx(ms_[b])], [YH+5.0, YM-5.0], lw='ref', color=INK, ls='ref')
F.marker(gx(hs[0]), YH-10.0, '^', 'small', fc=INK); F.text(gx(hs[0])+3.0, YH-13.6, 'gain', 6.8, color=MUTED)
F.marker(gx(hs[3]), YH-10.0, 'v', 'small', fc=INK); F.text(gx(hs[3])+3.0, YH-13.6, 'loss', 6.8, color=MUTED)
F.key_row(None, y+h-27.4, [('line', (INK, None, 'bold'), 'matched PAS'), ('line', (AXC, None, 'bold'), 'unmatched PAS')], center=x+w/2)
F.key_row(None, y+h-17.4, [('block', WIN, 'sequence window')], center=x+w/2)
DEFS = [('PAS', [('Driver', 'matched gain/loss drivers, 0–2'),
                 ('PAS fraction', 'matched / all discovery PASs')]),
        ('Architecture', [('A_g', '2 × matched / (human + mouse PASs)')]),
        ('Sequence', [('phastCons', 'mean, −100 to +50 nt of drivers'),
                      ('Motif rank', 'best hexamer, −35 to −9 nt (0–3)')])]
# DEFS (feature groups and definitions) are given in the caption, not drawn
F.values['a'] = dict(schematic=True, note='no measured coordinates')
F.notes.append('a: schematic only; drivers = PAS with the largest gain and the largest loss of usage (P minus reference); '
               'conserved = reciprocal one-to-one PAS match within 25 nt; model groups: PAS = driver + PAS fraction, Architecture = A_g, Sequence = phastCons + motif rank')

# ================================================================ b  model cohorts
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
cc = D.df('D19A', 'cis/descriptive_class_counts.tsv').set_index(['discovery', 'contrast'])
bx = Lin(0, 1000, x+46.0, x+w-14.0)
Bv = {}
for i, c in enumerate(CO):
    r = cc.loc[c]; yy = y+19.4+i*15.0
    F.text(x+0.4, yy-3.4, CL[c], 7.0)
    y1, n = int(r.model_Y1), int(r.model_n)
    F.rect(bx(0), yy-4.0, bx(y1)-bx(0), 8.0, fc=INK); F.rect(bx(y1), yy-4.0, bx(n)-bx(y1), 8.0, fc=NULLBAR)
    F.text(bx(n)+1.2, yy-3.2, f'{n:,}', 6.8, color=MUTED)
    Bv[CL[c]] = dict(n=n, Y1=y1)
F.xaxis(bx, y+19.4+3*15.0+6.0, [0, 500, 1000], ['0', '500', '1,000'], title='Model genes')
F.key_row(None, y+h-7.4, [('block', INK, 'Y = 1, same direction'), ('block', NULLBAR, 'Y = 0')], gap=5.0, center=x+w/2)
F.values['b'] = Bv

# ================================================================ c  descriptive distributions
x, y, w, h = F.panel('c', X0, Y2, RR-X0, H2)
fd = D.df('D19B', 'cis/feature_distribution_summary.tsv')
CLS = [('Shared', INK, 'o', INK), ('Recurrent-only', MUTED, 'o', 'white'), ('Matched background', NULLBAR, 'o', AXC)]
RNG = {'driver_conserved': (0, 2, [0, 1, 2]), 'frac_conserved_discovery_pas': (0, 1, [0, 0.5, 1]), 'A_g': (0, 1, [0, 0.5, 1]),
       'phastcons_driver_mean': (0, 1, [0, 0.5, 1]), 'motif_rank_driver_mean': (0, 3, [0, 1, 2, 3])}
LABW = 42.0; sg = 8.0
sw = (w-LABW-4*sg-2.0)/5
rp, gg = 6.2, 5.0
r0 = y+19.0
rows = {}
yy = r0
for c in CO:
    for k in range(3):
        rows[(c, k)] = yy; yy += rp
    yy += gg
ybot = yy-gg
Cv = {}
for i, c in enumerate(CO):
    F.text(x+0.4, rows[(c, 1)]-3.4, CL[c], 7.0)
for j, f in enumerate(FEAT):
    lo, hi, ticks = RNG[f]
    sx = Lin(lo, hi, x+LABW+j*(sw+sg)+3.0, x+LABW+j*(sw+sg)+sw-3.0)
    F.text((sx.p0+sx.p1)/2, y+3.0, FL[f], 7.5, ha='center')
    for c in CO:
        F.rect(sx.p0-3.0, rows[(c, 0)]-rp/2, sx.p1-sx.p0+6.0, 3*rp, fc=BAND)
    for c in CO:
        for k, (cls, col, m, fc) in enumerate(CLS):
            q = fd[(fd.discovery == c[0]) & (fd.contrast == c[1]) & (fd.feature == f) & (fd['class'] == cls)]
            assert len(q) == 1, (c, f, cls)
            q = q.iloc[0]; ry = rows[(c, k)]
            F.line([sx(q.q25), sx(q.q75)] if q.q75 > q.q25 else [sx(q.q25)-0.6, sx(q.q75)+0.6], [ry, ry], lw='bold', color=col)
            F.marker(sx(q['median']), ry, m, 'small', fc=fc, ec=INK if cls != 'Matched background' else AXC)
            Cv[f'{f}|{CL[c]}|{cls}'] = dict(n=int(q.n_observed), weight=float(q.weight_sum), q25=float(q.q25), median=float(q['median']), q75=float(q.q75))
    F.xaxis(sx, ybot+0.6, ticks, [num(t, 0) if float(t).is_integer() else f'{t:g}'.replace('0.', '.') for t in ticks])
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'Shared'), ('marker', ('o', 'white', INK), 'Recurrent-only'),
                          ('marker', ('o', AXC), 'Matched background'), ('line', (INK, None, 'bold'), 'weighted IQR'),
                          ('marker', ('o', 'white', INK), 'median (marker)')][:4], gap=7.0, center=x+w/2)
F.values['c'] = Cv
F.notes.append('c: weighted quartiles from the statistics table (match weights; background genes are weighted to the recurrent strata); '
               'marker = weighted median, bar = weighted IQR; descriptive, no test; supported layer, τ = .02')

# ================================================================ d  coefficients
WD = 290.0
x, y, w, h = F.panel('d', X0, Y3, WD, H3)
co = D.df('D19D', 'cis/recomputed_coefficients.tsv')
co = co[(co.layer == PRIMARY['layer']) & (abs(co.tau-PRIMARY['tau']) < 1e-9)]
LW_ = 46.0; dg = 6.0
dw = (w-LW_-3*dg-2.0)/4
rp = 15.0; d0 = y+24.0
Dv = {}
for k, f in enumerate(FEAT):
    F.text(x+0.4, d0+k*rp-3.4, FL[f], 7.0)
for j, c in enumerate(CO):
    sx = Lin(-1, 1, x+LW_+j*(dw+dg)+2.0, x+LW_+j*(dw+dg)+dw-2.0)
    F.text((sx.p0+sx.p1)/2, y+3.0, CL[c], 7.5, ha='center')
    z = co[(co.discovery == c[0]) & (co.contrast == c[1])].set_index('term')
    F.line([sx(0)]*2, [d0-6.0, d0+4*rp+6.0], lw='ref', color=AXC, ls='ref')
    for k, f in enumerate(FEAT):
        r = z.loc[f]; yy = d0+k*rp
        if np.isfinite(r.ci_low):
            F.line([sx.clip(r.ci_low), sx.clip(r.ci_high)], [yy, yy], lw='bold', color=NULLBAR)
        F.marker(sx(r.beta), yy, 'o', 'mark', fc=INK)
        Dv[f'{CL[c]}|{f}'] = dict(beta=float(r.beta), ci=[float(r.ci_low), float(r.ci_high)])
    yb = F.xaxis(sx, d0+4*rp+8.0, [-1, 0, 1], ['−1', '0', '1'])
F.text(x+LW_+(w-LW_)/2, yb+1.6, 'Adjusted logistic coefficient β, All-cis model', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), '95% interval')], gap=5.0, center=x+LW_+(w-LW_)/2)
F.values['d'] = Dv
F.notes.append('d: supported layer, τ = .02; intervals as stored in the statistics table; no interval clipped')

# ================================================================ e  feature correlations
x, y, w, h = F.panel('e', X0+WD+GAP, Y3, RR-(X0+WD+GAP), H3)
fc_ = D.df('D19E', 'cis/feature_correlations.tsv')
pairs = [(a, b) for i, a in enumerate(FEAT) for b in FEAT[i+1:]]
lw_e = max(textwidth(f'{FS[a]} · {FS[b]}', 7.0) for a, b in pairs)
cx0 = x+lw_e+3.0; cw = (x+w-1.0-cx0)/4; rh = 9.2; e0 = y+24.0
F.text(cx0+cw, y+3.0, 'Human', 7.5, ha='center'); F.text(cx0+3*cw, y+3.0, 'Mouse', 7.5, ha='center')
for j, c in enumerate(CO):
    F.text(cx0+cw*(j+0.5), y+13.0, c[1].replace('-', '–'), 7.0, ha='center', color=MUTED)
F.line([cx0+0.6, cx0+2*cw-0.6], [y+11.6]*2, lw='thin', color=INK); F.line([cx0+2*cw+0.6, cx0+4*cw-0.6], [y+11.6]*2, lw='thin', color=INK)
Ev = {}
for i, (a, b) in enumerate(pairs):
    yy = e0+i*rh
    F.text(cx0-2.0, yy+rh/2-3.4, f'{FS[a]} · {FS[b]}', 7.0, ha='right')
    for j, c in enumerate(CO):
        q = fc_[(fc_.discovery == c[0]) & (fc_.contrast == c[1]) & (fc_.feature1 == a) & (fc_.feature2 == b)]
        assert len(q) == 1
        v = float(q.spearman.iloc[0]); col = diverge(v, 1.0)
        F.rect(cx0+cw*j+0.4, yy+0.4, cw-0.8, rh-0.8, fc=col)
        F.text(cx0+cw*(j+0.5), yy+rh/2-3.2, f'{v:.2f}'.replace('0.', '.').replace('-', '−'), 6.8, ha='center', color=on(col))
        Ev[f'{a}|{b}|{CL[c]}'] = dict(rho=v, n=int(q.n_pairwise.iloc[0]))
F.cbar_key(x+w/2, y+h-7.4, 'Spearman ρ', '−1', '1', lambda f_: diverge(-1+2*f_, 1.0), wbar=40)
F.values['e'] = Ev
F.save(OUT, D)
