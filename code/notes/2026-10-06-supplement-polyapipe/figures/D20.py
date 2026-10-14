"""Supplementary Figure D20 (printed number from supdraw.SNUM) — predictive increment of the cis features (R4), polyApipe v2
primary. Takes over the four pre-specified primary increment tests from the former main Figure 4f.

a  the four pre-specified primary tests (supported layer, τ = .02): likelihood-ratio statistic of All cis over M1 against its
   conditional parametric-bootstrap null; bootstrap p and BH q within the four primary tests
b  the same All-cis-over-M1 test in all 32 planned settings (2 layers x 4 thresholds x 4 cohorts): q within the 32
c  out-of-fold gains of each model step at the primary setting (M1 over M0; PAS, Architecture, Sequence and All cis each over
   M1), mean of 10 paired repeats with conditional gene-bootstrap 95% intervals; log loss shown as reduction
d  the ten repeat-level All-cis-over-M1 gains, AUC against log-loss reduction, per cohort
e  All-cis-over-M1 repeat-mean gains across thresholds and PAS layers (all 32 settings)
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D20'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 96.0, 200.0, 168.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname(FID))
CO = [('human', 'P-S'), ('human', 'P-M'), ('mouse', 'P-S'), ('mouse', 'P-M')]
CL = {c: f'{c[0].capitalize()} {c[1].replace("-", "–")}' for c in CO}
CC = {c: (SC['S'] if c[1] == 'P-S' else SC['M']) for c in CO}          # contrast = colour of the reference state
CM = {c: ('o' if c[0] == 'human' else 's') for c in CO}                # species = marker
TAUS = [0.01, 0.02, 0.03, 0.05]
md = D.df('D20A', 'cis/recomputed_models.tsv')


def ctx_mark(xm, yy, c, ms='mark'):
    F.marker(xm, yy, CM[c], ms, fc=CC[c])


# ================================================================ a  primary tests
WA = 262.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
nl = D.df('D20A', 'cis/conditional_bootstrap_nulls.tsv.gz')
lx = Lin(0, 16, x+52.0, x+w-88.0)
CX = [x+w-64.0, x+w-29.0, x+w-8.5]
# 2026-10-12: panel title (LRT, All cis vs M1) moved to the caption
for cx_, lab in zip(CX, ['Y = 1 / n', 'p', 'q']):
    F.text(cx_, y+3.0, lab, 7.0, ha='center', color=MUTED)
Av = {}
for i, c in enumerate(CO):
    r = md[(md.layer == 'supported50') & (abs(md.tau-0.02) < 1e-9) & (md.discovery == c[0]) & (md.contrast == c[1])].iloc[0]
    assert bool(r.primary)
    z = nl[(nl.layer == r.layer) & (nl.contrast == r.contrast) & (nl.discovery == r.discovery) & (abs(nl.tau-r.tau) < 1e-9)].lrt_null.to_numpy()
    assert len(z) == r.bootstrap_draws
    assert abs((1+(z >= r.lrt).sum())/(1+len(z)) - r.p_boot) < 1e-9
    q025, q50, q975 = np.percentile(z, [2.5, 50, 97.5])
    yy = y+20.0+i*13.5
    ctx_mark(x+3.0, yy, c); F.text(x+8.4, yy-3.4, CL[c], 7.0)
    F.line([lx(q025), lx.clip(q975)], [yy, yy], lw='range', color=NULLBAR)
    F.line([lx(q50)]*2, [yy-3.0, yy+3.0], lw='bold', color=NULLW)
    sig = r.q_boot_primary_four < 0.05
    F.marker(lx.clip(r.lrt), yy, 'o', 'mark', fc=INK if sig else 'white', ec=INK)
    F.text(CX[0], yy-3.4, f'{int(r.positive)} / {int(r.n):,}', 7.0, ha='center')
    F.text(CX[1], yy-3.4, qfmt(r.p_boot), 7.0, ha='center')
    F.text(CX[2], yy-3.4, qfmt(r.q_boot_primary_four), 7.0, ha='center', font='Bol' if sig else 'Reg')
    Av[CL[c]] = dict(lrt=float(r.lrt), df=int(r.df), null_q025=float(q025), null_median=float(q50), null_q975=float(q975),
                     p_boot=float(r.p_boot), q_primary_four=float(r.q_boot_primary_four), n=int(r.n), Y1=int(r.positive), draws=int(len(z)))
F.xaxis(lx, y+20.0+3*13.5+7.0, [0, 4, 8, 12, 16], title='Likelihood-ratio statistic (5 df)')
F.key_row(None, y+h-7.4, [('marker', ('o', 'white', INK), 'Observed, q ≥ .05'), ('line', (NULLBAR, None, 'range'), 'Bootstrap null 95%'),
                          ('line', (NULLW, None, 'bold'), 'null median')], gap=6.0, center=x+w/2)
F.values['a'] = Av
F.notes.append('a: null range and median are display summaries of the stored 4,999 conditional bootstrap LRTs per test; the stored p_boot is reproduced exactly from them')

# ================================================================ b  32 planned settings
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
cx0 = x+50.0; cw = (x+w-1.0-cx0)/8; rh = 12.4; b0 = y+24.0
F.text(cx0+2*cw, y+3.0, 'Supported', 7.5, ha='center'); F.text(cx0+6*cw, y+3.0, 'Expanded', 7.5, ha='center')
for k in range(8):
    F.text(cx0+cw*(k+0.5), y+13.4, f'{TAUS[k % 4]:g}'.replace('0.', '.'), 7.0, ha='center', color=MUTED)
F.line([cx0+0.6, cx0+4*cw-0.6], [y+12.0]*2, lw='thin', color=INK); F.line([cx0+4*cw+0.6, cx0+8*cw-0.6], [y+12.0]*2, lw='thin', color=INK)
Bv = {}
qcol = lambda q: seq(min(-np.log10(q), 1.0), 0.0, 1.0)
for i, c in enumerate(CO):
    yy = b0+i*rh
    ctx_mark(x+3.0, yy+rh/2, c, 'small'); F.text(x+7.4, yy+rh/2-3.4, CL[c], 7.0)
    for k, (layer, tau) in enumerate([(l_, t_) for l_ in ('supported50', 'expanded_primary') for t_ in TAUS]):
        r = md[(md.layer == layer) & (abs(md.tau-tau) < 1e-9) & (md.discovery == c[0]) & (md.contrast == c[1])].iloc[0]
        fc = qcol(r.q_boot_planned_32)
        prim = layer == 'supported50' and abs(tau-0.02) < 1e-9
        F.rect(cx0+cw*k+0.4, yy+0.4, cw-0.8, rh-0.8, fc=fc, ec=INK if prim else 'none', lw='thin' if prim else 0)
        F.text(cx0+cw*(k+0.5), yy+rh/2-3.2, qfmt(r.q_boot_planned_32), 6.8, ha='center', color=on(fc))
        Bv[f'{CL[c]}|{layer}|{tau}'] = float(r.q_boot_planned_32)
F.text(cx0-2.0, y+13.4, 'τ', 7.0, ha='right', color=MUTED)   # outlined-cell and 32-test notes: caption
F.cbar_key(cx0+4*cw, y+h-7.4, 'q', '1', '≤ .1', lambda f_: seq(f_), wbar=40)
F.values['b'] = Bv

# ================================================================ c  paired CV gains, primary setting
WC = 292.0
x, y, w, h = F.panel('c', X0, Y2, WC, H2)
pu = D.df('D20B', 'cis/paired_cv_uncertainty.tsv')
pu = pu[(pu.layer == 'supported50') & (abs(pu.tau-0.02) < 1e-9) & (pu.model != 'M1')]
MODELS = [('M0', 'M1 over M0'), ('PAS', 'PAS'), ('Architecture', 'Architecture'), ('Sequence', 'Sequence'), ('All cis', 'All cis')]
LBW = 89.0; sg = 8.0
sw = (w-LBW-sg-2.0)/2
rp = 6.8; c0 = y+21.0
rows = {}; yy = c0
for c in CO:
    for m, _ in MODELS:
        rows[(c, m)] = yy; yy += rp
    yy += 4.0
ybot = yy-4.0
for c in CO:
    ctx_mark(x+3.0, rows[(c, 'M0')], c); F.text(x+8.4, rows[(c, 'M0')]-3.4, CL[c], 7.0)
    for m, lab in MODELS:
        F.text(x+LBW-2.0, rows[(c, m)]-3.2, lab, 6.8, ha='right', color=MUTED if m == 'M0' else INK)
Cv = {}
for j, (metric, title, lim, ticks, tl) in enumerate((('auc', 'AUC gain', (-0.08, 0.1), [-0.05, 0, 0.05], ['−.05', '0', '.05']),
                                                      ('logloss', 'Log-loss reduction', (-0.02, 0.03), [-0.02, 0, 0.02], ['−.02', '0', '.02']))):
    sx = Lin(lim[0], lim[1], x+LBW+j*(sw+sg)+2.0, x+LBW+j*(sw+sg)+sw-2.0)
    F.text((sx.p0+sx.p1)/2, y+3.0, title, 7.5, ha='center')
    for c in CO:
        F.rect(sx.p0-2.0, rows[(c, 'M0')]-rp/2, sx.p1-sx.p0+4.0, 5*rp, fc=BAND)
    F.line([sx(0)]*2, [c0-5.0, ybot+1.0], lw='ref', color=AXC, ls='ref')
    for c in CO:
        for m, lab in MODELS:
            r = pu[(pu.discovery == c[0]) & (pu.contrast == c[1]) & (pu.model == m) & (pu.metric == metric)].iloc[0]
            rev = (m == 'M0') != (metric == 'logloss')        # show M1 over M0, and log loss as a reduction
            d_, lo, hi = (-r.delta, -r.hi, -r.lo) if rev else (r.delta, r.lo, r.hi)
            ry = rows[(c, m)]
            F.line([sx.clip(lo), sx.clip(hi)], [ry, ry], lw='bold', color=NULLBAR)
            F.marker(sx(d_), ry, CM[c], 'small', fc=CC[c] if m == 'All cis' else INK)
            Cv[f'{metric}|{CL[c]}|{m}'] = dict(gain=float(d_), ci=[float(lo), float(hi)], repeats_mean=True, draws=int(r.draws))
    ybx = F.xaxis(sx, ybot+0.6, ticks, tl)
F.text(x+LBW+(w-LBW)/2, ybx+1.6, 'Higher = better; primary setting (supported layer, τ = .02)', 7.0, ha='center', color=MUTED)
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'Gene-bootstrap 95%')], gap=5.0, center=x+LBW+(w-LBW)/2)
F.values['c'] = Cv
F.notes.append('c: each step except M1-over-M0 is an addition to M1 (gains are not additive); log-loss signs and interval ends reversed to show reduction; '
               'mean of 10 fixed-fold repeats, intervals from the stored conditional gene bootstrap (2,000 draws)')

# ================================================================ d  repeat-level gains
x, y, w, h = F.panel('d', X0+WC+GAP, Y2, RR-(X0+WC+GAP), H2)
cv = D.df('D20C', 'cis/recomputed_cv.tsv')
p = cv[(cv.layer == 'supported50') & (abs(cv.tau-0.02) < 1e-9) & cv.model.isin(['M1', 'All cis'])]
wv = p.pivot_table(index=['discovery', 'contrast', 'repeat'], columns='model', values=['auc', 'logloss'])
g = (wv[('auc', 'All cis')]-wv[('auc', 'M1')]).rename('auc').to_frame()
g['ll'] = wv[('logloss', 'M1')]-wv[('logloss', 'All cis')]
g = g.reset_index()
side = (w-30.0-6.0-1.0)/2
sh_d = (h-16.0-30.0-24.0)/2
Dv = {}
for k, c in enumerate(CO):
    col, row = k % 2, k//2
    sx = Lin(-0.05, 0.05, x+30.0+col*(side+6.0), x+30.0+col*(side+6.0)+side)
    top = y+16.0+row*(sh_d+30.0)
    sy = Lin(-0.02, 0.01, top+sh_d, top)
    F.text((sx.p0+sx.p1)/2, top-13.0, CL[c], 7.5, ha='center')
    F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref'); F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    z = g[(g.discovery == c[0]) & (g.contrast == c[1])]
    assert len(z) == 10
    for a_, b_ in zip(z.auc, z.ll):
        F.marker(sx.clip(a_), sy.clip(b_), CM[c], 'small', fc=CC[c])
    F.yaxis(sy, sx.p0, [-0.02, -0.01, 0, 0.01], ['−.02', '−.01', '0', '.01'], 'Log-loss reduction' if col == 0 and row == 0 else None,
            labels_off=(col == 1), title_y=y+16.0+sh_d+15.0)
    F.xaxis(sx, sy.p0+1.0, [-0.04, 0, 0.04], ['−.04', '0', '.04'])
    Dv[CL[c]] = dict(auc_gain=[float(v) for v in z.auc], logloss_reduction=[float(v) for v in z.ll])
F.text(x+30.0+side+3.0, y+h-8.4, 'AUC gain, All cis over M1', 8.0, ha='center')
F.values['d'] = Dv
F.notes.append('d: repeat-level differences of the stored out-of-fold metrics (display transform); the ten repeats reuse one gene draw and are not independent')

# ================================================================ e  32 settings
x, y, w, h = F.panel('e', X0, Y3, RR-X0, H3)
cv2 = cv[cv.model.isin(['M1', 'All cis'])]
gm = cv2.groupby(['layer', 'discovery', 'contrast', 'tau', 'model'])[['auc', 'logloss']].mean().unstack('model')
gain = pd.DataFrame({'auc': gm[('auc', 'All cis')]-gm[('auc', 'M1')], 'logloss': gm[('logloss', 'M1')]-gm[('logloss', 'All cis')]}).reset_index()
LW_ = 30.0; sg = 8.0
sw = (w-LW_-3*sg-1.0)/4
sh = sw/1.8
Ev = {}
for row, (metric, lab, lim, ticks, tl) in enumerate((('auc', 'AUC gain', (-0.06, 0.08), [-0.04, 0, 0.04, 0.08], ['−.04', '0', '.04', '.08']),
                                                    ('logloss', 'Log-loss red.', (-0.04, 0.01), [-0.04, -0.02, 0], ['−.04', '−.02', '0']))):
    top = y+16.0+row*(sh+16.0)
    for k, c in enumerate(CO):
        sx = Lin(0.0, 0.06, x+LW_+k*(sw+sg), x+LW_+k*(sw+sg)+sw)
        sy = Lin(lim[0], lim[1], top+sh, top)
        F.lineplot(sx, sy, f'{metric}|{CL[c]}')
        if row == 0: F.text((sx.p0+sx.p1)/2, y+3.0, CL[c], 7.5, ha='center')
        F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
        for layer, ls in (('supported50', None), ('expanded_primary', 'dash')):
            z = gain[(gain.layer == layer) & (gain.discovery == c[0]) & (gain.contrast == c[1])].sort_values('tau')
            assert len(z) == 4
            F.line([sx(t) for t in z.tau], [sy(v) for v in z[metric]], lw='data', color=CC[c], ls=ls)
            for t, v in zip(z.tau, z[metric]):
                prim = layer == 'supported50' and abs(t-0.02) < 1e-9
                F.marker(sx(t), sy(v), CM[c], 'mark' if prim else 'small', fc=CC[c] if (layer == 'supported50') else 'white', ec=CC[c])
            Ev[f'{metric}|{CL[c]}|{layer}'] = {f'{t:g}': float(v) for t, v in zip(z.tau, z[metric])}
        F.yaxis(sy, sx.p0, ticks, tl, lab if k == 0 else None, labels_off=(k > 0))
        yb = F.xaxis(sx, sy.p0+1.0, [0.01, 0.02, 0.03, 0.05], ['.01', '.02', '.03', '.05'] if row == 1 else None, line=True) if row == 1 else \
            F.xaxis(sx, sy.p0+1.0, [0.01, 0.02, 0.03, 0.05], ['', '', '', ''])
F.text(x+LW_+(w-LW_)/2, yb+1.6, 'Effect threshold τ', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('line', (INK, None, 'data'), 'Supported layer'), ('line', (INK, 'dash', 'data'), 'Expanded layer'),
                          ('marker', ('o', INK), 'primary setting (larger)')], gap=7.0, center=x+LW_+(w-LW_)/2)
F.values['e'] = Ev
F.notes.append('e: repeat means of the stored out-of-fold metrics (display summary); no q values beyond b')
F.save(OUT, D)
