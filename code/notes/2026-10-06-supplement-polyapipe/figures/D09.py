"""Supplementary Figure D09 (printed number from supdraw.SNUM) — within-species recurrence and its robustness (R3),
polyApipe v2 primary.

a  recurrence enrichment (observed / null mean) under the five measurement-null models, two effect thresholds
b  proximal / distal composition of the recurrent genes
c  same-state pseudo-comparisons (pseudo-recurrent genes) against the observed recurrence on the same background
d  recurrent members kept, same direction, under six measurement variants
e  leave-one-study-out (human P–M, the only contrast with >= 2 remaining studies)
f  Guan series treated as separate pseudo-units (sensitivity only; human P–S)
All panels use the WPP metric.
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D09'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 204.0, 92.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
GROUPS = [('human', 'P_S', 'Human P–S'), ('human', 'P_M', 'Human P–M'), ('mouse', 'P_S', 'Mouse P–S'), ('mouse', 'P_M', 'Mouse P–M')]
TH = [(0.02, 'o', '|ΔDWUI| ≥ .02'), (0.05, 'D', '≥ .05')]

# ================================================================ a  null ladder
WA = 250.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
lad = D.df('D09A', 'recurrence_controls/recurrence_ladder.tsv')
MODELS = [('L0_feasible_series', 'L0 feasible series'), ('L1_measurement', 'L1 measurement'), ('L2_baseline', 'L2 baseline'),
          ('L2_PAS_first', 'L2 PAS-first'), ('L3_full', 'L3 full')]
labw = max(textwidth(l, 7.0) for _, l in MODELS)
sx = Lin(0.5, 4.2, x+labw+56, x+w-4)
F.text((sx.p0+sx.p1)/2, y+3.0, 'Recurrent genes, observed / null mean', 7.5, ha='center')
rh = 7.6; Av = {}; yy = y+20.0
F.line([sx(1.0)]*2, [yy-5, yy+rh*20+3*4.0-2], lw='ref', color=AXC, ls='ref')
for gi, (sp, ct, gl) in enumerate(GROUPS):
    F.text(x+0.4, yy+rh*2-3.4, gl, 7.5, font='Bol')
    if gi % 2 == 0: F.rect(x, yy-rh/2, w, rh*5, fc=BAND)
    for mk_, ml in MODELS:
        F.text(sx.p0-4.0, yy-3.4, ml, 7.0, ha='right')
        for th, mshape, _ in TH:
            q = lad[(lad.species == sp) & (lad.contrast == ct) & (lad.model == mk_) & (abs(lad.threshold-th) < 1e-9)].iloc[0]
            lo, hi = (q.expected-2*q.null_sd)/q.expected, (q.expected+2*q.null_sd)/q.expected
            F.line([sx(lo), sx(hi)], [yy+(-1.4 if th == 0.02 else 1.4)]*2, lw='thin', color=NULLBAR)
            F.marker(sx.clip(q.observed/q.expected), yy, mshape, 'small', fc=INK if th == 0.02 else 'white', ec=INK)
            Av[f'{gl}|{mk_}|{th}'] = dict(observed=int(q.observed), expected=float(q.expected), null_sd=float(q.null_sd),
                                           ratio=float(q.observed/q.expected), q_all_40=float(q.q_all_40), n_genes=int(q.n_genes))
        yy += rh
    yy += 4.0
ab = F.xaxis(sx, yy-rh/2+1.0, [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4], ['.5', '1', '1.5', '2', '2.5', '3', '3.5', '4'])
F.key_row(None, y+h-7.4, [('marker', ('o', INK), TH[0][2]), ('marker', ('D', 'white', INK), TH[1][2]),
                          ('line', (NULLBAR, None, 'thin'), 'Null mean ± 2 SD')], gap=7.0, center=x+w/2)
F.values['a'] = Av
F.notes.append('a: null band = (null mean ± 2 SD) / null mean from the recorded ladder null SD (display); all 40 configurations '
               'have q_all_40 < .001')

# ================================================================ b  direction composition (form of main Figure 3a, right)
xb = X0+WA+GAP; wb = RR-xb
H_B = 96.0
x, y, w, h = F.panel('b', xb, Y1, wb, H_B)
dc = D.df('D09B', 'panel_statistics/recurrence-direction-counts.tsv')
dc = dc[dc.metric == 'WPP']
b0 = x+48.0; bw_ = w-48.0-38.0
F.text(b0+bw_/2, y+3.0, 'Distal | proximal recurrent genes', 7.5, ha='center')
F.text(b0+bw_+19.0, y+3.0, '≥ .02', 7.0, ha='center', color=MUTED)
Bv = {}; yy = y+20.0
for sp, ct, gl in GROUPS:
    F.text(x+0.4, yy+2.8-3.4, gl, 7.0)
    for k, th in enumerate((0.02, 0.05)):
        q = dc[(dc.species == sp) & (dc.contrast == ct) & (abs(dc.threshold-th) < 1e-9)].set_index('direction').genes
        prox, dist = int(q['negative']), int(q['positive'])
        tot = max(prox+dist, 1); yk = yy+k*5.6; bh = 4.4 if th == 0.02 else 3.0
        F.rect(b0, yk-bh/2, bw_*dist/tot, bh, fc=DIST); F.rect(b0+bw_*dist/tot, yk-bh/2, bw_*prox/tot, bh, fc=PROX)
        if th == 0.02: F.text(b0+bw_+19.0, yk-3.0, f'{dist} | {prox}', 6.8, ha='center')
        Bv[f'{gl}|{th}'] = dict(proximal=prox, distal=dist)
    yy += 14.4
F.key_row(None, y+H_B-7.4, [('block', DIST, 'Distal'), ('block', PROX, 'Proximal')], gap=7.0,
          center=b0+bw_/2)
F.values['b'] = Bv
F.notes.append('b: form of main Figure 3a (right): bar = distal share (left, orange) and proximal share (right, blue) of the recurrent genes; '
               'upper bar |ΔDWUI| >= .02 with counts, lower thin bar >= .05')

# ================================================================ c  same-state pseudo-comparisons
x, y, w, h = F.panel('c', xb, Y1+H_B+GAP, wb, H1-H_B-GAP)
ss = D.df('D09C', 'measurement/same-state/recurrence_summary.tsv')
ss = ss[(ss.metric == 'WPP') & (abs(ss.threshold-0.02) < 1e-9)]
lx = Lin(0, 3, x+44, x+w-10.5)
F.text((lx.p0+lx.p1)/2, y+3.0, 'Recurrent genes, |ΔDWUI| ≥ .02 (log)', 7.5, ha='center')
Cv = {}; yy = y+19.0
for sp, ct, gl in GROUPS:
    F.text(x+0.4, yy+2.4-3.4, gl, 7.0)
    for k, kind in enumerate(('target', 'reference')):
        q = ss[(ss.species == sp) & (ss.contrast == ct) & (ss.pseudo_kind == kind)].iloc[0]
        yk = yy+k*5.0
        F.line([lx(math.log10(max(q.pseudo_median, 1))), lx(math.log10(max(q.pseudo_p95, 1)))], [yk, yk], lw='bold', color=NULLBAR)
        F.marker(lx(math.log10(max(q.pseudo_median, 1))), yk, 'o', 'small', fc='white', ec=MUTED)
        F.marker(lx(math.log10(q.observed_same_background_median)), yk, 'o', 'small', fc=INK)
        Cv[f'{gl}|{kind}'] = dict(pseudo_median=float(q.pseudo_median), pseudo_p95=float(q.pseudo_p95),
                                  observed=int(q.observed_same_background_median), eligible=int(q.eligible_median))
    yy += 15.0
F.xaxis(lx, yy-4.0, [0, 1, 2, 3], ['1', '10', '100', '1,000'])
F.key_row(None, y+h-7.4, [('marker', ('o', 'white', MUTED), 'Same-state halves'), ('marker', ('o', INK), 'Observed')], gap=6.0,
          center=(lx.p0+lx.p1)/2)
F.values['c'] = Cv
F.notes.append('c: two rows per contrast = target-state and reference-state pseudo-comparisons; values < 1 drawn at 1')

# ================================================================ d  member retention (form of main Figure 3b, members retained)
WD = 224.0
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
mr = D.df('D09D', 'measurement/nulls/member_retention.tsv')
mr = mr[(mr.metric == 'WPP') & (abs(mr.threshold-0.02) < 1e-9)]
VAR = [('native', 'Native'), ('interval_supported', 'Bootstrap interval'), ('depth20', '20 PAS UMIs'), ('depth40', '40 PAS UMIs'),
       ('cycle_matched', 'Cycle matched'), ('common', 'Common PAS axis')]
labw = max(textwidth(vl, 7.0) for _, vl in VAR)
cx0 = x+labw+5.0; cw = (x+w-cx0)/len(GROUPS); rh = 10.6; r0 = y+30.0; bwd = 12.0
F.text(cx0+cw*2, y+3.0, 'Members kept, same direction', 7.5, ha='center')
Dv = {}
for j_, (sp, ct, gl) in enumerate(GROUPS):
    sp_, ct_ = gl.split(' ')                                          # two-line header: species over contrast
    F.text(cx0+cw*j_+1.0, y+11.0, sp_, 7.0); F.text(cx0+cw*j_+1.0, y+18.6, ct_, 7.0, color=MUTED)
for i, (v, vl) in enumerate(VAR):
    yy = r0+rh*i
    if i % 2 == 0: F.rect(x+labw+2.0, yy-rh/2, w-labw-2.0, rh, fc=BAND)
    F.text(x+labw, yy-3.4, vl, 7.0, ha='right')
    for j_, (sp, ct, gl) in enumerate(GROUPS):
        q = mr[(mr.species == sp) & (mr.contrast == ct) & (mr.variant == v)].iloc[0]
        bx_ = cx0+cw*j_+1.0
        if q.status == 'tested':
            kept, orig = int(q.same_direction_intersection), int(q.original)
            F.rect(bx_, yy-2.2, bwd, 4.4, fc=TRACK); F.rect(bx_, yy-2.2, bwd*kept/orig, 4.4, fc=INK)
            F.text(bx_+bwd+1.6, yy-3.2, f'{kept}/{orig}', 6.8)
            Dv[f'{gl}|{v}'] = dict(kept=kept, original=orig, current=int(q.current))
        else:
            F.text(bx_+bwd/2, yy-3.2, '–', 6.8, ha='center', color=AXC)
            Dv[f'{gl}|{v}'] = None
F.values['d'] = Dv
F.notes.append('d: – = fewer than two studies under that variant (unavailable); |ΔDWUI| ≥ .02 (WPP in the source tables = DWUI)')

# ================================================================ e  leave one study out
WE = 128.0
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, WE, H2)
lo = D.df('D09E', 'recurrence/leave_one_study_out_nulls.tsv')
lo = lo[(lo.metric == 'WPP') & (~lo.series_matching) & (lo.species == 'human') & (lo.contrast == 'P_M') & (lo.status == 'tested')]
ex = Lin(0, 180, x+44, x+w-7)
F.text((ex.p0+ex.p1)/2, y+3.0, 'Human P–M', 7.5, ha='center')
Ev = {}; yy = y+22.0
for st in ['GSE118258', 'GSE150311', 'GSE221739']:
    F.rect(x, yy-6.0, w, 17.6, fc=BAND) if st != 'GSE150311' else None
    F.text(x+0.4, yy+2.6-3.4, 'w/o '+DATASET[st][0].replace('2020', '').replace('2023', ''), 7.0)
    for k, (th, mshape, _) in enumerate(TH):
        q = lo[(lo.omitted_study == st) & (abs(lo.threshold-th) < 1e-9)].iloc[0]
        yk = yy+k*5.6
        F.line([ex(q.null_p025), ex(q.null_p975)], [yk, yk], lw='bold', color=NULLBAR)
        F.marker(ex(q.observed), yk, mshape, 'small', fc=INK if th == 0.02 else 'white', ec=INK)
        Ev[f'{st}|{th}'] = dict(observed=int(q.observed), null=[float(q.null_p025), float(q.null_p975)], q_all_12=float(q.q_all_12))
    yy += 19.0
F.xaxis(ex, yy-8.0, [0, 60, 120, 180])
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'Null 95% range')], gap=6.0, center=(ex.p0+ex.p1)/2)
F.values['e'] = Ev
F.notes.append('e: P–S and both mouse contrasts keep a single study when one is omitted (not estimable); Guan2022 not evaluable for P–M')

# ================================================================ f  Guan pseudo-units
x, y, w, h = F.panel('f', X0+WD+GAP+WE+GAP, Y2, RR-(X0+WD+GAP+WE+GAP), H2)
gu = D.df('D09F', 'recurrence_controls/guan_series_unit_sensitivity.tsv')
gu = gu[(gu.species == 'human') & (gu.contrast == 'P_S')]
fx = Lin(0, 1400, x+20, x+w-10.5)
F.text((fx.p0+fx.p1)/2+4, y+3.0, 'Guan series', 7.5, ha='center')
Fv = {}; yy = y+28.0
for k, (th, mshape, lab) in enumerate(TH):
    q = gu[abs(gu.threshold-th) < 1e-9].iloc[0]
    yk = yy+k*24.0
    F.rect(x, yk-12.0, w, 24.0, fc=BAND) if k == 0 else None
    F.text(x+0.4, yk-3.4, ['.02', '.05'][k], 7.0)
    F.line([fx(q.null_p025), fx(q.null_p975)], [yk, yk], lw='bold', color=NULLBAR)
    F.marker(fx(q.observed), yk, mshape, 'small', fc=INK if th == 0.02 else 'white', ec=INK)
    Fv[str(th)] = dict(observed=int(q.observed), null=[float(q.null_p025), float(q.null_p975)], enrichment=float(q.enrichment))
F.xaxis(fx, yy+30.0, [0, 700, 1400], ['0', '700', '1,400'])
F.text((fx.p0+fx.p1)/2, y+h-17.4, 'Human P–S', 7.0, ha='center', color=MUTED)
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'Null 95%')], gap=6.0, center=(fx.p0+fx.p1)/2)
F.values['f'] = Fv
F.save(OUT, D)
