"""Supplementary Figure D13 (printed number from supdraw.SNUM) — intermediate-state exploration (R2), polyApipe v2 primary.

a  recurrent candidates -> endpoint measurable -> endpoint retained, then parallel R and E eligibility (human E unavailable)
b  metacells passing the frozen expression-distance support gate (q90, q95) per study and node
c  process classes of the R / E change relative to the P–S endpoint in both studies (primary: shared p08, requalified)
d  aligned intermediate effects, study 1 against study 2, coloured by class
e  nine settings: endpoint-aligned share after requalification, and class stability in the fixed primary cohort
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D13'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 112.0, 140.0, 104.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname(FID))
P = 'tables/intermediate/process/'
cf = D.df('D13A', P+'coverage_flow.tsv')
ss = D.df('D13C', P+'same_series_process_classes.tsv.gz')
CLS = [('endpoint_aligned', 'Aligned', SEQ_STOPS[2]), ('reverse', 'Reverse', DIST), ('within_threshold', 'Small', NULLBAR),
       ('direction_difference', 'Disagreement', DISAGREE), ('mixed_or_threshold_undetermined', 'Mixed', SEQ_STOPS[1]),
       ('endpoint_not_retained', 'Not retained', MUTED), ('unavailable', 'Unavailable', 'white')]
CC = {k: c for k, _, c in CLS}
GROUPS = [('human', 'R', 'Human R'), ('mouse', 'R', 'Mouse R'), ('mouse', 'E', 'Mouse E')]
prim = ss[(ss.variant == 'shared_p08') & (ss.scope == 'requalified') & (abs(ss.intermediate_threshold-0.02) < 1e-9)]

# ================================================================ a  eligibility flow
WA = 140.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
ax_ = Lin(0, 650, x+46, x+w-8)
F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Genes', 7.5, ha='center')
Av = {}; yy = y+20.0
for sp, lab, col in (('mouse', 'Mouse', STUDY['CIP'][1]), ('human', 'Human', STUDY['GSE118258'][1])):
    q = cf[(cf.species == sp) & (cf.variant == 'shared_p08') & (cf.scope == 'requalified')].iloc[0]
    F.text(x+0.4, yy-3.4, lab, 7.5, font='Bol')
    steps = [('original_candidates', 'Recurrent'), ('endpoint_measurable', 'Measurable'), ('endpoint_retained', 'Retained'),
             ('R', 'R eligible'), ('E', 'E eligible')]
    for k, (c_, cl) in enumerate(steps):
        yk = yy+7.6+k*7.0
        F.text(ax_.p0-3.0, yk-3.4, cl, 7.0, ha='right', color=INK if k < 3 else SC[c_[0]])
        v = q[c_]
        if v != v:
            F.text(ax_(0)+1.0, yk-3.4, 'not available', 7.0, color=MUTED); continue
        F.rect(ax_(0), yk-2.3, ax_(v)-ax_(0), 4.6, fc=NULLBAR if k < 3 else SC[c_[0]])
        F.text(ax_(v)+2.0, yk-3.4, str(int(v)), 6.8)
    Av[sp] = {c_: (None if q[c_] != q[c_] else int(q[c_])) for c_, _ in steps}
    yy += 44.0
F.xaxis(ax_, yy-7.0, [0, 300, 600])
F.values['a'] = Av

# ================================================================ b  expression support
WB = 160.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
es = D.df('D13B', P+'frozen_expression_support_summary.tsv')
NODES = ['S', 'R', 'E', 'P', 'M']
STD = [('GSE118258', 'Xing2020'), ('GSE178325', 'Guan2022'), ('GSE103221', 'Guo2019'), ('GSE115943', 'Schiebinger2019')]
labw = max(textwidth(n, 7.0) for _, n in STD)
hx = x+labw+22.0; cw = (x+w-hx)/5; rh = 8.6; r0 = y+14.0
F.text(hx+cw*2.5, y+3.0, 'Supported metacells', 7.5, ha='center')
for k, nd in enumerate(NODES):
    F.text(hx+cw*(k+0.5), r0-0.2-8.0+0.6, nd, 7.0, ha='center', color=SC[nd]) if False else None
Bv = {}
for i, (st, nm) in enumerate(STD):
    for g_, (col_, gl) in enumerate((('support_q90', 'q90'), ('support_q95', 'q95'))):
        yy = r0+rh*(2*i+g_)+2.0*i
        if g_ == 0: F.text(x+0.4, yy+rh-3.4, nm, 7.0)
        F.text(hx-2.0, yy+rh/2-3.2, gl, 6.8, ha='right', color=MUTED)
        for k, nd in enumerate(NODES):
            q = es[(es.study_id == st) & (es.node == nd)]
            if len(q) == 0:
                F.rect(hx+cw*k+0.4, yy+0.4, cw-0.8, rh-0.8, fc=BAND); continue
            q = q.iloc[0]; fr = q[col_]/q.metacells; fc = seq(fr, 0, 1)
            F.rect(hx+cw*k+0.4, yy+0.4, cw-0.8, rh-0.8, fc=fc)
            F.text(hx+cw*(k+0.5), yy+rh/2-3.2, str(int(q[col_])), 6.8, ha='center', color=on(fc))
            Bv[f'{st}|{nd}|{col_}'] = dict(supported=int(q[col_]), metacells=int(q.metacells))
yb = r0+rh*8+6.0
for k, nd in enumerate(NODES):
    F.text(hx+cw*(k+0.5), yb+0.4, nd, 7.5, ha='center', color=SC[nd])
F.cbar_key(x+w/2, y+h-7.4, 'Share of metacells', '0', '1', lambda f_: seq(f_, 0, 1), wbar=44)
F.values['b'] = Bv

# ================================================================ c  process classes
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
cx = Lin(0, 1, x+36, x+w-2)
F.text(x+w-48.5, y+3.0, 'Share of original P–S genes', 7.5, ha='center')
Cv = {}
for i, (sp, stt, gl) in enumerate(GROUPS):
    yy = y+23.0+i*19.5
    q = prim[(prim.species == sp) & (prim.state == stt)]
    tot = len(q); acc = 0; Cv[gl] = {}
    F.text(x+0.4, yy-3.4, gl, 7.0)
    for k, cl, cc in CLS:
        n = int((q.pair_class == k).sum())
        F.rect(cx(acc/tot), yy-5.8, cx((acc+n)/tot)-cx(acc/tot), 11.6, fc=cc, ec=AXC if cc == 'white' else 'white', lw='axis')
        acc += n; Cv[gl][k] = n
    Cv[gl]['total'] = tot
F.xaxis(cx, y+23.0+2*19.5+7.0, [0, 0.5, 1], ['0', '.5', '1'])
keys = [('block', cc, cl) for _, cl, cc in CLS]
F.key_row(None, y+h-27.4, keys[:3], gap=5.0, center=x+w/2)
F.key_row(None, y+h-17.4, keys[3:5], gap=5.0, center=x+w/2)
F.key_row(None, y+h-7.4, keys[5:], gap=5.0, center=x+w/2)
F.values['c'] = Cv

# ================================================================ d  aligned effects
x, y, w, h = F.panel('d', X0, Y2, RR-X0, H2)
pw = 92.0; gap_ = (w-30-3.5-3*pw)/2
Dv = {}
for j, (sp, stt, gl) in enumerate(GROUPS):
    q = prim[(prim.species == sp) & (prim.state == stt) & prim.pair_measurable]
    sx = Lin(-0.5, 0.8, x+30+(pw+gap_)*j, x+30+(pw+gap_)*j+pw)
    sy = Lin(-0.5, 0.8, y+16+pw, y+16)
    r1, r2 = q.route_1.iloc[0], q.route_2.iloc[0]
    F.text((sx.p0+sx.p1)/2, y+3.0, f'{gl} · {len(q)} genes', 7.5, ha='center')
    F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref'); F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    F.line([sx.p0, sx.p0], [sy.p0, sy.p1], lw='axis', color=AXC)
    F.yaxis(sy, sx.p0, [-0.4, 0, 0.4, 0.8], ['−.4', '0', '.4', '.8'], line=False)
    for k, cl, cc in CLS:
        z = q[q.pair_class == k]
        if len(z): F.points([sx.clip(v) for v in z.aligned_delta_1], [sy.clip(v) for v in z.aligned_delta_2], ms='small',
                            color=cc if cc != 'white' else AXC, alpha='secondary')
    F.xaxis(sx, sy.p0+1.0, [-0.4, 0, 0.4, 0.8], ['−.4', '0', '.4', '.8'])
    n1, n2 = STUDY[r1][0].replace('Schiebinger2019–', '').replace('Guo2019–', ''), STUDY[r2][0].replace('Schiebinger2019–', '')
    F.text((sx.p0+sx.p1)/2, sy.p0+12.6, f'{n1}, aligned Δ', 7.5, ha='center')
    F.vtext(sx.p0-19.0, (sy.p0+sy.p1)/2, f'{n2}, aligned Δ', 7.5)
    Dv[gl] = dict(pairs=int(len(q)), route_1=r1, route_2=r2)
F.key_row(None, y+h-7.4, [('marker', ('o', cc), cl) for k, cl, cc in CLS[:5]], gap=6.0, center=x+w/2)
F.values['d'] = Dv
F.notes.append('d: aligned Δ = intermediate-minus-source change multiplied by the endpoint direction; values beyond the axes drawn at the limit')

# ================================================================ e  settings
x, y, w, h = F.panel('e', X0, Y3, RR-X0, H3)
fo = D.df('D13E', P+'fixed_cohort_outcomes.tsv')
VARS = [('shared_p08', 'Purity ≥ .8 (primary)'), ('shared_p09', 'Purity ≥ .9'), ('shared_depth40', '40 PAS UMIs'),
        ('native_PAS_coordinates', 'Native PAS axis'), ('exclude_L00', 'Omit L00'), ('exclude_C012', 'Omit C012'), ('exclude_both', 'Omit both')]
labw = max(textwidth(l, 7.0) for _, l in VARS)
e1 = Lin(0, 0.6, x+labw+8, x+labw+8+150)
e2 = Lin(0, 1, e1.p1+24, x+w-80)
F.text((e1.p0+e1.p1)/2, y+3.0, 'Endpoint-aligned / measurable (requalified)', 7.5, ha='center')
F.text((e2.p0+e2.p1)/2, y+3.0, 'Stable class, fixed primary cohort', 7.5, ha='center')
MK = {'Human R': ('o', STUDY['GSE118258'][1]), 'Mouse R': ('s', STUDY['CIP'][1]), 'Mouse E': ('D', STUDY['WOT'][1])}
Ev = {}; rh = 9.0; r0 = y+21.0
for i, (v, vl) in enumerate(VARS):
    yy = r0+rh*i
    if i % 2 == 0: F.rect(x, yy-rh/2, w, rh, fc=BAND)
    F.text(x+0.4, yy-3.4, vl, 7.0, font='Bol' if v == 'shared_p08' else 'Reg')
    for sp, stt, gl in GROUPS:
        q = ss[(ss.variant == v) & (ss.scope == 'requalified') & (abs(ss.intermediate_threshold-0.02) < 1e-9) & (ss.species == sp) & (ss.state == stt)]
        if len(q) == 0: continue
        meas = int(q.pair_measurable.sum()); al = int((q.pair_class == 'endpoint_aligned').sum())
        m_, c_ = MK[gl]
        if meas: F.marker(e1(al/meas), yy, m_, 'small', fc=c_)
        f_ = fo[(fo.variant == v) & (fo.scope == 'fixed_primary_admission') & (fo.species == sp) & (fo.state == stt)]
        if len(f_):
            f_ = f_.iloc[0]
            F.marker(e2(f_.stable/f_.fixed_cohort), yy, m_, 'small', fc=c_)
            Ev[f'{v}|{gl}'] = dict(aligned=al, measurable=meas, stable=int(f_.stable), fixed_cohort=int(f_.fixed_cohort),
                                   unavailable=int(f_.unavailable))
F.xaxis(e1, r0+rh*len(VARS)-rh/2+1.0, [0, 0.2, 0.4, 0.6], ['0', '.2', '.4', '.6'])
F.xaxis(e2, r0+rh*len(VARS)-rh/2+1.0, [0, 0.5, 1], ['0', '.5', '1'])
F.text(x+w, r0-3.4, 'q90 / q95 support:', 7.0, ha='right', color=MUTED)
F.text(x+w, r0+5.6, 'no measurable genes', 7.0, ha='right', color=MUTED)
F.key_row(None, y+h-7.4, [('marker', MK[g][::-1][::-1], g) for _, _, g in GROUPS], gap=8.0, center=(e1.p0+e2.p1)/2)
F.values['e'] = Ev
F.save(OUT, D)
