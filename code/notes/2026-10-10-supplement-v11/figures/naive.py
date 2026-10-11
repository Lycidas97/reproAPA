"""Supplementary Figure S12 (v11): specificity of the shared proximal response to naive pluripotency.
Each row is one prespecified direct comparison: the human naive-related state (hN, N6/N8 pooled against P7) against one
human control contrast, both overlapped with the same mouse proximal set on the common one-to-one background of that
comparison (model: fine joint assay + mouse expression + human expression; 20,000 joint permutations).
a  overlap with the mouse proximal set, observed/expected, naive and control (|dDWUI| >= 0.02)
b  naive minus control fold, 0.02, against the paired null (95% range), q over the 14 comparisons
c  the same at 0.05
Data: notes/2026-10-10-proximal-module-specificity/direct_specificity_results.tsv (read-only)."""
import numpy as np
from common import *                                     # noqa
from common import supdraw, Reads, name, OUT

R = Reads()
MODEL = 'fine_joint_assay_mouse_expression_human_expression'
d = R.df('notes/2026-10-10-proximal-module-specificity/direct_specificity_results.tsv')
d = d[d.model == MODEL]
TARGETS = [('P_M_recurrent', 'Mouse P–M recurrent proximal set'), ('P_M_PS_same_contributor_dual', 'Mouse dual-reference proximal set')]
CTRL = [('primed_P7_S5_pool', 'Primed P–S, other library'), ('primed_P7_S7', 'Primed P–S'), ('M7_S7', 'M–S'), ('E7_S7', 'E–S'),
        ('primed_P7_M7', 'Primed P–M'), ('human_P_S_recurrent', 'P–S, recurrent'), ('human_P_M_recurrent', 'P–M, recurrent')]
SRC = {'primed_P7_S5_pool': 'Liu2020', 'primed_P7_S7': 'Liu2020', 'M7_S7': 'Liu2020', 'E7_S7': 'Liu2020', 'primed_P7_M7': 'Liu2020',
       'human_P_S_recurrent': 'all', 'human_P_M_recurrent': 'all'}
NAIVE = FAMILY_LOCAL['hN']
X0, XR = 4.5, PAGE_W-4.5
PITCH, GH = 9.0, 13.0
TOP = 4.5+24.0
rows = []
y = TOP
for tk, tl in TARGETS:
    rows.append(('head', tl, y)); y += GH
    for ck, cl in CTRL:
        rows.append((tk, ck, y)); y += PITCH
    y += 2.0
YB = y-2.0-PITCH/2+1.0
H = YB+24.6
F = supdraw.Fig(H, name('naive'))
WA, WB = 214.0, 120.0
xa, xb = X0, X0+WA+GAP
xc = xb+WB+GAP
for L, x0, w in (('a', xa, WA), ('b', xb, WB), ('c', xc, XR-xc)):
    F.panel(L, x0, 4.5, w, H-9.0)

def get(tk, ck, thr):
    r = d[(d.mouse_target == tk) & (d.control == ck) & np.isclose(d.threshold, thr)]
    assert len(r) == 1, (tk, ck, thr)
    return r.iloc[0]

# bands
for i, (tk, ck, yy) in enumerate(rows):
    if tk != 'head' and CTRL.index((ck, dict(CTRL)[ck])) % 2 == 0:
        for x0, w in ((xa, WA), (xb, WB), (xc, XR-xc)):
            F.use('a' if x0 == xa else 'b' if x0 == xb else 'c')
            F.rect(x0, yy-PITCH/2, w, PITCH, fc=BAND)
# ---------------------------------------------------------------- a
F.use('a')
sa = Lin(0.9, 1.9, xa+112.0, xa+WA-30.0)
F.text((sa.p0+sa.p1)/2, 4.5+0.6, 'Overlap with mouse set', 7.5, ha='center')
F.text((sa.p0+sa.p1)/2, 4.5+10.6, 'Obs./expected, |ΔDWUI| ≥ .02', 7.0, ha='center', color=MUTED)
F.text(xa+WA, 4.5+10.6, 'Pairs', 7.0, ha='right', color=MUTED)
F.line([sa(1)]*2, [TOP-4.5, YB], lw='ref', color=AXC, ls='ref')
for tk, ck, yy in rows:
    if tk == 'head':
        F.text(xa+0.2, yy-3.0, ck, 7.5); continue
    r = get(tk, ck, 0.02)
    nf, cf = r.naive_observed/r.naive_expected, r.control_observed/r.control_expected
    F.text(xa+6.0, yy-3.3, dict(CTRL)[ck], 7.0)
    F.line([sa(nf), sa(cf)], [yy, yy], lw='thin', color=MUTED)
    F.marker(sa(cf), yy, 'o', 'small', fc='white', ec=INK)
    F.marker(sa(nf), yy, 'o', 'small', fc=NAIVE, ec=NAIVE)
    F.text(xa+WA, yy-3.3, f'{int(r.n_orthologs)}', 7.0, ha='right', color=MUTED)
    F.values[f'a {tk}|{ck}'] = dict(n=int(r.n_orthologs), naive=[int(r.naive_observed), float(r.naive_expected)],
                                    control=[int(r.control_observed), float(r.control_expected)], naive_fold=float(nf), control_fold=float(cf))
ya = F.xaxis(sa, YB, [1, 1.5], ['1', '1.5'])
# ---------------------------------------------------------------- b, c
def forest(L, x0, w, thr, lim, ticks, labs):
    F.use(L)
    s = Lin(lim[0], lim[1], x0+5.0, x0+w-20.0)
    F.text(x0+16.0, 4.5+0.6, f'Naive − control', 7.5)
    F.text(x0+16.0, 4.5+10.6, f'Fold difference, {thr:.2f}'.replace('0.', '.'), 7.0, color=MUTED)
    F.text(x0+w, 4.5+10.6, 'q', 7.0, ha='right', color=MUTED)
    F.line([s(0)]*2, [TOP-4.5, YB], lw='ref', color=AXC, ls='ref')
    for tk, ck, yy in rows:
        if tk == 'head': continue
        r = get(tk, ck, thr)
        q = float(r.q_specificity14_within_model_threshold)
        F.line([s(r.null_p025), s(r.null_p975)], [yy, yy], lw='range', color=NULLBAR)
        F.marker(s(r.fold_difference), yy, 'o' if thr == 0.02 else 'D', 'small', fc=INK if q < 0.05 else 'white', ec=INK)
        F.text(x0+w, yy-3.3, f'{q:.2f}'.lstrip('0'), 7.0, ha='right')
        F.values[f'{L} {tk}|{ck}'] = dict(diff=float(r.fold_difference), null=[float(r.null_p025), float(r.null_p975)], p=float(r.p), q14=q)
    return F.xaxis(s, YB, ticks, labs)
yb = forest('b', xb, WB, 0.02, (-1.5, 1.5), [-1, 0, 1], ['−1', '0', '1'])
yc = forest('c', xc, XR-xc, 0.05, (-4.0, 2.0), [-4, -2, 0, 2], ['−4', '−2', '0', '2'])
yk = max(ya, yb, yc)+3.0
F.use('a')
F.key_row(0, yk, [('marker', ('o', NAIVE), 'Naive (hN–P)'), ('marker', ('o', 'white', INK), 'Control')], center=(sa.p0+sa.p1)/2)
F.use('b')
F.key_row(0, yk, [('line', (NULLBAR, None, 'range'), 'Null 95%'), ('marker', ('o', 'white', INK), 'q ≥ .05')], center=xb+(WB-20.0)/2)
F.use('c')
F.key_row(0, yk, [('marker', ('D', 'white', INK), 'q ≥ .05')], center=xc+(XR-xc-20.0)/2)
print('needed', yk+7.6+4.5, 'H', H)
F.save(OUT, R)
