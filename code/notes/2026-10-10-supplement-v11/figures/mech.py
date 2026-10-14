"""Supplementary Figure S19 (v11): explanations tested for the shared P–M proximal set (27 one-to-one orthologs, |dDWUI| >= 0.02).
a  S–M on the shared core (no P cells): cycle-score difference and core mean dDWUI per route (S above M = positive)
b  Schiebinger2019–2i within-date dose: correlation of the cycle-score increment with proximal strength across wells
c  candidate mechanisms and functions, observed/expected hits: ESCC seeds in full distal regions, external QKI-knockdown
   proximal and EMT distal targets, membrane GO categories (27 core), RNA-regulation GO categories (48 shared orthologs)
d  external U1 perturbation (Oh et al. 2020): fraction of measurable genes with the reciprocal response, core vs background\ne  PAS anchor architecture of the core against assay-matched expectation (q over 12 tests)
f–i former Supplementary Figure S15 a–d (cis-feature prediction of cross-species membership), vector crops, relettered;
   former S15 e–h are kept in Table S5e.
Data: notes/2026-10-10-cycle-source-attribution, shared-proximal-mechanisms, conserved-response (read-only);
former S15: notes/2026-10-06-supplement-polyapipe/consolidated/O/SuppFig.S15.pdf (= production Fig/SuppFig.S15.pdf)."""
import numpy as np
from common import *                                     # noqa
from common import supdraw, Reads, name, OUT, ROOT, crop, place, new_page, finish, sha

R = Reads()
X0, XR = 4.5, PAGE_W-4.5
import os, pdfplumber
# 2026-10-12: f–i are cropped from a composed page holding only former S15 a–d (consolidated O, bands D19 ab / D20 ab),
# rebuilt after the in-figure text cleanup (V11_S19_FI); letters are read from that page, the whole page is kept.
S15 = Path(os.environ['V11_S19_FI']) if os.environ.get('V11_S19_FI') else ROOT/'notes/2026-10-06-supplement-polyapipe/consolidated/O/SuppFig.S15.pdf'
with pdfplumber.open(S15) as _pdf:
    _pg = _pdf.pages[0]; H15 = float(_pg.height)
    _L = {w_['text']: (w_['x0'], w_['top']) for w_ in _pg.extract_words(extra_attrs=['size']) if w_['size'] > 10.5 and len(w_['text']) == 1}
OLD = {k_: (round(_L[k_][0], 1), round(_L[k_][1]-0.4, 1)) for k_ in 'abcd'}
H1, H2 = 97.0, 171.0   # 2026-10-12: room for the axis titles that replace the panel titles
Y_OLD = 4.5+H1+GAP+H2+GAP                  # top of the kept former-S15 rows
OLD_REG = (0.0, 0.0, PAGE_W, H15-4.5)        # whole composed page above its bottom margin
H = Y_OLD+(OLD_REG[3]-4.5)+4.5
F = supdraw.Fig(H, name('mech'))
y1, y2 = 4.5, 4.5+H1+GAP

# ---------------------------------------------------------------- a  S–M on the core
ROUTES = ['GSE118258', '178325_0618', 'GSE221739', 'GSE150311', 'CIP', 'WOT']
sm = R.df('notes/2026-10-10-cycle-source-attribution/sm-core/SM_route_summary.tsv')
sm = sm[(sm.design == 'main') & (sm.background == 'available_contributors') & (sm.gene_set == 'core27')]
WA = 262.0
x, y, w, h = F.panel('a', X0, y1, WA, H1)
# 2026-10-12: panel titles of a–e moved to axis titles or the caption
sc1 = Lin(-0.1, 0.4, x+86.0, x+150.0); sc2 = Lin(-0.06, 0.06, x+168.0, x+WA-34.0)
F.text((sc1.p0+sc1.p1)/2, y+3.0, 'Cycle score', 7.0, ha='center', color=MUTED)
F.text((sc2.p0+sc2.p1)/2, y+3.0, 'Mean ΔDWUI', 7.0, ha='center', color=MUTED)
F.text(XR if False else x+WA, y+3.0, 'Prox. | dist.', 7.0, ha='right', color=MUTED)
top, pitch = y+21.0, 9.8
for k, r in enumerate(ROUTES):
    g = sm[sm.route == r].iloc[0]
    yy = top+k*pitch
    if k % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
    F.rect(x+0.2, yy-2.8, 5.6, 5.6, fc=STUDY[r][1])
    F.text(x+7.8, yy-3.3, STUDY[r][0]+(' (obs. S)' if g.source_type == 'observed_S' else ''), 7.0)
    F.marker(sc1(g.mean_delta_Cycle_exact), yy, 'o', 'small', fc=INK)
    col = PROX if g.mean_DWUI < 0 else DIST
    F.marker(sc2(g.mean_DWUI), yy, 'o', 'small', fc=col, ec=col)
    F.text(x+WA, yy-3.3, f'{int(g.proximal)} | {int(g.distal)}', 7.0, ha='right')
    F.values[f'a {r}'] = dict(source=g.source_type, N=int(g.N), cycle_S_minus_M=float(g.mean_delta_Cycle_exact),
                              core_mean_dDWUI_S_minus_M=float(g.mean_DWUI), proximal=int(g.proximal), distal=int(g.distal))
yb = top+len(ROUTES)*pitch-pitch/2+1.0
for s_ in (sc1, sc2): F.line([s_(0)]*2, [top-pitch/2, yb], lw='ref', color=AXC, ls='ref')
F.xaxis(sc1, yb, [0, 0.2, 0.4], ['0', '.2', '.4'])
ya = F.xaxis(sc2, yb, [-0.05, 0, 0.05], ['−.05', '0', '.05'])
F.key_row(0, ya+2.0, [('marker', ('o', PROX), 'S more proximal'), ('marker', ('o', DIST), 'S more distal')], center=(sc2.p0+sc2.p1)/2)

# ---------------------------------------------------------------- b  WOT within-date dose
xb = X0+WA+GAP
x, y, w, h = F.panel('b', xb, y1, XR-xb, H1)
fp = R.df('notes/2026-10-10-cycle-source-attribution/wot-dose/fixed_panel_descriptive_statistics.tsv')
cf = R.df('notes/2026-10-10-cycle-source-attribution/wot-dose/core_formal_statistics_replay.tsv')
rows = []
for pan, lab in (('native85', '85 fixed genes'), ('ortholog45', '45 orthologs'), ('mouse_proximal8', '8 proximal recurrent')):
    g = fp[(fp.panel == pan) & (fp.predictor == 'delta_Cycle_exact_MC') & (fp.outcome == 'Y_mean') & (fp.batch_constant_centered == False)]
    if len(g): rows.append((lab, float(g.r.iloc[0]), float(g.leave_one_date_r_min.iloc[0]), float(g.leave_one_date_r_max.iloc[0]), None, int(g.n_dates.iloc[0])))
c0 = cf[cf.predictor == 'Cycle'].iloc[0]
rows.append(('Core subset', float(c0.r), float(c0.leave_one_date_r_min), float(c0.leave_one_date_r_max), float(c0.q2_BH), int(c0.n_dates)))
sb = Lin(-1, 1, x+66.0, XR-22.0)
F.text(XR-0.2, y+3.0, 'q', 7.0, ha='right', color=MUTED)
top, pitch = y+23.0, 12.8
for k, (lab, rr, lo, hi, q, nd) in enumerate(rows):
    yy = top+k*pitch
    if k % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
    F.text(x+0.2, yy-3.3, lab, 7.0)
    F.line([sb(lo), sb(hi)], [yy, yy], lw='bold', color=NULLBAR)
    F.marker(sb(rr), yy, 'o', 'small', fc=INK if (q is not None and q < 0.05) else 'white', ec=INK)
    F.text(XR-0.2, yy-3.3, '–' if q is None else f'{q:.2f}'.lstrip('0'), 7.0, ha='right', color=INK if q is not None else MUTED)
    F.values[f'b {lab}'] = dict(r=rr, leave_one_date=[lo, hi], q2=q, n_dates=nd)
yb = top+len(rows)*pitch-pitch/2+1.0
F.line([sb(0)]*2, [top-pitch/2, yb], lw='ref', color=AXC, ls='ref')
yb_ = F.xaxis(sb, yb, [-1, 0, 1], ['−1', '0', '1'], title='r, cycle vs proximal strength')
F.key_row(0, yb_+2.0, [('line', (NULLBAR, None, 'bold'), 'Leave one date out')], center=(sb.p0+sb.p1)/2)

# ---------------------------------------------------------------- c  candidate mechanisms
ms = R.df('notes/2026-10-10-shared-proximal-mechanisms/mechanism_summary.tsv')
ca = R.df('notes/2026-10-10-conserved-response/function/category_results.tsv')
ca = ca[(ca.primary == True) & (ca.set_name == 'main48')]
def mrow(branch, test, model):
    g = ms[(ms.branch == branch) & (ms.test == test) & (ms.model == model)]
    assert len(g) == 1, (branch, test, model); return g.iloc[0]
GROUPS = [
    ('Long-isoform clearance, 67 pairs', [('ESCC seeds, both species', mrow('decay', 'both_ESCC_ge7', 'coarse'))]),
    ('External targets, 605 genes', [('QKI knockdown, proximal', mrow('met-qki', 'QKI_KD_proximal', 'COARSE')),
                                     ('EMT, distal', mrow('met-qki', 'EMT_distal', 'COARSE'))]),
    ('Membrane functions, 751 pairs', [('ER membrane', mrow('membrane', 'ER_membrane', 'joint_measurement_PAS_capacity')),
                                       ('ER to Golgi transport', mrow('membrane', 'ER_to_Golgi_transport', 'joint_measurement_PAS_capacity')),
                                       ('Vesicle transport', mrow('membrane', 'vesicle_transport', 'joint_measurement_PAS_capacity')),
                                       ('Cytoskeleton organization', mrow('membrane', 'cytoskeleton_organization', 'joint_measurement_PAS_capacity'))]),
]
rna = []
for cat in ('RNA binding', 'mRNA processing', 'regulation of translation'):
    g = ca[ca.category_name == cat].iloc[0]
    rna.append(({'RNA binding': 'RNA binding', 'mRNA processing': 'mRNA processing', 'regulation of translation': 'Regulation of translation'}[cat], dict(observed=float(g.hits), expected=float(g.expected), q=float(g.q_primary6), N=int(g.N), n=int(g.cases))))
WC = 286.0
x, y, w, h = F.panel('c', X0, y2, WC, H2)
scc = Lin(0, 6, x+112.0, x+WC-46.0)
F.text(x+WC-18.0, y+3.0, 'n/exp.', 7.0, ha='right', color=MUTED); F.text(x+WC, y+3.0, 'q', 7.0, ha='right', color=MUTED)
yy = y+20.0; pitch = 8.85; k = 0; ylist = []
allg = GROUPS+[('RNA regulation, 48 shared orthologs', [(lab, v) for lab, v in rna])]
for gl, items in allg:
    F.text(x+0.2, yy-3.0, gl.split(',')[0], 7.0, color=MUTED); yy += pitch   # 2026-10-12: background sizes in the caption; values keys keep the full label
    for lab, r in items:
        o, e, q = (float(r['observed']), float(r['expected']), float(r['q'])) if isinstance(r, dict) else (float(r.observed), float(r.expected), float(r.q))
        if k % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
        k += 1
        F.text(x+6.0, yy-3.3, lab, 7.0)
        F.marker(scc(o/e), yy, 'o', 'small', fc=INK if q < 0.05 else 'white', ec=INK)
        F.text(x+WC-18.0, yy-3.3, f'{o:g}/{e:.1f}', 7.0, ha='right')
        F.text(x+WC, yy-3.3, f'{q:.2f}'.lstrip('0'), 7.0, ha='right')
        F.values[f'c {gl}|{lab}'] = dict(observed=o, expected=e, ratio=o/e, q=q)
        ylist.append(yy); yy += pitch
yb = yy-pitch/2+1.0
F.line([scc(1)]*2, [ylist[0]-pitch/2, yb], lw='ref', color=AXC, ls='ref')
yc = F.xaxis(scc, yb, [0, 1, 2, 4, 6], ['0', '1', '2', '4', '6'], title='Observed / expected')
F.key_row(0, yc+2.0, [('marker', ('o', 'white', INK), 'q ≥ .05')], center=(scc.p0+scc.p1)/2)

# ---------------------------------------------------------------- d  U1 reciprocal response, e  anchor architecture
xd = X0+WC+GAP
HD = 81.0
x, y, w, h = F.panel('d', xd, y2, XR-xd, HD)
u = ms[ms.branch == 'u1-kinetics']
sd = Lin(0, 1, x+24.0, XR-46.0)
top, pitch = y+19.0, 14.0
for k, (t, lab) in enumerate((('U1_reciprocal_log2_t0.0', '> 0'), ('U1_reciprocal_log2_t0.2', '≥ .2'), ('U1_reciprocal_log2_t0.5', '≥ .5'))):
    g = u[u.test == t].iloc[0]
    yy = top+k*pitch
    F.text(x+0.2, yy-1.4, lab, 7.0)
    for j, (num_, den, col, nm) in enumerate(((g.observed, g.n, INK, 'core'), (g.other_hits, g.other_n, NULLBAR, 'other'))):
        yb_ = yy-2.8+j*5.6
        F.rect(sd(0), yb_-1.8, sd(1)-sd(0), 3.6, fc=TRACK)
        F.rect(sd(0), yb_-1.8, sd(num_/den)-sd(0), 3.6, fc=col)
        if j == 0: F.text(XR-0.2, yy-3.3, f'{int(g.observed)}/{int(g.n)} | {int(g.other_hits)}/{int(g.other_n)}', 6.8, ha='right')
        F.values[f'd {t} {nm}'] = dict(hits=int(num_), n=int(den), fraction=float(num_/den))
yb = top+2*pitch+5.0
yd = F.xaxis(sd, yb, [0, 0.5, 1], ['0', '.5', '1'], title='Reciprocal response (fraction)')
F.key_row(0, yd+2.0, [('block', INK, 'Core'), ('block', NULLBAR, 'Other')], center=(sd.p0+sd.p1)/2)
print('d end', round(yd+2.0+7.6-y, 1))

ye0 = y2+HD+GAP
x, y, w, h = F.panel('e', xd, ye0, XR-xd, y2+H2-ye0)
arc = ms[(ms.branch == 'u1-kinetics') & (ms.model == 'assay_architecture')]
se = Lin(-0.6, 0.6, x+70.0, XR-18.0)
F.text(XR-0.2, y+3.0, 'q', 7.0, ha='right', color=MUTED)
top, pitch = y+20.0, 7.0
LAB = {'proximal_core_rank': 'proximal rank', 'core_rank_difference': 'rank difference', 'log_anchor_span': 'log span'}
k = 0
for sp, spl in (('human', 'Hs'), ('mouse', 'Mm')):
    for key, lab in LAB.items():
        g = arc[arc.test == f'{sp}_{key}'].iloc[0]
        yy = top+k*pitch
        if k % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
        F.text(x+0.2, yy-3.3, f'{spl} {lab}', 7.0)
        dv = float(g.observed-g.expected)
        F.marker(se.clip(dv), yy, 'o', 'small', fc=INK if g.q < 0.05 else 'white', ec=INK)
        F.text(XR-0.2, yy-3.3, f'{g.q:.2f}'.lstrip('0'), 7.0, ha='right')
        F.values[f'e {sp}_{key}'] = dict(observed=float(g.observed), expected=float(g.expected), q12=float(g.q), n_core=int(g.n), N=int(g.N))
        k += 1
yb = top+k*pitch-pitch/2+1.0
F.line([se(0)]*2, [top-pitch/2, yb], lw='ref', color=AXC, ls='ref')
yee = F.xaxis(se, yb, [-0.5, 0, 0.5], ['−.5', '0', '.5'], title='Core − expected')
print('e end', round(yee-y, 1), 'h', round(h, 1))

# ---------------------------------------------------------------- e–h  former S15 a–d, vector crops
dst = new_page(H)
LB = []
for L, (lx, ly) in OLD.items():
    LB.append((lx-0.5, ly-0.5, lx+8.5, ly+11.5))           # clear the old letter
blk = crop(S15, H15, OLD_REG, extra_clear=LB)
dy = Y_OLD-4.5
place(dst, H, blk, OLD_REG, 0.0, dy)
for (L, (lx, ly)), new in zip(OLD.items(), 'fghi'):
    F.t.text_top(lx+0.2, ly+0.4+dy, new, 11.0, font='Bol')
extra = dict(kept=dict(source=str(S15.relative_to(ROOT)), source_sha=sha(S15), region=OLD_REG, dy=dy,
                       letters={'a': 'f', 'b': 'g', 'c': 'h', 'd': 'i'}, dropped=['former S15 e–h (in Table S5e)']))
finish(dst, F, H, R, extra)
print('H', H)
