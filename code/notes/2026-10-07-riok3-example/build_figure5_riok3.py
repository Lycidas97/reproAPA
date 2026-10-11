"""Figure 5 · records adapter of the v5.1 production builder (layout unchanged; data from one result bundle).

Usage: build.py --records <bundle>/result-records.json --out <new dir>
Panel notes below describe the production figure drawn from the archived v2 data; their numbers are history, not checks.
Example genes, exact sites and the shared interval come from the bundle through display-selection.json.

Base: locked master Fig/Figure.5.pdf. Contract rows 5A–5F of main-panel-replacement-v5.tsv.
  A  observed interval direction agreement + descriptive gene-bootstrap CI + pairs/genes. No C0, no Cobs−C0.
  B  genes with net weighted predicted-site loss / gain (20 UMI), one shared zero-based scale for miRNA and RBP.
  C  equal-width screening flow 64 → 5 → 1 genes (mouse P–M), 13 interval pairs annotated beside the 5.
  D  partial retention (%) per route × contrast × category; P–S/P–M union relationship table; legend relabelled.
  E  Eif4ebp2 GGCAGTG two-site relation (changed / retained) and inclusion of the interval holding the changed site.
  F  Rad23b: family PAS usage M vs P, native MPRO / DWUI with n; interval inclusion M→P per paired library;
     site annotations counted, no single-seed highlight.
"""
import csv, gzip, json, sys, math
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'2026-10-06-main-figures-polyapipe/builders/tools-local'))
from v5draw import *  # noqa
from figrec import Run, MASTERS, pas_pos, pas_strand
HERE = Path(__file__).resolve().parent
R = Run('F5')
import numpy as np

MASTER = MASTERS/'Figure.5.pdf'      # locked master copy (Fig/ now holds promoted production files)
assert sha(MASTER) == 'd2b29ffa84a38374af32890e1ec3f43c0cdd6c331d3c66db76adfb8a634b1b3c'
OUT = R.out/'Figure.5.draft.pdf'
def tsv(p):
    op = gzip.open if str(p).endswith('.gz') else open
    with op(p, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))
use = R.rows

pl = pdfplumber.open(MASTER).pages[0]; W, H = pl.width, pl.height
pdf = pikepdf.open(MASTER)
CLEAR = {
    'A': [(55, 0, 334, 98), (6.5, 20.5, 53.0, 68.5)],
    'B': [(52, 113.5, 334, 190), (80, 103.5, 310, 113.2)],
    'C': [(343.5, 0, 454, 192), (334.5, 18, 343.5, 192), (333.0, 15.0, 449.0, 25.0), (333.0, 90.5, 449.0, 100.5)],
    'D': [(17, 213, 60, 247), (120, 211, 230, 253), (38, 262, 230, 322.5)],
    'E': [(249, 190, 454, 378), (239.5, 209.5, 249, 378)],
    'F': [(14, 379, 454, 527.2), (0, 391.5, 14, 527.2)],
}
removed = {k: len(clear(pdf, H, v)) for k, v in CLEAR.items()}
t = Top(W, H); ver = {'master_sha256': sha(MASTER), 'removed_objects': removed}
STUDY = {'CIP': ('#CF4C00', 's'), 'WOT': ('#207450', 'o'), 'GSE118258': ('#753399', 'o'), 'GSE150311': ('#648290', 'D')}
SC = dict(S='#E69F00', R='#0072B2', E='#009E73', P='#CC79A7', M='#8C6D31')

# ------------------------------------------------------------------ A
rc = use('F5A', 'panel_statistics/observed-interval-agreement.tsv')
XA = lambda v: 66.0+(v-0.5)*168.0
XR = lambda v: 167.0+(v+0.1)/1.1*50.0  # polyApipe redraw 2026-10-06: axis rescaled within the same span
t.marker(68.9, 8.8, 'o', 3.6, mew=0.3); t.text_top(78.4, 5.1, 'Observed', 8)
t.line([123.0, 133.0], [8.8, 8.8], lw=0.7, color='#263947'); t.text_top(136.8, 5.1, 'Gene-bootstrap 95% CI', 8)
_x = 136.8+textwidth('Gene-bootstrap 95% CI', 8)+12.0
t.marker(_x, 8.8, '^', 3.0, fc='#8A959C', ec='#8A959C', mew=0.3); t.text_top(_x+4.5, 5.1, 'miRNA', 8); _x += 4.5+textwidth('miRNA', 8)+7.0
t.marker(_x, 8.8, 'v', 3.0, fc='#8A959C', ec='#8A959C', mew=0.3); t.text_top(_x+4.5, 5.1, 'RBP', 8); _x += 4.5+textwidth('RBP', 8)
assert _x < 330.0, _x
t.text_top(222.0, 91.8, 'Agreeing / all', 8)
t.text_top(320.0, 91.8, 'Genes', 8, ha='right')
Ad = []
for (sp, con), y in zip([('mouse', 'SP'), ('mouse', 'MP'), ('human', 'SP'), ('human', 'MP')], [26.0, 41.0, 56.0, 71.0]):
    r = [r for r in rc if r['species'] == sp and r['contrast'] == con and r['category'] == 'all' and r['unit'] == 'region_pair']
    assert len(r) == 1; r = r[0]
    c, lo, hi, n, g = float(r['cobs']), float(r['cobs_ci_low']), float(r['cobs_ci_high']), int(r['n']), int(r['genes'])
    if y in (26.0, 56.0):
        t.rect(4.0, y-7.0, 317.5, 14.0, fc=BAND)
    t.text_top(7.1, y-3.9, {'SP': 'P–S', 'MP': 'P–M'}[con] and f"{sp.title()} {({'SP': 'P–S', 'MP': 'P–M'})[con]}", 8)
    if hi > lo:
        t.line([XA(lo), XA(hi)], [y, y], lw=0.7, color='#263947')
        for v in (lo, hi): t.line([XA(v)]*2, [y-1.35, y+1.35], lw=0.7, color='#263947')
    small = g < 20
    t.marker(XA(c), y, 'o', 3.6, fc='white' if small else INK, ec=INK, mew=0.8 if small else 0.3)
    rho, rlo, rhi = float(r['rho_mean']), float(r['rho_mean_ci_low']), float(r['rho_mean_ci_high'])
    if not (R.axis('F5A', f'{sp} {con} agreement', [c, lo, hi], 0.5, 1.0) & R.axis('F5A', f'{sp} {con} effect rho', [rho, rlo, rhi], -0.1, 1.0)):
        c, lo, hi = (min(max(v_, 0.5), 1.0) for v_ in (c, lo, hi)); rho, rlo, rhi = (min(max(v_, -0.1), 1.0) for v_ in (rho, rlo, rhi))
    t.line([XR(rlo), XR(rhi)], [y, y], lw=0.7, color='#263947')
    for v in (rlo, rhi): t.line([XR(v)]*2, [y-1.35, y+1.35], lw=0.7, color='#263947')
    t.marker(XR(rho), y, 's', 3.4, fc='white' if small else INK, ec=INK, mew=0.8 if small else 0.3)
    subs = {}
    for cat_, dy_, mk_ in [('miRNA', -4.3, '^'), ('RBP', 4.3, 'v')]:
        q_ = [q for q in rc if q['species'] == sp and q['contrast'] == con and q['category'] == cat_ and q['unit'] == 'region_pair']
        assert len(q_) == 1; q_ = q_[0]
        cc, clo, chi = float(q_['cobs']), float(q_['cobs_ci_low']), float(q_['cobs_ci_high'])
        rr_, rlo_, rhi_ = float(q_['rho_mean']), float(q_['rho_mean_ci_low']), float(q_['rho_mean_ci_high'])
        if not (R.axis('F5A', f'{sp} {con} {cat_} agreement', [cc, clo, chi], 0.5, 1.0) & R.axis('F5A', f'{sp} {con} {cat_} rho', [rr_, rlo_, rhi_], -0.1, 1.0)):
            cc, clo, chi = (min(max(v_, 0.5), 1.0) for v_ in (cc, clo, chi)); rr_, rlo_, rhi_ = (min(max(v_, -0.1), 1.0) for v_ in (rr_, rlo_, rhi_))
        for X_, v_, lo_, hi_ in [(XA, cc, clo, chi), (XR, rr_, rlo_, rhi_)]:
            if hi_ > lo_:
                t.line([X_(lo_), X_(hi_)], [y+dy_]*2, lw=0.5, color='#B5C0C7')
            t.marker(X_(v_), y+dy_, mk_, 3.0, fc='#8A959C', ec='#8A959C', mew=0.3, zorder=4)
        subs[cat_] = dict(cobs=cc, ci=[clo, chi], rho_mean=rr_, rho_ci=[rlo_, rhi_], pairs=int(q_['n']), genes=int(q_['genes']))
    c = float(r['cobs']); agree = int(round(c*n))
    if abs(agree-c*n) > 1e-6: R.flag('F5A', 'consistency', f'{sp} {con}: agreement fraction x pairs is not an integer')
    bx0_, bw2 = 222.0, 40.0
    n_ = max(n, 1)
    t.rect(bx0_, y-2.4, bw2*agree/n_, 4.8, fc='#263947'); t.rect(bx0_+bw2*agree/n_, y-2.4, bw2*(n-agree)/n_, 4.8, fc='#C9CFD3')
    t.text_top(bx0_+bw2+3.0, y-4.0, f'{agree}/{n}', 8)
    t.text_top(320.0, y-4.0, str(g), 8, ha='right')
    Ad.append(dict(agree=agree, species=sp, contrast=con, cobs=c, ci_low=lo, ci_high=hi, pairs=n, genes=g, open_marker=small,
                   rho_mean=rho, rho_ci=[rlo, rhi], site_classes=subs))
t.line([XA(0.5), XA(1)], [79.0, 79.0], lw=0.5, color='#637580')
for v, s_ in [(0.5, '50'), (0.75, '75'), (1, '100')]:
    t.line([XA(v)]*2, [79.0, 80.6], lw=0.5, color='#637580'); t.text_top(XA(v), 82.3, s_, 8, ha='center')
t.text_top(XA(0.75), 91.8, 'Direction agreement (%)', 8, ha='center')
t.line([XR(-0.1), XR(1)], [79.0, 79.0], lw=0.5, color='#637580')
for v, s_ in [(0, '0'), (0.5, '.5'), (1, '1')]:
    t.line([XR(v)]*2, [79.0, 80.6], lw=0.5, color='#637580'); t.text_top(XR(v), 82.3, s_, 8, ha='center')
t.text_top(XR(0.5), 91.8, 'Effect ρ', 8, ha='center')
assert XA(0.75)+textwidth('Direction agreement (%)', 8)/2 < XR(0.5)-textwidth('Effect ρ', 8)/2-4
ver['A'] = Ad

# ------------------------------------------------------------------ B
rs = use('F5B', 'weighted_retention/retention_summary.tsv')
gr = use('F5B', 'weighted_retention/gene_retention.tsv')
K = 0.55*80/250; Bd = []; SW = 30.0; SO = 83.0
panels = [('miRNA', 64.0), ('RBP', 200.0)]
# Bar scale: same 44 pt loss span as the locked layout, rescaled to 250 genes (ticks −200/−100/0/100) for the polyApipe counts
mx_loss, mx_gain = 250, 153  # polyApipe redraw 2026-10-06: axis rescaled within the same span
R.axis('F5B', 'genes with net loss', [int(r['negative']) for r in rs], 0, mx_loss)
fl = lambda v: math.copysign(math.log10(1+abs(v)), v)
for cat_, x0_ in panels:
    t.text_top(x0_+(SO+SW)/2-6.0, 106.4, {'miRNA': 'Predicted miRNA target sites', 'RBP': 'RBP motif matches'}[cat_], 7.5, ha='center')   # 2026-10-07 terminology
for cat, x0 in panels:
    z = x0+mx_loss*K; sx0 = x0+SO; SXg = lambda v, sx0=sx0: sx0+SW/2+fl(v)*SW/(2*math.log10(200))
    t.line([z, z], [116.0, 163.9], lw=0.5, color='#A6B2BA')
    t.line([SXg(0)]*2, [116.0, 163.9], lw=0.5, color='#A6B2BA', ls=(0, (1, 1.5)))
    for (sp, con), y in zip([('mouse', 'SP'), ('mouse', 'MP'), ('human', 'SP'), ('human', 'MP')], [121.65, 133.65, 145.6, 157.55]):
        r = [r for r in rs if r['species'] == sp and r['contrast'] == con and r['category'] == cat]; assert len(r) == 1; r = r[0]
        neg, pos, den = int(r['negative']), int(r['positive']), int(r['genes_with_complete_sites'])
        neg_txt, pos_txt = neg, pos
        neg = min(neg, mx_loss)
        if z+pos*K+1.5+textwidth(str(pos), 7.5) >= sx0-1.0:
            R.flag('F5B', 'axis_range', f'{sp} {con} {cat}: {pos} net-gain genes overrun the bar column', data_range=[0, pos])
            while pos > 0 and z+pos*K+1.5+textwidth(str(pos), 7.5) >= sx0-1.0: pos -= 1
        if y in (121.65, 145.6):
            t.rect(x0-12.0, y-5.9, SO+SW+14.0, 11.8, fc=BAND)
        t.rect(z-neg*K, y-3.55, neg*K, 7.1, fc='#263947'); t.rect(z, y-3.55, pos*K, 7.1, fc='#0072B2')
        t.text_top(z-neg*K-1.5, y-3.8, str(neg_txt), 7.5, ha='right'); t.text_top(z+pos*K+1.5, y-3.8, str(pos_txt), 7.5)
        g_ = [float(q['weighted_sites20']) for q in gr if q['species'] == sp and q['contrast'] == con and q['category'] == cat]
        nz = sorted(v for v in g_ if v != 0.0)
        if not (sum(v < 0 for v in nz) == int(r['negative']) and sum(v > 0 for v in nz) == int(r['positive']) and len(nz) == den):
            R.flag('F5B', 'consistency', f'{sp} {con} {cat}: per-gene weighted sites disagree with the summary counts')
        R.axis('F5B', f'{sp} {con} {cat} weighted sites', nz, -199, 199)
        nz = [min(max(v, -199), 199) for v in nz]
        for col_, sel in [('#263947', [v for v in nz if v < 0]), ('#0072B2', [v for v in nz if v > 0])]:
            if sel:
                t.ax.plot([SXg(v) for v in sel], [t.H-(y-3.0+6.0*(k_ % 5)/4) for k_ in range(len(sel))], ls='none', marker='o', ms=1.7,
                          mfc=col_, mec='none', alpha=0.6, zorder=4)
        Bd.append(dict(species=sp, contrast=con, category=cat, net_loss=neg_txt, net_gain=pos_txt, genes_with_complete_sites=den,
                       genes_zero=len(g_)-len(nz), per_gene_weighted_sites20=nz))
    t.line([z-mx_loss*K, z+mx_gain*K], [165.0, 165.0], lw=0.5, color='#63757F')
    for v in (-200, -100, 0, 100):
        t.line([z+v*K]*2, [165.0, 166.6], lw=0.5, color='#63757F'); t.text_top(z+v*K, 168.2, str(abs(v)), 7.5, ha='center')
    t.line([SXg(-199), SXg(199)], [165.0, 165.0], lw=0.5, color='#63757F')
    for v, s_ in [(-50, '−50'), (0, '0'), (50, '+50')]:
        t.line([SXg(v)]*2, [165.0, 166.6], lw=0.5, color='#63757F'); t.text_top(SXg(v), 168.2, s_, 7.5, ha='center')
    t.text_top(sx0+SW, 177.2, 'Weighted sites', 7.5, ha='right')
    assert SXg(50)+textwidth('+50', 7.5)/2 < 319.5, SXg(50)
_x = 7.1; t.text_top(_x, 177.2, 'Genes', 7.5); _x += textwidth('Genes', 7.5)+4.0
t.rect(_x, 178.4, 6.0, 5.0, fc='#263947'); _x += 9.0
t.text_top(_x, 177.2, 'net loss', 7.5); _x += textwidth('net loss', 7.5)+5.0
t.rect(_x, 178.4, 6.0, 5.0, fc='#0072B2'); _x += 9.0; t.text_top(_x, 177.2, 'net gain', 7.5); _x += textwidth('net gain', 7.5)
assert _x < 64.0+SO+SW-textwidth('Weighted sites', 7.5)-4.0, _x
ver['B'] = Bd

# ------------------------------------------------------------------ C
att = {r['stage']: r for r in use('F5C', 'candidate_screen/low_shift_attrition.tsv')}
full = use('F5C', 'candidate_screen/low_shift_candidates.tsv')
ST1, ST2, ST3 = 'low_abs_DWUI_lt002_MPRO_ge01', 'strict_recurrent_interval', 'both_abs_DWUI_lt001'
assert set(att) >= {ST1, ST2, ST3}, sorted(att)
ids_ = lambda s_: [g for g in att[s_]['gene_ids'].split(';') if g]
sym_ = {r['gene_id']: r['gene_symbol'] for r in full}
names = [sym_.get(g, g) for g in ids_(ST2)]; pairs13 = int(att[ST2]['interval_pairs']); one = [sym_.get(g, g) for g in ids_(ST3)]
n64 = len(full)
if n64 != int(att[ST1]['genes']): R.flag('F5C', 'consistency', f'candidate rows {n64} != first-stage genes {att[ST1]["genes"]}')
# the highlighted low-shift example (also 5F): example_selection low_shift_interval restricted to the final stage
exC = R.example('low_shift', 'F5C', use('F5C', 'candidate_screen/example_selection.tsv', lambda r: r['role'] == 'low_shift_interval' and r['gene_id'] in set(ids_(ST3))))
LOWID = exC['gene_id'] if exC else None; LOWSYM = exC['gene_symbol'] if exC else ''
t.text_top(347.0, 6.6, 'Mouse P–M', 8)
bx, bw = 345.0, 106.0
def box(y, h, lines):
    t.rect(bx, y, bw, h, fc='#F3F5F6', ec='#B5C0C7', lw=0.5)
    yy = y+2.6
    for s, size, font, col in lines:
        t.text_top(bx+bw/2, yy, s, size, font=font, ha='center', color=col); yy += size+1.3
def arrow(y0, y1):
    t.line([bx+bw/2]*2, [y0, y1-2.6], lw=0.7, color=MUTED)
    from matplotlib.patches import Polygon
    t.ax.add_patch(Polygon([(bx+bw/2-2.2, H-(y1-3.2)), (bx+bw/2+2.2, H-(y1-3.2)), (bx+bw/2, H-y1)], closed=True, fc=MUTED, ec='none'))
box(17.0, 12.5, [(f'{n64} genes', 8, 'Bol', INK)])
arrow(29.5, 35.5)
box(35.5, 12.5, [(f'{len(names)} strict-interval genes', 8, 'Bol', INK)])
arrow(48.0, 56.0)
# scatter of all 64 candidates: x = max |ΔDWUI| over the two studies, y = min MPRO
sx0, sx1, sy0, sy1 = 361.0, 448.0, 61.0, 160.0
SX = lambda v: sx0+v/0.02*(sx1-sx0); SY = lambda v: sy1-(min(v, 0.5)-0.1)/0.4*(sy1-sy0)
t.rect(SX(0), sy0, SX(0.01)-SX(0), sy1-sy0, fc='#EEF1F3')
pts = [(float(r['max_abs_DWUI']), float(r['min_MPRO']), r['gene_symbol'], int(r['strict_region_pairs'])) for r in full]
if not all(p[0] < 0.02 and p[1] >= 0.1 for p in pts): R.flag('F5C', 'consistency', 'candidates outside the |ΔDWUI| < .02, MPRO ≥ .1 screen')
R.axis('F5C', 'max |ΔDWUI|', [p[0] for p in pts], 0, 0.02)
R.axis('F5C', 'min MPRO (values above .5 sit on the top edge)', [p[1] for p in pts], 0.1, float('inf'))
for xv, yv, g, sp_ in pts:
    if sp_ > 0:
        t.marker(SX(xv), SY(yv), 'o', 3.6, fc=INK, ec=INK, mew=0.4, zorder=4)
    else:
        t.marker(SX(xv), SY(yv), 'o', 2.6, fc='white', ec='#8A959C', mew=0.5, zorder=3)
rx = next((p for p, r_ in zip(pts, full) if r_['gene_id'] == LOWID), None)
lx, ly = SX(0.0003), SY(0.47)
if rx is not None:
    t.marker(SX(rx[0]), SY(rx[1]), 'o', 6.4, fc='none', ec='#CC79A7', mew=1.0, zorder=5)
    t.text_top(lx, ly, LOWSYM, 7.5, font='RegIta', color='#CC79A7')
    t.line([lx+8.0, SX(rx[0])], [ly+8.6, SY(rx[1])-3.4], lw=0.5, color='#CC79A7')
elif exC:
    R.unavailable('F5C', f'{LOWSYM} is not among the plotted candidates')
if any(SX(0.0003)-1 < SX(p[0]) < SX(0.0003)+26 and ly-1 < SY(p[1]) < ly+9 for p in pts): R.flag('F5C', 'overlap', 'a candidate point lies under the example label')
t.line([sx0, sx1], [sy1+1.5]*2, lw=0.5, color=AXC); t.line([sx0-1.5]*2, [sy0, sy1], lw=0.5, color=AXC)
for v, s in [(0, '0'), (0.01, '.01'), (0.02, '.02')]:
    t.line([SX(v)]*2, [sy1+1.5, sy1+3.5], lw=0.5, color=AXC); t.text_top(SX(v), sy1+4.6, s, 7.5, ha='center')
for v, s in [(0.1, '.1'), (0.3, '.3'), (0.5, '.5')]:
    t.line([sx0-3.5, sx0-1.5], [SY(v)]*2, lw=0.5, color=AXC); t.text_top(sx0-4.5, SY(v)-3.8, s, 7.5, ha='right')
t.text_top((sx0+sx1)/2, sy1+12.4, 'Max |ΔDWUI|', 7.5, ha='center')
t.vtext(340.5, (sy0+sy1)/2, 'Min MPRO', 7.5)
t.marker(350.0, 185.2, 'o', 3.6, fc=INK, ec=INK, mew=0.4); t.text_top(354.0, 181.6, 'Strict interval', 7.5)
t.marker(410.0, 185.2, 'o', 2.6, fc='white', ec='#8A959C', mew=0.5); t.text_top(414.0, 181.6, 'Other', 7.5)
ver['C'] = dict(n64=n64, five=names, interval_pairs=pairs13, one=one)

# ------------------------------------------------------------------ D
pr = use('F5D', 'site_statistics/partial_retention_summary.tsv')
cs = use('F5D', 'site_statistics/cross_study_partial_retention.tsv', lambda r: r['exact_partial_recurrence'] == 'True')
R.derive('F5D', 'P–S/P–M union of exact partial relationships counted as distinct (gene, factor) and genes (GAP-D05)')
rel = {}
for r in cs:
    if r['exact_partial_recurrence'] == 'True' and r['contrast'] in ('SP', 'MP'):
        rel.setdefault((r['species'], r['category']), set()).add((r['gene_id'], r['factor']))
tab = {k: (len(v), len({g for g, f in v})) for k, v in rel.items()}
for txt, top in [('Constitutive', 225.4)]:
    t.text_top(17.8, top, txt, 8)
t.text_top(165.0, 198.1, 'miRNA', 8, ha='center') if False else None
for (sp, lab), y in [(('mouse', 'Mouse'), 215.1), (('human', 'Human'), 231.0)]:
    t.text_top(122.9, y, lab, 8)
    t.text_top(167.0, y, f'{tab.get((sp, "miRNA"), (0, 0))[0]} ({tab.get((sp, "miRNA"), (0, 0))[1]})', 8, ha='center')
    t.text_top(209.5, y, f'{tab.get((sp, "RBP"), (0, 0))[0]} ({tab.get((sp, "RBP"), (0, 0))[1]})', 8, ha='center')
t.text_top(174.8, 243.5, 'Relationships (genes)', 7.5, ha='center')
XD = {'miRNA': lambda v: 40.6+v*141.83, 'RBP': lambda v: 137.1+v*141.83}
rows = [('CIP', 'SP', 267.2), ('CIP', 'MP', 274.3), ('WOT', 'SP', 280.6), ('WOT', 'MP', 287.7),
        ('GSE118258', 'SP', 298.5), ('GSE118258', 'MP', 305.5), ('GSE150311', 'SP', 311.8), ('GSE150311', 'MP', 318.9)]
Dd = []
for cat in ('miRNA', 'RBP'):
    X = XD[cat]
    for v in (0, 0.3, 0.6):
        t.line([X(v)]*2, [264.6, 323.3], lw=0.35, color='#D9E0E4')
    for route, con, y in rows:
        r = [r for r in pr if r['route'] == route and r['contrast'] == con and r['category'] == cat]
        if len(r) != 1:
            R.unavailable('F5D', f'{route} {con} {cat}: no partial-retention estimate'); continue
        r = r[0]
        f, lo, hi = float(r['partial_fraction']), float(r['ci_low']), float(r['ci_high'])
        if not R.axis('F5D', f'{route} {con} {cat} partial retention (master axis 0–60 %)', [f, lo, hi], 0.0, 0.6): f, lo, hi = (min(max(v_, 0.0), 0.6) for v_ in (f, lo, hi))
        col, _ = STUDY[route]
        t.whisker(X(lo), X(hi), y, color=col, lw=0.7, cap=1.25)
        if con == 'SP': t.marker(X(f), y, 'o', 3.4, fc=col, ec=col, mew=0.6)
        else: t.marker(X(f), y, 'D', 3.2, fc='white', ec=col, mew=0.8)
        Dd.append(dict(route=route, contrast=con, category=cat, partial_fraction=f, ci_low=lo, ci_high=hi,
                       changed_multisite_pairs=int(r['changed_multisite_pairs'])))
ver['D'] = dict(points=Dd, relationships={f'{a}|{b}': v for (a, b), v in tab.items()})

# ------------------------------------------------------------------ E
specP = R.selection.get('partial_retention_mouse', {'mode': 'rank'})
selP = use('F5E', 'candidate_screen/example_selection.tsv', lambda r: r['role'] == 'partial_retention' and r['species'] == 'mouse')
_want = {specP.get('gene_id')} if specP.get('mode') == 'pin' else set()
_wsym = set() if specP.get('mode') == 'pin' else {r['gene_symbol'] for r in selP}
xs = use('F5E', 'site_statistics/exact_site_cross_study.tsv.gz', lambda r: r['species'] == 'mouse' and (r['gene_id'] in _want or r.get('gene_symbol_a') in _wsym))
T_ = lambda v: v == 'True'
def relation(gid, factor=None, contrast='SP'):
    """Exact two-site relation of one factor: recurrent strict changed sites and sites retained in both routes."""
    rows_ = [r for r in xs if r['gene_id'] == gid and r['contrast'] == contrast and (factor is None or r['factor'] == factor)]
    out = {}
    for r in rows_:
        d_ = out.setdefault((r['category'], r['factor']), dict(changed=set(), retained=set(), names=r.get('names_a', '')))
        if T_(r['recurrent_strict']): d_['changed'].add((int(r['start0']), int(r['end0'])))
        if T_(r['retained_a']) and T_(r['retained_b']): d_['retained'].add((int(r['start0']), int(r['end0'])))
    return {k: v for k, v in out.items() if v['changed'] and v['retained']}
GID5 = None; rel5 = None
if specP.get('mode') == 'pin':
    rels = relation(specP['gene_id'], specP.get('factor'), specP.get('contrast', 'SP'))
    if rels:
        GID5 = specP['gene_id']; GSYM5 = specP['gene_symbol']; (cat5, fac5), rel5 = sorted(rels.items())[0]
    else:
        R.flag('F5E', 'example_ineligible', f'pinned {specP.get("gene_symbol")} {specP.get("factor")} has no exact changed + retained relation in this source')
else:
    for r in sorted(selP, key=lambda r: float(r['rank'])):
        gid_ = next((x['gene_id'] for x in xs if x.get('gene_symbol_a') == r['gene_symbol']), None)
        rels = relation(gid_) if gid_ else {}
        if rels:
            GID5, GSYM5 = gid_, r['gene_symbol']; (cat5, fac5), rel5 = sorted(rels.items())[0]; break
R.chosen['partial_retention_mouse@F5E'] = dict(mode=specP.get('mode', 'rank'), gene_id=GID5, factor=None if rel5 is None else fac5)
E_OK = rel5 is not None
if not E_OK:
    R.unavailable('F5E', 'no eligible partial-retention example; panel left empty'); GSYM5 = ''
else:
    if len(rel5['changed']) != 1 or len(rel5['retained']) != 1:
        R.flag('F5E', 'row_capacity', f'{GSYM5} {fac5}: {len(rel5["changed"])} changed / {len(rel5["retained"])} retained sites; the diagram shows the first of each')
    CH, RT = sorted(rel5['changed'])[0], sorted(rel5['retained'])[0]
feat = use('F5E', 'examples/inclusion_features.tsv', lambda r: r['gene_id'] == GID5)
iss = use('F5E', 'examples/inclusion_state_summary.tsv', lambda r: r['gene_id'] == GID5)
MINUS5 = {r.get('strand') for r in feat} != {'+'}
_coords = [int(r[k]) for r in feat if r['kind'] == 'interval' for k in ('start0', 'end0')]+([*CH, *RT] if E_OK else [])
if specP.get('mode') == 'pin' and specP.get('window_end0'):
    ANCH = int(specP['window_end0'])
else:
    import math as _m
    ANCH = (_m.ceil((max(_coords)+100)/250)*250 if MINUS5 else _m.floor((min(_coords)-100)/250)*250) if _coords else 0
XG = (lambda c: 262.0+(ANCH-c)/1000*32.0) if MINUS5 else (lambda c: 262.0+(c-ANCH)/1000*32.0)
if _coords and max(abs(ANCH-c) for c in _coords) > 5600:
    R.flag('F5E', 'axis_range', f'{GSYM5}: drawn coordinates span more than the 5.6 kb track', data_range=[0, max(abs(ANCH-c) for c in _coords)], axis=[0, 5600])
ex = [r for r in xs if r['gene_id'] == GID5]
t.text_top(260.4, 197.7, GSYM5, 9, font='BolIta')
yt = 222.0
t.line([XG(ANCH), XG(ANCH-5600 if MINUS5 else ANCH+5600)], [yt, yt], lw=1.6, color='#C9D1D6') if E_OK else None
fe = [r for r in feat if r['route'] == 'CIP' and r['kind'] == 'interval']
cuts = sorted({int(r['start0']) for r in fe} | {int(r['end0']) for r in fe}, reverse=True)
civ = [r for r in fe if E_OK and int(r['start0']) <= CH[0] and int(r['end0']) >= CH[1]]
if E_OK and len(civ) != 1:
    R.unavailable('F5E', f'{GSYM5}: CIP interval holding the changed site not unique ({len(civ)})'); E_OK = False
if not E_OK:
    cuts, civ = [], [dict(start0=0, end0=0)]; CH = RT = (0, 0)
a, b = int(civ[0]['start0']), int(civ[0]['end0'])
t.rect(min(XG(a), XG(b)), yt-4.6, abs(XG(a)-XG(b)), 9.2, fc='#DCEAF3') if E_OK else None
for c in cuts:
    t.line([XG(c)]*2, [yt-3.6, yt+3.6], lw=0.8, color=INK)
for (s0, s1), lab, filled in ([(CH, 'Changed', True), (RT, 'Constitutive', False)] if E_OK else []):
    xc = XG((s0+s1)/2)
    t.marker(xc, yt-8.6, 'v', 4.4, fc='#0072B2' if filled else 'white', ec='#0072B2', mew=0.8)
    if lab == 'Constitutive':   # 2026-10-07 terminology: right edge kept where the centred 'Retained' label ended, clear of the scale bar
        t.text_top(xc+textwidth('Retained', 7.5)/2, 227.6, lab, 7.5, ha='right', color=MUTED)
    else:
        t.text_top(xc, 227.6, lab, 7.5, ha='center', color=MUTED)
LAB5 = (specP.get('label') if specP.get('mode') == 'pin' else None) or (f'{fac5}, {cat5} {"seed" if cat5 == "miRNA" else "motif"}' if E_OK else '')
t.text_top(448.0, 198.2, LAB5, 7.5, ha='right', color=MUTED)
if 448.0-textwidth(LAB5, 7.5) <= 305.0: R.flag('F5E', 'overlap', 'factor label runs into the gene name')
sb0 = 287.0
t.line([sb0, sb0+32.0], [233.0, 233.0], lw=0.9, color=INK)
t.text_top(sb0+35.0, 229.4, '1 kb', 7.5, color=MUTED)
if E_OK and not (sb0+35.0+textwidth('1 kb', 7.5) < min(XG(a), XG(b))-2.0 and sb0 > XG((RT[0]+RT[1])/2)+textwidth('Retained', 7.5)/2+2.0):
    R.flag('F5E', 'overlap', 'scale bar collides with the interval or the retained-site label')
lib = use('F5E', 'examples/inclusion_library_states.tsv.gz', lambda r: r['gene_id'] == GID5 and r['cutoff'] in ('20', '20.0'))
R.path('F5E', 'panel_statistics/inclusion-paired-means.tsv')    # registered; 5E draws per-state library values
Y0E, Y1E = 252.0, 291.0                       # inclusion row (100 % .. 0 %)
U0, U1 = 312.0, 351.0                         # PAS-usage row (100 % .. 0 %)
ys = lambda v: Y1E-v*(Y1E-Y0E)
yu = lambda v: U1-v/100*(U1-U0)
ee5 = use('F5E', 'examples/example_endpoint_effects.tsv.gz', lambda r: r['gene_id'] == GID5 and r['setting'] == 'supported50_p08')
SC5 = dict(S='#E69F00', P='#CC79A7', M='#8C6D31'); LS5 = dict(S=(0, (3, 1.6)), P='-', M=(0, (3, 1.2, 1, 1.2)))
Ed = {}; EU5 = {}
for k, (route, title, colr) in enumerate([('CIP', 'Guo2019–CIP', '#CF4C00'), ('WOT', 'Schiebinger2019–2i', '#207450')]):
    dx = 102.0*k
    t.text_top(305.6+dx, 240.6, title, 8, ha='center', color=colr)
    fr_ = [r for r in feat if r['route'] == route and r['kind'] == 'interval'
           and int(r['start0']) <= CH[0] and int(r['end0']) >= CH[1]]
    if not E_OK or len(fr_) != 1:
        if E_OK: R.unavailable('F5E', f'{route}: interval holding the changed site not unique ({len(fr_)})')
        continue
    fid = fr_[0]['feature_id']
    st = {r['state']: r for r in iss if r['feature_id'] == fid and r['cutoff'] == '20'}
    t.line([265.5+dx]*2, [ys(1.0), ys(0)], lw=0.5, color=AXC)
    for v in (1.0, 0.5, 0):
        t.line([263.5+dx, 265.5+dx], [ys(v)]*2, lw=0.5, color=AXC)
        if k == 0: t.text_top(262.8, ys(v)-4.0, str(int(v*100)), 8, ha='right')
    xs_ = {'S': 271+dx, 'R': 286.7+dx, 'E': 302.6+dx, 'P': 318.5+dx, 'M': 340.8+dx}
    t.line([330.0+dx]*2, [ys(1.0), ys(0)], lw=0.5, color='#AAB4BB', ls=(0, (1, 1.5)))
    path = [(xs_[s_], ys(float(st[s_]['equal_library']))) for s_ in 'SREP' if s_ in st]
    t.line([p[0] for p in path], [p[1] for p in path], lw=0.8, color='#ADB9C2')
    libs = {}
    for r in lib:
        if r['feature_id'] == fid:
            libs.setdefault(r['state'], []).append(float(r['inclusion']))
    for s_ in 'SREPM':
        if s_ not in st: continue
        r = st[s_]; v = float(r['equal_library']); lo, hi = float(r['minimum']), float(r['maximum'])
        vals = sorted(libs.get(s_, []))
        if not (vals and len(vals) == int(r['n_libraries']) and abs(min(vals)-lo) < 1e-9 and abs(max(vals)-hi) < 1e-9):
            R.flag('F5E', 'consistency', f'{route} {s_}: library inclusion values disagree with the state summary')
        for q_, vv in enumerate(vals):
            off = 0.0 if len(vals) == 1 else -3.6+7.2*q_/(len(vals)-1)
            t.marker(xs_[s_]+off, ys(vv), 'o', 2.3, fc='white', ec=SC[s_], mew=0.5, zorder=3)
        t.line([xs_[s_]]*2, [ys(hi), ys(lo)], lw=0.6, color=SC[s_])
        t.marker(xs_[s_], ys(v), 'o', 4.2, fc=SC[s_], ec='white', mew=0.5, zorder=5)
        t.text_top(xs_[s_], Y1E+2.4, s_, 8, ha='center', color=SC[s_])
    Ed[route] = {s_: dict(equal_library=float(st[s_]['equal_library']), minimum=float(st[s_]['minimum']), maximum=float(st[s_]['maximum']),
                         n_libraries=int(st[s_]['n_libraries']), libraries=sorted(libs.get(s_, []))) for s_ in st}
    Ed[route]['feature_id'] = fid
    # PAS-usage line plot (style of Fig. 3D): % of family usage at or distal to each position
    rps = [r for r in ee5 if r['route'] == route and r['contrast'] == 'P_S']; rpm = [r for r in ee5 if r['route'] == route and r['contrast'] == 'P_M']
    if len(rps) != 1 or len(rpm) != 1 or json.loads(rpm[0]['PAS_ids']) != json.loads(rps[0]['PAS_ids']):
        R.unavailable('F5E', f'{route}: native P–S and P–M usage on one PAS set not available; usage plot left empty'); continue
    ids = json.loads(rps[0]['PAS_ids']); pos = [pas_pos(i_) for i_ in ids]
    mn_ = pas_strand(ids[0]) == '-'
    prox = max(pos) if mn_ else min(pos)                             # minus strand: proximal = largest coordinate
    dd_ = (lambda p_: prox-p_) if mn_ else (lambda p_: p_-prox)
    us = {'S': json.loads(rps[0]['reference_usage']), 'P': json.loads(rps[0]['target_usage']), 'M': json.loads(rpm[0]['reference_usage'])}
    KB = 5.6; XU5 = lambda d, dx=dx: 267.0+dx+d/1000/KB*78.0
    ia, ib = int(fr_[0]['start0']), int(fr_[0]['end0'])
    t.rect(min(XU5(dd_(ia)), XU5(dd_(ib))), U0, abs(XU5(dd_(ia))-XU5(dd_(ib))), U1-U0, fc='#DCEAF3')
    t.line([265.5+dx]*2, [U0, U1], lw=0.5, color=AXC); t.line([267.0+dx, XU5(KB*1000)], [U1+1.5]*2, lw=0.5, color=AXC)
    for v in (100, 50, 0):
        t.line([263.5+dx, 265.5+dx], [yu(v)]*2, lw=0.5, color=AXC)
        if k == 0: t.text_top(262.8, yu(v)-4.0, str(v), 8, ha='right')
    for v in (0, 2, 4):
        t.line([XU5(v*1000)]*2, [U1+1.5, U1+3.5], lw=0.5, color=AXC); t.text_top(XU5(v*1000), U1+4.6, str(v), 8, ha='center')
    order = sorted(range(len(pos)), key=lambda i_: dd_(pos[i_]))       # proximal first
    if max(dd_(p_) for p_ in pos) > KB*1000: R.flag('F5E', 'axis_range', f'{route}: PAS span beyond {KB} kb', data_range=[0, max(dd_(p_) for p_ in pos)], axis=[0, KB*1000])
    for st_ in ('S', 'M', 'P'):
        u = np.array(us[st_], float); u = u/u.sum(); cur = 100.0; xs, yv = [XU5(0)], [yu(100)]
        for i_ in order:
            d = min(dd_(pos[i_]), KB*1000); xs += [XU5(d), XU5(d)]; yv += [yu(cur), yu(max(cur-100*u[i_], 0))]; cur -= 100*u[i_]
        xs.append(XU5(KB*1000)); yv.append(yu(max(cur, 0)))
        t.line(xs, yv, lw=1.0, color=SC5[st_], ls=LS5[st_])
        EU5[f'{route}|{st_}'] = dict(usage=[round(float(x), 4) for x in u], distance_nt=[dd_(p_) for p_ in pos])
t.vtext(244.7, (Y0E+Y1E)/2, 'I (%)', 8)
t.vtext(244.7, (U0+U1)/2, 'Distal (%)', 8)
t.text_top(356.6, U1+14.0, 'Distance from proximal PAS (kb)', 8, ha='center')
x = 318.0
for st_ in ('S', 'P', 'M'):
    t.line([x, x+10.0], [306.0, 306.0], lw=1.0, color=SC5[st_], ls=LS5[st_]); t.text_top(x+12.5, 302.4, st_, 7.5); x += 26.0
t.rect(x+4.0, 303.4, 6.0, 5.0, fc='#DCEAF3'); t.text_top(x+12.5, 302.4, 'interval', 7.5)
assert x+12.5+textwidth('interval', 7.5) < 449.0 and U1+14.0+8.0 < 376.8
ver['E_usage'] = EU5
ver['E'] = dict(changed=CH, retained=RT, interval=(a, b), states=Ed, pas_cuts_CIP=cuts, exact_rows=len(ex))

# ------------------------------------------------------------------ F  (2026-10-07: Riok3, cross-workflow example)
# Riok3 is not a low-shift screen hit (|ΔDWUI| > .01 in both routes with polyApipe); it is pinned here because its terminal-family
# redistribution and interval inclusion loss are supported with both PAS workflows. Its native, family and inclusion rows are not
# in this bundle's example/geometry tables (those cover the low-shift cohort only) and come from derive.py, which reproduces the
# bundle's Trip12 rows exactly with the same code paths (receipt.json, checks).
import hashlib as _hl
DER = Path(__file__).resolve().parent/'tables'
_dsha = {p.name: _hl.sha256(p.read_bytes()).hexdigest() for p in sorted(DER.glob('*.tsv'))}
GIDF, GSYMF = 'ENSMUSG00000024404', 'Riok3'
IV889 = (12288586, 12289475)
R.chosen['low_shift@F5F'] = dict(mode='pin', chosen=dict(gene_id=GIDF, gene_symbol=GSYMF), note='cross-workflow example; not a low-shift screen hit')
R.derive('F5F', f'Riok3 native endpoint, terminal-family decomposition and paired inclusion from derived tables {json.dumps(_dsha)} (notes/2026-10-07-riok3-example/derive.py; Trip12 reproduction exact)')
R.derive('F5F', 'shared interval fixed to chr18:[12288586,12289475) (889 nt), the interval supported by both PAS workflows; the polyApipe strict interval [12288550,12289754) has the same complete-retention mask [0,1]')
use('F5F', 'examples/example_endpoint_effects.tsv.gz', lambda r: False); use('F5F', 'examples/inclusion_paired_libraries.tsv.gz', lambda r: False)
use('F5F', 'geometry/axis-decomposition.tsv', lambda r: False); use('F5F', 'candidate_screen/recurrent_region_pairs.tsv.gz', lambda r: False)
use('F5F', 'panel_statistics/inclusion-paired-means.tsv', lambda r: False)
ee = tsv(DER/'example_endpoint_effects.tsv')
ip = [x for x in tsv(DER/'inclusion_paired_libraries.tsv') if (int(x['start0']), int(x['end0'])) == IV889]
pmF = [x for x in tsv(DER/'inclusion-paired-means.tsv') if x['cutoff'] in ('20', '20.0') and f'|{IV889[0]}-{IV889[1]}|' in x['feature_id']]
famF = {r['route']: r for r in tsv(DER/'axis-decomposition.tsv')}
inv = use('F5F', 'site_statistics/site_factor_inventory.tsv.gz', lambda r: r['gene_id'] == GIDF)
PSC_SEEDS = {r['seed'] for r in R.rows('F5F', 'external_mirna/external_seed_expression.tsv',
                                       lambda r: r['species'] == 'mouse' and r['expression_group'] == 'PSC')}
_pex = {}
for r_ in R.rows('F5F', 'rbp_expression/state_expression.tsv', lambda r: r['species'] == 'mouse' and r['state'] == 'P' and r['route'] in ('CIP', 'WOT')):
    _pex.setdefault(r_['factor'], set()).add(r_['route']) if r_['expressed_any_2CPM'] == 'True' else None
P_RBP = {f for f, rs in _pex.items() if rs == {'CIP', 'WOT'}}
RUG_ROWS = [
    ('miRNA', 'miRNA', lambda x: x['factor'] in PSC_SEEDS, 412.0, INK, 0.6),
    ('RBP_direct', 'RBP', lambda x: x['provenance'] == 'RBP', 421.5, '#8C6D31', 0.6),
    ('RBP_transfer', 'RBP', lambda x: x['provenance'] == 'RBP_transfer' and any(n in P_RBP for n in x['names'].split(';')), 431.0, '#C9B68F', 0.5)]
R.derive('F5F', 'transferred-motif RBP row: provenance RBP_transfer and any named RBP expressed >= 2 CPM in P of both CIP and WOT (rbp_expression/state_expression)')
R.derive('F5F', 'site rug limited to PSC-enriched miRNA seeds (external_seed_expression, mouse, expression_group PSC) and same-species RBP motifs (provenance RBP)')
RUG_LABEL = {'miRNA': 'PSC miRNA', 'RBP_direct': 'Direct RBP', 'RBP_transfer': 'Transferred RBP'}
RUG_SHORT = {'RBP_transfer': 'Transf. RBP'}
RUG_FS = 7.5
IVS = [IV889]
inany = lambda a, b: any(a >= s_ and b <= e_ for s_, e_ in IVS)
R.derive('F5F', 'family PAS usage renormalised within the terminal family of geometry/axis-decomposition (GAP-D06)')
t.text_top(24.9, 389.6, GSYMF, 9, font='BolIta')
_x = 92.0
t.line([_x, _x+10.0], [394.0, 394.0], lw=1.0, color='#8A959C', ls=(0, (3, 1.6))); t.text_top(_x+12.5, 390.0, 'M', 8); _x += 12.5+textwidth('M', 8)+9.0
t.line([_x, _x+10.0], [394.0, 394.0], lw=1.2, color=INK); t.text_top(_x+12.5, 390.0, 'P', 8); _x += 12.5+textwidth('P', 8)+12.0
_sl = 'Shared interval'
t.rect(_x, 391.4, 7.0, 6.0, fc='#E4DFEE'); t.text_top(_x+9.5, 390.0, _sl, 8); _x += 9.5+textwidth(_sl, 8)+12.0
for _c, _n in [('#CF4C00', 'Guo2019–CIP'), ('#207450', 'Schiebinger2019–2i')]:
    t.rect(_x, 391.4, 6.0, 6.0, fc=_c); t.text_top(_x+8.5, 390.0, _n, 8); _x += 8.5+textwidth(_n, 8)+12.0
assert _x < 452.0, _x
nsite = {}
for route in ('CIP', 'WOT'):
    rr = [r for r in inv if r['contrast'] == 'MP' and r['route'] == route and r['strict_changed'] == 'True' and inany(int(r['start0']), int(r['end0']))]
    nsite[route] = {c: sum(1 for r in rr if r['category'] == c) for c in ('miRNA', 'RBP')}
    nsite[route]['coords'] = {c: len({(r['start0'], r['end0']) for r in rr if r['category'] == c}) for c in ('miRNA', 'RBP')}
Fd = {}
for k, (route, title, colr) in enumerate([('CIP', 'Guo2019–CIP', '#CF4C00'), ('WOT', 'Schiebinger2019–2i', '#207450')]):
    dx = 230.2*k
    r = [r for r in ee if r['route'] == route]
    if len(r) != 1 or route not in famF:
        R.unavailable('F5F', f'{route}: native P–M usage or terminal family not available; column left empty'); continue
    r = r[0]
    ids = json.loads(r['PAS_ids']); pos = [pas_pos(i_) for i_ in ids]; strand = pas_strand(ids[0])
    tu = json.loads(r['target_usage']); ru = json.loads(r['reference_usage'])
    fs_, fe_ = int(float(famF[route]['family_start0'])), int(float(famF[route]['family_end0']))
    plusF = strand == '+'
    fam = sorted([i for i, p in enumerate(pos) if fs_ <= p <= fe_], key=lambda i: pos[i] if plusF else -pos[i])
    if len(fam) != int(float(famF[route]['terminal_native_pas'])):
        R.flag('F5F', 'consistency', f'{route}: {len(fam)} PAS inside family bounds != terminal_native_pas {famF[route]["terminal_native_pas"]}')
    if not fam:
        R.unavailable('F5F', f'{route}: no PAS inside the terminal family'); continue
    base = pos[fam[0]]
    DF = (lambda p: p-base) if plusF else (lambda p: base-p)
    FMAX = 2200.0
    R.axis('F5F', f'{route} complete family distance', [DF(pos[j]) for j in fam], 0, FMAX)
    XF = lambda nt, dx=dx: 34.6+dx+nt*0.0365*2500/FMAX
    YF = lambda v: 499.2-v*0.528
    t.rect(6.0+dx, 402.0, 212.0, 3.0, fc=colr)
    for iv in IVS:
        xa_, xb_ = sorted((XF(DF(iv[0])), XF(DF(iv[1]))))
        t.rect(xa_, YF(100), xb_-xa_, YF(0)-YF(100), fc='#E4DFEE')
    xl_ = min(XF(DF(v_)) for iv in IVS for v_ in iv); xr_ = max(XF(DF(v_)) for iv in IVS for v_ in iv)
    HDR = [f"MPRO {float(r['MPRO']):.3f}, n = {r['n_pairs']}", 'DWUI M → P',
           f"{float(r['reference_DWUI']):.3f} → {float(r['target_DWUI']):.3f}"]
    for k_, h_ in enumerate(HDR):
        t.text_top(216.0+dx, 408.8+9.6*k_, h_, 8, ha='right')
    OBST = [(216.0+dx-textwidth(h_, 8), 408.8+9.6*k_, 408.8+9.6*k_+8.0) for k_, h_ in enumerate(HDR)]
    OBST += [(153.2+dx-textwidth('100', 8), 440.2, 448.2)]
    rug = {}
    for key, cat, keep_, yr, cc, lw_ in RUG_ROWS:
        rr = [x for x in inv if x['contrast'] == 'MP' and x['route'] == route and x['strict_changed'] == 'True' and x['category'] == cat
              and inany(int(x['start0']), int(x['end0'])) and keep_(x)]
        ds = sorted({DF((int(x['start0'])+int(x['end0']))/2) for x in rr})
        if not ds: R.unavailable('F5F', f'{route} {key}: no site of this row inside the interval')
        for d in ds:
            t.line([XF(d)]*2, [yr, yr+4.2], lw=lw_, color=cc)
        lab_ = f'{RUG_LABEL[key]} ({len(ds)})'
        if xl_-2.0-textwidth(lab_, RUG_FS) > 30.6+dx+6.0:
            t.text_top(xl_-2.0, yr-2.1, lab_, RUG_FS, ha='right', color=MUTED)
        else:
            y0_, y1_ = yr-2.1, yr-2.1+RUG_FS
            lim_ = min([x_ for x_, a_, b_ in OBST if a_ < y1_+1.0 and b_ > y0_-1.0] + [218.0+dx]) - 4.0
            if xr_+2.5+textwidth(lab_, RUG_FS) >= lim_ and key in RUG_SHORT:
                lab_ = f'{RUG_SHORT[key]} ({len(ds)})'
            t.text_top(xr_+2.5, yr-2.1, lab_, RUG_FS, ha='left', color=MUTED)
            assert xr_+2.5+textwidth(lab_, RUG_FS) < lim_, lab_
        rug[key] = dict(rule=RUG_LABEL[key], annotations=len(rr), distinct_sites=len(ds), factors=sorted({x['factor'] for x in rr}),
                        named=sorted({n for x in rr for n in x['names'].split(';')}) if key == 'RBP_transfer' else None)
    for usage, style, col, lw in [(ru, (0, (3, 1.6)), '#8A959C', 1.0), (tu, '-', INK, 1.2)]:
        u = np.array([usage[i] for i in fam]); u = u/u.sum(); cur = 1.0; pts = [(XF(0), YF(100))]
        for i, j in enumerate(fam):
            d = DF(pos[j]); pts.append((XF(d), YF(cur*100))); cur -= u[i]; pts.append((XF(d), YF(max(cur, 0)*100)))
        pts.append((XF(FMAX), YF(max(cur, 0)*100)))
        t.line([p[0] for p in pts], [p[1] for p in pts], lw=lw, color=col, ls=style)
    t.line([30.6+dx, 30.6+dx], [YF(100), YF(0)], lw=0.5, color=AXC); t.line([30.6+dx, XF(FMAX)], [501.6, 501.6], lw=0.5, color=AXC)
    for v in (100, 50, 0):
        t.line([28.6+dx, 30.6+dx], [YF(v)]*2, lw=0.5, color=AXC); t.text_top(27.3+dx, YF(v)-4.0, str(v), 8, ha='right')
    for v in (0, 1000, 2000):
        t.line([XF(v)]*2, [501.6, 503.6], lw=0.5, color=AXC); t.text_top(XF(v), 504.0, str(v), 8, ha='center')
    t.text_top(XF(1000), 513.0, 'Distance (nt)', 8, ha='center')
    t.vtext(9.0+dx, 470.0, 'Family PAS usage (%)', 7.5)
    lo_, hi_ = 0.0, 1.0
    YI = lambda v: 496.2-(min(max(v, lo_), hi_)-lo_)/(hi_-lo_)*52.0
    xm, xp = 168.0+dx, 205.8+dx
    lib = [x for x in ip if x['route'] == route and x['cutoff'] in ('20', '20.0')]
    refs = np.array([float(x['reference_inclusion']) for x in lib]); tars = np.array([float(x['target_inclusion']) for x in lib])
    pm_ = [m for m in pmF if m['route'] == route]
    incl = []
    if not lib or len(pm_) != 1:
        R.unavailable('F5F', f'{route}: inclusion pairs / paired mean not available')
    else:
        pm_ = pm_[0]
        R.axis('F5F', f'{route} interval inclusion', list(refs)+list(tars), lo_, hi_)
        for a_, b_ in zip(refs, tars):
            t.line([xm, xp], [YI(a_), YI(b_)], lw=0.4, color=colr, alpha=0.55)
        t.line([xm, xp], [YI(float(pm_['reference_mean'])), YI(float(pm_['target_mean']))], lw=1.6, color=colr)
        t.text_top(xm, 504.0, 'M', 8, ha='center'); t.text_top(xp, 504.0, 'P', 8, ha='center')
        if int(float(pm_['n'])) != len(lib): R.flag('F5F', 'consistency', f'{route}: paired-mean n {pm_["n"]} != {len(lib)} libraries')
        incl.append(dict(interval=1, start0=IV889[0], end0=IV889[1], length_nt=IV889[1]-IV889[0], n=int(float(pm_['n'])),
                         source_feature_id=pm_['feature_id'], reference_mean=float(pm_['reference_mean']), target_mean=float(pm_['target_mean']),
                         delta_I_pp=float(pm_['delta_mean'])*100))
    t.line([156.5+dx]*2, [YI(1.0), YI(lo_)], lw=0.5, color=AXC)
    for v, s in [(1.0, '100'), (0.5, '50'), (0.0, '0')]:
        t.line([154.5+dx, 156.5+dx], [YI(v)]*2, lw=0.5, color=AXC); t.text_top(153.2+dx, YI(v)-4.0, s, 8, ha='right')
    if incl:
        lab = 'ΔI ' + ' / '.join(num(x['delta_I_pp'], 1) for x in incl) + ' pp; n=' + '/'.join(str(x['n']) for x in incl)
        fs_ = 8 if textwidth(lab, 8) < 216.0-112.0 else 7.5
        t.text_top(216.0+dx, 515.2, lab, fs_, ha='right')
        assert 216.0-textwidth(lab, fs_) > 108.0, lab
    t.vtext(135.0+dx, 470.0, 'Interval\ninclusion (%)', 7.5)
    Fd[route] = dict(MPRO=float(r['MPRO']), n_native_pairs=int(r['n_pairs']), DWUI_M=float(r['reference_DWUI']), DWUI_P=float(r['target_DWUI']),
                     delta_DWUI=float(r['delta_DWUI']), family_shift_nt=float(famF[route]['terminal_shift_nt']),
                     within_family_DWUI=float(famF[route]['within_family_DWUI']), nonterminal_DWUI=float(famF[route]['nonterminal_DWUI']),
                     intervals=incl, sites=nsite[route], site_rug=rug)
ver['F'] = Fd

t.save(R.out/'overlay.pdf'); overlay(pdf, R.out/'overlay.pdf', 0, 0); pdf.save(OUT)
render_png(OUT, R.out/'Figure.5.draft.png', 150)
ver['output_sha256'] = sha(OUT)
R.finish(5, OUT, json.loads(json.dumps(ver, default=str)))
