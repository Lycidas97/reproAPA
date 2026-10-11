"""Supplementary Figure D14 (printed number from supdraw.SNUM) — functional enrichment of the profile classes (R2),
polyApipe v2 primary.

a  test positions: class over-representation (ORA) and ranked-set tests per measure and PAS layer
b  q < .05 positions per test family (class ORA; ranked APA and GEX in the supported and expanded layers)
c  ranked APA positives by route in the two layers
d  complete corrected-q distributions: class ORA by species; ranked tests by measure and layer
e  ranked APA effect, supported against expanded layer, all terms
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D14'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 87.4, 126.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
rk = D.df('D14A', 'profile_gex/rank_enrichment_all_tests.tsv.gz')
ora = D.df('D14A', 'profile_gex/enrichment_all_tests.tsv.gz')
FAM = [('ORA', None, None, 'Class ORA'), ('rank', 'APA', 'supported50', 'APA, supported'), ('rank', 'APA', 'expanded_primary', 'APA, expanded'),
       ('rank', 'GEX', 'supported50', 'GEX, supported'), ('rank', 'GEX', 'expanded_primary', 'GEX, expanded')]
ROUTES = ['GSE118258', 'CIP', 'WOT']


def fam_rows(kind, meas, layer):
    if kind == 'ORA': return ora.q_value
    q = rk[(rk.measure == meas) & (rk.layer == layer)]
    return q.q_both_databases


labw = max(textwidth(l, 7.0) for *_, l in FAM)
# ================================================================ a  test positions
WA = 146.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
ax_ = Lin(0, 16000, x+labw+4, x+w-24)
F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Test positions', 7.5, ha='center')
Av = {}
for i, (kind, meas, layer, lab) in enumerate(FAM):
    yy = y+20.0+i*13.0
    n = len(fam_rows(kind, meas, layer))
    F.text(x+0.4, yy-3.4, lab, 7.0)
    F.rect(ax_(0), yy-3.6, ax_(n)-ax_(0), 7.2, fc=NULLBAR if kind == 'rank' else ORD3[2])
    F.text(ax_(n)+2.0, yy-3.4, f'{n:,}', 6.8, color=MUTED)
    Av[lab] = n
F.xaxis(ax_, y+20.0+4*13.0+5.0, [0, 8000, 16000], ['0', '8k', '16k'])
F.values['a'] = Av

# ================================================================ b  positives per family
WB = 146.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
bx = Lin(0, 1500, x+labw+4, x+w-7.5)
F.text((bx.p0+bx.p1)/2, y+3.0, 'q < .05 positions', 7.5, ha='center')
Bv = {}
for i, (kind, meas, layer, lab) in enumerate(FAM):
    yy = y+20.0+i*13.0
    q = fam_rows(kind, meas, layer); n = int((q < 0.05).sum())
    F.text(x+0.4, yy-3.4, lab, 7.0)
    F.line([bx(0), bx(n)], [yy, yy], lw='bold', color=NULLBAR)
    F.marker(bx(n), yy, 'o', 'mark', fc=INK if meas != 'GEX' else 'white', ec=INK)
    F.text(bx(n)+4.0, yy-3.4, f'{n:,}', 6.8)
    Bv[lab] = dict(q05=n, tests=int(len(q)), min_q=float(q.min()))
F.xaxis(bx, y+20.0+4*13.0+5.0, [0, 500, 1000, 1500], ["0", "500", "1k", "1.5k"])
F.values['b'] = Bv

# ================================================================ c  APA positives by route
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
cx = Lin(0, 45, x+16, x+w-3)
F.text((cx.p0+cx.p1)/2, y+3.0, 'APA q < .05, by route', 7.5, ha='center')
Cv = {}
for i, r in enumerate(ROUTES):
    yy = y+21.0+i*17.6
    F.marker(x+6.0, yy, STUDY[r][2], 'mark', fc=STUDY[r][1])
    Cv[r] = {}
    for k, (layer, cc) in enumerate((('supported50', ORD3[0]), ('expanded_primary', ORD3[2]))):
        q = rk[(rk.route == r) & (rk.measure == 'APA') & (rk.layer == layer)]
        n = int((q.q_both_databases < 0.05).sum())
        yk = yy+(k-0.5)*6.2
        F.rect(cx(0), yk-2.6, max(cx(n)-cx(0), 0.6), 5.2, fc=cc)
        F.text(cx(n)+2.0, yk-3.4, str(n), 6.8)
        Cv[r][layer] = dict(q05=n, tests=int(len(q)))
F.xaxis(cx, y+21.0+2*17.6+8.0, [0, 20, 40])
F.key_row(None, y+h-7.4, [('block', ORD3[0], 'Supported'), ('block', ORD3[2], 'Expanded')], gap=5.0, center=x+w/2)
F.values['c'] = Cv

# ================================================================ d  q distributions
ES = 80.0
WD = RR-X0-GAP-(ES+34)
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
cy = Lin(0, 1, y+90, y+16)
dw_ = (w-26-20-4)/2
d1 = Lin(0, 1, x+26, x+26+dw_); d2 = Lin(0, 1, x+26+dw_+20, x+w-4)
F.lineplot(d1, cy, 'ORA q'); F.lineplot(d2, cy, 'rank q')
F.text((d1.p0+d1.p1)/2, y+3.0, 'Class ORA', 7.5, ha='center'); F.text((d2.p0+d2.p1)/2, y+3.0, 'Ranked sets', 7.5, ha='center')
F.yaxis(cy, d1.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Cumulative fraction', title_y=(cy.p0+cy.p1)/2+4)
F.yaxis(cy, d2.p0, [0, 0.5, 1], ['0', '.5', '1'])


def ecdf(s_, v, **kw):
    v = np.sort(np.asarray(v, float)); n = len(v)
    xs_ = np.r_[0, np.repeat(v, 2), 1]; ys_ = np.r_[0, 0, np.repeat(np.arange(1, n+1)/n, 2)][:len(xs_)]
    F.line([s_(t) for t in xs_], [cy(t) for t in ys_], **kw)


Dv = {}
for sp, col in (('human', STUDY['GSE118258'][1]), ('mouse', STUDY['WOT'][1])):
    v = ora[ora.species == sp].q_value; ecdf(d1, v, lw='data', color=col); Dv[f'ORA|{sp}'] = dict(n=int(len(v)), min_q=float(v.min()))
for meas, layer, col, ls in (('APA', 'supported50', INK, None), ('APA', 'expanded_primary', INK, 'dash'),
                             ('GEX', 'supported50', MUTED, None), ('GEX', 'expanded_primary', MUTED, 'dash')):
    v = rk[(rk.measure == meas) & (rk.layer == layer)].q_both_databases
    ecdf(d2, v, lw='data', color=col, ls=ls); Dv[f'{meas}|{layer}'] = dict(n=int(len(v)), q05=int((v < 0.05).sum()))
for s_ in (d1, d2):
    F.line([s_(0.05)]*2, [cy.p0, cy.p1], lw='ref', color=AXC, ls='ref')
    F.xaxis(s_, cy.p0+1.0, [0, 0.5, 1], ['0', '.5', '1'])
F.text((d1.p0+d2.p1)/2, cy.p0+12.6, 'Corrected q', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('line', (STUDY['GSE118258'][1], None, 'data'), 'Human'), ('line', (STUDY['WOT'][1], None, 'data'), 'Mouse')],
          gap=6.0, center=(d1.p0+d1.p1)/2)
K2 = [('line', (INK, None, 'data'), 'APA'), ('line', (MUTED, None, 'data'), 'GEX'),
      ('line', (INK, None, 'data'), 'Supported'), ('line', (INK, 'dash', 'data'), 'Expanded')]
F.key_row(None, y+h-7.4, K2, gap=5.0, center=min((d2.p0+d2.p1)/2, x+w-0.4-F.key_width(K2, gap=5.0)/2))
F.values['d'] = Dv

# ================================================================ e  supported vs expanded APA effect
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, RR-(X0+WD+GAP), H2)
pv = rk[rk.measure == 'APA'].pivot_table(index=['route', 'database', 'term_id'], columns='layer', values='rank_effect').dropna().reset_index()
ex = Lin(-0.6, 0.6, x+30, x+30+ES); ey = Lin(-0.6, 0.6, y+16+ES, y+16)
F.text((ex.p0+ex.p1)/2, y+3.0, 'APA rank effect', 7.5, ha='center')
F.line([ex.p0, ex.p1], [ey.p0, ey.p1], lw='ref', color=AXC, ls='ref')
F.yaxis(ey, ex.p0, [-0.4, 0, 0.4], ['−.4', '0', '.4'], 'Expanded layer')
for r in ROUTES:
    z = pv[pv.route == r]
    F.points([ex.clip(v) for v in z.supported50], [ey.clip(v) for v in z.expanded_primary], color=STUDY[r][1], alpha='dense')
F.xaxis(ex, ey.p0+1.0, [-0.4, 0, 0.4], ['−.4', '0', '.4'], 'Supported layer')
F.key_row(None, y+h-7.4, [('marker', ('o', STUDY[r][1]), STUDY[r][0].replace('Schiebinger2019–', '').replace('Guo2019–', '').replace('2020', ''))
                          for r in ROUTES], gap=5.0, center=x+w/2)
F.values['e'] = dict(terms=int(len(pv)), r=float(np.corrcoef(pv.supported50, pv.expanded_primary)[0, 1]))
F.notes.append('b–d: ranked-set q = q_both_databases (GO and KEGG corrected together), as in the statistics counts; class ORA min q = '
               f'{ora.q_value.min():.2f}')
F.save(OUT, D)
