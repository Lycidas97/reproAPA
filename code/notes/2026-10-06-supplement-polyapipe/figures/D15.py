"""Supplementary Figure D15 (printed number from supdraw.SNUM) — GO/KEGG for the recurrence sets (R3), polyApipe v2 primary.

a  nested gene backgrounds: G0 expressed (annotated), G1 APA-testable, G2 recurrent, per species and contrast (GO annotation)
b  corrected tests with q < .05 for G1 vs G0, G2 vs G0 and G2 vs G1 (L3 full null, GO and KEGG)
c  the 14 G1-vs-G0 terms: expected -> observed genes, matched O/E and q; they describe measurability, not the recurrent set
d  conditional-null label mobility (expected changed fraction) under the five null models
e  detection power of the corrected G2 tests against an injected odds ratio
"""
import sys, math, textwrap
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D15'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 79.0, 170.0, 140.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname(FID))
fam = D.df('D15A', 'functional/functional_families.tsv')
GR = [('human', 'P_S', 'Human P–S'), ('human', 'P_M', 'Human P–M'), ('mouse', 'P_S', 'Mouse P–S'), ('mouse', 'P_M', 'Mouse P–M')]
GCOL = {'G0': ORD3[2], 'G1': ORD3[1], 'G2': ORD3[0]}

# ================================================================ a  nested backgrounds
WA = 224.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
ax_ = Lin(0, 15000, x+50, x+w-7)
F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Annotated genes (GO)', 7.5, ha='center')
Av = {}
for i, (sp, ct, gl) in enumerate(GR):
    yy = y+21.0+i*9.4
    a1 = fam[(fam.database == 'GO') & (fam.analysis_id == f'{sp}|{ct}|G1_vs_G0|all')].iloc[0]
    a2 = fam[(fam.database == 'GO') & (fam.analysis_id == f'{sp}|{ct}|G2_vs_G0|all')].iloc[0]
    g0, g1, g2 = int(a1.N), int(a1.n), int(a2.n)
    F.text(x+0.4, yy-3.4, gl, 7.0)
    F.rect(ax_(0), yy-3.4, ax_(g0)-ax_(0), 6.8, fc=GCOL['G0'])
    F.rect(ax_(0), yy-3.4, ax_(g1)-ax_(0), 6.8, fc=GCOL['G1'])
    F.rect(ax_(0), yy-3.4, ax_(g2)-ax_(0), 6.8, fc=GCOL['G2'])
    Av[gl] = dict(G0=g0, G1=g1, G2=g2)
F.xaxis(ax_, y+21.0+3*9.4+4.6, [0, 5000, 10000, 15000], ['0', '5k', '10k', '15k'])
F.key_row(None, y+h-7.4, [('block', GCOL['G0'], 'G0 expressed'), ('block', GCOL['G1'], 'G1 APA-testable'), ('block', GCOL['G2'], 'G2 recurrent')],
          gap=6.0, center=x+w/2)
F.values['a'] = Av

# ================================================================ b  corrected positives
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
ex_ = {db: D.df('D15B', f'functional/{db.lower()}_exact.tsv.gz') for db in ('GO', 'KEGG')}
bx = Lin(0, 16, x+36, x+w-30)
F.text((bx.p0+bx.p1)/2, y+3.0, 'q < .05 positions', 7.5, ha='center')
F.text(x+w, y+3.0, 'Tests', 7.0, ha='right', color=MUTED)
DBC = {'GO': '#2D6A96', 'KEGG': '#9CC3DD'}
Bv = {}
for i, (cmp_, cl) in enumerate([('G1_vs_G0', 'G1 vs G0'), ('G2_vs_G0', 'G2 vs G0'), ('G2_vs_G1', 'G2 vs G1')]):
    yy = y+21.0+i*13.6
    F.text(x+0.4, yy-3.4, cl, 7.0)
    acc = 0; tests = 0; Bv[cl] = {}
    for db in ('GO', 'KEGG'):
        q = ex_[db]; q = q[(q.model == 'L3_full') & (q.comparison == cmp_)]
        n = int((q.q_database_model < 0.05).sum()); tests += len(q)
        if n: F.rect(bx(acc), yy-3.6, bx(acc+n)-bx(acc), 7.2, fc=DBC[db])
        acc += n; Bv[cl][db] = dict(q05=n, tests=int(len(q)))
    F.text(bx(acc)+2.5, yy-3.4, str(acc), 7.0)
    F.text(x+w, yy-3.4, f'{tests:,}', 6.8, ha='right', color=MUTED)
F.xaxis(bx, y+21.0+2*13.6+6.4, [0, 8, 16])
F.key_row(None, y+h-7.4, [('block', DBC['GO'], 'GO'), ('block', DBC['KEGG'], 'KEGG')], gap=6.0, center=(bx.p0+bx.p1)/2)
F.values['b'] = Bv

# ================================================================ c  14 terms
x, y, w, h = F.panel('c', X0, Y2, RR-X0, H2)
td = D.df('D15C', 'supplement_statistics/endpoint-GO-KEGG-display.tsv').sort_values('display_order')
SCOL = {'P_M': SC['M'], 'P_S': SC['S']}                    # supplement rule: contrast = colour of the reference state
SMK = {'human': 'o', 'mouse': 's'}                           # species = marker shape
TL = {r.term_name: r.term_name+(' (KEGG)' if r.database == 'KEGG' else '') for r in td.itertuples()}
labw = max(textwidth(t, 7.0) for t in TL.values())+2.0
c1 = Lin(0, 240, x+labw+12, x+labw+12+100)
c2 = Lin(1.0, 2.0, c1.p1+22, x+w-34)
F.text((c1.p0+c1.p1)/2, y+3.0, 'Expected → observed genes', 7.5, ha='center')
F.text((c2.p0+c2.p1)/2, y+3.0, 'Matched O/E', 7.5, ha='center')
F.text(x+w, y+3.0, 'q', 7.0, ha='right', color=MUTED)
rh = 9.6; r0 = y+21.0
Cv = []
for i, r in enumerate(td.itertuples()):
    yy = r0+rh*i
    if i % 2 == 0: F.rect(x, yy-rh/2, w, rh, fc=BAND)
    col = SCOL[r.contrast]; mk = SMK[r.species]
    nm = TL[r.term_name]
    F.text(x+labw, yy-3.4, nm, 7.0, ha='right')
    F.line([c1(r.expected), c1(r.hits)], [yy, yy], lw='thin', color=col)
    F.marker(c1(r.expected), yy, mk, 'small', fc='white', ec=col)
    F.marker(c1(r.hits), yy, mk, 'small', fc=col)
    F.line([c2(1.0), c2(r.matched_fold)], [yy, yy], lw='thin', color=col)
    F.marker(c2(r.matched_fold), yy, mk, 'small', fc=col)
    F.text(x+w, yy-3.4, qfmt(r.q_database_model), 6.8, ha='right')
    Cv.append(dict(term=r.term_name, db=r.database, species=r.species, contrast=r.contrast, hits=int(r.hits),
                   expected=float(r.expected), matched_OE=float(r.matched_fold), q=float(r.q_database_model)))
ab = r0+rh*len(td)-rh/2+1.0
F.line([c2(1.0)]*2, [r0-rh/2, ab-1.0], lw='ref', color=AXC, ls='ref')
F.xaxis(c1, ab, [0, 80, 160, 240]); F.xaxis(c2, ab, [1, 1.5, 2])
F.key_row(None, y+h-7.4, [('marker', ('o', 'white', INK), 'Expected'), ('marker', ('o', INK), 'Observed'),
                          ('block', SC['S'], 'P–S'), ('block', SC['M'], 'P–M'), ('marker', ('o', INK), 'Human'), ('marker', ('s', INK), 'Mouse')],
          gap=7.0, center=x+w/2)
F.values['c'] = Cv
F.notes.append('c: all terms are G1 (measurable) vs G0 (eligible), L3 full null, q_database_model; none are G2 (recurrent) discoveries')

# ================================================================ d  mobility
WD = 224.0
x, y, w, h = F.panel('d', X0, Y3, WD, H3)
mo = D.df('D15D', 'functional/functional_mobility.tsv')
MODELS = [('L0_eligibility', 'L0', NULLW, 'dot'), ('L1_measurement', 'L1', MUTED, 'dash'), ('L2_baseline', 'L2 baseline', AXC, None),
          ('L2_PAS_first', 'L2 PAS-first', '#9CC3DD', None), ('L3_full', 'L3 full', INK, None)]
dx = Lin(0, 0.4, x+30, x+w-3); dy = Lin(0, 1, y+110, y+12)
F.lineplot(dx, dy, 'mobility')
F.text((dx.p0+dx.p1)/2, y+3.0, f'{mo.analysis_id.nunique()} analyses', 7.5, ha='center')
F.yaxis(dy, dx.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Cumulative fraction', title_y=(dy.p0+dy.p1)/2+6, title_x=x+5.4)
Dv = {}
for m_, ml, col, ls in MODELS:
    v = np.sort(mo[mo.model == m_].expected_changed_fraction.to_numpy()); n = len(v)
    xs_ = np.r_[0, np.repeat(v, 2)]; ys_ = np.r_[0, 0, np.repeat(np.arange(1, n+1)/n, 2)][:len(xs_)]
    F.line([dx.clip(t) for t in xs_], [dy(t) for t in ys_], lw='data' if m_ == 'L3_full' else 'thin', color=col, ls=ls)
    Dv[m_] = dict(n=int(n), median=float(np.median(v)))
F.xaxis(dx, dy.p0+1.0, [0, 0.1, 0.2, 0.3, 0.4], ['0', '.1', '.2', '.3', '.4'], 'Expected changed-label fraction')
F.key_row(None, y+h-7.4, [('line', (c, l_, 'data' if m == 'L3_full' else 'thin'), ml) for m, ml, c, l_ in MODELS], gap=5.0, center=x+w/2)
F.values['d'] = Dv

# ================================================================ e  power
x, y, w, h = F.panel('e', X0+WD+GAP, Y3, RR-(X0+WD+GAP), H3)
pw_ = D.df('D15E', 'functional/functional_power.tsv.gz')
pw_ = pw_[(pw_.model == 'L3_full') & pw_.analysis_id.str.contains('G2_vs')]
px = Lin(1, 8, x+24, x+w-3); py = Lin(0, 1, y+110, y+12)
F.lineplot(px, py, 'power')
F.text((px.p0+px.p1)/2, y+3.0, 'G2 tests, L3 full null', 7.5, ha='center')
F.yaxis(py, px.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Detection probability', title_y=(py.p0+py.p1)/2+6, title_x=x+5.4)
Ev = {}
for db, col in (('GO', DBC['GO']), ('KEGG', DBC['KEGG'])):              # same database colours as b
    for th, ls in (('nominal', None), ('bonferroni_terms', 'dash')):
        q = pw_[(pw_.database == db) & (pw_.threshold == th)].groupby('odds_ratio').rejection_probability.median()
        F.line([px(v) for v in q.index], [py(v) for v in q.values], lw='data', color=col, ls=ls)
        for o_, v in q.items(): F.marker(px(o_), py(v), 'o', 'small', fc=col)
        Ev[f'{db}|{th}'] = {int(k): float(v) for k, v in q.items()}
F.xaxis(px, py.p0+1.0, [1, 2, 4, 8], title='Injected odds ratio')
F.key_row(None, y+h-7.4, [('line', (DBC['GO'], None, 'data'), 'GO'), ('line', (DBC['KEGG'], None, 'data'), 'KEGG'),
                          ('line', (INK, None, 'data'), 'Nominal'), ('line', (INK, 'dash', 'data'), 'Corrected')], gap=5.0, center=x+w/2)
F.values['e'] = Ev
F.notes.append('e: median rejection probability over terms (display summary); corrected = Bonferroni over terms')
F.save(OUT, D)
