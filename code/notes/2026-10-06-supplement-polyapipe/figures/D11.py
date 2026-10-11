"""Supplementary Figure D11 (printed number from supdraw.SNUM) — dual-reference candidates and PDLIM7 (R3), polyApipe v2.

a  recurrent candidates -> matched-measurable -> dual-qualified (both P–S and P–M), human Xing2020 and mouse CIP/WOT
b  DWUI-only / shared / MPRO-only recurrent membership (|ΔDWUI| >= .02, MPRO >= .1)
c  PDLIM7 matched ΔDWUI for P–S, P–M and M–S in the three routes
d  matched P–S against P–M ΔDWUI for every measurable gene; dual-qualified filled
e  PDLIM7 cumulative PAS usage (same contributors for S, P and M) in the three routes
"""
import sys, math, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D11'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 108.0, 117.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
RTS = [('GSE118258', 'human'), ('CIP', 'mouse'), ('WOT', 'mouse')]
SLS = {'S': 'dash', 'P': None, 'M': 'dashdot'}      # line styles of main Figures 3d, 4e, 5e

# ================================================================ a  admission ladder
WA = 150.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
lad = D.df('D11A', 'candidate_screen/dual_reference_ladder.tsv').set_index('species')
ax_ = Lin(0, 450, x+42, x+w-9)
F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Genes', 7.5, ha='center')
STEPS = [('recurrent_candidates', 'Recurrent'), ('matched_measurable', 'Measurable'), ('dual_reference', 'Aligned')]
Av = {}; yy = y+20.0
for sp, lab, cols in [('human', 'Human', [STUDY['GSE118258'][1]]), ('mouse', 'Mouse', [STUDY['CIP'][1], STUDY['WOT'][1]])]:
    q = lad.loc[sp]
    F.text(x+0.4, yy-3.4, lab, 7.5, font='Bol')
    for k, (col, sl) in enumerate(STEPS):
        yk = yy+8.6+k*8.6
        F.text(ax_.p0-3.0, yk-3.4, sl, 7.0, ha='right')
        v = int(q[col])
        if col == 'dual_reference':
            F.rect(ax_(0), yk-2.6, ax_(q.distal)-ax_(0), 5.2, fc=DIST)
            F.rect(ax_(q.distal), yk-2.6, ax_(v)-ax_(q.distal), 5.2, fc=PROX)
        else:
            F.rect(ax_(0), yk-2.6, ax_(v)-ax_(0), 5.2, fc=NULLBAR)
        F.text(ax_(v)+2.0, yk-3.4, str(v), 7.0)
    Av[sp] = {c: int(q[c]) for c in ['recurrent_candidates', 'matched_measurable', 'dual_reference', 'distal', 'proximal']}
    yy += 35.0
F.xaxis(ax_, yy-6.6, [0, 200, 400])
F.key_row(None, y+h-7.4, [('block', DIST, 'Distal'), ('block', PROX, 'Proximal')], gap=6.0, center=x+w/2)
F.values['a'] = Av

# ================================================================ b  metric membership
WB = 150.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
mm = D.df('D11B', 'matched_pms/metric_membership.tsv.gz')
mm = mm[(abs(mm.DWUI_threshold-0.02) < 1e-9) & (abs(mm.MPRO_threshold-0.1) < 1e-9)]
CO = [('DWUI_only', 'DWUI only', ORD3[0]), ('shared', 'Both', ORD3[1]), ('MPRO_only', 'MPRO only', ORD3[2])]
bx = Lin(0, 1100, x+44, x+w-1)
F.text((bx.p0+bx.p1)/2, y+3.0, 'Recurrent genes', 7.5, ha='center')
Bv = {}
for i, (sp, ct, gl) in enumerate([('human', 'P_S', 'Human P–S'), ('human', 'P_M', 'Human P–M'), ('mouse', 'P_S', 'Mouse P–S'),
                                  ('mouse', 'P_M', 'Mouse P–M')]):
    yy = y+22.0+i*19.4
    F.text(x+0.4, yy-3.4, gl, 7.0)
    acc = 0; Bv[gl] = {}
    for co, cl, cc in CO:
        n = int(((mm.species == sp) & (mm.contrast == ct) & (mm.cohort == co)).sum())
        F.rect(bx(acc), yy-5.8, bx(acc+n)-bx(acc), 11.6, fc=cc, ec='white', lw='axis')
        if bx(acc+n)-bx(acc) > 16: F.text((bx(acc)+bx(acc+n))/2, yy-3.2, str(n), 6.8, ha='center', color=on(cc))
        acc += n; Bv[gl][co] = n
F.xaxis(bx, y+22.0+3*19.4+6.0, [0, 500, 1000])
F.key_row(None, y+h-7.4, [('block', cc, cl) for _, cl, cc in CO], gap=5.0, center=x+w/2)
F.values['b'] = Bv

# ================================================================ c  PDLIM7 matched effects
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
dm = D.df('D11D', 'supplement_statistics/dual-reference-all-measurable.tsv.gz')
pd7 = dm[dm.gene_symbol.str.upper() == 'PDLIM7'].set_index('route')
cx = Lin(-0.7, 0.1, x+30, x+w-1)
F.text(x+15.0+textwidth('c', 11, 'Bol'), y+2.4, 'PDLIM7/Pdlim7', 8.0, font='BolIta')   # 2026-10-08: both species
F.text(x+w, y+3.0, 'ΔDWUI', 7.5, ha='right')   # 2026-10-08: shorter header (matched in caption) to make room for the two-species title
Cv = {}; r0 = y+20.0; rh = 7.4
F.line([cx(0)]*2, [r0-5, r0+rh*8+2*3.0+5], lw='ref', color=AXC, ls='ref')
for i, (r, sp) in enumerate(RTS):
    q = pd7.loc[r]; Cv[r] = {}
    for k, (c_, cl) in enumerate([('DWUI_PS', 'P–S'), ('DWUI_PM', 'P–M'), ('DWUI_MS', 'M–S')]):
        yy = r0+rh*(3*i+k)+3.0*i
        if k == 0: F.marker(x+4.0, yy+rh, STUDY[r][2], 'mark', fc=STUDY[r][1])
        F.text(x+10.0, yy-3.4, cl, 7.0)
        F.line([cx(0), cx(q[c_])], [yy, yy], lw='thin', color=STUDY[r][1])
        F.marker(cx(q[c_]), yy, 'o', 'mark', fc=STUDY[r][1])
        Cv[r][c_] = float(q[c_])
F.xaxis(cx, r0+rh*8+6.0+4.5, [-0.6, -0.4, -0.2, 0], ['−.6', '−.4', '−.2', '0'])
F.key_row(None, y+h-7.4, [('marker', (STUDY[r][2], STUDY[r][1]), STUDY[r][0].split('–')[-1] if r != 'GSE118258' else 'Xing')
                          for r, _ in RTS], gap=5.0, center=x+w/2)
F.values['c'] = Cv

# ================================================================ d  P–S vs P–M, all measurable genes
WD = 186.0
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
pw = (w-30-10-4)/2
Dv = {}
for j, (sp, routes, title) in enumerate([('human', ['GSE118258'], 'Human'), ('mouse', ['CIP', 'WOT'], 'Mouse')]):
    sx = Lin(-0.6, 0.6, x+30+(pw+10)*j, x+30+(pw+10)*j+pw)
    sy = Lin(-0.6, 0.6, y+16+pw, y+16)
    q = dm[dm.species == sp]
    F.text((sx.p0+sx.p1)/2, y+3.0, f'{title} · {q.gene_id.nunique()} genes', 7.5, ha='center')
    F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref'); F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    F.line([sx.p0, sx.p0], [sy.p0, sy.p1], lw='axis', color=AXC)
    if j == 0: F.yaxis(sy, sx.p0, [-0.6, -0.3, 0, 0.3, 0.6], ['−.6', '−.3', '0', '.3', '.6'], 'P–M ΔDWUI', line=False)
    else: F.yaxis(sy, sx.p0, [-0.6, -0.3, 0, 0.3, 0.6], line=False, labels_off=True)
    for dual in (False, True):                                   # encoding of main Figure 3c: open = P–M cohort, ink = aligned
        for r in routes:
            z = q[q.route == r]; zz = z[z.dual_qualified_gene == dual]
            F.points([sx.clip(v) for v in zz.DWUI_PS], [sy.clip(v) for v in zz.DWUI_PM], ms='small', m=STUDY[r][2],
                     color=INK if dual else MUTED, alpha='full' if dual else 'secondary', fc=INK if dual else 'white', mew=0.0 if dual else 0.5)
            if dual: Dv[r] = dict(genes=int(len(z)), dual=int(z.dual_qualified_gene.sum()))
    if sp == 'mouse':                                           # study = marker shape (empty upper-left corner)
        F.key_row(sx.p0+2.0, sy.p1+1.0, [('marker', ('s', INK), 'CIP'), ('marker', ('o', INK), '2i')], gap=4.0, size=7.0)
    F.xaxis(sx, sy.p0+1.0, [-0.6, -0.3, 0, 0.3, 0.6], ['−.6', '−.3', '0', '.3', '.6'])
F.text(x+30+pw+5, y+16+pw+14.0, 'P–S ΔDWUI', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'P–S/P–M aligned'), ('marker', ('o', 'white', MUTED), 'P–M cohort'),
                          ], gap=6.0, center=x+30+pw+5)
F.values['d'] = Dv
F.notes.append('d: values beyond ±.6 drawn at the axis limit')

# ================================================================ e  PDLIM7 cumulative usage
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, RR-(X0+WD+GAP), H2)
cu = D.df('D11E', 'panel_statistics/example-cumulative-PAS-curves.tsv.gz')
cu = cu[(cu.gene_symbol.str.upper() == 'PDLIM7') & (cu.source == 'same_contributors_PMS')]
pw = (w-26-2*8-4)/3
ey = Lin(0, 1, y+20+min(pw, 62.0), y+20)
Ev = {}
for j, (r, sp) in enumerate(RTS):
    q = cu[cu.route == r]
    xmax = q.oriented_position_nt.max()/1000
    ex = Lin(0, 10, x+26+(pw+8)*j, x+26+(pw+8)*j+pw)
    F.lineplot(ex, ey, r)
    F.rect(ex.p0, y+3.6, 5.6, 5.6, fc=STUDY[r][1]) if False else None
    F.marker(ex.p0+2.8, y+6.0, STUDY[r][2], 'mark', fc=STUDY[r][1])
    F.text(ex.p0+7.6, y+3.0, STUDY[r][0].split('–')[-1] if r != 'GSE118258' else 'Xing2020', 7.0)
    F.line([ex.p0, ex.p0], [ey.p0, ey.p1], lw='axis', color=AXC)
    if j == 0: F.yaxis(ey, ex.p0, [0, 0.5, 1], ['0', '50', '100'], 'Distal (%)', line=False, title_y=(ey.p0+ey.p1)/2+4)
    else: F.yaxis(ey, ex.p0, [0, 0.5, 1], line=False, labels_off=True)
    Ev[r] = {}
    for st in ['S', 'P', 'M']:
        z = q[q.state == st].sort_values('oriented_position_nt')
        xs_ = [0.0]+list(z.oriented_position_nt/1000)+[10.0]; ys_ = [0.0]+list(z.cumulative)+[1.0]
        px, py = [], []
        for i_ in range(len(xs_)-1):
            px += [xs_[i_], xs_[i_+1]]; py += [ys_[i_], ys_[i_]]
        F.line([ex(v) for v in px], [ey(1-v) for v in py], lw='data', color=SC[st], ls=SLS[st])   # distal fraction = 1 − cumulative (main Figure 3d)
        Ev[r][st] = [[int(p_), float(c_)] for p_, c_ in zip(z.oriented_position_nt, z.cumulative)]
    F.xaxis(ex, ey.p0+1.0, [0, 5, 10])
F.text(x+26+(3*pw+16)/2, ey.p0+12.6, 'Distance from proximal PAS (kb)', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('line', (SC[s], SLS[s], 'data'), s) for s in ['S', 'P', 'M']], gap=8.0, center=x+26+(3*pw+16)/2)
F.values['e'] = Ev
F.save(OUT, D)
