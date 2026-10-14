"""Supplementary Figure D27 (printed number from supdraw.SNUM) — axis scale and local terminal-family change of the 51
low-net-shift genes (R5), polyApipe v2 primary.

a  native DWUI normalisation span in the two mouse studies for all 51 genes; colour = whether a terminal-family estimate
   exists in both studies (36) or neither (15); ring = gene with a strict recurrent interval (8)
b  terminal-family shift (P − M, nt) in the two studies for the 36 genes with both estimates; strict-interval genes named
c  the eight strict-interval genes: terminal-family shift per study, and the exact decomposition of native ΔDWUI into
   within-family position, family mass and outside-family terms (bold names: |ΔDWUI| < .01 in both studies, the final 4)
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D27'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 153.0, 160.0   # 2026-10-12: panel titles moved to axis titles or the caption
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
RT = ('CIP', 'WOT')
ad = D.df('D27A', 'geometry/axis-decomposition.tsv')
assert ad.gene_id.nunique() == 51
ok = ad.assign(ok=ad.family_status.eq('measurable_native_family')).groupby('gene_id').ok.sum()
strict = ad.groupby('gene_id').strict_region_pairs.first().gt(0)
sym = ad.groupby('gene_id').gene_symbol.first()
at = D.df('D27C', 'candidate_screen/low_shift_attrition.tsv').set_index('stage') if ('tables/candidate_screen/low_shift_attrition.tsv' in D.rb.allowed.get('D27C', set())) else None
FINAL = {'Trip12', 'Rhoa', 'Capzb', 'Rad23b'}
if at is not None:
    FINAL = {sym[g] for g in at.loc['both_abs_DWUI_lt001', 'gene_ids'].split(';')}
OKC = {2: INK, 1: DIST, 0: NULLBAR}
WA = (RR-X0-GAP)/2

# ================================================================ a  axis spans
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
sp = ad.pivot(index='gene_id', columns='route', values='axis_span_nt')/1000
lx = lambda v: np.log10(v)
side = w-34.0-66.0
sx = Lin(0, 2.6, x+30.0, x+30.0+side); sy = Lin(0, 2.6, y+16.0+side, y+16.0)
F.line([sx.p0, sx.p1], [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref')
Av = {}
for k in (0, 1, 2):
    z = sp[ok.reindex(sp.index) == k]
    for g_, r in z.iterrows():
        F.marker(sx.clip(lx(r.CIP)), sy.clip(lx(r.WOT)), 'o', 'small', fc=OKC[k])
    Av[f'family_estimates_{k}_routes'] = int(len(z))
for g_, r in sp[strict.reindex(sp.index)].iterrows():
    F.marker(sx.clip(lx(r.CIP)), sy.clip(lx(r.WOT)), 'o', 'mark', fc='none', ec=INK)
Av['strict'] = int(strict.sum())
T = [0, 1, 2]; TL = ['1', '10', '100']
F.yaxis(sy, sx.p0, T, TL, 'Span, Schiebinger2019–2i (kb)')
F.xaxis(sx, sy.p0+1.0, T, TL, 'Span, Guo2019–CIP (kb)')
kx = sx.p1+6.0
F.key_row(kx, y+20.0, [('marker', ('o', INK), f'family in both ({Av["family_estimates_2_routes"]})')], size=7.0)
F.key_row(kx, y+30.0, [('marker', ('o', NULLBAR), f'in neither ({Av["family_estimates_0_routes"]})')], size=7.0)
F.key_row(kx, y+40.0, [('marker', ('o', 'white', INK), f'strict interval ({Av["strict"]})')], size=7.0)
F.values['a'] = Av
assert Av['family_estimates_1_routes'] == 0

# ================================================================ b  terminal shifts
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WA, H1)
ts = ad.pivot(index='gene_id', columns='route', values='terminal_shift_nt').dropna()
bx = Lin(-600, 150, x+34.0, x+w-3.0); by = Lin(-600, 150, y+16.0+side, y+16.0)
F.line([bx(0)]*2, [by.p0, by.p1], lw='ref', color=AXC, ls='ref'); F.line([bx.p0, bx.p1], [by(0)]*2, lw='ref', color=AXC, ls='ref')
F.line([bx.p0, bx.p1], [by.p0, by.p1], lw='ref', color=AXC, ls='dot')
Bv = {}
for g_, r in ts.iterrows():
    st_ = bool(strict.get(g_, False))
    F.marker(bx(r.CIP), by(r.WOT), 'o', 'mark' if st_ else 'small', fc=INK if st_ else NULLBAR)
    if st_: Bv[sym[g_]] = dict(CIP=float(r.CIP), WOT=float(r.WOT))
OFF = {'Trip12': (6, -2, 'left'), 'Cisd2': (6, -6, 'left'), 'Prpf38b': (-6, 2, 'right'), 'Nxn': (-6, -9, 'right'),
       'Rad23b': (6, 3, 'left'), 'Capzb': (-6, -12, 'right'), 'Rhoa': (6, -8, 'left')}
for name, v in Bv.items():
    dx_, dy_, ha = OFF.get(name, (6, -3, 'left'))
    F.text(bx(v['CIP'])+dx_, by(v['WOT'])+dy_, name, 7.0, font='BolIta' if name in FINAL else 'RegIta', ha=ha)
T = [-600, -300, 0]; TL = ['−600', '−300', '0']
F.yaxis(by, bx.p0, T, TL, 'Shift, Schiebinger2019–2i (nt)')
F.xaxis(bx, by.p0+1.0, T, TL, 'Shift, Guo2019–CIP (nt)')
F.values['b'] = dict(n=int(len(ts)), strict=Bv)

# ================================================================ c  strict genes
x, y, w, h = F.panel('c', X0, Y2, RR-X0, H2)
z = ad[strict.reindex(ad.gene_id).values]
order = sorted(z.gene_symbol.unique(), key=lambda s_: (s_ not in FINAL, s_))
assert len(order) == 8
LBW = 30.0
cx1 = Lin(-600, 150, x+LBW+10.0, x+LBW+150.0)
cx2 = Lin(-30, 45, x+LBW+176.0, x+w-4.0)
COMP = [('within_family_DWUI', SEQ_STOPS[2], 'Family position'), ('family_mass_DWUI', SEQ_STOPS[1], 'Family mass'), ('nonterminal_DWUI', OUTSIDE, 'Outside family')]
Cv = {}; rp = 13.0; c0 = y+20.0
for i, g in enumerate(order):
    yy = c0+i*rp
    F.rect(x+LBW+6.0, yy-rp/2+0.6, w-LBW-6.0, rp-1.2, fc=BAND) if i % 2 == 0 else None
    F.text(x+LBW, yy-3.4, g, 7.0, font='BolIta' if g in FINAL else 'RegIta', ha='right')
    Cv[g] = {}
    for k, route in enumerate(RT):
        r = z[(z.gene_symbol == g) & (z.route == route)].iloc[0]
        ry = yy+(k-0.5)*5.0
        if np.isfinite(r.terminal_shift_nt):
            F.marker(cx1(r.terminal_shift_nt), ry, STUDY[route][2], 'small', fc=STUDY[route][1])
        elif k == 0:
            F.text(cx1.p1-2.0, yy-3.2, 'NA, both studies', 6.8, ha='right', color=MUTED)
        pos = neg = 0.0
        if r.family_status == 'measurable_native_family':
            for col, cc, lab in COMP:
                v = float(r[col])*1000
                left = pos if v >= 0 else neg+v
                F.rect(cx2.clip(left), ry-2.0, cx2.clip(left+abs(v))-cx2.clip(left), 4.0, fc=cc)
                if v >= 0: pos += v
                else: neg += v
            assert abs(r.native_delta_DWUI-sum(r[c] for c, _, _ in COMP)) < 1e-9
        F.marker(cx2(r.native_delta_DWUI*1000), ry, 'D', 'small', fc='white', ec=STUDY[route][1])
        Cv[g][route] = dict(terminal_shift_nt=None if not np.isfinite(r.terminal_shift_nt) else float(r.terminal_shift_nt),
                            native=float(r.native_delta_DWUI), family_status=r.family_status,
                            **({c: float(r[c]) for c, _, _ in COMP} if r.family_status == 'measurable_native_family' else {}))
yend = c0+7*rp+rp/2
F.line([cx1(0)]*2, [c0-6.0, yend], lw='ref', color=AXC, ls='ref'); F.line([cx2(0)]*2, [c0-6.0, yend], lw='ref', color=AXC, ls='ref')
F.xaxis(cx1, yend+1.0, [-600, -300, 0], ['−600', '−300', '0'], title='Terminal-family shift (nt)')
F.xaxis(cx2, yend+1.0, [-30, 0, 30], ['−30', '0', '30'], title='Contribution to native ΔDWUI × 1,000')
F.key_row(None, y+h-17.4, [('marker', ('s', STUDY['CIP'][1]), 'Guo2019–CIP (upper)'), ('marker', ('o', STUDY['WOT'][1]), 'Schiebinger2019–2i (lower)'),
                           ('marker', ('D', 'white', INK), 'native net ΔDWUI')], gap=7.0, center=x+w/2)
F.key_row(None, y+h-7.4, [('block', cc, lab) for _, cc, lab in COMP]+[('block', 'white', 'NA = no measurable terminal family')], gap=7.0, center=x+w/2)
F.values['c'] = Cv
F.notes.append('c: bold = the four genes with |ΔDWUI| < .01 in both studies (Figure 5 screen); unavailable terms are not replaced by zero')
F.save(OUT, D)
