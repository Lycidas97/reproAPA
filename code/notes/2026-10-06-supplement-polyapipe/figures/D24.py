"""Supplementary Figure D24 (printed number from supdraw.SNUM) — external miRNA seed groups in lost 3' UTR windows (R5),
polyApipe v2 primary. Panel a is also the source of the 'PSC miRNA' rug of Figure 5f.

a  mouse external small-RNA expression of each seed (PSC against somatic, log10 RPM + 1) and its expression group
b  share of all lost 3' UTR base pairs that have a matched control window, per species x contrast
c  gene-level lost − control seed-site density for the four seed groups in mouse P–S and P–M
d  the eight primary tests (external somatic − PSC; marker barrier − ESCC), estimate, interval and q over the eight
e  mouse external somatic − PSC contrast under the alternative expression mappings and the position-randomisation nulls
f  mouse P–S: mean after deleting each of the eight genes, against the full-cohort mean
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D24'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 138.0, 140.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
CX = [('human', 'SP'), ('human', 'MP'), ('mouse', 'SP'), ('mouse', 'MP')]
CL = {c: f'{c[0].capitalize()} {c[1][1]}–{c[1][0]}' for c in CX}
GRP = {'PSC': SC['P'], 'somatic': SC['S'], 'other': NULLBAR, 'ESCC': INK, 'barrier': MUTED}

# ================================================================ a  seed expression
WA = 138.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
se = D.df('D24A', 'external_mirna/external_seed_expression.tsv')
m = se[se.species == 'mouse']
ax_ = Lin(0, 6, x+26.0, x+w-3.0); ay = Lin(0, 6, y+16.0+(ax_.p1-ax_.p0)*0.86, y+16.0)
# 2026-10-12: panel titles moved to axis titles or the caption throughout this figure
F.line([ax_.p0, ax_.p1], [ay.p0, ay(6)], lw='ref', color=AXC, ls='ref')
Av = {}
for g in ('other', 'somatic', 'PSC'):
    z = m[m.expression_group == g]
    F.points([ax_.clip(v) for v in np.log10(z.somatic+1)], [ay.clip(v) for v in np.log10(z.PSC+1)], ms='small', color=GRP[g], alpha='secondary')
    Av[g] = int(len(z))
F.yaxis(ay, ax_.p0, [0, 2, 4, 6], None, 'PSC log10(RPM + 1)')
F.xaxis(ax_, ay.p0+1.0, [0, 2, 4, 6], None, 'Somatic log10(RPM + 1)')
F.key_row(None, y+h-7.4, [('marker', ('o', GRP[g]), f'{g} {Av[g]}') for g in ('PSC', 'somatic', 'other')], gap=5.0, center=x+w/2)
F.values['a'] = dict(groups=Av, scope=sorted(m.scope.unique()))

# ================================================================ b  matched coverage
WB = 120.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
at = D.df('D24B', 'external_mirna/attrition.tsv').set_index(['species', 'contrast'])
bx = Lin(0, 0.1, x+44.0, x+w-6.0)
Bv = {}
for i, c in enumerate(CX):
    r = at.loc[c]; yy = y+22.0+i*25.0; fr = r.bp_matched/r.bp_all
    F.text(x+0.4, yy-3.4, CL[c], 7.0)
    F.rect(bx(0), yy-5.0, bx(fr)-bx(0), 10.0, fc=INK)
    F.text(bx(fr)+1.4, yy-3.2, f'{fr*100:.1f}%', 6.8, color=MUTED)
    F.text(bx(0), yy+6.6, f'{int(r.matched_windows)}/{int(r.all_lost_windows)} windows', 6.8, color=MUTED)
    Bv[CL[c]] = dict(bp_fraction=float(fr), matched_windows=int(r.matched_windows), all_windows=int(r.all_lost_windows),
                     matched_genes=int(r.matched_genes), all_genes=int(r.all_lost_genes))
F.xaxis(bx, y+24.0+3*25.0+15.0, [0, 0.05, 0.1], ['0', '.05', '.1'], 'Matched / all lost bp')
F.values['b'] = Bv

# ================================================================ c  gene-level densities
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
gf = D.df('D24C', 'external_mirna/gene_family_densities.tsv')
cx_ = Lin(-30, 15, x+56.0, x+w-7.0)
Cv = {}; yy = y+20.0; rp = 11.0
y_top = yy-5
rng = np.random.default_rng(24)
for c in (('mouse', 'SP'), ('mouse', 'MP')):
    for k, g in enumerate(('PSC', 'somatic', 'ESCC', 'barrier')):
        z = gf[(gf.species == c[0]) & (gf.contrast == c[1]) & (gf.group == g)].difference.to_numpy()
        if k == 0: F.text(x+0.4, yy-3.4, CL[c].replace('Mouse ', 'Mm '), 7.0)
        F.text(x+54.0, yy-3.4, g, 7.0, ha='right', color=MUTED)
        F.points([cx_.clip(v) for v in z], yy+rng.uniform(-2.4, 2.4, len(z)), ms='small', color=GRP[g], alpha='secondary')
        F.line([cx_(np.median(z))]*2, [yy-3.6, yy+3.6], lw='bold', color=INK)
        Cv[f'{CL[c]}|{g}'] = dict(n=int(len(z)), median=float(np.median(z)), clipped=int(((z < -30) | (z > 15)).sum()))
        yy += rp
    yy += 4.0
F.line([cx_(0)]*2, [y_top, yy-8], lw='ref', color=AXC, ls='ref')
F.xaxis(cx_, yy-6.0, [-30, -15, 0, 15], ['−30', '−15', '0', '15'], title='Lost − control sites/kb')
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'gene'), ('line', (INK, None, 'bold'), 'median')], gap=6.0, center=(cx_.p0+cx_.p1)/2)
F.values['c'] = Cv
F.notes.append('c: one point per gene (8 mouse P–S, 24 mouse P–M); median is a display summary; points beyond the axis drawn at the limit')

# ================================================================ d  primary tests
WD = 176.0
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
pt = D.df('D24D', 'external_mirna/primary_tests.tsv')
dx = Lin(-1.0, 1.1, x+80.0, x+w-26.0)
F.text(x+w-16.0, y+3.0, 'n', 7.0, ha='center', color=MUTED); F.text(x+w-4.5, y+3.0, 'q', 7.0, ha='center', color=MUTED)
TL = {'external_somatic_minus_PSC': 'somatic − PSC', 'marker_barrier_minus_ESCC': 'barrier − ESCC'}
Dv = {}; yy = y+20.0; rp = 11.0
y_top = yy-5
for c in CX:
    for k, t_ in enumerate(('external_somatic_minus_PSC', 'marker_barrier_minus_ESCC')):
        r = pt[(pt.species == c[0]) & (pt.contrast == c[1]) & (pt.test == t_)].iloc[0]
        if k == 0: F.text(x+0.4, yy+rp/2-3.4, CL[c].replace('Human ', 'Hs ').replace('Mouse ', 'Mm '), 7.0)
        F.text(x+76.0, yy-3.4, TL[t_], 6.8, ha='right', color=MUTED)
        desc = r.status.startswith('descriptive')
        F.line([dx.clip(r.low), dx.clip(r.high)], [yy, yy], lw='bold', color=NULLBAR)
        F.marker(dx(r.estimate), yy, 'o', 'small', fc='white' if desc else INK, ec=INK)
        F.text(x+w-16.0, yy-3.2, str(int(r.n_genes)), 6.8, ha='center', color=MUTED)
        F.text(x+w-4.5, yy-3.2, '–' if desc else ('1' if r.q_BH8 > 0.999 else qfmt(r.q_BH8)), 6.8, ha='center', color=MUTED if desc else INK)
        Dv[f'{CL[c]}|{t_}'] = dict(estimate=float(r.estimate), interval=[float(r.low), float(r.high)], n=int(r.n_genes), status=r.status, q8=float(r.q_BH8))
        yy += rp
    yy += 2.0
F.line([dx(0)]*2, [y_top, yy-6], lw='ref', color=AXC, ls='ref')
F.xaxis(dx, yy-4.0, [-1, 0, 1], ['−1', '0', '1'], title='Density contrast')
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'tested'), ('marker', ('o', 'white', INK), '< 10 genes')], gap=6.0, center=x+w/2)
F.values['d'] = Dv

# ================================================================ e  sensitivity and position nulls
WE = 160.0
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, WE, H2)
ms_ = D.df('D24E', 'external_mirna/expression_mapping_sensitivity_tests.tsv')
pn = D.df('D24E', 'external_mirna/position_null_tests.tsv')
D.df('D24E', 'external_mirna/complete_test_family.tsv')     # complete 30-test family (receipt)
pn = pn[(pn.species == 'mouse') & (pn.test == 'external_somatic_minus_PSC')]
ex = Lin(-0.2, 0.1, x+60.0, x+w-5.0)
SET = [('exact_totalRPM', 'Exact, total'), ('isomir_mappedRPM', 'isomiR, mapped'), ('isomir_totalRPM', 'isomiR, total')]
Ev = {}; yy = y+20.0; rp = 7.6
y_top = yy-5
for co in ('SP', 'MP'):
    F.text(x+0.4, yy-3.4, f'{co[1]}–{co[0]}', 7.0, font='Bol')
    yy += rp
    for st_, lab in SET:
        r = ms_[(ms_.contrast == co) & (ms_.setting == st_)].iloc[0]
        F.text(x+57.0, yy-3.4, lab, 6.8, ha='right', color=MUTED)
        F.line([ex.clip(r.low), ex.clip(r.high)], [yy, yy], lw='bold', color=NULLBAR)
        F.marker(ex(r.estimate), yy, 'o', 'small', fc='white' if r.status.startswith('descriptive') else INK, ec=INK)
        Ev[f'{co}|{st_}'] = dict(estimate=float(r.estimate), interval=[float(r.low), float(r.high)], n=int(r.n_genes))
        yy += rp
    for mdl, lab in (('coupled', 'Pos. coupled'), ('independent', 'Pos. indep.')):
        r = pn[(pn.contrast == co) & (pn.model == mdl)].iloc[0]
        F.text(x+57.0, yy-3.4, lab, 6.8, ha='right', color=MUTED)
        F.line([ex.clip(r.null_low), ex.clip(r.null_high)], [yy, yy], lw='range', color=NULLBAR)
        F.marker(ex(r.observed_density_contrast), yy, 'D', 'small', fc=INK)
        Ev[f'{co}|position_{mdl}'] = dict(observed=float(r.observed_density_contrast), null=[float(r.null_low), float(r.null_high)], p=float(r.p))
        yy += rp
    yy += 2.0
F.line([ex(0)]*2, [y_top, yy-6], lw='ref', color=AXC, ls='ref')
F.xaxis(ex, yy-4.0, [-0.2, -0.1, 0, 0.1], ['−.2', '−.1', '0', '.1'], title='Somatic − PSC')
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'gene 95%'), ('line', (NULLBAR, None, 'range'), 'null 95%'), ('marker', ('D', INK), 'observed')],
          gap=4.0, center=x+w/2)
F.values['e'] = Ev
F.notes.append('e: all 30 planned controls have q = 1 (complete family in the source table)')

# ================================================================ f  leave-one-out
x, y, w, h = F.panel('f', X0+WD+GAP+WE+GAP, Y2, RR-(X0+WD+GAP+WE+GAP), H2)
inf = D.df('D24F', 'external_mirna/mouse_SP_influence.tsv')
full = float(inf.difference.mean())
assert np.allclose(inf.leave_one_out_mean, (inf.difference.sum()-inf.difference)/(len(inf)-1))
fx = Lin(-0.05, 0.03, x+10.0, x+w-4.0)
yy0 = y+21.0; rp = 12.2
F.line([fx(0)]*2, [yy0-5, yy0+7*rp+5], lw='ref', color=AXC, ls='ref')
F.line([fx(full)]*2, [yy0-5, yy0+7*rp+5], lw='thin', color=INK, ls='dash')
Fv = []
for i, r in enumerate(inf.itertuples()):
    yy = yy0+i*rp
    F.text(x+0.4, yy-3.4, str(i+1), 7.0, color=MUTED)
    F.line([fx(full), fx(r.leave_one_out_mean)], [yy, yy], lw='thin', color=AXC)
    F.marker(fx(r.leave_one_out_mean), yy, 'o', 'small', fc=INK)
    Fv.append(dict(gene_id=r.gene_id, difference=float(r.difference), leave_one_out_mean=float(r.leave_one_out_mean)))
F.xaxis(fx, yy0+7*rp+6.0, [-0.04, 0, 0.02], ['−.04', '0', '.02'], 'Mean after deletion')
F.key_row(None, y+h-7.4, [('line', (INK, 'dash', 'thin'), 'all 8')], gap=5.0, center=(fx.p0+fx.p1)/2)
F.values['f'] = dict(full_mean=full, genes=Fv, crosses_zero=bool((inf.leave_one_out_mean > 0).any() and (inf.leave_one_out_mean < 0).any()))
F.notes.append('f: rows = omitted gene in source-table order; deletion means fall on both sides of zero')
F.save(OUT, D)
