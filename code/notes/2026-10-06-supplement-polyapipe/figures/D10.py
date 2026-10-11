"""Supplementary Figure D10 (printed number from supdraw.SNUM) — well-level inference in Schiebinger2019–2i and external
PAS-seq (R3), polyApipe v2 primary.

a  DWUI vs MPRO joint-family q for the 614 genes tested at well level (mouse P–M members highlighted)
b  genes with joint q < .05 for five metrics: all 614 tested genes and the 143 tested recurrent mouse P–M members
c  external PAS-seq whole-library contrasts (TVD and MPRO) with source-choice ranges and the matching design gate
d  source-variance fraction (two weightings) and the minimum joint q over 100 six-unit downsampling draws
e  exact minimum two-sided p attainable with 4–6 sign units
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D10'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 136.0, 113.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
WOTC = STUDY['WOT'][1]
wt = D.df('D10D', 'supplement_statistics/well-metric-tests-audited.tsv.gz')
rm = D.df('D10B', 'measurement/nulls/recurrence_members.tsv.gz')      # native variant = the primary recurrent members
REC = set(rm[(rm.variant == 'native') & (rm.metric == 'WPP') & (rm.species == 'mouse') & (rm.contrast == 'P_M')
             & (abs(rm.threshold-0.02) < 1e-9)].gene_id)
tested = set(wt.gene_id)
RT = REC & tested
METRICS = ['DWUI', 'RWUI', 'MPRO', 'PDUI', 'PPUI']

# ================================================================ a  DWUI vs MPRO q
WA = 150.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
p = wt.pivot_table(index='gene_id', columns='metric', values='q_BH_joint')
lx = Lin(0, 4, x+24, x+w-2); ly = Lin(0, 4, y+102, y+12)
F.yaxis(ly, lx.p0, [0, 1, 2, 3, 4], title='−log10 q, MPRO')
F.xaxis(lx, ly.p0+1.0, [0, 1, 2, 3, 4], title='−log10 q, DWUI')
th = -math.log10(0.05)
F.line([lx(th)]*2, [ly.p0, ly.p1], lw='ref', color=AXC, ls='ref'); F.line([lx.p0, lx.p1], [ly(th)]*2, lw='ref', color=AXC, ls='ref')
vx = -np.log10(p['DWUI'].clip(lower=1e-4)); vy = -np.log10(p['MPRO'].clip(lower=1e-4))
both = (p['DWUI'] < 0.05) & (p['MPRO'] < 0.05)
isrec = p.index.isin(RT)
for mask, col, al in [(~both & ~isrec, NULLW, 'secondary'), (both & ~isrec, INK, 'secondary'), (isrec, WOTC, 'full')]:
    F.points([lx.clip(v) for v in vx[mask]], [ly.clip(v) for v in vy[mask]], ms='small' if col == WOTC else 'point', color=col, alpha=al)
F.text((lx.p0+lx.p1)/2+6, y+3.0, f'{len(p)} genes, both q < .05: {int(both.sum())}', 7.5, ha='center')
F.key_row(None, y+h-7.4, [('marker', ('o', WOTC), 'Recurrent P–M'), ('marker', ('o', INK), 'Both q < .05')], gap=6.0, center=(lx.p0+lx.p1)/2)
F.values['a'] = dict(genes=int(len(p)), both_q05=int(both.sum()), recurrent_tested=len(RT))

# ================================================================ b  five metrics
WB = 142.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
bx = Lin(0, 1, x+30, x+w-36)
F.text((bx.p0+bx.p1)/2, y+3.0, 'Joint q < .05, fraction', 7.5, ha='center')
F.text(x+w, y+3.0, 'n', 7.0, ha='right', color=MUTED)
Bv = {}; yy = y+20.0
for mtr in METRICS:
    q = wt[wt.metric == mtr]
    a_ = int(q.joint_q05.sum()); r_ = int(q[q.gene_id.isin(RT)].joint_q05.sum())
    F.text(x+0.4, yy+1.4-3.4, mtr, 7.0, font='Bol' if mtr == 'DWUI' else 'Reg')
    for k, (num_, den, col, bh) in enumerate([(a_, len(tested), INK, 4.4), (r_, len(RT), WOTC, 3.0)]):
        yk = yy+k*7.2
        F.rect(bx(0), yk-bh/2, bx(1)-bx(0), bh, fc=BAND)
        F.rect(bx(0), yk-bh/2, bx(num_/den)-bx(0), bh, fc=col)
        F.text(x+w, yk-3.2, f'{num_}/{den}', 6.8, ha='right', color=INK if k == 0 else WOTC)
    Bv[mtr] = dict(all=[a_, len(tested)], recurrent=[r_, len(RT)])
    yy += 19.8
F.xaxis(bx, yy-8.0, [0, 0.5, 1], ['0', '.5', '1'])
F.key_row(None, y+h-7.4, [('block', INK, 'All tested'), ('block', WOTC, 'Recurrent')], gap=6.0, center=x+w/2)
F.values['b'] = Bv
F.notes.append(f'b: recurrent = mouse P–M recurrent members (WPP, |Δ| ≥ .02) that were tested at well level ({len(RT)} of {len(REC)})')

# ================================================================ c  external PAS-seq
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
vs = D.df('D10C', 'external_pas/vector_summary.tsv')
CONF = [('primary', 'Primary'), ('count50', 'Count ≥ 50'), ('count100', 'Count ≥ 100'), ('controls2', '2 controls')]
tx = Lin(0, 0.02, x+46, x+w-30)
mx_ = None
F.text((tx.p0+tx.p1)/2, y+3.0, 'TVD, whole library', 7.5, ha='center')
F.text(x+w, y+3.0, 'Gate', 7.0, ha='right', color=MUTED)
Cv = {}; yy = y+27.0
for i, (cf, cl) in enumerate(CONF):
    q = vs[(vs.configuration == cf) & (vs.statistic == 'tvd')].iloc[0]
    if i % 2 == 0: F.rect(x, yy-11.0, w, 22.0, fc=BAND)
    F.text(x+0.4, yy-5.6, cl, 7.0)
    F.text(x+0.4, yy+1.0, f'{int(q.matched_n)}/{int(q.measurable_n)}', 6.8, color=MUTED)
    F.line([tx(q.sensitivity_min), tx(q.sensitivity_max)], [yy, yy], lw='bold', color=NULLBAR)
    F.marker(tx(q.estimate), yy, 'o', 'mark', fc=INK)
    F.text(x+w-0.4, yy-3.4, 'pass' if q.design_gate_pass else 'fail', 7.0, ha='right', color=INK if q.design_gate_pass else MUTED)
    Cv[cf] = dict(tvd=float(q.estimate), range=[float(q.sensitivity_min), float(q.sensitivity_max)], gate=bool(q.design_gate_pass),
                  matched=int(q.matched_n), measurable=int(q.measurable_n), max_abs_smd=float(q.max_abs_smd),
                  label_assignments=int(q.label_assignments))
    yy += 22.0
F.xaxis(tx, yy-10.0, [0, 0.01, 0.02], ['0', '.01', '.02'])
F.key_row(None, y+h-17.4, [('line', (NULLBAR, None, 'bold'), 'Source range')], gap=6.0, center=(tx.p0+tx.p1)/2)
F.text(x+w/2, y+h-7.4, 'Matched / measurable genes', 7.0, ha='center', color=MUTED)
F.values['c'] = Cv
F.notes.append('c: six whole-library label assignments; descriptive diagnostic, not an external significance test')

# ================================================================ d  source variance and six-unit draws
WD = 296.0
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
ec = D.df('D10D', 'panel_statistics/well-source-variance-source_variance_fraction-ecdf.tsv.gz')
ds = D.df('D10D', 'supplement_statistics/well_six_unit_downsample_summary.tsv')
cy = Lin(0, 1, y+78, y+16)
d1 = Lin(0, 1, x+26, x+142)
d2 = Lin(0.04, 0.09, x+167, x+w-7)
F.lineplot(d1, cy, 'source variance'); F.lineplot(d2, cy, 'min q')
F.text((d1.p0+d1.p1)/2, y+3.0, 'Source-variance fraction', 7.5, ha='center')
F.text((d2.p0+d2.p1)/2, y+3.0, 'Six-unit draws (100)', 7.5, ha='center')
F.yaxis(cy, d1.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Cumulative fraction', title_y=(cy.p0+cy.p1)/2+4)
F.yaxis(cy, d2.p0, [0, 0.5, 1], ['0', '.5', '1'])
Dv = {}
for mt, col, ls in [('count_WPP', INK, None), ('equal_MC_WPP', MUTED, 'dash')]:
    q = ec[ec.metric == mt].sort_values('order')
    xs_ = np.r_[q.value.iloc[0], q.value.to_numpy()]; ys_ = np.r_[0, q.cumulative.to_numpy()]
    F.line([d1(v) for v in np.repeat(xs_, 2)[1:]], [cy(v) for v in np.repeat(ys_, 2)[:-1]], lw='data', color=col, ls=ls)
    Dv[mt] = dict(n=int(len(q)), median=float(q.value.iloc[int(np.searchsorted(q.cumulative.to_numpy(), 0.5))]))
for mtr in METRICS:
    v = np.sort(ds[ds.metric == mtr].min_joint_q.to_numpy()); n = len(v)
    xs_ = np.repeat(v, 2); ys_ = np.repeat(np.arange(n+1)/n, 2)[1:-1]
    F.line([d2.clip(t) for t in xs_], [cy(t) for t in ys_], lw='data' if mtr == 'DWUI' else 'thin',
           color=INK if mtr == 'DWUI' else NULLW, alpha='full' if mtr == 'DWUI' else 'secondary')
    Dv[f'draws|{mtr}'] = dict(draws=int(n), min=float(v.min()), median=float(np.median(v)), max=float(v.max()),
                              downsample_significant=int(ds[ds.metric == mtr].downsample_significant.max()))
F.xaxis(d1, cy.p0+1.0, [0, 0.5, 1], ['0', '.5', '1'])
F.line([d2(0.05)]*2, [cy.p0, cy.p1], lw='ref', color=AXC, ls='ref')
F.xaxis(d2, cy.p0+1.0, [0.04, 0.05, 0.06, 0.07, 0.08, 0.09], ['.04', '.05', '.06', '.07', '.08', '.09'], 'Minimum joint q')
F.key_row(None, y+h-7.4, [('line', (INK, None, 'data'), 'Count weights'), ('line', (MUTED, 'dash', 'data'), 'Equal metacells')],
          gap=6.0, center=(d1.p0+d1.p1)/2)
F.key_row(None, y+h-7.4, [('line', (INK, None, 'data'), 'DWUI'), ('line', (NULLW, None, 'thin'), 'Other metrics')], gap=6.0,
          center=(d2.p0+d2.p1)/2)
F.values['d'] = Dv

# ================================================================ e  attainable p
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, RR-(X0+WD+GAP), H2)
pr = D.df('D10E', 'axis_extensions/pool_resolution.tsv')
ex = Lin(3.5, 6.5, x+28, x+w-1)
ey = Lin(0, 0.14, y+78, y+16)
F.text((ex.p0+ex.p1)/2, y+3.0, 'Exact sign test, minimum', 7.5, ha='center')
F.yaxis(ey, ex.p0, [0, 0.05, 0.1], ['0', '.05', '.1'], 'Two-sided p')
F.line([ex.p0, ex.p1], [ey(0.05)]*2, lw='ref', color=AXC, ls='ref')
Ev = {}
for _, q in pr.iterrows():
    F.marker(ex(q.n_units), ey(q.exact_two_sided_sign_min_p), 'o', 'mark', fc=INK)
    F.text(ex(q.n_units), ey(q.exact_two_sided_sign_min_p)-10.0, f'{q.exact_two_sided_sign_min_p:g}'.replace('0.', '.'), 7.0, ha='center')
    Ev[int(q.n_units)] = dict(min_p=float(q.exact_two_sided_sign_min_p), genes=int(q.genes))
F.xaxis(ex, ey.p0+1.0, [4, 5, 6], title='Independent sign units')
F.text((ex.p0+ex.p1)/2, y+h-7.4, 'Genes: ' + ' / '.join(f'{int(g):,}' for g in pr.genes), 7.0, ha='center', color=MUTED)
F.values['e'] = Ev
F.save(OUT, D)
