"""Supplementary Figure D18 (printed number from supdraw.SNUM) — same-PAS cross-species comparison (R4), polyApipe v2
primary only (v3 points, columns and legend removed as the plot contract requires).

a  qualified gene pairs per design at the 25 and 50 nt reciprocal matching tolerances (P–M–S design: Xing x CIP, Xing x 2i
   and their union used by the pair mean; P–E–M control: Liu x CIP)
b  human against mouse matched-rank ΔPAS position for every qualified gene pair (25 nt), P–M and P–S contrasts
c  per-pair and pair-mean Spearman ρ (rank axis, 25 nt) with gene-resampling 95% intervals
d  paired correlation differences on the rank and distance axes; matched-null 95% range and the two correction families
   (conditional null, capacity null) for the rank summary
e  tolerance sensitivity (25 vs 50 nt) of the pair-mean summaries, rank and distance axes
f, g  PDLIM7 / Pdlim7: annotated short- and long-CDS transcript models, the two matched PASs and their P and M usage in
   Xing x CIP
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D18'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 132.0, 108.0, 116.0   # 2026-10-12: panel titles moved to axis titles or the caption
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname(FID))
TOL = 25
PAIRS = {'H1_C': ('GSE118258', 'CIP'), 'H1_W': ('GSE118258', 'WOT'), 'H2_C': ('GSE150311', 'CIP')}
CS = {'PM': SC['M'], 'PS': SC['S']}            # a contrast is coloured by its reference state
EXON, CDSC = NULLBAR, INK

cs = D.df('D18A', 'same_pas/v2_primary_correlation_summary.tsv')
assert set(cs.variant) == {'v2_primary'}
nt = D.df('D18A', 'same_pas/v2_primary_null_strata_audit.tsv')


def ngenes(design, pair, tol):
    v = cs[(cs.design == design) & (cs.pair == pair) & (cs.tolerance == tol) & (cs.axis == 'rank')].n_genes.unique()
    assert len(v) == 1, (design, pair, tol, v)
    return int(v[0])


# ================================================================ a  qualified genes
WA = 100.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
bx = Lin(0, 260, x+28, x+w-12)
T25, T50 = SEQ_STOPS[2], SEQ_STOPS[1]
ROWS_A = [('PMS', 'H1_C', y+34.0), ('PMS', 'H1_W', y+50.0), ('PMS', 'pair_mean', y+66.0), ('PEM', 'H2_C', y+92.0)]
F.text(x+0.4, y+16.0, 'P, M and S measured', 7.0, color=MUTED)
F.text(x+0.4, y+76.0, 'P, E and M measured', 7.0, color=MUTED)
Av = {}
for design, pair, yy in ROWS_A:
    if pair == 'pair_mean':
        F.rect(x+0.4, yy-2.8, 5.6, 5.6, fc=STUDY['GSE118258'][1]); F.text(x+7.6, yy-3.4, '× any', 7.0, color=MUTED)
    else:
        F.study_pair(x+0.4, yy, *PAIRS[pair])
    for k, (tol, col) in enumerate(((25, T25), (50, T50))):
        n = ngenes(design, pair, tol)
        F.rect(bx(0), yy-6.6+k*6.8, bx(n)-bx(0), 6.4, fc=col)
        F.text(bx(n)+1.2, yy-6.6+k*6.8+0.1, str(n), 6.8, color=MUTED)
        Av[f'{design}|{pair}|{tol}'] = n
for design in ('PMS', 'PEM'):                       # the union equals the null-strata gene count
    for tol in (25, 50):
        a_ = nt[(nt.design == design) & (nt.tolerance == tol)].iloc[0]
        assert int(a_.n_genes) == ngenes(design, 'pair_mean', tol), (design, tol)
F.xaxis(bx, y+99.0, [0, 100, 200], title='Qualified gene pairs')
F.key_row(None, y+h-7.4, [('block', T25, '25 nt'), ('block', T50, '50 nt')], gap=6.0, center=(bx.p0+bx.p1)/2)
F.values['a'] = Av
F.notes.append('a: qualified = reciprocal PAS match within the tolerance, same P/M/S (or P/E/M) contributors, matched mass >= 0.8 in every state; '
               'the union row is the pair-mean denominator and equals the null-strata gene count (12 PMS / 4 PEM strata, no singletons)')

# ================================================================ b  matched-rank scatter
WB = 200.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
ge = D.df('D18B', 'same_pas/v2_primary_same_axis_gene_effects.tsv.gz')
ge = ge[(ge.tolerance == TOL) & (ge.design == 'PMS') & (ge.mass_min >= 0.8)]
LIM = 0.4
gx = 8.0
xs0 = x+27.0
side = (x+w-1.0-xs0-gx)/2
Bv = {}
for j, pair in enumerate(('H1_C', 'H1_W')):
    z = ge[ge.pair == pair]
    assert len(z) == ngenes('PMS', pair, TOL), (pair, len(z))
    sx = Lin(-LIM, LIM, xs0+j*(side+gx), xs0+j*(side+gx)+side); sy = Lin(-LIM, LIM, y+18.0+side, y+18.0)
    hr, mr = PAIRS[pair]
    cx = (sx.p0+sx.p1)/2
    F.study_pair(cx-9.1, y+6.4, hr, mr)   # gene count (len(z)) moved to the caption
    F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref'); F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    ticks = [-0.3, 0, 0.3]; tl = ['−.3', '0', '.3']
    F.yaxis(sy, sx.p0, ticks, tl, 'Mouse rank effect' if j == 0 else None, labels_off=(j == 1))
    nclip = 0
    for ct in ('PS', 'PM'):
        hx, my = z[f'rank_{ct}_h'].to_numpy(), z[f'rank_{ct}_m'].to_numpy()
        nclip += int(((np.abs(hx) > LIM) | (np.abs(my) > LIM)).sum())
        F.points([sx.clip(v) for v in hx], [sy.clip(v) for v in my], ms='point', color=CS[ct], alpha='secondary')
    yb = F.xaxis(sx, sy.p0+1.0, ticks, tl)
    Bv[pair] = dict(n=int(len(z)), clipped_points=nclip)
F.text(xs0+side+gx/2, yb+1.6, 'Human rank effect', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('marker', ('o', CS['PM']), 'P–M'), ('marker', ('o', CS['PS']), 'P–S')], gap=6.0, center=xs0+side+gx/2)
F.values['b'] = Bv
F.notes.append(f'b: points beyond ±{LIM} drawn at the limit (counts in values); one point per gene pair and contrast')

# ================================================================ c  per-pair rho
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
cx_ = Lin(-0.3, 0.45, x+44.0, x+w-3.0)
q = cs[(cs.tolerance == TOL) & (cs.design == 'PMS') & (cs.axis == 'rank')]
Cv = {}; yy = y+21.0; rh = 8.4
y_top = yy-5
for metric, lab, col in (('PM_rho', 'P–M', CS['PM']), ('PS_rho', 'P–S', CS['PS']), ('MS_rho', 'M–S', MUTED)):
    F.text(x+0.4, yy+rh-3.4, lab, 7.5, font='Bol', color=col)
    for pair in ('H1_C', 'H1_W', 'pair_mean'):
        r = q[(q.pair == pair) & (q.metric == metric)].iloc[0]
        if pair == 'pair_mean':
            F.text(x+21.0, yy-3.4, 'Mean', 7.0, color=MUTED)
        else:
            F.study_pair(x+21.0, yy, *PAIRS[pair])
        F.whisker(cx_.clip(r.low), cx_.clip(r.high), yy, color=INK, lw='thin', cap=1.3)       # CI style of main Figure 4d
        F.marker(cx_(r.estimate), yy, 'o', 'small' if pair != 'pair_mean' else 'mark', fc=col)
        Cv[f'{metric}|{pair}'] = dict(n=int(r.n_genes), rho=float(r.estimate), ci=[float(r.low), float(r.high)])
        yy += rh
    yy += 4.0
F.line([cx_(0)]*2, [y_top, yy-6], lw='ref', color=AXC, ls='ref')
F.xaxis(cx_, yy-4.0, [-0.2, 0, 0.2, 0.4], ['−.2', '0', '.2', '.4'], 'Spearman ρ')
F.key_row(None, y+h-7.4, [('line', (INK, None, 'thin'), 'Bootstrap 95%')], gap=5.0, center=(cx_.p0+cx_.p1)/2)
F.values['c'] = Cv

# ================================================================ d  paired differences and their tests
WD = 262.0
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
cn = D.df('D18C', 'same_pas/v2_primary_conditional_null_tests.tsv')
cap = D.df('D18C', 'same_pas/capacity_sensitivity_tests.tsv')
dx = Lin(-0.4, 0.45, x+76.0, x+w-66.0)
QX1, QX2 = x+w-50.0, x+w-17.5
F.text(QX1, y+3.0, 'Coarse', 7.0, ha='center'); F.text(QX2, y+3.0, 'Capacity', 7.0, ha='center')   # q8 tiles of main Figure 4d
ROWS_D = [('PMS', 'pair_mean', 'PM_minus_PS_rho', 'ρ(P–M) − ρ(P–S)'), ('PMS', 'pair_mean', 'PM_minus_MS_rho', 'ρ(P–M) − ρ(M–S)'),
          ('PEM', 'H2_C', 'PM_minus_EM_rho', 'ρ(P–M) − ρ(E–M)')]
Dv = {}; yy = y+20.0
y_top = yy-6
for design, pair, metric, lab in ROWS_D:
    F.text(x+0.4, yy+1.0, lab, 7.0)
    if pair == 'pair_mean':
        F.text(x+0.4, yy+9.6-3.4+1.6, 'Mean of Xing pairs', 7.0, color=MUTED)
    else:
        F.study_pair(x+0.4, yy+9.6+1.6, *PAIRS[pair])
    nl = cn[(cn.tolerance == TOL) & (cn.design == design) & (cn.metric == metric)].iloc[0]
    cp = cap[(cap.design == design) & (cap.metric == metric)].iloc[0]
    assert abs(nl.estimate - cs[(cs.tolerance == TOL) & (cs.design == design) & (cs.pair == pair) & (cs.axis == 'rank') & (cs.metric == metric)].estimate.iloc[0]) < 1e-9
    for k, (axis, m) in enumerate((('rank', 'o'), ('distance', 's'))):
        r = cs[(cs.tolerance == TOL) & (cs.design == design) & (cs.pair == pair) & (cs.axis == axis) & (cs.metric == metric)].iloc[0]
        ry = yy+4.2+k*7.0
        if axis == 'rank':
            F.line([dx.clip(nl.null_low), dx.clip(nl.null_high)], [ry, ry], lw='range', color=NULLBAR)
        F.whisker(dx.clip(r.low), dx.clip(r.high), ry, color=INK, lw='thin', cap=1.3)
        F.marker(dx(r.estimate), ry, m, 'mark', fc=INK if axis == 'rank' else 'white', ec=INK)
        Dv[f'{metric}|{axis}'] = dict(n=int(r.n_genes), estimate=float(r.estimate), ci=[float(r.low), float(r.high)])
    sig1, sig2 = nl.q_BH_family < 0.05, cp.q_BH_family8 < 0.05
    for qx_, sig, qv_ in ((QX1, sig1, nl.q_BH_family), (QX2, sig2, cp.q_BH_family8)):
        F.rect(qx_-5.0, yy+4.2-4.0, 10.0, 8.0, fc=INK if sig else 'white', ec=INK if sig else NULLW, lw='axis')
        F.text(qx_, yy+11.2-3.4, qfmt(qv_), 6.8, ha='center', color=MUTED)
    Dv[f'{metric}|tests'] = dict(null_95=[float(nl.null_low), float(nl.null_high)], q_conditional=float(nl.q_BH_family),
                                 q_capacity=float(cp.q_BH_family8), family_size=int(nl.family_size))
    yy += 22.0
F.line([dx(0)]*2, [y_top, yy-8], lw='ref', color=AXC, ls='ref')
F.xaxis(dx, yy-6.0, [-0.4, -0.2, 0, 0.2, 0.4], ['−.4', '−.2', '0', '.2', '.4'], title='Correlation difference')
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'Rank'), ('marker', ('s', 'white', INK), 'Distance'),
                          ('line', (INK, None, 'thin'), '95% interval'), ('line', (NULLBAR, None, 'range'), 'Null 95%'),
                          ('block', INK, 'q8 < .05')], gap=5.0, size=7.0, center=x+w/2)
F.values['d'] = Dv
F.notes.append('d: q tiles as in main Figure 4d (coarse = conditional matched null, capacity null; each within its family of 8; filled q < .05), '
               'q printed under each tile; rank-axis summary only')

# ================================================================ e  tolerance sensitivity
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, RR-(X0+WD+GAP), H2)
q = cs[(cs.design == 'PMS') & (cs.pair == 'pair_mean')]
sub_w = (w-26.0-6.0-2.0)/2
Ev = {}
for j, (metric, lab) in enumerate((('PM_rho', 'ρ(P–M)'), ('PM_minus_PS_rho', 'ρ(P–M) − ρ(P–S)'))):
    sx = Lin(20, 55, x+26.0+j*(sub_w+6.0), x+26.0+j*(sub_w+6.0)+sub_w)
    sy = Lin(0, 0.3, y+18.0+sub_w/1.2, y+18.0)
    F.lineplot(sx, sy, metric)
    F.text((sx.p0+sx.p1)/2, y+3.0, lab, 7.5, ha='center')
    F.yaxis(sy, sx.p0, [0, 0.1, 0.2, 0.3], ['0', '.1', '.2', '.3'], 'Pair mean' if j == 0 else None, labels_off=(j == 1))
    for axis, ls, m in (('rank', None, 'o'), ('distance', 'dash', 's')):
        z = q[(q.axis == axis) & (q.metric == metric)].sort_values('tolerance')
        F.line([sx(v) for v in z.tolerance], [sy(v) for v in z.estimate], lw='data', color=INK, ls=ls)
        for t_, v in zip(z.tolerance, z.estimate):
            F.marker(sx(t_), sy(v), m, 'small', fc=INK if axis == 'rank' else 'white', ec=INK)
        Ev[f'{metric}|{axis}'] = {int(t_): float(v) for t_, v in zip(z.tolerance, z.estimate)}
    yb = F.xaxis(sx, sy.p0+1.0, [25, 50], ['25', '50'])
F.text(x+26.0+sub_w+3.0, yb+1.6, 'Matching tolerance (nt)', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('line', (INK, None, 'data'), 'Rank'), ('line', (INK, 'dash', 'data'), 'Distance')], gap=6.0, center=x+w/2)
F.values['e'] = Ev
F.notes.append('e: 25 nt n = 216, 50 nt n = 251 pair-mean genes')

# ================================================================ f, g  PDLIM7 structure and matched PAS usage
tx = D.df('D18H', 'geometry/transcript-models.tsv.gz')
ex_ = D.df('D18H', 'supplement_statistics/same-PAS-example-exon-CDS.tsv.gz')
hp = D.df('D18E', 'axis_extensions/homolog_site_pairs.tsv.gz')
hp = hp[(hp.human_gene_id == 'ENSG00000196923') & (hp.tolerance == TOL) & (hp.design == 'qualified_P_M') & (hp.pair == 'H1_C')]
assert len(hp) == 2
IDS = {'human': ['ENST00000355572', 'ENST00000355841'], 'mouse': ['ENSMUST00000069968', 'ENSMUST00000046246']}
GENE = {'human': ('PDLIM7', 'chr5'), 'mouse': ('Pdlim7', 'chr13')}
WF = (RR-X0-GAP)/2
Gv = {}
for j, (L, sp) in enumerate((('f', 'human'), ('g', 'mouse'))):
    x, y, w, h = F.panel(L, X0+j*(WF+GAP), Y3, WF, H3)
    z = tx[tx.transcript_id.isin(IDS[sp]) & (tx.kind == 'transcript')].set_index('transcript_id').loc[IDS[sp]].reset_index()
    origin = int(z.end0.max())                               # minus strand: 5' end is the largest coordinate
    span = (origin-int(z.start0.min()))/1000
    gx = Lin(0, span+1.4, x+46.0, x+w-3.0)
    sym, chrom = GENE[sp]
    F.text(x+16.0, y+0.4, sym, 8.0, font='BolIta' if sp == 'human' else 'BolIta')
    # chromosome, strand and Ensembl release: caption
    rows = {}
    for i, a in enumerate(z.itertuples()):
        ry = y+26.0+i*18.0; rows[a.transcript_id] = ry
        F.text(x+0.4, ry-3.4, 'Short CDS' if i == 0 else 'Long CDS', 7.0)
        F.text(gx((origin-a.end0)/1000), ry-11.0, a.transcript_id, 6.8, color=MUTED)
        F.line([gx((origin-a.end0)/1000), gx((origin-a.start0)/1000)], [ry, ry], lw='thin', color=AXC)
        seg = ex_[ex_.transcript_id == a.transcript_id]
        assert len(seg[seg.kind == 'exon']) > 0 and len(seg[seg.kind == 'CDS']) > 0
        for kind, hh, col in (('exon', 4.4, EXON), ('CDS', 7.6, CDSC)):
            for v in seg[seg.kind == kind].itertuples():
                lo, hi = gx((origin-v.end0)/1000), gx((origin-v.start0)/1000)
                F.rect(lo, ry-hh/2, max(hi-lo, 0.4), hh, fc=col)
    # matched PAS usage, same genomic x
    cuts = np.abs(hp['human_signed_cut' if sp == 'human' else 'mouse_signed_cut'].to_numpy())
    pos = (origin-cuts)/1000
    sh = 's' if sp == 'human' else 'm'
    uy = Lin(0, 1, y+86.0, y+57.0)
    F.yaxis(uy, gx.p0-4.0, [0, 0.5, 1], ['0', '.5', '1'], None)
    F.vtext(x+3.6, (uy.p0+uy.p1)/2, 'Usage', 7.5)
    for k, (pp, r) in enumerate(sorted(zip(pos, hp.itertuples()), key=lambda t: t[0])):
        F.line([gx(pp)]*2, [y+22.0, uy.p0], lw='ref', color=INK, ls='ref')
        pv, mv = (r.target_h, r.reference_h) if sp == 'human' else (r.target_m, r.reference_m)
        F.line([gx(pp)]*2, [uy(pv), uy(mv)], lw='thin', color=AXC)
        F.marker(gx(pp), uy(pv), 'o', 'mark', fc=SC['P'])
        F.marker(gx(pp), uy(mv), 's', 'mark', fc=SC['M'])
        F.text(gx(pp)-3.6, uy(pv)-3.2, f'{pv:.2f}'.replace('0.', '.'), 6.8, color=SC['P'], ha='right')
        F.text(gx(pp)+3.6, uy(mv)-3.2, f'{mv:.2f}'.replace('0.', '.'), 6.8, color=SC['M'], ha='left')
    sep = abs(pos[0]-pos[1])
    ticks = list(range(0, int(span)+1, 4))
    F.xaxis(gx, uy.p0+1.0, ticks, [str(t) for t in ticks], '5′ → 3′ distance from transcript start (kb)')
    Gv[sp] = dict(symbol=sym, transcripts=IDS[sp], matched_cuts=[int(c) for c in cuts], separation_kb=round(float(sep), 3),
                  P_usage=[float(v) for v in (hp.target_h if sp == 'human' else hp.target_m)],
                  M_usage=[float(v) for v in (hp.reference_h if sp == 'human' else hp.reference_m)])
    F.line([gx(min(pos)), gx(max(pos))], [y+53.0]*2, lw='thin', color=MUTED)
    F.text(gx((pos[0]+pos[1])/2), y+46.2, f'{sep:.2f} kb'.replace('-', '−'), 7.0, ha='center', color=MUTED)
    if L == 'f':
        F.key_row(None, y+h-7.4, [('block', EXON, 'Exon / UTR'), ('block', CDSC, 'CDS'), ('line', (INK, 'ref', 'ref'), 'Matched PAS')],
                  gap=6.0, center=(gx.p0+gx.p1)/2)
    else:
        F.key_row(None, y+h-7.4, [('marker', ('o', SC['P']), 'P usage'), ('marker', ('s', SC['M']), 'M usage')], gap=6.0, center=(gx.p0+gx.p1)/2)
F.values['fg'] = Gv
F.notes.append('f, g: Xing x CIP P–M pair; raw matched-site fractions (unassigned mass not renormalised); keys under f (structure) and g (usage) apply to both')
F.save(OUT, D)
