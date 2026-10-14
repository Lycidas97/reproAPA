"""Supplementary Figure D21 (printed number from supdraw.SNUM) — recurrent interval changes and weighted predicted-site
content (R5), polyApipe v2 primary.

a  gene denominators per species x contrast cohort: union pool, genes measurable in both routes, genes with a recurrent
   interval pair (and the number of interval pairs)
b  the low-net-shift screen of Figure 5: 51 -> 8 -> 4 genes (all counts are genes; interval pairs in brackets)
c  matched-interval inclusion change (ΔI) in the two studies of each species, all qualified interval pairs (P–S and P–M)
d  observed sign agreement of the interval pairs with gene-bootstrap 95% intervals (descriptive)
e  weighted predicted-site change per gene (20 UMI), by cohort and site category, with loss / gain gene counts and genes without a complete site
f  sensitivity of the weighted estimate to the 50 UMI library threshold (50 − 20 UMI), cumulative over genes
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D21'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 99.0, 176.0, 142.0   # 2026-10-12: panel titles moved to axis titles or the caption
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname(FID))
CX = [('human', 'SP'), ('human', 'MP'), ('mouse', 'SP'), ('mouse', 'MP')]
CL = {c: f'{c[0].capitalize()} {c[1][1]}–{c[1][0]}' for c in CX}          # SP -> P–S
ROUTES = {'human': ('GSE118258', 'GSE150311'), 'mouse': ('CIP', 'WOT')}
st = D.df('D21A', 'site_statistics/region_consistency_statistics.tsv')
st = st[(st.category == 'all') & (st.unit == 'region_pair')].set_index(['species', 'contrast'])

# ================================================================ a  denominators
WA = 305.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
claims = D.df('D21A', 'union_pool/figure5_union_claims.tsv')
me = D.df('D21A', 'intervals/model_effects.tsv')
D.df('D21A', 'union_pool/figure5_MPRO_admission.tsv')      # admission audit (receipt only)
ax_ = Lin(0, 2200, x+64.0, x+w-50.0)
F.text(x+w-22.0, y+3.0, 'Interval pairs', 7.0, ha='center', color=MUTED)
Av = {}
for i, c in enumerate(CX):
    sp, co = c; yy = y+21.0+i*14.0
    hr, mr = ROUTES[sp]
    F.rect(x+0.4, yy-2.8, 5.6, 5.6, fc=STUDY[hr][1]); F.rect(x+7.4, yy-2.8, 5.6, 5.6, fc=STUDY[mr][1])
    F.text(x+16.0, yy-3.4, CL[c], 7.0)
    genes = set(claims.loc[(claims.species == sp) & (claims.contrast == co), 'gene_id'])
    q = me[(me.species == sp) & (me.contrast == co) & me.eligible & (me.annotated_pas >= 2) & (me.n_ref50 >= 1) & (me.n_target50 >= 1)]
    nr = q.groupby('gene_id').route.nunique(); meas = len(set(nr[nr == 2].index) & genes)
    s_ = st.loc[c]; ig, ip = int(s_.genes), int(s_.n)
    F.line([ax_(ig), ax_(len(genes))], [yy, yy], lw='thin', color=AXC)
    F.marker(ax_(len(genes)), yy, 'o', 'mark', fc='white', ec=INK)
    F.marker(ax_(meas), yy, 'o', 'mark', fc=INK)
    F.marker(ax_(ig), yy, 'D', 'mark', fc=DIST)
    F.text(x+w-22.0, yy-3.4, f'{ip}', 7.0, ha='center')
    Av[CL[c]] = dict(union=len(genes), measurable_both_routes=meas, interval_genes=ig, interval_pairs=ip)
F.xaxis(ax_, y+21.0+3*14.0+7.0, [0, 1000, 2000], ['0', '1,000', '2,000'], title='Genes')
F.key_row(None, y+h-7.4, [('marker', ('o', 'white', INK), 'Union pool'), ('marker', ('o', INK), 'Measurable in both studies'),
                          ('marker', ('D', DIST), 'With interval pair')], gap=6.0, center=x+w/2)
F.values['a'] = Av
F.notes.append('a: measurable = eligible gene model with >= 2 annotated PASs and >= 1 reference and target library at 50 UMI in both studies '
               '(stored eligibility columns, no new selection); interval pairs are not a gene attrition stage')

# ================================================================ b  low-shift screen
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
at = D.df('D21B', 'candidate_screen/low_shift_attrition.tsv').set_index('stage')
ls = D.df('D21B', 'geometry/low-shift-cohort.tsv').set_index('gene_id')
STG = [('low_abs_DWUI_lt002_MPRO_ge01', 'Low net shift', '|ΔDWUI| < .02, MPRO ≥ .1'), ('strict_recurrent_interval', 'Strict recurrent interval', 'in both mouse studies'),
       ('both_abs_DWUI_lt001', 'Both |ΔDWUI| < .01', 'CIP and 2i')]
Bv = {}
bxw = 30.0; BX = x+0.4
for i, (k, lab, sub) in enumerate(STG):
    yy = y+18.0+i*27.0
    n, npair = int(at.loc[k, 'genes']), int(at.loc[k, 'interval_pairs'])
    F.rect(BX, yy, bxw, 15.0, fc=BAND, ec=NULLBAR, lw='axis')                 # box and arrow of main Figure 5c
    F.text(BX+bxw/2, yy+3.8, str(n), 8.0, ha='center', font='Bol')
    F.text(BX+bxw+5.0, yy+0.2, lab, 7.0)
    F.text(BX+bxw+5.0, yy+8.0, sub, 6.8, color=MUTED)
    if i < 2:
        F.arrow_down(BX+bxw/2, yy+15.6, yy+27.0)
    Bv[lab] = dict(genes=n, interval_pairs=npair, gene_ids=at.loc[k, 'gene_ids'].split(';'))
fin = [ls.loc[g_, 'gene_symbol'] for g_ in at.loc['both_abs_DWUI_lt001', 'gene_ids'].split(';')]
F.text(BX+bxw+5.0, y+18.0+2*27.0+16.6, ', '.join(fin), 7.0, font='RegIta')
F.values['b'] = Bv
assert [Bv[s[1]]['genes'] for s in STG] == [51, 8, 4]

# ================================================================ c  ΔI scatter
WC = 300.0
x, y, w, h = F.panel('c', X0, Y2, WC, H2)
rg = D.df('D21C', 'site_statistics/all_changed_region_pairs.tsv.gz')
rg = rg[rg.contrast.isin(['SP', 'MP'])]
side = (w-30.0-12.0-1.0)/2
Cv = {}
for j, sp in enumerate(('human', 'mouse')):
    sx = Lin(-100, 100, x+30.0+j*(side+12.0), x+30.0+j*(side+12.0)+side); sy = Lin(-100, 100, y+18.0+side, y+18.0)
    a_, b_ = ROUTES[sp]
    cx_ = (sx.p0+sx.p1)/2
    F.text(cx_, y+3.0, sp.capitalize(), 7.5, ha='center')
    F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref'); F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    for co, m in (('SP', 'o'), ('MP', '^')):
        z = rg[(rg.species == sp) & (rg.contrast == co)]
        assert len(z) == int(st.loc[(sp, co), 'n']), (sp, co, len(z))
        ea, eb = z.effect_mean_a.to_numpy()*100, z.effect_mean_b.to_numpy()*100
        col = np.where(ea*eb <= 0, 'grey', np.where(ea < 0, 'loss', 'gain'))
        for key, cc in (('grey', NULLBAR), ('loss', PROX), ('gain', DIST)):
            k_ = col == key
            F.points([sx.clip(v) for v in ea[k_]], [sy.clip(v) for v in eb[k_]], ms='small', color=cc, alpha='secondary', m=m)
        Cv[f'{sp}|{co}'] = dict(pairs=int(len(z)), both_loss=int((col == 'loss').sum()), both_gain=int((col == 'gain').sum()),
                                other=int((col == 'grey').sum()))
    F.yaxis(sy, sx.p0, [-50, 0, 50], ['−50', '0', '50'], 'ΔI in second study (pp)' if j == 0 else None,
            labels_off=(j == 1))
    F.rect(sx.p0+3.0, sy.p1+2.0, 5.6, 5.6, fc=STUDY[b_][1]); F.text(sx.p0+10.6, sy.p1+1.4, STUDY[b_][0], 7.0, color=MUTED)
    yb = F.xaxis(sx, sy.p0+1.0, [-50, 0, 50], ['−50', '0', '50'])
    F.rect(sx.p1-5.6-1.0-textwidth(STUDY[a_][0], 7.0)-2.0, sy.p0-9.0, 5.6, 5.6, fc=STUDY[a_][1])
    F.text(sx.p1-1.0, sy.p0-9.6, STUDY[a_][0], 7.0, ha='right', color=MUTED)
F.text(x+30.0+side+6.0, yb+1.6, 'ΔI in first study (pp)', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('marker', ('o', PROX), 'Both loss'), ('marker', ('o', DIST), 'Both gain'), ('marker', ('o', NULLBAR), 'Other signs'),
                          ('marker', ('o', INK), 'P–S'), ('marker', ('^', INK), 'P–M')], gap=6.0, center=x+30.0+side+6.0)
F.values['c'] = Cv
F.notes.append('c: one marker per interval pair (exact ties overplot); first study = Xing2020 / Guo2019–CIP, second = Liu2020 / Schiebinger2019–2i')

# ================================================================ d  agreement (form of main Figure 5a, left)
x, y, w, h = F.panel('d', X0+WC+GAP, Y2, RR-(X0+WC+GAP), H2)
sa = D.df('D21D', 'site_statistics/region_consistency_statistics.tsv')
sa = sa[sa.unit == 'region_pair'].set_index(['species', 'contrast', 'category'])
dx = Lin(0.5, 1.0, x+50.0, x+w-22.0)
F.text(x+w-9.0, y+3.0, 'n', 7.0, ha='center', color=MUTED)
Dv = {}; pitch = 28.5; y0d = y+32.5
for i, c in enumerate(CX):
    yy = y0d+pitch*i
    if i % 2 == 0: F.rect(x, yy-pitch/2+1.0, w, pitch-2.0, fc=BAND)
    F.text(x+0.4, yy-3.4, CL[c], 7.0)
    for cat, dy, mk in (('miRNA', -8.0, '^'), ('all', 0.0, 'o'), ('RBP', 8.0, 'v')):
        s_ = sa.loc[(c[0], c[1], cat)]; lo, hi, v = float(s_.cobs_ci_low), float(s_.cobs_ci_high), float(s_.cobs)
        if cat == 'all':
            F.whisker(dx(lo), dx(hi), yy, color=INK, lw='thin', cap=1.3)
            F.marker(dx(v), yy, mk, 'mark', fc='white' if int(s_.genes) < 20 else INK, ec=INK)
        else:
            if hi > lo: F.line([dx(lo), dx(hi)], [yy+dy]*2, lw='axis', color=NULLBAR)
            F.marker(dx(v), yy+dy, mk, 'small', fc=AXC)
        F.text(x+w-9.0, yy+dy-3.2, str(int(s_.n)), 6.8, ha='center', color=MUTED if cat != 'all' else INK)
        Dv[f'{CL[c]}|{cat}'] = dict(n=int(s_.n), genes=int(s_.genes), agreement=float(v), ci=[float(lo), float(hi)])
yend = y0d+pitch*3+pitch/2
F.line([dx(1.0)]*2, [y0d-pitch/2+1.0, yend-1.0], lw='ref', color=AXC, ls='ref')
F.xaxis(dx, yend+1.0, [0.5, 0.75, 1.0], ['50', '75', '100'], 'Direction agreement (%)')
F.key_row(None, y+h-17.4, [('marker', ('o', INK), 'Observed'), ('line', (INK, None, 'thin'), 'Gene-bootstrap 95%')], gap=5.0, size=7.0,
          center=x+w/2)
F.key_row(None, y+h-7.4, [('marker', ('^', AXC), 'miRNA'), ('marker', ('v', AXC), 'RBP')],   # every cohort has >= 20 genes (23-246)
          gap=5.0, size=7.0, center=x+w/2)
Dv = {k: Dv[k] for k in [f'{CL[c]}|{cat}' for c in CX for cat in ('all', 'miRNA', 'RBP')]}
F.values['d'] = Dv
F.notes.append('d: fraction of interval pairs with the same ΔI sign; all sites and per predicted-site category; descriptive agreement within the selected cohorts, no expectation or observed-minus-expected claim')

# ================================================================ e  weighted site change
WE = 300.0
x, y, w, h = F.panel('e', X0, Y3, WE, H3)
gr = D.df('D21E', 'weighted_retention/gene_retention.tsv')
LBL = 66.0
fl = lambda v: math.copysign(math.log10(1+abs(v)), v)                 # symmetric log axis of main Figure 5b
ex = Lin(fl(-250), fl(250), x+LBL, x+LBL+132.0)
cxw = Lin(-250, 250, x+LBL+132.0+14.0, x+w-32.0)
F.text(x+w-14.0, y+3.0, 'No sites', 7.0, ha='center', color=MUTED)
Ev = {}; yy = y+20.0; rh = 10.4
y_top = yy-6
for c in CX:
    for cat in ('miRNA', 'RBP'):
        z = gr[(gr.species == c[0]) & (gr.contrast == c[1]) & (gr.category == cat)].weighted_sites20.to_numpy()
        if cat == 'miRNA': F.text(x+0.4, yy-3.4, CL[c], 7.0)
        F.text(x+LBL-3.0, yy-3.4, cat, 7.0, ha='right', color=MUTED)
        q10, q25, q50, q75, q90 = np.percentile(z, [10, 25, 50, 75, 90])
        nz = z[z != 0]                                                     # per-gene strip of main Figure 5b (zeros counted at right)
        jit = (np.arange(len(nz)) % 5)/4*6.0-3.0
        for sel, cc in ((nz < 0, INK), (nz > 0, PROX)):
            F.points([ex(fl(v)) for v in nz[sel]], yy+jit[sel], ms='point', color=cc, alpha='secondary')
        nl, n0, ng = int((z < 0).sum()), int((z == 0).sum()), int((z > 0).sum())
        # every gene with W = 0 has no complete site (S6c.1: zero = 0 among genes with complete sites)
        assert n0 == int((gr[(gr.species == c[0]) & (gr.contrast == c[1]) & (gr.category == cat)].sites.astype(float) == 0).sum()), 'no-site genes'
        F.rect(cxw(-nl), yy-3.4, cxw(0)-cxw(-nl), 6.8, fc=INK); F.rect(cxw(0), yy-3.4, cxw(ng)-cxw(0), 6.8, fc=PROX)
        F.text(x+w-14.0, yy-3.4, str(n0), 7.0, ha='center')
        Ev[f'{CL[c]}|{cat}'] = dict(n=int(len(z)), q10=float(q10), q25=float(q25), median=float(q50), q75=float(q75), q90=float(q90),
                                    loss=nl, zero=n0, gain=ng, clipped=int(((z < -60) | (z > 40)).sum()))
        yy += rh
    yy += 3.0
F.line([ex(0)]*2, [y_top, yy-4], lw='ref', color=AXC, ls='ref'); F.line([cxw(0)]*2, [y_top, yy-4], lw='axis', color=AXC)
F.xaxis(ex, yy-2.0, [fl(-200), fl(-5), 0, fl(5), fl(200)], ['−200', '−5', '0', '5', '200'], title='Weighted predicted-site change')
F.xaxis(cxw, yy-2.0, [-200, 0, 200], ['200', '0', '200'], title='Genes')
# legend: one entry per direction, gene point and gene-count bar side by side (point = one gene, bar = number of genes)
kent = [(INK, 'Net loss'), (PROX, 'Net gain')]; kgap = 10.0
kw_ = sum(7.6+7.6+textwidth(t, 7.0) for _, t in kent)+kgap*(len(kent)-1)
kx = x+w/2-kw_/2; ky = y+h-7.4
for cc, t in kent:
    F.marker(kx+2.8, ky+3.8, 'o', 'mark', fc=cc); kx += 7.6
    F.rect(kx, ky+1.0, 5.6, 5.6, fc=cc); kx += 7.6
    kx += F.text(kx, ky, t, 7.0)+kgap
F.values['e'] = Ev
F.notes.append('e: form and colours of main Figure 5b: per-gene net weighted sites at 20 UMI on a symmetric log axis (sign · log10(1 + |x|)), '
               'net loss ink, net gain blue, genes without a complete site (W = 0 by construction; values key zero) counted at right; quantiles kept in values; predicted sites do not establish binding')

# ================================================================ f  50 vs 20 UMI
x, y, w, h = F.panel('f', X0+WE+GAP, Y3, RR-(X0+WE+GAP), H3)
d50 = (gr.weighted_sites50-gr.weighted_sites20)
fx = Lin(-4, 4, x+21.0, x+w-3.0); fy = Lin(0, 1, y+16.0+(fx.p1-fx.p0)/1.25, y+16.0)
F.lineplot(fx, fy, '50-20')
F.line([fx(0)]*2, [fy.p0, fy.p1], lw='ref', color=AXC, ls='ref')
Fv = {}
for cat, col, lsd in (('miRNA', INK, None), ('RBP', INK, 'dash')):
    v = np.sort(d50[gr.category == cat].to_numpy()); n = len(v)
    xs, ys = [fx.p0], [fy(0)]
    for k, t in enumerate(v):
        xx = fx.clip(t); xs += [xx, xx]; ys += [ys[-1], fy((k+1)/n)]
    xs.append(fx.p1); ys.append(ys[-1])
    F.line(xs, ys, lw='data', color=col, ls=lsd)
    Fv[cat] = dict(n=n, median=float(np.median(v)), q05=float(np.percentile(v, 5)), q95=float(np.percentile(v, 95)),
                   clipped=int((np.abs(v) > 4).sum()))
F.yaxis(fy, fx.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Cumulative fraction')
F.xaxis(fx, fy.p0+1.0, [-4, 0, 4], ['−4', '0', '4'], 'Difference in weighted sites')
F.key_row(None, y+h-7.4, [('line', (INK, None, 'data'), 'miRNA'), ('line', (INK, 'dash', 'data'), 'RBP')], gap=6.0, center=(fx.p0+fx.p1)/2)
F.values['f'] = Fv
F.save(OUT, D)
