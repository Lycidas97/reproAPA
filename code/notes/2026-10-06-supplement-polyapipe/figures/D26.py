"""Supplementary Figure D26 (printed number from supdraw.SNUM) — Trip12, the low-net-shift example of Figure 5 (R5),
polyApipe v2 primary. Intervals are numbered as in Figure 5f: 1 = 842 nt (nearer the proximal PAS), 2 = 706 nt.

a  inclusion of each strict interval in M and P of every paired library (20 UMI), both studies
b  native P–M ΔDWUI and MPRO (same-library estimator; contributors differ from the interval contributors)
c  interval ΔI at 20 and 50 UMI with the leave-one-library range (influence, not a confidence interval)
d  external long-read terminal-exon models (GSE282319) with the same complete intron chain: interval 2 is fully excluded
   by the short model (one strict pair); the short model keeps 2 nt of interval 1 (no strict pair); no sample detects both
e  archived AGO2 CLIP clusters overlapping the two intervals and the 68 miRNA seed annotations in them
f  every native PAS of Trip12 with its P and M usage: whole axis and the terminal-family zoom
g  exact decomposition of native ΔDWUI into the within-family position, family mass and outside-family terms
h  every predicted site in the two intervals by evidence tier (the three upper tiers are the Figure 5f rug)
"""
import sys, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D26'
D = Data(); OUT = outdir(FID)
X = ExtData(Path(__file__).resolve().parent.parent/'tools'/'external-sources-trip12.json')
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3, H4 = 122.0, 90.0, 120.0, 91.0   # 2026-10-12: panel titles and notes moved to axis titles or the caption
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; Y4 = Y3+H3+GAP; H = Y4+H4+MARGIN
F = Fig(H, figname(FID))
GID, SYM = 'ENSMUSG00000026219', 'Trip12'
RT = ('CIP', 'WOT')
IV = {1: (84701184, 84702026), 2: (84700478, 84701184)}          # Figure 5f numbering
IVC = {1: '#9CC3DD', 2: '#2D6A96'}
O_TE = 84702535                                                 # 5' end of the terminal exon (minus strand)
kb = lambda c: (O_TE-c)/1000
ivl = lambda k: IV[k][1]-IV[k][0]

rr = D.df('D26A', 'examples/recurrent_region_pairs.tsv.gz')
got = sorted({(int(a), int(b)) for a, b in rr[(rr.gene_symbol_a == SYM) & (rr.contrast == 'MP')][['shared_start0', 'shared_end0']].values})
assert got == sorted(IV.values()), got


def rmark(xm, yy, route, ms='mark'):
    F.marker(xm, yy, STUDY[route][2], ms, fc=STUDY[route][1])


# ================================================================ a  paired libraries
WA = 176.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
pl = D.df('D26A', 'examples/inclusion_paired_libraries.tsv.gz')
pl = pl[(pl.gene_id == GID) & (pl.kind == 'interval') & (pl.cutoff == 20)]
ay = Lin(0, 0.7, y+103.0, y+20.0)
gw = (w-24.0)/4
Av = {}
F.yaxis(ay, x+22.0, [0, 0.2, 0.4, 0.6], ['0', '20', '40', '60'], 'Interval inclusion (%)')
for j, (route, k) in enumerate([(r_, k_) for r_ in RT for k_ in (1, 2)]):
    z = pl[(pl.route == route) & (pl.start0 == IV[k][0]) & (pl.end0 == IV[k][1])]
    assert len(z) > 0
    x0 = x+24.0+j*gw; xm, xp = x0+gw*0.2, x0+gw*0.72
    for r in z.itertuples():
        F.line([xm, xp], [ay(r.reference_inclusion), ay(r.target_inclusion)], lw='thin', color=STUDY[route][1], alpha='secondary')
    mm, mp = float(z.reference_inclusion.mean()), float(z.target_inclusion.mean())
    F.line([xm, xp], [ay(mm), ay(mp)], lw='bold', color=INK)
    F.marker(xm, ay(mm), 'o', 'small', fc=SC['M']); F.marker(xp, ay(mp), 'o', 'small', fc=SC['P'])
    F.text(xm, ay.p0+2.0, 'M', 7.0, ha='center', color=SC['M']); F.text(xp, ay.p0+2.0, 'P', 7.0, ha='center', color=SC['P'])
    F.text(x0+gw/2, ay.p0+11.0, f'{k}', 7.5, ha='center', font='Bol')
    Av[f'{route}|I{k}'] = dict(n=int(len(z)), mean_M=mm, mean_P=mp)
F.line([x+24.0+2*gw]*2, [y+12.0, ay.p0+18.0], lw='ref', color=AXC, ls='ref')
for j, route in enumerate(RT):
    cx = x+24.0+gw*(2*j+1)
    F.text(cx, y+3.0, STUDY[route][0], 7.5, ha='center')
F.values['a'] = Av

# ================================================================ b  native metrics
WB = 92.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
ee = D.df('D26B', 'examples/example_endpoint_effects.tsv.gz')
ee = ee[(ee.gene_id == GID) & (ee.setting == 'supported50_p08') & (ee.contrast == 'P_M') & ee.route.isin(RT)].set_index('route')
Bv = {}
for row, (col, lab, lim, ticks, tl) in enumerate((('delta_DWUI', 'Native ΔDWUI', (-0.02, 0.02), [-0.02, 0, 0.02], ['−.02', '0', '.02']),
                                                   ('MPRO', 'MPRO', (0, 0.4), [0, 0.2, 0.4], ['0', '.2', '.4']))):
    top = y+3.0+row*63.0
    sx = Lin(lim[0], lim[1], x+14.0, x+w-6.0)
    for i, route in enumerate(RT):
        yy = top+15.0+i*13.0
        rmark(x+4.0, yy, route)
        r = ee.loc[route]
        F.line([sx(0) if col == 'delta_DWUI' else sx.p0, sx(r[col])], [yy, yy], lw='thin', color=AXC)
        F.marker(sx(r[col]), yy, 'o', 'mark', fc=INK)
        Bv[f'{route}|{col}'] = float(r[col]); Bv[f'{route}|pairs'] = int(r.n_pairs)
    if col == 'delta_DWUI':
        F.line([sx(0)]*2, [top+10.0, top+29.0], lw='ref', color=AXC, ls='ref')
    F.xaxis(sx, top+36.0, ticks, tl, title=lab)
F.values['b'] = Bv

# ================================================================ c  interval ΔI sensitivity
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
es = D.df('D26C', 'examples/inclusion_effect_sensitivity.tsv')
es = es[(es.gene_id == GID) & (es.kind == 'interval') & (es.contrast == 'MP')]
cx_ = Lin(-30, 0, x+44.0, x+w-4.0)
Cv = {}; yy = y+20.0; rp = 8.6
for k in (1, 2):
    for route in RT:
        for cut in (20, 50):
            r = es[(es.route == route) & (es.start0 == IV[k][0]) & (es.end0 == IV[k][1]) & (es.cutoff == cut)].iloc[0]
            if route == 'CIP' and cut == 20: F.text(x+0.4, yy-3.4, f'Interval {k}', 7.0, font='Bol')
            if cut == 20: rmark(x+30.0, yy+rp/2, route, 'small')
            F.text(x+42.0, yy-3.2, str(cut), 6.8, ha='right', color=MUTED)
            F.line([cx_.clip(r.leave_one_min*100), cx_.clip(r.leave_one_max*100)], [yy, yy], lw='bold', color=NULLBAR)
            F.marker(cx_(r.equal_library*100), yy, 'o', 'small', fc=INK)
            Cv[f'I{k}|{route}|{cut}'] = dict(dI_pp=float(r.equal_library*100), leave_one=[float(r.leave_one_min*100), float(r.leave_one_max*100)], libraries=int(r.n_libraries))
            yy += rp
        yy += 1.0
    yy += 3.0
F.xaxis(cx_, yy-4.0, [-30, -20, -10, 0], ['−30', '−20', '−10', '0'], title='Interval ΔI, P − M (pp)')
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'leave-one-library range')], gap=5.0, center=(cx_.p0+cx_.p1)/2)
F.values['c'] = Cv

# ================================================================ shared genomic frame for d and h
GX0, GX1 = -0.05, 2.25


def iv_shade(gx, ytop, ybot, labels=True):
    for k in (1, 2):
        a_, b_ = sorted((kb(IV[k][1]), kb(IV[k][0])))
        F.rect(gx(a_), ytop, gx(b_)-gx(a_), ybot-ytop, fc=BAND)
        if labels: F.text(gx((a_+b_)/2), ytop-8.4, f'{k}  ({ivl(k)} nt)', 7.0, ha='center', color=MUTED)
    F.line([gx(kb(IV[1][0]))]*2, [ytop, ybot], lw='ref', color=AXC, ls='ref')


# ================================================================ d  long-read models
WD = 300.0
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
ta = D.df('D26D', 'longread/transcript_audit.tsv')
cp = D.df('D26D', 'longread/compatible_isoform_pairs.tsv')
D.df('D26D', 'longread/structural_summary.tsv')
ba = X.df('longread_boundary')
LONG, SHORT = 'Trip12.G1389.16', 'Trip12.G1389.149'
gx = Lin(GX0, GX1, x+70.0, x+w-1.0)
iv_shade(gx, y+16.0, y+60.0)
Dv = {}
for i, tid in enumerate((LONG, SHORT)):
    a = ta[ta.transcript_id == tid].iloc[0]
    ex_ = json.loads(a.exons_mm39); lo, hi = min(ex_, key=lambda e: e[0])
    ry = y+28.0+i*20.0
    F.text(x+0.4, ry-3.4, ('Long model ' if i == 0 else 'Short model ')+tid.split('.', 1)[1], 7.0)
    a0, a1 = sorted((kb(hi), kb(lo)))
    F.rect(gx.clip(a0), ry-4.0, gx.clip(a1)-gx.clip(a0), 8.0, fc=MUTED)
    Dv[tid] = dict(terminal_exon=[int(lo), int(hi)], max_TPM=float(a.max_TPM), samples_TPM_positive=int(a.n_samples_TPM_positive))
for k in (1, 2):
    r = ba[ba.feature_id.str.contains(f'|{IV[k][0]}|{IV[k][1]}|', regex=False)].iloc[0]
    a_, b_ = sorted((kb(IV[k][1]), kb(IV[k][0])))
    txt = f'kept {int(r.proximal_overlap_nt)} nt' if not r.strict_full_exclusion else 'fully excluded'
    F.text(gx((a_+b_)/2), y+56.0-3.4+8.0, txt, 6.8, ha='center', color=INK if r.strict_full_exclusion else MUTED)
    Dv[f'I{k}'] = dict(length=int(r.length_nt), short_model_overlap_nt=int(r.proximal_overlap_nt), strict=bool(r.strict_full_exclusion),
                       both_TPM_positive_samples=int(r.both_TPM_positive_samples))
assert int(ba.both_TPM_positive_samples.max()) == 0
yb = F.xaxis(gx, y+70.0, [0, 0.5, 1, 1.5, 2], ['0', '.5', '1', '1.5', '2'])
F.text((gx.p0+gx.p1)/2, yb+1.6, 'Distance from the terminal-exon start (kb)', 8.0, ha='center')
F.values['d'] = Dv

# ================================================================ e  AGO2
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, RR-(X0+WD+GAP), H2)
ov = X.df('ago2_overlap'); X.df('ago2_archive'); X.df('ago2_units')
ov['kind'] = ['interval' if '|shared_region|' in f else ('mirna' if '|site|miRNA|' in f else 'other') for f in ov.feature_id]
TR = [('WT1A', 'WT1 A'), ('WT1B', 'WT1 B'), ('WT2', 'WT2'), ('WT_merged', 'WT merged'), ('KO_merged', 'Dicer-KO merged')]
cx0 = x+58.0; cw = (x+w-1.0-cx0)/3; rh = 14.0; e0 = y+14.0
for j, lab in enumerate(('Int. 1', 'Int. 2', 'Seeds')):
    F.text(cx0+cw*(j+0.5), y+3.0, lab, 7.0, ha='center')
nmi = int((ov[ov.label == 'WT2'].kind == 'mirna').sum())
Ev = {}
for i, (lab, name) in enumerate(TR):
    yy = e0+i*rh+(3.0 if i >= 3 else 0.0)
    F.text(cx0-2.0, yy+rh/2-3.4, name, 7.0, ha='right')
    z = ov[ov.label == lab]
    vals = []
    for k in (1, 2):
        q = z[(z.kind == 'interval') & z.feature_id.str.contains(f'|{IV[k][0]}|{IV[k][1]}|', regex=False)].iloc[0]
        vals.append((int(q.overlap_bp_sum), f'{int(q.overlap_bp_sum)} nt' if q.overlap_bp_sum else '0'))
    mi = z[z.kind == 'mirna']; hits = int((mi.overlap_clusters > 0).sum())
    vals.append((hits, f'{hits}/{len(mi)}'))
    for j, (v, t_) in enumerate(vals):
        fc = SHADE if v > 0 else 'white'
        F.rect(cx0+cw*j+0.4, yy+0.4, cw-0.8, rh-0.8, fc=fc, ec=AXC, lw='axis')
        F.text(cx0+cw*(j+0.5), yy+rh/2-3.2, t_, 6.8, ha='center', color=INK if v > 0 else MUTED)
    Ev[lab] = dict(I1_bp=vals[0][0], I2_bp=vals[1][0], mirna_annotations_hit=hits, mirna_annotations=int(len(mi)))
F.values['e'] = Ev
assert nmi == 68

# ================================================================ f  native PAS
WF = 300.0
x, y, w, h = F.panel('f', X0, Y3, WF, H3)
ns = D.df('D26F', 'geometry/native-sites.tsv.gz')
ns = ns[ns.gene_id == GID]
ad = D.df('D26F', 'geometry/axis-decomposition.tsv')
ad = ad[ad.gene_id == GID].set_index('route')
oN = int(ns.cut_interbase.max())
fy = Lin(0, 0.7, y+78.0, y+18.0)
fx1 = Lin(-1, 64, x+24.0, x+108.0)
fx2 = Lin(GX0, 3.8, x+124.0, x+w-3.0)
F.text((fx1.p0+fx1.p1)/2, y+3.0, 'Whole PAS axis', 7.5, ha='center')
F.text((fx2.p0+fx2.p1)/2, y+3.0, 'Terminal family', 7.5, ha='center')
fam = (ad.family_start0.iloc[0], ad.family_end0.iloc[0])
F.rect(fx1((oN-fam[1])/1000), fy.p1, max(fx1((oN-fam[0])/1000)-fx1((oN-fam[1])/1000), 1.0), fy.p0-fy.p1, fc=BAND)
for k in (1, 2):
    a_, b_ = sorted((kb(IV[k][1]), kb(IV[k][0])))
    F.rect(fx2(a_), fy.p1, fx2(b_)-fx2(a_), fy.p0-fy.p1, fc=BAND)
    F.text(fx2((a_+b_)/2), fy.p1-1.0+0.0, str(k), 7.0, ha='center', color=MUTED)
F.line([fx2(kb(IV[1][0]))]*2, [fy.p1, fy.p0], lw='data', color='white')
Fv = {}
for route, dxo in (('CIP', -0.6), ('WOT', 0.6)):
    z = ns[ns.route == route]
    for r in z.itertuples():
        for sx, xv in ((fx1, (oN-r.cut_interbase)/1000), (fx2, kb(r.cut_interbase))):
            if sx is fx2 and not r.terminal_family: continue
            xx = sx.clip(xv)+dxo
            F.line([xx, xx], [fy(r.M_usage), fy(r.P_usage)], lw='thin', color=STUDY[route][1])
            F.marker(xx, fy(r.P_usage), STUDY[route][2], 'small', fc=STUDY[route][1])
            F.marker(xx, fy(r.M_usage), STUDY[route][2], 'small', fc='white', ec=STUDY[route][1])
    Fv[route] = [dict(cut=int(r.cut_interbase), P=float(r.P_usage), M=float(r.M_usage), family=bool(r.terminal_family)) for r in z.itertuples()]
F.yaxis(fy, fx1.p0, [0, 0.2, 0.4, 0.6], ['0', '.2', '.4', '.6'], 'PAS usage')
F.yaxis(fy, fx2.p0, [0, 0.2, 0.4, 0.6], None, labels_off=True)
yb = F.xaxis(fx1, fy.p0+1.0, [0, 30, 60], ['0', '30', '60'], 'from proximal PAS (kb)')
F.xaxis(fx2, fy.p0+1.0, [0, 1, 2, 3], None, 'from terminal-exon start (kb)')
F.key_row(None, y+h-7.4, [('marker', ('s', STUDY['CIP'][1]), 'Guo2019–CIP'), ('marker', ('o', STUDY['WOT'][1]), 'Schiebinger2019–2i'),
                          ('marker', ('o', INK), 'P'), ('marker', ('o', 'white', INK), 'M')], gap=6.0, center=x+w/2)
F.values['f'] = Fv
F.notes.append('f: every current native PAS, zero-usage sites included; shade = frozen terminal family (left) and the two strict intervals (right)')

# ================================================================ g  decomposition
x, y, w, h = F.panel('g', X0+WF+GAP, Y3, RR-(X0+WF+GAP), H3)
COMP = [('within_family_DWUI', SEQ_STOPS[2], 'Family position'), ('family_mass_DWUI', SEQ_STOPS[1], 'Family mass'), ('nonterminal_DWUI', OUTSIDE, 'Outside family')]
gx_ = Lin(-10, 6, x+20.0, x+w-4.0)
Gv = {}
for i, route in enumerate(RT):
    r = ad.loc[route]; yy = y+26.0+i*22.0
    rmark(x+4.0, yy, route)
    pos = neg = 0.0
    for col, cc, lab in COMP:
        v = float(r[col])*1000
        left = pos if v >= 0 else neg+v
        F.rect(gx_(left), yy-5.0, gx_(left+abs(v))-gx_(left), 10.0, fc=cc)
        if v >= 0: pos += v
        else: neg += v
    F.marker(gx_(r.native_delta_DWUI*1000), yy, 'D', 'mark', fc='white', ec=INK)
    resid = float(r.native_delta_DWUI-sum(r[c] for c, _, _ in COMP))
    assert abs(resid) < 1e-9, resid
    Gv[route] = {c: float(r[c]) for c, _, _ in COMP} | dict(native=float(r.native_delta_DWUI), terminal_shift_nt=float(r.terminal_shift_nt))
F.line([gx_(0)]*2, [y+18.0, y+56.0], lw='ref', color=AXC, ls='ref')
F.xaxis(gx_, y+62.0, [-8, -4, 0, 4], ['−8', '−4', '0', '4'], title='ΔDWUI × 1,000')
for k, (col, cc, lab) in enumerate(COMP):
    F.key_row(x+6.0, y+86.0+k*11.4, [('block', cc, lab)], size=7.0)
F.key_row(x+82.0, y+86.0, [('marker', ('D', 'white', INK), 'native net')], size=7.0)
F.values['g'] = Gv

# ================================================================ h  site tiers
x, y, w, h = F.panel('h', X0, Y4, RR-X0, H4)
sv = X.df('site_evidence'); sm = X.df('site_summary')
sv = sv[sv.route == 'CIP']
TIERS = [('PSC_miRNA', 'PSC miRNA', INK, True), ('direct_RBP_motif', 'Direct RBP', SC['M'], True),
         ('transferred_RBP_P_expressed', 'Transferred RBP', OUTSIDE, True), ('other_miRNA', 'Other miRNA', AXC, False),
         ('transferred_RBP_not_P_expressed_both_routes', 'Transf. RBP, not P-expressed', AXC, False)]
hx = Lin(GX0, GX1, x+120.0, x+w-44.0)
iv_shade(hx, y+16.0, y+68.0, labels=False)
for k in (1, 2):
    a_, b_ = sorted((kb(IV[k][1]), kb(IV[k][0])))
    F.text(hx((a_+b_)/2), y+3.0, f'Interval {k} ({ivl(k)} nt)', 7.5, ha='center')
F.text(x+w-20.5, y+3.0, 'Annot. / pos.', 7.0, ha='center', color=MUTED)
Hv = {}
for i, (tier, lab, cc, shown) in enumerate(TIERS):
    yy = y+21.0+i*10.4
    F.text(x+0.4, yy-3.4, lab, 7.0, color=INK if shown else MUTED)
    z = sv[sv.evidence_tier == tier]
    pos_ = sorted({(int(a), int(b)) for a, b in z[['start0', 'end0']].values})
    for a, b in pos_:
        xm = hx(kb((a+b)/2))
        F.line([xm, xm], [yy-3.6, yy+3.6], lw='thin' if shown else 'axis', color=cc)
    s_ = sm[(sm.route == 'CIP') & (sm.evidence_tier == tier)].iloc[0]
    assert int(s_.distinct_coordinates) == len(pos_) and int(s_.annotations) == len(z)
    F.text(x+w-20.5, yy-3.4, f'{int(s_.annotations)} / {len(pos_)}', 7.0, ha='center')
    Hv[tier] = dict(annotations=int(len(z)), positions=len(pos_), in_figure5f=shown)
yb = F.xaxis(hx, y+21.0+4*10.4+7.0, [0.5, 1, 1.5, 2], ['.5', '1', '1.5', '2'])
F.text((hx.p0+hx.p1)/2, yb+1.6, 'Distance from the terminal-exon start (kb)', 8.0, ha='center')
# 2026-10-12: tier note (Figure 4d rug; predicted, not binding) moved to the caption
F.values['h'] = Hv
F.notes.append('h: Guo2019–CIP rows; Schiebinger2019–2i has identical site coordinates and tiers (checked)')
F.save(OUT, [D, X])
