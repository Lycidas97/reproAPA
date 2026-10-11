"""Supplementary Figure D16 (printed number from supdraw.SNUM) — cross-species sharing: sensitivity, direction and continuous
correspondence (R4), polyApipe v2 primary.

a  shared proximal candidates (P–S, P–M) against the joint measurement/PAS null at four thresholds in the supported and
   expanded layers (the primary .02 / supported result is Figure 4c)
b  same-direction genes for the four human x mouse study pairs, split into shared proximal and shared distal (Q5 null)
c  all 80 direction-difference tests (v2): complete-family q
d  human against mouse observed-state ΔDWUI for every measurable gene pair, by state
e  Spearman ρ with gene-resampling 95% intervals for all 24 study pairs, R, E, P and M (Figure 4b adds M–S)
f  conditional-injection power of the sharing test (L1 capacity and L3 full nulls)
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D16'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 112.0, 272.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
HUM = {'H1': 'GSE118258', 'H2': 'GSE150311', 'H3': '178325_0618', 'H4': 'GSE221739'}
MOU = {'W': 'WOT', 'C': 'CIP'}
PAIRS = [('GSE118258', 'GSE103221', 'GSE118258', 'CIP'), ('GSE118258', 'GSE115943', 'GSE118258', 'WOT'),
         ('GSE178325', 'GSE103221', '178325_0618', 'CIP'), ('GSE178325', 'GSE115943', '178325_0618', 'WOT')]


# ================================================================ a  threshold and layer sensitivity
WA = 214.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
ot = D.df('D16A', 'cross_overlap/overlap_threshold_layer_sensitivity.tsv')
pw = (w-48-10-4)/2
Av = {}
for j, (ct, cl) in enumerate((('PS', 'P–S'), ('PM', 'P–M'))):
    sx = Lin(0, 70, x+48+(pw+10)*j, x+48+(pw+10)*j+pw)
    F.text((sx.p0+sx.p1)/2, y+3.0, f'{cl}, shared proximal', 7.5, ha='center')
    r = 0
    for li, (layer, ll) in enumerate((('supported50', 'Supported'), ('expanded_primary', 'Expanded'))):
        for th in (0.01, 0.02, 0.03, 0.05):
            yy = y+20.0+r*8.6+li*3.0
            q = ot[(ot.contrast == ct) & (ot.layer == layer) & (abs(ot.threshold-th) < 1e-9)].iloc[0]
            if j == 0:
                if th == 0.01: F.text(x+0.4, yy-3.4, ll, 7.0)
                F.text(sx.p0-3.0, yy-3.4, f'{th:g}'.replace('0.', '.'), 7.0, ha='right', color=INK if (th == 0.02 and layer == 'supported50') else MUTED,
                       font='Bol' if (th == 0.02 and layer == 'supported50') else 'Reg')
            F.line([sx(q.null_q025), sx(q.null_q975)], [yy, yy], lw='bold', color=NULLBAR)
            sig = q.q_16 < 0.05
            F.marker(sx(q.observed), yy, 'o', 'mark', fc=INK if sig else 'white', ec=INK)
            Av[f'{ct}|{layer}|{th}'] = dict(observed=int(q.observed), null=[float(q.null_q025), float(q.null_q975)], q16=float(q.q_16))
            r += 1
    F.xaxis(sx, y+20.0+7*8.6+3.0+5.4, [0, 35, 70])
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'q16 < .05'), ('marker', ('o', 'white', INK), 'q16 ≥ .05'),
                          ('line', (NULLBAR, None, 'bold'), 'Null 95%')], gap=6.0, center=x+w/2)
F.values['a'] = Av

# ================================================================ b  same-direction genes by pair
WB = 108.0
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, WB, H1)
dd = D.df('D16B', 'direction_difference/v2_original_direction_decomposition.tsv')
dd = dd[(dd.metric == 'WPP') & (abs(dd.threshold-0.02) < 1e-9) & (dd.model == 'Q5')]
bx = Lin(0, 220, x+26, x+w-1.5)
F.text((bx.p0+bx.p1)/2, y+3.0, 'Same direction', 7.5, ha='center')
Bv = {}
for i, (hs, ms, hr, mr) in enumerate(PAIRS):
    yy = y+21.0+i*16.0
    F.study_pair(x+0.4, yy, hr, mr)
    q = dd[(dd.human_study == hs) & (dd.mouse_study == ms)].set_index('test')
    neg, pos = int(q.loc['shared_negative'].observed), int(q.loc['shared_positive'].observed)
    F.rect(bx(0), yy-4.2, bx(neg)-bx(0), 8.4, fc=PROX); F.rect(bx(neg), yy-4.2, bx(neg+pos)-bx(neg), 8.4, fc=DIST)
    s_ = q.loc['same_direction']
    F.line([bx(s_.null_p025), bx(s_.null_p975)], [yy+6.4]*2, lw='bold', color=NULLBAR)
    Bv[f'{hr}|{mr}'] = dict(proximal=neg, distal=pos, same=int(s_.observed), null=[float(s_.null_p025), float(s_.null_p975)],
                            q_joint_all_56=float(s_.q_joint_all_56))
F.xaxis(bx, y+21.0+3*16.0+9.0, [0, 100, 200])
F.key_row(None, y+h-17.4, [('block', PROX, 'Proximal'), ('block', DIST, 'Distal')], gap=5.0, center=x+w/2)
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'Null 95%')], gap=5.0, center=x+w/2)
F.values['b'] = Bv

# ================================================================ c  direction tests
x, y, w, h = F.panel('c', X0+WA+GAP+WB+GAP, Y1, RR-(X0+WA+GAP+WB+GAP), H1)
dt = D.df('D16C', 'direction_difference/direction_contrast_tests.tsv')
op = D.df('D16C', 'direction_difference/original_pair_direction_sensitivity.tsv')
FAMS = [('pair', 'distal_minus_proximal_excess_rate', dt, 'Excess rate'), ('pair', 'log_distal_vs_proximal_relative_enrichment', dt, 'Log enrichment'),
        ('unique_gene_pooled', None, dt, 'Pooled genes'), ('pair_heterogeneity', None, dt, 'Between pairs'), ('pair', None, op, 'Pair-specific null')]
qx = Lin(0, 1, x+50, x+w-3)
F.text(x+w/2+6, y+3.0, 'q, 80 tests', 7.5, ha='center')
Cv = {}; rng = np.random.default_rng(16); tot = 0
for i, (sc_, tst, src, lab) in enumerate(FAMS):
    yy = y+21.0+i*14.4
    q = src[src.scope == sc_]
    if tst: q = q[q.test == tst]
    F.text(x+0.4, yy-3.4, lab, 7.0)
    F.points([qx(v) for v in q.q_complete_80], yy+rng.uniform(-3.2, 3.2, len(q)), ms='small', color=INK, alpha='secondary')
    Cv[lab] = dict(n=int(len(q)), min_q=float(q.q_complete_80.min())); tot += len(q)
F.line([qx(0.05)]*2, [y+15, y+21.0+4*14.4+5], lw='ref', color=AXC, ls='ref')
F.xaxis(qx, y+21.0+4*14.4+7.0, [0, 0.5, 1], ['0', '.5', '1'])
F.text(x+w/2, y+h-7.4, f'{tot} tests; minimum q = {min(v["min_q"] for v in Cv.values()):.2f}'.replace('0.', '.'), 7.0, ha='center', color=MUTED)
assert tot == 80, tot
F.values['c'] = Cv

# ================================================================ d  gene-pair scatter
WD = 150.0; HD = 162.0
x, y, w, h = F.panel('d', X0, Y2, WD, HD)
gp = D.df('D16D', 'supplement_statistics/continuous-gene-pairs.tsv.gz')
gp = gp[gp.layer == 'supported50']
sx = Lin(-0.6, 0.6, x+30, x+w-3); sy = Lin(-0.6, 0.6, y+16+(sx.p1-sx.p0), y+16)
F.text((sx.p0+sx.p1)/2, y+3.0, f'{len(gp):,} gene pairs', 7.5, ha='center')
F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref'); F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
F.yaxis(sy, sx.p0, [-0.6, -0.3, 0, 0.3, 0.6], ['−.6', '−.3', '0', '.3', '.6'], 'Mouse ΔDWUI vs S')
for st in ['R', 'E', 'P', 'M']:
    z = gp[gp.state == st]
    F.points([sx.clip(v) for v in z.human_effect], [sy.clip(v) for v in z.mouse_effect], color=SC[st], alpha='dense', rasterized=True)
F.xaxis(sx, sy.p0+1.0, [-0.6, -0.3, 0, 0.3, 0.6], ['−.6', '−.3', '0', '.3', '.6'], 'Human ΔDWUI vs S')
F.key_row(None, y+h-7.4, [('marker', ('o', SC[s]), f'{s}–S') for s in 'REPM' if (gp.state == s).any()], gap=5.0, center=(sx.p0+sx.p1)/2)
F.values['d'] = dict(pairs=int(len(gp)), by_state={s: int((gp.state == s).sum()) for s in 'REPM'})
F.notes.append('d: supported layer; points beyond ±.6 drawn at the limit; no pooled correlation is reported for this mixed pool')

# ================================================================ e  all study pairs
x, y, w, h = F.panel('e', X0+WD+GAP, Y2, RR-(X0+WD+GAP), H2)
cp = D.df('D16E', 'cross_continuous/continuous_pair_summary.tsv')
cp = cp[cp.layer == 'supported50']
ex = Lin(-0.15, 0.25, x+46, x+w-23)
F.text(x+w, y+3.0, 'Genes', 7.0, ha='right', color=MUTED)
F.text((ex.p0+ex.p1)/2, y+3.0, 'Human × mouse study pairs, supported layer', 7.5, ha='center')
Ev = {}; yy = y+20.0; rh = 9.0
for st in ['R', 'E', 'P', 'M']:
    q = cp[cp.state == st].sort_values(['human', 'mouse'])
    F.text(x+0.4, yy+rh*(len(q)-1)/2-3.4, f'{st}–S', 7.5, color=SC[st], font='Bol')
    for r in q.itertuples():
        F.study_pair(x+20.0, yy, HUM[r.human], MOU[r.mouse])
        F.line([ex.clip(r.lo), ex.clip(r.hi)], [yy, yy], lw='bold', color=NULLBAR)
        F.marker(ex(r.rho), yy, 'o', 'small', fc=INK)
        F.text(x+w, yy-3.4, f'{int(r.n):,}', 6.8, ha='right', color=MUTED)
        Ev[f'{st}|{HUM[r.human]}|{MOU[r.mouse]}'] = dict(n=int(r.n), rho=float(r.rho), ci=[float(r.lo), float(r.hi)])
        yy += rh
    yy += 3.0
F.line([ex(0)]*2, [y+15, yy-3], lw='ref', color=AXC, ls='ref')
F.xaxis(ex, yy-2.0, [-0.1, 0, 0.1, 0.2], ['−.1', '0', '.1', '.2'])
F.key_row(None, y+h-7.4, [('block', STUDY[r][1], STUDY[r][0].replace('Schiebinger2019–', '').replace('Guo2019–', ''))
                          for r in ['GSE118258', 'GSE150311', '178325_0618', 'GSE221739', 'CIP', 'WOT']], gap=5.0, center=x+w/2)
F.values['e'] = Ev

# ================================================================ f  power
x, y, w, h = F.panel('f', X0, Y2+HD+GAP, WD, H2-HD-GAP)
pwr = D.df('D16F', 'cross_controls/cross_species_power.tsv')
pwr = pwr[(pwr.experiment == 'conditional_injection') & (pwr.generator == 'L3_full') & (abs(pwr.alpha-0.05) < 1e-9)]
fx = Lin(1, 8, x+26, x+w-4); fy = Lin(0, 1, y+16+(fx.p1-fx.p0)/1.98, y+16)
F.lineplot(fx, fy, 'power')
F.text((fx.p0+fx.p1)/2, y+3.0, 'Power by injected odds ratio', 7.5, ha='center')
F.yaxis(fy, fx.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Power', title_y=(fy.p0+fy.p1)/2+3)
Fv = {}
for ct, col in (('PS', INK), ('PM', MUTED)):
    for tm, ls in (('L1_capacity', 'dash'), ('L3_full', None)):
        q = pwr[(pwr.contrast == ct) & (pwr.tested_model == tm)].sort_values('odds_ratio')
        F.line([fx(v) for v in q.odds_ratio], [fy(v) for v in q.rejection_rate], lw='data', color=col, ls=ls)
        for o_, v in zip(q.odds_ratio, q.rejection_rate): F.marker(fx(o_), fy(v), 'o', 'small', fc=col)
        Fv[f'{ct}|{tm}'] = {int(o_): float(v) for o_, v in zip(q.odds_ratio, q.rejection_rate)}
F.xaxis(fx, fy.p0+1.0, [1, 2, 4, 8])
F.key_row(None, y+h-7.4, [('line', (INK, None, 'data'), 'P–S'), ('line', (MUTED, None, 'data'), 'P–M'), ('line', (INK, 'dash', 'data'), 'L1')],
          gap=5.0, center=x+w/2)
F.values['f'] = Fv
F.save(OUT, D)
