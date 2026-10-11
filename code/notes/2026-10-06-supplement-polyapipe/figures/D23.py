"""Supplementary Figure D23 (printed number from supdraw.SNUM) — predicted-site density in changed versus matched stable
3' UTR windows and the detection limit of this test (R5), polyApipe v2 primary.

a  per-gene density difference (changed − matched stable, sites/kb) for each species x contrast x site-category cohort,
   gene-equal mean, gene-bootstrap 95% interval where >= 10 genes, and permutation-test q over the 8 tests
b  detection power of the matched-gene test against injected site density (mouse P–M, the only cohort with >= 10 genes)
c  GC balance of the matched windows (lost window against its control)
d  lost windows by matching outcome
e  background ladder: mean density difference under nine background definitions for all eight cohorts
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D23'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 132.0, 150.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
CX = [('human', 'SP'), ('human', 'MP'), ('mouse', 'SP'), ('mouse', 'MP')]
CL = {c: f'{c[0].capitalize()} {c[1][1]}–{c[1][0]}' for c in CX}
CC = {c: (SC['S'] if c[1] == 'SP' else SC['M']) for c in CX}
CM = {c: ('o' if c[0] == 'human' else 's') for c in CX}
CAT = {'miRNA': INK, 'RBP': SEQ_STOPS[1]}

# ================================================================ a  density differences
WA = 262.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
gd = D.df('D23A', 'motif_background/gene_density.tsv')
ds = D.df('D23A', 'motif_background/density_summary.tsv')
ax_ = Lin(-80, 130, x+73.0, x+w-46.0)
F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Changed − stable sites per kb', 7.5, ha='center')
F.text(x+w-32.0, y+3.0, 'Genes', 7.0, ha='center', color=MUTED); F.text(x+w-9.0, y+3.0, 'q', 7.0, ha='center', color=MUTED)
Av = {}; yy = y+20.0; rp = 11.0
y_top = yy-6
rng = np.random.default_rng(23)
for c in CX:
    for k, cat in enumerate(('miRNA', 'RBP')):
        if k == 0:
            F.marker(x+3.0, yy, CM[c], 'mark', fc=CC[c]); F.text(x+8.4, yy-3.4, CL[c], 7.0)
        F.text(x+70.0, yy-3.4, cat, 7.0, ha='right', color=MUTED)
        z = gd[(gd.species == c[0]) & (gd.contrast == c[1]) & (gd.category == cat)].density_difference.to_numpy()
        s_ = ds[(ds.species == c[0]) & (ds.contrast == c[1]) & (ds.category == cat)].iloc[0]
        assert len(z) == int(s_.n_genes)
        F.points([ax_.clip(v) for v in z], yy+rng.uniform(-2.6, 2.6, len(z)), ms='small', color=CAT[cat], alpha='secondary')
        if np.isfinite(s_.lo):
            F.line([ax_(s_.lo), ax_(s_.hi)], [yy, yy], lw='bold', color=INK)
        F.line([ax_(s_.mean_difference)]*2, [yy-4.0, yy+4.0], lw='bold', color=INK)
        F.text(x+w-32.0, yy-3.4, str(int(s_.n_genes)), 7.0, ha='center')
        F.text(x+w-9.0, yy-3.4, ('1' if s_.q_8_tests > 0.999 else qfmt(s_.q_8_tests)) if s_.status == 'matched_gene_permutation' else '–', 7.0, ha='center',
               color=INK if s_.status == 'matched_gene_permutation' else MUTED)
        Av[f'{CL[c]}|{cat}'] = dict(n=int(s_.n_genes), mean=float(s_.mean_difference), ci=[float(s_.lo), float(s_.hi)] if np.isfinite(s_.lo) else None,
                                   q8=float(s_.q_8_tests), status=s_.status, clipped=int(((z < -80) | (z > 130)).sum()))
        yy += rp
    yy += 2.0
F.line([ax_(0)]*2, [y_top, yy-6], lw='ref', color=AXC, ls='ref')
F.xaxis(ax_, yy-4.0, [-50, 0, 50, 100], ['−50', '0', '50', '100'])
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'gene'), ('line', (INK, None, 'bold'), 'mean | 95% interval'),
                          ('line', (INK, None, 'ref'), '– fewer than 10 genes')][:2], gap=7.0, center=(ax_.p0+ax_.p1)/2)
F.values['a'] = Av
F.notes.append('a: q over the 8 cohort tests; cohorts with fewer than 10 genes are descriptive only (no interval, q shown as –); points jittered vertically')

# ================================================================ b  power
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
pw = D.df('D23B', 'motif_controls/motif_power.tsv')
pw = pw[(abs(pw.alpha-0.05) < 1e-9) & (pw.species == 'mouse') & (pw.contrast == 'MP')]
bx = Lin(0, 20, x+24.0, x+w-4.0); by = Lin(0, 1, y+16.0+(bx.p1-bx.p0)/1.9, y+16.0)
F.lineplot(bx, by, 'power')
F.text((bx.p0+bx.p1)/2, y+3.0, 'Mouse P–M, α = .05', 7.5, ha='center')
Bv = {}
for cat in ('miRNA', 'RBP'):
    for prev, lsd in (('all', None), ('quarter', 'dash')):
        z = pw[(pw.category == cat) & (pw.prevalence == prev)].sort_values('injected_sites_per_kb')
        F.line([bx(v) for v in z.injected_sites_per_kb], [by(v) for v in z.rejection_rate], lw='bold', color=CAT[cat], ls=lsd)
        for v, r in zip(z.injected_sites_per_kb, z.rejection_rate):
            F.marker(bx(v), by(r), 'o', 'small', fc=CAT[cat] if prev == 'all' else 'white', ec=CAT[cat])
        Bv[f'{cat}|{prev}'] = {int(v): float(r) for v, r in zip(z.injected_sites_per_kb, z.rejection_rate)}
F.yaxis(by, bx.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Detection probability')
F.xaxis(bx, by.p0+1.0, [0, 5, 10, 20], title='Injected sites/kb per gene')
F.key_row(None, y+h-17.4, [('line', (CAT['miRNA'], None, 'bold'), 'miRNA'), ('line', (CAT['RBP'], None, 'bold'), 'RBP')], gap=6.0, center=(bx.p0+bx.p1)/2)
F.key_row(None, y+h-7.4, [('line', (INK, None, 'bold'), 'all genes'), ('line', (INK, 'dash', 'bold'), 'quarter')], gap=6.0, center=(bx.p0+bx.p1)/2)
F.values['b'] = Bv
F.notes.append('b: the other three cohorts have fewer than 10 genes, so the test is not run and power is 0 by design; 1,000 simulations per point')

# ================================================================ c  GC balance
WC = 128.0
x, y, w, h = F.panel('c', X0, Y2, WC, H2)
mw = D.df('D23C', 'motif_background/matched_windows.tsv.gz').drop_duplicates(['species', 'contrast', 'gene_id', 'start0', 'end0'])
gx = Lin(0.2, 0.65, x+24.0, x+w-3.0); gy = Lin(0.2, 0.65, y+16.0+88.0, y+16.0)
F.text((gx.p0+gx.p1)/2, y+3.0, f'{len(mw)} matched windows', 7.5, ha='center')
F.line([gx.p0, gx.p1], [gy.p0, gy.p1], lw='ref', color=AXC, ls='ref')
Cv = {}
for c in CX:
    z = mw[(mw.species == c[0]) & (mw.contrast == c[1])]
    for a_, b_ in zip(z.GC_change, z.GC_control):
        F.marker(gx.clip(a_), gy.clip(b_), CM[c], 'small', fc=CC[c])
    Cv[CL[c]] = dict(n=int(len(z)), max_abs_GC_difference=float((z.GC_change-z.GC_control).abs().max()) if len(z) else None)
F.yaxis(gy, gx.p0, [0.2, 0.4, 0.6], ['.2', '.4', '.6'], 'Control GC')
F.xaxis(gx, gy.p0+1.0, [0.2, 0.4, 0.6], ['.2', '.4', '.6'], 'Lost-window GC')
F.key_row(None, y+h-17.4, [('marker', ('o', CC[('human', 'SP')]), 'P–S'), ('marker', ('o', CC[('human', 'MP')]), 'P–M')], gap=6.0, center=(gx.p0+gx.p1)/2)
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'Human'), ('marker', ('s', INK), 'Mouse')], gap=6.0, center=(gx.p0+gx.p1)/2)
F.values['c'] = Cv

# ================================================================ d  matching outcome
WD = 130.0
x, y, w, h = F.panel('d', X0+WC+GAP, Y2, WD, H2)
at = D.df('D23D', 'motif_background/motif_attrition.tsv')
ORD = [('no_equal_length_window', MISSING, 'No equal-length window'), ('no_GC_match', NULLBAR, 'No GC match'),
       ('no_joint_position_match', NULLW, 'No position match'), ('matched', INK, 'Matched')]
dx = Lin(0, 680, x+44.0, x+w-14.0)
F.text((dx.p0+dx.p1)/2, y+3.0, 'Lost windows', 7.5, ha='center')
Dv = {}
for i, c in enumerate(CX):
    yy = y+24.0+i*22.0
    F.text(x+0.4, yy-3.4, CL[c], 7.0)
    z = at[(at.species == c[0]) & (at.contrast == c[1])].status.value_counts()
    left = 0
    for st_, col, lab in ORD:
        n = int(z.get(st_, 0))
        if n: F.rect(dx(left), yy-7.5, dx(left+n)-dx(left), 15.0, fc=col)
        left += n
    assert left == int(z.sum())
    F.text(dx(left)+1.4, yy-3.2, str(left), 6.8, color=MUTED)
    Dv[CL[c]] = {st_: int(z.get(st_, 0)) for st_, _, _ in ORD}
F.xaxis(dx, y+24.0+3*22.0+7.0, [0, 300, 600])
ky = y+24.0+3*22.0+22.0
for k, (st_, col, lab) in enumerate(ORD):
    F.key_row(x+4.0, ky+k*9.0, [('block', col, lab)], size=7.0)
F.values['d'] = Dv

# ================================================================ e  background ladder
x, y, w, h = F.panel('e', X0+WC+GAP+WD+GAP, Y2, RR-(X0+WC+GAP+WD+GAP), H2)
bl = D.df('D23E', 'motif_controls/background_ladder.tsv')
D.df('D23E', 'motif_controls/grid_summary.tsv')         # full 240-test grid (receipt; linked in the caption)
ROWS = [('same_gene_original_fixed', 'same_gene', 'Fixed'), ('six_background_common', 'same_gene', 'Same gene'),
        ('six_background_common', 'cross_gene', 'Cross gene')]
LEV = [('length', 'length'), ('GC', 'GC'), ('GC_position', 'GC + pos.')]
COLS = [(c, cat) for c in CX for cat in ('miRNA', 'RBP')]
lab_w = 62.0
cx0 = x+lab_w; cw = (x+w-1.0-cx0)/8; rh = 11.2; e0 = y+22.0
LIM = 12.0
for j, c in enumerate(CX):
    F.text(cx0+cw*(2*j+1), y+3.0, CL[c].replace('Human ', 'H ').replace('Mouse ', 'M '), 7.0, ha='center')
    F.line([cx0+cw*2*j+0.6, cx0+cw*(2*j+2)-0.6], [y+11.0]*2, lw='thin', color=INK)
    for k, cat in enumerate(('mi', 'R')):
        F.text(cx0+cw*(2*j+k+0.5), y+13.4, cat, 6.8, ha='center', color=MUTED)
Ev = {}; yy = e0
for cohort, src, glab in ROWS:
    for lev, ll in LEV:
        F.text(cx0-2.0, yy+rh/2-3.4, f'{glab}, {ll}' if lev == 'length' else ll, 7.0, ha='right', color=INK if lev == 'length' else MUTED)
        for j, (c, cat) in enumerate(COLS):
            r = bl[(bl.cohort == cohort) & (bl.source == src) & (bl.level == lev) & (bl.species == c[0]) & (bl.contrast == c[1]) & (bl.category == cat)]
            assert len(r) == 1
            r = r.iloc[0]; fc = diverge(float(r.mean_difference), LIM)
            tested = r.status == 'paired_gene_signs'
            F.rect(cx0+cw*j+0.4, yy+0.4, cw-0.8, rh-0.8, fc=fc, ec=INK if tested else 'none', lw='thin' if tested else 0)
            F.text(cx0+cw*(j+0.5), yy+rh/2-3.2, '0' if abs(r.mean_difference) < 0.5 else f'{r.mean_difference:.0f}'.replace('-', '−'), 6.8, ha='center', color=on(fc))
            Ev[f'{cohort}|{src}|{lev}|{CL[c]}|{cat}'] = dict(mean=float(r.mean_difference), n=int(r.n_genes), tested=bool(tested), q_all=float(r.q_all_backgrounds))
        yy += rh
    yy += 2.0
F.text(x+w/2, yy+1.0, 'outlined: tested (≥ 10 genes); all q = 1', 7.0, ha='center', color=MUTED)
F.cbar_key(x+w/2, y+h-7.4, 'sites/kb', '−12', '+12', lambda f_: diverge(-LIM+2*LIM*f_, LIM), wbar=40)
assert (bl.q_all_backgrounds >= 0.999).all()
F.values['e'] = Ev
F.notes.append('e: mi = miRNA, R = RBP; Same gene rows of the six-background cohort coincide with the fixed cohort (same windows); colour clipped at ±12 sites/kb')
F.save(OUT, D)
