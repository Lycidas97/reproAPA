"""Supplementary Figure D25 (printed number from supdraw.SNUM) — expression of the candidate RNA-binding proteins (R5),
polyApipe v2 primary.

a  state-resolved expression (mean log2(CPM + 1)) of the six frozen display RBPs per species in both studies; circle =
   below 2 CPM in every library, dash = state not sampled
b  expression filter of the third rug row of Figure 5f: every RBP named by a transferred motif in the two Trip12 strict
   intervals, P-state expression in Guo2019–CIP and Schiebinger2019–2i; an RBP is kept when it reaches 2 CPM in P of both
   studies (Rbm24 and Celf3 removed)
c, d  P–S and P–M expression change of every candidate RBP in the two studies of each species, and the direction counts
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D25'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 102.0, 124.0, 150.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname(FID))
ROUTES = {'human': ('GSE118258', 'GSE150311'), 'mouse': ('CIP', 'WOT')}
CS = {'SP': SC['S'], 'MP': SC['M']}
cand = D.df('D25A', 'rbp_expression/candidates.tsv')
sx_ = D.df('D25A', 'rbp_expression/state_expression.tsv')
LMAX = 11.0

# ================================================================ a  display RBPs
x, y, w, h = F.panel('a', X0, Y1, RR-X0, H1)
hw = (w-12.0)/2
Av = {}
for j, sp in enumerate(('human', 'mouse')):
    x0 = x+j*(hw+12.0)
    facs = sorted(cand.loc[(cand.species == sp) & cand.display_selected, 'factor'])
    assert len(facs) == 6
    lw_ = max(textwidth(f, 7.0, 'RegIta') for f in facs)+5.0
    cx0 = x0+lw_; cw = (x0+hw-cx0-6.0)/10; rh = 10.0; r0 = y+28.0
    F.text(x0+16.0 if j == 0 else x0, y+0.4, sp.capitalize(), 8.0, font='Bol')
    for k, route in enumerate(ROUTES[sp]):
        xm = cx0+(5*k+2.5)*cw+(6.0 if k else 0.0)
        F.rect(xm-textwidth(STUDY[route][0], 7.0)/2-7.6, y+12.6, 5.6, 5.6, fc=STUDY[route][1])
        F.text(xm, y+12.0, STUDY[route][0], 7.0, ha='center')
        for i, st_ in enumerate('SREPM'):
            F.text(cx0+(5*k+i+0.5)*cw+(6.0 if k else 0.0), y+20.6, st_, 6.8, ha='center', color=SC[st_])
    for i, f in enumerate(facs):
        yy = r0+i*rh
        F.text(cx0-2.0, yy+rh/2-3.4, f, 7.0, font='RegIta', ha='right')
        for k, route in enumerate(ROUTES[sp]):
            for q, st_ in enumerate('SREPM'):
                r = sx_[(sx_.species == sp) & (sx_.factor == f) & (sx_.route == route) & (sx_.state == st_)].iloc[0]
                cx = cx0+(5*k+q)*cw+(6.0 if k else 0.0)
                if r.status != 'available' or not np.isfinite(r.mean_log2CPM1):
                    F.rect(cx+0.4, yy+0.4, cw-0.8, rh-0.8, fc='white', ec=AXC, lw='axis')
                    F.text(cx+cw/2, yy+rh/2-3.2, '–', 6.8, ha='center', color=MUTED)
                    Av[f'{sp}|{f}|{route}|{st_}'] = None
                    continue
                fc = seq(r.mean_log2CPM1, 0, LMAX)
                F.rect(cx+0.4, yy+0.4, cw-0.8, rh-0.8, fc=fc)
                if not bool(r.expressed_any_2CPM):
                    F.marker(cx+cw/2, yy+rh/2, 'o', 'small', fc='white', ec=on(fc) if on(fc) == 'white' else INK)
                Av[f'{sp}|{f}|{route}|{st_}'] = dict(log2CPM1=float(r.mean_log2CPM1), expressed_2CPM=bool(r.expressed_any_2CPM))
F.cbar_key(x+w/2-60.0, y+h-7.4, 'Mean log2(CPM + 1)', '0', '11', lambda f_: seq(f_*LMAX, 0, LMAX), wbar=60)
F.key_row(x+w/2+70.0, y+h-7.4, [('marker', ('o', 'white', INK), '< 2 CPM in all libraries'), ('block', 'white', '– not sampled')], gap=6.0)
F.values['a'] = Av

# ================================================================ b  Trip12 transferred-motif filter
x, y, w, h = F.panel('b', X0, Y2, RR-X0, H2)
GID = 'ENSMUSG00000026219'
inv = D.df('D25B', 'site_statistics/site_factor_inventory.tsv.gz')
rr = D.df('D25B', 'candidate_screen/recurrent_region_pairs.tsv.gz')
ivs = sorted({(int(a), int(b)) for a, b in rr[(rr.gene_id == GID) & (rr.contrast == 'MP')][['shared_start0', 'shared_end0']].values})
t = inv[(inv.gene_id == GID) & (inv.provenance == 'RBP_transfer')]
t = t[[any(r.start0 >= s0 and r.end0 <= e0 for s0, e0 in ivs) for r in t.itertuples()]].drop_duplicates(['factor', 'start0', 'end0', 'names'])
pP = sx_[(sx_.species == 'mouse') & (sx_.state == 'P') & sx_.route.isin(['CIP', 'WOT'])]
okr = {f_: set(g_.loc[g_.expressed_any_2CPM.astype(bool), 'route']) for f_, g_ in pP.groupby('factor')}
P_RBP = {f_ for f_, rs in okr.items() if rs == {'CIP', 'WOT'}}
t = t.assign(keep=[any(n in P_RBP for n in s_.split(';')) for s_ in t.names])
names = sorted({n for s_ in t.names for n in s_.split(';')})
npos = {n: int(t[[n in s_.split(';') for s_ in t.names]][['start0', 'end0']].drop_duplicates().shape[0]) for n in names}
lv = {(f_, r_): float(v) for f_, r_, v in pP[['factor', 'route', 'mean_log2CPM1']].itertuples(index=False)}
kept = sorted([n for n in names if n in P_RBP], key=lambda n: -np.mean([lv[(n, 'CIP')], lv[(n, 'WOT')]]))
drop = sorted([n for n in names if n not in P_RBP])
assert drop == ['Celf3', 'Rbm24'], drop
assert int(t.keep.sum()) == 29 and t[t.keep][['start0', 'end0']].drop_duplicates().shape[0] == 28
order = kept+drop
by = Lin(0, 9, y+72.0, y+20.5)
x0b = x+22.0; cwb = (x+w-1.0-x0b-6.0)/len(order)
F.text(x0b+cwb*len(kept)/2, y+3.0, 'Kept', 7.5, ha='center')   # 2026-10-12: keep rule and counts moved to the caption
F.text(x0b+cwb*(len(kept)+1)+6.0, y+3.0, 'Removed', 7.5, ha='center', color=MUTED)
F.yaxis(by, x0b-2.0, [0, 3, 6, 9], None, 'log2(CPM + 1)')
F.line([x0b+cwb*len(kept)+3.0]*2, [by.p1-4.0, by.p0+4.0], lw='ref', color=AXC, ls='ref')
Bv = {}
for i, n in enumerate(order):
    cx = x0b+cwb*(i+0.5)+(6.0 if n in drop else 0.0)
    F.line([cx]*2, [by.p0, by.p1], lw='ref', color=BAND)
    for k, (route, dxo) in enumerate((('CIP', -2.2), ('WOT', 2.2))):
        v = lv[(n, route)]
        good = route in okr.get(n, set())
        F.marker(cx+dxo, by.clip(v), STUDY[route][2], 'small', fc=STUDY[route][1] if good else 'white', ec=STUDY[route][1])
    F.text(cx, by.p0+5.0, str(npos[n]), 6.8, ha='center', color=MUTED)
    F.vlabel(cx, by.p0+14.0+textwidth(n, 7.0, 'RegIta'), n, 7.0, font='RegIta', color=INK if n in kept else MUTED)
    Bv[n] = dict(CIP=lv[(n, 'CIP')], WOT=lv[(n, 'WOT')], CIP_2CPM='CIP' in okr.get(n, set()), WOT_2CPM='WOT' in okr.get(n, set()),
                 kept=n in kept, positions=npos[n])
F.text(x0b-4.0, by.p0+5.0, 'sites', 6.8, ha='right', color=MUTED)
F.key_row(None, y+h-7.4, [('marker', ('s', STUDY['CIP'][1]), 'Guo2019–CIP'), ('marker', ('o', STUDY['WOT'][1]), 'Schiebinger2019–2i'),
                          ('marker', ('o', 'white', INK), 'below 2 CPM in all P libraries')], gap=7.0, center=x+w/2)
F.values['b'] = dict(intervals=ivs, kept=kept, removed=drop, rbps=Bv)
F.notes.append('b: same table, interval and filter as Figure 5f (site_factor_inventory, provenance RBP_transfer, inside the shared strict P–M intervals); '
               'numbers under the columns are distinct site positions named by that RBP; 29 kept annotations at 28 positions')

# ================================================================ c, d  cross-study changes
cs = D.df('D25C', 'rbp_expression/cross_study_summary.tsv')
WC = (RR-X0-GAP)/2
DIRS = [('both_increased', 'both up'), ('both_decreased', 'both down'), ('opposite', 'opposite'), ('zero_or_tied', 'zero/tied')]
Cv = {}
for j, (L, sp) in enumerate((('c', 'human'), ('d', 'mouse'))):
    x, y, w, h = F.panel(L, X0+j*(WC+GAP), Y3, WC, H3)
    a_, b_ = ROUTES[sp]
    q = cs[(cs.species == sp) & cs.contrast.isin(['SP', 'MP'])]
    side = 104.0
    sx = Lin(-6, 7, x+28.0, x+28.0+side); sy = Lin(-6, 7, y+16.0+side, y+16.0)
    F.text((sx.p0+sx.p1)/2, y+3.0, f'{sp.capitalize()}, {q.factor.nunique()} RBPs', 7.5, ha='center')
    F.line([sx(0)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref'); F.line([sx.p0, sx.p1], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    for co, mk in (('SP', 'o'), ('MP', '^')):
        z = q[q.contrast == co]
        for r in z.itertuples():
            F.marker(sx.clip(r.effect_a), sy.clip(r.effect_b), mk, 'small', fc=CS[co] if r.expressed_in_both_studies else 'white', ec=CS[co])
    F.yaxis(sy, sx.p0, [-6, -3, 0, 3, 6], ['−6', '−3', '0', '3', '6'], f'Δ, {STUDY[b_][0]}')
    F.xaxis(sx, sy.p0+1.0, [-6, -3, 0, 3, 6], ['−6', '−3', '0', '3', '6'], f'Δ, {STUDY[a_][0]}')
    # direction counts
    bx0 = sx.p1+14.0; bxx = Lin(0, 0.66, bx0+40.0, x+w-10.0)
    F.text((bx0+x+w)/2, y+3.0, 'Direction', 7.5, ha='center')
    yy = y+24.0
    for co in ('SP', 'MP'):
        z = q[q.contrast == co]; n = len(z)
        F.text(bx0, yy-3.4, f'{co[1]}–{co[0]} ({n})', 7.0, color=CS[co] if False else INK)
        yy += 10.0
        for d_, lab in DIRS:
            k = int((z.direction == d_).sum())
            if d_ == 'zero_or_tied' and k == 0: continue
            F.text(bx0+38.0, yy-3.4, lab, 6.8, ha='right', color=MUTED)
            F.rect(bxx(0), yy-3.2, bxx(k/n)-bxx(0), 6.4, fc=CS[co])
            F.text(bxx(k/n)+1.2, yy-3.2, str(k), 6.8, color=MUTED)
            Cv[f'{sp}|{co}|{d_}'] = k
            yy += 9.0
        yy += 6.0
    F.key_row(None, y+h-7.4, [('marker', ('o', CS['SP']), 'P–S'), ('marker', ('^', CS['MP']), 'P–M'),
                              ('marker', ('o', 'white', INK), 'not expressed in both')], gap=6.0, center=x+w/2)
F.values['cd'] = Cv
F.notes.append('c, d: Δ = mean log2(CPM + 1) of P minus the reference state per study; every candidate RBP and contrast shown; '
               'expression change does not establish binding or activity')
F.save(OUT, D)
