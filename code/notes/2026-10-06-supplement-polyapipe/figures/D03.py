"""Supplementary Figure D03 (printed number from supdraw.SNUM) — frozen expression framework and design usability
(expression only; unchanged by the PAS reanalysis).

a, b  frozen joint expression embedding (display subsample), human and mouse, coloured by state family
c  routes usable for each APA design: source-qualified P–S, same-library P–M
d  human -> mouse state correspondence: row-normalised weight of the 62 x 49 actual state-pair edges, ordered by family
e  native programme scores (qualified library x state, >= 20 cells) per species: scale and composition differ
f  held-out expression correlation excess for the 12 human x mouse route pairs: original fit, its null range, and the
   full refit variants
"""
import sys, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D03'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 182.0, 262.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
PAL = json.loads(D.path('D03A', 'frozen/palette.json').read_text())['state_groups']
ST = D.df('D03D', 'frozen/states.tsv')
CODE = dict(shared_s='S', shared_r='R', shared_e='E', shared_tp='P', shared_m='M', shared_k='K', human_n='hN', human_x='hX',
            human_u='hU', mouse_b='mB', mouse_n='mN', mouse_i='mI', mouse_y='mY', mouse_u='mU')
FORD = list(CODE)

# ================================================================ a, b  embeddings
dc = D.df('D03A', 'frozen/display_cells.tsv.gz')
WE = 150.0
Av = {}
for j, (L, sp) in enumerate((('a', 'human'), ('b', 'mouse'))):
    x, y, w, h = F.panel(L, X0+j*(WE+GAP), Y1, WE, H1)
    q = dc[dc.species_key == sp].sample(frac=1, random_state=27)
    side = min(w-2.0, h-12.0-24.0)
    ex = Lin(-4.8, 14.6, x+(w-side)/2, x+(w-side)/2+side); ey = Lin(-3.8, 15.6, y+12.0+side, y+12.0)
    F.text(x+16.0, y+0.4, sp.capitalize(), 8.0, font='Bol')
    F.text(x+w-1.0, y+1.0, f'{len(q):,} cells shown', 7.0, ha='right', color=MUTED)
    for fam in FORD:
        z = q[q.reviewed_group == fam]
        if len(z): F.points([ex.clip(v) for v in z.UMAP1], [ey.clip(v) for v in z.UMAP2], ms='point', color=PAL[fam], alpha='secondary', rasterized=True)
    fams = [f for f in FORD if (q.reviewed_group == f).any()]
    half = (len(fams)+1)//2
    F.key_row(None, y+h-17.4, [('marker', ('o', PAL[f]), CODE[f]) for f in fams[:half]], size=7.0, gap=5.0, center=x+w/2)
    F.key_row(None, y+h-7.4, [('marker', ('o', PAL[f]), CODE[f]) for f in fams[half:]], size=7.0, gap=5.0, center=x+w/2)
    Av[sp] = {CODE[f]: int((q.reviewed_group == f).sum()) for f in fams}
F.values['ab'] = Av
F.notes.append('a, b: frozen UMAP coordinates of the stored 60,000-cell display subsample; no embedding or annotation refitted; draw order shuffled (seed 27)')

# ================================================================ c  design usability
x, y, w, h = F.panel('c', X0+2*(WE+GAP), Y1, RR-(X0+2*(WE+GAP)), H1)
pc = D.df('D03C', 'global/planned_route_contrasts.tsv')
ROUTES = ['GSE118258', 'GSE150311', '178325_0618', '178325_0809', 'GSE221739', 'CIP', 'WOT']
c1, c2 = x+w-40.0, x+w-14.0
F.text(c1, y+3.0, 'Source', 7.0, ha='center'); F.text(c1, y+10.0, 'P–S', 7.0, ha='center')
F.text(c2, y+3.0, 'Paired', 7.0, ha='center'); F.text(c2, y+10.0, 'P–M', 7.0, ha='center')
Cv = {}
for i, r in enumerate(ROUTES):
    yy = y+26.0+i*21.0
    F.rect(x+0.4, yy-2.8, 5.6, 5.6, fc=STUDY[r][1]); F.text(x+8.4, yy-3.4, STUDY[r][0], 7.0)
    q = pc[pc.route == r]
    a_ = bool(((q.contrast == 'P_S') & (q.design == 'true_source')).any()); b_ = bool(((q.contrast == 'P_M') & (q.design == 'same_library')).any())
    for cc, ok in ((c1, a_), (c2, b_)):
        F.marker(cc, yy, 'o', 'mark', fc=INK if ok else 'white', ec=INK if ok else AXC)
    Cv[r] = dict(source_P_S=a_, paired_P_M=b_)
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'used'), ('marker', ('o', 'white', AXC), 'not used')], gap=6.0, center=x+w/2)
F.values['c'] = Cv

# ================================================================ d  state correspondence
WD = 268.0
x, y, w, h = F.panel('d', X0, Y2, WD, H2)
ed = D.df('D03D', 'frozen/actual_state_pair_edges.tsv.gz')
hs = ST[ST.species == 'human'].copy(); ms = ST[ST.species == 'mouse'].copy()
for t_ in (hs, ms):
    t_['fo'] = t_.family.map({f: i for i, f in enumerate(FORD)})
hs = hs.sort_values(['fo', 'figure_row']); ms = ms.sort_values(['fo', 'figure_row'])
M_ = ed.pivot(index='human_state', columns='mouse_state', values='human_row_share').reindex(index=hs.state_key, columns=ms.state_key)
ml = 30.0
mx0 = x+ml; mw = x+w-1.0-mx0
cw = mw/len(ms); rh = (h-28.0-17.0)/len(hs)
my0 = y+28.0
for i, hk in enumerate(hs.state_key):
    for k, mk in enumerate(ms.state_key):
        v = M_.loc[hk, mk]
        fc = '#F3F5F6' if not (v == v) else seq(float(v))
        F.rect(mx0+k*cw, my0+i*rh, cw+0.05, rh+0.05, fc=fc)
# family strips and labels
for frame, horiz in ((hs, False), (ms, True)):
    codes = frame.family.tolist()
    pos = 0
    for f in dict.fromkeys(codes):
        n = codes.count(f)
        if horiz:
            F.rect(mx0+pos*cw, my0-6.0, n*cw, 4.0, fc=PAL[f])
            if n*cw > 6.0: F.text(mx0+(pos+n/2)*cw, my0-15.0, CODE[f][-1] if CODE[f] in ('S', 'R', 'E', 'P', 'M', 'K') else CODE[f], 6.8, ha='center', color=MUTED)
            if pos: F.line([mx0+pos*cw]*2, [my0, my0+len(hs)*rh], lw='ref', color=AXC)
        else:
            F.rect(mx0-6.0, my0+pos*rh, 4.0, n*rh, fc=PAL[f])
            if n*rh > 6.0: F.text(mx0-8.0, my0+(pos+n/2)*rh-3.2, CODE[f], 6.8, ha='right', color=MUTED)
            if pos: F.line([mx0, mx0+len(ms)*cw], [my0+pos*rh]*2, lw='ref', color=AXC)
        pos += n
F.text(mx0+mw/2, y+3.0, 'Mouse display families', 7.5, ha='center')
F.vtext(x+3.6, my0+len(hs)*rh/2, 'Human display families', 7.5)
F.cbar_key(mx0+mw/2, y+h-7.4, 'Row share of edge weight', '0', '1', lambda f_: seq(f_), wbar=60)
F.values['d'] = dict(human_states=hs.state_key.tolist(), mouse_states=ms.state_key.tolist(),
                     max_share=float(np.nanmax(M_.to_numpy())), edges=int(len(ed)), empty_rows=int((M_.fillna(0).sum(axis=1) == 0).sum()))
F.notes.append('d: rows/columns ordered by frozen family then figure row; light grey = no estimated edge (not zero weight); no lineage identity is inferred from a link')

# ================================================================ e  programme scales
x0e = X0+WD+GAP; we = RR-x0e
HE = 88.0
x, y, w, h = F.panel('e', x0e, Y2, we, HE)
pg = D.df('D03E', 'frozen/program_by_state_library.tsv.gz')
pg = pg[pg.n_cells >= 20]
sx = Lin(0, 3.5, x+62.0, x+w-4.0)
F.text((sx.p0+sx.p1)/2, y+3.0, 'Native program score', 7.5, ha='center')
Ev = {}; yy = y+20.0
for sp in ('human', 'mouse'):
    for k, prog in enumerate(('ECM', 'Epithelial', 'Pluripotency')):
        v = pg[(pg.species == sp) & (pg.program == prog)].mean_log1p_cp10k.to_numpy()
        if k == 0: F.text(x+0.4, yy-3.4, sp.capitalize(), 7.0)
        F.text(x+60.0, yy-3.4, prog, 7.0, ha='right', color=MUTED)
        q10, q25, q50, q75, q90 = np.percentile(v, [10, 25, 50, 75, 90])
        F.line([sx.clip(q10), sx.clip(q90)], [yy, yy], lw='thin', color=AXC)
        F.line([sx(q25), sx(q75)], [yy, yy], lw='range', color=NULLBAR)
        F.line([sx(q50)]*2, [yy-3.0, yy+3.0], lw='bold', color=INK)
        Ev[f'{sp}|{prog}'] = dict(n=int(len(v)), q10=float(q10), q25=float(q25), median=float(q50), q75=float(q75), q90=float(q90))
        yy += 9.0
    yy += 3.0
F.xaxis(sx, yy-4.0, [0, 1, 2, 3])
F.values['e'] = Ev

# ================================================================ f  held-out expression excess
x, y, w, h = F.panel('f', x0e, Y2+HE+GAP, we, H2-HE-GAP)
re_ = D.df('D03F', 'frozen/S2C_route_expression.tsv')
HU = {'GSE118258': 'GSE118258', 'GSE178325::donor_0605': '178325_0605', 'GSE178325::donor_0618': '178325_0618', 'GSE221739': 'GSE221739'}
MO = {'gse103221_cip': 'CIP', 'gse103221_osk': 'OSK', 'wot_common_to_2i': 'WOT'}
orig = re_[re_.setting == 'Original']
# row labels (2026-10-07): text instead of colour swatches; human route on the first row of each block, mouse route per row
HLAB = {'GSE118258': 'Xing2020', 'GSE178325::donor_0605': 'Guan 0605', 'GSE178325::donor_0618': 'Guan 0618', 'GSE221739': 'Panariello'}
MLAB = {'gse103221_cip': 'CIP', 'gse103221_osk': 'OSK', 'wot_common_to_2i': '2i'}
LX = x+0.4+max(textwidth(v, 6.8) for v in HLAB.values())+3.0
fx = Lin(-0.005, 0.055, LX+max(textwidth(v, 6.8) for v in MLAB.values())+5.0, x+w-4.0)
F.text((fx.p0+fx.p1)/2, y+3.0, 'Held-out correlation excess', 7.5, ha='center')
Fv = {}; rp = 9.0; y0 = y+21.0
for i, r in enumerate(orig.itertuples()):
    yy = y0+i*rp
    if i == 0 or r.human != orig.human.iloc[i-1]:
        F.text(x+0.4, yy-3.4, HLAB[r.human], 6.8)
    F.text(LX, yy-3.4, MLAB[r.mouse], 6.8, color=MUTED)
    F.line([fx.clip(r.null_q025-r.null_mean), fx.clip(r.null_q975-r.null_mean)], [yy, yy], lw='range', color=NULLBAR)
    z = re_[(re_.human == r.human) & (re_.mouse == r.mouse) & (re_.setting == 'Full variants') & (re_.status == 'OK')]
    if len(z):
        F.points([fx.clip(v) for v in z.excess], yy+np.linspace(-2.4, 2.4, len(z)), ms='point', color=PROX, alpha='full')
    F.marker(fx(r.excess), yy, 'D', 'small', fc=INK)
    Fv[f'{r.human}|{r.mouse}'] = dict(original=float(r.excess), null=[float(r.null_q025-r.null_mean), float(r.null_q975-r.null_mean)],
                                      refits_ok=int(len(z)), refits=[float(v) for v in z.excess])
F.line([fx(0)]*2, [y0-5.0, y0+11*rp+5.0], lw='ref', color=AXC, ls='ref')
F.xaxis(fx, y0+11*rp+6.0, [0, 0.02, 0.04], ['0', '.02', '.04'])
F.key_row(None, y+h-17.4, [('marker', ('D', INK), 'original'), ('marker', ('o', PROX), 'refit variants')], gap=6.0, center=x+w/2)
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'range'), 'original null 95%, centered')], gap=6.0, center=x+w/2)
F.values['f'] = Fv
F.notes.append('f: up to 12 full refit variants per pair; non-evaluable refits retained in the source table; the null range belongs to the original fit only; Guo2019–OSK included')
F.save(OUT, D)
