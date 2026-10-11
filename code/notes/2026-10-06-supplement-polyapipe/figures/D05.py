"""Supplementary Figure D05 (printed number from supdraw.SNUM) — route means, tandem-UTR restriction and state-reference
sensitivity (R1), polyApipe v2 primary.

a  P–S route mean ΔDWUI over 16 settings (4 PAS settings x 4 gene backgrounds); supported -> expanded change split into
   re-estimation on common genes and genes admitted only by the expanded layer
b  gene-level P–S ΔDWUI distribution per route (primary setting) and the route mean range over documented-source choices
c  mean P–S ΔDWUI for all eligible genes and for genes whose primary PAS pair lies in one annotated terminal 3' UTR
d  hN − P in the single RSeT library, same genes, four PAS settings (v2 only)
e  the three edges of the same RSeT library (hN−S, P−S, hN−P) at the primary setting
f  R/E/P/M minus observed S per route: available genes and the route-specific fixed common background
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D05'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2, H3, H4 = 110.0, 100.0, 86.0, 150.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; Y4 = Y3+H3+GAP; H = Y4+H4+MARGIN
F = Fig(H, figname(FID))
ROUTES = ['GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT']
PRIMARY = ('supported50_p08', 'available')
SETTINGS = ['supported50_p08', 'supported50_p09', 'supported50_p08_mc5', 'expanded_primary_p08']
SLAB = {'supported50_p08': 'Primary setting', 'supported50_p09': 'Purity ≥ .9', 'supported50_p08_mc5': '≥ 5 metacells',
        'expanded_primary_p08': 'Expanded layer'}    # primary = reference-supported layer, state purity >= 0.8, >= 3 metacells
COMMON, ADDED = ORD3[0], ORD3[2]


def route_rows(y0, pitch):
    return {r: y0+pitch*i for i, r in enumerate(ROUTES)}


def route_markers(xm, rows):
    for r, yy in rows.items():
        F.marker(xm, yy, STUDY[r][2], 'mark', fc=STUDY[r][1])


def quant_row(sx, yy, v, col, half=3.2):
    """gene distribution summary: q10–q90 thin, IQR bold, median tick, mean marker (display summary of the gene table)"""
    q10, q25, q50, q75, q90 = np.percentile(v, [10, 25, 50, 75, 90])
    F.line([sx.clip(q10), sx.clip(q90)], [yy, yy], lw='thin', color=AXC)
    F.line([sx.clip(q25), sx.clip(q75)], [yy, yy], lw='range', color=NULLBAR)
    F.line([sx.clip(q50)]*2, [yy-half, yy+half], lw='bold', color=INK)
    F.marker(sx.clip(float(np.mean(v))), yy, 'o', 'small', fc=col)
    return dict(n=int(len(v)), mean=float(np.mean(v)), q10=float(q10), q25=float(q25), median=float(q50), q75=float(q75), q90=float(q90))


# ================================================================ a  16 settings and layer decomposition
x, y, w, h = F.panel('a', X0, Y1, RR-X0, H1)
g = D.df('D05A', 'global/global_direction_summary.tsv')
g = g[(g.contrast == 'P_S') & (g.design == 'true_source')]
dec = D.df('D05A', 'global/pas_layer_coverage_decomposition.tsv')
dec = dec[(dec.contrast == 'P_S') & (dec.design == 'true_source')].set_index(['route', 'subset'])
rows = route_rows(y+21.0, 12.0)
xa = Lin(-0.04, 0.02, x+16, x+268)
xd = Lin(-0.009, 0.006, x+300, x+w-0.2)
F.text((xa.p0+xa.p1)/2, y+3.0, 'P–S route mean ΔDWUI, 16 settings', 7.5, ha='center')
F.text((xd.p0+xd.p1)/2, y+3.0, 'Supported → expanded change', 7.5, ha='center')
route_markers(x+5.0, rows)
for s_, v in ((xa, 0.0), (xd, 0.0)):
    F.line([s_(v)]*2, [rows[ROUTES[0]]-6, rows[ROUTES[-1]]+6], lw='ref', color=AXC, ls='ref')
Av = {}
for r, yy in rows.items():
    q = g[g.route == r]; assert len(q) == 16 and (q.status == 'measured').all()
    F.line([xa(q['mean'].min()), xa(q['mean'].max())], [yy, yy], lw='range', color=NULLBAR)
    oth = q[~((q.setting == PRIMARY[0]) & (q.background == PRIMARY[1]))]
    F.points([xa(v) for v in oth['mean']], [yy]*len(oth), ms='small', color=AXC, alpha='secondary')   # grey as in the key; study = row marker
    pm = float(q[(q.setting == PRIMARY[0]) & (q.background == PRIMARY[1])]['mean'].iloc[0])
    F.marker(xa(pm), yy, 'o', 'mark', fc='white', ec=DIST if pm > 0 else PROX)        # primary setting as in main Figure 1d
    # decomposition (display arithmetic on the recorded subset means and n)
    sc, ec, eo = dec.loc[(r, 'supported_common')], dec.loc[(r, 'expanded_common')], dec.loc[(r, 'expanded_only')]
    N = int(ec.n+eo.n)
    common = ec.n/N*(ec['mean']-sc['mean']); added = eo.n/N*(eo['mean']-sc['mean'])
    tot = q[(q.setting == 'expanded_primary_p08') & (q.background == 'available')]['mean'].iloc[0]-pm
    assert abs(common+added-tot) < 1e-9, (r, common+added, tot)
    for off, val, col in ((-2.4, common, COMMON), (2.4, added, ADDED)):
        F.rect(min(xd(0), xd(val)), yy+off-2.1, abs(xd(val)-xd(0)), 4.2, fc=col)
    Av[r] = dict(settings={f'{a}|{b}': float(m) for a, b, m in zip(q.setting, q.background, q['mean'])}, primary=pm,
                 range=[float(q['mean'].min()), float(q['mean'].max())], change_supported_to_expanded=float(tot),
                 common_reestimation=float(common), added_genes=float(added), n_common=int(ec.n), n_added=int(eo.n))
ay = rows[ROUTES[-1]]+7.0
F.xaxis(xa, ay, [-0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02], ['−.04', '−.03', '−.02', '−.01', '0', '.01', '.02'])
F.xaxis(xd, ay, [-0.008, -0.004, 0, 0.004], ['−.008', '−.004', '0', '.004'])
F.key_row(None, y+h-17.4, [('marker', (STUDY[r][2], STUDY[r][1]), STUDY[r][0]) for r in ROUTES], gap=7.0, center=x+w/2)
F.key_row(None, y+h-7.4, [('marker', ('o', 'white', INK), 'Primary setting'), ('marker', ('o', AXC), 'Other 15 settings'),
                          ('line', (NULLBAR, None, 'range'), 'Range')], gap=7.0, center=(xa.p0+xa.p1)/2)
F.key_row(None, y+h-7.4, [('block', COMMON, 'Common genes'), ('block', ADDED, 'Added genes')], gap=7.0, center=(xd.p0+xd.p1)/2)
F.values['a'] = Av
F.notes.append('a right: change = n_common/N (expanded_common − supported_common) + n_added/N (expanded_only − supported_common), '
               'N = n_common + n_added; display arithmetic, checked to equal the recorded expanded − supported route mean')

# ================================================================ b  gene distributions and source choice
WB = 266.0
x, y, w, h = F.panel('b', X0, Y2, WB, H2)
re_ = D.df('D05B', 'state_vectors/route_effects.primary.tsv.gz', usecols=['route', 'contrast', 'design', 'setting', 'delta_DWUI'])
re_ = re_[(re_.contrast == 'P_S') & (re_.design == 'true_source') & (re_.setting == 'supported50_p08')]
scr = D.df('D05B', 'panel_statistics/source-choice-ranges.tsv').set_index('route')
rows = route_rows(y+21.0, 12.0)
xb = Lin(-0.15, 0.10, x+16, x+152)
xs = Lin(-0.045, 0.025, x+160, x+w-8)
F.text((xb.p0+xb.p1)/2, y+3.0, 'Gene P–S ΔDWUI', 7.5, ha='center')
F.text((xs.p0+xs.p1)/2, y+3.0, 'Route mean by source', 7.5, ha='center')
F.text(x+w, y+3.0, 'n', 7.0, ha='right', color=MUTED)
route_markers(x+5.0, rows)
for s_ in (xb, xs):
    F.line([s_(0)]*2, [rows[ROUTES[0]]-6, rows[ROUTES[-1]]+6], lw='ref', color=AXC, ls='ref')
Bv = {}
for r, yy in rows.items():
    v = re_[re_.route == r].delta_DWUI.dropna().to_numpy()
    Bv[r] = quant_row(xb, yy, v, STUDY[r][1])
    c = scr.loc[r]
    if c.maximum > c.minimum:
        F.line([xs(c.minimum), xs(c.maximum)], [yy, yy], lw='range', color=NULLBAR)
    F.marker(xs(c['median']), yy, STUDY[r][2], 'small', fc=STUDY[r][1])
    F.text(x+w, yy-3.4, f'{int(c.n_settings)}', 7.0, ha='right', color=MUTED)
    Bv[r]['source_choice'] = dict(median=float(c['median']), min=float(c.minimum), max=float(c.maximum), n_settings=int(c.n_settings))
ay = rows[ROUTES[-1]]+7.0
F.xaxis(xb, ay, [-0.15, -0.1, -0.05, 0, 0.05, 0.1], ['−.15', '−.1', '−.05', '0', '.05', '.1'])
F.xaxis(xs, ay, [-0.04, -0.02, 0, 0.02], ['−.04', '−.02', '0', '.02'])
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'range'), 'IQR'), ('line', (AXC, None, 'thin'), '10–90%'),
                          ('marker', ('o', AXC), 'Mean')], gap=6.0, center=(xb.p0+xb.p1)/2)
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'range'), 'Min–max')], gap=6.0, center=(xs.p0+xs.p1)/2)
F.values['b'] = Bv

# ================================================================ c  tandem-UTR restriction
x, y, w, h = F.panel('c', X0+WB+GAP, Y2, RR-(X0+WB+GAP), H2)
tc = D.df('D05C', 'panel_statistics/tandem-comparison.tsv').set_index(['route', 'subset'])
rows = route_rows(y+21.0, 12.0)
xc = Lin(-0.06, 0.03, x+16, x+w-46)
F.text((xc.p0+xc.p1)/2, y+3.0, 'Mean P–S ΔDWUI', 7.5, ha='center')
F.text(x+w, y+3.0, 'All / tandem', 7.0, ha='right', color=MUTED)
route_markers(x+5.0, rows)
F.line([xc(0)]*2, [rows[ROUTES[0]]-6, rows[ROUTES[-1]]+6], lw='ref', color=AXC, ls='ref')
Cv = {}
for r, yy in rows.items():
    a_, t_ = tc.loc[(r, 'All genes')], tc.loc[(r, 'Tandem-compatible')]
    F.line([xc(a_['mean']), xc(t_['mean'])], [yy, yy], lw='thin', color=STUDY[r][1])
    F.marker(xc(a_['mean']), yy, 'o', 'mark', fc=STUDY[r][1])
    F.marker(xc(t_['mean']), yy, 'o', 'mark', fc='white', ec=STUDY[r][1])
    F.text(x+w, yy-3.4, f'{int(a_.n):,} / {int(t_.n):,}'.replace(',', ','), 7.0, ha='right')
    Cv[r] = dict(all=dict(mean=float(a_['mean']), n=int(a_.n)), tandem=dict(mean=float(t_['mean']), n=int(t_.n)))
ay = rows[ROUTES[-1]]+7.0
F.xaxis(xc, ay, [-0.06, -0.04, -0.02, 0, 0.02], ['−.06', '−.04', '−.02', '0', '.02'])
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'All genes'), ('marker', ('o', 'white', INK), 'Tandem 3′ UTR')],
          gap=7.0, center=x+w/2)
F.values['c'] = Cv

# ================================================================ d, e  RSeT same-library contrasts (v2 only)
t = D.df('D05D', 'naive_state/RSeT_same_library_triangle_gene_effects.tsv.gz')
t = t[(t.metric == 'DWUI') & (t.version == 'v2_primary')]
WD = 226.0
for L, x0, w0 in (('d', X0, WD), ('e', X0+WD+GAP, RR-(X0+WD+GAP))):
    x, y, w, h = F.panel(L, x0, Y3, w0, H3)
    sx = Lin(-0.06, 0.04, x+(66 if L == 'd' else 36), x+w-20)
    F.text(x+w, y+3.0, 'n', 7.0, ha='right', color=MUTED)
    if L == 'd':
        F.text((sx.p0+sx.p1)/2, y+3.0, 'hN − P ΔDWUI, one RSeT library', 7.5, ha='center')
        items = [(SLAB[s], t[t.setting == s].hN_P.dropna().to_numpy()) for s in SETTINGS]
    else:
        F.text((sx.p0+sx.p1)/2, y+3.0, 'Same library, primary setting', 7.5, ha='center')
        p = t[t.setting == 'supported50_p08']
        items = [('hN − S', p.hN_S.dropna().to_numpy()), ('P − S', p.P_S.dropna().to_numpy()), ('hN − P', p.hN_P.dropna().to_numpy())]
    pitch = 12.0 if L == 'd' else 18.0
    r0 = y+20.0
    F.line([sx(0)]*2, [r0-6, r0+pitch*(len(items)-1)+6], lw='ref', color=AXC, ls='ref')
    vals = {}
    for i, (lab, v) in enumerate(items):
        yy = r0+pitch*i
        F.text(x+0.4 if i or True else x, yy-3.4, lab, 7.0)
        vals[lab] = quant_row(sx, yy, v, INK)
        F.text(x+w, yy-3.4, str(len(v)), 7.0, ha='right', color=MUTED)
    ay = y+20.0+pitch*3+7.0 if L == 'd' else r0+pitch*(len(items)-1)+7.0
    ay = y+63.0
    F.xaxis(sx, ay, [-0.06, -0.04, -0.02, 0, 0.02, 0.04], ['−.06', '−.04', '−.02', '0', '.02', '.04'])
    F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'range'), 'IQR'), ('line', (AXC, None, 'thin'), '10–90%'),
                              ('marker', ('o', INK), 'Mean')], gap=6.0, center=(sx.p0+sx.p1)/2)
    F.values[L] = vals
F.notes.append('d, e: v2 annotation only (contract scope: remove cross-version matching and the v3 line); hN − S = hN − P + P − S per gene')

# ================================================================ f  state reference: available and fixed common genes (cells of main Figure 1d)
x, y, w, h = F.panel('f', X0, Y4, RR-X0, H4)
sr = D.df('D05F', 'geometry/state-reference-summary.tsv')
FR = ['GSE118258', 'GSE150311', '178325_0330', '178325_0605', '178325_0618', '178325_0809', 'GSE221739',
      'CIP', 'WOT', 'OSK', '171127_GETM', '171127_OSKM']
STS = ['R', 'E', 'P', 'M']
BGS = [('available', 'Available genes'), ('fixed_common_observed_states_this_route', 'Fixed common genes per route')]
labw = max(textwidth(STUDY[r][0], 7.0) for r in FR)
xl = x+13.0
NW = 24.0                                            # gene-count column of the fixed common background
hx0 = xl+7.6+labw+3.0
hw = (x+w-hx0-14.0-NW)/2
cw = hw/4; rh = 9.2; r0 = y+24.0
F.vtext(x+4.6, r0+rh*3.5, 'Human', 7.5); F.vtext(x+4.6, r0+rh*9.5+3.0, 'Mouse', 7.5)
LIM = 0.03; K_ = (cw/2-1.5)/LIM
Fv = {}
for b_, (bg, title) in enumerate(BGS):
    hx = hx0+b_*(hw+14.0)
    F.text(hx+hw/2, y+3.0, title, 7.5, ha='center')
    for k, s_ in enumerate(STS):
        F.text(hx+cw*(k+0.5), y+12.6, s_, 7.5, ha='center', color=SC[s_])
        F.line([hx+cw*k+1.5, hx+cw*(k+1)-1.5], [y+21.0]*2, lw='bold', color=SC[s_])
    for i, r in enumerate(FR):
        yy = r0+rh*i+(3.0 if i >= 7 else 0)
        if b_ == 0:
            F.rect(xl, yy+rh/2-2.8, 5.6, 5.6, fc=STUDY[r][1]); F.text(xl+7.6, yy+rh/2-3.4, STUDY[r][0], 7.0)
        nfix = None
        for k, s_ in enumerate(STS):
            q = sr[(sr.route == r) & (sr.background == bg) & (sr.state == s_)]
            assert len(q) == 1
            q = q.iloc[0]
            cxm, cym = hx+cw*(k+0.5), yy+rh/2
            if q.status == 'measured':
                v = float(q['mean']); vc = max(-LIM, min(LIM, v))
                F.rect(hx+cw*k+0.6, cym-3.6, cw-1.2, 7.2, fc=BAND)
                F.rect(min(cxm, cxm+vc*K_), cym-2.7, abs(vc)*K_, 5.4, fc=DIST if v > 0 else PROX)
                F.line([cxm]*2, [cym-4.0, cym+4.0], lw='ref', color=AXC)
                Fv[f'{bg}|{r}|{s_}'] = dict(mean=v, n=int(q.n)); nfix = int(q.n)
            else:
                F.text(cxm, cym-3.2, '–', 6.8, ha='center', color=AXC)
                Fv[f'{bg}|{r}|{s_}'] = None
        if b_ == 1 and nfix is not None:
            F.text(x+w-0.4, yy+rh/2-3.2, f'{nfix:,}', 6.8, ha='right', color=MUTED)
F.text(x+w-0.4, y+12.6, 'Genes', 7.0, ha='right', color=MUTED)
# key of main Figure 1d: 'Proximal −.03 [bar] +.03 Distal', bar half-length = ±.03 = the cell's bar scale
kx = hx0+hw+7.0; hb = LIM*K_; ky = y+h-7.4
kparts = [('Mean ΔDWUI vs observed S', 7.0, 8.0), ('Proximal', 7.0, 3.0), ('−.03', 7.0, 1.8), (None, 2*hb, 1.8), ('+.03', 7.0, 3.0), ('Distal', 7.0, 0.0)]
kw_ = sum((textwidth(t, sz) if t else sz)+g for t, sz, g in kparts)
kx0 = kx-kw_/2
for t, sz, g in kparts:
    if t is None:
        F.rect(kx0, ky+1.0, hb, 5.6, fc=PROX); F.rect(kx0+hb, ky+1.0, hb, 5.6, fc=DIST); kx0 += sz+g
    else:
        F.text(kx0, ky, t, sz, color=MUTED if t in ('−.03', '+.03') else INK); kx0 += textwidth(t, sz)+g
F.values['f'] = Fv
F.notes.append('f: cells drawn as in main Figure 1d (bar from the cell centre, proximal blue / distal orange, ±.03 = half cell, clipped); '
               '– = state not observed or no observed S in the route (unavailable, not zero); gene counts per cell in values, '
               'the fixed common background count printed at right')
F.save(OUT, D)
