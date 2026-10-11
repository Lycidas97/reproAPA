"""Supplementary Figure D22 (printed number from supdraw.SNUM) — partial retention: the same factor's exact predicted sites
change and stay (R5), polyApipe v2 primary.

a  fraction of changed multisite gene–factor relationships that are partial (some exact sites of the same factor changed,
   others retained), per study, contrast and site category, with 95% intervals
b  cross-study partial relationships and their genes per species and category (P–S and P–M pooled)
c, d  the two frozen partial-retention examples (P–M, one RBP factor each): mouse Snx30 / Carhsp1 sites and human SH3PXD2A /
   BOLL sites; left, every native PAS of the gene model on the genomic axis, filled when the 3' end at that PAS keeps the
   site; right, equal-library inclusion of the changed and the retained site by state (20 and 50 UMI), library values and
   empirical IQR
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D22'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 101.0, 106.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H2+MARGIN
F = Fig(H, figname(FID))
CAT = {'miRNA': INK, 'RBP': SEQ_STOPS[1]}
ROLE = {'changed': INK, 'retained': SEQ_STOPS[1]}
ROUTE_ORDER = [('human', 'GSE118258'), ('human', 'GSE150311'), ('mouse', 'CIP'), ('mouse', 'WOT')]

# ================================================================ a  partial fraction (form of main Figure 5d)
WA = 300.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
pr = D.df('D22A', 'site_statistics/partial_retention_summary.tsv')
pr = pr[pr.contrast.isin(['SP', 'MP'])]
Av = {}
for sp, route in ROUTE_ORDER:                     # values in the original order
    for co in ('SP', 'MP'):
        for cat in ('miRNA', 'RBP'):
            r = pr[(pr.route == route) & (pr.contrast == co) & (pr.category == cat)].iloc[0]
            Av[f'{route}|{co}|{cat}'] = dict(fraction=float(r.partial_fraction), ci=[float(r.ci_low), float(r.ci_high)],
                                             partial=int(r.partial_including_mixed), changed=int(r.changed_multisite_pairs), changed_genes=int(r.changed_genes))
ORD_A = [('CIP', 'SP'), ('CIP', 'MP'), ('WOT', 'SP'), ('WOT', 'MP'), ('GSE118258', 'SP'), ('GSE118258', 'MP'), ('GSE150311', 'SP'), ('GSE150311', 'MP')]
ry_, yy = {}, y+21.0
for k, (route, co) in enumerate(ORD_A):
    if k: yy += 6.6 if co == 'MP' else (11.0 if k == 4 else 8.4)
    ry_[(route, co)] = yy
LBA, GW = 34.0, 22.0
sw_ = (w-LBA-GW-12.0)/2
for route, co in ORD_A:
    yy = ry_[(route, co)]
    if co == 'SP': F.marker(x+3.0, yy+3.3, STUDY[route][2], 'mark', fc=STUDY[route][1])
    F.text(x+8.4, yy-3.4, f'{co[1]}–{co[0]}', 7.0)
    F.text(x+w-GW/2, yy-3.4, str(Av[f'{route}|{co}|miRNA']['changed_genes']), 6.8, ha='center', color=MUTED)
F.text(x+w-GW/2, y+3.0, 'Genes', 7.0, ha='center', color=MUTED)
for j_, cat in enumerate(('miRNA', 'RBP')):
    px = Lin(0, 0.6, x+LBA+j_*(sw_+12.0), x+LBA+j_*(sw_+12.0)+sw_)
    F.text((px.p0+px.p1)/2, y+3.0, f'{cat}, partial retention (%)', 7.5, ha='center')
    for v in (0, 0.3, 0.6): F.line([px(v)]*2, [y+16.0, ry_[ORD_A[-1]]+4.0], lw='ref', color=BAND)
    for route, co in ORD_A:
        a_ = Av[f'{route}|{co}|{cat}']; yy = ry_[(route, co)]; col = STUDY[route][1]
        F.whisker(px(a_['ci'][0]), px(a_['ci'][1]), yy, color=col, lw='thin', cap=1.2)
        if co == 'SP': F.marker(px(a_['fraction']), yy, 'o', 'small', fc=col)
        else: F.marker(px(a_['fraction']), yy, 'D', 'small', fc='white', ec=col)
    F.xaxis(px, ry_[ORD_A[-1]]+5.0, [0, 0.3, 0.6], ['0', '30', '60'])
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'P–S'), ('marker', ('D', 'white', INK), 'P–M'), ('line', (INK, None, 'thin'), '95% interval'),
                          ('marker', ('s', STUDY['CIP'][1]), 'CIP'), ('marker', ('o', STUDY['WOT'][1]), '2i'),
                          ('marker', ('o', STUDY['GSE118258'][1]), 'Xing2020'), ('marker', ('s', STUDY['GSE150311'][1]), 'Liu2020')],
          gap=5.0, size=7.0, center=x+w/2)
F.values['a'] = Av
lo, hi = min(v['fraction'] for v in Av.values()), max(v['fraction'] for v in Av.values())
F.notes.append(f'a: partial fraction range {lo:.3f}–{hi:.3f}; genes = genes with a changed multisite relationship (miRNA row) in that study')

# ================================================================ b  relationship counts
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
rc = D.df('D22B', 'supplement_statistics/partial-relationship-counts.tsv')
bx = Lin(0, 1200, x+30.0, x+w-4.0)
F.text((bx.p0+bx.p1)/2, y+3.0, 'Partial in both studies', 7.5, ha='center')
Bv = {}; yy = y+22.0
for sp in ('human', 'mouse'):
    for k, cat in enumerate(('miRNA', 'RBP')):
        r = rc[(rc.species == sp) & (rc.category == cat)].iloc[0]
        ry = yy+k*13.0
        if k == 0: F.text(x+0.4, ry-3.4, sp.capitalize(), 7.0)
        F.rect(bx(0), ry-4.5, bx(r.relationships)-bx(0), 9.0, fc=CAT[cat])
        F.line([bx(r.genes)]*2, [ry-4.5, ry+4.5], lw='bold', color='white' if cat == 'miRNA' else INK)
        F.text(bx(r.relationships)+1.4, ry-3.2, f'{int(r.relationships)}', 6.8, color=MUTED)
        Bv[f'{sp}|{cat}'] = dict(relationships=int(r.relationships), genes=int(r.genes), scope=r.contrast_scope)
    yy += 28.0
F.xaxis(bx, yy-12.0, [0, 500, 1000], ['0', '500', '1,000'], 'Relationships')
F.key_row(None, y+h-7.4, [('line', (INK, None, 'bold'), 'genes')], gap=5.0, center=x+w/2)
F.values['b'] = Bv
F.notes.append('b: relationship = same gene and factor partial in both studies of the species; P–S and P–M pooled; tick = number of genes')

# ================================================================ c, d  examples
sp_sites = D.df('D22C', 'supplement_statistics/partial-selected-site-pairs.tsv')
axs = D.df('D22D', 'supplement_statistics/partial-site-PAS-axis.tsv')
masks = D.df('D22C', 'supplement_statistics/partial-site-PAS-masks.tsv')
summ = D.df('D22E', 'supplement_statistics/partial-site-state-summary.tsv')
libs = D.df('D22E', 'supplement_statistics/partial-site-library-inclusion.tsv.gz')
EX = [('c', 'mouse', 'Snx30', ('CIP', 'WOT')), ('d', 'human', 'SH3PXD2A', ('GSE118258', 'GSE150311'))]
Gv = {}
for L, sp, sym, routes in EX:
    yP = Y2 if L == 'c' else Y3
    x, y, w, h = F.panel(L, X0, yP, RR-X0, H2)
    s_ = sp_sites[sp_sites.gene_symbol == sym]
    gid = s_.gene_id.iloc[0]; fac = s_.factor.iloc[0]; strand = s_.strand.iloc[0]; chrom = s_.chrom.iloc[0]
    assert s_.gene_id.nunique() == 1 and s_.factor.nunique() == 1
    ax_ = axs[axs.gene_id == gid]
    allc = list(ax_.cut_interbase)+list(s_.start0)+list(s_.end0)
    origin = min(allc) if strand == '+' else max(allc)
    tr = (lambda c: (c-origin)/1000) if strand == '+' else (lambda c: (origin-c)/1000)
    span = max(tr(c) for c in allc)
    WG = 196.0
    gx = Lin(-0.2, span+0.3, x+70.0, x+WG)
    F.text(x+16.0, y+0.4, sym, 8.0, font='BolIta')
    F.text(x+16.0+textwidth(sym, 8.0, 'BolIta')+3.0, y+1.0, f'{fac.upper()} sites, P–M, {chrom} ({strand.replace("-", "−")})', 7.0, color=MUTED)
    rows_y = {}
    Gl = {}
    for i, route in enumerate(routes):
        for k, role in enumerate(('changed', 'retained')):
            ry = y+24.0+i*26.0+k*10.0; rows_y[(route, role)] = ry
            if k == 0:
                F.marker(x+3.0, ry+5.0, STUDY[route][2], 'mark', fc=STUDY[route][1])
                F.text(x+8.4, ry+5.0-3.4, STUDY[route][0], 7.0)
            m = masks[(masks.gene_id == gid) & (masks.route == route) & (masks.site_role == role)].sort_values('pas_order')
            F.line([gx(tr(m.cut_interbase.min()) if strand == '+' else tr(m.cut_interbase.max())), gx(max(tr(c) for c in m.cut_interbase))],
                   [ry, ry], lw='thin', color=AXC)
            for c_, keep in zip(m.cut_interbase, m.retains_site):
                F.marker(gx(tr(c_)), ry, 'o', 'small', fc=ROLE[role] if keep else 'white', ec=ROLE[role])
            site = s_[(s_.route == route) & (s_.site_role == role)].iloc[0]
            a0, a1 = sorted((tr(site.start0), tr(site.end0)))
            F.rect(gx(a0)-1.0, ry-4.6, max(gx(a1)-gx(a0), 0.6)+2.0, 2.2, fc=ROLE[role])
            Gl[f'{route}|{role}'] = dict(site=[int(site.start0), int(site.end0)], pas=[int(c) for c in m.cut_interbase],
                                         retains=[int(v) for v in m.retains_site])
    yb = F.xaxis(gx, y+24.0+26.0+10.0+8.0, [t for t in range(0, int(span)+1)], [str(t) for t in range(0, int(span)+1)])
    F.text((gx.p0+gx.p1)/2, yb+1.6, '5′ → 3′ distance (kb)', 8.0, ha='center')
    # inclusion by state
    sub_w = (x+w-(x+WG+8.0)-8.0-1.0)/2
    for j, route in enumerate(routes):
        sx0 = x+WG+8.0+j*(sub_w+8.0)
        sx = Lin(-0.5, 4.8, sx0+(18.0 if j == 0 else 0.0), sx0+sub_w)
        sy = Lin(0, 1, y+80.0, y+16.0)
        F.marker(sx.p0+4.0, sy.p1-7.0, STUDY[route][2], 'mark', fc=STUDY[route][1])
        F.text(sx.p0+9.0, sy.p1-10.4, STUDY[route][0], 7.0)
        F.line([sx(3.75)]*2, [sy.p0, sy.p1], lw='ref', color=AXC, ls='ref')
        for role in ('changed', 'retained'):
            for cut, lsd in ((20, None), (50, 'dash')):
                z = summ[(summ.gene_id == gid) & (summ.route == route) & (summ.site_role == role) & (summ.cutoff == cut)].set_index('state').reindex(list('SREPM'))
                assert len(z) == 5
                F.line([sx(i) for i in range(4)], [sy(v) for v in z.equal_library.iloc[:4]], lw='data', color=ROLE[role], ls=lsd)
                if cut == 20:
                    for i, st_ in enumerate('SREPM'):
                        xp = i if st_ != 'M' else 4.5
                        v = libs[(libs.gene_id == gid) & (libs.route == route) & (libs.site_role == role) & (libs.cutoff == 20) & (libs.state == st_)].inclusion.to_numpy()
                        off = np.linspace(-0.18, 0.18, len(v)) if len(v) > 1 else np.zeros(len(v))
                        F.points([sx(xp+o) for o in off], [sy(t) for t in v], ms='point', color=SC[st_] if role == 'changed' else ROLE[role],
                                 alpha='secondary')   # changed-site libraries in state colour (main Figure 5e)
                        F.line([sx(xp)]*2, [sy(z.loc[st_, 'q25']), sy(z.loc[st_, 'q75'])], lw='range', color=ROLE[role], alpha='dense')
                for i, st_ in enumerate('SREPM'):
                    xp = i if st_ != 'M' else 4.5
                    F.marker(sx(xp), sy(z.loc[st_, 'equal_library']), 'o', 'small', fc=ROLE[role] if cut == 20 else 'white', ec=ROLE[role])
                Gl[f'{route}|{role}|{cut}'] = {st_: float(z.loc[st_, 'equal_library']) for st_ in 'SREPM'}
        F.yaxis(sy, sx.p0, [0, 0.5, 1], ['0', '50', '100'], 'Interval inclusion (%)' if j == 0 else None, labels_off=(j == 1))
        F.xaxis(sx, sy.p0+1.0, [0, 1, 2, 3, 4.5], ['']*5)
        for v_, st_ in zip([0, 1, 2, 3, 4.5], 'SREPM'):                    # state-coloured state labels (main Figure 5e)
            F.text(sx(v_), sy.p0+1.0+TICK+1.0, st_, 7.5, ha='center', color=SC[st_])
    kx = x+WG+8.0+18.0+(x+w-(x+WG+8.0+18.0))/2
    F.key_row(None, y+h-7.4, [('line', (ROLE['changed'], None, 'data'), 'Changed site'), ('line', (ROLE['retained'], None, 'data'), 'Constitutive site'),
                              ('line', (INK, None, 'data'), '20 UMI'), ('line', (INK, 'dash', 'data'), '50 UMI')], gap=6.0, center=kx)
    F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'PAS end keeps site'), ('marker', ('o', 'white', INK), 'site lost')], gap=6.0,
              center=(gx.p0+gx.p1)/2-10.0)
    Gv[sym] = Gl
F.values['cd'] = Gv
F.notes.append('c, d: both frozen partial-retention examples of the statistics package; predicted sites do not establish binding; '
               'IQR is empirical over libraries, not a confidence interval; M is a separate reference state, not a time step')
F.save(OUT, D)
