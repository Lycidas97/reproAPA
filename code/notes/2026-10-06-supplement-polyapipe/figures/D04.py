"""Supplementary Figure D04 (printed number from supdraw.SNUM) — PAS measurement support under the integrated v2 primary reference (R1.2), polyApipe.

A  per dataset: primary PAS split by reference support (±50 nt); PAS and UMI fractions within ±25/50/100 nt windows
B  canonical poly(A)-signal start frequency relative to the PAS, supported vs unsupported primary PAS vs matched background
C  read-end positions around supported and unsupported primary PAS (mean per-site fraction)
D  genes with complete S–R–E–P profiles in the three fitting routes at 20/40/80 PAS-UMI depth, two PAS layers
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

D = Data(); OUT = outdir('D04')
X0, R = MARGIN, PAGE_W-MARGIN
H1, H2, H3 = 115.0, 99.0, 101.5   # 2026-10-12: titles -> axis titles (a), title rows removed (c, d)
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP; H = Y3+H3+MARGIN
F = Fig(H, figname('D04'))
SUPC, UNSC = ORD3[0], ORD3[2]
DS = DATASET_ORDER

# ================================================================ A  catalogue and reference support
x, y, w, h = F.panel('A', X0, Y1, R-X0, H1)
sup = D.df('D04A', 'panel_statistics/primary-PAS-and-UMI-support.tsv').set_index('study_id')
c1 = Lin(0, 120, x+64, x+176); c2 = Lin(0.6, 1.0, x+192, x+308); c3 = Lin(0.88, 1.0, x+324, x+w-3)
rows = []; yy = y+20.8          # first band starts below the letter clear zone
for i, g in enumerate(DS):
    if i == 4: yy += 4.0
    rows.append((g, yy)); yy += 9.3
F.rect(x, rows[0][1]-4.65, w, rows[3][1]-rows[0][1]+9.3, fc=BAND)
for s_ in (c2, c3):
    F.line([s_.p1]*2, [rows[0][1]-4.65, rows[-1][1]+4.65], lw='ref', color=AXC, ls='ref')
Av = {}
for g, yy in rows:
    nm, col = DATASET[g]; r_ = sup.loc[g]
    F.rect(x+1.0, yy-2.8, 5.6, 5.6, fc=col); F.text(x+9.0, yy-3.4, nm, 7.0)
    s50, tot = r_.supported50_sites/1000, r_.sites/1000
    F.rect(c1(0), yy-3.2, c1(s50)-c1(0), 6.4, fc=SUPC); F.rect(c1(s50), yy-3.2, c1(tot)-c1(s50), 6.4, fc=UNSC)
    F.text(c1(tot)+1.6, yy-3.4, f'{tot:.1f}', 6.8, color=MUTED)
    fr = {wdw: (r_[f'supported{wdw}_sites']/r_.sites, r_[f'supported{wdw}_UMI']/r_.UMI) for wdw in (25, 50, 100)}
    for k_, s_ in ((0, c2), (1, c3)):
        F.line([s_.clip(fr[25][k_]), s_.clip(fr[100][k_])], [yy, yy], lw='range', color=NULLBAR)
        F.marker(s_.clip(fr[50][k_]), yy, 'o', 'small', fc=INK)
    Av[g] = dict(primary_sites=int(r_.sites), supported50_sites=int(r_.supported50_sites), UMI=int(r_.UMI),
                 PAS_fraction={k: float(v[0]) for k, v in fr.items()}, UMI_fraction={k: float(v[1]) for k, v in fr.items()})
ay = rows[-1][1]+5.7
F.xaxis(c1, ay, [0, 40, 80, 120], title='Primary PAS (thousands)')
F.xaxis(c2, ay, [0.6, 0.7, 0.8, 0.9, 1.0], ['.6', '.7', '.8', '.9', '1'], title='Fraction of PAS within window')
F.xaxis(c3, ay, [0.88, 0.92, 0.96, 1.0], ['.88', '.92', '.96', '1'], title='Fraction of UMIs within window')
F.key_row(None, y+H1-7.4, [('block', SUPC, 'Within 50 nt'), ('block', UNSC, 'Other primary PAS')], gap=7.0, center=(c1.p0+c1.p1)/2)
F.key_row(None, y+H1-7.4, [('marker', ('o', INK), '±50 nt'), ('line', (NULLBAR, None, 'range'), '±25 to ±100 nt'),
                            ('line', (AXC, 'ref', 'ref'), 'all')], gap=8.0, center=(c2.p0+c3.p1)/2)
F.values['A'] = Av

# ================================================================ B  canonical signal position
x, y, w, h = F.panel('B', X0, Y2, R-X0, H2)
CATS = [('supported50', INK, None, 'data', 'Supported primary PAS'), ('unsupported50', MUTED, 'dash', 'data', 'Unsupported primary PAS'),
        ('matched_genomic_background', NULLW, 'dot', 'data', 'Matched genomic background')]
by = Lin(0, 0.09, y+66, y+17)
y0x = x+26; pw = (x+w-y0x-6*6.0)/7
Bv = {}
for i, g in enumerate(DS):
    nm, col = DATASET[g]
    s_ = Lin(-100, 100, y0x+(pw+6)*i, y0x+(pw+6)*i+pw)
    F.text((s_.p0+s_.p1)/2, y+3.0, nm, 7.0, ha='center')     # curves are not dataset-coloured, so no colour block
    F.lineplot(s_, by, g)
    F.line([s_.p0, s_.p0], [by.p0, by.p1], lw='axis', color=AXC)
    if i == 0: F.yaxis(by, s_.p0, [0, 0.04, 0.08], ['0', '.04', '.08'], 'Frequency', line=False)
    F.line([s_(0)]*2, [by.p0, by.p1], lw='ref', color=AXC, ls='ref')
    sc = D.df('D04B', f'pas_qc/{DATASET_DIR.get(g, g)}/sequence_curves.tsv.gz')
    sc = sc[(sc.alignment == 'representative') & (sc.feature == 'canonical_signal_start')]
    Bv[g] = {}
    for cat, cc, ls, lw, _ in CATS:
        m = sc[sc.category == cat].sort_values('position_nt')
        F.line([s_(v) for v in m.position_nt], [by.clip(v) for v in m.frequency], lw=lw, color=cc, ls=ls)
        pk = m.loc[m.frequency.idxmax()]
        Bv[g][cat] = dict(peak_position=int(pk.position_nt), peak_frequency=float(pk.frequency), denominator=int(pk.denominator))
    F.xaxis(s_, by.p0+2.0, [-50, 0, 50], ['−50', '0', '50'])
F.text(y0x+(7*pw+36)/2, y+81.0, 'Canonical signal start relative to PAS (nt)', 8.0, ha='center')
F.key_row(None, y+H2-7.4, [('line', (cc, ls, lw), lab) for _, cc, ls, lw, lab in CATS], gap=9.0, center=y0x+(7*pw+36)/2)
F.values['B'] = Bv

# ================================================================ C  read-end positions
# Row 3: two header rows (panel title, subplot titles), two line charts side by side per panel (aspect <= 2).
WC = 222.0
x, y, w, h = F.panel('C', X0, Y3, WC, H3)
PY0, PY1, AXY = y+13.0, y+71.0, y+73.0          # data area top/bottom and x-axis line, shared by C and D
LEG3 = y+H3-7.4
gp = 9.0; pw_ = (x+w-0.3-(x+26)-gp)/2
CX = [Lin(-400, 400, x+26+(pw_+gp)*j, x+26+(pw_+gp)*j+pw_) for j in range(2)]
cy = Lin(0, 0.045, PY1, PY0)
Cv = {}
for j, (cat, lab) in enumerate([('supported50', 'Supported'), ('unsupported50', 'Unsupported')]):
    cx = CX[j]
    F.lineplot(cx, cy, cat)
    F.text((cx.p0+cx.p1)/2, y+3.0, lab, 7.5, ha='center')
    F.line([cx.p0, cx.p0], [cy.p0, cy.p1], lw='axis', color=AXC)
    F.yaxis(cy, cx.p0, [0, 0.02, 0.04], ['0', '.02', '.04'], line=False, labels_off=(j == 1))
    F.line([cx(0)]*2, [cy.p0, cy.p1], lw='ref', color=AXC, ls='ref')
    for g in DS:
        rp = D.df('D04C', f'pas_qc/{DATASET_DIR.get(g, g)}/read_profile_curves.tsv.gz')
        m = rp[rp.category == cat].sort_values('position_nt')
        F.line([cx(v) for v in m.position_nt], [cy.clip(v) for v in m.mean_site_fraction], lw='thin', color=DATASET[g][1])
        Cv.setdefault(g, {})[cat] = dict(max_fraction=float(m.mean_site_fraction.max()),
                                          at_nt=int(m.loc[m.mean_site_fraction.idxmax()].position_nt), sites=int(m.selected_sites.max()))
    F.xaxis(cx, AXY, [-200, 0, 200], ['−200', '0', '200'])
F.vtext(x+4.4, (PY0+PY1)/2+12.6, 'Mean fraction per site', 8.0)
F.text((CX[0].p0+CX[1].p1)/2, AXY+12.0, 'Read end relative to PAS (nt)', 8.0, ha='center')
# 2026-10-12: 'Line colors: studies as in a' moved to the caption
F.values['C'] = Cv

# ================================================================ D  depth and layer sensitivity of complete profiles
xd = X0+WC+GAP; wd = R-xd
x, y, w, h = F.panel('D', xd, Y3, wd, H3)
pos = D.df('D04D', 'profiles/all_route_node_positions.tsv.gz', usecols=['route', 'setting', 'depth', 'gene_id', 'node'])
dxs = {20: 0, 40: 1, 80: 2}
pw_ = (x+w-0.3-(x+24)-gp)/2
DX = [Lin(-0.25, 2.25, x+24+(pw_+gp)*j, x+24+(pw_+gp)*j+pw_) for j in range(2)]
dyy = Lin(0, 5, PY1, PY0)
Dv = {}
for j, (setting, lab) in enumerate((('supported50_p08', 'Supported layer'), ('expanded_primary_p08', 'Expanded layer'))):
    dx = DX[j]
    F.lineplot(dx, dyy, setting)
    F.text((dx.p0+dx.p1)/2, y+3.0, lab, 7.5, ha='center')
    F.line([dx.p0, dx.p0], [dyy.p0, dyy.p1], lw='axis', color=AXC)
    F.yaxis(dyy, dx.p0, [0, 1, 2, 3, 4, 5], line=False, labels_off=(j == 1),
            title='Genes (×1,000)' if j == 0 else None, title_y=(PY0+PY1)/2+3.0)
    for r in FIT_ROUTES:
        nm, col, mk = STUDY[r]
        sub = pos[(pos.route == r) & (pos.setting == setting) & pos.node.isin(list('SREP'))]
        cnt = sub.groupby(['depth', 'gene_id']).node.nunique()
        comp = {int(d_): int((cnt.loc[d_] == 4).sum()) for d_ in (20, 40, 80)}
        F.line([dx(dxs[d_]) for d_ in comp], [dyy(v/1000) for v in comp.values()], lw='thin', color=col)
        for d_, v in comp.items():
            F.marker(dx(dxs[d_]), dyy(v/1000), mk, 'small', fc=col)
        Dv[f'{r}|{setting}'] = comp
    F.xaxis(dx, AXY, [0, 1, 2], ['20', '40', '80'])
F.text((DX[0].p0+DX[1].p1)/2, AXY+12.0, 'Minimum PAS UMIs per metacell', 8.0, ha='center')
F.key_row(None, LEG3, [('marker', (STUDY[r][2], STUDY[r][1]), STUDY[r][0]) for r in FIT_ROUTES], gap=6.0, center=(DX[0].p0+DX[1].p1)/2)
F.notes.append('D: complete-profile counts are the number of genes with all four S/R/E/P nodes per route, setting and depth (display count)')
F.notes.append('C, D: line charts split into side-by-side subplots so every data area has width/height <= 2 (house rule)')
F.values['D'] = Dv

F.save(OUT, D)
