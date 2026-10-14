"""Supplementary Figure D12 (printed number from supdraw.SNUM) — APA–GEX response identity and direction (R3), polyApipe v2.

a  the prespecified eight-test family: tested class (ΔAPA same direction in both studies, ΔGEX opposite) observed / null mean
   under the response-identity and direction-specific nulls, 5 baseline bins (3 bins shown as open markers), both species
b  all four ΔAPA x ΔGEX relation classes for the shared genes, with comparable / shuffleable denominators
c  PAS usage (P vs S, cumulative from the proximal PAS) for AP1S2 (human) and Rtf2 (mouse), each in both studies
d  matched ΔDWUI and ΔGEX (log2 FC) for the same genes and studies
"""
import sys, math, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'D12'
D = Data(); OUT = outdir(FID)
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 132.0, 112.0
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
fam = D.df('D12A', 'apa_gex/prespecified_eight_test_family.tsv')
nul = D.df('D12A', 'apa_gex/four_class_nulls.tsv.gz')
SPC = {'human': STUDY['GSE118258'][1], 'mouse': STUDY['WOT'][1]}

# ================================================================ a  eight-test family
WA = 262.0
x, y, w, h = F.panel('a', X0, Y1, WA, H1)
sx = Lin(0.6, 2.1, x+92, x+w-50)
# 2026-10-12: panel titles moved to axis titles or the caption throughout this figure
F.text(x+w-27, y+3.0, 'q, 5 bins', 7.0, ha='center', color=MUTED)
F.text(x+w-27, y+11.0, '(3 bins)', 7.0, ha='center', color=MUTED)
rows = [(sp, th, nm) for sp in ('human', 'mouse') for th in (0.02, 0.05) for nm in ('response_identity', 'direction_specific')]
NM = {'response_identity': 'Identity', 'direction_specific': 'Direction'}
rh = 10.4; r0 = y+21.0
F.line([sx(1)]*2, [r0-6, r0+rh*7+4+6], lw='ref', color=AXC, ls='ref')
Av = {}
for i, (sp, th, nm) in enumerate(rows):
    yy = r0+rh*i+(4.0 if i >= 4 else 0)
    if i % 4 == 0:
        F.rect(x, yy-rh/2, w, rh*4, fc=BAND) if sp == 'human' else None
        F.text(x+0.4, yy+rh*1.5-3.4, sp.capitalize(), 7.5, font='Bol')
    F.text(sx.p0-4.0, yy-3.4, f'|Δ| ≥ {th:g}'.replace('0.', '.')+' '+NM[nm], 7.0, ha='right')
    for k, ss in enumerate((1, 2)):
        q = nul[(nul.species == sp) & (abs(nul.threshold-th) < 1e-9) & (nul.null_model == nm) & (nul.baseline_bins == 5)
                & (nul.shuffled_study == ss) & (nul.cell == 'B')].iloc[0]
        yk = yy+(k-0.5)*2.6
        F.line([sx(q.null_p025/q.null_mean), sx(q.null_p975/q.null_mean)], [yk, yk], lw='thin', color=NULLBAR)
    f5 = fam[(fam.species == sp) & (abs(fam.threshold-th) < 1e-9) & (fam.null_model == nm) & (fam.baseline_bins == 5)].iloc[0]
    f3 = fam[(fam.species == sp) & (abs(fam.threshold-th) < 1e-9) & (fam.null_model == nm) & (fam.baseline_bins == 3)].iloc[0]
    F.marker(sx(f3.enrichment_min), yy, 'o', 'small', fc='white', ec=SPC[sp])
    F.marker(sx(f5.enrichment_min), yy, 'o', 'mark', fc=SPC[sp])
    F.text(x+w-27, yy-3.4, qfmt(f5.q_family)+f' ({qfmt(f3.q_family)})', 6.8, ha='center')
    Av[f'{sp}|{th}|{nm}'] = dict(observed=int(f5.observed), enrichment_min_B5=float(f5.enrichment_min), q_B5=float(f5.q_family),
                                  enrichment_min_B3=float(f3.enrichment_min), q_B3=float(f3.q_family))
ab = F.xaxis(sx, r0+rh*7+4+rh/2+1.0, [1, 1.5, 2], ['1', '1.5', '2'], title='Observed / null mean')
F.key_row(None, y+h-7.4, [('marker', ('o', INK), '5 bins'), ('marker', ('o', 'white', INK), '3 bins'),
                          ('line', (NULLBAR, None, 'thin'), 'Null 95%, each shuffled study')], gap=6.0, center=(sx.p0+x+w)/2-10)
F.values['a'] = Av
F.notes.append('a: marker = enrichment_min (the smaller of the two shuffled-study enrichments) as in the test family; '
               'bars = null 2.5–97.5% / null mean for each shuffled study, 5 baseline bins')

# ================================================================ b  four relation classes (2 x 2 tiles of main Figure 3e)
x, y, w, h = F.panel('b', X0+WA+GAP, Y1, RR-(X0+WA+GAP), H1)
gx0 = x+52.0; cw, ch, sg = 27.0, 17.0, 8.0
Bv = {}
LIM_E = 0.8                                            # colour = log2(observed / null mean), clipped as in Figure 3e
for k, sp in enumerate(('human', 'mouse')):
    sx0 = gx0+k*(2*cw+sg)
    comp = int(nul[(nul.species == sp) & (nul.baseline_bins == 5)].comparable_genes.iloc[0])
    F.text(sx0+cw, y+3.0, f'{sp.capitalize()} {comp:,}', 7.0, ha='center')
    F.text(sx0+cw/2, y+11.4, 'Same', 6.8, ha='center', color=MUTED); F.text(sx0+1.5*cw, y+11.4, 'Opp.', 6.8, ha='center', color=MUTED)
F.text(gx0-3.0, y+11.4, 'ΔGEX', 6.8, ha='right', color=MUTED)
for b_, th in enumerate((0.02, 0.05)):
    by0 = y+19.4+b_*(2*ch+6.0)
    F.text(x+0.4, by0+ch-3.4, f'{th:g}'.replace('0.', '.'), 7.0, font='Bol')
    for rr, lab in ((0, 'ΔAPA same'), (1, 'ΔAPA opp.')):
        F.text(gx0-3.0, by0+rr*ch+ch/2-3.4, lab, 6.8, ha='right', color=MUTED)
    for k, sp in enumerate(('human', 'mouse')):
        sx0 = gx0+k*(2*cw+sg)
        q = nul[(nul.species == sp) & (abs(nul.threshold-th) < 1e-9) & (nul.null_model == 'response_identity') & (nul.baseline_bins == 5)]
        shuf = min(set(q.shuffled_study), key=lambda s_: float(q[(q.shuffled_study == s_) & (q.cell == 'B')].enrichment.iloc[0]))
        c1 = q[q.shuffled_study == 1].set_index('cell'); cs = q[q.shuffled_study == shuf].set_index('cell')
        assert (c1.observed == cs.observed).all()
        for c_, rr, cc in (('A', 0, 0), ('B', 0, 1), ('C', 1, 0), ('D', 1, 1)):
            e = float(cs.loc[c_].enrichment); o = int(cs.loc[c_].observed)
            fc = diverge(math.log2(e), LIM_E); tx = sx0+cc*cw; ty = by0+rr*ch
            F.rect(tx+0.5, ty+0.5, cw-1.0, ch-1.0, fc=fc, ec=INK if c_ == 'B' else 'none', lw='data' if c_ == 'B' else 0)
            F.text(tx+cw/2, ty+1.6, str(o), 7.0, ha='center', font='Bol' if c_ == 'B' else 'Reg', color=on(fc))
            F.text(tx+cw/2, ty+8.8, f'×{e:.2f}', 6.8, ha='center', color=on(fc))
        Bv[f'{sp}|{th}'] = dict(cells={c_: int(c1.loc[c_].observed) for c_ in 'ABCD'}, comparable=int(c1.comparable_genes.iloc[0]),
                                shuffleable=int(c1.shuffleable_genes.iloc[0]), shuffled_study_shown=int(shuf),
                                enrichment={c_: float(cs.loc[c_].enrichment) for c_ in 'ABCD'})
F.cbar_key(gx0+cw*2+sg/2, y+h-17.4, 'Obs. / null', '×0.57', '×1.74', lambda f_: diverge(-LIM_E+2*LIM_E*f_, LIM_E), wbar=44)
F.key_row(None, y+h-7.4, [('block', 'white', 'Tested class (outlined)')], gap=6.0, center=gx0+cw*2+sg/2)
F.values['b'] = Bv
F.notes.append('b: 2 x 2 tiles as in main Figure 3e: response-identity null, five baseline bins, the shuffled study with the smaller '
               'enrichment of the tested class (same rule as Figure 3e); colour = log2(observed / null mean) clipped at ±0.8; '
               'comparable genes 3,606 (human) / 2,312 (mouse); the tested class is ΔAPA same · ΔGEX opposite')

# ================================================================ c  example PAS usage
re_ = D.df('D12C', 'state_vectors/route_effects.primary.tsv.gz', usecols=['gene_symbol', 'route', 'contrast', 'design', 'setting',
                                                                            'PAS_ids', 'target_usage', 'reference_usage', 'delta_DWUI'])
re_ = re_[(re_.contrast == 'P_S') & (re_.design == 'true_source') & (re_.setting == 'supported50_p08') & re_.gene_symbol.isin(['AP1S2', 'Rtf2'])]
EX = [('AP1S2', ['GSE118258', '178325_0618', '178325_0809']), ('Rtf2', ['CIP', 'WOT'])]
x, y, w, h = F.panel('c', X0, Y2, WA, H2)
pw = (w-29-12-4)/2
cy = Lin(0, 1, y+76, y+20)
Cv = {}


def oriented(ids):
    pos = [int(s.split(':')[-1].split('_')[-2]) for s in ids]; strand = ids[0].split('_')[-1].strip('"')
    p0 = min(pos) if strand == 'f' else max(pos)
    return [abs(p-p0)/1000 for p in pos]


for j, (g, routes) in enumerate(EX):
    q = re_[re_.gene_symbol == g]
    span = max(max(oriented(json.loads(s))) for s in q.PAS_ids)
    ex = Lin(0, math.ceil(span), x+29+(pw+12)*j, x+29+(pw+12)*j+pw)
    F.lineplot(ex, cy, g)
    F.text((ex.p0+ex.p1)/2, y+3.0, g, 8.0, ha='center', font='BolIta')
    F.line([ex.p0, ex.p0], [cy.p0, cy.p1], lw='axis', color=AXC)
    if j == 0: F.yaxis(cy, ex.p0, [0, 0.5, 1], ['0', '50', '100'], 'Distal (%)', line=False, title_y=(cy.p0+cy.p1)/2+4)
    else: F.yaxis(cy, ex.p0, [0, 0.5, 1], line=False, labels_off=True)
    Cv[g] = {}
    for r in routes:
        z = q[q.route == r].iloc[0]; ids = json.loads(z.PAS_ids); o = oriented(ids)
        for st, col_ in (('reference_usage', 'S'), ('target_usage', 'P')):
            u = json.loads(z[st]); order = np.argsort(o); xs_ = np.array(o)[order]; cs = np.cumsum(np.array(u)[order])
            px = [0.0]; py = [0.0]
            for xi, ci in zip(xs_, cs):
                px += [xi, xi]; py += [py[-1], ci]
            px.append(ex.v1); py.append(py[-1])
            F.line([ex(v) for v in px], [cy(1-v) for v in py], lw='data' if col_ == 'P' else 'thin', color=STUDY[r][1],   # distal = 1 − cumulative
                   ls=None if col_ == 'P' else 'dash')
        Cv[g][r] = dict(PAS=len(ids), delta_DWUI=float(z.delta_DWUI))
    F.xaxis(ex, cy.p0+1.0, list(range(0, int(math.ceil(span))+1, 2 if span > 4 else 1)))
F.text(x+29+pw+6, cy.p0+12.6, 'Distance from proximal PAS (kb)', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('line', (INK, 'dash', 'thin'), 'S'), ('line', (INK, None, 'data'), 'P')] +
          [('marker', (STUDY[r][2], STUDY[r][1]), STUDY[r][0].replace('Schiebinger2019–', '').replace('Guo2019–', '').replace('Guan2022–', ''))
           for r in ['GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT']], gap=4.0, size=7.0, center=x+w/2)
F.values['c'] = Cv

# ================================================================ d  matched effects
x, y, w, h = F.panel('d', X0+WA+GAP, Y2, RR-(X0+WA+GAP), H2)
pi = D.df('D12D', 'apa_gex/permutation_input.tsv.gz')
DR = [('AP1S2', 'GSE118258__GSE178325', [('1', 'GSE118258'), ('2', '178325')]), ('Rtf2', 'GSE103221__GSE115943', [('1', 'CIP'), ('2', 'WOT')])]
ax_ = Lin(-0.35, 0.05, x+38, x+90)
gx = Lin(-1.5, 3.5, x+102, x+w-1)
_gex_label = 'ΔGEX, log2(CPM + 1)'   # production label (2026-10-07 label sync); script synced 2026-10-07 terminology pass
Dv = {}; yy = y+18.0
for g, axis, studies in DR:
    q = pi[(pi.gene_symbol == g) & (pi.axis_id == axis) & (pi.setting == 'shared_p08')].iloc[0]
    F.text(x+0.4, yy+4.6-3.4, g, 7.5, font='RegIta')
    for k, (sfx, r) in enumerate(studies):
        yk = yy+k*10.0
        col = STUDY.get(r, ('', DATASET.get('GSE178325')[1], 'X'))[1] if r != '178325' else '#A82C60'
        mk = STUDY.get(r, ('', '', 'X'))[2] if r != '178325' else 'X'
        F.marker(x+31.0, yk, mk, 'small', fc=col)
        for s_, v, lo, hi in ((ax_, q[f'APA_{sfx}'], q[f'alo_{sfx}'], q[f'ahi_{sfx}']), (gx, q[f'GEX_{sfx}'], q[f'glo_{sfx}'], q[f'ghi_{sfx}'])):
            if hi > lo: F.line([s_(lo), s_(hi)], [yk, yk], lw='bold', color=NULLBAR)
            F.marker(s_.clip(v), yk, 'o', 'mark', fc=col)
        Dv[f'{g}|{r}'] = dict(APA=float(q[f'APA_{sfx}']), APA_range=[float(q[f'alo_{sfx}']), float(q[f'ahi_{sfx}'])],
                              GEX=float(q[f'GEX_{sfx}']), GEX_range=[float(q[f'glo_{sfx}']), float(q[f'ghi_{sfx}'])])
    yy += 36.0
for s_ in (ax_, gx):
    F.line([s_(0)]*2, [y+14, yy-8], lw='ref', color=AXC, ls='ref')
F.xaxis(ax_, yy-7.0, [-0.3, -0.15, 0], ['−.3', '−.15', '0'], title='ΔDWUI')
gb = F.xaxis(gx, yy-7.0, [-1, 0, 1, 2, 3])
F.text(min((gx.p0+gx.p1)/2, x+w-0.3-textwidth(_gex_label, 8.0)/2), gb+1.6, _gex_label, 8.0, ha='center')   # axis title kept inside the panel
F.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'Guan2022 series range')], gap=6.0, center=x+w/2)
F.values['d'] = Dv
F.notes.append('d: shared_p08 setting of the APA–GEX analysis; Guan2022 = 0618 and 0809 series (point = study value, bar = series range)')
F.save(OUT, D)
