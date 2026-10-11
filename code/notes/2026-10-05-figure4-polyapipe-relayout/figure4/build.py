"""Figure 4 · polyApipe relayout (2026-10-05). Data from ONE result bundle; v2 primary only; no cis panel.

Usage: build.py --records <bundle>/result-records.json --out <new dir>

Narrative (Result 4): cross-species sharing depends on the study pair (A) and genome-wide continuous
correspondence is small (B); among selected proximal candidates P–M sharing exceeds three nulls (C); on
homologous PAS, P–M correspondence is positive and exceeds P–S (D, the full eight-test same-PAS family with
the binned trends that show it per study pair); PDLIM7 illustrates one gene that sits in both C and D (E).
The cis prediction panel (old 4F) leaves the main figure; its full model belongs to the supplement.

Layout (draft slots; tools-local/retile.py sets the 4.0 mm gutters and equal heights A|B and C+E|D):
  A (0, 0, 302.8, 172)      conditional same-direction sharing, DWUI and RWUI columns; log obs/null axis
  B (302.8, 0, 453.5, 172)  continuous correspondence R/E/P/M relative to observed S, with pair-ρ tiles
  C (0, 172, 240.7, 330)    shared proximal candidates P–S / P–M against the Joint, Quintile and Range nulls, q row
  E (0, 330, 240.7, 593.7)  PDLIM7 two homologous PAS: usage step curves and native ΔDWUI strip
  D (240.7, 172, 453.5, 593.7)  same-PAS forest (PMS: P–S, M–S, P–M ρ, P–M − P–S, P–M − M–S; PEM: E–M, P–M ρ,
                            P–M − E–M), coarse and capacity q tiles, then binned trends P–M, P–S, M–S
Old -> new: A->A, B->B (+M–S row), C->C (+q row), D->D (v3 removed, all eight family tests), E->E, F->supplement.
Typography, colours and helpers are those of the production builder (main-figures-v5/figure4, records adapter).
"""
import csv, gzip, json, sys, math
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools-local'))
from v5draw import *  # noqa
from figrec import Run, MASTERS
R = Run('F4')

MASTER = MASTERS/'Figure.4.pdf'      # locked master copy: canvas size and nothing else (every object is removed)
LOCK_SHA = '6ba591b8e8210808d06520d1da7c3af0c4f19b6d5e0eac3c0a85855c0ba33386'
assert sha(MASTER) == LOCK_SHA, 'master changed'
OUT = R.out/'Figure.4.draft.pdf'
use = R.rows

pl = pdfplumber.open(MASTER).pages[0]; W, H = pl.width, pl.height
pdf = pikepdf.open(MASTER)
removed = {'all': len(clear(pdf, H, [(-5, -5, W+5, H+5)]))}
t = Top(W, H)
ver = {'master_sha256': LOCK_SHA, 'removed_objects': removed, 'info': {}}
INFO = ver['info']
def info(panel, mark, n, channels):
    INFO.setdefault(panel, []).append(dict(mark=mark, n=int(n), channels=int(channels)))
PANELS4 = dict(A=(0, 0, 302.8, 172.0), B=(302.8, 0, 453.5, 172.0), C=(0, 172.0, 240.7, 330.0),
               E=(0, 330.0, 240.7, 593.7), D=(240.7, 172.0, 453.5, 593.7))
ver['panels'] = PANELS4
for L, (x, top) in dict(A=(4.2, 9.1), B=(304.3, 7.3), C=(4.2, 173.3), D=(242.3, 173.3), E=(4.2, 331.3)).items():
    t.text_top(x, top, L, 11, font='Bol')
SCOL = {'GSE118258': '#753399', 'GSE150311': '#648290', '178325_0618': '#A82C60', '178325_0809': '#C47A2C',
        'GSE221739': '#30363D', 'CIP': '#CF4C00', 'WOT': '#207450'}
SNAME = {'GSE118258': 'Xing2020', 'GSE150311': 'Liu2020', '178325_0618': 'Guan2022–0618', '178325_0809': 'Guan2022–0809',
         'GSE221739': 'Panariello2023', 'CIP': 'Guo2019–CIP', 'WOT': 'Schiebinger2019–2i'}
def block(x, ytop, key, w=6.0, h=6.0):
    t.rect(x, ytop, w, h, fc=SCOL[key])
def key_rows(x0, ytop, keys, xmax, size=7.5, gap=6.0, pitch=10.0, names=None):
    """Study colour key wrapped onto as many rows as needed between x0 and xmax; returns the last row's top."""
    x, y = x0, ytop
    for k in keys:
        nm = (names or {}).get(k, SNAME[k]); w_ = 7.6+textwidth(nm, size)
        if x > x0 and x+w_ > xmax:
            x, y = x0, y+pitch
        block(x, y+1.0, k, 5.6, 5.6); t.text_top(x+7.6, y, nm, size); x += w_+gap
    return y
def pair_blocks(x, ytop, h, m):
    block(x, ytop, h); t.text_top(x+7.2, ytop-1.2, '×', 7.5, color=MUTED); block(x+13.2, ytop, m); return x+19.2

# ---------------------------------------------------------------- A
rows = use('F4A', 'cross_sharing/sharing_summary.tsv')
pairs = [('GSE118258', 'GSE103221'), ('GSE118258', 'GSE115943'), ('GSE178325', 'GSE103221'), ('GSE178325', 'GSE115943')]
PK = {'GSE118258': 'GSE118258', 'GSE178325': '178325_0618', 'GSE103221': 'CIP', 'GSE115943': 'WOT'}
ys = [34.8, 46.9, 59.0, 71.1, 83.2, 95.4, 107.5, 119.6]
BLK = [('WPP', 'DWUI', 30.0), ('RWUI', 'RWUI', 171.0)]
RLO, RHI = 0.6, 2.0                       # obs/null on a log axis (polyApipe Guan ratios reach 1.94)
PROXC, DISTC = '#0072B2', '#D55E00'
Ad = []
shared = {(r['metric'], r['human_study'], r['mouse_study'], round(float(r['threshold']), 3), r['test']): int(r['observed']) for r in rows}
for i, (h, m) in enumerate(pairs):
    yc = (ys[2*i]+ys[2*i+1])/2
    t.rect(4.0, ys[2*i]-5.5, 294.0, ys[2*i+1]-ys[2*i]+11.0, fc=BAND) if i % 2 == 0 else None
    block(8.0, yc-3.0, PK[h]); block(16.0, yc-3.0, PK[m])
for metric, head, b in BLK:
    X = lambda v, b=b: b+(math.log(v)-math.log(RLO))/(math.log(RHI)-math.log(RLO))*52.0
    xq = b+52.0+3.0+15.0; xb = b+90.0; bw = 18.0
    t.text_top(b+62.0, 5.0, head, 8, ha='center')
    t.text_top(X(1.0), 20.1, 'Obs./null', 7.5, ha='center')
    t.text_top(xq-3.0, 19.1, 'q', 7.5); t.text_top(xq+1.2, 21.2, '8', 6.8)
    t.text_top(xb+bw/2+5.0, 20.1, 'Shared', 7.5, ha='center')
    t.line([X(1.0), X(1.0)], [28.8, 125.6], lw=0.6, color='#AAB4BB')
    t.line([X(RLO), X(RHI)], [125.6, 125.6], lw=0.5, color='#AAB4BB')
    for v, lab in [(0.7, '0.7'), (1.0, '1'), (1.4, '1.4'), (2.0, '2')]:
        t.line([X(v)]*2, [125.6, 127.8], lw=0.6, color=INK)
        t.text_top(X(v), 129.6, lab, 7.5, ha='center')
    k = 0
    for h, m in pairs:
        for thr in ('0.02', '0.05'):
            r = [r for r in rows if r['metric'] == metric and r['test'] == 'same_direction' and r['human_study'] == h
                 and r['mouse_study'] == m and abs(float(r['threshold'])-float(thr)) < 1e-9]
            assert len(r) == 1; r = r[0]
            mean = float(r['null_mean']); lo = float(r['null_p025'])/mean; hi = float(r['null_p975'])/mean
            enr = float(r['enrichment']); q8 = float(r['q_primary_8']); q24 = float(r['q_joint_primary_24'])
            if abs(enr-float(r['observed'])/mean) > 1e-9: R.flag('F4A', 'consistency', f'{metric} {h}x{m} {thr}: enrichment != observed/null mean')
            if not R.axis('F4A', f'{metric} {h}x{m} {thr} obs/null', [lo, hi, enr], RLO, RHI):
                lo, hi, enr = (min(max(v_, RLO), RHI) for v_ in (lo, hi, enr))
            y = ys[k]; k += 1
            t.line([X(lo), X(hi)], [y, y], lw=1.4, color=NULLBAR)
            filled = q8 < 0.05
            t.marker(X(enr), y, 'o' if thr == '0.02' else 'D', 3.2 if thr == '0.02' else 3.6, fc=INK if filled else 'white', ec=INK, mew=0.65)
            lab = f'{q8:.4f}'.lstrip('0') + (' †' if q24 < 0.05 else '')   # 2026-10-08: q without leading zero, as elsewhere
            if q24 < 0.05:
                t.text_top(xq-15.0, y-4.2, lab, 8); assert xq-15.0+textwidth(lab, 8) < xb-1.0
            else:
                t.text_top(xq, y-4.2, lab, 8, ha='center')
            key = (metric, h, m, round(float(thr), 3))
            nneg, npos, nsame = shared[key+('shared_negative',)], shared[key+('shared_positive',)], shared[key+('same_direction',)]
            if not nneg+npos == nsame == int(r['observed']): R.flag('F4A', 'consistency', f'{metric} {h}x{m} {thr}: shared counts {nneg}+{npos} vs {nsame}/{r["observed"]}')
            den_ = max(nneg+npos, 1)
            t.rect(xb, y-2.1, bw*nneg/den_, 4.2, fc=PROXC); t.rect(xb+bw*nneg/den_, y-2.1, bw*npos/den_, 4.2, fc=DISTC)
            t.text_top(xb+bw+1.6, y-3.8, str(nsame), 7.5, color=MUTED)
            assert xb+bw+1.6+textwidth(str(nsame), 7.5) < 300.0
            Ad.append(dict(metric=metric, human=h, mouse=m, threshold=thr, observed=int(r['observed']), null_mean=mean,
                           ratio=enr, null_lo=lo, null_hi=hi, q8=q8, q24=q24, filled=filled, dagger=q24 < 0.05,
                           shared_proximal=nneg, shared_distal=npos))
info('A', 'forest: observed ratio, null interval, fill = q8 < .05', 16, 4)
info('A', 'q8 values (table)', 16, 1)
info('A', 'stacked bar: proximal / distal shares + count', 16, 2)
x = 8.0
t.text_top(x, 141.9, 'Threshold', 7.5); x += textwidth('Threshold', 7.5)+4.0
t.marker(x+2, 145.7, 'o', 3.2, fc=INK, ec=INK, mew=0.65); t.text_top(x+6.0, 141.9, '0.02', 7.5); x += 6.0+textwidth('0.02', 7.5)+7.0
t.marker(x+2, 145.7, 'D', 3.6, fc=INK, ec=INK, mew=0.65); t.text_top(x+6.0, 141.9, '0.05', 7.5); x += 6.0+textwidth('0.05', 7.5)+9.0
t.rect(x, 143.6, 5.0, 4.4, fc=PROXC); t.text_top(x+7.0, 141.9, 'Proximal', 7.5); x += 7.0+textwidth('Proximal', 7.5)+6.0
t.rect(x, 143.6, 5.0, 4.4, fc=DISTC); t.text_top(x+7.0, 141.9, 'Distal', 7.5); x += 7.0+textwidth('Distal', 7.5)+8.0
t.text_top(x, 141.9, '† cross-metric q', 7.5); x2 = x+textwidth('† cross-metric q', 7.5)
t.text_top(x2+0.4, 143.8, '24', 6.8); t.text_top(x2+0.8+textwidth('24', 6.8), 141.9, ' < .05', 7.5)
assert x2+0.8+textwidth('24 < .05', 7.5) < 292.0
key_rows(8.0, 152.0, ['GSE118258', '178325_0618', 'CIP', 'WOT'], xmax=292.0, names={'178325_0618': 'Guan2022'})
ver['A'] = Ad

# ---------------------------------------------------------------- B
cont = {r['state']: r for r in use('F4B', 'cross_continuous/continuous_gene_block_summary.tsv') if r['layer'] == 'supported50'}
pairs_ = use('F4B', 'cross_continuous/continuous_pair_summary.tsv')
seen = set(); pairrows = []
for r in pairs_:
    k_ = (r['layer'], r['state'], r['mouse'], r['human'])
    if r['layer'] == 'supported50' and k_ not in seen:
        seen.add(k_); pairrows.append(r)
b0, b1 = -0.06, 0.10
XB = lambda v: 330.0+(v-b0)/(b1-b0)*58.0
STATES = [('R', 'R–S'), ('E', 'E–S'), ('P', 'P–S'), ('M', 'M–S')]
STATES = [s_ for s_ in STATES if s_[0] in cont]
for s_ in set(cont)-{s for s, _ in STATES}: R.flag('F4B', 'row_capacity', f'state {s_} in the bundle has no row')
YB = {st: 44.0+15.0*i for i, (st, _) in enumerate(STATES)}
ybot = YB[STATES[-1][0]]+8.0
t.line([XB(0)]*2, [37.0, ybot], lw=0.5, color=AXC, ls=(0, (1, 1.5)))
Bd = {}
HUM = [('H1', 'GSE118258'), ('H2', 'GSE150311'), ('H3', '178325_0618'), ('H4', 'GSE221739')]
PCOL = [('C', h) for h, _ in HUM]+[('W', h) for h, _ in HUM]
gx0, gcw, GL = 396.0, 6.6, 0.20
GX = lambda k: gx0+k*gcw+(1.4 if k >= 4 else 0.0)
for i, (st, lab) in enumerate(STATES):
    y = YB[st]
    t.rect(305.0, y-7.5, 147.0, 15.0, fc=BAND) if i % 2 == 0 else None
    t.text_top(306.5, y-3.8, lab, 7.5)
    r = cont[st]; rho, lo, hi = float(r['rho']), float(r['lo']), float(r['hi'])
    pr_ = {(x['mouse'], x['human']): float(x['rho']) for x in pairrows if x['state'] == st}
    if len(pr_) != int(r['pairs']): R.flag('F4B', 'consistency', f'{st}: {len(pr_)} pair rows vs pairs={r["pairs"]}')
    if not R.axis('F4B', f'{st} pooled rho', [lo, hi, rho], b0, b1): lo, hi, rho = (min(max(v_, b0), b1) for v_ in (lo, hi, rho))
    R.axis('F4B', f'{st} pair rho (colour)', list(pr_.values()), -GL, GL)
    for pc_ in set(pr_)-set(PCOL): R.flag('F4B', 'row_capacity', f'{st}: pair {pc_} has no tile column')
    t.whisker(XB(lo), XB(hi), y, color=INK, lw=0.8, cap=1.3); t.marker(XB(rho), y, 'o', 3.3, zorder=4)
    for k, pc in enumerate(PCOL):
        v = pr_.get(pc)
        t.rect(GX(k), y-6.5, gcw-0.7, 13.0, fc=MISSING if v is None else diverge(v, GL))
    Bd[st] = dict(rho=rho, lo=lo, hi=hi, union_genes=int(r['union_genes']), pairs=int(r['pairs']),
                  pair_rhos={f'{a}|{b}': v for (a, b), v in pr_.items()})
info('B', 'forest: pooled rho with CI', len(STATES), 3)
info('B', 'tile grid: per-pair rho (colour)', len(STATES)*8, 1)
for k_, (mk, x0_) in enumerate([('CIP', GX(0)), ('WOT', GX(4))]):
    t.rect(x0_, 20.0, 4*gcw-0.7, 4.5, fc=SCOL[mk])
for k, (mo, hu) in enumerate(PCOL):
    t.rect(GX(k), 27.0, gcw-0.7, 4.5, fc=SCOL[dict(HUM)[hu]])
t.text_top(GX(0)-2.0, 19.0, 'Mouse', 7.0, ha='right', color=MUTED); t.text_top(GX(0)-2.0, 26.0, 'Human', 7.0, ha='right', color=MUTED)
ya = ybot+3.5
t.line([XB(b0), XB(b1)], [ya, ya], lw=0.5, color=AXC)
for v, lab in [(-0.04, '−.04'), (0, '0'), (0.08, '+.08')]:
    t.line([XB(v)]*2, [ya, ya+2], lw=0.5, color=AXC); t.text_top(XB(v), ya+3.8, lab, 7.5, ha='center')
t.text_top(XB(0.02), ya+13.0, 'Mean Spearman ρ', 7.5, ha='center')
CBW = 48.0
for i_ in range(40):
    v = -GL+2*GL*(i_+0.5)/40
    t.rect(gx0+i_*CBW/40, ya, CBW/40+0.1, 3.0, fc=diverge(v, GL))
for v, lab in [(-GL, '−.2'), (0, '0'), (GL, '+.2')]:
    t.text_top(gx0+(v+GL)/(2*GL)*CBW, ya+4.6, lab, 7.0, ha='center', color=MUTED)
t.text_top(gx0+CBW/2, ya+13.0, 'Pair ρ', 7.5, ha='center')
keys = ['GSE118258', 'GSE150311', '178325_0618', 'GSE221739', 'CIP', 'WOT']
yk0 = ya+28.0
for i_, kk in enumerate(keys):
    col, row = i_ // 3, i_ % 3
    x = 306.0+col*66.0; yt = yk0+row*11.0
    block(x, yt+1.0, kk, 5.6, 5.6); t.text_top(x+7.6, yt, SNAME[kk], 7.5)
    assert x+7.6+textwidth(SNAME[kk], 7.5) < 452.0
assert yk0+2*11.0+9.0 < 172.0, yk0
ver['B'] = Bd

# ---------------------------------------------------------------- C
ov = use('F4C', 'cross_overlap/overlap_null_results.tsv')
NUL = [('joint_measurement_PAS', 'Joint', '^'), ('independent_q5_reference_capacity', 'Quintile', 'o'),
       ('independent_range10_reference_capacity', 'Range', 's')]
y0c, y1c, gmax = 199.0, 271.0, 36
YC = lambda g: y1c-g/gmax*(y1c-y0c)
Cd = {}
t.vtext(22.0, (y0c+y1c)/2, 'Shared proximal', 7.5)
for j, (con, head) in enumerate([('PS', 'P–S'), ('PM', 'P–M')]):
    xl = 46.0+j*96.0; w = 80.0
    t.text_top(xl+w/2, 175.0, head, 8, ha='center')
    t.line([xl, xl], [y0c, y1c], lw=0.5, color=AXC); t.line([xl, xl+w], [y1c, y1c], lw=0.5, color=AXC)
    for g in (0, 10, 20, 30):
        t.line([xl-2.1, xl], [YC(g)]*2, lw=0.5, color=AXC)
        t.text_top(xl-3.4, YC(g)-3.6, str(g), 7.5, ha='right')
    d = [r for r in ov if r['contrast'] == con]
    obs = {int(r['observed']) for r in d}; assert len(obs) == 1; obs = obs.pop()
    t.line([xl, xl+w], [YC(obs)]*2, lw=0.9, color=INK)
    t.text_top(xl+w/2, 184.6, f'{obs} observed', 7.5, ha='center', color=MUTED)
    nsig = 0; Cd[con] = {'observed': obs}
    for k, (model, lab, mk) in enumerate(NUL):
        r = [r for r in d if r['model'] == model]; assert len(r) == 1; r = r[0]
        xc = xl+w*(1+2*k)/6; m = float(r['null_mean']); lo = float(r['null_q025']); hi = float(r['null_q975'])
        R.axis('F4C', f'{con} {model} null range / observed', [hi, lo, obs], 0, gmax)
        q = float(r['q_within_model_two_contrasts'])
        t.vwhisker(xc, YC(hi), YC(lo)); t.marker(xc, YC(m), mk, 3.6, fc=INK if q < 0.05 else 'white', ec=NULLW if q >= 0.05 else INK, mew=1.0)
        t.text_top(xc, y1c+3.0, lab, 7.5, ha='center')
        t.text_top(xc, y1c+12.5, qfmt(q), 7.0, ha='center', color=MUTED)
        Cd[con][lab] = dict(null_mean=m, q025=lo, q975=hi, p=float(r['empirical_p']), q=q)
        nsig += q < 0.05
    if j == 0:
        t.text_top(xl-3.4, y1c+12.5, 'q', 7.0, ha='right', color=MUTED)
    Cd[con]['nulls_q_lt_05'] = nsig
info('C', 'null mean + 95% range per named null, fill = q < .05', 6, 4)
info('C', 'observed count line', 2, 1)
info('C', 'q per null (table)', 6, 1)
ver['C'] = Cd

# ---------------------------------------------------------------- D
bank = use('F4D', 'same_pas/v2_primary_correlation_summary.tsv', lambda r: r['tolerance'] in ('25', '25.0') and r['axis'] == 'rank' and r['pair'] == 'pair_mean')
coarse = use('F4D', 'same_pas/v2_primary_conditional_null_tests.tsv')
capac = use('F4D', 'same_pas/capacity_sensitivity_tests.tsv')
if any(r['variant'] != 'v2_primary' for r in bank+coarse+capac):
    R.flag('F4D', 'contract', 'bundle carries a non-v2 annotation variant; this figure draws v2 primary only')
def est(design, metric):
    r = [r for r in bank if r['variant'] == 'v2_primary' and r['design'] == design and r['metric'] == metric]
    assert len(r) == 1, (design, metric); r = r[0]; return float(r['estimate']), float(r['low']), float(r['high']), int(r['n_genes'])
def qv(tab, design, metric):
    col = 'q_BH_family' if 'q_BH_family' in tab[0] else 'q_BH_family8'
    r = [r for r in tab if r['variant'] == 'v2_primary' and r['design'] == design and r['metric'] == metric
         and r.get('tolerance', '25') in ('25', '25.0')]
    assert len(r) == 1, (design, metric, len(r)); return float(r[0][col])
# the eight tests of the same-PAS family, in two cohorts
# row labels name the correlation explicitly (review 2026-10-07): differences are rho(P-M) - rho(reference)
BLOCKS = [('PMS', [('PS_rho', 'ρ(P–S)'), ('MS_rho', 'ρ(M–S)'), ('PM_rho', 'ρ(P–M)'), ('PM_minus_PS_rho', 'ρ(P–M) − ρ(P–S)'),
                   ('PM_minus_MS_rho', 'ρ(P–M) − ρ(M–S)')]),
          ('PEM', [('EM_rho', 'ρ(E–M)'), ('PM_rho', 'ρ(P–M)'), ('PM_minus_EM_rho', 'ρ(P–M) − ρ(E–M)')])]
fam = {(r['design'], r['metric']) for r in coarse if r['variant'] == 'v2_primary' and r.get('tolerance') in ('25', '25.0')}
drawn = {(d_, m_) for d_, rr in BLOCKS for m_, _ in rr}
if fam != drawn: R.flag('F4D', 'row_capacity', f'coarse family {sorted(fam)} differs from the drawn rows {sorted(drawn)}')
xlab = 256.0; vmin, vmax = -0.4, 0.5
LABW = max(textwidth(l_, 8) for _, rr_ in BLOCKS for _, l_ in rr_)
xl0, xl1 = xlab+LABW+6.0, 386.0   # forest starts after the longest row label
assert xl1-xl0 >= 60.0, (LABW, xl0)
XD = lambda v: xl0+(v-vmin)/(vmax-vmin)*(xl1-xl0)
QX = {'coarse': 396.0, 'capacity': 429.0}; QW = 10.0
for nm, lab in [('coarse', 'Coarse'), ('capacity', 'Capacity')]:
    t.text_top(QX[nm]+QW/2, 175.6, lab, 7.5, ha='center')
assert QX['capacity']+QW/2+textwidth('Capacity', 7.5)/2 < 452.0
PITCH = 13.0
Dd = []; y = 176.0
for bi, (design, rr) in enumerate(BLOCKS):
    if design == 'PMS':
        x = pair_blocks(262.0, y+0.8, 'GSE118258', 'CIP'); x = pair_blocks(x+6.0, y+0.8, 'GSE118258', 'WOT')
    else:
        x = pair_blocks(262.0, y+0.8, 'GSE150311', 'CIP')
    n_ = est(design, rr[0][0])[3]
    t.text_top(x+5.0, y, f'n = {n_}', 7.5, color=MUTED)
    ytop = y+16.0
    for ri, (metric, lab) in enumerate(rr):
        yy = ytop+ri*PITCH
        if ri % 2 == 0:
            t.rect(xlab-2.0, yy-6.0, 196.0, 12.0, fc=BAND)
        t.text_top(xlab, yy-4.0, lab, 8)
        e, lo, hi, n = est(design, metric)
        if n != n_: R.flag('F4D', 'consistency', f'{design} {metric}: n {n} vs block n {n_}')
        if not R.axis('F4D', f'{design} {metric}', [lo, hi, e], vmin, vmax): e, lo, hi = (min(max(v_, vmin), vmax) for v_ in (e, lo, hi))
        t.whisker(XD(lo), XD(hi), yy, color=INK, lw=0.8, cap=1.2)
        t.marker(XD(e), yy, 'o', 3.0, fc=INK, ec=INK, mew=0.8)
        qc, qk = qv(coarse, design, metric), qv(capac, design, metric)
        for nm, q_ in [('coarse', qc), ('capacity', qk)]:
            t.rect(QX[nm], yy-4.0, QW, 8.0, fc=INK if q_ < 0.05 else 'white', ec=INK if q_ < 0.05 else NULLW, lw=0.6)
        Dd.append(dict(design=design, metric=metric, estimate=e, low=lo, high=hi, n_genes=n, q_coarse=qc, q_capacity=qk))
    yend = ytop+(len(rr)-1)*PITCH
    t.line([XD(0)]*2, [ytop-7.0, yend+7.0], lw=0.5, color=AXC, ls=(0, (1, 1.5)))
    y = yend+12.0
info('D', 'forest: estimate with CI (eight-test family)', 8, 3)
info('D', 'q tiles (coarse, capacity): q8 < .05', 16, 1)
ya = y-3.0
t.line([XD(vmin), XD(vmax)], [ya, ya], lw=0.5, color=AXC)
for v, lab in [(-0.2, '−.2'), (0, '0'), (0.2, '+.2'), (0.4, '+.4')]:
    t.line([XD(v)]*2, [ya, ya+2], lw=0.5, color=AXC); t.text_top(XD(v), ya+3.6, lab, 7.5, ha='center')
_tw = textwidth('Spearman ρ or paired Δρ', 8)   # keep the title clear of the q legend squares
t.text_top(min((xl0+xl1)/2, QX['coarse']-4.0-_tw/2), ya+12.6, 'Spearman ρ or paired Δρ', 8, ha='center')
x = (QX['coarse']+QX['capacity']+QW)/2-22.0; yq = ya+3.4
t.rect(QX['coarse'], yq+1.0, 6.0, 5.6, fc=INK, ec=INK, lw=0.6); x = QX['coarse']+8.0
t.text_top(x, yq, 'q', 7.5); x2 = x+textwidth('q', 7.5)
t.text_top(x2+0.4, yq+2.1, '8', 6.8); t.text_top(x2+0.8+textwidth('8', 6.8), yq, ' < .05', 7.5)
t.rect(QX['coarse'], yq+11.0, 6.0, 5.6, fc='white', ec=NULLW, lw=0.6); t.text_top(QX['coarse']+8.0, yq+10.0, '≥ .05', 7.5)
assert x2+0.8+textwidth('8 < .05', 7.5) < 452.0
ylast = key_rows(xlab, ya+25.0, ['GSE118258', 'GSE150311', 'CIP', 'WOT'], xmax=452.0, gap=5.0)
# binned trend: mean mouse rank effect within equal-count bins of the human rank effect, per study pair and contrast
import pandas as pd
rbins = R.df('F4D', 'panel_statistics/same-pas-rank-bins.tsv'); rsum = R.df('F4D', 'panel_statistics/same-pas-rank-bin-summary.tsv')
COLS = [('PM', 'P–M'), ('PS', 'P–S'), ('MS', 'M–S')]
PAIRS = [('H1_C', 'CIP', 's', -1.2), ('H1_W', 'WOT', 'o', 1.2)]
NB = 8; YLO, YHI = -0.18, 0.10
gx0, gw, ggap = 278.0, 52.0, 7.0
gy0, gh = ylast+10.0+22.0, 70.0
Gd = []
for ci, (k, klab) in enumerate(COLS):
    cx = gx0+ci*(gw+ggap)
    X = lambda b_, cx=cx: cx+(b_+0.5)/NB*gw
    Y = lambda v: gy0+gh-(np.clip(v, YLO, YHI)-YLO)/(YHI-YLO)*gh
    t.text_top(cx+gw/2, gy0-20.5, klab, 8, ha='center')
    rx_ = []
    t.rect(cx, gy0, gw, gh, fc='white', ec='#C9D1D6', lw=0.5)
    t.line([cx, cx+gw], [Y(0)]*2, lw=0.4, color='#C9D1D6', ls=(0, (1, 1.5)))
    for pi, (pair, mk, mshape, off) in enumerate(PAIRS):
        bb = rbins[(rbins.pair == pair) & (rbins.contrast == k) & (rbins.status == 'estimated')].sort_values('bin')
        ss_ = rsum[(rsum.pair == pair) & (rsum.contrast == k)]
        if len(bb) != NB or len(ss_) != 1:
            R.unavailable('F4D', f'trend {k} {pair}: {len(bb)} of {NB} bins estimated; series not drawn'); continue
        rho = float(ss_.rho.iloc[0]); nR = int(ss_.n.iloc[0])
        bx_, by_, be_ = bb.x_mean.astype(float).tolist(), bb.y_mean.astype(float).tolist(), bb.y_se.astype(float).tolist()
        R.axis('F4D', f'trend {k} {pair} bin mean ± s.e.', [v-e for v, e in zip(by_, be_)]+[v+e for v, e in zip(by_, be_)], YLO, YHI)
        col = SCOL[mk]
        t.line([X(b_)+off for b_ in range(NB)], [Y(v) for v in by_], lw=0.8, color=col, alpha=0.9)
        for b_, (yv, ev) in enumerate(zip(by_, be_)):
            t.line([X(b_)+off]*2, [Y(yv-ev), Y(yv+ev)], lw=0.6, color=col)
            t.marker(X(b_)+off, Y(yv), mshape, 3.0, fc=col, ec='white', mew=0.3, zorder=5)
        rx_.append((num(rho, 2, plus=True), col))
        Gd.append(dict(contrast=k, pair=pair, n=nR, rho=rho, bins=NB, bin_x=bx_, bin_mean=by_, bin_se=be_))
    wr = textwidth('ρ ', 7.5)+sum(textwidth(a_, 7.5) for a_, _ in rx_)+4.0; xr = cx+gw/2-wr/2
    t.text_top(xr, gy0-10.5, 'ρ', 7.5); xr += textwidth('ρ ', 7.5)
    for a_, c_ in rx_:
        t.text_top(xr, gy0-10.5, a_, 7.5, color=c_); xr += textwidth(a_, 7.5)+4.0
    assert wr < gw+6.0, wr
    for b_ in (0, 3, 7):
        t.line([X(b_)]*2, [gy0+gh, gy0+gh+2], lw=0.5, color=AXC); t.text_top(X(b_), gy0+gh+3.4, str(b_+1), 7.5, ha='center')
    if ci == 0:
        for v, s_ in [(-0.1, '−.1'), (0, '0'), (0.1, '+.1')]:
            t.line([cx-2.0, cx], [Y(v)]*2, lw=0.5, color=AXC); t.text_top(cx-2.8, Y(v)-3.8, s_, 7.5, ha='right')
assert gx0+3*gw+2*ggap < 449.0
info('D', 'binned trend: bin means with s.e. (y, s.e.)', 3*2*NB, 2)
info('D', 'per-pair rho annotation', 6, 1)
t.vtext(254.0, gy0+gh/2, 'Mouse rank effect', 7.5)
gxe = gx0+3*gw+2*ggap
t.text_top((gx0+gxe)/2, gy0+gh+13.0, 'Human rank-effect bin', 8, ha='center')
yk = gy0+gh+24.0
x = xlab; t.text_top(x, yk, 'Xing2020 ×', 7.5); x += textwidth('Xing2020 ×', 7.5)+5.0
for mk_, shp_, lab_ in [('CIP', 's', 'Guo2019–CIP'), ('WOT', 'o', 'Schiebinger2019–2i')]:
    w_ = 12.0+textwidth(lab_, 7.5)
    if x+w_ > 452.0:
        x, yk = xlab+textwidth('Xing2020 ×', 7.5)+5.0, yk+10.0
    t.line([x, x+9.0], [yk+3.8]*2, lw=0.8, color=SCOL[mk_]); t.marker(x+4.5, yk+3.8, shp_, 3.0, fc=SCOL[mk_], ec='white', mew=0.3, zorder=5)
    t.text_top(x+12.0, yk, lab_, 7.5); x += w_+8.0
ver['D_trend'] = Gd
ver['D'] = Dd
ver['D_bottom'] = yk+9.0

# ---------------------------------------------------------------- E
R.optional('F4E', 'geometry/PAS-transcript-coordinate-check.tsv.gz')   # registered structural annotation; not drawn
specE = R.selection.get('homolog_example', {'mode': 'rank'})
hp_ = R.optional('F4E', 'axis_extensions/homolog_site_pairs.tsv.gz')
hom = use('F4E', 'axis_extensions/homolog_site_pairs.tsv.gz', lambda r: r['tolerance'] in ('25', '25.0')) if hp_ else []
NEED = {('H1_C', 'qualified_P_S'), ('H1_C', 'qualified_P_M'), ('H1_W', 'qualified_P_S'), ('H1_W', 'qualified_P_M')}
exE = None
if hp_:
    _have = {}
    for r in hom:
        _have.setdefault((r['human_gene_id'], r['mouse_gene_id']), set()).add((r['pair'], r['design']))
    eligE = [dict(gene_id=h, mouse_gene_id=m) for (h, m), s_ in _have.items() if NEED <= s_]
    if specE.get('mode') == 'pin':
        exE = R.example('homolog_example', 'F4E', eligE)
        if exE and specE.get('mouse_gene_id') not in (None, exE['mouse_gene_id']):
            R.flag('F4E', 'example_ineligible', 'pinned mouse ortholog differs from the homologous pair'); exE = None
    else:
        R.unavailable('F4E', 'rank rule needs a statistics-side homologous example export (GAP-D04); pin a pair in display-selection.json')
HG, MG = (exE['gene_id'], exE['mouse_gene_id']) if exE else (specE.get('gene_id'), specE.get('mouse_gene_id'))
mem = use('F4E', 'cross_trace/member_endpoint_effects.tsv.gz', lambda r: r['human_gene_id'] == HG and r['setting'] == 'supported50_p08')
eff = {(r['contrast'], r['route']): r for r in mem}
h3b = [r for r in hom if r['human_gene_id'] == HG and r['design'] == 'qualified_P_S' and r['pair'] == 'H3b_C']
_strand = {r['human_id'].split('|')[-2] for r in hom if r['human_gene_id'] == HG}
sites = sorted({(int(r['human_id'].split('|')[-1]), int(r['mouse_id'].split('|')[-1])) for r in hom if r['human_gene_id'] == HG}, reverse=_strand != {'+'})
E_OK = exE is not None and len(sites) == 2
if exE is not None and len(sites) != 2:
    R.unavailable('F4E', f'{specE.get("gene_symbol", HG)}: {len(sites)} homologous PAS (layout shows exactly two); usage curves left empty')
if E_OK:
    (hp, mp), (hd, md) = sites
    hkb, mkb = abs(hp-hd)/1000, abs(mp-md)/1000
    R.derive('F4E', 'distal share of the two homologous PAS renormalised within the pair (GAP-D04)')
else:
    hp = mp = hd = md = hkb = mkb = None
_g = specE.get('gene_symbol', ''); t.text_top(21.0, 332.0, f'{_g}/{_g.capitalize()}' if _g else '', 8, font='BolIta')   # 2026-10-08: human/mouse symbols
def usage(pair, design, side):
    out = {}
    for r in hom:
        if r['human_gene_id'] == HG and r['pair'] == pair and r['design'] == design:
            pos = int(r[f'{"human" if side == "h" else "mouse"}_id'].split('|')[-1])
            out[pos] = (float(r[f'reference_{side}']), float(r[f'target_{side}']))
    return out
SC_ = dict(S='#E69F00', P='#CC79A7', M='#8C6D31')
LS_ = dict(S=(0, (3, 1.6)), P='-', M=(0, (3, 1.2, 1, 1.2)))
PLOTS = [('GSE118258', 'H1_C', 'h', hp, hd, hkb), ('CIP', 'H1_C', 'm', mp, md, mkb), ('WOT', 'H1_W', 'm', mp, md, mkb)]
ex0, ew, egap = 30.0, 58.0, 10.0
ey0, eh = 355.0, 46.0
EU = {}
for i, (skey, pair, side, pprox, pdist, kb) in enumerate(PLOTS):
    x0 = ex0+i*(ew+egap)
    XK = lambda v, x0=x0: x0+v/9.5*ew
    YK = lambda v: ey0+eh-v/100*eh
    t.rect(x0, 346.2, ew, 4.5, fc=SCOL[skey])
    t.line([x0-1.5]*2, [ey0, ey0+eh], lw=0.5, color=AXC); t.line([x0, x0+ew], [ey0+eh+1.5]*2, lw=0.5, color=AXC)
    if E_OK:
        ps, pm = usage(pair, 'qualified_P_S', side), usage(pair, 'qualified_P_M', side)
        cur = {'S': ps[pprox][0], 'P': ps[pprox][1], 'M': pm[pprox][0]}
        dis = {'S': ps[pdist][0], 'P': ps[pdist][1], 'M': pm[pdist][0]}
        for st in ('S', 'M', 'P'):
            tot = cur[st]+dis[st]; frac = 100*dis[st]/tot
            xs = [XK(0), XK(0), XK(kb), XK(kb), XK(9.5)]; ys_ = [YK(100), YK(frac), YK(frac), YK(0), YK(0)]
            t.line(xs, ys_, lw=1.0, color=SC_[st], ls=LS_[st])
            EU[f'{SNAME[skey]}|{st}'] = dict(proximal=cur[st], distal=dis[st], distal_fraction_pct=frac, shared_mass=tot)
    else:
        t.rect(x0+1.0, ey0+1.0, ew-2.0, eh-2.0, fc='none', ec='#C9D1D6', lw=0.5, ls=(0, (2, 2)))
    for v in (0, 4, 8):
        t.line([XK(v)]*2, [ey0+eh+1.5, ey0+eh+3.5], lw=0.5, color=AXC); t.text_top(XK(v), ey0+eh+4.8, str(v), 7.5, ha='center')
    for v in (0, 50, 100):
        t.line([x0-3.5, x0-1.5], [YK(v)]*2, lw=0.5, color=AXC)
        if i == 0: t.text_top(x0-4.3, YK(v)-3.8, str(v), 7.5, ha='right')
if not E_OK:
    t.text_top(ex0+(3*ew+2*egap)/2, ey0+eh/2-4.0, 'awaiting homolog_site_pairs.tsv.gz', 7.0, ha='center', color=MUTED)
info('E', 'PAS-usage step curves (3 states x 3 studies, 4 vertices)', 9 if E_OK else 0, 4)
t.vtext(15.0, ey0+eh/2, 'Distal (%)', 7.5)
t.text_top(ex0+(3*ew+2*egap)/2, ey0+eh+14.0, 'Distance from proximal PAS (kb)', 7.5, ha='center')
x = 120.0
for st in ('S', 'P', 'M'):
    t.line([x, x+10.0], [335.6, 335.6], lw=1.0, color=SC_[st], ls=LS_[st]); t.text_top(x+12.5, 332.0, st, 7.5); x += 25.0
COLS_E = ['GSE118258', '178325_0618', '178325_0809', 'GSE150311', 'GSE221739', 'CIP', 'WOT']
cx0, cw = 40.0, 22.0
CX = lambda i: cx0+i*(cw+1.0)+(3.0 if i >= 5 else 0.0)
hb = 430.0
for i, route in enumerate(COLS_E):
    t.rect(CX(i), hb, cw, 4.5, fc=SCOL[route])
Ed = []
for j, (con, lab) in enumerate([('P_S', 'P–S'), ('P_M', 'P–M')]):
    ytop = 437.0+j*15.5
    t.text_top(12.0, ytop+3.5, lab, 8)
    for i, route in enumerate(COLS_E):
        x = CX(i); v = None; how = 'not in design'
        if route == '178325_0809':
            if con == 'P_S' and not hp_:
                how = 'pending: homolog_site_pairs.tsv.gz not in this bundle copy'
            elif con == 'P_S' and h3b:
                v = float(h3b[0]['native_h']); how = 'homolog_site_pairs H3b qualified_P_S native_h'
        elif (con, route) in eff:
            r = eff[(con, route)]
            if r['status'] == 'measured':
                v = float(r['delta']); how = 'member_endpoint_effects supported50_p08'
            else:
                how = r['status']
        if v is None:
            t.rect(x, ytop, cw, 14.5, fc=MISSING)
            if how.startswith('pending'):
                t.rect(x+0.4, ytop+0.4, cw-0.8, 13.7, fc='none', ec=MUTED, lw=0.5, ls=(0, (1.5, 1.5)))
        else:
            t.rect(x, ytop, cw, 14.5, fc=diverge(v, 0.3))
            t.text_top(x+cw/2, ytop+3.3, '−' if v < 0 else '+', 8, ha='center', color='white' if abs(v) > 0.2 else INK)
        Ed.append(dict(contrast=con, route=route, value=v, source=how))
R.claim('F4E', 'every available native ΔDWUI of the example is negative', all(e['value'] < 0 for e in Ed if e['value'] is not None))
R.axis('F4E', 'native ΔDWUI (colour; values below −0.3 use the extend triangle)', [e['value'] for e in Ed if e['value'] is not None], float('-inf'), 0.3)
info('E', 'native ΔDWUI tiles (study x contrast)', 14, 1)
cby = 475.5
t.text_top(12.0, cby-2.4, 'Native ΔDWUI', 7.5)
cbx0 = 12.0+textwidth('Native ΔDWUI', 7.5)+12.0; cbx1 = cbx0+56.0
for i in range(80):
    v = -0.3+0.6*(i+0.5)/80
    t.rect(cbx0+i*(cbx1-cbx0)/80, cby, (cbx1-cbx0)/80+0.15, 4.0, fc=diverge(v, 0.3))
from matplotlib.patches import Polygon
if min([e['value'] for e in Ed if e['value'] is not None] or [0]) < -0.3:
    t.ax.add_patch(Polygon([(cbx0-0.1, H-cby), (cbx0-0.1, H-cby-4.0), (cbx0-3.6, H-cby-2.0)], closed=True, fc=DIV[0], ec='none'))
for v, lab in [(-0.3, '−.3'), (0, '0'), (0.3, '+.3')]:
    x = cbx0+(v+0.3)/0.6*(cbx1-cbx0); t.text_top(x, cby+5.6, lab, 7.5, ha='center')
xm = cbx1+14.0
t.rect(xm, cby, 7.0, 4.0, fc=MISSING); t.text_top(xm+9.5, cby-2.4, 'Not available', 7.5)
assert xm+9.5+textwidth('Not available', 7.5) < 236.0
key_rows(12.0, cby+19.0, ['GSE118258', '178325_0618', '178325_0809', 'GSE150311', 'GSE221739', 'CIP', 'WOT'], xmax=236.0)
ver['E'] = dict(gene=specE.get('gene_symbol'), human_gene_id=HG, mouse_gene_id=MG, curves_drawn=E_OK, usage=EU,
                human_sites=[hp, hd], mouse_sites=[mp, md], human_kb=hkb, mouse_kb=mkb, cells=Ed)

# ---------------------------------------------------------------- save
t.save(R.out/'overlay.pdf'); overlay(pdf, R.out/'overlay.pdf', 0, 0); pdf.save(OUT)
render_png(OUT, R.out/'Figure.4.draft.png', 150)
ver['output_sha256'] = sha(OUT)
R.finish(4, OUT, json.loads(json.dumps(ver, default=str)))
