"""Supplementary Figures D01 (human) and D02 (mouse) — frozen expression states: marker evidence, sampling support and
programme scores in one state order (printed numbers from supdraw.SNUM). Expression only; unchanged by the PAS reanalysis.

a  frozen marker panel across every state (columns, grouped by study/route unit): colour = relative expression (0–1 per
   gene), area = detected-cell fraction; ring = small-stratum fallback (no library with >= 20 cells)
b  sampling support of the same states: libraries with >= 20 cells, and the share of the largest library
c  native ECM, epithelial and pluripotency programme scores of every qualified library x state observation (>= 20 cells)
Usage: D_markers.py human|mouse
"""
import sys, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

SP = sys.argv[1]
FID = {'human': 'D01', 'mouse': 'D02'}[SP]
PA = FID+'A'
D = Data(); OUT = outdir(FID)
ST = D.df(PA, 'frozen/states.tsv'); ST = ST[ST.species == SP].reset_index(drop=True)
MB = D.df(PA, 'frozen/marker_by_state.tsv'); MB = MB[MB.species == SP].set_index(['state_key', 'gene'])
SEL = json.loads(D.path(PA, 'frozen/marker_selection.json').read_text())['markers'][SP]
PAL = json.loads(D.path(PA, 'frozen/palette.json').read_text())['state_groups']
PG = D.df(FID+'C', 'frozen/program_by_state_library.tsv.gz'); PG = PG[(PG.species == SP) & (PG.n_cells >= 20)]
CODE = dict(shared_s='S', shared_r='R', shared_e='E', shared_tp='P', shared_m='M', shared_k='K', human_n='hN', human_x='hX',
            human_u='hU', mouse_b='mB', mouse_n='mN', mouse_i='mI', mouse_y='mY', mouse_u='mU')
NAME = dict(shared_s='source', shared_r='remodelling', shared_e='epithelial', shared_tp='pluripotency', shared_m='mesenchymal',
            shared_k='keratin', human_n='naive-related', human_x='endoderm-secretory', human_u='unresolved', mouse_b='trophoblast-like',
            mouse_n='neural', mouse_i='IFN-high', mouse_y='myeloid', mouse_u='unresolved')
UNIT = {'GSE118258': ('Xing2020', STUDY['GSE118258'][1]), 'GSE150311': ('Liu2020', STUDY['GSE150311'][1]),
        'GSE178325': ('Guan2022', DATASET['GSE178325'][1]), 'GSE221739': ('Panariello2023', STUDY['GSE221739'][1]),
        'gse103221_osk': ('Guo2019–OSK', STUDY['OSK'][1]), 'gse103221_cip': ('Guo2019–CIP', STUDY['CIP'][1]),
        'gse171127_d3_getm': ('Jaber D3 GETM', STUDY['171127_GETM'][1]), 'gse171127_d3_oskm': ('D3 OSKM', STUDY['171127_OSKM'][1]),
        'gse171127_d6_pooled': ('D6', DATASET['GSE171127'][1]), 'gse171127_d12_pooled': ('D12', DATASET['GSE171127'][1]),
        'wot_common_to_2i': ('Schiebinger2019–2i', STUDY['WOT'][1])}
GENES = [(grp, g) for grp, gs in SEL.items() for g in gs]
NS, NG = len(ST), len(GENES)
assert NS == (62 if SP == 'human' else 49) and NG == (31 if SP == 'human' else 30)

X0, RR = MARGIN, PAGE_W-MARGIN
gw = max(textwidth(g, 7.0, 'RegIta') for _, g in GENES)
grw = max(textwidth(grp, 6.8) for grp in SEL)
CX0 = X0+grw+4.0+gw+4.0
CW = (RR-1.0-CX0)/NS
RP = 7.0
HDR_A = 34.0
HA = HDR_A+NG*RP+3.0+12.0+12.0+10.0+5.0
HB, HC = 62.0, 96.0
Y1 = MARGIN; Y2 = Y1+HA+GAP; Y3 = Y2+HB+GAP; H = Y3+HC+MARGIN
F = Fig(H, figname(FID))
cx = lambda i: CX0+CW*(i+0.5)
units = list(dict.fromkeys(ST.unit))
bounds = []
for u in units[:-1]:
    bounds.append(CX0+CW*(ST.index[ST.unit == u].max()+1))


def unit_rules(ytop, ybot):
    for b in bounds:
        F.line([b, b], [ytop, ybot], lw='ref', color=AXC)


# ================================================================ a  marker matrix
x, y, w, h = F.panel('a', X0, Y1, RR-X0, HA)
HEAD = []
for u in units:
    key = 'GSE171127' if u.startswith('gse171127') else u
    if HEAD and HEAD[-1][0] == key: HEAD[-1][1].append(u)
    else: HEAD.append((key, [u]))
for key, us in HEAD:
    idx = ST.index[ST.unit.isin(us)]
    name, col = (('Jaber2022', DATASET['GSE171127'][1]) if key == 'GSE171127' else UNIT[key])
    xm = (cx(idx.min())+cx(idx.max()))/2; span = CW*len(idx)
    lab = name if textwidth(name, 7.0)+7.6 < span+4.0 else name.split('–')[-1].replace('Jaber ', '')
    tw = textwidth(lab, 7.0)
    x0u = min(max(xm-(tw+7.6)/2, CX0), RR-1.0-(tw+7.6))
    F.rect(x0u, y+4.2, 5.6, 5.6, fc=col)
    F.text(x0u+7.6, y+3.6, lab, 7.0)
yf = y+14.0
for i, r in ST.iterrows():
    F.rect(CX0+CW*i+0.2, yf, CW-0.4, 7.0, fc=PAL[r.family])
    code = CODE[r.family]
    F.text(cx(i), yf+7.8, code[-1] if len(code) == 2 else code, 6.8, ha='center', color=INK)
y0 = y+HDR_A
unit_rules(yf, y0+NG*RP)
dmax = min(CW, RP)-0.7
prev = None
V = {}
for k, (grp, g) in enumerate(GENES):
    yy = y0+RP*(k+0.5)
    if grp != prev:
        F.text(X0+0.4, yy-3.3, grp, 6.8, color=MUTED); prev = grp
        if k: F.line([X0+grw+4.0, RR-1.0], [y0+RP*k]*2, lw='ref', color=BAND)
    F.text(CX0-3.0, yy-3.4, g, 7.0, font='RegIta', ha='right')
    for i, r in ST.iterrows():
        v = MB.loc[(r.state_key, g)]
        F.dot(cx(i), yy, float(v.detected_fraction), seq(float(v.relative_expression)), dmax)
        if bool(v.fallback_small):
            F.marker(cx(i), yy, 'o', 'mark', fc='none', ec=AXC)
        V[f'{r.state_key}|{g}'] = [round(float(v.relative_expression), 4), round(float(v.detected_fraction), 4), bool(v.fallback_small)]
yb = y0+NG*RP+2.0
for i in range(NS):
    if (i+1) % 5 == 0 or i == 0:
        F.text(cx(i), yb, str(i+1), 6.8, ha='center', color=MUTED)
F.text(CX0-3.0, yb, 'State', 6.8, ha='right', color=MUTED)
ky = yb+11.0
parts_x = CX0+(RR-CX0)/2
F.cbar_key(CX0+(RR-CX0)*0.24, ky, 'Relative expression', '0', '1', lambda f_: seq(f_), wbar=50)
xd = parts_x-6.0
xd += F.text(xd, ky, 'Detected cells', 7.5)+6.0
for fr, lab in ((0.1, '10%'), (0.5, '50%'), (1.0, '100%')):
    F.dot(xd+dmax/2, ky+3.8, fr, seq(0.6), dmax); xd += dmax+2.0
    xd += F.text(xd, ky, lab, 7.0)+6.0
F.marker(xd+2.8, ky+3.8, 'o', 'mark', fc='none', ec=AXC); F.text(xd+7.6, ky, 'fallback', 7.0)
fams = [f for f in CODE if f in set(ST.family)]
sh = [f for f in fams if f.startswith('shared')]; lo = [f for f in fams if not f.startswith('shared')]
F.key_row(None, ky+12.0, [('block', PAL[f], f'{CODE[f]} {NAME[f]}') for f in sh], size=7.0, gap=6.0, center=(CX0+RR)/2)
F.key_row(None, ky+22.0, [('block', PAL[f], f'{CODE[f]} {NAME[f]}') for f in lo], size=7.0, gap=6.0, center=(CX0+RR)/2)
F.values['a'] = dict(states=ST.state_key.tolist(), genes=[g for _, g in GENES], values=V)

# ================================================================ b  sampling support
x, y, w, h = F.panel('b', X0, Y2, RR-X0, HB)
nmax = int(ST.n_libraries_ge20.max())
sy1 = Lin(0, nmax, y+27.0, y+4.0)
sy2 = Lin(0, 1, y+HB-4.0, y+35.0)
F.text(X0+16.0, (sy1.p0+sy1.p1)/2+3.0, 'Libraries', 7.0); F.text(X0+16.0, (sy1.p0+sy1.p1)/2+10.0, '≥ 20 cells', 7.0, color=MUTED)
F.text(X0+0.4, (sy2.p0+sy2.p1)/2-6.0, 'Largest', 7.0); F.text(X0+0.4, (sy2.p0+sy2.p1)/2+1.0, 'library share', 7.0, color=MUTED)
F.yaxis(sy1, CX0, [0, nmax], ['0', str(nmax)])
F.yaxis(sy2, CX0, [0, 1], ['0', '1'])
for i, r in ST.iterrows():
    n = r.n_libraries_ge20
    if n == n and n > 0:
        F.rect(CX0+CW*i+0.6, sy1(n), CW-1.2, sy1(0)-sy1(n), fc=PAL[r.family])
    F.rect(CX0+CW*i+0.6, sy2(r.max_library_fraction), CW-1.2, sy2(0)-sy2(r.max_library_fraction), fc=MUTED)
unit_rules(y+2.0, y+HB-2.0)
F.values['b'] = {r.state_key: dict(libraries_ge20=None if r.n_libraries_ge20 != r.n_libraries_ge20 else int(r.n_libraries_ge20),
                                   max_library_fraction=float(r.max_library_fraction)) for r in ST.itertuples()}

# ================================================================ c  programme scores
x, y, w, h = F.panel('c', X0, Y3, RR-X0, HC)
PROG = ['ECM', 'Epithelial', 'Pluripotency']
top = max(float(PG.mean_log1p_cp10k.max()), 0.5)
lim = float(np.ceil(top*2)/2)
sub = (HC-4.0-2*5.0-14.0)/3
Cv = {}
idx = {k: i for i, k in enumerate(ST.state_key)}
assert (PG.family == PG.state_key.map(ST.set_index('state_key').family)).all()
for j, prog in enumerate(PROG):
    sy = Lin(0, lim, y+4.0+j*(sub+5.0)+sub, y+4.0+j*(sub+5.0))
    F.text(X0+(16.0 if j == 0 else 0.4), (sy.p0+sy.p1)/2-3.4, prog, 7.0)
    F.yaxis(sy, CX0, [0, lim], ['0' if j == 2 else '', f'{lim:g}'])
    q = PG[PG.program == prog]
    for sk, g in q.groupby('state_key'):
        i = idx[sk]; n = len(g)
        offs = np.linspace(-CW*0.28, CW*0.28, n) if n > 1 else [0.0]
        F.points([cx(i)+o for o in offs], [sy.clip(v) for v in g.mean_log1p_cp10k], ms='point', color=PAL[g.family.iloc[0]], alpha='full')
    Cv[prog] = dict(observations=int(len(q)), states=int(q.state_key.nunique()))
unit_rules(y+2.0, y+HC-14.0)
F.text((CX0+RR)/2, y+HC-7.4, 'Native mean log1p(CP10k) of each qualified library × state (≥ 20 cells); color = display family', 7.0, ha='center', color=MUTED)
F.values['c'] = Cv
F.notes.append('frozen upstream annotation and native expression; columns in the state order of frozen/states.tsv (state key in values); '
               'vertical rules separate study/route units (Jaber2022: D3 GETM, D3 OSKM, D6 pooled, D12 pooled)')
F.save(OUT, D)
