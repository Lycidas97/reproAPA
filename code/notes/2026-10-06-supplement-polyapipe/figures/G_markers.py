"""Supplementary Figures G1 (human) and G2 (mouse) — native marker evidence for every frozen expression state, by study,
with the stability of each state under nine neighbouring Leiden partitions (printed numbers from supdraw.SNUM).

Re-laid-out from the 2026-09-14 state-grouping supplement (SMC1/SMC2, notes/2026-09-14-state-grouping-robustness/
submission-v3/render.py): same states, same gene panels (10 core markers per species + the study-local markers), same J/F
values. Changes: house layout; families and review flags from the current frozen annotation in the supplement bundle
(frozen/states.tsv; two GSE171127 D3 states are now mU, previously mY); marker values use the bundle's library weighting
(equal weight over libraries with >= 20 cells in the state, all libraries when none has 20 — identical to the bundle's
frozen/marker_by_state.tsv for the 61 genes it carries, checked below).
Usage: G_markers.py human|mouse
"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

SP = sys.argv[1]
FID = {'human': 'G1', 'mouse': 'G2'}[SP]
D = Data(); X = ExtData(panel=FID+'A'); OUT = outdir(FID)
ST = D.df(FID+'A', 'frozen/states.tsv')
ST = ST[ST.species == SP].copy()
BM = D.df(FID+'A', 'frozen/marker_by_state.tsv')
PAL = __import__('json').loads(D.path(FID+'A', 'frozen/palette.json').read_text())['state_groups']
M = X.df('marker_by_state')
PS = X.df('state_partition_summary').set_index('state_id')
LIB = X.df('state_library_inventory')

# ---- display values: the bundle's weighting rule, reproduced for the study-local genes and checked on the bundle genes
M['mean_disp'] = M.mean_log1p_cp10k_ge20.where(M.mean_log1p_cp10k_ge20.notna(), M.mean_log1p_cp10k)
M['det_disp'] = M.detected_fraction_ge20.where(M.detected_fraction_ge20.notna(), M.detected_fraction)
chk = BM.merge(M, left_on=['state_key', 'gene'], right_on=['state_id', 'gene'], suffixes=('_b', ''))
assert len(chk) == len(BM[BM.species == SP]) or len(chk) == len(BM)
dmax_mean = float(np.abs(chk.mean_log1p_cp10k_b-chk.mean_disp).max()); dmax_det = float(np.abs(chk.detected_fraction_b-chk.det_disp).max())
assert dmax_mean < 1e-9 and dmax_det < 1e-9, (dmax_mean, dmax_det)
M = M.set_index(['state_id', 'gene'])

CODE = dict(shared_s='S', shared_r='R', shared_e='E', shared_tp='P', shared_m='M', shared_k='K', human_n='hN', human_x='hX',
            human_u='hU', mouse_b='mB', mouse_n='mN', mouse_i='mI', mouse_y='mY', mouse_u='mU')
FORDER = {f: i for i, f in enumerate(['shared_s', 'shared_r', 'shared_e', 'shared_tp', 'shared_m', 'shared_k', 'human_n',
                                      'human_x', 'human_u', 'mouse_b', 'mouse_n', 'mouse_i', 'mouse_y', 'mouse_u'])}
CORE = {'human': ['DCN', 'COL1A1', 'EPCAM', 'KRT8', 'KRT18', 'DPPA4', 'ESRG', 'TOP2A', 'FOS', 'ISG15'],
        'mouse': ['Dcn', 'Col1a1', 'Epcam', 'Krt8', 'Krt18', 'Dppa4', 'Sall4', 'Top2a', 'Fos', 'Isg15']}
EXTRA = {'GSE118258': ['CD248', 'LOX', 'PRRX2', 'CAV1', 'H19', 'MFAP5'], 'GSE150311': ['DNMT3L', 'DPPA3', 'KRT17', 'TUBB3'],
         'GSE178325': ['MMP3', 'KRT19', 'RARRES1', 'RBP1', 'TTR', 'APOA1'], 'GSE221739': ['FN1', 'ACTA2'],
         'wot_common_to_2i': ['Timp1', 'Pclaf', 'Dmrtc2', 'Aldh3a1', 'Tubb3', 'Elavl3', 'Peg10', 'Phlda2'],
         'gse103221_cip': ['Pcolce', 'Tgfbi', 'Flt1', 'Pdpn', 'Cldn7', 'Cldn4', 'Sprr1a'],
         'gse103221_osk': ['Fst', 'Fbln5', 'Tm4sf1', 'Ube2c', 'Cdca8', 'Cenpf'],
         'GSE171127': ['Cldn4', 'Sprr1a', 'Tubb3', 'Elavl3', 'Lyz2', 'Msr1', 'Peg10', 'Phlda2']}
BOUNDARY = {'GSE118258::studywide::C005', 'gse103221_cip|C012', 'wot_common_to_2i|WOT::common_to_2i::L00'}
ALIAS = {'gse171127_d3_getm': 'G3', 'gse171127_d3_oskm': 'O3', 'gse171127_d6_pooled': 'D6', 'gse171127_d12_pooled': 'D12'}
UORDER = {u: i for i, u in enumerate(ALIAS)}
if SP == 'human':
    PANELS = [('a', 'GSE118258', 'Xing2020', '#753399'), ('b', 'GSE150311', 'Liu2020', '#648290'),
              ('c', 'GSE178325', 'Guan2022', '#A82C60'), ('d', 'GSE221739', 'Panariello2023', '#30363D')]
else:
    PANELS = [('a', 'wot_common_to_2i', 'Schiebinger2019–2i', '#207450'), ('b', 'gse103221_osk', 'Guo2019–OSK', '#255FBC'),
              ('c', 'gse103221_cip', 'Guo2019–CIP', '#CF4C00'), ('d', 'GSE171127', 'Jaber2022', '#766A8A')]
LOWJ = 0.80
SEQ_MAX = 5.0


def rows_of(key):
    q = ST[ST.study_id == key] if key.startswith('GSE') else ST[ST.unit == key]
    q = q.assign(fo=q.family.map(FORDER), uo=q.unit.map(UORDER).fillna(0))
    return q.sort_values(['uo', 'fo', 'short_id'])


R = {L: rows_of(k) for L, k, *_ in PANELS}
GENES = {L: CORE[SP]+EXTRA[k] for L, k, *_ in PANELS}
LMAX = max(textwidth(g, 7.0, 'RegIta') for gs in GENES.values() for g in gs)
HDR = 10.6+LMAX+2.0                    # title row + rotated gene labels + gap
LEGH = 21.0                            # two legend lines under the bottom panels
X0, RR = MARGIN, PAGE_W-MARGIN
WL = (RR-X0-GAP)/2; XR = X0+WL+GAP
nL = len(R['a'])+len(R['c']); nR = len(R['b'])+len(R['d'])
PITCH_LONG = 8.0
nmax = max(nL, nR)
# 2026-10-12: J/F and bold/flag notes moved to the caption, so only the left column (size/colour key) reserves LEGH
TOT = 2*HDR+GAP+0.8+max(nL*PITCH_LONG+LEGH, nR*PITCH_LONG)
pitch = {'left': (TOT-2*HDR-GAP-LEGH-0.8)/nL, 'right': (TOT-2*HDR-GAP-0.8)/nR}
H = MARGIN*2+TOT
F = Fig(H, figname(FID))
VALS = {}


def panel(L, key, name, col, x, y, w, pitch_, legend_lines):
    rows = R[L]; genes = GENES[L]; n = len(rows)
    h = HDR+n*pitch_+0.4+(LEGH if legend_lines else 0)
    F.panel(L, x, y, w, h)
    tx = x+0.2+textwidth(L, 11, 'Bol')+F.CLEAR_R+0.6
    F.rect(tx, y+3.2, 5.6, 5.6, fc=col)
    tw = F.text(tx+7.6, y+2.4, name, 8.0)
    F.text(tx+7.6+tw+3.0, y+3.1, key.split('_')[0].upper() if not key.startswith('GSE') else key, 7.0, color=MUTED) \
        if not key.startswith('wot') else F.text(tx+7.6+tw+3.0, y+3.1, 'GSE115943', 7.0, color=MUTED)
    # row labels
    labs = []
    for r in rows.itertuples():
        s = (ALIAS[r.unit]+' ' if r.unit in ALIAS else '')+r.short_id
        sparse = not (LIB[LIB.state_id == r.state_key].n_cells >= 20).any()
        labs.append((s+('†' if sparse else ''), r, sparse))
    lw_ = max(textwidth(s, 7.0, 'Bol') for s, *_ in labs)
    xl = x+0.4+lw_; xf = xl+8.2
    jw = 14.0; xj = x+w-2*jw+1.6
    gx0 = xf+7.5; cw = (xj-1.0-gx0)/len(genes)
    dmax = min(cw, pitch_, F.DOT_MAX)-0.9
    y0 = y+HDR
    F.text(xj+jw/2, y0-9.2, 'J', 7.5, ha='center'); F.text(xj+1.5*jw, y0-9.2, 'F', 7.5, ha='center')
    for i, g in enumerate(genes):
        gxc = gx0+cw*(i+0.5)
        F.vlabel(gxc, y0-2.0, g, 7.0, font='RegIta')
    F.line([gx0+cw*len(CORE[SP])]*2, [y0-0.6, y0+n*pitch_-0.2], lw='ref', color=AXC, ls='ref')
    VALS[L] = {}
    for j, (s, r, sparse) in enumerate(labs):
        yc = y0+pitch_*(j+0.5)
        bold = r.state_key in BOUNDARY
        F.text(xl, yc-3.4, s, 7.0, font='Bol' if bold else 'Reg', ha='right')
        code = CODE[r.family]+('*' if r.review_flag in ('M', 'C') else '')
        F.text(xl+2.2, yc-3.4, code, 7.0, color=PAL[r.family])
        for i, g in enumerate(genes):
            v = M.loc[(r.state_key, g)]
            assert bool(v.all_libraries_feature_available)
            F.dot(gx0+cw*(i+0.5), yc, float(v.det_disp), seq(float(v.mean_disp), 0, SEQ_MAX), dmax)
        ps = PS.loc[r.state_key]; assert int(ps.runs) == 9
        for k, val in enumerate([float(ps.median_jaccard), float(ps.median_family_neighbor_fraction)]):
            F.text(xj+jw*(k+0.5), yc-3.3, f'{val:.2f}'.replace('0.', '.').replace('1.00', '1'), 6.8, ha='center',
                   font='Bol' if val < LOWJ else 'Reg', color=INK if val < LOWJ else MUTED)
        VALS[L][r.state_key] = dict(label=s, family=CODE[r.family], review_flag=None if r.review_flag != r.review_flag else r.review_flag,
                                    sparse=bool(sparse), J=float(ps.median_jaccard), F=float(ps.median_family_neighbor_fraction),
                                    markers={g: [round(float(M.loc[(r.state_key, g)].mean_disp), 4), round(float(M.loc[(r.state_key, g)].det_disp), 4)] for g in genes})
    if legend_lines:
        yl = y0+n*pitch_+0.4+3.6
        legend_lines(x, w, gx0, xj, yl)
    return y+h


def centred(xc, yl, parts):
    """parts: list of (kind, payload, width); drawn left to right, centred on xc."""
    tot = sum(w_ for *_, w_ in parts); xs = xc-tot/2
    for kind, pay, w_ in parts:
        if kind == 'text': F.text(xs, yl, pay[0], pay[1], color=pay[2] if len(pay) > 2 else INK)
        elif kind == 'dot': F.dot(xs+3.6, yl+3.8, pay, seq(3.5, 0, SEQ_MAX), F.DOT_MAX-0.9)
        elif kind == 'cbar': F.cbar(xs+textwidth('0', 7.0)+1.8, yl+1.4, pay, '0', '5', lambda f: seq(f*SEQ_MAX, 0, SEQ_MAX))
        xs += w_


def leg_size_colour(x, w, gx0, xj, yl):
    xc = (gx0+xj)/2
    parts = [('text', ('Detected cells', 7.5), textwidth('Detected cells', 7.5)+5)]
    for fr, lab in [(0.1, '10%'), (0.5, '50%'), (1.0, '100%')]:
        parts += [('dot', fr, 8.6), ('text', (lab, 7.0), textwidth(lab, 7.0)+6)]
    centred(xc, yl, parts[:-1]+[('text', ('100%', 7.0), textwidth('100%', 7.0))])
    cbw = 50.0
    F.cbar_key(xc, yl+10.0, 'Mean log1p(CP10k)', '0', '5', lambda f: seq(f*SEQ_MAX, 0, SEQ_MAX), wbar=cbw)





yb = panel('a', PANELS[0][1], PANELS[0][2], PANELS[0][3], X0, MARGIN, WL, pitch['left'], None)
panel('c', PANELS[2][1], PANELS[2][2], PANELS[2][3], X0, yb+GAP, WL, pitch['left'], leg_size_colour)
yb = panel('b', PANELS[1][1], PANELS[1][2], PANELS[1][3], XR, MARGIN, WL, pitch['right'], None)
panel('d', PANELS[3][1], PANELS[3][2], PANELS[3][3], XR, yb+GAP, WL, pitch['right'], None)   # J/F/bold/flag notes: caption
F.values = VALS
F.notes.append(f'marker display rule reproduced and checked against bundle frozen/marker_by_state.tsv: max |Δ| mean {dmax_mean:.1e}, detected {dmax_det:.1e}')
F.notes.append('row pitch differs between the two columns so that both columns end at the same height: '
               + f"{pitch['left']:.2f} / {pitch['right']:.2f} pt")
F.save(OUT, [D, X])
