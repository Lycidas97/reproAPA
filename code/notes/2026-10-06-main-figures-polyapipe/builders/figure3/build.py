"""Figure 3 · records adapter of the v5.1 production builder (layout unchanged; data from one result bundle).

Usage: build.py --records <bundle>/result-records.json --out <new dir>
Panel notes below describe the production figure drawn from the archived v2 data; their numbers are history, not checks.
Example genes come from display-selection.json (pinned or rank rule) and are validated against the bundle.

Base: locked master Fig/Figure.3.pdf. Contract: main-panel-replacement-v5.tsv rows 3A–3E.
  A  four candidate-recurrence sets at |ΔDWUI| ≥ 0.02, series_matching = True; fold = observed / null mean,
     counts printed against the null median; Studies column kept.
  B  set-level controls (P–S, 0.02; six-rule row removed) + retained members (bootstrap, independent wells)
     + same-state pseudo vs observed medians (reference S and target P halves, named separately).
  C  dual-reference display subset + all-gene scatters; human 223→199→112 (Xing only), mouse 257→235→128.
  D  Pdlim7: state DWUI (matched P/M/S contributors) and PAS distal-fraction curves. Metacell quantiles are
     NOT available locally (export needed on the server) — state means only, flagged in verification.
  E  APA–GEX: same-direction APA with opposite GEX, identity and fixed-responder direction tests by species.
"""
import json, sys, math
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools-local'))
from v5draw import *  # noqa
from figrec import Run, MASTER_F3, pas_pos, pas_strand
HERE = Path(__file__).resolve().parent
R = Run('F3')
import numpy as np

# Fig/Figure.3.pdf is now the promoted v5 final; the locked master is archived here (same SHA-256).
MASTER = MASTER_F3
assert sha(MASTER) == '0c51339a2a124cfddda4c7aea6c1f666bedcc0ad23d058c01ab7d373ed5d5010'
OUT = R.out/'Figure.3.draft.pdf'
use = R.rows
f2 = lambda s: abs(float(s)-0.02) < 1e-9

pl = pdfplumber.open(MASTER).pages[0]; W, H = pl.width, pl.height
pdf = pikepdf.open(MASTER)
CLEAR = {
    'A': [(34, 15, 236.5, 129), (4, 129, 236.5, 145.5)],
    'B': [(246.5, 0, 454, 146.2), (236.5, 18, 246.5, 146.2)],
    'C': [(14, 146.2, 454, 437.4), (0, 165, 14, 437.4)],
    'D': [(18, 436, 233.5, 641), (0, 451, 18, 641)],
    'E': [(243, 436, 454, 641), (233.5, 451, 243, 641)],
}
removed = {k: len(clear(pdf, H, v)) for k, v in CLEAR.items()}
t = Top(W, H); ver = {'master_sha256': sha(MASTER), 'removed_objects': removed}
STUDY = {'CIP': ('#CF4C00', 's'), 'WOT': ('#207450', 'o'), 'GSE118258': ('#753399', 'o'),
         'GSE150311': ('#648290', 'D'), 'GSE221739': ('#30363D', '^')}
SC = dict(S='#E69F00', R='#0072B2', E='#009E73', P='#CC79A7', M='#8C6D31')

# ------------------------------------------------------------------ A
rec = use('F3A', 'recurrence/recurrence_summary.tsv')
dirc = use('F3A', 'panel_statistics/recurrence-direction-counts.tsv', lambda r: r['metric'] == 'WPP' and f2(r['threshold']))
clipA = lambda v: min(max(v, 0.6), 3.4)
XA = lambda v: 58.0+(v-0.6)*20.0
f5 = lambda s: abs(float(s)-0.05) < 1e-9
t.text_top(88.0, 18.5, 'Obs./null', 8, ha='center')
t.text_top(130.5, 18.5, 'Obs./median', 8, ha='center'); _hd = [('Dist.', '#D55E00'), (' / ', INK), ('prox.', '#0072B2')]; _hx = 174.0-sum(textwidth(a_, 7.5) for a_, _ in _hd)/2
for a_, c_ in _hd:
    t.text_top(_hx, 19.0, a_, 7.5, color=c_); _hx += textwidth(a_, 7.5)
t.text_top(213.5, 18.5, 'Studies', 8, ha='center')
t.line([XA(1), XA(1)], [30.7, 114.8], lw=0.5, color='#98A1A8', ls=(0, (1, 1.5)))
t.line([XA(0.6), XA(3.4)], [114.8, 114.8], lw=0.5, color=INK)
for v in (1, 2, 3):
    t.line([XA(v)]*2, [114.8, 116.8], lw=0.5, color=INK); t.text_top(XA(v), 118.8, str(v), 8, ha='center')
Ad = []
for (sp, con), y in zip([('human', 'P_S'), ('mouse', 'P_S'), ('human', 'P_M'), ('mouse', 'P_M')], [41.2, 62.9, 84.6, 106.2]):
    r = [r for r in rec if r['metric'] == 'WPP' and r['species'] == sp and r['contrast'] == con and f2(r['threshold'])
         and r['series_matching'] == 'True' and r['variant'] == 'primary']
    assert len(r) == 1; r = r[0]
    obs = int(r['observed']); mean = float(r['null_mean']); med = float(r['null_median'])
    lo, hi = float(r['null_p025'])/mean, float(r['null_p975'])/mean; fold = float(r['enrichment'])
    if abs(fold-obs/mean) > 1e-9: R.flag('F3A', 'consistency', f'{sp} {con}: enrichment != observed/null mean')
    if not R.axis('F3A', f'{sp} {con} obs/null (0.02)', [lo, hi, fold], 0.6, 3.4): lo, hi, fold = clipA(lo), clipA(hi), clipA(fold)
    t.line([XA(lo), XA(hi)], [y-3.0, y-3.0], lw=1.7, color='#C9CFD3')
    t.marker(XA(fold), y-3.0, 'o', 4.0, mew=0.3)
    t.text_top(130.5, y-8.2, f'{obs}/{int(round(med))}', 8, ha='center')
    r5 = [r_ for r_ in rec if r_['metric'] == 'WPP' and r_['species'] == sp and r_['contrast'] == con and f5(r_['threshold'])
          and r_['series_matching'] == 'True' and r_['variant'] == 'primary']
    assert len(r5) == 1; r5 = r5[0]
    if r5['n_studies'] != r['n_studies']: R.flag('F3A', 'consistency', f'{sp} {con}: n_studies differs between thresholds')
    obs5, mean5, med5 = int(r5['observed']), float(r5['null_mean']), float(r5['null_median'])
    lo5, hi5, fold5 = float(r5['null_p025'])/mean5, float(r5['null_p975'])/mean5, float(r5['enrichment'])
    if not R.axis('F3A', f'{sp} {con} obs/null (0.05)', [lo5, hi5, fold5], 0.6, 3.4): lo5, hi5, fold5 = clipA(lo5), clipA(hi5), clipA(fold5)
    if XA(fold5)+2.5 >= 130.5-textwidth(f'{obs}/{int(round(med))}', 8)/2: R.flag('F3A', 'overlap', f'{sp} {con}: 0.05 marker reaches the count column')
    t.line([XA(lo5), XA(hi5)], [y+3.5, y+3.5], lw=1.7, color='#C9CFD3')
    t.marker(XA(fold5), y+3.5, 'D', 3.6, mew=0.3)
    t.text_top(130.5, y+0.6, f'{obs5}/{int(round(med5))}', 8, ha='center')
    dc = {m['direction']: int(m['genes']) for m in dirc if m['species'] == sp and m['contrast'] == con}
    dist, prox = dc.get('positive', 0), dc.get('negative', 0)
    if dist+prox != obs: R.flag('F3A', 'consistency', f'{sp} {con}: direction counts {dist}+{prox} != observed {obs}')
    tot_ = max(dist+prox, 1)
    bx0_, bw_ = 156.0, 36.0
    t.rect(bx0_, y-5.0, bw_*dist/tot_, 4.4, fc='#D55E00'); t.rect(bx0_+bw_*dist/tot_, y-5.0, bw_*prox/tot_, 4.4, fc='#0072B2')
    t.text_top(bx0_+bw_/2, y+0.8, f'{dist} | {prox}', 7.5, ha='center')
    t.text_top(213.5, y-6.5, r['n_studies'], 8, ha='center')
    ex3 = r.get('exact_3_observed', '')
    if ex3 not in ('', 'nan'):
        pass                                     # exact-3 count -> caption
    Ad.append(dict(species=sp, contrast=con, observed=obs, null_mean=mean, null_median=med, fold=fold, null_lo=lo, null_hi=hi,
                   n_studies=int(r['n_studies']), p=float(r['empirical_p']), distal=dist, proximal=prox,
                   exact_3=None if ex3 in ('', 'nan') else int(float(ex3)),
                   t005=dict(observed=obs5, null_mean=mean5, null_median=med5, fold=fold5, null_lo=lo5, null_hi=hi5, p=float(r5['empirical_p']))))
_x = 6.0; t.text_top(_x, 131.2, '|ΔDWUI| ≥', 8); _x += textwidth('|ΔDWUI| ≥', 8)+5.0
t.marker(_x, 135.0, 'o', 4.0, mew=0.3); t.text_top(_x+4.5, 131.2, '0.02', 8); _x += 4.5+textwidth('0.02', 8)+7.0
t.marker(_x, 135.0, 'D', 3.6, mew=0.3); t.text_top(_x+4.5, 131.2, '0.05', 8); _x += 4.5+textwidth('0.05', 8)+12.0
t.line([_x, _x+8.7], [135.0, 135.0], lw=1.7, color='#C9CFD3'); t.text_top(_x+12.0, 131.2, '95% perm. range', 8)
assert _x+12.0+textwidth('95% perm. range', 8) < 231.0, _x
ver['A'] = Ad

# ------------------------------------------------------------------ B
mn = use('F3B', 'measurement/nulls/recurrence_summary.tsv')
lad = use('F3B', 'recurrence_controls/recurrence_ladder.tsv')
ret = use('F3B', 'measurement/nulls/member_retention.tsv')
ss = use('F3B', 'measurement/same-state/recurrence_summary.tsv')
well = {r['background']: r for r in use('F3B', 'panel_statistics/well-denominators.tsv') if r['metric'] == 'DWUI'}
clipB = lambda v: min(max(v, 0.75), 2.0)
XB = lambda v: 339.4+(v-1)*107.7
CTRL = [('interval_supported', 'Bootstrap interval'), ('depth20', '20 PAS UMIs'), ('cycle_matched', 'Cycle-score matched'),
        ('common', 'Common PAS axis'), ('L3_full', 'Nested backgrounds')]
Bd = {'controls': []}
for i, (var, lab) in enumerate(CTRL):
    yh = 24.1+9.1*i; ym = yh+3.8
    t.text_top(312.2, yh-1.0, lab, 7.5, ha='right')
    for sp, y in [('human', yh), ('mouse', ym)]:
        if var == 'L3_full':
            r = [r for r in lad if r['species'] == sp and r['contrast'] == 'P_S' and f2(r['threshold']) and r['model'] == 'L3_full'][0]
            obs, mean, sd = int(r['observed']), float(r['expected']), float(r['null_sd'])
            lo, hi = (mean-1.96*sd)/mean, (mean+1.96*sd)/mean; kind = 'expected ± 1.96 SD'
        else:
            r = [r for r in mn if r['metric'] == 'WPP' and r['species'] == sp and r['contrast'] == 'P_S' and f2(r['threshold'])
                 and r['series_matching'] == 'True' and r['variant'] == var][0]
            obs, mean = int(r['observed']), float(r['null_mean'])
            lo, hi = float(r['null_p025'])/mean, float(r['null_p975'])/mean; kind = 'permutation 2.5–97.5%'
        ratio = obs/mean
        if not R.axis('F3B', f'{var} {sp} obs/null', [lo, hi, ratio], 0.75, 2.0): lo, hi, ratio = clipB(lo), clipB(hi), clipB(ratio)
        t.line([XB(lo), XB(hi)], [y, y], lw=1.7, color='#C9CFD3')
        t.marker(XB(ratio), y, 'o', 3.1, fc=INK if sp == 'human' else 'white', ec=INK, mew=0.7)
        Bd['controls'].append(dict(control=var, species=sp, observed=obs, null_mean=mean, ratio=ratio, lo=lo, hi=hi, interval=kind))
ya = 24.1+9.1*4+3.8+5.0
t.line([XB(1), XB(1)], [20.9, ya-1.5], lw=0.5, color='#989FA8', ls=(0, (1, 1.5)))
t.line([317.8, XB(2)], [ya, ya], lw=0.5, color='#989FA8')
for v, lab in [(1, '1'), (1.5, '1.5'), (2, '2')]:
    t.line([XB(v)]*2, [ya, ya+2], lw=0.5, color='#989FA8'); t.text_top(XB(v), ya+3.6, lab, 7.5, ha='center')
t.marker(255.2, ya+7.4, 'o', 3.1, mew=0.7); t.text_top(259.0, ya+3.6, 'Human', 7.5)
t.marker(287.5, ya+7.4, 'o', 3.1, fc='white', mew=0.7); t.text_top(291.3, ya+3.6, 'Mouse', 7.5)
y0 = ya+13.4
t.text_top(251.2, y0, 'Members retained', 8)
t.text_top(345.8, y0+0.3, 'Human', 7.5, ha='center'); t.text_top(412.6, y0+0.3, 'Mouse', 7.5, ha='center')
def bar(x, y, num_, den):
    t.rect(x, y-2.2, 32.3, 4.4, fc='#E6EAED'); t.rect(x, y-2.2, 32.3*num_/den, 4.4, fc='#263947')
    t.text_top(x+34.1, y-3.7, f'{num_}/{den}', 7.5)
rr = {(r['species'], r['contrast']): r for r in ret if r['variant'] == 'interval_supported' and r['metric'] == 'WPP' and f2(r['threshold'])}
yb = y0+13.0
assert 312.2-textwidth('Bootstrap interval', 7.5) > 238.5
t.text_top(312.2, yb-3.7, 'Bootstrap interval', 7.5, ha='right')   # 2026-10-08: same name as the upper rows and S10d
bh = (int(rr[('human', 'P_S')]['same_direction_intersection']), int(rr[('human', 'P_S')]['original']))
bm = (int(rr[('mouse', 'P_S')]['same_direction_intersection']), int(rr[('mouse', 'P_S')]['original']))
bar(317.8, yb, *bh); bar(384.5, yb, *bm)
# Well-level test (Schiebinger2019-2i, mouse P-M): its own gene background, shown apart from recurrent members.
# Denominators come from the statistics side (recurrent members with series_matching=True, the 3A/3B cohort).
nw, dw = int(well['all_tested']['joint_q05']), int(well['all_tested']['tested'])
nr, n_tested = int(well['recurrent_P_M_intersection']['joint_q05']), int(well['recurrent_P_M_intersection']['tested'])
yw = yb+9.1
t.text_top(251.2, yw-3.7, 'Wells, mouse P–M (DWUI q < .05)', 7.5, color=MUTED)
yw1 = yw+8.6; yw2 = yw1+9.1
for yy, lab, (a, b) in [(yw1, 'All tested genes', (nw, dw)), (yw2, 'Recurrent members', (nr, n_tested))]:
    assert 312.2-textwidth(lab, 7.5) > 238.5
    t.text_top(312.2, yy-3.7, lab, 7.5, ha='right')
    if b: bar(384.5, yy, a, b)
    else: R.unavailable('F3B', f'wells: {lab} has no tested genes')
SS = {}
for k, (kind, lab) in enumerate([('reference', 'S halves, same state'), ('target', 'P halves, same state')]):
    yy = yw2+9.1*(k+1)
    assert 312.2-textwidth(lab, 7.5) > 238.5, lab
    t.text_top(312.2, yy-3.7, lab, 7.5, ha='right')
    for sp, xc in [('human', 333.9), ('mouse', 400.6)]:
        r = [r for r in ss if r['pseudo_kind'] == kind and r['species'] == sp and r['contrast'] == 'P_S' and r['metric'] == 'WPP' and f2(r['threshold'])][0]
        a, b = float(r['pseudo_median']), float(r['observed_same_background_median'])
        fmt = lambda v: str(int(v)) if v == int(v) else f'{v:.1f}'
        t.text_top(xc, yy-3.7, f'{fmt(a)} / {fmt(b)}', 7.5, ha='center')
        SS[f'{kind}|{sp}'] = dict(pseudo_median=a, observed_median=b, eligible_median=float(r['eligible_median']))
assert yw2+9.1*2+4 < 145.0, yw2
Bd.update(bootstrap=dict(human=bh, mouse=bm), wells=(nw, dw),
          wells_recurrent=dict(tested=n_tested, dwui_q_lt_005=nr), same_state=SS)
ver['B'] = Bd

# ------------------------------------------------------------------ C
cnt = {r['species']: r for r in use('F3C', 'candidate_screen/dual_reference_ladder.tsv')}
cand = use('F3C', 'candidate_screen/dual_reference_candidates.tsv')
ev = use('F3C', 'recurrence/study_event_covariates.tsv.gz')
msv = use('F3C', 'matched_pms/matched_state_vectors.tsv.gz', lambda r: r['setting'] == 'supported50_p08')
memC = use('F3C', 'recurrence/recurrence_members.tsv.gz', lambda r: r['metric'] == 'WPP' and r['contrast'] == 'P_M' and f2(r['threshold'])
           and r['variant'] == 'primary' and r['series_matching'] == 'True')
valsC = {(r['route'], r['gene_id']): r for r in msv}
evk = {(r['unit'], r['gene_id']): int(r['event']) for r in ev if r['contrast'] == 'P_M' and r['metric'] == 'WPP' and f2(r['threshold'])}
LEAD = {}
def pick(sp, routes, slot, n):
    """Display rows: the example slot (pinned or rank rule) first, then the statistics-side dual-reference
    candidates by descending minimum absolute effect; values are the matched P/M/S state vectors."""
    cs = [dict(c) for c in cand if c['species'] == sp]
    ok = sorted([c for c in cs if all((x, c['gene_id']) in valsC for x in routes)], key=lambda c: -float(c['minimum_abs_effect']))
    if len(ok) != len(cs): R.flag('F3C', 'consistency', f'{sp}: {len(cs)-len(ok)} dual-reference candidates without matched vectors')
    for i_, c in enumerate(ok): c['rank'] = i_+1
    first = R.example(slot, 'F3C', ok)
    order = ([first] if first else [])+[c for c in ok if first is None or c['gene_id'] != first['gene_id']]
    if first: LEAD[sp] = first['gene_id']
    n_ = R.capacity('F3C', f'{sp} display rows', len(order), n, subset=True)
    return [(c['gene_id'], {x: valsC[(x, c['gene_id'])] for x in routes}) for c in order[:n_]], len(cs)
mrows, mdual = pick('mouse', ['CIP', 'WOT'], 'dual_reference_mouse', 10)
hrows, hdual = pick('human', ['GSE118258'], 'dual_reference_human', 8)
for sp_, d_ in [('mouse', mdual), ('human', hdual)]:
    if d_ != int(cnt[sp_]['dual_reference']): R.flag('F3C', 'consistency', f'{sp_}: candidate rows {d_} != ladder {cnt[sp_]["dual_reference"]}')
allv = [float(v[x][k]) for _, v in mrows for x in ('CIP', 'WOT') for k in ('DWUI_PS', 'DWUI_PM')] + \
       [float(v['GSE118258'][k]) for _, v in hrows for k in ('DWUI_PS', 'DWUI_PM')]
smin, smax = -0.75, 0.5
R.axis('F3C', 'strip ΔDWUI', allv, smin, smax)
STRIP = [('DWUI_PS', 'P–S ΔDWUI', 167.0), ('DWUI_PM', 'P–M ΔDWUI', 249.0)]; SW = 66.3
XS = lambda x0, v: x0+(v-smin)/(smax-smin)*SW
t.text_top(121.5, 157.1, 'Native P–M', 8, ha='center')
for k, lab, x0 in STRIP:
    t.text_top(x0+SW/2, 157.2, lab, 8, ha='center')
def axis_strip(ytop):
    for k, lab, x0 in STRIP:
        for v, s in [(-0.5, '−0.5'), (0, '0'), (0.5, '+0.5')]:
            t.text_top(XS(x0, v), ytop, s, 8, ha='center')
def sign(e):
    return {1: '+', -1: '−', 0: '0'}.get(e, '·')
Cm = []
def block(head, ytop, rows, units, routes):
    t.text_top(31.0, ytop, head, 9)
    for j, u in enumerate(units):
        c, m = STUDY[u if u in STUDY else u]; t.marker(104.0+j*17.5, ytop+4.6, m, 4.0, fc=c, ec=c, mew=0.5)
    axis_strip(ytop+5.2)
    y = ytop+18.5
    for g, v in rows:
        sym = v[routes[0]]['gene_symbol']; bold = g in LEAD.values()
        t.text_top(4.8, y-1.2, sym, 8.5, font='BolIta' if bold else 'RegIta')
        unit_ids = {'CIP': 'GSE103221', 'WOT': 'GSE115943'}
        for j, u in enumerate(units):
            e = evk.get((unit_ids.get(u, u), g))
            fc = {1: '#EDBB94', -1: '#9CC8E2', 0: '#F2F2F2'}.get(e, MISSING)
            t.rect(104.0+j*17.5-7.6, y-2.6, 15.2, 8.0, fc=fc)
            if e is not None:
                t.text_top(104.0+j*17.5, y-1.1, sign(e), 8, ha='center')
        for k, lab, x0 in STRIP:
            t.line([XS(x0, 0)]*2, [y-3.0, y+6.0], lw=0.5, color='#B8C2C8')
            for x in routes:
                c, m = STUDY[x]; val = min(max(float(v[x][k]), smin), smax)
                t.marker(XS(x0, val), y+1.6, m, 3.4, fc=c, ec=c, mew=0.4)
        Cm.append(dict(block=head, gene=sym, gene_id=g, signs={u: evk.get((unit_ids.get(u, u), g)) for u in units},
                       values={x: dict(PS=float(v[x]['DWUI_PS']), PM=float(v[x]['DWUI_PM'])) for x in routes}))
        y += 9.05
    return y
yend = block('Mouse', 168.9, mrows, ['CIP', 'WOT'], ['CIP', 'WOT'])
yend = block('Human', yend+6.0, hrows, ['GSE118258', 'GSE150311', 'GSE221739'], ['GSE118258'])
# scatters
dual_m = {c['gene_id'] for c in cand if c['species'] == 'mouse'}
_meas = lambda route: {r['gene_id'] for r in msv if r['route'] == route and r['DWUI_PS'] not in ('', 'nan') and r['DWUI_PM'] not in ('', 'nan')}
coh_m = {r['gene_id'] for r in memC if r['species'] == 'mouse'} & _meas('CIP') & _meas('WOT')   # matched-measurable in both routes
if len(coh_m) != int(cnt['mouse']['matched_measurable']): R.flag('F3C', 'consistency', f'P–M cohort {len(coh_m)} != ladder {cnt["mouse"]["matched_measurable"]}')
gm = {route: {g: g in dual_m for g in coh_m} for route in ('CIP', 'WOT')}
XSC = lambda v: 397.35+v*48.65
Cs = {}
for route, title, ytop in [('CIP', 'Guo2019–CIP', 165.9), ('WOT', 'Schiebinger2019–2i', 282.7)]:
    t.text_top(397.35, ytop, title, 8, ha='center')
    y0s = ytop+45.1                                   # value 0 on the y axis
    YSC = lambda v, y0s=y0s: y0s-v*58.35
    allp = [(float(r['DWUI_PS']), float(r['DWUI_PM']), r['gene_id']) for r in msv if r['route'] == route and r['setting'] == 'supported50_p08'
            and r['DWUI_PS'] not in ('', 'nan') and r['DWUI_PM'] not in ('', 'nan')]
    xs = np.array([p[0] for p in allp]); ys = np.array([p[1] for p in allp])
    keep = (np.abs(xs) <= 1) & (ys >= -0.5) & (ys <= 0.3)
    t.ax.scatter([XSC(v) for v in xs[keep]], [H-YSC(v) for v in ys[keep]], s=1.2, c='#C9CFD3', lw=0, rasterized=True, zorder=1)
    coh = [(a, b) for a, b, g in allp if g in gm[route]]; dual = [(a, b) for a, b, g in allp if gm[route].get(g)]
    t.ax.scatter([XSC(a) for a, b in coh], [H-YSC(b) for a, b in coh], s=4.0, facecolors='white', edgecolors='#5F6F78', lw=0.35, zorder=2)
    t.ax.scatter([XSC(a) for a, b in dual], [H-YSC(b) for a, b in dual], s=4.0, c=INK, lw=0, zorder=3)
    t.line([XSC(-1), XSC(1)], [YSC(0)]*2, lw=0.4, color='#B8C2C8', ls=(0, (1, 1.5)))
    t.line([XSC(0)]*2, [YSC(0.3), YSC(-0.5)], lw=0.4, color='#B8C2C8', ls=(0, (1, 1.5)))
    t.line([345.6, 345.6], [YSC(0.3), YSC(-0.5)], lw=0.5, color=AXC); t.line([XSC(-1), XSC(1)], [YSC(-0.5)+2.3]*2, lw=0.5, color=AXC)
    for v, s in [(0.3, '+0.3'), (0, '0'), (-0.5, '−0.5')]:
        t.line([343.6, 345.6], [YSC(v)]*2, lw=0.5, color=AXC); t.text_top(342.8, YSC(v)-4.0, s, 8, ha='right')
    for v, s in [(-1, '−1'), (0, '0'), (1, '+1')]:
        t.line([XSC(v)]*2, [YSC(-0.5)+2.3, YSC(-0.5)+4.3], lw=0.5, color=AXC); t.text_top(XSC(v), YSC(-0.5)+5.8, s, 8, ha='center')
    t.text_top(397.35, YSC(-0.5)+15.2, 'P–S ΔDWUI', 8, ha='center')
    t.vtext(322.0, y0s, 'P–M ΔDWUI', 8)
    n_out = int((~keep).sum())
    Cs[route] = dict(all_measurable=len(allp), outside_axes=n_out, cohort=len(coh), dual=len(dual))
    t.text_top(449.0, ytop+9.5, f'n = {len(allp)}', 7.5, ha='right', color=MUTED)
if not Cs['CIP']['dual'] == Cs['WOT']['dual'] == mdual: R.flag('F3C', 'consistency', f'dual genes in scatters {Cs["CIP"]["dual"]}/{Cs["WOT"]["dual"]} != {mdual}')
yl = 379.6
for x, u, lab in [(9.0, 'CIP', 'Guo2019–CIP'), (111.5, 'WOT', 'Schiebinger2019–2i')]:
    c, m = STUDY[u]; t.marker(x, yl+4.0, m, 4.0, fc=c, ec=c, mew=0.5); t.text_top(x+6.7, yl, lab, 8)
for x, u, lab in [(9.0, 'GSE118258', 'Xing2020'), (97.0, 'GSE150311', 'Liu2020'), (176.0, 'GSE221739', 'Panariello2023')]:
    c, m = STUDY[u]; t.marker(x, yl+17.1, m, 4.0, fc=c, ec=c, mew=0.5); t.text_top(x+6.7, yl+13.1, lab, 8)
for k, (sp, lab, tag) in enumerate([('mouse', 'Mouse', '(2 studies)'), ('human', 'Human', '(Xing2020 only)')]):
    c = cnt[sp]; yy = 405.5+12.3*k
    a, b, d = int(c['recurrent_candidates']), int(c['matched_measurable']), int(c['dual_reference'])
    t.text_top(4.8, yy, lab, 8); t.text_top(51.1, yy, f'{b}/{a} measurable', 8)
    t.rect(128.0, yy+1.6, 78.0, 4.6, fc='#E6EAED'); t.rect(128.0, yy+1.6, 78.0*d/max(b, 1), 4.6, fc='#263947')
    t.text_top(209.4, yy, f'{d}/{b} aligned', 8)
for k, (lab, kind) in enumerate([('All measurable', 'all'), (f'P–M cohort ({int(cnt["mouse"]["matched_measurable"])})', 'coh'),
                                 (f'P–S/P–M aligned ({mdual})', 'dual')]):
    yy = 391.2+12.55*k; x = 343.0
    if kind == 'all': t.marker(x, yy+4.0, 'o', 2.6, fc='#C9CFD3', ec='#C9CFD3', mew=0)
    elif kind == 'coh': t.marker(x, yy+4.0, 'o', 3.0, fc='white', ec='#5F6F78', mew=0.5)
    else: t.marker(x, yy+4.0, 'o', 3.0, fc=INK, ec=INK, mew=0)
    t.text_top(349.2, yy, lab, 8)
ver['C'] = dict(display_rule='example slot first (display-selection.json), then dual-reference candidates by descending minimum_abs_effect',
                rows=Cm, scatter=Cs, counts={k: dict(v) for k, v in cnt.items()})

# ------------------------------------------------------------------ D
exD = R.example('dual_reference_mouse', 'F3D', use('F3D', 'candidate_screen/example_selection.tsv', lambda r: r['role'] == 'dual_reference' and r['species'] == 'mouse'))
GID = exD['gene_id'] if exD else None; GSYM = exD['gene_symbol'] if exD else ''
if exD and LEAD.get('mouse') != GID: R.flag('F3D', 'consistency', '3D example differs from the first 3C mouse row')
pd7 = {r['route']: r for r in use('F3D', 'matched_pms/matched_state_vectors.tsv.gz', lambda r: r['gene_id'] == GID and r['setting'] == 'supported50_p08')}
R.path('F3D', 'panel_statistics/example-cumulative-PAS-curves.tsv.gz')   # registered; does not hold every selected gene (GAP-D03)
R.derive('F3D', 'state and per-triplet cumulative PAS-usage steps and per-triplet state DWUI drawn from matched_state_vectors / weighted_triplets usage (GAP-D03)')
if set(pd7) != {'CIP', 'WOT'}:
    R.unavailable('F3D', f'{GSYM or "no example"}: matched P/M/S vectors for CIP and WOT not both available; panel left empty')
t.text_top(31.6, 439.3, GSYM, 9, font='BolIta')
YD = lambda v: 470.6+(0.95-v)/(0.95-0.15)*35.0
Dd = {}
import pandas as _pd
trip = R.df('F3D', 'matched_pms/weighted_triplets.tsv.gz', usecols=['gene_id', 'route', 'setting', 'contrast', 'design', 'library_key_a', 'library_key_b', 'delta_DWUI', 'triplet_weight',
                                          'PAS_ids', 'M_PAS_ids', 'target_usage', 'reference_usage', 'M_usage'])
trip = trip[(trip.gene_id == GID) & (trip.setting == 'supported50_p08') & (trip.contrast == 'P_S') & (trip.design == 'true_source')]
SXs = lambda v, dx: 46.0+dx+(v+0.7)/0.5*58.0
for k, (route, title) in enumerate([('CIP', 'Guo2019–CIP'), ('WOT', 'Schiebinger2019–2i')] if set(pd7) == {'CIP', 'WOT'} else []):
    dx = 105.6*k; r = pd7[route]
    _nt = int((trip.route == route).sum()); _tt = f'{title}, n = {_nt}'
    t.text_top(74.3+dx, 455.9, _tt, 7.5, ha='center'); assert 74.3+dx+textwidth(_tt, 7.5)/2 < 231.0
    t.line([34.3+dx, 34.3+dx], [YD(0.95), YD(0.15)], lw=0.5, color=AXC)
    for v, s in [(0.9, '0.9'), (0.55, '0.55'), (0.2, '0.2')]:
        t.line([32.3+dx, 34.3+dx], [YD(v)]*2, lw=0.8, color=INK); t.text_top(30.7+dx, YD(v)-3.8, s, 7.5, ha='right')
    xs_ = {'S': 45.8+dx, 'P': 74.3+dx, 'M': 102.9+dx}
    t.line([88.6+dx]*2, [YD(0.95), YD(0.15)], lw=0.5, color='#AAB4BB', ls=(0, (1, 1.5)))
    vals = {s: float(r[f'{s}_DWUI']) for s in 'SPM'}
    R.axis('F3D', f'{route} state DWUI', list(vals.values()), 0.15, 0.95)
    vals_raw = dict(vals); vals = {s: min(max(v_, 0.15), 0.95) for s, v_ in vals.items()}
    t.line([xs_['S'], xs_['P']], [YD(vals['S']), YD(vals['P'])], lw=0.8, color='#ADB9C2')
    for s, mk in [('S', '^'), ('P', 'o'), ('M', 'D')]:
        t.marker(xs_[s], YD(vals[s]), mk, 4.2 if s != 'M' else 3.6, fc=SC[s], ec=SC[s], mew=0.5)
        t.text_top(xs_[s], 509.2, s, 7.5, ha='center', color=SC[s])
    tq = trip[trip.route == route]
    # per-triplet state DWUI (usage . normalised distal distance; reproduces delta_DWUI and the state means exactly)
    _ids = json.loads(r['PAS_ids']); _pos = np.array([pas_pos(i) for i in _ids]); _minus = pas_strand(_ids[0]) == '-'
    _dn = ((_pos.max()-_pos) if _minus else (_pos-_pos.min()))/(_pos.max()-_pos.min())
    for s_ in 'SPM':
        if abs(float(np.dot(json.loads(r[f'{s_}_usage']), _dn))-vals_raw[s_]) > 1e-9:
            R.flag('F3D', 'consistency', f'{route} {s_}: usage does not reproduce the state DWUI')
    tdw = []
    for _, tr in tq.iterrows():
        dS, dP, dM = (float(np.dot(json.loads(tr[c_]), _dn)) for c_ in ('reference_usage', 'target_usage', 'M_usage'))
        if abs((dP-dS)-float(tr.delta_DWUI)) > 1e-9: R.flag('F3D', 'consistency', f'{route}: triplet usage does not reproduce delta_DWUI')
        R.axis('F3D', f'{route} triplet state DWUI', [dS, dP, dM], 0.15, 0.95)
        tdw.append(tuple(min(max(v_, 0.15), 0.95) for v_ in (dS, dP, dM)))
    for dS, dP, dM in sorted(set(tdw)):
        t.line([xs_['S'], xs_['P']], [YD(dS), YD(dP)], lw=0.4, color='#C9D1D6', zorder=1)
    for x_, k_ in [(xs_['S'], 0), (xs_['P'], 1), (xs_['M'], 2)]:
        vv = sorted({v_[k_] for v_ in tdw})
        t.ax.plot([x_]*len(vv), [t.H-YD(v_) for v_ in vv], ls='none', marker='o', ms=2.0, mfc='#AAB5BC', mec='none', zorder=2)
    wm = float(np.average(tq.delta_DWUI, weights=tq.triplet_weight))
    if abs(wm-float(r['DWUI_PS'])) > 1e-6: R.flag('F3D', 'consistency', f'{route}: weighted triplet mean != DWUI_PS')
    if r.get('n_triplets', '') not in ('', 'nan') and int(float(r['n_triplets'])) != len(tq): R.flag('F3D', 'consistency', f'{route}: n_triplets != weighted triplets')
    vals_t = sorted(tq.delta_DWUI.astype(float))
    clipD = lambda v: min(max(v, -0.7), -0.2)
    R.axis('F3D', f'{route} ΔDWUI strip', vals_t+[float(r['DWUI_PS']), float(r['DWUI_PM'])], -0.7, -0.2)
    vals_t = [clipD(v_) for v_ in vals_t]
    t.text_top(43.0+dx, 521.5-3.8, 'P–S', 7.5, ha='right'); t.text_top(43.0+dx, 528.8-3.8, 'P–M', 7.5, ha='right')
    for k_, v_ in enumerate(vals_t):
        t.marker(SXs(v_, dx), 521.5+(-1.4+2.8*(k_ % 3)/2), 'o', 2.0, fc='#AAB5BC', ec='none', mew=0, zorder=3)
    t.marker(SXs(clipD(float(r['DWUI_PS'])), dx), 521.5, 'o', 3.6, mew=0.3, zorder=5)
    t.marker(SXs(clipD(float(r['DWUI_PM'])), dx), 528.8, 'D', 3.4, mew=0.3, zorder=5)
    t.line([SXs(-0.7, dx), SXs(-0.2, dx)], [532.6]*2, lw=0.5, color=AXC)
    for v_, s_ in [(-0.6, '−0.6'), (-0.4, '−0.4'), (-0.2, '−0.2')]:
        t.line([SXs(v_, dx)]*2, [532.6, 534.4], lw=0.5, color=AXC); t.text_top(SXs(v_, dx), 535.4, s_, 7.5, ha='center')
    t.text_top(43.0+dx, 535.4, 'ΔDWUI', 7.5, ha='right')
    ids = json.loads(r['PAS_ids']); pos = [pas_pos(i) for i in ids]; strand = pas_strand(ids[0])
    prox = max(pos) if strand == '-' else min(pos); dist = [abs(prox-p)/1000 for p in pos]
    if max(dist) > 9.6: R.flag('F3D', 'axis_range', f'{route}: PAS span {max(dist):.1f} kb beyond the 9.6 kb axis', data_range=[0, max(dist)], axis=[0, 9.6])
    XK = lambda kb, dx=dx: 38.7+dx+kb*8.6
    for kb_ in sorted(set(round(v_, 3) for v_ in dist)):
        t.marker(XK(kb_), 560.8, 'v', 3.2, fc=INK, ec=INK, mew=0.3)
    if dx == 0:
        t.text_top(30.7, 557.0, 'PAS', 7.5, ha='right')
    t.line([34.3+dx, 34.3+dx], [567.1, 611.2], lw=0.5, color=AXC); t.line([34.3+dx, XK(9.6)], [613.9, 613.9], lw=0.5, color=AXC)
    for v, s in [(100, '100'), (50, '50'), (0, '0')]:
        yv = 611.2-v/100*44.1; t.line([32.3+dx, 34.3+dx], [yv]*2, lw=0.8, color=INK); t.text_top(30.7+dx, yv-3.8, s, 7.5, ha='right')
    for v in (0, 4, 8):
        t.line([XK(v)]*2, [613.9, 615.9], lw=0.8, color=INK); t.text_top(XK(v), 617.7, str(v), 7.5, ha='center')
    curves = {}
    order_ = sorted(range(len(pos)), key=lambda i: dist[i]); thin = {}
    for s, col_ in [('S', 'reference_usage'), ('P', 'target_usage'), ('M', 'M_usage')]:
        uniq = []
        for _, tr in tq.iterrows():
            assert json.loads(tr.PAS_ids) == ids and json.loads(tr.M_PAS_ids) == ids
            u_ = tuple(round(float(z), 12) for z in json.loads(tr[col_]))
            if u_ not in uniq:
                uniq.append(u_)
        for u_ in uniq:
            cur = 1.0; pts = [(XK(0), 611.2-cur*44.1)]
            for i in order_:
                pts.append((XK(dist[i]), 611.2-cur*44.1)); cur -= u_[i]; pts.append((XK(dist[i]), 611.2-max(cur, 0)*44.1))
            pts.append((XK(9.6), 611.2-max(cur, 0)*44.1))
            t.line([p[0] for p in pts], [p[1] for p in pts], lw=0.45, color=SC[s], alpha=0.4, zorder=1)
        thin[s] = [list(u_) for u_ in uniq]
    pu = np.average(np.array([json.loads(z) for z in tq.target_usage]), axis=0, weights=tq.triplet_weight)
    assert np.allclose(pu, json.loads(r['P_usage']), atol=1e-6), (route, pu)
    for s, ls in [('S', (0, (3, 1.6))), ('P', '-'), ('M', (0, (4, 1.5, 1, 1.5)))]:
        u = json.loads(r[f'{s}_usage']); order = sorted(range(len(pos)), key=lambda i: dist[i])
        xs_c, ys_c = [XK(0)], []; cur = 1.0
        pts = [(XK(0), 611.2-cur*44.1)]
        for i in order:
            pts.append((XK(dist[i]), 611.2-cur*44.1)); cur -= u[i]; pts.append((XK(dist[i]), 611.2-max(cur, 0)*44.1))
        pts.append((XK(9.6), 611.2-max(cur, 0)*44.1))
        t.line([p[0] for p in pts], [p[1] for p in pts], lw=1.05, color=SC[s], ls=ls)
        curves[s] = dict(usage=u, distances_kb=dist)
    Dd[route] = dict(S_DWUI=vals_raw['S'], P_DWUI=vals_raw['P'], M_DWUI=vals_raw['M'], DWUI_PS=float(r['DWUI_PS']), DWUI_PM=float(r['DWUI_PM']),
                     n_target_libraries=int(r['n_target_libraries']), n_source_libraries=int(r['n_source_libraries']), pas=ids, curves=curves,
                     triplet_usage_unique=thin, triplet_state_DWUI=[dict(S=a_, P=b_, M=c_) for a_, b_, c_ in tdw],
                     PS_triplets=[dict(target=a_, source=b_, delta_DWUI=float(c_), weight=float(w_)) for a_, b_, c_, w_ in
                                  zip(tq.library_key_a, tq.library_key_b, tq.delta_DWUI, tq.triplet_weight)])
t.vtext(9.0, 487.8, 'DWUI', 7.5); t.vtext(11.1, 589.1, 'Distal (%)', 7.5)   # 2026-10-08: as Figures 4e, 5e and S12e,h
for x0, s, ls in [(45.8, 'S', (0, (3, 1.6))), (97.2, 'P', '-'), (148.5, 'M', (0, (4, 1.5, 1, 1.5)))]:
    t.line([x0, x0+14.2], [552.8]*2, lw=1.05, color=SC[s], ls=ls); t.text_top(x0+17.1, 549.0, s, 7.5)
t.text_top(116.2, 629.4, 'Distance from proximal PAS (kb)', 7.5, ha='center')
t.marker(186.0, 552.8, 'o', 2.0, fc='#AAB5BC', ec='none', mew=0); t.text_top(190.0, 549.0, 'Triplet', 7.5)
assert 190.0+textwidth('Triplet', 7.5) < 230.0
ver['D'] = dict(values=Dd, note='State values are matched-contributor state means (matched_state_vectors); metacell quantile ranges need a server export and are not drawn.')

# ------------------------------------------------------------------ E
fam = use('F3E', 'apa_gex/prespecified_eight_test_family.tsv')
nul = use('F3E', 'apa_gex/four_class_nulls.tsv.gz')
# 2 x 2 class tiles per species: rows APA same/opposite, columns GEX same/opposite (cells A, B / C, D)
from matplotlib.colors import LinearSegmentedColormap
cmr = LinearSegmentedColormap.from_list('r', DIV)
lr = lambda e: cmr((np.clip(np.log2(e), -0.8, 0.8)+0.8)/1.6)
Tiles = {}
comp_ = {sp: int([r for r in fam if r['species'] == sp and f2(r['threshold']) and r['baseline_bins'] == '5'][0]['comparable_genes'])
         for sp in ('human', 'mouse')}
cw_, ch_ = 38.0, 22.0
for k, (sp, lab) in enumerate([('human', 'Human'), ('mouse', 'Mouse')]):
    gx = 286.0+k*84.0
    nn = [r for r in nul if r['species'] == sp and r['null_model'] == 'response_identity' and f2(r['threshold']) and r['baseline_bins'] == '5']
    shuf = min({r['shuffled_study'] for r in nn}, key=lambda s: float([r for r in nn if r['shuffled_study'] == s and r['cell'] == 'B'][0]['enrichment']))
    cells = {r['cell']: r for r in nn if r['shuffled_study'] == shuf}
    hl, nl = lab, f' {comp_[sp]:,} genes'
    wtot = textwidth(hl, 8, 'Bol')+textwidth(nl, 7.5); xs0 = gx+cw_-wtot/2
    t.text_top(xs0, 452.0, hl, 8, font='Bol'); t.text_top(xs0+textwidth(hl, 8, 'Bol'), 452.5, nl, 7.5, color=MUTED)
    t.text_top(gx+cw_/2, 462.5, 'Same', 7.5, ha='center', color=MUTED); t.text_top(gx+1.5*cw_, 462.5, 'Opp.', 7.5, ha='center', color=MUTED)
    for ci_, (cell, rr, cc) in enumerate([('A', 0, 0), ('B', 0, 1), ('C', 1, 0), ('D', 1, 1)]):
        r = cells[cell]; e = float(r['enrichment']); o = int(r['observed'])
        x = gx+cc*cw_; y = 472.0+rr*ch_
        t.rect(x, y, cw_-1.0, ch_-1.0, fc=lr(e), ec=INK if cell == 'B' else 'none', lw=1.0 if cell == 'B' else 0)
        t.text_top(x+(cw_-1)/2, y+2.6, str(o), 8, ha='center', font='Bol' if cell == 'B' else 'Reg')
        t.text_top(x+(cw_-1)/2, y+11.4, f'×{e:.2f}', 7.0, ha='center', color=MUTED)
        Tiles[f'{sp}|{cell}'] = dict(observed=o, null_mean=float(r['null_mean']), enrichment=e, shuffled_study=shuf, comparable_genes=comp_[sp])
    if k == 0:
        # Factor names carry a delta: the tiles classify direction of change, not levels.
        t.text_top(gx-3.0, 472.0+6.6, 'Same', 7.5, ha='right', color=MUTED); t.text_top(gx-3.0, 472.0+ch_+6.6, 'Opp.', 7.5, ha='right', color=MUTED)
        assert 472.0+6.6+7.5 < 472.0+ch_/2+6.6 and 472.0+ch_/2+6.6+7.5 < 472.0+ch_+6.6  # sits between the two row labels
        t.text_top(248.0, 472.0+ch_/2+6.6, 'ΔAPA', 7.5, color=MUTED)
t.text_top(248.0, 464.0, 'ΔGEX', 7.5, color=MUTED)
for i_ in range(40):
    v = -0.8+1.6*(i_+0.5)/40
    t.rect(286.0+i_*1.6, 521.0, 1.7, 3.6, fc=cmr((v+0.8)/1.6))
for v, s_ in [(-0.8, '×0.57'), (0, '×1'), (0.8, '×1.74')]:
    t.text_top(286.0+(v+0.8)/1.6*64.0, 526.6, s_, 7.0, ha='center', color=MUTED)
t.rect(372.0, 520.6, 6.0, 4.4, fc='white', ec=INK, lw=1.0); t.text_top(381.0, 519.2, 'Tested class', 7.5, color=MUTED)
XE = lambda v: 318.0+(v-0.7)/1.0*96.0
t.text_top(366.0, 542.0, 'Obs. / null mean', 8, ha='center'); t.text_top(438.0, 542.0, 'q', 8, ha='center')
Ed = []; y = 556.0
for sp, lab in [('human', 'Human'), ('mouse', 'Mouse')]:
    fr = [r for r in fam if r['species'] == sp and f2(r['threshold']) and r['baseline_bins'] == '5']
    assert len(fr) == 2
    obs = int(fr[0]['observed']); comp = int(fr[0]['comparable_genes'])
    for model, mlab in [('response_identity', 'identity'), ('direction_specific', 'direction')]:
        f = [r for r in fr if r['null_model'] == model][0]
        nn = [r for r in nul if r['species'] == sp and r['null_model'] == model and f2(r['threshold']) and r['baseline_bins'] == '5' and r['cell'] == 'B']
        assert len(nn) == 2
        n = min(nn, key=lambda r: float(r['enrichment']))
        mean = float(n['null_mean']); lo, hi = float(n['null_p025'])/mean, float(n['null_p975'])/mean; e = float(n['enrichment'])
        assert abs(e-float(f['enrichment_min'])) < 1e-9 and int(n['observed']) == obs
        t.rect(246.0, y-7.0, 203.0, 14.0, fc=BAND) if mlab == 'identity' else None
        t.text_top(248.0, y-4.0, f'{lab} {mlab}', 8)
        t.line([XE(lo), XE(hi)], [y-2.7, y-2.7], lw=1.7, color='#C9CFD3')
        q = float(f['q_family'])
        t.marker(XE(e), y-2.7, 'o', 4.0, fc=INK if q < 0.05 else 'white', ec=INK, mew=0.7)
        t.text_top(438.0, y-7.4, qfmt(q), 8, ha='center')
        f5_ = [r_ for r_ in fam if r_['species'] == sp and f5(r_['threshold']) and r_['baseline_bins'] == '5' and r_['null_model'] == model]
        assert len(f5_) == 1; f5_ = f5_[0]
        nn5 = [r_ for r_ in nul if r_['species'] == sp and r_['null_model'] == model and f5(r_['threshold']) and r_['baseline_bins'] == '5' and r_['cell'] == 'B']
        assert len(nn5) == 2; n5 = min(nn5, key=lambda r_: float(r_['enrichment']))
        m5 = float(n5['null_mean']); lo5, hi5, e5 = float(n5['null_p025'])/m5, float(n5['null_p975'])/m5, float(n5['enrichment'])
        assert abs(e5-float(f5_['enrichment_min'])) < 1e-9 and int(n5['observed']) == int(f5_['observed'])
        q5 = float(f5_['q_family'])
        t.line([XE(lo5), XE(hi5)], [y+2.9, y+2.9], lw=1.7, color='#C9CFD3')
        t.marker(XE(e5), y+2.9, 'D', 3.6, fc=INK if q5 < 0.05 else 'white', ec=INK, mew=0.7)
        t.text_top(438.0, y+0.4, qfmt(q5), 8, ha='center')
        Ed.append(dict(species=sp, model=model, observed=obs, comparable=comp, null_mean=mean, enrichment=e, lo=lo, hi=hi, q=q,
                       shuffled_study=n['shuffled_study'],
                       t005=dict(observed=int(f5_['observed']), null_mean=m5, enrichment=e5, lo=lo5, hi=hi5, q=q5, shuffled_study=n5['shuffled_study'])))
        y += 14.5
    y += 4.0
R.axis('F3E', 'obs/null mean', [v_ for d in Ed for v_ in (d['lo'], d['hi'], d['enrichment'], d['t005']['lo'], d['t005']['enrichment'])],
       0.7+(312-318.0)/96.0, 0.7+(424-318.0)/96.0)
ya = y-9.0
t.line([XE(1)]*2, [549.0, ya-3], lw=0.5, color='#98A1A8', ls=(0, (1, 1.5)))
t.line([XE(0.8), XE(1.8)], [ya, ya], lw=0.5, color=INK)
for v, s_ in [(1, '1'), (1.5, '1.5')]:
    t.line([XE(v)]*2, [ya, ya+2], lw=0.5, color=INK); t.text_top(XE(v), ya+3.6, s_, 7.5, ha='center')
yl_ = ya+14.0
t.marker(251.0, yl_+3.8, 'o', 4.0, fc=INK, ec=INK, mew=0.7); t.text_top(256.0, yl_, 'q < .05', 7.5, color=MUTED)
t.marker(292.0, yl_+3.8, 'o', 4.0, fc='white', ec=INK, mew=0.7); t.text_top(297.0, yl_, 'q ≥ .05', 7.5, color=MUTED)
_x = 337.0; t.text_top(_x, yl_, '|ΔDWUI| ≥', 7.5, color=MUTED); _x += textwidth('|ΔDWUI| ≥', 7.5)+5.0
t.marker(_x, yl_+3.8, 'o', 4.0, fc=INK, ec=INK, mew=0.7); t.text_top(_x+4.5, yl_, '0.02', 7.5, color=MUTED); _x += 4.5+textwidth('0.02', 7.5)+7.0
t.marker(_x, yl_+3.8, 'D', 3.6, fc=INK, ec=INK, mew=0.7); t.text_top(_x+4.5, yl_, '0.05', 7.5, color=MUTED)
assert _x+4.5+textwidth('0.05', 7.5) < 450.0, _x
assert yl_+8.0 < 639.0, yl_
ver['E'] = dict(tests=Ed, tiles=Tiles)

t.save(R.out/'overlay.pdf'); overlay(pdf, R.out/'overlay.pdf', 0, 0); pdf.save(OUT)
render_png(OUT, R.out/'Figure.3.draft.png', 150)
ver['output_sha256'] = sha(OUT)
R.finish(3, OUT, json.loads(json.dumps(ver, default=str)))
