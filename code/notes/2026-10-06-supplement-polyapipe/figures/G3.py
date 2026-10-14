"""Supplementary Figure G3 — stability of the frozen expression partitions and within-library programme contrasts of
the three boundary states (printed number from supdraw.SNUM).

Re-laid-out from the 2026-09-14 state-grouping supplement (SMC4 A–B and SMC3, submission-v3/render.py), same inputs and
values. a: 11 original partitions, nine neighbouring Leiden runs each (adjacent resolution x three seeds): original-family
purity (B-cubed precision) and ARI, median and range, with original / alternative cluster counts. b–d: target minus
reference programme score in each library where both groups have >= 20 cells. Not shown: SMC4 C–F (APA membership after
excluding the boundary states), which was computed on the earlier C1–C4 profile classes and does not apply to the current
polyApipe k = 3 classes.
"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'tools'))
from supdraw import *   # noqa

FID = 'G3'
D = Data(); X = ExtData(panel='G3A'); OUT = outdir(FID)
RUN = X.df('partition_runs'); DIFF = X.df('within_library_program_differences')
DIFF = DIFF[DIFF.both_groups_at_least20]
X0, RR = MARGIN, PAGE_W-MARGIN
H1, H2 = 136.0, 111.5
Y1 = MARGIN; Y2 = Y1+H1+GAP; H = Y2+H2+MARGIN
F = Fig(H, figname(FID))
LIBC = AXC

# ================================================================ a  partition stability
UNITS = [('gse118258', 'Xing2020', '#753399', ''), ('gse150311', 'Liu2020', '#648290', ''),
         ('gse178325', 'Guan2022', '#A82C60', ''), ('gse221739', 'Panariello2023', '#30363D', ''),
         ('wot_common_to_2i', 'Schiebinger2019–2i', '#207450', ''), ('gse103221_cip', 'Guo2019–CIP', '#CF4C00', ''),
         ('gse103221_osk', 'Guo2019–OSK', '#255FBC', ''), ('gse171127_d3_getm', 'Jaber2022', '#766A8A', 'D3 GETM'),
         ('gse171127_d3_oskm', 'Jaber2022', '#B3A4C4', 'D3 OSKM'), ('gse171127_d6_pooled', 'Jaber2022', '#766A8A', 'D6'),
         ('gse171127_d12_pooled', 'Jaber2022', '#766A8A', 'D12')]
x, y, w, h = F.panel('a', X0, Y1, RR-X0, H1)
labw = max(textwidth(n, 7.0)+(textwidth(' '+s, 7.0, 'Reg') if s else 0) for _, n, _, s in UNITS)
xl = x+0.4                                    # rows start below the letter clear zone
cnt_w = 52.0
pa0 = xl+7.6+labw+10.0; vw = 16.0
pw_ = (x+w-cnt_w-pa0-2*vw-12.0)/2
pa = Lin(0, 1, pa0, pa0+pw_)
pb = Lin(0, 1, pa.p1+vw+12.0, pa.p1+vw+12.0+pw_)
F.text((pa.p0+pa.p1)/2, y+3.0, 'Original-family purity (B-cubed precision)', 7.5, ha='center')
F.text((pb.p0+pb.p1)/2, y+3.0, 'ARI with the original partition', 7.5, ha='center')
xc = x+w-cnt_w/2
F.text(xc, y+3.0, 'Clusters', 7.5, ha='center'); F.text(xc, y+11.0, 'original / runs', 7.0, ha='center', color=MUTED)
pitch = 8.6; r0 = y+24.0
for s_ in (pa, pb):
    F.line([s_(1.0)]*2, [r0-5.0, r0+pitch*10+5.0], lw='ref', color=AXC, ls='ref')
Av = {}
for j, (u, nm, col, sub) in enumerate(UNITS):
    yy = r0+pitch*j
    if j % 2 == 0: F.rect(x, yy-pitch/2, w, pitch, fc=BAND)
    F.rect(xl, yy-2.8, 5.6, 5.6, fc=col)
    tw = F.text(xl+7.6, yy-3.4, nm, 7.0)
    if sub: F.text(xl+7.6+tw+2.0, yy-3.4, sub, 7.0, color=MUTED)
    q = RUN[(RUN.unit == u) & (RUN.role == 'sensitivity')]; assert len(q) == 9
    base = RUN[(RUN.unit == u) & (RUN.role == 'baseline_replay')].iloc[0]
    assert float(base.ARI) == 1.0
    Av[u] = {}
    for s_, met in ((pa, 'family_Bcubed_precision'), (pb, 'ARI')):
        v = q[met]
        F.line([s_(v.min()), s_(v.max())], [yy, yy], lw='bold', color=NULLBAR)
        F.marker(s_(v.median()), yy, 'o', 'small', fc=INK)
        F.text(s_.p1+3.0, yy-3.4, f'{v.median():.2f}'.replace('0.', '.').replace('1.00', '1'), 6.8, color=MUTED)
        Av[u][met] = dict(min=float(v.min()), median=float(v.median()), max=float(v.max()), runs=9)
    F.text(xc, yy-3.4, f'{int(base.original_states)} / {int(q.alternative_clusters.min())}–{int(q.alternative_clusters.max())}', 7.0, ha='center')
    Av[u]['clusters'] = dict(original=int(base.original_states), runs_min=int(q.alternative_clusters.min()),
                             runs_max=int(q.alternative_clusters.max()))
ay = r0+pitch*10+pitch/2+1.0
for s_ in (pa, pb):
    F.xaxis(s_, ay, [0, 0.25, 0.5, 0.75, 1.0], ['0', '.25', '.5', '.75', '1'])
F.key_row(None, y+h-7.4, [('marker', ('o', INK), 'Median of 9 runs'), ('line', (NULLBAR, None, 'bold'), 'Range')],
          gap=8.0, center=(pa.p0+pb.p1)/2)
F.values['a'] = Av

# ================================================================ b–d  boundary programme contrasts
CASES = [('b', 'WOT L00 vs other R', 'Schiebinger2019–2i', '#207450', 'L00 vs L03, L06',
          [('fibroblast_ecm_retention', 'Source / ECM'), ('cell_cycle_g2m', 'G2/M')], 'Δ AUCell', (-0.2, 0.6), [-0.2, 0, 0.2, 0.4, 0.6]),
         ('c', 'CIP C012 vs other E', 'Guo2019–CIP', '#CF4C00', 'C012 vs C008, C009, C011',
          [('fibroblast_ecm_retention', 'Source / ECM'), ('met_epithelialization', 'Epithelial'),
           ('pluripotency_acquisition', 'Pluripotency'), ('stress_ap1', 'AP1 / stress')], 'Δ AUCell', (-0.05, 0.15), [-0.05, 0, 0.05, 0.1, 0.15]),
         ('d', 'Human C005 vs C004/C006', 'Xing2020', '#753399', 'C005 vs C004, C006',
          [('Source / ECM', 'ECM'), ('Epithelial', 'Epithelial'), ('Pluripotency', 'Pluripotency')],
          'Δ mean log1p(CP10k)', (-0.4, 0.4), [-0.4, -0.2, 0, 0.2, 0.4])]
WC = (RR-X0-2*GAP)/3
PROGW = max(textwidth(lab, 7.0) for *_, progs, _, _, _ in [(c[0], c[5], 0, 0, 0) for c in CASES] for _, lab in progs)
Vb = {}
for i, (L, comp, study, col, title, progs, metric, rng, ticks) in enumerate(CASES):
    x, y, w, h = F.panel(L, X0+(WC+GAP)*i, Y2, WC, H2)
    tx = x+0.2+textwidth(L, 11, 'Bol')+F.CLEAR_R+0.6
    F.rect(tx, y+3.2, 5.6, 5.6, fc=col); F.text(tx+7.6, y+2.4, study, 8.0)
    # 2026-10-12: comparison subtitle (title) moved to the caption
    cntw = textwidth('29 / 29', 7.0)
    sx = Lin(rng[0], rng[1], x+PROGW+5.0, x+w-cntw-6.0)
    top, bot = y+18.0, y+71.0
    pitch_ = (bot-top)/len(progs)
    F.text(x+w, top-6.6, 'Δ > 0', 7.0, ha='right', color=MUTED)
    F.line([sx(0)]*2, [top, bot], lw='ref', color=AXC, ls='ref')
    Vb[L] = {}
    for j, (pk, lab) in enumerate(progs):
        yy = top+pitch_*(j+0.5)
        q = DIFF[(DIFF.comparison == comp) & (DIFF.program == pk)].sort_values('library_key')
        assert len(q)
        F.text(x+0.4, yy-3.4, lab, 7.0)
        jit = np.linspace(-1, 1, len(q))*min(0.3*pitch_, 6.0) if len(q) > 1 else np.zeros(1)
        F.points([sx.clip(v) for v in q.delta], yy+jit, ms='small', color=LIBC, alpha='full' if len(q) < 5 else 'secondary')
        med = float(q.delta.median())
        if len(q) > 1:
            F.line([sx(med)]*2, [yy-min(0.36*pitch_, 7.0), yy+min(0.36*pitch_, 7.0)], lw='bold', color=INK)
        F.text(x+w, yy-3.4, f'{int((q.delta > 0).sum())} / {len(q)}', 7.0, ha='right')
        Vb[L][pk] = dict(n_libraries=int(len(q)), median=med, positive=int((q.delta > 0).sum()),
                         values=[float(v) for v in q.delta], metric=str(q.metric.iloc[0]))
        assert q.delta.between(rng[0], rng[1]).all(), (L, pk)
    F.xaxis(sx, bot+1.0, ticks, [(f'{t:g}'.replace('-0.', '−.').replace('0.', '.') if t else '0') for t in ticks], metric)
    items = [('marker', ('o', LIBC), 'One library')] + ([('line', (INK, None, 'bold'), 'Median')] if L != 'd' else [])
    F.key_row(None, y+h-7.4, items, gap=8.0, center=(sx.p0+sx.p1)/2)
F.values['bcd'] = Vb
F.notes.append('b–d: libraries in which both target and reference groups have >= 20 cells; programme scores as in the '
               '2026-09-13 marker analysis (existing AUCell for mouse, fixed-panel mean log1p(CP10k) for human)')
F.notes.append('SMC4 C–F (APA membership after excluding L00/C012/C005) not drawn: computed on the earlier C1–C4 classes')
F.save(OUT, [X, D])
