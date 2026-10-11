"""Supplementary Figure S10 (v11): event classes of the shared P–M proximal response and corresponding PAS switches.
a  event-class matrix of the 751-pair P–M background (human x mouse class), shared proximal pairs at |dDWUI| >= 0.02
b  composition of the shared proximal set at nine thresholds (both tandem / both ALE-compatible / other or mixed)
c  direction pairing within the tandem and ALE classes under both null models at 0.02 (observed count, 95% null range)
d  DYNLL2/Dynll2: corresponding PAS switch from source to P (segment, homologous seed matches, inclusion per study)
e  attrition of the 48 shared orthologs to corresponding PAS switches at 25 and 50 nt
Data: notes/2026-10-10-utr-event-stratification, notes/2026-10-10-conserved-response/sites (read-only)."""
import math
import numpy as np
from common import *                                     # noqa
from common import supdraw, Reads, name, OUT

R = Reads()
CLS = [('tandem', 'Tandem'), ('ALE', 'ALE'), ('mixed', 'Mixed'), ('other', 'Other')]
OTHER = {'intronic_compatible', 'terminal_exon_nonUTR'}
cls4 = lambda c: 'other' if c in OTHER else c
PAIRCOL = {'both_tandem': ORD3[0], 'both_ALE': ORD3[1], 'other_mixed': ORD3[2]}
H = 4.5+128.0+11.3+89.0+11.3+132.4+4.5
F = supdraw.Fig(H, name('strat'))
X0, XR = 4.5, PAGE_W-4.5

# ------------------------------------------------------------------ a  class matrix
def panel_a(y0, h):
    x, y, w, _ = F.panel('a', X0, y0, 170.0, h)
    m = R.df('notes/2026-10-10-utr-event-stratification/classification/species_class_matrix.tsv')
    core = R.df('notes/2026-10-10-utr-event-stratification/classification/shared_proximal_core_all_thresholds.tsv')
    core = core[np.isclose(core.threshold, 0.02)]
    M = {}; C = {}
    for r in m.itertuples(): M[(cls4(r.human_event_class), cls4(r.mouse_event_class))] = M.get((cls4(r.human_event_class), cls4(r.mouse_event_class)), 0)+int(r.n_pairs)
    for r in core.itertuples(): C[(cls4(r.human_event_class), cls4(r.mouse_event_class))] = C.get((cls4(r.human_event_class), cls4(r.mouse_event_class)), 0)+1
    assert sum(M.values()) == 751 and sum(C.values()) == 27
    cx0, cw, ch, top = x+42.0, 32.0, 21.0, y+22.0
    F.text(cx0+2*cw, y+0.6, 'Mouse class', 7.5, ha='center')
    for j, (k, lab) in enumerate(CLS):
        F.text(cx0+j*cw+cw/2, y+11.0, lab, 7.0, ha='center')
        F.text(x+38.5, top+j*ch+ch/2-3.3, lab, 7.0, ha='right')
    F.vtext(x+4.0, top+2*ch, 'Human class', 7.5)
    lmax = math.log10(315)
    for i, (hk, _) in enumerate(CLS):
        for j, (mk, _) in enumerate(CLS):
            n = M.get((hk, mk), 0); c = C.get((hk, mk), 0)
            fc = supdraw.seq(math.log10(n)/lmax if n else 0.0)
            F.rect(cx0+j*cw, top+i*ch, cw-1.0, ch-1.0, fc=fc)
            lab = f'{n}' + (f' ({c})' if c else '')
            F.text(cx0+j*cw+(cw-1.0)/2, top+i*ch+ch/2-3.8, lab, 7.0, ha='center', color=supdraw.on(fc))
            F.values[f'a {hk}|{mk}'] = dict(pairs=n, shared_proximal_002=c)
    yk = top+4*ch+4.0
    F.cbar_key(cx0+2*cw-0.5, yk, 'Pairs', '1', '315', lambda u: supdraw.seq(u), wbar=44.0)
    F.text(cx0+2*cw-0.5, yk+10.5, '(n) shared proximal, |ΔDWUI| ≥ .02', 7.0, ha='center', color=MUTED)
    return yk+10.5+0.93*7.0

# ------------------------------------------------------------------ b  composition across thresholds
def panel_b(y0, h):
    x0 = X0+170.0+GAP
    x, y, w, _ = F.panel('b', x0, y0, XR-x0, h)
    core = R.df('notes/2026-10-10-utr-event-stratification/classification/shared_proximal_core_all_thresholds.tsv')
    T = sorted(core.threshold.unique())
    bx, sc = x+30.0, Lin(0, 60, x+30.0, x+30.0+168.0)
    F.text(x+20.0, y+0.6, 'Shared proximal P–M pairs by event class', 7.5)
    F.text(XR-0.2, y+11.0, 'Strict', 7.0, ha='right', color=MUTED)
    ytop, pitch = y+22.0, 9.0
    for i, t in enumerate(T):
        g = core[np.isclose(core.threshold, t)]
        yy = ytop+i*pitch
        if np.isclose(t, 0.02): F.rect(x+0.0, yy-4.3, XR-x, pitch, fc=BAND)
        F.text(x+25.0, yy-3.3, f'{t:.3f}'.lstrip('0'), 7.0, ha='right')
        cnt = {k: int((g.pair_class == k).sum()) for k in ('both_tandem', 'both_ALE')}
        cnt['other_mixed'] = len(g)-cnt['both_tandem']-cnt['both_ALE']
        acc = 0
        for k in ('both_tandem', 'both_ALE', 'other_mixed'):
            if cnt[k]: F.rect(sc(acc), yy-2.6, sc(acc+cnt[k])-sc(acc), 5.2, fc=PAIRCOL[k])
            acc += cnt[k]
        F.text(sc(acc)+1.6, yy-3.3, str(acc), 7.0)
        strict = int((g.strict_full_axis_sensitivity == True).sum())
        F.text(XR-0.2, yy-3.3, str(strict), 7.0, ha='right', color=MUTED)
        F.values[f'b t={t:.3f}'] = dict(cnt, total=len(g), strict_full_axis=strict)
    yb = ytop+(len(T)-1)*pitch+4.6
    yy = F.xaxis(sc, yb, [0, 20, 40, 60], title='Pairs')
    F.yaxis(Lin(0, 1, yb, ytop-4.6), x+29.0, [], line=False)
    F.vtext(x+4.0, (ytop+yb)/2-4.0, '|ΔDWUI| threshold', 7.5)
    F.key_row(0, yy+2.5, [('block', PAIRCOL['both_tandem'], 'Both tandem'), ('block', PAIRCOL['both_ALE'], 'Both ALE'),
                          ('block', PAIRCOL['other_mixed'], 'Other or mixed')], center=(sc.p0+sc.p1)/2)
    return yy+2.5+7.6

# ------------------------------------------------------------------ c  direction pairing within classes, two nulls
COMBOS = [('proximal', 'proximal'), ('distal', 'distal'), ('distal', 'proximal'), ('proximal', 'distal')]
DC = dict(proximal=PROX, distal=DIST)
def panel_c(y0, h):
    x, y, w, _ = F.panel('c', X0, y0, XR-X0, h)
    d = R.df('notes/2026-10-10-utr-event-stratification/overlap/directional_tests.tsv')
    d = d[np.isclose(d.threshold, 0.02)]
    ytop, pitch = y+30.0, 11.0
    F.text(x+4.5, y+20.0, 'Hs', 7.0, color=SPECIES['human'], ha='center')
    F.text(x+15.5, y+20.0, 'Mm', 7.0, color=SPECIES['mouse'], ha='center')
    for k, (hd, md) in enumerate(COMBOS):
        yy = ytop+k*pitch
        if k % 2 == 0: F.rect(x, yy-pitch/2, XR-x, pitch, fc=BAND)
        F.rect(x+1.7, yy-2.8, 5.6, 5.6, fc=DC[hd]); F.rect(x+12.7, yy-2.8, 5.6, 5.6, fc=DC[md])
    blocks = [('both_tandem', 'full_assay'), ('both_tandem', 'measurement_PAS_capacity'), ('both_ALE', 'full_assay'), ('both_ALE', 'measurement_PAS_capacity')]
    b0, gin, gout = x+30.0, 8.0, 16.0
    bw = (XR-b0-2*gin-gout)/4
    PX = [b0, b0+bw+gin, b0+2*bw+gin+gout, b0+3*bw+2*gin+gout]
    for bi, (cls, model) in enumerate(blocks):
        px = PX[bi]
        sc = Lin(0, 20, px, px+bw-34.0)
        g = d[(d.pair_class == cls) & (d.model == model)]
        if bi % 2 == 0:
            n = int(g.measured_orthologs.iloc[0])
            F.text(px+bw+gin/2, y+0.6, f'{"Both tandem-3′ UTR" if cls == "both_tandem" else "Both ALE-compatible"}, {n} pairs', 7.5, ha='center')
        F.text((sc.p0+sc.p1)/2, y+12.0, 'Full assay' if model == 'full_assay' else 'Capacity', 7.0, ha='center', color=MUTED)
        F.text(px+bw, y+12.0, 'n/exp.', 7.0, ha='right', color=MUTED)
        for k, (hd, md) in enumerate(COMBOS):
            r = g[(g.human_direction == hd) & (g.mouse_direction == md)].iloc[0]
            yy = ytop+k*pitch
            q = r.q_primary4 if (hd, md) == ('proximal', 'proximal') else r.q_secondary16_per_threshold
            F.line([sc(r.null_p025), sc(r.null_p975)], [yy, yy], lw='range', color=NULLBAR)
            F.line([sc(r.analytic_expected)]*2, [yy-2.6, yy+2.6], lw='thin', color=MUTED)
            F.marker(sc(r.observed), yy, 'o', 'mark', fc=INK if q < 0.05 else 'white', ec=INK)
            F.text(px+bw, yy-3.3, f'{int(r.observed)}/{r.analytic_expected:.1f}', 7.0, ha='right')
            F.values[f'c {cls}|{model}|{hd}-{md}'] = dict(observed=int(r.observed), expected=float(r.analytic_expected),
                                                         null=[float(r.null_p025), float(r.null_p975)], q=float(q),
                                                         q_family='q_primary4' if (hd, md) == ('proximal', 'proximal') else 'q_secondary16_per_threshold')
        yb = ytop+3*pitch+pitch/2
        yy = F.xaxis(sc, yb, [0, 10, 20], title='Pairs' if bi == 0 else None)
    F.key_row(0, yy+2.5, [('marker', ('o', INK), 'q < .05'), ('marker', ('o', 'white', INK), 'q ≥ .05'),
                          ('line', (NULLBAR, None, 'range'), 'Null 95%'), ('line', (MUTED, None, 'thin'), 'Expected'),
                          ('block', PROX, 'Proximal'), ('block', DIST, 'Distal')], center=(b0+XR)/2)
    return yy+2.5+7.6

# ------------------------------------------------------------------ d  DYNLL2 / Dynll2
HG = 'ENSG00000264364'
STUDIES_D = [('GSE118258', 'GSE118258'), ('GSE178325', '178325_0618'), ('GSE178325', '178325_0809'), ('GSE103221', 'CIP'), ('GSE115943', 'WOT')]
def panel_d(y0, h, wd):
    x, y, w, _ = F.panel('d', X0, y0, wd, h)
    F.text(x+16.0, y+0.4, 'DYNLL2/Dynll2', 8.0, font='BolIta')
    F.text(x+16.0+supdraw.textwidth('DYNLL2/Dynll2', 8.0, 'BolIta')+4.0, y+0.9, 'source to P', 7.5, color=MUTED)
    iv = R.df('notes/2026-10-10-conserved-response/sites/strict_case_intervals.tsv')
    iv = iv[(iv.human_gene_id == HG) & (iv.tolerance == 25)]
    sp = R.df('notes/2026-10-10-conserved-response/sites/reciprocal_predicted_seed_pairs.tsv')
    sp = sp[(sp.human_gene_id == HG) & (sp.tolerance == 25)]
    inc = R.df('notes/2026-10-10-conserved-response/sites/strict_case_interval_inclusion.tsv')
    inc = inc[(inc.human_gene_id == HG) & (inc.tolerance == 25)]
    hu = iv[iv.species == 'human'].iloc[0]; mo = iv[iv.species == 'mouse'].iloc[0]
    assert hu.strand == '+' and mo.strand == '-'
    Lh, Lm = int(hu.end0-hu.start0), int(mo.end0-mo.start0)
    hpos = lambda a, b: ((a+b)/2-hu.start0)/1000            # + strand: distance from the proximal cut
    mpos = lambda a, b: (mo.end0-(a+b)/2)/1000               # - strand
    sx = Lin(0, 1.8, x+30.0, x+150.0)
    yh, ym = y+30.0, y+62.0
    F.text(x+0.2, yh-3.3, 'Human', 7.0); F.text(x+0.2, ym-3.3, 'Mouse', 7.0)
    seeds = sorted({(int(r.human_start0), int(r.human_end0), int(r.mouse_start0), int(r.mouse_end0)) for r in sp.itertuples()})
    assert len(sp) == 13
    for hs, he, ms_, me in seeds:
        F.line([sx(hpos(hs, he)), sx(mpos(ms_, me))], [yh+3.0, ym-3.0], lw='ref', color=AXC)
    for yy, L in ((yh, Lh), (ym, Lm)):
        F.line([sx(-0.0), sx(1.8)], [yy, yy], lw='thin', color=AXC)
        F.rect(sx(0), yy-2.5, sx(L/1000)-sx(0), 5.0, fc=SHADE)
        F.line([sx(0)]*2, [yy-5.0, yy+5.0], lw='data', color=PROX)
        F.line([sx(L/1000)]*2, [yy-5.0, yy+5.0], lw='data', color=DIST)
    for hs, he, ms_, me in seeds:
        F.line([sx(hpos(hs, he))]*2, [yh-2.5, yh+2.5], lw='thin', color=INK)
        F.line([sx(mpos(ms_, me))]*2, [ym-2.5, ym+2.5], lw='thin', color=INK)
    F.text(sx(Lh/1000), yh-13.6, f'{Lh:,} nt', 7.0, ha='right', color=MUTED)
    F.text(sx(Lm/1000), ym+6.6, f'{Lm:,} nt', 7.0, ha='right', color=MUTED)
    yax = ym+19.6
    yy = F.xaxis(sx, yax, [0, 0.5, 1, 1.5], ['0', '0.5', '1', '1.5'], title='Distance from proximal PAS (kb)')
    F.values['d segment'] = dict(human=[int(hu.start0), int(hu.end0)], mouse=[int(mo.start0), int(mo.end0)], human_nt=Lh, mouse_nt=Lm,
                                 pooled_inclusion=dict(human=[float(hu.reference_inclusion), float(hu.target_inclusion)],
                                                       mouse=[float(mo.reference_inclusion), float(mo.target_inclusion)]),
                                 seed_pairs=13, seed_positions=len(seeds))
    # inclusion S -> P per study
    ix0 = x+wd-92.0
    si = Lin(0, 100, ix0+10.0, x+wd-7.0)
    F.text((si.p0+si.p1)/2, y+11.0, 'Segment inclusion', 7.5, ha='center')
    rt, rp = y+23.0, 10.2
    for k, (study, series) in enumerate(STUDIES_D):
        r = inc[(inc.study_id == study) & (inc.series == series)].iloc[0]
        yy_ = rt+k*rp+(3.0 if k >= 3 else 0.0)
        F.rect(ix0, yy_-2.8, 5.6, 5.6, fc=STUDY[series][1])
        a, b = 100*r.reference_inclusion, 100*r.target_inclusion
        F.line([si(a), si(b)], [yy_, yy_], lw='thin', color=MUTED)
        F.marker(si(a), yy_, 'o', 'small', fc=SC['S']); F.marker(si(b), yy_, 'o', 'small', fc=SC['P'])
        F.values[f'd inclusion {series}'] = dict(S=a, P=b, delta_pp=b-a)
    yb = rt+4*rp+3.0+5.0
    F.line([si(50)]*2, [rt-4.0, yb], lw='ref', color=AXC, ls='ref')
    yi = F.xaxis(si, yb, [0, 50, 100], title='Inclusion (%)')
    yk = max(yy, yi)+3.0
    F.key_row(0, yk, [('line', (PROX, None, 'data'), 'Proximal PAS'), ('line', (DIST, None, 'data'), 'Distal PAS'),
                      ('line', (INK, None, 'thin'), 'Seed match'), ('marker', ('o', SC['S']), 'S'), ('marker', ('o', SC['P']), 'P')],
              center=x+wd/2)
    keys = [('block', STUDY[s][1], STUDY[s][0]) for _, s in STUDIES_D]
    F.key_row(0, yk+10.5, keys[:3], center=x+wd/2)
    F.key_row(0, yk+21.0, keys[3:], center=x+wd/2)
    return yk+21.0+7.6

# ------------------------------------------------------------------ e  attrition
STAGES = [('all_contexts', 'Shared orthologs'), ('strict_UTR', 'Strict tandem-3′ UTR'), ('order_eligible', '≥ 2 matched PASs in order'),
          ('mass_eligible', 'Matched usage ≥ .8'), ('driver_judgeable', 'Gain/loss PASs judgeable'),
          ('double_driver_correspondence', 'Gain and loss PASs match'), ('two_study_exact_event_both', '≥ 2 studies per species')]
def panel_e(y0, h, x0):
    x, y, w, _ = F.panel('e', x0, y0, XR-x0, h)
    a = R.df('notes/2026-10-10-conserved-response/sites/event_attrition.tsv')
    a = a[(a.source_population == 'main48') & np.isclose(a.gene_threshold, 0.02) & (a.context_role == 'shared_recurrent')]
    F.text(x+16.0, y+4.0, 'Corresponding switches', 7.5)
    sc = Lin(0, 48, x+0.4, XR-24.0)
    F.text(XR-0.2, y+4.6, '25 | 50 nt', 7.0, ha='right', color=MUTED)
    yt, pitch = y+15.8, 15.4
    for k, (st, lab) in enumerate(STAGES):
        yy = yt+k*pitch
        F.text(x+0.2, yy, lab, 7.0)
        for j, (tol, col) in enumerate(((25, INK), (50, ALT_OTHER))):
            n = int(a[(a.tolerance == tol) & (a.stage == st)].unique_orthologs.iloc[0])
            yb = yy+8.6+j*4.0
            F.rect(sc(0), yb-1.5, sc(48)-sc(0), 3.0, fc=TRACK)
            if n: F.rect(sc(0), yb-1.5, sc(n)-sc(0), 3.0, fc=col)
            F.values[f'e {tol} {st}'] = n
        n25 = int(a[(a.tolerance == 25) & (a.stage == st)].unique_orthologs.iloc[0]); n50 = int(a[(a.tolerance == 50) & (a.stage == st)].unique_orthologs.iloc[0])
        F.text(XR-0.2, yy+7.0, f'{n25} | {n50}', 7.0, ha='right')
    yk = yt+(len(STAGES)-1)*pitch+8.6+4.0+4.0
    F.key_row(0, yk, [('block', INK, '25 nt'), ('block', ALT_OTHER, '50 nt'), ('block', TRACK, 'Shared orthologs')], center=x+w/2)
    return yk+7.6

# ------------------------------------------------------------------ layout
H1, H2, H3 = 128.0, 89.0, 132.4
y1 = 4.5; y2 = y1+H1+GAP; y3 = y2+H2+GAP
ea = panel_a(y1, H1); eb = panel_b(y1, H1)
ec = panel_c(y2, H2)
WD = 285.0
ed = panel_d(y3, H3, WD)
ee = panel_e(y3, H3, X0+WD+GAP)
print('bottoms', round(ea-y1, 1), round(eb-y1, 1), round(ec-y2, 1), round(ed-y3, 1), round(ee-y3, 1), 'page', y3+H3+4.5)
F.save(OUT, R)
