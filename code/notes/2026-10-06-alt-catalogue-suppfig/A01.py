"""Supplementary Figure S31 (working ID A01) — the main analyses repeated with alternative PAS calling.

Narrative: alternative calling (same libraries, expression states, reference and statistics; separate calling, merging,
counting, point internal-priming filter and 50-nt terminal protection) reproduces the route- and test-level results
(a–c), whereas the overlap of individual members shrinks with selection depth (d–f).

Layout v2 (2026-10-06, user request: fewer redundant columns, less empty space, more varied chart forms):
a  slopegraph per contrast (P–S documented source, P–M same library): route mean ΔDWUI, polyApipe -> alternative calling,
   16-setting range behind each end; gene-level agreement between the workflows (Spearman rho, same direction among genes
   with |ΔDWUI| >= .02 in both) printed with each route (former panel b)
b  dodged dot plot of the count-based tests: observed / null mean, polyApipe (upper) and alternative (lower) sub-lines,
   each on its own null 95% range, fill = q < .05 (2026-10-06 respaced: rows 3.6 pt apart within, ~6 pt between)
c  dodged forest plot of the same-PAS rank correlations (95% CIs, same sub-lines as b) and arrow plot of the four cis
   increment tests (q scale)
d  alluvial plot of the k = 3 profile classes: polyApipe class -> label-matched alternative class, plus members selected
   by one workflow only (former 3 x 3 tables and the profile-candidate rows)
e  overlap bars (polyApipe only | both | alternative only, as fractions of the union, centred on the shared part) for
   recurrent genes, site-retention genes and exact predicted sites
f  waffle of the low-shift mouse P–M screen: every gene selected by either workflow, strict-interval genes marked per
   workflow, |ΔDWUI| < .01 genes named

Reads only the pinned comparison tables listed in sources.json (no new tests; display-side summaries only).
"""
import sys, itertools, math
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'2026-10-06-supplement-polyapipe/tools'))
from supdraw import *   # noqa
from v5draw import qfmt  # noqa
BASE = Path(__file__).resolve().parent     # supdraw exports its own HERE

NAME = 'SuppFig.S31'
OUT = BASE/'candidates'/(sys.argv[1] if len(sys.argv) > 1 else 'run1')


class PinnedTables:
    """Reads the comparison tables listed in sources.json directly from the repository, after checking each pinned SHA-256.
    S31 is not part of supplement contract r2 (its inputs are the dual-main concordance tables), so it keeps its own
    reader instead of supdraw.ExtData, which since contract r2 reads through the bundle. Receipt keys are unchanged."""

    def __init__(self, manifest):
        import json
        self.m = json.loads(Path(manifest).read_text()); self.reads = {}

    def path(self, key):
        f = self.m['files'][key]; p = ROOT/f['path']; got = digest(p)
        if got != f['sha256']:
            raise RuntimeError(f'{key}: {p} sha256 {got} differs from the pinned {f["sha256"]}')
        self.reads[key] = dict(path=f['path'], sha256=got)
        return p

    def df(self, key, **kw):
        return pd.read_csv(self.path(key), sep='\t', **kw)

    def receipt(self):
        return dict(external_sources=self.reads, external_manifest_pinned=self.m['pinned'], external_note=self.m['note'])


E = PinnedTables(BASE/'sources.json')
PS, AL = 'polyapipe_common', 'reject_original_point_terminal_p50'
from matplotlib.colors import to_hex as _hex
# polyApipe ink; alternative calling = Figure 2e ramp at 0.55 (#5489AF). The ramp stop #2D6A96 was too close to ink
# (dE 26, L* 43 vs 24) for the 3-pt markers and the keys; 0.55 gives dE 34 (L* 55) and avoids #0072B2 (proximal / R).
CP, CA = INK, _hex(seq(0.55)).upper()
assert CA == '#5489AF', CA                    # house colours only (2026-10-06: teal replaced, then lightened, user request)
AM = '^'                                      # alternative-calling marker: triangle (diamond read too close to the circle)
MK = {PS: ('o', CP), AL: (AM, CA)}
SIG = 0.05
X0, RR = MARGIN, PAGE_W-MARGIN
# overlap colours (d, e, f): one workflow only = light tint of its colour; both = dark blend
# polyApipe only = light slate (supplement overlap greys, e.g. S11); alternative only = light ramp blue; both = ramp end;
# both but other class (d) = mid slate. All already used in S1–S30 / Figure 2e.
# polyApipe only lightened and alternative only darkened along the ramp so the two one-workflow tints separate (dE 16 -> see receipt).
ONLY_P, ONLY_A, BOTH, BOTH_X = ALT_ONE, _hex(seq(0.40)).upper(), SEQ_STOPS[3], ALT_OTHER   # registered tokens, same colours


def mk(F, x, y, src, sig=True, ms='small'):
    m, c = MK[src]
    F.marker(x, y, m, ms, fc=c if sig else 'white', ec=c)


def arrow(F, x0, x1, y, color=MUTED, gap0=2.4, gap1=2.8, head=3.2):
    """Thin arrow from the polyApipe marker (x0) to the alternative-calling marker (x1); omitted when they nearly coincide."""
    if abs(x1-x0) < gap0+gap1+1.5: return False
    s = 1 if x1 > x0 else -1
    a, b = x0+s*gap0, x1-s*gap1
    hl = min(head, abs(b-a))
    if abs(b-a) > hl+0.3:
        F.line([a, b-s*hl], [y, y], lw='thin', color=color)
    F.poly([(b, y), (b-s*hl, y-1.6), (b-s*hl, y+1.6)], fc=color)
    return True


def pct(v, d=1):
    return f'{100*v:.{d}f}%'


def head(F, x, y, s):
    F.text(x, y, s, 7.0, font='Bol')


def dodge(ys, sep, lo, hi):
    """Label positions: keep order, at least sep apart, within [lo, hi]."""
    idx = np.argsort(ys); out = np.array(ys, float)[idx]
    for i in range(1, len(out)): out[i] = max(out[i], out[i-1]+sep)
    if out[-1] > hi:
        out -= out[-1]-hi
        for i in range(len(out)-2, -1, -1): out[i] = min(out[i], out[i+1]-sep)
    out = np.maximum(out, lo)
    res = np.empty_like(out); res[idx] = out
    return res


# ------------------------------------------------------------------------------------------------ data
glob = {s: E.df('global_poly' if s == PS else 'global_alt') for s in (PS, AL)}
eff = E.df('effect')
eff = eff[(eff.metric == 'delta_DWUI') & (eff.threshold == 0.02) & (eff.setting == 'supported50_p08')]
rec = E.df('recurrence'); rec = rec[(rec.metric == 'WPP') & (rec.threshold == 0.02) & (rec.series_matching == True) & (rec.variant == 'primary')]
agx = E.df('apa_gex'); agx = agx[(agx.threshold == 0.02) & (agx.baseline_bins == 5)]
prox = E.df('proximal')
spc = E.df('samepas_corr'); spc = spc[(spc.pair == 'pair_mean') & (spc.axis == 'rank') & (spc.tolerance == 25)]
spq = E.df('samepas_q')
cis = E.df('cis'); cis = cis[cis.primary == True]
rconc = E.df('rec_conc'); rconc = rconc[rconc.scope == 'common_evaluable']
tconc = E.df('traj_conc'); ret = E.df('retention'); sites = E.df('sites'); low = E.df('lowshift')
lowg = E.df('lowshift_genes')
mem = {s: E.df('member_poly' if s == PS else 'member_alt') for s in (PS, AL)}

# ------------------------------------------------------------------------------------------------ layout
H1, H2, H3 = 132.0, 208.4, 166.5                 # page 538.5 pt = 189.97 mm (<= 190 mm series convention)
WB = 228.0; WC = RR-X0-WB-GAP
WD, WE = 172.0, 136.0; WF = RR-X0-WD-WE-2*GAP
TOP0 = 15.8                                    # first content row below the panel-letter clear zone
Y1 = MARGIN; Y2 = Y1+H1+GAP; Y3 = Y2+H2+GAP
F = Fig(Y3+H3+MARGIN, NAME)
LEG = 7.5*Fig.TEXT_H+1.0                       # height of one key row
GHEAD = 9.5                                    # group header row
DOD = 1.8                                      # b, c: polyApipe (upper) / alternative (lower) sub-line offset from the row centre
AXIS_H = TICK+1.0+7.5*Fig.TEXT_H+1.6+8.0*Fig.TEXT_H   # x axis with ticks, labels and title


def key_y(y, h):
    return y+h-LEG


# a ---------------------------------------------------------------------------------------- slopegraphs
x, y, w, h = F.panel('a', X0, Y1, RR-X0, H1)
FACETS = [('P_S', 'true_source', 'P–S, documented source', ['GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT']),
          ('P_M', 'same_library', 'P–M, paired', ['GSE118258', 'GSE150311', 'GSE221739', 'CIP', 'WOT'])]
FG = 16.0
fw = (w-FG)/2
top = y+9.0                                    # column headers share the title row (was y+13 with headers one row lower)
ybot = y+h-LEG-2.0-8.0*Fig.TEXT_H-7.0*Fig.TEXT_H-4.0     # room for the end labels and the ordering line
sy = Lin(0.022, -0.035, top+2.0, ybot)        # upper limit just above the largest 16-setting range (.019)
va = {}
for k, (con, des, title, routes) in enumerate(FACETS):
    fx = x+k*(fw+FG)
    ax = fx+30.0 if k == 0 else fx+6.0
    pcol = x+(k+1)*fw+k*FG-1.5                          # right edge of the same-direction column
    rcol = pcol-textwidth('same dir.', 7.0)-6.0         # right edge of the rho column
    lab0 = rcol-textwidth('.00', 7.0)-6.0-max(textwidth(STUDY[r][0], 7.0) for r in routes)-7.6
    xl = ax+12.0; xr = lab0-16.0
    F.text((ax+pcol)/2, y+0.4, title, 7.5, ha='center')
    F.text(rcol, y+0.9, 'ρ', 7.0, ha='right', color=MUTED)            # bottom-aligned with the 7.5-pt title
    F.text(pcol, y+0.9, 'same dir.', 7.0, ha='right', color=MUTED)
    F.yaxis(sy, ax, [-0.03, -0.02, -0.01, 0, 0.01, 0.02], ['−.03', '−.02', '−.01', '0', '.01', '.02'],
            title='Route mean ΔDWUI' if k == 0 else None, labels_off=(k == 1))
    F.line([ax, xr+2.0], [sy(0)]*2, lw='ref', color=AXC, ls='ref')
    pts = {}
    for r in routes:
        pts[r] = {}
        for s in (PS, AL):
            g = glob[s]; g = g[(g.route == r) & (g.contrast == con) & (g.design == des) & (g.status == 'measured')]
            assert len(g) == 16, (s, r, con, len(g))
            prim = float(g[(g.setting == 'supported50_p08') & (g.background == 'available')]['mean'].iloc[0])
            pts[r][s] = dict(primary=prim, min16=float(g['mean'].min()), max16=float(g['mean'].max()))
    for r in routes:                                    # ranges first (behind), then lines and ends
        col = STUDY[r][1]
        for s, xx in ((PS, xl), (AL, xr)):
            d = pts[r][s]
            F.line([xx]*2, [sy(d['max16']), sy(d['min16'])], lw='range', color=NULLBAR)
    for r in routes:
        col = STUDY[r][1]
        F.line([xl, xr], [sy(pts[r][PS]['primary']), sy(pts[r][AL]['primary'])], lw='data', color=col)
        F.marker(xl, sy(pts[r][PS]['primary']), 'o', 'small', fc=CP)      # ends coloured by workflow (key as in b, c);
        F.marker(xr, sy(pts[r][AL]['primary']), AM, 'small', fc=CA)       # the study is the line colour and label block
    ly = dodge([sy(pts[r][AL]['primary']) for r in routes], 8.6, top, ybot)
    for r, yl in zip(routes, ly):
        ye = sy(pts[r][AL]['primary'])
        F.line([xr+3.0, lab0-2.0], [ye, yl], lw='ref', color=AXC)
        F.study_block(lab0, yl-2.8, r)
        F.text(lab0+7.6, yl-3.4, STUDY[r][0], 7.0)
        e_ = eff[(eff.route == r) & (eff.contrast == con)].iloc[0]
        F.text(rcol, yl-3.4, f'{e_.spearman:.2f}'.replace('0.', '.'), 7.0, ha='right')
        F.text(pcol, yl-3.4, pct(float(e_.event_both_direction_agreement)), 7.0, ha='right')
        pts[r]['agreement'] = dict(common_genes=int(e_.common), spearman=float(e_.spearman), both_responders=int(e_.event_both),
                                   same_direction=float(e_.event_both_direction_agreement))
    yl_ = ybot+3.0
    F.text(xl, yl_, 'polyApipe', 7.0, ha='center'); F.text(xr, yl_, 'Alternative', 7.0, ha='center')
    ordv = {}
    for s in (PS, AL):
        g = glob[s]; g = g[(g.contrast == con) & (g.design == des) & (g.status == 'measured') & (g.route.isin(routes))]
        gaps = [grp[grp.species == 'human']['mean'].min()-grp[grp.species == 'mouse']['mean'].max() for _, grp in g.groupby(['setting', 'background'])]
        ordv[s] = dict(settings=len(gaps), human_above_mouse=int(sum(v > 0 for v in gaps)), min_gap=float(min(gaps)))
    pts['ordering'] = ordv
    hm = [ordv[s]['human_above_mouse'] for s in (PS, AL)]
    msg = ('Human above mouse in all 16 settings, both workflows' if hm == [16, 16] else
           f'Human above mouse in {hm[0]}/16 (polyApipe) and {hm[1]}/16 (alternative) settings')
    F.text((ax+pcol)/2, yl_+8.4, msg, 7.0, ha='center', color=MUTED)
    va[con] = pts
F.key_row(None, key_y(y, h), [('marker', ('o', CP), 'polyApipe'), ('marker', (AM, CA), 'Alternative calling'),
                               ('line', (NULLBAR, None, 'range'), '16 settings (range)')], center=x+w/2)
F.values['a'] = va

# b ---------------------------------------------------------------------------------------- count-based tests (arrows)
x, y, w, h = F.panel('b', X0, Y2, WB, H2)
rows = []
for sp, c in (('human', 'P_S'), ('human', 'P_M'), ('mouse', 'P_S'), ('mouse', 'P_M')):
    d = {}
    for s in (PS, AL):
        r_ = rec[(rec.source == s) & (rec.species == sp) & (rec.contrast == c)].iloc[0]
        d[s] = dict(ratio=r_.observed/r_.null_mean, lo=r_.null_p025/r_.null_mean, hi=r_.null_p975/r_.null_mean, q=float(r_.q_all_24),
                    observed=int(r_.observed), null_mean=float(r_.null_mean))
    rows.append(('Within-species recurrence', f'{sp.capitalize()} {c.replace("_", "–")}', d))
for model, lab in (('response_identity', 'Identity'), ('direction_specific', 'Direction')):
    for sp in ('human', 'mouse'):
        d = {}
        for s in (PS, AL):
            r_ = agx[(agx.source == s) & (agx.species == sp) & (agx.null_model == model)].iloc[0]
            d[s] = dict(ratio=float(r_.enrichment_min), lo=None, hi=None, q=float(r_.q_family), observed=int(r_.observed), null_mean=None)
        rows.append(('APA recurrence, opposite GEX change', f'{lab}, {sp}', d))
for c, cl in (('PM', 'P–M'), ('PS', 'P–S')):
    for model, ml in (('joint_measurement_PAS', 'joint'), ('independent_q5_reference_capacity', 'quintile'),
                      ('independent_range10_reference_capacity', 'range')):
        d = {}
        for s in (PS, AL):
            r_ = prox[(prox.source == s) & (prox.contrast == c) & (prox.model == model)].iloc[0]
            d[s] = dict(ratio=r_.observed/r_.null_mean, lo=r_.null_q025/r_.null_mean, hi=r_.null_q975/r_.null_mean,
                        q=float(r_.q_within_model_two_contrasts), observed=int(r_.observed), null_mean=float(r_.null_mean))
        rows.append(('Shared proximal candidates, by null', f'{cl}, {ml}', d))
labw = max(textwidth(r[1], 7.0) for r in rows)
sc = Lin(np.log(0.25), np.log(3.2), x+labw+8.0, x+w-1.6)
L = lambda v: sc(np.log(max(min(v, 3.2), 0.25)))
F.text((sc.p0+sc.p1)/2, y+5.0, 'Observed / null mean', 7.5, ha='center')
ngr = len({r[0] for r in rows})
pitch = (H2-TOP0-ngr*GHEAD-1.0-AXIS_H-2.0-LEG)/len(rows)
yy = y+TOP0; cur = None; vb = {}; segs = []; first = None
for g, lab, d in rows:
    if g != cur:
        if first is not None: segs.append((first, yy))
        head(F, x+0.2, yy+0.6, g); yy += GHEAD; cur = g; first = yy
    ym = yy+pitch/2
    F.text(x+labw+0.2, ym-3.4, lab, 7.0, ha='right')
    for s, off in ((PS, -DOD), (AL, DOD)):
        if d[s]['lo'] is not None:
            F.line([L(d[s]['lo']), L(d[s]['hi'])], [ym+off]*2, lw='bold', color=NULLBAR)
    for s, off in ((PS, -DOD), (AL, DOD)):
        mk(F, L(d[s]['ratio']), ym+off, s, sig=d[s]['q'] < SIG)
    vb[f'{g} | {lab}'] = d
    yy += pitch
segs.append((first, yy))
for y0_, y1_ in segs:
    F.line([L(1.0)]*2, [y0_, y1_], lw='ref', color=AXC, ls='ref')
F.xaxis(sc, yy+1.0, [np.log(v) for v in (0.25, 0.5, 1, 2)], ['0.25', '0.5', '1', '2'], title='Ratio (log scale)')
F.key_row(None, key_y(y, h), [('marker', ('o', CP), 'polyApipe'), ('marker', (AM, CA), 'Alternative'),
                               ('marker', ('o', 'white', INK), 'q ≥ .05'), ('line', (NULLBAR, None, 'bold'), 'null 95%')],
          center=x+w/2, gap=6.0)
F.values['b'] = vb

# c ---------------------------------------------------------------------------------------- same-PAS and cis (arrows)
x, y, w, h = F.panel('c', X0+WB+GAP, Y2, WC, H2)
ROWS_D = [('PMS', 'PS_rho', 'P–S'), ('PMS', 'MS_rho', 'M–S'), ('PMS', 'PM_rho', 'P–M'), ('PMS', 'PM_minus_PS_rho', 'P–M − P–S'),
          ('PMS', 'PM_minus_MS_rho', 'P–M − M–S'), ('PEM', 'PM_rho', 'P–M'), ('PEM', 'EM_rho', 'E–M'), ('PEM', 'PM_minus_EM_rho', 'P–M − E–M')]
GD = {'PMS': 'Same PAS, Xing2020 pairs (P/M/S)', 'PEM': 'Same PAS, Liu2020–Guo2019–CIP'}
CISR = [('P-S', 'human', 'P–S, h→m'), ('P-S', 'mouse', 'P–S, m→h'), ('P-M', 'human', 'P–M, h→m'), ('P-M', 'mouse', 'P–M, m→h')]
labd = max(textwidth(r[2], 7.0) for r in ROWS_D+[(0, 0, c_[2]) for c_ in CISR])
sd = Lin(-0.3, 0.6, x+labd+8.0, x+w-5.0)
F.text((sd.p0+sd.p1)/2, y+5.0, 'Rank correlation ρ', 7.5, ha='center')
nrow = len(ROWS_D)+len(CISR)
pd_ = (H2-TOP0-3*GHEAD-2*(1.0+AXIS_H)-2.0-2.0-LEG)/nrow
yy = y+TOP0; cur = None; vc = {}; segs = []; first = None
for des, met, lab in ROWS_D:
    if des != cur:
        if first is not None: segs.append((first, yy))
        head(F, x+0.2, yy+0.6, GD[des]); yy += GHEAD; cur = des; first = yy
    ym = yy+pd_/2
    F.text(x+labd+0.2, ym-3.4, lab, 7.0, ha='right')
    dd = {}
    for s, off in ((PS, -DOD), (AL, DOD)):
        r_ = spc[(spc.source == s) & (spc.design == des) & (spc.metric == met)].iloc[0]
        q_ = float(spq[(spq.source == s) & (spq.design == des) & (spq.metric == met)].q_BH_family8.iloc[0])
        F.line([sd.clip(r_.low), sd.clip(r_.high)], [ym+off]*2, lw='thin', color=MK[s][1])
        dd[s] = dict(estimate=float(r_.estimate), low=float(r_.low), high=float(r_.high), n_genes=int(r_.n_genes), q_capacity=q_)
    for s, off in ((PS, -DOD), (AL, DOD)):
        mk(F, sd.clip(dd[s]['estimate']), ym+off, s, sig=dd[s]['q_capacity'] < SIG)
    vc[f'{GD[des]} | {lab}'] = dd
    yy += pd_
segs.append((first, yy))
for y0_, y1_ in segs:
    F.line([sd(0)]*2, [y0_, y1_], lw='ref', color=AXC, ls='ref')
yy = F.xaxis(sd, yy+1.0, [-0.2, 0, 0.2, 0.4, 0.6], ['−.2', '0', '.2', '.4', '.6'], title='ρ (thin lines, 95% CI)')+2.0
head(F, x+0.2, yy+0.6, 'cis increment tests (four prespecified)'); yy += GHEAD
sq = Lin(0.0, 1.5, sd.p0, sd.p1)                      # -log10 q
first = yy
for c, sp, lab in CISR:
    ym = yy+pd_/2
    F.text(x+labd+0.2, ym-3.4, lab, 7.0, ha='right')
    qq = {s: float(cis[(cis.source == s) & (cis.contrast == c) & (cis.discovery == sp)].q_boot_primary_four.iloc[0]) for s in (PS, AL)}
    a_ = arrow(F, sq.clip(-math.log10(qq[PS])), sq.clip(-math.log10(qq[AL])), ym)
    for s in (PS, AL):
        mk(F, sq.clip(-math.log10(qq[s])), ym, s, sig=qq[s] < SIG)
    vc[f'cis | {lab}'] = dict(q={s: qq[s] for s in (PS, AL)}, arrow=a_)
    yy += pd_
F.line([sq(-math.log10(0.05))]*2, [first, yy], lw='ref', color=AXC, ls='ref')
F.xaxis(sq, yy+1.0, [0, -math.log10(0.3), 1.0, -math.log10(0.05)], ['1', '.3', '.1', '.05'], title='q (log scale)')
F.key_row(None, key_y(y, h), [('marker', ('o', CP), 'polyApipe'), ('marker', (AM, CA), 'Alternative'),
                               ('marker', ('o', 'white', INK), 'q ≥ .05')], center=x+w/2, gap=6.0)
F.values['c'] = vc

# d ---------------------------------------------------------------------------------------- alluvial: profile classes
x, y, w, h = F.panel('d', X0, Y3, WD, H3)
CLS = ['C1', 'C2', 'C3']
ve = {}
NODE = 4.0
nl = x+textwidth('none 000', 7.0)+3.0                  # left node x (labels to its left)
nr = x+w-1.5-textwidth('000 none', 7.0)-1.5-NODE        # right node x
F.text(nl+NODE/2, y+0.4, 'polyApipe', 7.0, ha='center', color=MUTED)
F.text(nr+NODE/2, y+0.4, 'Alternative', 7.0, ha='center', color=MUTED)
blk = (h-10.0-LEG-2.0)/2
for k, sp in enumerate(('human', 'mouse')):
    a = mem[PS][mem[PS].species == sp].set_index('gene_id').cluster
    b = mem[AL][mem[AL].species == sp].set_index('gene_id').cluster
    com = a.index.intersection(b.index)
    tab = pd.crosstab(a[com], b[com]).reindex(index=CLS, columns=CLS, fill_value=0)
    best = max(itertools.permutations(CLS), key=lambda p: sum(tab.loc[CLS[i], p[i]] for i in range(3)))
    relab = {best[i]: CLS[i] for i in range(3)}         # alternative label -> matched polyApipe label
    bm = b.map(relab)
    tab = pd.crosstab(a[com], bm[com]).reindex(index=CLS, columns=CLS, fill_value=0)
    agree = sum(tab.iloc[i, i] for i in range(3))/len(com)
    tc = tconc[tconc.species == sp].iloc[0]
    assert len(com) == int(tc.common_selected) and abs(agree-float(tc.optimal_label_matched_agreement)) < 1e-9, (sp, len(com), agree)
    onlyP = a.drop(com).value_counts().reindex(CLS, fill_value=0)
    onlyA = bm.drop(com).value_counts().reindex(CLS, fill_value=0)
    nP = {c: int(tab.loc[c].sum()+onlyP[c]) for c in CLS}; nA = {c: int(tab[c].sum()+onlyA[c]) for c in CLS}
    tot_l = sum(nP.values())+int(onlyA.sum()); tot_r = sum(nA.values())+int(onlyP.sum())
    by0 = y+10.0+k*(blk+1.0)
    hd = 9.0
    F.text(x+w-1.5, by0, f'{sp.capitalize()} · {len(com)} in both · {pct(agree)} same class · ARI {float(tc.ARI):.2f}'.replace(' 0.', ' .'), 7.0, ha='right')
    gapn = 2.0
    avail = blk-hd-3*gapn-1.0
    u = avail/max(tot_l, tot_r)
    # left stack: C1..C3 (polyApipe), then 'alternative only'; right stack: C1..C3 (alternative, matched), then 'polyApipe only'
    yl = {}; yc = by0+hd
    for c in CLS: yl[c] = yc; yc += nP[c]*u+gapn
    yl['none'] = yc
    yr = {}; yc = by0+hd
    for c in CLS: yr[c] = yc; yc += nA[c]*u+gapn
    yr['none'] = yc
    flows = []   # (left key, right key, n, colour)
    for c in CLS:
        for c2 in CLS:
            flows.append((c, c2, int(tab.loc[c, c2]), BOTH if c == c2 else BOTH_X))
        flows.append((c, 'none', int(onlyP[c]), ONLY_P))
    for c2 in CLS:
        flows.append(('none', c2, int(onlyA[c2]), ONLY_A))
    curL = dict(yl); curR = dict(yr)
    # order the outgoing/incoming flows so ribbons do not cross more than needed
    for lk, rk, n, col in flows:
        if n == 0: continue
        hgt = n*u
        a0, b0 = curL[lk], curR[rk]; curL[lk] += hgt; curR[rk] += hgt
        xs = np.linspace(nl+NODE, nr, 24); t = (xs-xs[0])/(xs[-1]-xs[0]); sm = t*t*(3-2*t)
        top_ = a0+(b0-a0)*sm; bot_ = top_+hgt
        F.poly(list(zip(xs, top_))+list(zip(xs[::-1], bot_[::-1])), fc=col, alpha=ALPHA['full'])
    for side, ys_, nx, cnt_ in (('L', yl, nl, {**nP, 'none': int(onlyA.sum())}), ('R', yr, nr, {**nA, 'none': int(onlyP.sum())})):
        for key in CLS+['none']:
            hh = cnt_[key]*u
            F.rect(nx, ys_[key], NODE, max(hh, 0.6), fc=INK if key != 'none' else MUTED)
            lab = key if key != 'none' else 'none'
            txt = f'{lab} {cnt_[key]}'
            if side == 'L': F.text(nx-1.5, ys_[key]+hh/2-3.3, txt, 7.0, ha='right', color=INK if key != 'none' else MUTED)
            else: F.text(nx+NODE+1.5, ys_[key]+hh/2-3.3, f'{cnt_[key]} {lab}', 7.0, color=INK if key != 'none' else MUTED)
    ve[sp] = dict(common=len(com), polyapipe_candidates=int(tc.selected_polyapipe), alternative_candidates=int(tc.selected_apadetect),
                  matched_alternative_labels=list(best), table=tab.values.tolist(), polyapipe_only=onlyP.to_dict(),
                  alternative_only=onlyA.to_dict(), agreement=agree, ARI=float(tc.ARI))
    assert sum(nP.values()) == int(tc.selected_polyapipe) and sum(nA.values()) == int(tc.selected_apadetect)
F.key_row(None, key_y(y, h), [('block', BOTH, 'same class'), ('block', BOTH_X, 'other class')], center=x+w/2, gap=6.0)
F.values['d'] = ve

# e ---------------------------------------------------------------------------------------- overlap bars
x, y, w, h = F.panel('e', X0+WD+GAP, Y3, WE, H3)
rows_e = []
for sp, c in (('human', 'P_S'), ('human', 'P_M'), ('mouse', 'P_S'), ('mouse', 'P_M')):
    r_ = rconc[(rconc.species == sp) & (rconc.contrast == c)].iloc[0]
    rows_e.append(('Recurrent genes', f'{sp.capitalize()} {c.replace("_", "–")}', int(r_.polyapipe), int(r_.apadetect), int(r_.shared_genes)))
for sp, c in (('human', 'SP'), ('human', 'MP'), ('mouse', 'SP'), ('mouse', 'MP')):
    r_ = ret[(ret.species == sp) & (ret.contrast == c) & (ret.category == 'miRNA')].iloc[0]
    r2 = ret[(ret.species == sp) & (ret.contrast == c) & (ret.category == 'RBP')].iloc[0]
    assert (r_.polyapipe, r_.apadetect_p50, r_.shared) == (r2.polyapipe, r2.apadetect_p50, r2.shared)
    rows_e.append(('Site-retention genes', f'{sp.capitalize()} P–{c[0]}', int(r_.polyapipe), int(r_.apadetect_p50), int(r_.shared)))
for sp, c in (('human', 'SP'), ('human', 'MP'), ('mouse', 'SP'), ('mouse', 'MP')):
    r_ = sites[(sites.species == sp) & (sites.contrast == c)].iloc[0]
    rows_e.append(('Exact predicted sites', f'{sp.capitalize()} P–{c[0]}', int(r_.polyapipe), int(r_.apadetect_p50), int(r_.shared)))
labe = max(textwidth(r[1], 7.0) for r in rows_e)
jw = textwidth('.00', 7.0)
se = Lin(-1.0, 1.0, x+labe+5.0, x+w-1.5-jw-5.0)
F.text((se.p0+se.p1)/2, y+5.0, 'Share of the union', 7.5, ha='center')
F.text(x+w-1.5, y+5.0, 'J', 7.5, ha='right')
ng = len({r[0] for r in rows_e})
pe = (h-TOP0-ng*GHEAD-1.0-AXIS_H+8.0*Fig.TEXT_H+1.6-2.0-LEG)/len(rows_e)    # axis without title
yy = y+TOP0; cur = None; vv = {}
for g, lab, nPo, nAo, nS in rows_e:
    if g != cur:
        head(F, x+0.2, yy+0.6, g); yy += GHEAD; cur = g
    U = nPo+nAo-nS
    fp, fs, fa = (nPo-nS)/U, nS/U, (nAo-nS)/U
    ym = yy+pe/2; bh = pe-2.2
    F.text(x+labe+0.2, ym-3.4, lab, 7.0, ha='right')
    F.rect(se(-fs/2-fp), ym-bh/2, se(-fs/2)-se(-fs/2-fp), bh, fc=ONLY_P)
    F.rect(se(-fs/2), ym-bh/2, se(fs/2)-se(-fs/2), bh, fc=BOTH)
    F.rect(se(fs/2), ym-bh/2, se(fs/2+fa)-se(fs/2), bh, fc=ONLY_A)
    F.text(x+w-1.5, ym-3.4, f'{nS/U:.2f}'.replace('0.', '.'), 7.0, ha='right')
    vv[f'{g} | {lab}'] = dict(polyapipe=nPo, alternative=nAo, shared=nS, union=U, jaccard=nS/U)
    yy += pe
F.line([se(0)]*2, [y+TOP0+GHEAD, yy], lw='ref', color=AXC, ls='ref')
F.xaxis(se, yy+1.0, [-1, -0.5, 0, 0.5, 1], ['1', '.5', '0', '.5', '1'])
F.key_row(None, key_y(y, h), [('block', ONLY_P, 'polyApipe'), ('block', BOTH, 'both'), ('block', ONLY_A, 'alternative')],
          center=x+w/2, gap=5.0)
F.values['e'] = vv

# f ---------------------------------------------------------------------------------------- waffle: low-shift screen
x, y, w, h = F.panel('f', X0+WD+WE+2*GAP, Y3, WF, H3)
gP = lowg[lowg.source == PS].set_index('gene_id'); gA = lowg[lowg.source == AL].set_index('gene_id')
assert len(gP) == int(low[low.stage == 'low_abs_DWUI_lt002_MPRO_ge01'].polyapipe.iloc[0])
assert len(gA) == int(low[low.stage == 'low_abs_DWUI_lt002_MPRO_ge01'].apadetect_p50.iloc[0])
both_ids = gP.index.intersection(gA.index)
assert len(both_ids) == int(low[low.stage == 'low_abs_DWUI_lt002_MPRO_ge01'].shared.iloc[0])


def flags(gid):
    sp_ = gP.loc[gid] if gid in gP.index else None; sa_ = gA.loc[gid] if gid in gA.index else None
    strictP = sp_ is not None and sp_.strict_region_pairs > 0; strictA = sa_ is not None and sa_.strict_region_pairs > 0
    finP = strictP and sp_.max_abs_DWUI < 0.01; finA = strictA and sa_.max_abs_DWUI < 0.01
    sym = (sp_ if sp_ is not None else sa_).gene_symbol
    return dict(symbol=sym, strictP=bool(strictP), strictA=bool(strictA), finalP=bool(finP), finalA=bool(finA))


FL = {g: flags(g) for g in gP.index.union(gA.index)}
nsP = sum(v['strictP'] for v in FL.values()); nsA = sum(v['strictA'] for v in FL.values())
nfP = sum(v['finalP'] for v in FL.values()); nfA = sum(v['finalA'] for v in FL.values())
for stage, col in (('strict_recurrent_interval', None), ('both_abs_DWUI_lt001', None)):
    r_ = low[low.stage == stage].iloc[0]
    key = 'strict' if stage.startswith('strict') else 'final'
    nP_ = nsP if key == 'strict' else nfP; nA_ = nsA if key == 'strict' else nfA
    ns_ = sum(v['strictP'] and v['strictA'] for v in FL.values()) if key == 'strict' else sum(v['finalP'] and v['finalA'] for v in FL.values())
    assert (nP_, nA_, ns_) == (int(r_.polyapipe), int(r_.apadetect_p50), int(r_.shared)), (stage, nP_, nA_, ns_)
BLOCKS = [('polyApipe only', [g for g in gP.index if g not in both_ids], ONLY_P), ('both', list(both_ids), BOTH),
          ('alternative only', [g for g in gA.index if g not in both_ids], ONLY_A)]
NC, CELL = 4, 8.4                              # 4 columns: larger cells, the 43-gene block is 11 rows
bw_ = NC*CELL
bgap = (w-1.0-3*bw_)/2
F.text(x+w/2, y+5.0, 'Low-net-shift P–M screen', 7.5, ha='center')   # same height as the e title
HEADS = [('polyApipe', 'only'), ('Both', ''), ('Alternative', 'only')]
hy = y+TOP0                                    # block headers start below the letter clear zone
wy0 = hy+2*7.4+2.4
vf = dict(blocks={}, strict=dict(polyapipe=nsP, alternative=nsA), final=dict(polyapipe=nfP, alternative=nfA))
named = []; bottoms = []
for b_, (bn, ids, col) in enumerate(BLOCKS):
    bx = x+0.5+b_*(bw_+bgap)
    l1, l2 = HEADS[b_]
    F.text(bx+bw_/2, hy, l1, 7.0, ha='center')
    F.text(bx+bw_/2, hy+7.4, f'{l2} · {len(ids)}' if l2 else f'{len(ids)}', 7.0, ha='center', color=MUTED)
    ids = sorted(ids, key=lambda g: (-(FL[g]['finalP'] or FL[g]['finalA']), -(FL[g]['strictP'] or FL[g]['strictA']), FL[g]['symbol']))
    mc = 'white' if col in (BOTH, ONLY_A) else INK   # markers on coloured cells are white
    for i, g in enumerate(ids):
        cx_ = bx+(i % NC)*CELL; cy_ = wy0+(i//NC)*CELL
        f_ = FL[g]
        F.rect(cx_+0.4, cy_+0.4, CELL-0.8, CELL-0.8, fc=col)
        if f_['strictP'] and f_['strictA']:
            F.marker(cx_+CELL/2-2.1, cy_+CELL/2, 'o', 'small', fc=mc)
            F.marker(cx_+CELL/2+2.1, cy_+CELL/2, AM, 'small', fc=mc)
        elif f_['strictP']:
            F.marker(cx_+CELL/2, cy_+CELL/2, 'o', 'small', fc=mc)
        elif f_['strictA']:
            F.marker(cx_+CELL/2, cy_+CELL/2, AM, 'small', fc=mc)
        if f_['finalP'] or f_['finalA']: named.append((b_, g))
    vf['blocks'][bn] = [FL[g]['symbol'] for g in ids]
    bottoms.append(wy0+math.ceil(len(ids)/NC)*CELL)
# names of the strict |ΔDWUI| < .01 genes directly under their own block, with the workflow's marker
for b_, (bn, ids, col) in enumerate(BLOCKS):
    bx = x+0.5+b_*(bw_+bgap)
    nm = sorted([FL[g] for bb, g in named if bb == b_], key=lambda v: v['symbol'])
    for j, f_ in enumerate(nm):
        ms_ = [m for m, k in (('o', 'finalP'), (AM, 'finalA')) if f_[k]]
        tw = textwidth(f_['symbol'], 7.0); gw = 4.2*len(ms_)
        nx = bx+bw_/2-(gw+tw)/2; yy_ = bottoms[b_]+2.6+j*7.4
        for k_, m in enumerate(ms_):
            F.marker(nx+1.5+4.2*k_, yy_+3.4, m, 'small', fc=CP if m == 'o' else CA)
        F.text(nx+gw, yy_, f_['symbol'], 7.0, font='RegIta')
vf['named_final'] = [dict(symbol=FL[g]['symbol'], polyapipe=FL[g]['finalP'], alternative=FL[g]['finalA']) for _, g in named]
STAGES = [('low_abs_DWUI_lt002_MPRO_ge01', '|ΔDWUI| < .02'), ('strict_recurrent_interval', 'strict interval'), ('both_abs_DWUI_lt001', '|ΔDWUI| < .01')]
vf['stages'] = {}                              # stage counts kept in values (readable from the waffle; table dropped)
for st, lab_ in STAGES:
    r_ = low[low.stage == st].iloc[0]
    vf['stages'][lab_] = dict(polyapipe=int(r_.polyapipe), alternative=int(r_.apadetect_p50), shared=int(r_.shared))
F.text(x+w/2, key_y(y, h)-19.0, 'markers, strict interval', 7.0, ha='center', color=MUTED)
F.text(x+w/2, key_y(y, h)-10.6, 'names, also |ΔDWUI| < .01', 7.0, ha='center', color=MUTED)
F.key_row(None, key_y(y, h), [('marker', ('o', CP), 'polyApipe'), ('marker', (AM, CA), 'Alternative')], center=x+w/2, gap=6.0)
F.values['f'] = vf

F.notes += ['Alternative calling = reject_original_point_terminal_p50 (point internal-priming filter + 50-nt terminal protection); '
            'polyApipe = polyapipe_common v2_primary. No new tests; ratios, ranges, overlaps and label matching are display-side.',
            'a: primary end = supported50_p08, available background; vertical bars = range over 4 settings x 4 backgrounds; '
            'labels dodged at the alternative end; rho and same direction from effect-concordance (supported50_p08, |ΔDWUI| >= .02).',
            'b: recurrence uses series-matched nulls (q_all_24); APA-GEX uses the 5-bin baseline (q over the 8-test family) and has no '
            'null range; proximal overlap q within model over the two contrasts. Null ranges (upper bar polyApipe, lower bar alternative) '
            'are null 2.5/97.5 percentiles divided by the null mean. Each row is dodged: polyApipe on the upper sub-line, alternative on '
            'the lower, each marker on its own null range (b) or CI (c, same-PAS rows); cis rows keep the arrow form '
            '(arrows omitted where the two markers nearly coincide).',
            'c: capacity-conditioned null q (8-test family); CI = gene-bootstrap 95% CI of the two-pair mean (25 nt, rank axis), clipped '
            'at -0.3; cis rows: q over the four prespecified primary increment tests, on a log q axis.',
            'd: alternative-calling labels permuted to the matching that maximises agreement on the genes selected by both; ribbon '
            'heights proportional to genes; n and agreement checked against trajectory-concordance.tsv.',
            'e: bars are fractions of the union of the two sets, centred on the shared part; J = Jaccard index. Recurrent genes on the '
            'common measurable background; site-retention miRNA and RBP gene sets identical (checked).',
            'f: one square per gene selected by either workflow at the first screen step (|ΔDWUI| < .02, MPRO >= .10); markers = strict '
            'recurrent interval in that workflow (circle polyApipe, triangle alternative; white on coloured cells); names under their block = strict genes with |ΔDWUI| < .01 in the workflow(s) marked before the name (Figure 5c); counts checked against '
            'low-shift-candidate-concordance.tsv.']
# display-side summaries quoted in the legend (computed here, not hard-coded)
_ag = [d_['agreement'] for con_ in F.values['a'].values() for r_, d_ in con_.items() if r_ != 'ordering']
F.notes.append(f"Legend summaries: medians over the {len(_ag)} route x contrast pairs, Spearman rho "
               f"{np.median([g_['spearman'] for g_ in _ag]):.2f} and same direction {100*np.median([g_['same_direction'] for g_ in _ag]):.1f}%; "
               f"same-class agreement human {100*F.values['d']['human']['agreement']:.1f}% and mouse {100*F.values['d']['mouse']['agreement']:.1f}%; "
               f"filled markers q < {SIG:.2f}.".replace('0.05', '.05'))
F.save(OUT, E)
