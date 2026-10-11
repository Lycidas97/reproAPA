"""Cross-species main figure for the v11 numbering (Figure 3), 2026-10-10.
a = direction pairing at |dDWUI| >= 0.02 (former candidate 4b, moved left, unchanged content)
b = threshold scan of the same pairing (notes/2026-10-10-ortholog-threshold-scan), full_assay model, q72 per model
c-e = production Figure 4 c-e, vector crop at the original position (letters unchanged)
Former 4a lives in Supplementary Figure S13a (10-08 decision). Output: out-v11/ only; Fig/ is never written.
Also assembles the renumbered main-figure set in out-v11/ (copies, no redraw).
Usage: python3 compose_v11.py"""
import sys, json, shutil
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compose as C
from compose import supdraw, Lin, PROX, DIST, SPECIES, INK, MUTED, AXC, NULLBAR, BAND, W, FIG, ROOT, read, crop, place, new_page, sha
from supdraw import STUDY, SC, MISSING

OUT2 = HERE/'out-v11'; OUT2.mkdir(exist_ok=True)
C.OUT = OUT2
H = 481.12
TOP = 4.4                      # top of the first-row panels (as in the production figure)
YAX = 129.5                    # shared x-axis line of the first row
GAP = 11.3


def pairing(F, x0):
    """Former candidate 4b, identical layout shifted to x0."""
    dx = x0-303.6
    y0, w, h = TOP, 449.04-303.6, 154.9
    F.panel('A', x0, y0, w, h)
    d = read('notes/2026-10-10-recurrence-utr-ortholog/ortholog_3x3.tsv')
    d = d[(d.model == 'full_assay') & np.isclose(d.threshold, 0.02) & d.tested]
    sf = Lin(0, 3.5, 346+dx, 393+dx)
    combos = [('proximal', 'proximal'), ('distal', 'distal'), ('distal', 'proximal'), ('proximal', 'distal')]
    hdr = y0+3.0
    F.text(323.3+dx, hdr, 'Hs', 7.0, color=SPECIES['human'], ha='center')
    F.text(334.3+dx, hdr, 'Mm', 7.0, color=SPECIES['mouse'], ha='center')
    F.text((sf.p0+sf.p1)/2, hdr, 'Obs./exp.', 7.5, ha='center')
    F.text(428.0+dx, hdr, 'n/exp.', 7.0, ha='right', color=MUTED)
    F.text(449.0+dx, hdr, 'q', 7.0, ha='right', color=MUTED)
    rh, gap, ytop = 12.5, 6.0, 22.0
    y = ytop
    for gi, ct in enumerate(('P_S', 'P_M')):
        gy = y
        if gi == 0: F.rect(x0, gy, w, 4*rh, fc=BAND)
        F.text(x0, gy+2*rh-3.3, ct.replace('_', '–'), 7.0)
        for hd, md in combos:
            yc = y+rh/2
            F.rect(320.5+dx, yc-2.8, 5.6, 5.6, fc=PROX if hd == 'proximal' else DIST)
            F.rect(331.5+dx, yc-2.8, 5.6, 5.6, fc=PROX if md == 'proximal' else DIST)
            r = d[(d.contrast == ct) & (d.human_direction == hd) & (d.mouse_direction == md)].iloc[0]
            e = r.analytic_expected; qv = r.q_eight_per_threshold_model
            F.line([sf(r.null_p025/e), sf(r.null_p975/e)], [yc, yc], lw='range', color=NULLBAR)
            F.marker(sf(r.observed/e), yc, 'o', 'mark', fc=INK if qv < 0.05 else 'white', ec=INK)
            F.text(428.0+dx, yc-3.3, f'{int(r.observed)}/{e:.1f}', 7.0, ha='right')
            F.text(449.0+dx, yc-3.3, f'{qv:.3f}'.lstrip('0'), 7.0, ha='right')
            F.values[f'a {ct} human {hd} mouse {md}'] = dict(obs=int(r.observed), expected=e, range=[r.null_p025, r.null_p975], q8=qv,
                                                            orthologs=int(r.measured_orthologs))
            y += rh
        y += gap
    yb = y-gap
    assert abs(yb+1.5-YAX) < 1e-9
    F.line([sf(1)]*2, [ytop, yb], lw='ref', color=AXC, ls='ref')
    yy = F.xaxis(sf, YAX, [0, 1, 2, 3], title='Obs./expected')
    ky = yy+3.5
    F.key_row(0, ky, [('block', PROX, 'Proximal'), ('block', DIST, 'Distal'), ('marker', ('o', INK), 'q < 0.05')], center=x0+w/2)
    return ky


SERIES = [(('proximal', 'proximal'), PROX, None, 'data'),
          (('distal', 'distal'), DIST, None, 'data'),
          (('distal', 'proximal'), MUTED, 'dash', 'thin')]
TICKS = [0.01, 0.03, 0.05]


def scan(F, x0, x1, ky):
    y0 = TOP; h = 154.9
    F.panel('B', x0, y0, x1-x0, h)
    d = read('notes/2026-10-10-ortholog-threshold-scan/ortholog_3x3_threshold_scan.tsv')
    d = d[(d.model == 'full_assay') & d.tested]
    left, right, sgap = x0+18.0, x1-8.0, 20.0
    ws = (right-left-sgap)/2
    sy = Lin(0, 4.0, YAX, y0+20.0)
    F.values['b rule'] = 'points with analytic expected count < 1 omitted; filled = q72 < 0.05 (full_assay, 72 tests per model)'
    for i, ct in enumerate(('P_S', 'P_M')):
        px0 = left+i*(ws+sgap)
        sx = Lin(0.01, 0.05, px0, px0+ws)
        F.lineplot(sx, sy, label=ct)
        g = d[d.contrast == ct]
        n = int(g.measured_orthologs.iloc[0])
        F.text(px0+ws/2, y0+3.0, f'{ct.replace("_", "–")}, {n:,} pairs', 7.5, ha='center')
        F.line([sx.p0, sx.p1], [sy(1)]*2, lw='ref', color=AXC, ls='ref')
        for (hd, md), col, ls, lw in SERIES:
            s = g[(g.human_direction == hd) & (g.mouse_direction == md)].sort_values('threshold')
            s = s[s.analytic_expected >= 1.0]
            xs = [sx(t) for t in s.threshold]; ys = [sy(f) for f in s.fold]
            F.line(xs, ys, lw=lw, color=col, ls=ls)
            for t, f, q, o, e in zip(s.threshold, s.fold, s.q_scan_72_per_model, s.observed, s.analytic_expected):
                F.marker(sx(t), sy(f), 'o', 'small', fc=col if q < 0.05 else 'white', ec=col)
                F.values[f'b {ct} human {hd} mouse {md} t={t:.3f}'] = dict(obs=int(o), expected=float(e), fold=float(f), q72=float(q))
        F.yaxis(sy, px0, [1, 2, 3, 4], title='Obs./expected' if i == 0 else None, labels_off=i > 0)
        F.xaxis(sx, YAX, TICKS, ['0.01', '0.03', '0.05'], title='|ΔDWUI| threshold')
    F.key_row(0, ky, [('line', (PROX, None, 'data'), 'Same proximal'), ('line', (DIST, None, 'data'), 'Same distal'),
                      ('line', (MUTED, 'dash', 'thin'), 'Hs distal, Mm proximal')], center=x0+(x1-x0)/2)


LEFT_R = 229.1                 # right edge of the production left column (d starts at 240.4)
C_TOP, E_TOP, E_BOT = 170.0, 296.8, 477.2
STRAT = 'notes/2026-10-10-utr-event-stratification/overlap/primary_proximal_overlap.tsv'
CLASSES = [('both_tandem', 'Tandem 3′ UTR', PROX, None), ('both_ALE', 'Alt. last exon', PROX, 'dash')]


def panel_c(F):
    """Shared P–M proximal pairs by event class (both species tandem-3′ UTR-compatible, both ALE-compatible), full-assay null."""
    x0, y0 = 4.5, C_TOP
    h = E_TOP-GAP-C_TOP
    F.panel('C', x0, y0, LEFT_R-x0, h)
    d = read(STRAT); d = d[d.model == 'full_assay']
    sf = Lin(0, 4.0, 62.0, 122.0)
    hdr = y0+3.0
    F.text((sf.p0+sf.p1)/2, hdr, 'Obs./exp.', 7.5, ha='center')
    F.text(147.0, hdr, 'n/exp.', 7.0, ha='right', color=MUTED)
    rh, ytop = 31.0, y0+22.0
    for i, (cls, lab, col, ls) in enumerate(CLASSES):
        g = d[d.pair_class == cls]
        yc = ytop+i*rh+rh/2
        if i == 0: F.rect(x0, ytop, 147.0-x0+1.0, rh, fc=BAND)
        F.text(x0, yc-7.4, lab, 7.0)
        F.text(x0, yc+0.6, f'{int(g.measured_orthologs.iloc[0])} pairs', 7.0, color=MUTED)
        r2 = g[np.isclose(g.threshold, 0.02)].iloc[0]; r5 = g[np.isclose(g.threshold, 0.05)].iloc[0]
        e = r2.analytic_expected
        F.line([sf(r2.null_p025/e), sf(r2.null_p975/e)], [yc, yc], lw='range', color=NULLBAR)
        if r5.analytic_expected >= 1.0:
            F.marker(sf(r5.observed/r5.analytic_expected), yc, 'D', 'mark', fc=INK if r5.q_scan36 < 0.05 else 'white', ec=INK)
        F.marker(sf(r2.observed/e), yc, 'o', 'mark', fc=INK if r2.q_primary4 < 0.05 else 'white', ec=INK)
        F.text(147.0, yc-3.3, f'{int(r2.observed)}/{e:.1f}', 7.0, ha='right')
        F.values[f'c forest {cls}'] = dict(obs02=int(r2.observed), exp02=float(e), range02=[float(r2.null_p025), float(r2.null_p975)],
                                          q4=float(r2.q_primary4), obs05=int(r5.observed), exp05=float(r5.analytic_expected), q36_05=float(r5.q_scan36))
    yb = ytop+2*rh
    F.line([sf(1)]*2, [ytop, yb], lw='ref', color=AXC, ls='ref')
    yy = F.xaxis(sf, yb+1.5, [0, 1, 2, 3, 4], title='Obs./expected')
    px0, px1 = 162.0, LEFT_R-8.0
    sx = Lin(0.01, 0.05, px0, px1); sy = Lin(0, 4.0, yb+1.5, y0+20.0)
    F.lineplot(sx, sy, label='class scan', standalone=False)
    F.text((px0+px1)/2, hdr, 'Same proximal', 7.5, ha='center')
    F.line([sx.p0, sx.p1], [sy(1)]*2, lw='ref', color=AXC, ls='ref')
    for cls, lab, col, ls in CLASSES:
        g = d[(d.pair_class == cls) & (d.analytic_expected >= 1.0)].sort_values('threshold')
        F.line([sx(t) for t in g.threshold], [sy(o/e) for o, e in zip(g.observed, g.analytic_expected)], lw='data', color=col, ls=ls)
        for t, o, e, q in zip(g.threshold, g.observed, g.analytic_expected, g.q_scan36):
            F.marker(sx(t), sy(o/e), 'o', 'small', fc=col if q < 0.05 else 'white', ec=col)
            F.values[f'c scan {cls} t={t:.3f}'] = dict(obs=int(o), expected=float(e), fold=float(o/e), q36=float(q))
    F.yaxis(sy, px0, [1, 2, 3, 4], title=None)
    F.xaxis(sx, yb+1.5, TICKS, ['0.01', '0.03', '0.05'], title='|ΔDWUI| threshold')
    F.key_row(0, yy+3.5, [('marker', ('o', 'white', INK), '0.02'), ('marker', ('D', 'white', INK), '0.05'),
                          ('line', (PROX, None, 'data'), 'Tandem'), ('line', (PROX, 'dash', 'data'), 'Alt. last exon'),
                          ('marker', ('o', INK), 'q < 0.05')],
              center=x0+(LEFT_R-x0)/2)


# ------------------------------------------------------------------------------------------------- e  PDLIM7 / Pdlim7
HG, MG = 'ENSG00000196923', 'ENSMUSG00000021493'
TXS = {'human': ('ENST00000355572', 'ENST00000355841'), 'mouse': ('ENSMUST00000069968', 'ENSMUST00000046246')}
ROUTES_E = ['GSE118258', '178325_0618', '178325_0809', 'GSE150311', 'GSE221739', 'CIP', 'WOT']
LS_STATE = dict(S='dash', P=None, M='dashdot')


def panel_e(F, B):
    x0, y0 = 4.5, E_TOP
    F.panel('E', x0, y0, LEFT_R-x0, E_BOT-y0)
    F.text(x0+16.0, y0+0.4, 'PDLIM7/Pdlim7', 8.0, font='BolIta')
    hp = B.df('F4E', 'axis_extensions/homolog_site_pairs.tsv.gz')
    hp = hp[(hp.human_gene_id == HG) & np.isclose(hp.tolerance, 25)]
    tx = B.df('D18H', 'geometry/transcript-models.tsv.gz'); ex_ = B.df('D18H', 'supplement_statistics/same-PAS-example-exon-CDS.tsv.gz')
    sites = sorted({(int(a.split('|')[-1]), int(b.split('|')[-1])) for a, b in zip(hp.human_id, hp.mouse_id)}, reverse=True)
    assert len(sites) == 2 and all(r.split('|')[3] == '-' for r in hp.human_id) and all(r.split('|')[3] == '-' for r in hp.mouse_id)
    (hprox, mprox), (hdist, mdist) = sites          # minus strand: the proximal cut has the larger coordinate
    prox = dict(human=hprox, mouse=mprox); dist = dict(human=hdist, mouse=mdist)

    def usage(pair, design, side):
        g = hp[(hp.pair == pair) & (hp.design == design)]
        return {int(r[f'{"human" if side == "h" else "mouse"}_id'].split('|')[-1]): (float(r[f'reference_{side}']), float(r[f'target_{side}']))
                for _, r in g.iterrows()}

    PL = [('GSE118258', 'H1_C', 'h', 'human'), ('CIP', 'H1_C', 'm', 'mouse'), ('WOT', 'H1_W', 'm', 'mouse')]
    ex0, egap = x0+28.0, 8.0
    ew = (LEFT_R-1.5-ex0-2*egap)/3
    XR = (-1.0, 9.5)
    ybar, yshort, ylong = y0+11.0, y0+21.5, y0+30.5
    sy = Lin(0, 100, y0+83.0, y0+37.5)
    F.text(x0+0.2, yshort-3.4, 'Short', 7.0); F.text(x0+0.2, ylong-3.4, 'Long', 7.0)
    EU = {}
    for i, (route, pair, side, sp) in enumerate(PL):
        px = ex0+i*(ew+egap)
        sx = Lin(XR[0], XR[1], px, px+ew)
        F.lineplot(sx, sy, label=route, standalone=False)
        F.rect(px, ybar, ew, 4.5, fc=STUDY[route][1])
        P0 = prox[sp]
        dk = lambda c, P0=P0: (P0-c)/1000
        for tid, yr in zip(TXS[sp], (yshort, ylong)):
            t = tx[(tx.transcript_id == tid) & (tx.kind == 'transcript')].iloc[0]
            lo, hi = max(XR[0], dk(t.end0)), min(XR[1], dk(t.start0))
            F.line([sx(lo), sx(hi)], [yr, yr], lw='thin', color=AXC)
            seg = ex_[ex_.transcript_id == tid]
            lx = seg[seg.kind == 'exon'].sort_values('start0').iloc[0]      # minus strand: the 3′-most exon
            tcol = PROX if tid == TXS[sp][0] else DIST
            for kind, hh, col in (('exon', 3.0, NULLBAR), ('CDS', 5.4, INK)):
                for v in seg[seg.kind == kind].itertuples():
                    a, b = dk(v.end0), dk(v.start0)
                    if b < XR[0] or a > XR[1]: continue
                    a, b = max(a, XR[0]), min(b, XR[1])
                    term = v.start0 >= lx.start0 and v.end0 <= lx.end0
                    F.rect(sx(a), yr-hh/2, max(sx(b)-sx(a), 0.6), hh, fc=tcol if term else col)
            F.values.setdefault('e last exons', {})[tid] = dict(start0=int(lx.start0), end0=int(lx.end0), kb=[dk(lx.end0), dk(lx.start0)])
        kb = (prox[sp]-dist[sp])/1000
        for v in (0.0, kb):
            F.line([sx(v)]*2, [yshort-4.0, sy.p0], lw='ref', color=INK, ls='ref')
        ps, pm = usage(pair, 'qualified_P_S', side), usage(pair, 'qualified_P_M', side)
        cur = dict(S=ps[prox[sp]][0], P=ps[prox[sp]][1], M=pm[prox[sp]][0]); dis = dict(S=ps[dist[sp]][0], P=ps[dist[sp]][1], M=pm[dist[sp]][0])
        for st in ('S', 'M', 'P'):
            frac = 100*dis[st]/(cur[st]+dis[st])
            F.line([sx(XR[0]), sx(0), sx(0), sx(kb), sx(kb), sx(XR[1])], [sy(100), sy(100), sy(frac), sy(frac), sy(0), sy(0)],
                   lw='data', color=SC[st], ls=LS_STATE[st])
            EU[f'{route}|{st}'] = dict(proximal=cur[st], distal=dis[st], distal_pct=frac)
        F.yaxis(sy, px-1.5, [0, 50, 100], title='Distal (%)' if i == 0 else None, labels_off=i > 0)
        F.xaxis(sx, sy.p0+1.5, [0, 4, 8], ['0', '4', '8'])
    yy = sy.p0+1.5+2+1+0.93*7.5
    F.text(ex0+(3*ew+2*egap)/2, yy+1.6, 'Distance from proximal PAS (kb)', 8.0, ha='center')
    yy += 1.6+0.93*8.0
    mem = B.df('F4E', 'cross_trace/member_endpoint_effects.tsv.gz')
    mem = mem[(mem.human_gene_id == HG) & (mem.setting == 'supported50_p08')]
    h3b = hp[(hp.pair == 'H3b_C') & (hp.design == 'qualified_P_S')]
    cw = (LEFT_R-ex0-6*1.0-3.0)/7
    CX = lambda j: ex0+j*(cw+1.0)+(3.0 if j >= 5 else 0.0)
    hb = yy+4.0
    for j, r in enumerate(ROUTES_E): F.rect(CX(j), hb, cw, 4.0, fc=STUDY[r][1])
    cells = []
    for k, (con, lab) in enumerate((('P_S', 'P–S'), ('P_M', 'P–M'))):
        yt = hb+5.5+k*12.5
        F.text(x0+0.2, yt+2.4, lab, 7.0)
        for j, route in enumerate(ROUTES_E):
            v = None
            if route == '178325_0809':
                if con == 'P_S' and len(h3b): v = float(h3b.native_h.iloc[0])
            else:
                m = mem[(mem.contrast == con) & (mem.route == route)]
                if len(m) and m.status.iloc[0] == 'measured': v = float(m.delta.iloc[0])
            if v is None:
                F.rect(CX(j), yt, cw, 11.5, fc=MISSING)
            else:
                fc = supdraw.diverge(v, 0.3); F.rect(CX(j), yt, cw, 11.5, fc=fc)
                F.text(CX(j)+cw/2, yt+2.2, '−' if v < 0 else '+', 7.5, ha='center', color=supdraw.on(fc))
            cells.append(dict(contrast=con, route=route, value=v))
    yk = hb+5.5+2*12.5+3.0
    F.key_row(0, yk, [('line', (SC['S'], 'dash', 'data'), 'S'), ('line', (SC['P'], None, 'data'), 'P'),
                      ('line', (SC['M'], 'dashdot', 'data'), 'M'), ('block', PROX, 'Short last exon'), ('block', DIST, 'Long last exon')],
              center=x0+(LEFT_R-x0)/2)
    yk2 = yk+10.5
    tw = supdraw.textwidth('Native ΔDWUI', 7.5)
    wbar, gap = 50.0, 6.0
    tot = tw+gap+supdraw.textwidth('−.3', 7.0)+1.8+wbar+1.8+supdraw.textwidth('+.3', 7.0)+14.0+7.6+supdraw.textwidth('Not available', 7.5)
    xs0 = x0+(LEFT_R-x0)/2-tot/2
    F.text(xs0, yk2, 'Native ΔDWUI', 7.5)
    cb0 = xs0+tw+gap+supdraw.textwidth('−.3', 7.0)+1.8
    F.cbar(cb0, yk2+1.4, wbar, '−.3', '+.3', lambda u: supdraw.diverge(-0.3+0.6*u, 0.3))
    xm = cb0+wbar+1.8+supdraw.textwidth('+.3', 7.0)+14.0
    F.rect(xm, yk2+1.0, 5.6, 5.6, fc=MISSING); F.text(xm+7.6, yk2, 'Not available', 7.5)
    keys = [('block', STUDY[r][1], STUDY[r][0]) for r in ROUTES_E]
    F.key_row(0, yk2+10.5, keys[:4], center=x0+(LEFT_R-x0)/2)
    F.key_row(0, yk2+21.0, keys[4:], center=x0+(LEFT_R-x0)/2)
    F.values['e'] = dict(human_sites=[hprox, hdist], mouse_sites=[mprox, mdist], transcripts=TXS, usage=EU, cells=cells,
                         note='usage renormalised within the two matched PASs (as production 4e); native ΔDWUI below −0.3 drawn at the ramp end')


def figure3():
    src = FIG/'Figure.4.pdf'
    DREG = (240.0, 166.0, W, H)
    dst = new_page(H)
    place(dst, H, crop(src, H, DREG), DREG, DREG[0], DREG[1])
    F = supdraw.Fig(H, 'Figure.3', W=W)
    B = supdraw.Data()
    xa = 4.5; wa = 449.04-303.6
    ky = pairing(F, xa)
    scan(F, xa+wa+GAP, 449.04, ky)
    panel_c(F)
    panel_e(F, B)
    rec = dict(source=str(src.relative_to(ROOT)), source_sha=sha(src), kept=[dict(region=DREG, dy=0)], new=['a', 'b', 'c', 'e'],
               former={'a': 'candidate 4b (out/Figure.4.pdf)', 'b': 'new: threshold scan', 'c': 'new: event-class stratification',
                       'd': 'production 4d (crop)', 'e': 'production 4e redrawn with transcript models',
                       'production 4a': 'Supplementary Figure S13a', 'production 4c': 'to the Supplement'},
               reads=dict(C.READS), bundle=B.receipt())
    C.finish(dst, F, 'Figure.3', H, rec)


def assemble():
    """Renumbered set for the v11 narrative (copies only)."""
    m = {'Figure.1.pdf': FIG/'Figure.1.pdf',
         'Figure.2.pdf': HERE/'out'/'Figure.3.pdf',
         'Figure.4.pdf': ROOT/'notes/2026-10-08-narrative-figures/out/Figure.4.pdf'}
    for k, s in m.items():
        shutil.copyfile(s, OUT2/k)
        C.RECEIPT.setdefault('set', {})[k] = dict(source=str(Path(s).relative_to(ROOT)), sha256=sha(s))
    C.RECEIPT['set']['Figure.3.pdf'] = dict(source='composed here', sha256=sha(OUT2/'Figure.3.pdf'))


if __name__ == '__main__':
    figure3(); assemble()
    C.RECEIPT['production_replaced'] = False
    (OUT2/'receipt.json').write_text(json.dumps(C.RECEIPT, indent=1, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)))
