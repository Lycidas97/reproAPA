"""Divergence-narrative main-figure candidates (2026-10-10). Spec: DESIGN.zh.md in this folder.
Kept panels are vector crops of the production PDFs in Fig/ (moved, never scaled, no builder rerun); replaced panels are
drawn with supdraw design tokens and checked. Output: out/ only; Fig/ is never written.
Usage: python3 compose.py [3] [4]"""
import sys, json, hashlib, subprocess, tempfile
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'notes/2026-10-06-main-figures-polyapipe/builders/tools-local'))
sys.path.insert(0, str(ROOT/'notes/2026-10-06-supplement-polyapipe/tools'))
from v5draw import clear, prune_xobjects, INK, MUTED, AXC, NULLBAR, BAND, textwidth   # noqa
import pikepdf, pdfplumber
import supdraw
from supdraw import Lin, PROX, DIST, STUDY, SPECIES, on
OUT = HERE/'out'; OUT.mkdir(exist_ok=True)
FIG = ROOT/'revision-backups/20261011-before-v11-apply/Fig'
W = 453.5433
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
READS, RECEIPT = {}, {}


def read(rel, **kw):
    p = ROOT/rel; READS[rel] = sha(p)
    return pd.read_csv(p, sep='\t', **kw)


def crop(src, H, region, extra_clear=()):
    x0, y0, x1, y1 = region
    blk = pikepdf.open(src)
    if extra_clear: clear(blk, H, list(extra_clear))
    clear(blk, H, [(-20, -20, x0-0.05, H+20), (x1+0.05, -20, W+20, H+20), (-20, -20, W+20, y0-0.05), (-20, y1+0.05, W+20, H+20)])
    prune_xobjects(blk)
    p = blk.pages[0]; box = [x0, H-y1, x1, H-y0]; p.MediaBox = box; p.CropBox = box
    return blk


def place(dst, Hd, blk, region, nx, ny):
    x0, y0, x1, y1 = region
    dst.pages[0].add_overlay(blk.pages[0], pikepdf.Rectangle(nx, Hd-(ny+y1-y0), nx+x1-x0, Hd-ny))


def new_page(H):
    p = pikepdf.new(); p.add_blank_page(page_size=(W, H)); return p


def finish(dst, F, name, H, rec):
    chk = F.check()
    tmp = Path(tempfile.mkdtemp()); F.t.save(tmp/'overlay.pdf')
    ov = pikepdf.open(tmp/'overlay.pdf')
    dst.pages[0].add_overlay(ov.pages[0], pikepdf.Rectangle(0, 0, W, H))
    out = OUT/f'{name}.pdf'; dst.save(out)
    try:
        from pdfclean import finalize
        finalize(out, label=name, log=OUT/'sizes.tsv')
    except Exception as e:
        rec['pdfclean_error'] = repr(e)
    subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(out), str(OUT/name)], check=True, capture_output=True)
    rec.update(pdf=str(out.relative_to(ROOT)), sha256=sha(out), page_pt=[W, H], panels_new=F.panels,
               checks={k: chk[k] for k in ('outside', 'text_overlaps', 'letter_zone', 'colour_unnamed', 'line_aspect_bad', 'sizes',
                                             'slack', 'internal_vertical_strip', 'internal_horizontal_strip')},
               tokens=supdraw.TOKENS_USED, values=F.values)
    RECEIPT[name] = rec
    print(name, 'H', H, '| outside', len(chk['outside']), 'overlaps', len(chk['text_overlaps']), 'letter', len(chk['letter_zone']),
          'colour', len(chk['colour_unnamed']), 'sizes', chk['sizes'])
    print('  slack', chk['slack']); print('  vstrip', chk['internal_vertical_strip']); print('  hstrip', chk['internal_horizontal_strip'])
    for k in ('outside', 'text_overlaps', 'letter_zone', 'colour_unnamed'):
        if chk[k]: print('  ', k, json.dumps(chk[k], default=str)[:900])


def fmt(v, d=2, sign=False):
    s = f'{v:+.{d}f}' if sign else f'{v:.{d}f}'
    return s.replace('-', '−')


def name_of(route):
    return STUDY[route][0]


# =============================================================================================== Figure 3
def fig3a(F):
    x0, y0, w, h = 4.5, 5.1, 222.3, 123.0
    F.panel('A', x0, y0, w, h)
    rec = read('notes/2026-10-10-directional-recurrence/directional_recurrence.tsv')
    comp = read('notes/2026-10-10-directional-recurrence/direction_composition.tsv')
    groups = [('human', 'P_S', 'Human P–S'), ('mouse', 'P_S', 'Mouse P–S'), ('human', 'P_M', 'Human P–M'), ('mouse', 'P_M', 'Mouse P–M')]
    sf = Lin(-1, 3, 80, 140)                     # log2 obs/null
    L2 = lambda v: sf(np.log2(max(v, 0.5)))
    sc = Lin(0, 100, 180, 219)
    rh, gap, ytop = 9.6, 3.0, 22.1
    F.text((sf.p0+sf.p1)/2, y0+3.0, 'Obs./null', 7.5, ha='center')
    F.text(172.5, y0+3.0, 'Obs./mean', 7.5, ha='right')
    F.text((sc.p0+sc.p1)/2, y0+3.0, 'Distal, %', 7.5, ha='center')
    yb = ytop+4*2*rh+3*gap
    for gi, (sp, ct, lab) in enumerate(groups):
        gy = ytop+gi*(2*rh+gap)
        if gi % 2 == 0: F.rect(x0, gy, w, 2*rh, fc=BAND)
        F.text(x0, gy+rh-3.3, lab, 7.0)
        for ri, d in enumerate(('proximal', 'distal')):
            yc = gy+ri*rh+rh/2
            F.text(46.0, yc-3.3, d, 7.0, color=PROX if d == 'proximal' else DIST)
            r2 = rec[(rec.species == sp) & (rec.contrast == ct) & (rec.direction == d) & np.isclose(rec.threshold, 0.02)].iloc[0]
            r5 = rec[(rec.species == sp) & (rec.contrast == ct) & (rec.direction == d) & np.isclose(rec.threshold, 0.05)].iloc[0]
            F.line([L2(r2.null_p025/r2.null_mean), L2(r2.null_p975/r2.null_mean)], [yc, yc], lw='range', color=NULLBAR)
            F.marker(L2(r5.fold_enrichment), yc, 'D', 'mark', fc=INK)
            F.marker(L2(r2.fold_enrichment), yc, 'o', 'mark', fc=INK)
            F.text(172.5, yc-3.3, f'{int(r2.observed)}/{r2.null_mean:.0f}', 7.0, ha='right')
            F.values[f'3a {sp} {ct} {d}'] = dict(obs02=int(r2.observed), null_mean02=r2.null_mean, fold02=r2.fold_enrichment,
                                                 range02=[r2.null_p025, r2.null_p975], fold05=r5.fold_enrichment, q02=r2.q_threshold_8)
        c = comp[(comp.species == sp) & (comp.contrast == ct) & np.isclose(comp.threshold, 0.02)].iloc[0]
        cy = gy+rh
        F.line([sc(100*c.null_distal_fraction_p025), sc(100*c.null_distal_fraction_p975)], [cy, cy], lw='range', color=NULLBAR)
        F.marker(sc(100*c.observed_distal_fraction), cy, 'o', 'mark', fc=INK)
        F.values[f'3a composition {sp} {ct}'] = dict(obs=c.observed_distal_fraction, null=[c.null_distal_fraction_p025, c.null_distal_fraction_p975], q=c.q_threshold_4)
    F.line([sf(0)]*2, [ytop, yb], lw='ref', color=AXC, ls='ref')
    ya = yb+1.5
    F.xaxis(sf, ya, [-1, 0, 1, 2, 3], ['0.5', '1', '2', '4', '8'])
    F.xaxis(sc, ya, [0, 50, 100], ['0', '50', '100'])
    F.key_row(0, ya+11.2, [('marker', ('o', INK), '|ΔDWUI| ≥ 0.02'), ('marker', ('D', INK), '0.05'),
                           ('line', (NULLBAR, None, 'range'), '95% perm. range')], center=x0+w/2)


def p_s_effects():
    q = read('notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/global/qualified_date_effects.tsv.gz')
    q = q[(q.setting == 'supported50_p08') & (q.contrast == 'P_S')]
    g = q.groupby(['route', 'gene_id']).APA_delta.mean().reset_index()
    b = read('notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/global/background_gene_membership.tsv.gz')
    b = b[(b.setting == 'supported50_p08') & (b.contrast == 'P_S') & (b.background == 'common_one_to_one_ortholog_all_settings')]
    s = read('notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/global/global_direction_summary.tsv')
    s = s[(s.setting == 'supported50_p08') & (s.contrast == 'P_S') & (s.background == 'common_one_to_one_ortholog_all_settings')]
    R = {}
    for r in ['GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT']:
        x = g[(g.route == r) & g.gene_id.isin(b[b.route == r].gene_id)].APA_delta.values
        row = s[s.route == r].iloc[0]
        assert len(x) == row.n and abs(x.mean()-row['mean']) < 1e-9 and (x <= -0.02).sum() == row.proximal_002 and (x >= 0.02).sum() == row.distal_002
        pr, di = x[x <= -0.02], x[x >= 0.02]
        R[r] = dict(n=len(x), mean=x.mean(), pos=x[x > 0].sum()/len(x), neg=x[x < 0].sum()/len(x), fp=len(pr)/len(x), fd=len(di)/len(x),
                    np_=len(pr), nd=len(di), mp=pr.mean(), md=di.mean(), pq=np.percentile(pr, [10, 25, 50, 75, 90]).tolist(), dq=np.percentile(di, [10, 25, 50, 75, 90]).tolist())
        assert abs(R[r]['pos']-row.positive_contribution) < 1e-9 and abs(R[r]['neg']-row.negative_contribution) < 1e-9
    hum = {k: (R['GSE118258'][k]+(R['178325_0618'][k]+R['178325_0809'][k])/2)/2 for k in ('mean', 'pos', 'neg', 'fp', 'fd', 'mp', 'md')}
    mou = {k: (R['CIP'][k]+R['WOT'][k])/2 for k in ('mean', 'pos', 'neg', 'fp', 'fd', 'mp', 'md')}
    # study-equal quantile (Vincent) averages: each percentile averaged across studies, Guan2022 series averaged first
    for k in ('pq', 'dq'):
        hum[k] = ((np.array(R['GSE118258'][k])+(np.array(R['178325_0618'][k])+np.array(R['178325_0809'][k]))/2)/2).tolist()
        mou[k] = ((np.array(R['CIP'][k])+np.array(R['WOT'][k]))/2).tolist()
    return R, hum, mou


def fig3c(F, y0):
    x0, w = 4.5, 444.5
    R, hum, mou = p_s_effects()
    s1, s2, s3 = Lin(-0.04, 0.04, 98, 206), Lin(-45, 45, 228, 326), Lin(-0.2, 0.2, 342, 441)
    rh = 10.5; ytop = y0+17.0
    ys = {}; y = ytop
    order = [('GSE118258', 'h'), ('178325_0618', 'h'), ('178325_0809', 'h'), ('gap', 3.5), ('CIP', 'm'), ('WOT', 'm'), ('gap', 6.0),
             ('human', 's'), ('mouse', 's')]
    for k, t in order:
        if k == 'gap': y += t; continue
        ys[k] = y+rh/2; y += rh
    yend = y
    hdr = y0+3.0
    F.text(22.5, hdr, '1,306 orthologs', 7.0, color=MUTED)
    F.text((s1.p0+s1.p1)/2, hdr, 'Mean ΔDWUI split by sign', 7.5, ha='center')
    F.text((s2.p0+s2.p1)/2, hdr, 'Genes with |ΔDWUI| ≥ 0.02, %', 7.5, ha='center')
    F.text((s3.p0+s3.p1)/2, hdr, 'ΔDWUI of these genes', 7.5, ha='center')
    F.vtext(8.6, (ys['GSE118258']+ys['178325_0809'])/2, 'Human', 7.0)
    F.vtext(8.6, (ys['CIP']+ys['WOT'])/2, 'Mouse', 7.0)
    for r in ('GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT'):
        F.study_block(14.5, ys[r]-2.8, r); F.text(22.5, ys[r]-3.3, name_of(r), 7.0)
    for sp in ('human', 'mouse'):
        F.rect(14.5, ys[sp]-2.8, 5.6, 5.6, fc=SPECIES[sp]); F.text(22.5, ys[sp]-3.3, f'{sp.capitalize()}, study mean', 7.0)
    for sc in (s1, s2, s3):
        F.line([sc(0)]*2, [ytop, yend], lw='ref', color=AXC, ls='ref')
    bh = 5.6

    def arms(sc, y, neg, pos, mean=None):
        F.rect(sc(neg), y-bh/2, sc(0)-sc(neg), bh, fc=PROX)
        F.rect(sc(0), y-bh/2, sc(pos)-sc(0), bh, fc=DIST)
        if mean is not None: F.marker(sc(mean), y, 'D', 'small', fc=INK)
    for r in ('GSE118258', '178325_0618', '178325_0809', 'CIP', 'WOT'):
        v = R[r]; y = ys[r]
        arms(s1, y, v['neg'], v['pos'], v['mean'])
        arms(s2, y, -100*v['fp'], 100*v['fd'])
        for q, col in ((v['pq'], PROX), (v['dq'], DIST)):
            F.line([s3.clip(q[0]), s3.clip(q[4])], [y, y], lw='thin', color=col)
            F.line([s3(q[1]), s3(q[3])], [y, y], lw='bold', color=col)
            F.marker(s3(q[2]), y, 'o', 'small', fc='white', ec=col)
        F.marker(s3(v['mp']), y, 'o', 'small', fc=PROX); F.marker(s3(v['md']), y, 'o', 'small', fc=DIST)
        F.values[f'3c {r}'] = v
    for sp, v in (('human', hum), ('mouse', mou)):
        arms(s1, ys[sp], v['neg'], v['pos'], v['mean']); arms(s2, ys[sp], -100*v['fp'], 100*v['fd'])
        for q, col in ((v['pq'], PROX), (v['dq'], DIST)):
            F.line([s3.clip(q[0]), s3.clip(q[4])], [ys[sp]]*2, lw='thin', color=col)
            F.line([s3(q[1]), s3(q[3])], [ys[sp]]*2, lw='bold', color=col)
            F.marker(s3(q[2]), ys[sp], 'o', 'small', fc='white', ec=col)
        F.marker(s3(v['mp']), ys[sp], 'o', 'small', fc=PROX); F.marker(s3(v['md']), ys[sp], 'o', 'small', fc=DIST)
        F.values[f'3c {sp} study mean'] = v
    dpos, dneg = hum['pos']-mou['pos'], hum['neg']-mou['neg']; tot = dpos+dneg
    F.values['3c human-mouse'] = dict(diff=tot, distal_arm=dpos, proximal_arm=dneg, distal_share=dpos/tot, proximal_share=dneg/tot)
    ya = yend+1.5
    F.xaxis(s1, ya, [-0.04, -0.02, 0, 0.02, 0.04], ['−0.04', '−0.02', '0', '+0.02', '+0.04'])
    F.xaxis(s2, ya, [-40, -20, 0, 20, 40], ['40', '20', '0', '20', '40'])
    F.xaxis(s3, ya, [-0.2, -0.1, 0, 0.1, 0.2], ['−0.2', '−0.1', '0', '+0.1', '+0.2'])
    yk = ya+11.2
    F.key_row(0, yk, [('block', PROX, 'Proximal'), ('block', DIST, 'Distal'), ('marker', ('D', INK), 'Net mean'),
                      ('line', (INK, None, 'bold'), 'Quartiles'), ('line', (INK, None, 'thin'), '10th–90th percentile'),
                      ('marker', ('o', INK), 'Mean'), ('marker', ('o', 'white', INK), 'Median')], center=x0+w/2)
    h = yk+7.6-y0
    if 'C' in F.panels: F.panels['C'] = (x0, y0, w, h)
    return h


def fig3d(F, y0, h):
    x0, w = 4.5, 222.3
    F.panel('D', x0, y0, w, h)
    t = read('notes/2026-10-10-recurrence-utr-ortholog/tandem_directional_recurrence.tsv')
    L = read('notes/2026-10-10-recurrence-utr-ortholog/tandem_recurrent_length.tsv')
    t = t[(t.contrast == 'P_S') & np.isclose(t.threshold, 0.02)]
    L = L[(L.contrast == 'P_S') & np.isclose(L.threshold, 0.02) & (L.level == 'route')]
    rh = 12.0
    sf = Lin(0, 3.5, 84, 172)
    F.text(22.5, y0+3.0, 'Strict tandem 3′UTR, recurrent P–S genes', 7.5)
    ytop = y0+21.0; y = ytop; rows = []
    for sp in ('human', 'mouse'):
        for d, word in (('proximal', 'shortening'), ('distal', 'lengthening')):
            rows.append((sp, d, word, y+rh/2)); y += rh
        y += 3.0
    yb1 = y-3.0
    F.text(224.6, ytop-9.0, 'n/exp.', 7.0, ha='right', color=MUTED)
    F.rect(x0, ytop, w, 2*rh, fc=BAND)
    for sp, d, word, yc in rows:
        if d == 'proximal': F.text(x0, yc+rh/2-3.3, sp.capitalize(), 7.0)
        F.text(42.5, yc-3.3, word, 7.0, color=PROX if d == 'proximal' else DIST)
        r = t[(t.species == sp) & (t.direction == d)].iloc[0]; e = r.analytic_expected
        F.line([sf(r.null_p025/e), sf(r.null_p975/e)], [yc, yc], lw='range', color=NULLBAR)
        F.marker(sf(r.observed/e), yc, 'o', 'mark', fc=INK)
        F.text(224.6, yc-3.3, f'{int(r.observed)}/{e:.0f}', 7.0, ha='right')
        F.values[f'3d {sp} {d}'] = dict(obs=int(r.observed), expected=e, range=[r.null_p025, r.null_p975], fold=r.observed/e, q=r.q_eight_per_threshold)
    F.line([sf(1)]*2, [ytop, yb1], lw='ref', color=AXC, ls='ref')
    yy = F.xaxis(sf, yb1+1.5, [0, 1, 2, 3], title='Obs./expected')
    y2 = yy+6.0
    F.text(x0, y2, 'Usage-weighted retained length change, nt', 7.5)
    sl = Lin(-300, 420, 92, 222)
    ytop2 = y2+12.0; y = ytop2; ys = {}
    for r in ('GSE118258', '178325_0618', '178325_0809', 'gap', 'CIP', 'WOT'):
        if r == 'gap': y += 3.0; continue
        ys[r] = y+rh/2; y += rh
    yb2 = y
    F.rect(x0, ytop2, w, 3*rh, fc=BAND)
    for r, yc in ys.items():
        F.study_block(x0, yc-2.8, r); F.text(x0+8.0, yc-3.3, name_of(r), 7.0)
        for d, col in (('distal', DIST), ('proximal', PROX)):
            v = L[(L.route == r) & (L.direction == d)].iloc[0]
            F.line([sl(v.median_delta_nt), sl(v.mean_delta_nt)], [yc, yc], lw='thin', color=col)
            F.marker(sl(v.median_delta_nt), yc, 'o', 'mark', fc='white', ec=col)
            F.marker(sl(v.mean_delta_nt), yc, 'o', 'mark', fc=col)
            F.values[f'3d length {r} {d}'] = dict(n=int(v.n_genes), mean_nt=v.mean_delta_nt, median_nt=v.median_delta_nt)
    F.line([sl(0)]*2, [ytop2, yb2], lw='ref', color=AXC, ls='ref')
    yy = F.xaxis(sl, yb2+1.5, [-200, 0, 200, 400], ['−200', '0', '+200', '+400'])
    F.key_row(0, yy+3.0, [('marker', ('o', DIST), 'Lengthening'), ('marker', ('o', PROX), 'Shortening'),
                          ('marker', ('o', INK), 'Mean'), ('marker', ('o', 'white', INK), 'Median')], center=x0+w/2)
    return yy+3.0+7.6-y0


def fig3():
    src = FIG/'Figure.3.pdf'; Hs = 618.4
    B = (231.0, 0, W, 137.0); E = (231.0, 416.0, W, Hs)
    F0 = supdraw.Fig(Hs+200, 'probe', W=W); hc = fig3c(F0, 139.0)
    for k in supdraw.TOKENS_USED: supdraw.TOKENS_USED[k] = {}
    yd = 139.0+hc+11.3
    dy = yd+0.4-419.4
    hd = 614.3-418.6
    H = round(yd+hd+4.1, 2)
    dst = new_page(H)
    place(dst, H, crop(src, Hs, B), B, B[0], B[1])
    # 3e axis labels: 'ΔGEX'/'ΔAPA' + 'Same'/'Opp.' -> explicit 'GEX same' / 'APA same' (direction agreement of the two studies)
    pg = pdfplumber.open(src).pages[0]
    words = [w for w in pg.extract_words() if w['x0'] > 231 and 425 < w['top'] < 480 and w['text'] in ('ΔGEX', 'ΔAPA', 'Same', 'Opp.')]
    boxes = [(w['x0']-0.3, w['top']-0.3, w['x1']+0.3, w['bottom']+0.3) for w in words]
    place(dst, H, crop(src, Hs, E, boxes), E, E[0], E[1]+dy)
    F = supdraw.Fig(H, 'Figure.3', W=W)
    relab = []
    for w in words:
        if w['text'] in ('ΔGEX', 'ΔAPA'): continue
        col = w['top'] < 440                     # column headers sit above the cells, row labels left of them
        new_s = ('GEX ' if col else 'APA ')+('same' if w['text'] == 'Same' else 'opp.')
        x, ha = ((w['x0']+w['x1'])/2, 'center') if col else (w['x1'], 'right')
        supdraw.TOKENS_USED['size'][7.5] = supdraw.TOKENS_USED['size'].get(7.5, 0)+1
        F.t.text_top(x, w['top']+dy, new_s, 7.5, color=INK, ha=ha)
        relab.append(dict(old=w['text'], new=new_s, x=round(x, 2), top=round(w['top']+dy, 2), ha=ha))
    fig3a(F)
    F.panel('C', 4.5, 139.0, 444.5, hc); fig3c(F, 139.0)
    used = fig3d(F, yd, hd)
    rec = dict(source=str(src.relative_to(ROOT)), source_sha=sha(src), kept=[dict(region=B, dy=0), dict(region=E, dy=round(dy, 2))],
               replaced=['a', 'c', 'd'], d_content_height=round(used, 1), reads=dict(READS),
               e_relabelled=dict(cleared=[[round(v, 2) for v in b] for b in boxes], drawn=relab))
    finish(dst, F, 'Figure.3', H, rec)


# =============================================================================================== Figure 4
def fig4b(F):
    x0, y0, w, h = 303.6, 4.4, 449.04-303.6, 154.9
    F.panel('B', x0, y0, w, h)
    d = read('notes/2026-10-10-recurrence-utr-ortholog/ortholog_3x3.tsv')
    d = d[(d.model == 'full_assay') & np.isclose(d.threshold, 0.02) & d.tested]
    sf = Lin(0, 3.5, 346, 393)
    combos = [('proximal', 'proximal'), ('distal', 'distal'), ('distal', 'proximal'), ('proximal', 'distal')]
    hdr = y0+3.0
    F.text(323.3, hdr, 'Hs', 7.0, color=SPECIES['human'], ha='center')
    F.text(334.3, hdr, 'Mm', 7.0, color=SPECIES['mouse'], ha='center')
    F.text((sf.p0+sf.p1)/2, hdr, 'Obs./exp.', 7.5, ha='center')
    F.text(428.0, hdr, 'n/exp.', 7.0, ha='right', color=MUTED)
    F.text(449.0, hdr, 'q', 7.0, ha='right', color=MUTED)
    rh, gap, ytop = 12.5, 6.0, 22.0
    y = ytop
    for gi, ct in enumerate(('P_S', 'P_M')):
        gy = y
        if gi == 0: F.rect(x0, gy, w, 4*rh, fc=BAND)
        F.text(x0, gy+2*rh-3.3, ct.replace('_', '–'), 7.0)
        for hd, md in combos:
            yc = y+rh/2
            F.rect(320.5, yc-2.8, 5.6, 5.6, fc=PROX if hd == 'proximal' else DIST)
            F.rect(331.5, yc-2.8, 5.6, 5.6, fc=PROX if md == 'proximal' else DIST)
            r = d[(d.contrast == ct) & (d.human_direction == hd) & (d.mouse_direction == md)].iloc[0]
            e = r.analytic_expected; qv = r.q_eight_per_threshold_model
            F.line([sf(r.null_p025/e), sf(r.null_p975/e)], [yc, yc], lw='range', color=NULLBAR)
            F.marker(sf(r.observed/e), yc, 'o', 'mark', fc=INK if qv < 0.05 else 'white', ec=INK)
            F.text(428.0, yc-3.3, f'{int(r.observed)}/{e:.1f}', 7.0, ha='right')
            F.text(449.0, yc-3.3, f'{qv:.3f}'.lstrip('0'), 7.0, ha='right')
            F.values[f'4b {ct} human {hd} mouse {md}'] = dict(obs=int(r.observed), expected=e, range=[r.null_p025, r.null_p975], q8=qv,
                                                             orthologs=int(r.measured_orthologs))
            y += rh
        y += gap
    yb = y-gap
    F.line([sf(1)]*2, [ytop, yb], lw='ref', color=AXC, ls='ref')
    yy = F.xaxis(sf, yb+1.5, [0, 1, 2, 3], title='Obs./expected')
    F.key_row(0, yy+3.5, [('block', PROX, 'Proximal'), ('block', DIST, 'Distal'), ('marker', ('o', INK), 'q < 0.05')], center=x0+w/2)


def fig4():
    src = FIG/'Figure.4.pdf'; H = 481.12
    A = (0, 0, 300.5, 166.0); CDE = (0, 166.0, W, H)
    dst = new_page(H)
    place(dst, H, crop(src, H, A), A, 0, 0); place(dst, H, crop(src, H, CDE), CDE, 0, CDE[1])
    F = supdraw.Fig(H, 'Figure.4', W=W)
    fig4b(F)
    rec = dict(source=str(src.relative_to(ROOT)), source_sha=sha(src), kept=[dict(region=A), dict(region=CDE)], replaced=['b'], reads=dict(READS))
    finish(dst, F, 'Figure.4', H, rec)


if __name__ == '__main__':
    which = sys.argv[1:] or ['3', '4']
    if '3' in which: fig3()
    if '4' in which: fig4()
    prev = json.loads((OUT/'receipt.json').read_text()) if (OUT/'receipt.json').exists() else {}
    prev.update(RECEIPT); prev['production_replaced'] = False
    (OUT/'receipt.json').write_text(json.dumps(prev, indent=1, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)))
