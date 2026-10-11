"""Narrative-lock main-figure candidates (2026-10-08): Figure 3 (compact 3c), new Figure 4 (old 5, new panel b),
new Figure 5 (old 4 c-e). Kept panels are vector crops of the production PDFs in Fig/ (no builder rerun, values and
labels unchanged); new content is drawn with v5draw.Top in the house style. Spec: reviews/2026-10-08-narrative-revision-plan/FIGURE-SPEC.md.
Usage: python3 compose_main.py [3] [4] [5]"""
import sys, json, hashlib, shutil, subprocess, tempfile
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TL = ROOT/'notes/2026-10-06-main-figures-polyapipe/builders/tools-local'
sys.path.insert(0, str(TL))
from v5draw import Top, clear, INK, MUTED, AXC, NULLBAR, NULLW, BAND, textwidth, prune_xobjects   # noqa
import pikepdf, pdfplumber
OUT = HERE/'out'; OUT.mkdir(exist_ok=True)
FIG = ROOT/'Fig'
W = 453.5433
ORANGE, BLUE = '#D55E00', '#0072B2'
CIP_C, WOT_C = '#CF4C00', '#207450'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
RECEIPT = {}


def letters(src):
    pg = pdfplumber.open(src).pages[0]; H = float(pg.height)
    out = {}
    for c in pg.chars:
        if c['size'] >= 9.5 and 'Bol' in c['fontname'] and len(c['text']) == 1 and c['text'].isalpha():
            out[c['text']] = dict(box=(c['x0'], c['top'], c['x1'], c['bottom']), ox=float(c['matrix'][4]), base=H-float(c['matrix'][5]))
    return out


def crop(src, H, region, extra_clear=()):
    """Copy of src keeping only objects that intersect region (top-origin x0,y0,x1,y1); extra_clear boxes removed too."""
    x0, y0, x1, y1 = region
    blk = pikepdf.open(src)
    if extra_clear:
        clear(blk, H, list(extra_clear))
    clear(blk, H, [(-20, -20, x0-0.05, H+20), (x1+0.05, -20, W+20, H+20), (-20, -20, W+20, y0-0.05), (-20, y1+0.05, W+20, H+20)])
    prune_xobjects(blk)
    p = blk.pages[0]; box = [x0, H-y1, x1, H-y0]; p.MediaBox = box; p.CropBox = box
    return blk


def place(dst, Hd, blk, region, nx, ny):
    x0, y0, x1, y1 = region
    dst.pages[0].add_overlay(blk.pages[0], pikepdf.Rectangle(nx, Hd-(ny+y1-y0), nx+x1-x0, Hd-ny))


def new_page(H):
    p = pikepdf.new(); p.add_blank_page(page_size=(W, H)); return p


def finish(dst, top, name, H):
    tmp = Path(tempfile.mkdtemp())
    top.save(tmp/'overlay.pdf')
    ov = pikepdf.open(tmp/'overlay.pdf')
    dst.pages[0].add_overlay(ov.pages[0], pikepdf.Rectangle(0, 0, W, H))
    out = OUT/f'{name}.pdf'; dst.save(out)
    try:
        sys.path.insert(0, str(TL)); from pdfclean import finalize
        finalize(out, label=name, log=OUT/'sizes.tsv')
    except Exception as e:
        print('pdfclean skipped:', repr(e))
    subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(out), str(OUT/name)], check=True, capture_output=True)
    RECEIPT[name] = dict(pdf=str(out.relative_to(ROOT)), sha256=sha(out), page_pt=[W, H])
    return out


def letter(top, L, ox, base):
    top.text(ox, base, L, 11.0, font='Bol')


# ===================================================================================================== house style
# Rules: notes/2026-10-02-supplement-major-revision/VISUAL-PARADIGMS.zh.md (sec. 4-5) and the supdraw design tokens
# (sizes 6.8/7.0/7.5/8.0/11, LW/MS/MEW tokens, named colours, ~3 mm ink gap). New content is drawn with supdraw.Fig so
# the token and layout checks (outside / text overlap / letter zone) apply; kept panels are moved, never scaled.
sys.path.insert(0, str(ROOT/'notes/2026-10-06-supplement-polyapipe/tools'))
import supdraw
from supdraw import Lin, PROX, DIST, SEQ_STOPS, OUTSIDE, TRACK
INKGAP = 8.5          # visible gap between panel inks (~3 mm)
M_ = 4.5              # outer margin


def ink(src, R):
    pg = pdfplumber.open(src).pages[0]; x0, y0, x1, y1 = R; bb = [1e9, 1e9, -1e9, -1e9]
    for L in (pg.chars, pg.rects, pg.lines, pg.curves):
        for o in L:
            if o['x1']-o['x0'] > 400 and o['bottom']-o['top'] > 400: continue
            cx, cy = (o['x0']+o['x1'])/2, (o['top']+o['bottom'])/2
            if x0 <= cx <= x1 and y0 <= cy <= y1:
                bb = [min(bb[0], o['x0']), min(bb[1], o['top']), max(bb[2], o['x1']), max(bb[3], o['bottom'])]
    return bb


def finish_fig(dst, F, name, H, rec):
    chk = F.check()
    bad = {k: chk[k] for k in ('outside', 'text_overlaps', 'letter_zone') if chk.get(k)}
    rec['new_content_checks'] = dict(outside=len(chk['outside']), text_overlaps=len(chk['text_overlaps']), letter_zone=len(chk['letter_zone']),
                                     sizes=chk['sizes'], tokens=supdraw.TOKENS_USED)
    if bad: print(name, 'NEW-CONTENT CHECK', json.dumps(bad, default=str)[:1500])
    finish(dst, F.t, name, H); RECEIPT[name].update(rec)


def swarm(xs, d, half):
    """Deterministic 1-D beeswarm: offsets along y so that markers of diameter d do not overlap (within +-half)."""
    order = np.argsort(xs); placed = []; off = np.zeros(len(xs))
    cand = [0.0]+[s*k*d*0.5 for k in range(1, 60) for s in (1, -1)]
    for i in order:
        best = None
        for c in cand:
            if abs(c) > half: break
            if all((xs[i]-x)**2+(c-y)**2 >= (1.05*d)**2 for x, y in placed if abs(xs[i]-x) < 1.05*d): best = c; break
        if best is None: best = min(cand[:int(2*half/(0.5*d))+1], key=lambda c: sum(1 for x, y in placed if abs(xs[i]-x) < d and abs(c-y) < d))
        off[i] = best; placed.append((xs[i], best))
    return off


# ===================================================================================================== Figure 3
def fig3():
    src = FIG/'Figure.3.pdf'; Hs = 618.4; Ls = letters(src)
    AB = (0, 0, W, 137.0); S1 = (315.6, 145.0, W, 246.5); S2 = (315.6, 262.0, W, 363.5); DE = (0, 414.0, W, Hs)
    i1, i2, ide = ink(src, S1), ink(src, S2), ink(src, DE)
    sw = i1[2]-i1[0]
    s2x = i2[0]; s1x = s2x-sw-INKGAP                       # two scatters side by side, right edge kept
    sy = i1[1]                                              # scatter tops stay at the old 3c scatter top
    KEY_TOP = i1[3]+4.0; C_BOT = KEY_TOP+7.0                # one shared key row closes panel c
    dy_de = (C_BOT+INKGAP+0.5) - ide[1]
    H = round(Hs+dy_de, 2)
    dst = new_page(H)
    rec = dict(source=str(src.relative_to(ROOT)), source_sha=sha(src), blocks=[])
    place(dst, H, crop(src, Hs, AB), AB, 0, 0); rec['blocks'].append(dict(region=AB, to=[0, 0]))
    plus = [(305.0, 155.0, 324.5, 164.6), (305.0, 270.0, 324.5, 279.6)]
    for R, ii, nx in ((S1, i1, s1x), (S2, i2, s2x)):
        place(dst, H, crop(src, Hs, R, plus), R, R[0]+(nx-ii[0]), R[1]+(sy-ii[1]))
        rec['blocks'].append(dict(region=R, to=[R[0]+(nx-ii[0]), R[1]+(sy-ii[1])]))
    rec['blocks'].append(dict(cleared=plus))
    place(dst, H, crop(src, Hs, DE), DE, 0, DE[1]+dy_de); rec['blocks'].append(dict(region=DE, dy=dy_de))
    F = supdraw.Fig(H, 'Figure.3')
    # letters of the moved d, e (in the crop) stay; panel c is registered for the new content
    cx0, cy0 = M_, Ls['c']['box'][1]
    F.panel('c', cx0, cy0, s1x-INKGAP-cx0, C_BOT-cy0)
    t = pd.read_csv(ROOT/'reviews/2026-10-08-plan-a-redesign/background/results/raw-tests.tsv', sep='\t')
    q = pd.read_csv(ROOT/'reviews/2026-10-08-plan-a-redesign/corrected-tests.tsv', sep='\t')
    def row(sp, test):
        r = t[(t.pas_source == 'polyapipe') & (t.species == sp) & (t.test == test)].iloc[0]
        return dict(obs=int(r.observed), mean=float(r.analytic_null_mean), lo=float(r.null_p025), hi=float(r.null_p975),
                    q=float(q[q.test_id == r.test_id]['q_BH_family'].iloc[0]), cand=int(r.candidates))
    V = {sp: {k: row(sp, k) for k in ('bounded_MS', 'S_vector')} for sp in ('mouse', 'human')}
    MEAS = dict(mouse=(395, 420), human=(250, 268))
    assert (V['mouse']['bounded_MS']['obs'], V['human']['bounded_MS']['obs']) == (222, 128)
    assert (V['mouse']['bounded_MS']['cand'], V['human']['bounded_MS']['cand']) == (395, 250)
    QX = s1x-INKGAP                                          # right edge of the q column = panel-c left block edge
    sc = Lin(0.8, 1.3, 72.0, QX-26.0)
    F.text(QX, cy0+0.6, 'q', 8.0, ha='right')
    F.text((sc.p0+sc.p1)/2, cy0+0.6, 'Obs./null mean', 8.0, ha='center')
    LAB = dict(bounded_MS='Residual', S_vector='S profile')
    ROWS = dict(mouse=(cy0+24.0, cy0+39.5), human=(cy0+64.0, cy0+79.5))
    yax = cy0+91.0
    F.line([sc(1.0)]*2, [cy0+12.0, yax], lw='ref', color=AXC, ls='ref')
    for sp, ys in ROWS.items():
        F.text(M_, ys[0]-3.9, sp.capitalize(), 8.0)
        F.text(M_, ys[1]-3.6, f"{V[sp]['bounded_MS']['obs']}/{MEAS[sp][0]}", 7.5, color=MUTED)
        for k, yy in zip(('bounded_MS', 'S_vector'), ys):
            v = V[sp][k]
            F.text(36.0, yy-3.6, LAB[k], 7.5)
            F.line([sc(v['lo']/v['mean']), sc(v['hi']/v['mean'])], [yy, yy], lw='range', color=NULLBAR)
            F.marker(sc(v['obs']/v['mean']), yy, 'o', 'mark', fc=INK)
            qs = (f"{v['q']:.2f}" if v['q'] >= .01 else f"{v['q']:.4f}").lstrip('0')
            F.text(QX, yy-3.6, qs, 7.5, ha='right')
    F.xaxis(sc, yax, [0.8, 1.0, 1.2], labels=['0.8', '1', '1.2'], size=7.5)
    # shared key row: null rows (left block) and the two scatters (right block)
    F.key_row(None, KEY_TOP-0.6, [('marker', ('o', INK), 'Observed'), ('line', (NULLBAR, None, 'range'), '95% null')], size=7.5,
              center=(M_+QX)/2)
    items = [('All measurable', '#C9CFD3', '#C9CFD3', 0.0, 2.6), ('P–M cohort (395)', 'white', '#5F6F78', 0.5, 3.0), ('Aligned (222)', INK, INK, 0.0, 3.0)]
    wk = sum(7.6+textwidth(l, 7.5)+8.0 for l, *_ in items)-8.0
    kx = (s1x+i2[2])/2-wk/2
    for lab, fc, ec, mew, ms in items:                      # exact marker styles of the scatter (not size tokens)
        F.t.marker(kx+2.8, KEY_TOP+3.2, m='o', ms=ms, fc=fc, ec=ec, mew=mew)
        kx += 7.6; kx += textwidth(lab, 7.5); F.t.text_top(kx-textwidth(lab, 7.5), KEY_TOP-0.6, lab, 7.5); kx += 8.0
    rec.update(values_3c=V, measurable=MEAS, layout=dict(scatters=[s1x, s2x, sy], c_bottom=C_BOT, dy_de=dy_de))
    finish_fig(dst, F, 'Figure.3', H, rec)


# ===================================================================================================== new Figure 4 (old 5)
def fig4():
    src = FIG/'Figure.5.pdf'; Hs = 523.38; Ls = letters(src)
    RA, RB, RC, RF = (0, 0, 326.0, 105.0), (0, 105.0, 326.0, 197.0), (326.0, 0, W, 197.0), (0, 380.0, W, Hs)
    iA, iB, iC, iF = ink(src, RA), ink(src, RB), ink(src, RC), ink(src, RF)
    # rows: [a | c] with new b under a; d = old f (Riok3); [new e (Fgfr1op2) | f = old b (sites)]
    b_top, b_bot = iA[3]+INKGAP, iC[3]
    dyF = (b_bot+INKGAP) - iF[1]
    r4_top = iF[3]+dyF+INKGAP
    dxB, dyB = iC[2]-iB[2], r4_top-iB[1]
    H = round(iB[3]+dyB+M_, 2)
    dst = new_page(H)
    rec = dict(source=str(src.relative_to(ROOT)), source_sha=sha(src), blocks=[])
    lbox = lambda L: (Ls[L]['box'][0]-0.6, Ls[L]['box'][1]-0.6, Ls[L]['box'][2]+0.6, Ls[L]['box'][3]+0.6)
    place(dst, H, crop(src, Hs, RA), RA, 0, 0); rec['blocks'].append(dict(panel='a', region=RA))
    place(dst, H, crop(src, Hs, RC), RC, RC[0], 0); rec['blocks'].append(dict(panel='c', region=RC))
    place(dst, H, crop(src, Hs, RF, [lbox('f')]), RF, 0, RF[1]+dyF); rec['blocks'].append(dict(panel='old f -> d', region=RF, dy=dyF))
    place(dst, H, crop(src, Hs, RB, [lbox('b')]), RB, dxB, RB[1]+dyB); rec['blocks'].append(dict(panel='old b -> f', region=RB, dx=dxB, dy=dyB))
    F = supdraw.Fig(H, 'Figure.4')
    letter(F.t, 'd', Ls['f']['ox'], Ls['f']['base']+dyF)
    letter(F.t, 'f', Ls['b']['ox']+dxB, Ls['b']['base']+dyB)
    # ---------------- b: gene-wide shift against recurrent local interval change (mouse P–M), P2 strip with count columns
    NS = ROOT/'reviews/2026-10-08-narrative-logic-chain/netshift-local'
    def load(w):
        g = pd.read_csv(NS/w/'gene_classes.tsv', sep='\t')
        g = g[(g.species == 'mouse') & (g.contrast == 'MP') & (g.cls != 'gene_wide_unavailable')].copy()
        sgn = np.where(g.local_direction.eq('loss'), -1.0, 1.0)
        g['agree'] = np.minimum(g.dDWUI_a*sgn, g.dDWUI_b*sgn)        # gene-wide shift along the interval change, smaller study
        g['k'] = np.where(g.agree >= 0.02, 'cap', np.where(g.agree > 0, 'low', 'opp'))
        cls = g.cls.map(lambda c: 'cap' if c.startswith('captured') else ('low' if c.startswith('same_sign') else 'opp'))
        assert (cls.values == g.k.values).all()                      # the strip position reproduces the class exactly
        return g
    g, g5 = load('polyapipe'), load('p50')
    cnt = {k: int((g.k == k).sum()) for k in ('cap', 'low', 'opp')}; cnt5 = {k: int((g5.k == k).sum()) for k in ('cap', 'low', 'opp')}
    assert len(g) == 122 and cnt == dict(cap=97, low=20, opp=5) and len(g5) == 19 and cnt5 == dict(cap=13, low=5, opp=1), (cnt, cnt5)
    bx0, by0 = M_, b_top
    F.panel('b', bx0, by0, iA[2]-bx0, b_bot-by0)
    LIM = (-0.1, 0.6)
    tr = lambda v: np.sign(v)*np.sqrt(np.abs(v))
    SX0, SX1 = bx0+16.0, iA[2]-6.0                         # one strip across the whole panel width
    sx = lambda v: SX0+(tr(np.clip(v, *LIM))-tr(LIM[0]))/(tr(LIM[1])-tr(LIM[0]))*(SX1-SX0)
    yc, half, yax = by0+38.0, 13.0, by0+56.0
    # header: what is plotted (left, after the letter zone) and the colour key (right)
    F.text(bx0+16.0, by0+0.6, f'Mouse P–M, {len(g)} genes with a strict recurrent interval', 8.0)
    kw = F.key_width([('marker', None, 'Lost'), ('marker', None, 'Gained')], size=7.5)
    F.key_row(SX1-kw, by0+0.6, [('marker', ('o', PROX), 'Lost'), ('marker', ('o', DIST), 'Gained')], size=7.5)
    # class brackets with counts directly above the points they count (house rule 9: counts next to their objects)
    by_ = by0+21.0
    for lo, hi, k, lab in ((LIM[0], 0.0, 'opp', 'Opposite'), (0.0, 0.02, 'low', '<0.02'), (0.02, LIM[1], 'cap', 'Both ≥0.02')):
        x0_, x1_ = sx(lo)+(1.2 if lo > LIM[0] else 0), sx(hi)-(1.2 if hi < LIM[1] else 0)
        F.line([x0_, x0_, x1_, x1_], [by_+2.0, by_, by_, by_+2.0], lw='axis', color=AXC)
        F.text((x0_+x1_)/2, by_-9.4, f'{lab} ({cnt[k]})', 7.5, ha='center')
    F.line([sx(0)]*2, [by_+3.0, yax], lw='ref', color=AXC, ls='ref')
    F.line([sx(0.02)]*2, [by_+3.0, yax], lw='ref', color=AXC, ls='dash')
    pos = {}; lab0 = 'polyApipe'
    xs = np.array([sx(v) for v in g.agree]); off = swarm(xs, 3.6, half)
    for (i, r), x_, o in zip(g.reset_index(drop=True).iterrows(), xs, off):
        F.marker(x_, yc+o, 'o', 'mark', fc=PROX if r.local_direction == 'loss' else DIST)
        pos[(lab0, r.gene_symbol)] = (x_, yc+o)
    F.line([SX0, SX1], [yax, yax], lw='axis', color=AXC)
    for v, lab in [(-0.05, '−0.05'), (0, '0'), (0.02, '0.02'), (0.1, '0.1'), (0.3, '0.3'), (0.5, '0.5')]:
        F.line([sx(v)]*2, [yax, yax+2.0], lw='axis', color=AXC); F.text(sx(v), yax+3.0, lab, 7.5, ha='center')
    F.text((SX0+SX1)/2, yax+12.0, 'Gene-wide ΔDWUI along the interval change (smaller of two studies)', 8.0, ha='center')
    for sym, lx, ly in (('Fgfr1op2', 8.0, -9.5), ('Riok3', -8.0, 8.5)):       # short leaders (house rule 9)
        x_, y_ = pos[(lab0, sym)]
        sx_ = np.sign(lx)
        F.line([x_+sx_*1.6, x_+lx-sx_*0.6], [y_+np.sign(ly)*1.2, y_+ly-np.sign(ly)*0.2], lw='ref', color=MUTED)
        F.text(x_+lx, y_+ly-(7.2 if ly < 0 else 0.4), sym, 7.5, font='RegIta', ha='left' if lx > 0 else 'right')
    # ---------------- e: Fgfr1op2 decomposition of the gene-wide shift (same terms and colours as Figures S19/S20)
    D = json.loads((HERE/'opposite-sign/opposite-sign.json').read_text())['Fgfr1op2']
    dec = {r['route']: r for r in D['decomposition']}
    iv = [i for i in D['strict_intervals'] if i['length'] == 970][0]
    ex0, ey0 = M_, r4_top
    ew = (iB[0]+dxB)-INKGAP-ex0; eh = iB[3]+dyB-ey0
    F.panel('e', ex0, ey0, ew, eh)
    F.t.text_top(20.0, ey0+0.2, 'Fgfr1op2', 9.0, font='BolIta')
    ex = Lin(-0.03, 0.06, 16.0, ex0+ew-30.0)
    F.text(ex0+ew, ey0+13.0, 'ΔI', 8.0, ha='right')
    rows = [('CIP', ey0+24.5), ('WOT', ey0+38.5)]
    F.line([ex(0)]*2, [ey0+17.0, ey0+46.0], lw='ref', color=AXC, ls='ref')
    for rt, yy in rows:
        r = dec[rt]; F.study_block(M_, yy-2.8, rt)
        F.rect(ex(r['within']), yy-3.8, ex(0)-ex(r['within']), 7.6, fc=SEQ_STOPS[2])
        F.rect(ex(0), yy-3.8, ex(r['mass'])-ex(0), 7.6, fc=SEQ_STOPS[1])
        F.rect(ex(r['mass']), yy-3.8, ex(r['mass']+r['outside'])-ex(r['mass']), 7.6, fc=OUTSIDE)
        F.marker(ex(r['dDWUI']), yy, 'D', 'small', fc=INK)
        F.text(ex0+ew, yy-3.6, f"{iv[rt]['20']['dI_pp']:.0f} pp".replace('-', '−'), 7.5, ha='right')
    ybot = F.xaxis(ex, ey0+47.0, [-0.02, 0, 0.04], labels=['−0.02', '0', '+0.04'], size=7.5)
    F.key_row(ex0, ybot+1.0, [('block', SEQ_STOPS[2], 'Family position'), ('block', SEQ_STOPS[1], 'Mass')], size=7.5)
    F.key_row(ex0, ybot+10.0, [('block', OUTSIDE, 'Outside family'), ('marker', ('D', INK), 'ΔDWUI')], size=7.5)
    rec['new_b'] = dict(source=[str((NS/'polyapipe/gene_classes.tsv').relative_to(ROOT)), str((NS/'p50/gene_classes.tsv').relative_to(ROOT))],
                        sha=[sha(NS/'polyapipe/gene_classes.tsv'), sha(NS/'p50/gene_classes.tsv')], counts=cnt, counts_alt_not_drawn=cnt5,
                        x='min over the two studies of ΔDWUI signed by the interval change; signed square-root axis', limits=LIM,
                        clipped=int(((g.agree < LIM[0]) | (g.agree > LIM[1])).sum()+((g5.agree < LIM[0]) | (g5.agree > LIM[1])).sum()))
    rec['new_e'] = dict(source='opposite-sign/opposite-sign.json', sha=sha(HERE/'opposite-sign/opposite-sign.json'), decomposition=dec,
                        interval=iv['interval'], dI_20=dict(CIP=iv['CIP']['20'], WOT=iv['WOT']['20']))
    rec['letters'] = dict(a='kept', b='new', c='kept', d=f'old f moved {dyF:+.2f}', e='new', f=f'old b moved ({dxB:+.2f}, {dyB:+.2f})')
    rec['dropped'] = 'old d (partial retention) and e (Eif4ebp2): same content as Supplementary Figure S17c-f'
    finish_fig(dst, F, 'Figure.4', H, rec)


# ===================================================================================================== new Figure 5 (old 4)
def fig5():
    src = FIG/'Figure.4.pdf'; Hs = 481.12; Ls = letters(src)
    REG = (0, 164.0, W, Hs)
    dy = Ls['a']['box'][1] - Ls['c']['box'][1]
    H = round(Hs+dy, 2)
    dst = new_page(H); top = Top(W, H)
    lb = lambda L: (Ls[L]['box'][0]-0.6, Ls[L]['box'][1]-0.6, Ls[L]['box'][2]+0.6, Ls[L]['box'][3]+0.6)
    place(dst, H, crop(src, Hs, REG, [lb('c'), lb('d'), lb('e')]), REG, 0, REG[1]+dy)
    for old, new in (('c', 'a'), ('d', 'b'), ('e', 'c')):
        letter(top, new, Ls[old]['ox'], Ls[old]['base']+dy)
    RECEIPT_ = dict(source=str(src.relative_to(ROOT)), source_sha=sha(src), region=REG, dy=dy, letters='c,d,e -> a,b,c',
                    dropped='old a (study-pair sharing) and b (genome-wide correlation): to Supplementary Figure S13')
    finish(dst, top, 'Figure.5', H); RECEIPT['Figure.5'].update(RECEIPT_)


if __name__ == '__main__':
    which = sys.argv[1:] or ['3', '4', '5']
    for w in which:
        dict(**{'3': fig3, '4': fig4, '5': fig5})[w]()
    old = json.loads((OUT/'receipt.json').read_text()) if (OUT/'receipt.json').exists() else {}
    old.update(RECEIPT); (OUT/'receipt.json').write_text(json.dumps(old, indent=1, default=float, ensure_ascii=False))
    print(json.dumps({k: v.get('page_pt') for k, v in RECEIPT.items()}))
