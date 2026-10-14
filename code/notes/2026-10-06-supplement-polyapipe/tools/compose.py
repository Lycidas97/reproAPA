"""Consolidated supplementary figures (2026-10-06 tightening, 31 -> 21; plan in notes/2026-10-06-suppfig-consolidation).

Each consolidated figure is assembled from the existing per-figure scripts (figures/<work id>.py), not from their PDFs:
a source script runs with supdraw.Fig replaced by a Recorder that stores every drawing primitive per panel in the
script's own coordinates; the kept panels are then placed in bands (rows of blocks, a block = panels of one source moved
together) and replayed into one supdraw.Fig, which runs the usual layout checks, pdfclean, PNG and receipt.
Panel letters are reassigned in reading order (bands, then blocks, then the source's own letter order).

Geometry changes needed to fit a 190 mm page are source-text patches listed in the spec (each must match exactly once);
the standalone scripts and their candidates/ renders are untouched. Values of every kept panel are compared with the
standalone render (candidates/<work id>/...values.json) and the result is written to the receipt.

Usage: python3 tools/compose.py A B ...   (plan letters; no argument = all)
"""
import sys, json, re, shutil, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(HERE))
import supdraw
from supdraw import (PAGE_W, GAP, MARGIN, SIZES, INK, Lin, textwidth, SNUM, outdir as src_outdir, TOKENS_USED)
RealFig = supdraw.Fig
OUT = BASE/'consolidated'
ALT = BASE.parent/'2026-10-06-alt-catalogue-suppfig/final'

# ------------------------------------------------------------------------------------------------------------ spec
# Print order = first citation in the main text (main/results R1-R5, Discussion), checked 2026-10-06:
# R1 S1-S6 (A,B,C), S7 (D), S8 (E); R2 S9 (F), S10 (G), S17 (I), S11 (H); R3 S12 (J), S16 (L), S14 (K), S21 (N);
# R4 S19 (M), S22/S23 (O); R5 S24 (P), S26 (Q), S27/S28 (R), S30 (T), S29 (S); Discussion S31 (U).
PRINT = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'L', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'T', 'S', 'U']
PNUM = {k: i+1 for i, k in enumerate(PRINT)}
SCRIPT = {'D01': ('D_markers.py', ['human']), 'D02': ('D_markers.py', ['mouse']),
          'G1': ('G_markers.py', ['human']), 'G2': ('G_markers.py', ['mouse'])}

MARKERS_C = ('section', '# ================================================================ c  programme scores', "F.values['c'] = Cv", r'''# ================================================================ c  programme scores
x, y, w, h = F.panel('c', X0, Y3, RR-X0, HC)
# Column IDs link the program scores to the retained marker panels.
# Recover space formerly occupied by marker-gene labels so 6.8-pt IDs clear one another.
CX0 = X0+49.0
CW = (RR-1.0-CX0)/NS
bounds = [CX0+CW*(ST.index[ST.unit == u].max()+1) for u in units[:-1]]
for key, us in HEAD:
    idx_ = ST.index[ST.unit.isin(us)]
    name, col = (('Jaber2022', DATASET['GSE171127'][1]) if key == 'GSE171127' else UNIT[key])
    xm = (cx(idx_.min())+cx(idx_.max()))/2; span = CW*len(idx_)
    lab = name if textwidth(name, 7.0)+7.6 < span+12.0 else name.split('–')[-1].replace('Jaber ', '')
    tw = textwidth(lab, 7.0)
    x0u = min(max(xm-(tw+7.6)/2, CX0), RR-1.0-(tw+7.6))
    F.rect(x0u, y+4.2, 5.6, 5.6, fc=col)
    F.text(x0u+7.6, y+3.6, lab, 7.0)
yf = y+13.6
for i, r in ST.iterrows():
    fcol = PAL[r.family]
    F.rect(CX0+CW*i+0.2, yf, CW-0.4, 7.6, fc=fcol)
    code = CODE[r.family]
    F.text(cx(i), yf+0.6, code[-1] if len(code) == 2 else code, 6.8, ha='center', color=on(fcol))
aliases = {'gse171127_d3_getm': 'G3', 'gse171127_d3_oskm': 'O3',
           'gse171127_d6_pooled': 'D6', 'gse171127_d12_pooled': 'D12'}
for i, r in ST.iterrows():
    label = (aliases[r.unit]+' ' if r.unit in aliases else '')+r.short_id
    if r.n_libraries_ge20 != r.n_libraries_ge20 or r.n_libraries_ge20 == 0:
        label += '†'
    F.vlabel(cx(i), y+53.0, label, 6.8)
HDRC = 56.0
PROG = ['ECM', 'Epithelial', 'Pluripotency']
top = max(float(PG.mean_log1p_cp10k.max()), 0.5)
lim = float(np.ceil(top*2)/2)
sub = (HC-HDRC-2*5.0-4.2)/3
Cv = {}
idx = {k: i for i, k in enumerate(ST.state_key)}
assert (PG.family == PG.state_key.map(ST.set_index('state_key').family)).all()
for j, prog in enumerate(PROG):
    sy = Lin(0, lim, y+HDRC+j*(sub+5.0)+sub, y+HDRC+j*(sub+5.0))
    F.text(X0+9.4, (sy.p0+sy.p1)/2-3.4, prog, 7.0)
    F.yaxis(sy, CX0, [0, lim], ['0' if j == 2 else '', f'{lim:g}'])
    q = PG[PG.program == prog]
    for sk, g in q.groupby('state_key'):
        i = idx[sk]; n = len(g)
        offs = np.linspace(-CW*0.28, CW*0.28, n) if n > 1 else [0.0]
        F.points([cx(i)+o for o in offs], [sy.clip(v) for v in g.mean_log1p_cp10k], ms='point', color=PAL[g.family.iloc[0]], alpha='full')
    Cv[prog] = dict(observations=int(len(q)), states=int(q.state_key.nunique()))
unit_rules(yf, y+HC-4.2)
# 2026-10-12: the explanatory line under the panel moved to the caption; the quantity is a y-axis title
F.vtext(X0+3.8, y+HDRC+(HC-HDRC-4.2)/2, 'log1p(CP10k)', 8.0)
F.values['c'] = Cv''')
MARKERS_H = ('HB, HC = 62.0, 96.0', 'HB, HC = 62.0, 122.0')

POWER_E = ('section', '# ================================================================ e  power', "F.values['e'] = Ev", r'''# ================================================================ e  power
# consolidated figure: widened to share a row with the APA rank-effect panel, so one line chart per database
x, y, w, h = F.panel('e', X0+WD+GAP, Y3, 305.2, H3)
pw_ = D.df('D15E', 'functional/functional_power.tsv.gz')
pw_ = pw_[(pw_.model == 'L3_full') & pw_.analysis_id.str.contains('G2_vs')]
py = Lin(0, 1, y+110, y+12)
sg_ = 10.0; xl_ = x+24.0; sw_ = (x+w-3.0-xl_-sg_)/2
Ev = {}
for k, (db, col) in enumerate((('GO', INK), ('KEGG', INK))):          # database named in each subplot title; colour not needed
    px = Lin(1, 8, xl_+k*(sw_+sg_), xl_+k*(sw_+sg_)+sw_)
    F.lineplot(px, py, 'power '+db)
    F.text((px.p0+px.p1)/2, y+3.0, db, 7.5, ha='center')   # 2026-10-12: facet label only; G2 tests / L3 full null in the caption
    F.yaxis(py, px.p0, [0, 0.5, 1], ['0', '.5', '1'], 'Detection probability' if k == 0 else None, title_y=(py.p0+py.p1)/2+6,
            title_x=x+5.4, labels_off=(k == 1))
    for th, ls in (('nominal', None), ('bonferroni_terms', 'dash')):
        q = pw_[(pw_.database == db) & (pw_.threshold == th)].groupby('odds_ratio').rejection_probability.median()
        F.line([px(v) for v in q.index], [py(v) for v in q.values], lw='data', color=col, ls=ls)
        for o_, v in q.items(): F.marker(px(o_), py(v), 'o', 'small', fc=col)
        Ev[f'{db}|{th}'] = {int(k_): float(v) for k_, v in q.items()}
    yb_ = F.xaxis(px, py.p0+1.0, [1, 2, 4, 8])
F.text((xl_+x+w-3.0)/2, yb_+1.6, 'Injected odds ratio', 8.0, ha='center')
F.key_row(None, y+h-7.4, [('line', (INK, None, 'data'), 'Nominal'), ('line', (INK, 'dash', 'data'), 'Bonferroni-corrected over terms')],
          gap=6.0, center=(xl_+x+w-3.0)/2)
F.values['e'] = Ev''')

# bands: list of rows; a row is a list of blocks (source work id, its panel letters in source order)
SPEC = {
    'A': dict(title='human states', bands=[[('G1', 'abcd')], [('D01', 'c')]], patches={'D01': [MARKERS_C, MARKERS_H]}),
    'B': dict(title='mouse states', bands=[[('G2', 'abcd')], [('D02', 'c')], [('G3', 'bcd')]],
              patches={'D02': [MARKERS_C, MARKERS_H],
                       'G3': [('H1, H2 = 136.0, 111.5', 'H1, H2 = 136.0, 88.5'), ('top, bot = y+18.0, y+71.0', 'top, bot = y+18.0, y+57.0'),
                              # new letters f-h: the study block must clear the widest of them, not the old b-d
                              ("tx = x+0.2+textwidth(L, 11, 'Bol')+F.CLEAR_R+0.6", "tx = x+0.2+max(textwidth(c_, 11, 'Bol') for c_ in 'fgh')+F.CLEAR_R+0.6")]}),
    'C': dict(bands=[[('D03', 'abcdef')]]),
    'D': dict(bands=[[('D04', 'ABCD')]]),
    'E': dict(bands=[[('D05', 'abcdef')]]),
    'F': dict(bands=[[('D06', 'ABCDE')]]),
    'G': dict(bands=[[('D07', 'ABCDE')]]),
    'H': dict(bands=[[('D08', 'abcd')]]),
    'I': dict(title='functional enrichment', bands=[[('D14', 'abc')], [('D15', 'ab')], [('D15', 'c')], [('D14', 'e'), ('D15', 'e')]],
              patches={'D14': [('H1, H2 = 87.4, 126.0', 'H1, H2 = 87.4, 140.0'), ('ES = 80.0', 'ES = 94.0')], 'D15': [POWER_E]}),
    'J': dict(bands=[[('D09', 'abcdef')], [('D10', 'abc')]]),
    'K': dict(bands=[[('D11', 'abcde')], [('D12', 'abcd')]]),
    'L': dict(bands=[[('D13', 'abcde')]]),
    'M': dict(bands=[[('D16', 'abcdef')]]),
    'N': dict(bands=[[('D18', 'abcdefg')]]),
    'O': dict(title='cis features and prediction', bands=[[('D19', 'ab')], [('D20', 'ab')], [('D20', 'cd')], [('D19', 'de')]],
              patches={'D20': [('H1, H2, H3 = 96.0, 200.0, 168.0', 'H1, H2, H3 = 96.0, 173.0, 168.0'),
                               ('rp = 6.8; c0 = y+21.0', 'rp = 6.2; c0 = y+20.0'),
                               ('    yy += 4.0\nybot = yy-4.0', '    yy += 3.0\nybot = yy-3.0'),
                               ("F.text(x+LBW+(w-LBW)/2, ybx+1.6, 'Higher = better; primary setting (supported layer, τ = .02)', 7.0, ha='center', color=MUTED)\n",
                                '# consolidated figure: the primary-setting note moved to the legend\n')]}),
    'P': dict(bands=[[('D21', 'abcdef')]]),
    # Q and R differ from the plan (Q = S25 + S28a,b; R = S26a,b + S27): that Q needs 209 mm (the plan's 194 mm omitted the
    # row gaps). Same kept panels, regrouped: Q = site density (S26a,b) + partial retention (S25); R = miRNA (S27) and RBP
    # (S28a,b) expression context, both defining Figure 5f rows.
    'Q': dict(title='site density and partial retention', bands=[[('D23', 'ab')], [('D22', 'abcd')]],
              patches={'D23': [('ax_ = Lin(-80, 130, x+73.0, x+w-46.0)', 'ax_ = Lin(-80, 130, x+79.0, x+w-46.0)'),   # cohort and
                               ("F.text(x+70.0, yy-3.4, cat, 7.0, ha='right', color=MUTED)",                         # category labels
                                "F.text(x+76.0, yy-3.4, cat, 7.0, ha='right', color=MUTED)")]}),                     # touched
    'R': dict(title='miRNA and RBP expression', bands=[[('D24', 'abcdef')], [('D25', 'a')], [('D25', 'b')]],
              patches={'D25': [('H1, H2, H3 = 102.0, 124.0, 150.0', 'H1, H2, H3 = 96.0, 120.0, 150.0'),
                               ('rh = 10.0; r0 = y+28.0', 'rh = 9.0; r0 = y+28.0'),
                               ('by = Lin(0, 9, y+72.0, y+20.5)', 'by = Lin(0, 9, y+68.0, y+20.5)')]}),
    'S': dict(bands=[[('D26', 'abcdefgh')]]),
    'T': dict(bands=[[('D27', 'abc')]]),
    'U': dict(external='A01'),
}
# old printed number of each work id (candidates/ and captions/ use these)
OLD = dict(SNUM); OLD['A01'] = 31


def figname(cid):
    return f'SuppFig.S{PNUM[cid]}'


def cdir(cid):
    return OUT/cid


# ------------------------------------------------------------------------------------------------------------ recorder
class Recorder:
    """Stand-in for supdraw.Fig during a source script: records primitives per panel (script coordinates)."""
    TEXT_H, CLEAR_R, CLEAR_B, DOT_MAX, KEY_GLYPH, ASPECT = (RealFig.TEXT_H, RealFig.CLEAR_R, RealFig.CLEAR_B, RealFig.DOT_MAX,
                                                           RealFig.KEY_GLYPH, RealFig.ASPECT)
    xaxis, yaxis, study_block, study_pair, key_width, key_row, cbar, cbar_key, letter_zone, whisker, arrow_down = (
        RealFig.xaxis, RealFig.yaxis, RealFig.study_block, RealFig.study_pair, RealFig.key_width, RealFig.key_row, RealFig.cbar, RealFig.cbar_key,
        RealFig.letter_zone, RealFig.whisker, RealFig.arrow_down)
    LAST = None

    def __init__(self, H, name, W=PAGE_W):
        self.W, self.H, self.name = W, H, name
        self.panels, self.calls, self.order, self.cur = {}, {}, [], None
        self.values, self.notes, self.saved = {}, [], None
        Recorder.LAST = self

    def panel(self, L, x, y, w, h):
        self.panels[L] = (x, y, w, h); self.cur = L; self.calls.setdefault(L, [])
        if L not in self.order: self.order.append(L)
        return x, y, w, h

    def use(self, L):
        self.cur = L

    def _rec(self, *call):
        if self.cur is None: raise RuntimeError(f'{self.name}: drawing before the first panel: {call[0]}')
        self.calls[self.cur].append(call)

    def text(self, x, top, s, size=7.5, font='Reg', color=INK, ha='left'):
        assert size in SIZES, (s, size); self._rec('text', (x, top, s, size, font, color, ha), {}); return textwidth(s, size, font)

    def vtext(self, x, yc, s, size=8.0, color=INK):
        self._rec('vtext', (x, yc, s, size, color), {})

    def vlabel(self, x, ybot, s, size=7.0, font='Reg', color=INK):
        self._rec('vlabel', (x, ybot, s, size, font, color), {}); return textwidth(s, size, font)

    def dot(self, x, y, frac, fc, dmax):
        self._rec('dot', (x, y, frac, fc, dmax), {})

    def line(self, xs, ys, lw='axis', color=INK, ls=None, alpha='full', **kw):
        self._rec('line', (list(xs), list(ys), lw, color, ls, alpha), kw)

    def rect(self, x, y, w, h, fc='white', ec='none', lw=0, **kw):
        self._rec('rect', (x, y, w, h, fc, ec, lw), kw)

    def marker(self, x, y, m='o', ms='mark', fc=INK, ec=None, mew=None, zorder=3):
        self._rec('marker', (x, y, m, ms, fc, ec, mew, zorder), {})

    def poly(self, pts, fc, ec='none', lw=0, alpha=1.0):
        self._rec('poly', (list(pts), fc, ec, lw, alpha), {})

    def points(self, xs, ys, ms='point', color=INK, alpha='dense', rasterized=False, m='o', mew=0.0, fc=None):
        self._rec('points', (np.array(xs, float), np.array(ys, float), ms, color, alpha, rasterized, m, mew, fc), {})

    def lineplot(self, sx, sy, label='', standalone=True):
        self._rec('lineplot', (Lin(sx.v0, sx.v1, sx.p0, sx.p1), Lin(sy.v0, sy.v1, sy.p0, sy.p1), label, standalone), {})

    def save(self, outdir, data, extra=None):
        self.saved = dict(data=list(data) if isinstance(data, (list, tuple)) else [data], extra=extra)


def replay(M, call, dx, dy):
    name, a, k = call
    if name in ('text', 'vtext', 'vlabel', 'dot', 'rect', 'marker'):
        getattr(M, name)(a[0]+dx, a[1]+dy, *a[2:], **k)
    elif name == 'line':
        M.line([v+dx for v in a[0]], [v+dy for v in a[1]], *a[2:], **k)
    elif name == 'points':
        M.points(a[0]+dx, a[1]+dy, *a[2:], **k)
    elif name == 'poly':
        M.poly([(px+dx, py+dy) for px, py in a[0]], *a[1:], **k)
    elif name == 'lineplot':
        sx, sy = a[0], a[1]
        M.lineplot(Lin(sx.v0, sx.v1, sx.p0+dx, sx.p1+dx), Lin(sy.v0, sy.v1, sy.p0+dy, sy.p1+dy), *a[2:])
    else:
        raise ValueError(name)


# ------------------------------------------------------------------------------------------------------------ sources
def patched_source(fid, patches):
    script, args = SCRIPT.get(fid, (f'{fid}.py', []))
    p = BASE/'figures'/script
    src = p.read_text(); log = []
    for pt in patches:
        if pt[0] == 'section':
            _, start, end, new = pt
            assert src.count(start) == 1 and src.count(end) == 1, (fid, start[:40], end)
            i = src.index(start); j = src.index(end, i)+len(end)
            src = src[:i]+new+src[j:]; log.append(dict(kind='section', start=start, end=end, new_sha256=hashlib.sha256(new.encode()).hexdigest()))
        else:
            old, new = pt
            assert src.count(old) == 1, (fid, old)
            src = src.replace(old, new); log.append(dict(kind='replace', old=old, new=new))
    return p, args, src, log


def run_source(fid, patches):
    p, args, src, log = patched_source(fid, patches)
    supdraw.Fig = Recorder; Recorder.LAST = None
    argv = sys.argv
    try:
        sys.argv = [str(p)]+args
        ns = {'__name__': '__compose__', '__file__': str(p)}
        exec(compile(src, str(p), 'exec'), ns)
    finally:
        sys.argv = argv; supdraw.Fig = RealFig
    R = Recorder.LAST
    assert R is not None and R.saved is not None, f'{fid}: script did not save'
    return R, dict(script=str(p.relative_to(BASE)), args=args, script_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                   patched_sha256=hashlib.sha256(src.encode()).hexdigest(), patches=log)


class Merged:
    def __init__(self, recs):
        self.recs = recs

    def receipt(self):
        out = {}
        for r in self.recs:
            for k, v in r.items():
                if k in ('tables', 'external_sources'):
                    d = out.setdefault(k, {})
                    for t, m in v.items():
                        if t in d and 'panels' in m:
                            d[t]['panels'] = sorted(set(d[t]['panels']) | set(m['panels']))
                        else:
                            d[t] = json.loads(json.dumps(m, default=str))
                elif k in out and out[k] != v:
                    raise RuntimeError(f'receipts disagree on {k}: {out[k]} / {v}')
                else:
                    out[k] = v
        return out


def bbox(rects):
    return (min(r[0] for r in rects), min(r[1] for r in rects), max(r[0]+r[2] for r in rects), max(r[1]+r[3] for r in rects))


def compose(cid):
    sp = SPEC[cid]
    if 'external' in sp:
        return external(cid, sp['external'])
    fids = list(dict.fromkeys(f for band in sp['bands'] for f, _ in band))
    recs, prov = {}, {}
    for f in fids:
        recs[f], prov[f] = run_source(f, sp.get('patches', {}).get(f, []))
    # layout ------------------------------------------------------------------
    place, lmap, warn = [], {f: {} for f in fids}, []
    used = {f: set() for f in fids}
    ycur, nxt = MARGIN, 0
    for band in sp['bands']:
        boxes = []
        for f, Ls in band:
            for L in Ls:
                assert L in recs[f].panels, (cid, f, L, list(recs[f].panels)); assert L not in used[f]; used[f].add(L)
            boxes.append(bbox([recs[f].panels[L] for L in Ls]))
        width = sum(b[2]-b[0] for b in boxes)+GAP*(len(boxes)-1)
        if abs(width-(PAGE_W-2*MARGIN)) > 0.6: warn.append(f'band {[b for b in band]} width {width:.1f} pt')
        hs = [b[3]-b[1] for b in boxes]
        if max(hs)-min(hs) > 0.6: warn.append(f'band {[b for b in band]} heights {[round(v, 1) for v in hs]}')
        xcur = MARGIN
        for (f, Ls), b in zip(band, boxes):
            dx, dy = xcur-b[0], ycur-b[1]
            for L in Ls:
                newL = 'abcdefghijklmnopqrstuvwxyz'[nxt]; nxt += 1
                lmap[f][L] = newL; place.append((f, L, newL, dx, dy))
            xcur += b[2]-b[0]+GAP
        ycur += max(hs)+GAP
    H = ycur-GAP+MARGIN
    # replay ------------------------------------------------------------------
    for d in TOKENS_USED.values(): d.clear()
    M = RealFig(H, figname(cid))
    for f, L, newL, dx, dy in place:
        x, y, w, h = recs[f].panels[L]
        M.panel(newL, x+dx, y+dy, w, h)
        for call in recs[f].calls[L]:
            replay(M, call, dx, dy)
    # values, notes, provenance --------------------------------------------------
    srcvals, mapping, dropped = {}, {}, {}
    for f in fids:
        R = recs[f]; lm = {L.lower(): v for L, v in lmap[f].items()}
        dropped[f] = [L.lower() for L in R.order if L not in lmap[f]]
        old = src_outdir(f)/f'SuppFig.S{OLD[f]}.values.json'
        oldv = json.loads(old.read_text()) if old.exists() else None
        same = []
        for k, v in R.values.items():
            ks = k.lower()
            if all(c in lm or c in dropped[f] for c in ks):
                kept = [c for c in ks if c in lm]
                if not kept: continue
                assert len(kept) == len(ks), (cid, f, k, 'partially dropped values key')
                M.values[''.join(lm[c] for c in ks)] = v
            else:
                M.values[f'{f}:{k}'] = v
            if oldv is not None:
                same.append(json.dumps(oldv.get(k), sort_keys=True, default=str) == json.dumps(json.loads(json.dumps(v, default=str)), sort_keys=True))
        srcvals[f] = dict(standalone_values=str(old.relative_to(BASE)) if old.exists() else None,
                          kept_values_identical=bool(same) and all(same))
        mapping[f] = {L.lower(): newL for L, newL in lmap[f].items()}
        tag = ', '.join(f'{a}→{b}' for a, b in mapping[f].items())
        for n in R.notes:
            M.notes.append(f'[{f} (old S{OLD[f]}) {tag}] {n}')
    rec = Merged([d.receipt() for f in fids for d in recs[f].saved['data']])
    extra = dict(consolidated=dict(plan_letter=cid, title=sp.get('title'), print_order=PNUM[cid], bands=[[[f, Ls] for f, Ls in band] for band in sp['bands']],
                                   letter_map={f'{f} (old S{OLD[f]})': mapping[f] for f in fids},
                                   dropped_panels={f'{f} (old S{OLD[f]})': dropped[f] for f in fids if dropped[f]},
                                   sources=prov, values_check=srcvals, layout_warnings=warn))
    chk = M.save(cdir(cid), rec, extra)
    print(f'{cid} -> {figname(cid)}  {H/72*25.4:.1f} mm  map {mapping}  values identical {[v["kept_values_identical"] for v in srcvals.values()]}')
    for w_ in warn: print('  LAYOUT', w_)
    return chk


def external(cid, fid):
    """S31 (A01) is drawn in notes/2026-10-06-alt-catalogue-suppfig from its own pinned tables; copy it under the new number."""
    d = cdir(cid); d.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png', 'values.json'):
        shutil.copyfile(ALT/f'SuppFig.S31.{ext}', d/f'{figname(cid)}.{ext}')
    r = json.loads((ALT/'SuppFig.S31.receipt.json').read_text())
    r['consolidated'] = dict(plan_letter=cid, print_order=PNUM[cid], copied_from='notes/2026-10-06-alt-catalogue-suppfig/final/SuppFig.S31.pdf',
                             source_pdf_sha256=hashlib.sha256((ALT/'SuppFig.S31.pdf').read_bytes()).hexdigest(),
                             letter_map={'A01 (old S31)': {L.lower(): L.lower() for L in r['panels']}})
    r['figure'] = figname(cid)
    (d/f'{figname(cid)}.receipt.json').write_text(json.dumps(r, indent=1, ensure_ascii=False))
    print(f'{cid} -> {figname(cid)}  copied from A01 (old S31)')


if __name__ == '__main__':
    for cid in (sys.argv[1:] or PRINT):
        compose(cid)
