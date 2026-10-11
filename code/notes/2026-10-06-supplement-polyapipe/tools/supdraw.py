"""Supplementary-figure toolkit (polyApipe v2-primary redraw, 2026-10-06).

House style is taken from the main figures (v5 builders): Nimbus Sans L, type sizes 6.8/7.0/7.5/8.0 and 11 (panel
letters), ink #263D4A, axis grey #87969E at 0.5 pt, 4.0 mm (11.3 pt) gaps between panel rectangles, ~4.5 pt outer
margin, panel letter in a fixed top-left slot inside the rectangle, study colour blocks instead of dataset names.
All coordinates are top-origin PDF points.

Data come from ONE result bundle (the supplement package).  A panel reads a table either through its registered
component (full ResultBundle checks) or, when the bundle registers no component for that panel, through the plot
contract: the table must be bound to that panel and its file must match the table-manifest SHA.  Every read is
recorded in the receipt.
"""
import sys, json, math, subprocess
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'Fig').is_dir() and (p/'main.tex').exists())
sys.path.insert(0, str(ROOT/'notes/2026-10-06-main-figures-polyapipe/builders/tools-local'))
sys.path.insert(0, str(ROOT/'analysis_code/pas_rerun'))
from v5draw import Top, textwidth, INK, MUTED, AXC, NULLBAR, NULLW, MISSING, BAND, num, qfmt, sha  # noqa
# Diverging scale = orange/blue DIV #0072B2 / #F7F7F7 / #D55E00 (supplement house style). Defined here on purpose: the
# main-figure v5draw switched to Crameri vik on 2026-10-06, and a vik trial for the supplement was rejected by the user the
# same day, so the supplement no longer follows v5draw.diverge.
from matplotlib.colors import LinearSegmentedColormap as _LSC
DIV = ['#0072B2', '#F7F7F7', '#D55E00']
DIV_CMAP = _LSC.from_list('d', DIV)
DIV_NAME = 'orange/blue #0072B2/#F7F7F7/#D55E00'


def diverge(v, lim):
    """Signed value -> diverging colour; clipped at ±lim (blue = negative / proximal, orange = positive / distal)."""
    return DIV_CMAP((np.clip(v, -lim, lim)+lim)/(2*lim))
from result_records import ResultBundle, contained   # noqa
from pipeline import digest                           # noqa

BUNDLE_R1 = ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/supplement'   # sealed statistics bundle (contract r1)
BUNDLE = Path(__file__).resolve().parent.parent/'contract-r2'/'bundle'                # same tables, contract r2
PAGE_W = 453.54
GAP = 11.3
MARGIN = 4.5
SIZES = (6.8, 7.0, 7.5, 8.0, 11.0)
TICK = 2.0

# ---------------------------------------------------------------- design tokens (the only values a figure may use)
# Line widths by role (pt). Dominant values of the main figures, collapsed to one value per role.
LW = dict(ref=0.4,      # reference / guide lines: dotted zero or expectation lines, separators
          axis=0.5,     # axes, ticks, cell outlines, colour-bar frames
          thin=0.8,     # secondary data: non-highlighted series, connectors, hollow-marker edges
          data=1.0,     # primary data lines and curves, observed-value lines
          bold=1.6,     # emphasis: highlighted series, interval bars (IQR / CI)
          range=2.2)    # wide range bars: null or setting ranges
MS = dict(point=1.2, small=3.0, mark=3.6)          # marker sizes: point clouds, secondary marks, primary marks
MEW = dict(none=0.0, fill=0.5, hollow=0.8)          # marker edge widths: filled markers, hollow markers
ALPHA = dict(full=1.0, secondary=0.6, dense=0.35)   # opaque; secondary layers (draws, per-library lines); dense clouds
DASH = dict(ref=(0, (1, 1.5)), dash=(0, (3.2, 1.6)), dot=(0, (0.9, 1.3)),
            dashdot=(0, (4, 1.5, 1, 1.5)))   # dashdot = state M in PAS-usage curves (main Figures 3d, 4e, 5e)
TOKENS_USED = {'lw': {}, 'ms': {}, 'mew': {}, 'alpha': {}, 'size': {}}


def _tok(table, name, v):
    """Resolve a token name (or an exact token value) and count its use; anything else is an error."""
    tab = {'lw': LW, 'ms': MS, 'mew': MEW, 'alpha': ALPHA}[table]
    if isinstance(v, str):
        if v not in tab: raise ValueError(f'unknown {table} token {v!r}')
        key, val = v, tab[v]
    else:
        hit = [k for k, x in tab.items() if abs(x-v) < 1e-9]
        if not hit: raise ValueError(f'{table}={v} is not a design token {tab}')
        key, val = hit[0], v
    TOKENS_USED[table][key] = TOKENS_USED[table].get(key, 0)+1
    return val

# Figure numbering. Scripts and candidate folders use working IDs (D01–D27 = statistics-contract figures, G1–G3 = the
# 2026-09-14 state-grouping material added on 2026-10-06); the printed S number comes only from this table, so a
# later reordering changes one line. Current order: S1–S3 = D01–D03, S4–S6 = G1–G3, S7–S30 = D04–D27.
ORDER = ['D01', 'D02', 'D03', 'G1', 'G2', 'G3'] + [f'D{i:02d}' for i in range(4, 28)]
SNUM = {fid: i+1 for i, fid in enumerate(ORDER)}


def figname(fid):
    return f'SuppFig.S{SNUM[fid]}'


def outdir(fid):
    return Path(__file__).resolve().parent.parent/'candidates'/fid


# state colours (main figures)
SC = dict(S='#E69F00', R='#0072B2', E='#009E73', P='#CC79A7', M='#8C6D31')
PROX, DIST = '#0072B2', '#D55E00'

# Sequential heatmaps use the Figure 2e ramp (blue family of the diverging DIV #0072B2 / #F7F7F7 / #D55E00, calmer
# than DIV's blue arm for large filled areas), so main and supplement share one sequential and one diverging scale.
# SEQ_STOPS can be overridden (SUP_SEQ env) for side-by-side comparisons; 'slate' is the superseded ramp.
import os as _os
SEQ_SETS = {
    'divblue': ['#F7F7F7', '#0072B2', '#004A75'],            # blue arm of DIV, extended to a darker end
    'main':    ['#EEF3F7', '#9CC3DD', '#2D6A96', '#163A57'],  # Figure 2e sequential (builders/figure2/build.py SEQ)
    'slate':   ['#F3F5F6', '#9DB3BE', '#2E4A58'],             # superseded supplement ramp
}
SEQ_NAME = _os.environ.get('SUP_SEQ', 'main')
SEQ_STOPS = SEQ_SETS[SEQ_NAME]


# ------------------------------------------------------------------------------- named colour tokens (2026-10-06)
# Every colour a figure draws must be one of these names or lie on the sequential / diverging ramp; Fig.check() lists any
# other colour under checks['colour_unnamed'] and tools/audit.py counts it as a hard failure.
ORD3 = (SEQ_STOPS[3], SEQ_STOPS[2], SEQ_STOPS[1])   # ordinal categories, dark -> light (Figure 2e ramp stops)
TRACK = '#E6EAED'     # empty part of a part-of-whole bar (main Figures 3b, 3c)
SHADE = '#DCEAF3'     # interval / window shading (main Figure 5e)
OUTSIDE = '#C9B68F'   # outside the terminal family (main Figure 5f)
DISAGREE = '#E3B04B'  # the two studies disagree in direction (Supplementary Figure S11 relation classes); used nowhere else
ALT_ONE = '#DDE3E7'   # set found with one PAS-calling workflow only (Supplementary Figure S21 d-f, user-approved v4 design)
ALT_OTHER = '#8297A4' # matched label, other class (Supplementary Figure S21 d)
SPECIES = {'mouse': '#146B4A', 'human': '#6F2C91'}   # species markers when a panel pools routes (main Figure 2f)
# local state families of the frozen palette (frozen/palette.json; same colours as main Figure 1a)
FAMILY_LOCAL = {'K': '#A65628', 'hN': '#BCBD22', 'hX': '#A65628', 'hU': '#4B5563', 'mB': '#7B61A8', 'mN': '#4C5A9E',
                'mI': '#17BECF', 'mY': '#E15759', 'mU': '#9AA0A6'}


def seq(v, lo=0.0, hi=1.0):
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list('seq_'+SEQ_NAME, SEQ_STOPS)
    return cm(min(max((v-lo)/(hi-lo), 0.0), 1.0))


def on(fc):
    """Text colour for a label inside a filled cell: white on dark fills, ink otherwise (WCAG relative luminance)."""
    from matplotlib.colors import to_rgb
    lin = [c/12.92 if c <= 0.04045 else ((c+0.055)/1.055)**2.4 for c in to_rgb(fc)]
    return 'white' if 0.2126*lin[0]+0.7152*lin[1]+0.0722*lin[2] < 0.30 else INK
# study / route identity: name, colour, marker.  Figure 4/5 colours where they exist, Figure 2 otherwise.
STUDY = {
    'GSE118258': ('Xing2020', '#753399', 'o'), 'GSE150311': ('Liu2020', '#648290', 's'),
    'GSE221739': ('Panariello2023', '#30363D', '^'),
    '178325_0330': ('Guan2022–0330', '#73501A', 'D'), '178325_0605': ('Guan2022–0605', '#A66B00', 'P'),
    '178325_0618': ('Guan2022–0618', '#A82C60', 'X'), '178325_0809': ('Guan2022–0809', '#C47A2C', '>'),
    'CIP': ('Guo2019–CIP', '#CF4C00', 's'), 'WOT': ('Schiebinger2019–2i', '#207450', 'o'),
    'OSK': ('Guo2019–OSK', '#255FBC', '^'), '171127_GETM': ('Jaber2022–GETM D3', '#766A8A', 'D'),
    '171127_OSKM': ('Jaber2022–OSKM D3', '#B3A4C4', 'v')}
FIT_ROUTES = ['WOT', 'CIP', 'GSE118258']
# dataset (GEO series) identity for per-dataset QC panels; colour of the study's main route
DATASET = {'GSE118258': ('Xing2020', '#753399'), 'GSE150311': ('Liu2020', '#648290'), 'GSE178325': ('Guan2022', '#A82C60'),
           'GSE221739': ('Panariello2023', '#30363D'), 'GSE103221': ('Guo2019', '#CF4C00'), 'GSE115943': ('Schiebinger2019', '#207450'),
           'GSE171127': ('Jaber2022', '#766A8A')}
DATASET_ORDER = ['GSE118258', 'GSE150311', 'GSE178325', 'GSE221739', 'GSE103221', 'GSE115943', 'GSE171127']
DATASET_DIR = {'GSE115943': 'WOT_GSE115943'}          # the three S–R–E–P fitting routes, display order of Figure 2


# --------------------------------------------------------------------------------------------- data
# Tables a redesigned panel needs but the bundle's plot contract does not bind to that panel.  Each entry is read only
# after the table-manifest SHA check and is reported in the receipt; the list is the to-do for the contract update.
# Contract r2 (contract-r2/, 2026-10-06) binds every table the house-style supplement reads; nothing is read outside the
# interface any more. Add an entry here only for a new out-of-interface read, and register it in the next revision.
EXCEPTIONS = {}


class Data:
    def __init__(self, bundle=BUNDLE):
        self.rb = ResultBundle(bundle); self.root = self.rb.root; self.reads = {}
        self.revision = json.loads((self.root/'plot-contract.json').read_text()).get('contract_revision', 'r1')
        ExtData._D = self

    def path(self, panel, table):
        rel = table if table.startswith('tables/') else 'tables/'+table
        key = f'{panel}:{rel[7:]}'
        if key in self.rb.components:
            p = self.rb.table(panel, key); mode = 'component'
        else:
            exc = EXCEPTIONS.get((panel, rel))
            if rel not in self.rb.allowed.get(panel, set()) and exc is None:
                raise SystemExit(f'{panel}: {rel} is not bound to this panel in plot-contract.json')
            meta = self.rb.manifest.get(rel)
            p = contained(self.root, rel)
            if meta is None or digest(p) != meta['sha256']:
                raise SystemExit(f'{rel}: file does not match the table manifest')
            mode = 'plot_contract' if exc is None else 'interface_exception: '+exc
        self.reads.setdefault(rel, dict(sha256=self.rb.manifest[rel]['sha256'], panels=[], mode=mode))
        if panel not in self.reads[rel]['panels']: self.reads[rel]['panels'].append(panel)
        return p

    def df(self, panel, table, **kw):
        import pandas as pd
        return pd.read_csv(self.path(panel, table), sep='\t', **kw)

    def receipt(self):
        return dict(source=self.rb.identity, result_records_sha256=digest(self.root/'result-records.json'), contract_revision=self.revision,
                    bundle=str(self.root).split('/notes/')[-1],
                    tables=self.reads, statistics_recomputed_by_renderer=False, production_replaced=False)


class ExtData:
    """Formerly external, now contract-r2 tables (tools/external-sources*.json). Each key is read through the bundle at its
    `bundle_table`, bound to `panel` (or the first listed panel), after checking the pinned origin SHA; reads are listed in
    the receipt under external_sources with their origin path."""
    _D = None

    def __init__(self, manifest=Path(__file__).resolve().parent/'external-sources.json', panel=None, data=None):
        self.m = json.loads(Path(manifest).read_text()); self.reads = {}; self.panel = panel
        self.D = data or ExtData._D or Data(); ExtData._D = self.D

    def path(self, key):
        f = self.m['files'][key]
        panel = self.panel if self.panel in f['panels'] else f['panels'][0]
        p = self.D.path(panel, f['bundle_table'])
        got = digest(p)
        if got != f['sha256']:
            raise RuntimeError(f'{key}: {p} sha256 {got} differs from the pinned {f["sha256"]}')
        self.reads[key] = dict(bundle_table=f['bundle_table'], panel=panel, origin=f['path'], sha256=got)
        return p

    def df(self, key, **kw):
        import pandas as pd
        return pd.read_csv(self.path(key), sep='\t', **kw)

    def receipt(self):
        return dict(external_sources=self.reads, external_manifest_pinned=self.m['pinned'], external_note=self.m['note'])


# --------------------------------------------------------------------------------------------- scales
class Lin:
    def __init__(self, v0, v1, p0, p1):
        self.v0, self.v1, self.p0, self.p1 = v0, v1, p0, p1
    def __call__(self, v):
        return self.p0+(v-self.v0)/(self.v1-self.v0)*(self.p1-self.p0)
    def clip(self, v):
        lo, hi = min(self.v0, self.v1), max(self.v0, self.v1)
        return self(min(max(v, lo), hi))


# --------------------------------------------------------------------------------------------- colour registry
def _named():
    from matplotlib.colors import to_hex
    d = {}
    def add(n, c): d.setdefault(to_hex(c).upper(), n)
    for n, c in dict(INK=INK, MUTED=MUTED, AXC=AXC, NULLBAR=NULLBAR, NULLW=NULLW, MISSING=MISSING, BAND=BAND, PROX=PROX, DIST=DIST,
                     TRACK=TRACK, SHADE=SHADE, OUTSIDE=OUTSIDE, DISAGREE=DISAGREE, ALT_ONE=ALT_ONE, ALT_OTHER=ALT_OTHER, white='white').items(): add(n, c)
    for k, c in SC.items(): add('state '+k, c)
    for k, c in FAMILY_LOCAL.items(): add('family '+k, c)
    for k, c in SPECIES.items(): add('species '+k, c)
    for k, (n, c, m) in STUDY.items(): add('study '+n, c)
    for k, (n, c) in DATASET.items(): add('dataset '+n, c)
    for i, c in enumerate(SEQ_STOPS): add(f'seq stop {i}', c)
    for i, c in enumerate(DIV): add(f'div stop {i}', c)
    return d


def colour_unnamed(colours):
    """Colours that are neither a named token nor a point of the sequential / diverging ramp."""
    from matplotlib.colors import to_rgba, to_hex
    named = _named()
    ramps = [[to_rgba(seq(t)) for t in np.linspace(0, 1, 401)], [to_rgba(DIV_CMAP(t)) for t in np.linspace(0, 1, 401)]]
    bad = []
    for c, uses in colours.items():
        rgba = to_rgba(c)
        if to_hex(rgba).upper() in named: continue
        if any(min(sum((a-b)**2 for a, b in zip(rgba[:3], r[:3])) for r in ramp) < 2e-5 for ramp in ramps): continue
        bad.append([to_hex(rgba).upper(), sorted(f'{p}:{k}' for p, k in uses)[:6]])
    return bad


# --------------------------------------------------------------------------------------------- figure
class Fig:
    """Top-origin canvas with panel rectangles, a bbox registry and layout checks."""
    TEXT_H = 0.93            # text box height / size (Nimbus Sans L ascender+descender, measured)
    CLEAR_R, CLEAR_B = 8.0, 5.0   # letter clear zone: pt to the right of / below the letter box

    def __init__(self, H, name, W=PAGE_W):
        self.W, self.H, self.name = W, H, name
        self.t = Top(W, H); self.panels = {}; self.items = []; self.cur = None; self.values = {}; self.notes = []
        self.lineplots = []      # (panel, label, x0, y0, x1, y1) data areas of line charts, aspect checked
        self.colours = {}        # colour -> list of (panel, primitive), checked against the named tokens

    # panels -----------------------------------------------------------------
    def panel(self, L, x, y, w, h):
        self.panels[L] = (x, y, w, h); self.cur = L
        self.text(x+0.2, y+0.4, L.lower(), 11, font='Bol')          # panel letters are lowercase (house style)
        return x, y, w, h

    def letter_zone(self, L):
        """Clear zone around the panel letter: the letter box grown by CLEAR_R to the right and CLEAR_B below.
        Nothing else may enter it (house rule: no dense ink next to the letter)."""
        x, y, w, h = self.panels[L]
        lw = textwidth(L.lower(), 11, 'Bol')
        return (x, y, x+0.2+lw+self.CLEAR_R, y+0.4+self.TEXT_H*11+self.CLEAR_B)

    def use(self, L):
        self.cur = L

    # registry ---------------------------------------------------------------
    def _reg(self, kind, x0, y0, x1, y1, text=None, size=None):
        self.items.append(dict(panel=self.cur, kind=kind, box=(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)),
                               text=text, size=size))

    def _c(self, kind, *cols):
        for c in cols:
            if c is None or (isinstance(c, str) and c == 'none'): continue
            key = c if isinstance(c, str) else tuple(round(float(v), 4) for v in np.asarray(c).ravel()[:4])
            self.colours.setdefault(key, set()).add((self.cur, kind))

    # primitives ---------------------------------------------------------------
    def text(self, x, top, s, size=7.5, font='Reg', color=INK, ha='left'):
        assert size in SIZES, (s, size)
        TOKENS_USED['size'][size] = TOKENS_USED['size'].get(size, 0)+1
        self._c('text', color)
        w = textwidth(s, size, font)
        x0 = x if ha == 'left' else (x-w/2 if ha == 'center' else x-w)
        self._reg('text', x0, top, x0+w, top+self.TEXT_H*size, s, size)
        self.t.text_top(x, top, s, size, font=font, color=color, ha=ha)
        return w

    def vtext(self, x, yc, s, size=8.0, color=INK):
        assert size in SIZES
        self._c('text', color)
        w = textwidth(s, size)
        self._reg('text', x-self.TEXT_H*size/2, yc-w/2, x+self.TEXT_H*size/2, yc+w/2, s, size)
        self.t.vtext(x, yc, s, size, color=color)

    def vlabel(self, x, ybot, s, size=7.0, font='Reg', color=INK):
        """Rotated (bottom-to-top) label centred on x whose lower end sits at ybot (column headers)."""
        assert size in SIZES
        TOKENS_USED['size'][size] = TOKENS_USED['size'].get(size, 0)+1
        w = textwidth(s, size, font)
        self._c('text', color)
        self._reg('text', x-self.TEXT_H*size/2, ybot-w, x+self.TEXT_H*size/2, ybot, s, size)
        self.t.s.text(x, self.H-ybot, s, size, font, color, 'center', 'bottom', rotation=90)
        return w

    DOT_MAX = 7.2   # largest dot diameter (pt) of an area-encoded dot plot

    def dot(self, x, y, frac, fc, dmax):
        """Dot-plot cell: area proportional to frac (diameter dmax*sqrt(frac)); data-encoded, not a size token."""
        assert dmax <= self.DOT_MAX+1e-9
        self._c('dot', fc)
        d = dmax*math.sqrt(max(frac, 0.0))
        if d <= 0: return
        self._reg('marker', x-d/2, y-d/2, x+d/2, y+d/2)
        self.t.ax.plot([x], [self.H-y], ls='none', marker='o', ms=d, mfc=fc, mec='none', mew=0, zorder=3)

    def line(self, xs, ys, lw='axis', color=INK, ls=None, alpha='full', **kw):
        xs, ys = list(xs), list(ys)
        lw = _tok('lw', 'lw', lw); kw['alpha'] = _tok('alpha', 'alpha', alpha); self._c('line', color)
        if ls is not None: kw['ls'] = DASH.get(ls, ls) if isinstance(ls, str) and ls in DASH else ls
        self._reg('line', min(xs), min(ys), max(xs), max(ys))
        self.t.line(xs, ys, lw=lw, color=color, **kw)

    def rect(self, x, y, w, h, fc='white', ec='none', lw=0, **kw):
        if lw: lw = _tok('lw', 'lw', lw)
        self._c('rect', fc, ec if lw else None)
        self._reg('rect', x, y, x+w, y+h)
        self.t.rect(x, y, w, h, fc=fc, ec=ec, lw=lw, **kw)

    def marker(self, x, y, m='o', ms='mark', fc=INK, ec=None, mew=None, zorder=3):
        ms = _tok('ms', 'ms', ms); self._c('marker', fc, ec)
        if mew is None: mew = 'hollow' if fc in ('white', 'none') else 'fill'
        if m == '|': mew = 'bold'; mew_v = LW['bold']
        else: mew_v = _tok('mew', 'mew', mew)
        mew = mew_v
        r = ms/2+mew/2
        self._reg('marker', x-r, y-r, x+r, y+r)
        self.t.marker(x, y, m, ms, fc=fc, ec=fc if ec is None else ec, mew=mew, zorder=zorder)

    def poly(self, pts, fc, ec='none', lw=0, alpha=1.0):
        from matplotlib.patches import Polygon
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        self._c('poly', fc, ec if lw else None)
        self._reg('poly', min(xs), min(ys), max(xs), max(ys))
        self.t.ax.add_patch(Polygon([(x, self.H-y) for x, y in pts], closed=True, fc=fc, ec=ec, lw=lw, alpha=alpha))

    def points(self, xs, ys, ms='point', color=INK, alpha='dense', rasterized=False, m='o', mew=0.0, fc=None):
        xs, ys = np.asarray(xs, float), np.asarray(ys, float)
        ms = _tok('ms', 'ms', ms); alpha = _tok('alpha', 'alpha', alpha); self._c('points', color, fc)
        if len(xs):
            r = ms/2
            self._reg('points', xs.min()-r, ys.min()-r, xs.max()+r, ys.max()+r)
            self.t.ax.plot(xs, self.H-ys, ls='none', marker=m, ms=ms, mfc=color if fc is None else fc, mec=color,
                           mew=mew, alpha=alpha, rasterized=rasterized, zorder=2)

    # composite helpers ------------------------------------------------------
    def xaxis(self, sc, y, ticks, labels=None, title=None, size=7.5, tsize=8.0, line=True, title_gap=1.6):
        if line:
            self.line([sc.p0, sc.p1], [y, y], lw='axis', color=AXC)
        for i, v in enumerate(ticks):
            self.line([sc(v)]*2, [y, y+TICK], lw='axis', color=AXC)
            lab = labels[i] if labels else num(v, 0 if float(v).is_integer() else 2)
            self.text(sc(v), y+TICK+1.0, lab, size, ha='center')
        if title:
            self.text((sc.p0+sc.p1)/2, y+TICK+1.0+self.TEXT_H*size+title_gap, title, tsize, ha='center')
        return y+TICK+1.0+self.TEXT_H*size+(title_gap+self.TEXT_H*tsize if title else 0)

    def yaxis(self, sc, x, ticks, labels=None, title=None, size=7.5, tsize=8.0, line=True, title_x=None, title_y=None,
              labels_off=False):
        if line:
            self.line([x, x], [sc.p0, sc.p1], lw='axis', color=AXC)
        wmax = 0
        for i, v in enumerate(ticks):
            self.line([x-TICK, x], [sc(v)]*2, lw='axis', color=AXC)
            if labels_off: continue
            lab = labels[i] if labels else num(v, 0 if float(v).is_integer() else 2)
            wmax = max(wmax, self.text(x-TICK-1.2, sc(v)-self.TEXT_H*size/2+0.3, lab, size, ha='right'))
        if title:
            tx = title_x if title_x is not None else x-TICK-1.2-wmax-2.0-self.TEXT_H*tsize/2
            self.vtext(tx, (sc.p0+sc.p1)/2 if title_y is None else title_y, title, tsize)
        return wmax

    def whisker(self, x0, x1, y, color=INK, lw='thin', cap=1.3):
        """Interval with end caps (main-figure CI style, e.g. Figures 4d, 5a, 5d)."""
        Top.whisker(self, x0, x1, y, color=color, lw=lw, cap=cap)

    def arrow_down(self, x, y0, y1, color=MUTED):
        """Vertical flow arrow (main Figure 5c)."""
        self.line([x, x], [y0, y1-2.6], lw='thin', color=color)
        self.poly([(x-2.2, y1-3.2), (x+2.2, y1-3.2), (x, y1)], fc=color)

    def study_block(self, x, y, key, s=5.6):
        self.rect(x, y, s, s, fc=STUDY[key][1])

    def study_pair(self, xm, yy, hr, mr):
        self.rect(xm, yy-2.8, 5.6, 5.6, fc=STUDY[hr][1]); self.text(xm+7.2, yy-3.4, '×', 7.0, color=MUTED)
        self.rect(xm+12.6, yy-2.8, 5.6, 5.6, fc=STUDY[mr][1])

    KEY_GLYPH = dict(block=7.6, marker=7.6, line=12.5)

    def key_width(self, items, size=7.5, gap=8.0):
        return sum(self.KEY_GLYPH[k]+textwidth(lab, size)+gap for k, _, lab in items)-gap

    def key_row(self, x, y, items, size=7.5, gap=8.0, center=None):
        """items: (kind, spec, label); kind 'block' (colour), 'marker' ((m, colour)), 'line' ((colour, ls, lw)).
        center: x of the subplot centre the key refers to; the whole row is centred on it (house style)."""
        if center is not None:
            x = center-self.key_width(items, size, gap)/2
        for kind, spec, lab in items:
            if kind == 'block':
                self.rect(x, y+1.0, 5.6, 5.6, fc=spec, ec=AXC if spec == 'white' else 'none', lw='axis' if spec == 'white' else 0); x += 7.6
            elif kind == 'marker':                      # (m, fc) or (m, fc, ec)
                self.marker(x+2.8, y+3.8, spec[0], 'mark', fc=spec[1], ec=spec[2] if len(spec) > 2 else None); x += 7.6
            elif kind == 'line':                        # (colour, dash or None, lw token)
                self.line([x, x+10.0], [y+3.8]*2, lw=spec[2], color=spec[0], ls=spec[1]); x += 12.5
            x += self.text(x, y, lab, size)+gap
        return x-gap

    def cbar(self, x0, y0, wbar, lo_lab, hi_lab, fn, n=40, hbar=4.2):
        """Horizontal colour bar with end labels (7.0 pt) left and right of the bar."""
        for i in range(n):
            self.rect(x0+wbar*i/n, y0, wbar/n+0.05, hbar, fc=fn(i/(n-1)))
        self.text(x0-1.8, y0-1.6, lo_lab, 7.0, ha='right'); self.text(x0+wbar+1.8, y0-1.6, hi_lab, 7.0)

    def cbar_key(self, center, y, title, lo_lab, hi_lab, fn, wbar=60.0, size=7.5, gap=6.0):
        """Legend line 'title  lo [bar] hi', centred on the subplot it explains."""
        tw = textwidth(title, size); lw_ = textwidth(lo_lab, 7.0); hw = textwidth(hi_lab, 7.0)
        x0 = center-(tw+gap+lw_+1.8+wbar+1.8+hw)/2
        self.text(x0, y, title, size)
        self.cbar(x0+tw+gap+lw_+1.8, y+1.4, wbar, lo_lab, hi_lab, fn)

    ASPECT = (1.0, 2.0)          # house rule: a panel's own line-chart subplot has data width/height within [1, 2]

    def lineplot(self, sx, sy, label='', standalone=True):
        """Register a line chart's data area (x scale, y scale). A chart that is the panel's own subplot (the panel is
        made of line charts only) must have ASPECT[0] <= w/h <= ASPECT[1]; a line chart that is one component of a
        composite panel (e.g. beside a heatmap or forest plot) is recorded but exempt (standalone=False)."""
        x0, x1 = sorted((sx.p0, sx.p1)); y0, y1 = sorted((sy.p0, sy.p1))
        self.lineplots.append((self.cur, label, x0, y0, x1, y1, standalone))

    # checks -----------------------------------------------------------------
    def check(self, tol=0.35):
        out = dict(outside=[], text_overlaps=[], sizes=sorted({i['size'] for i in self.items if i['size']}))
        for it in self.items:
            P = self.panels.get(it['panel'])
            if P is None: out['outside'].append(dict(item=it, why='no panel')); continue
            x, y, w, h = P; b = it['box']
            if b[0] < x-tol or b[1] < y-tol or b[2] > x+w+tol or b[3] > y+h+tol:
                out['outside'].append(dict(panel=it['panel'], kind=it['kind'], text=it['text'], box=[round(v, 1) for v in b],
                                           rect=[round(x, 1), round(y, 1), round(x+w, 1), round(y+h, 1)]))
        tx = [i for i in self.items if i['kind'] == 'text']
        for i in range(len(tx)):
            for j in range(i+1, len(tx)):
                a, b = tx[i]['box'], tx[j]['box']
                ox = min(a[2], b[2])-max(a[0], b[0]); oy = min(a[3], b[3])-max(a[1], b[1])
                if ox > 0.4 and oy > 0.6:
                    out['text_overlaps'].append([tx[i]['text'], tx[j]['text'], round(ox, 1), round(oy, 1)])
        out.update(self.density())
        out['line_aspect'] = [[L, lab, round((x1-x0)/(y1-y0), 2), 'standalone' if st else 'component']
                              for L, lab, x0, y0, x1, y1, st in self.lineplots]
        out['line_aspect_bad'] = [v for v in out['line_aspect'] if v[3] == 'standalone' and
                                  not (self.ASPECT[0]-0.005 <= v[2] <= self.ASPECT[1]+0.005)]
        out['letter_zone'] = []
        for L in self.panels:
            z = self.letter_zone(L)
            for it in self.items:
                if it['panel'] != L or (it['text'] == L.lower() and it['size'] == 11.0): continue
                b = it['box']
                if b[0] < z[2] and b[2] > z[0] and b[1] < z[3] and b[3] > z[1]:
                    out['letter_zone'].append([L, it['kind'], it['text'], [round(v, 1) for v in b]])
        out['colour_unnamed'] = colour_unnamed(self.colours)
        # panel rectangles: gaps and margins
        R = self.panels
        out['margins'] = dict(left=round(min(v[0] for v in R.values()), 2), top=round(min(v[1] for v in R.values()), 2),
                              right=round(self.W-max(v[0]+v[2] for v in R.values()), 2),
                              bottom=round(self.H-max(v[1]+v[3] for v in R.values()), 2))
        return out

    def density(self, letter=(12.0, 11.0)):
        """Ink geometry per panel: slack between the ink box and the rectangle, the widest empty strip that crosses
        the whole panel (vertical: full height; horizontal: full width), and ink gaps between adjacent panels.
        House rule: slack <= 1.5 pt, internal strips < 8 pt (clearly below the 11.3 pt panel gap)."""
        res = dict(slack={}, internal_vertical_strip={}, internal_horizontal_strip={}, ink_gaps={})
        ink = {}
        for L, (x, y, w, h) in self.panels.items():
            bx = [i['box'] for i in self.items if i['panel'] == L]
            x0 = min(b[0] for b in bx); y0 = min(b[1] for b in bx); x1 = max(b[2] for b in bx); y1 = max(b[3] for b in bx)
            ink[L] = (x0, y0, x1, y1)
            res['slack'][L] = dict(left=round(x0-x, 1), top=round(y0-y, 1), right=round(x+w-x1, 1), bottom=round(y+h-y1, 1))
            def widest(intervals, lo, hi):
                iv = sorted((max(a, lo), min(b, hi)) for a, b in intervals if b > lo and a < hi)
                best, at, cur = 0.0, None, lo
                for a, b in iv:
                    if a-cur > best: best, at = a-cur, (round(cur-(x if lo == x0 else y), 1), round(a-(x if lo == x0 else y), 1))
                    cur = max(cur, b)
                return [round(best, 1), at]
            # exclude the letter slot so the letter does not hide an empty column/row next to it
            body = [i['box'] for i in self.items if i['panel'] == L and not (i['text'] == L.lower() and i['size'] == 11.0)]
            res['internal_vertical_strip'][L] = widest([(b[0], b[2]) for b in body], x0, x1)
            res['internal_horizontal_strip'][L] = widest([(b[1], b[3]) for b in body], y0, y1)
        Ls = list(self.panels); P = self.panels
        for i, A in enumerate(Ls):
            for B in Ls[i+1:]:
                (ax, ay, aw, ah), (bx_, by, bw, bh) = P[A], P[B]; a, b = ink[A], ink[B]
                if min(ay+ah, by+bh)-max(ay, by) > 5 and abs(abs(bx_-(ax+aw)) - GAP) < 1.0 or abs(abs(ax-(bx_+bw)) - GAP) < 1.0 and min(ay+ah, by+bh)-max(ay, by) > 5:
                    res['ink_gaps'][f'{A}-{B} (h)'] = round(max(b[0]-a[2], a[0]-b[2]), 1)
                elif min(ax+aw, bx_+bw)-max(ax, bx_) > 5 and (abs(by-(ay+ah)-GAP) < 1.0 or abs(ay-(by+bh)-GAP) < 1.0):
                    res['ink_gaps'][f'{A}-{B} (v)'] = round(max(b[1]-a[3], a[1]-b[3]), 1)
        return res

    def save(self, outdir, data, extra=None):
        outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
        chk = self.check()
        pdf = outdir/f'{self.name}.pdf'
        self.t.save(pdf)
        try:
            from pdfclean import finalize
            finalize(pdf, label=self.name, log=outdir/'sizes.tsv')
        except Exception as e:                                    # clean-up is cosmetic; record and continue
            chk['pdfclean_error'] = repr(e)
        subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(pdf), str(outdir/self.name)], check=True,
                       capture_output=True)
        srcs = data if isinstance(data, (list, tuple)) else [data]
        rec = {}
        for d_ in srcs: rec.update(d_.receipt())
        data = srcs[0]
        rec.update(figure=self.name, page_pt=[self.W, self.H], panels=self.panels, checks=chk,
                                         design_tokens=dict(LW=LW, MS=MS, MEW=MEW, ALPHA=ALPHA, used=TOKENS_USED),
                                         notes=self.notes, pdf_sha256=digest(pdf),
                                         colour_scales=dict(diverging=DIV_NAME, sequential=SEQ_NAME+' '+'/'.join(SEQ_STOPS)))
        if extra: rec.update(extra)
        (outdir/f'{self.name}.receipt.json').write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
        (outdir/f'{self.name}.values.json').write_text(json.dumps(self.values, indent=1, ensure_ascii=False, default=str))
        bad = len(chk['outside'])+len(chk['text_overlaps'])+len(chk['letter_zone'])+len(chk['line_aspect_bad'])+len(chk['colour_unnamed'])
        print(self.name, 'outside', len(chk['outside']), 'overlaps', len(chk['text_overlaps']), 'letter-zone', len(chk['letter_zone']), 'line-aspect', len(chk['line_aspect_bad']), 'sizes', chk['sizes'],
              'margins', chk['margins'])
        print('  slack', chk['slack']); print('  v-strip', chk['internal_vertical_strip'], 'h-strip', chk['internal_horizontal_strip'])
        print('  ink gaps', chk['ink_gaps'])
        if chk['outside']: print('  OUTSIDE', [(o.get('panel'), o.get('text'), o.get('box')) for o in chk['outside']][:8])
        if chk['text_overlaps']: print('  OVERLAPS', chk['text_overlaps'][:8])
        if chk['letter_zone']: print('  LETTER-ZONE', chk['letter_zone'][:12])
        if chk['colour_unnamed']: print('  COLOUR-UNNAMED', chk['colour_unnamed'][:12])
        print('  line aspect', chk['line_aspect'])
        return chk
