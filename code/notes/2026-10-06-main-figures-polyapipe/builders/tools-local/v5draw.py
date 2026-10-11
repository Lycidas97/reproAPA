"""Shared helpers for the v5 in-place main-figure candidates (records adapter copy; no legacy data roots).

Same drawing path as notes/2026-10-02-main-figure-revision: remove the replaced panel
regions from the locked master PDF with pdfedit.remove_regions and overlay new vector
content drawn with sheet.Sheet in PDF points. All coordinates here are TOP-ORIGIN
(pdfplumber convention, the same as layout-lock.json); conversion happens inside.
"""
import json, sys, hashlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'Fig').is_dir() and (p/'main.tex').exists())   # .../reproAPA, wherever this folder sits
TOOLS = ROOT/'notes/2026-10-02-main-figure-revision/tools'
sys.path.insert(0, str(TOOLS))
from sheet import Sheet, textwidth, FONTS   # noqa: E402
from pdfedit import remove_regions, overlay  # noqa: E402
import pikepdf, pdfplumber                  # noqa: E402

# Legacy statistics roots (IR/V2/V3) are intentionally absent: builders read result records only.

INK = '#263D4A'; MUTED = '#586870'; AXC = '#87969E'; NULLBAR = '#B5C0C7'; NULLW = '#758791'
MISSING = '#DEDEDE'; BAND = '#F3F5F6'
DIV = ['#0072B2', '#F7F7F7', '#D55E00']


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(8 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def num(v, d=3, plus=False):
    s = f'{v:+.{d}f}' if plus else f'{v:.{d}f}'
    return s.replace('-', '−')


def qfmt(q):
    """Two significant digits, leading zero dropped (house style '.0032', '.22')."""
    if q >= 0.1:
        s = f'{q:.2f}'
    elif q >= 0.01:
        s = f'{q:.3f}'
    else:
        s = f'{q:.4f}'
    return s.replace('0.', '.', 1)


def diverge(v, lim):
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list('d', DIV)
    return cm((np.clip(v, -lim, lim)+lim)/(2*lim))


_CAL = {}
def _calibrate():
    """Measured offset between a baseline and pdfplumber's text 'top' for each size."""
    if _CAL:
        return _CAL
    import tempfile
    s = Sheet(200, 100)
    sizes = [6.8, 7.0, 7.2, 7.5, 8.0, 11.0]
    for i, z in enumerate(sizes):
        s.text(10+i*35, 50, 'Xg0', z)
    p = Path(tempfile.mkdtemp())/'cal.pdf'; s.save(p)
    words = pdfplumber.open(p).pages[0].extract_words(extra_attrs=['size'])
    for w in words:
        _CAL[round(w['size'], 1)] = w['top']-(100-50)   # top-origin baseline is 50
    return _CAL


class Top:
    """Sheet wrapper in top-origin coordinates. text_top places by pdfplumber bbox top."""
    def __init__(self, W, H):
        self.W, self.H = W, H; self.s = Sheet(W, H); self.ax = self.s.ax
        self.cal = _calibrate(); self.log = []

    def text(self, x, ybase, s, size=8, font='Reg', color=INK, ha='left', **kw):
        self.log.append(dict(text=s, x=x, ybase=ybase, size=size))
        return self.s.text(x, self.H-ybase, s, size, font, color, ha, 'baseline', **kw)

    def text_top(self, x, top, s, size=8, font='Reg', color=INK, ha='left', **kw):
        off = self.cal.get(round(size, 1), -0.95*size)
        return self.text(x, top-off, s, size, font, color, ha, **kw)

    def vtext(self, x, ycenter, s, size=7.5, color=INK):
        """Rotated (bottom-to-top) text centred on (x, ycenter)."""
        return self.s.text(x, self.H-ycenter, s, size, 'Reg', color, 'center', 'center', rotation=90)

    def line(self, xs, ys, lw=0.5, color=INK, **kw):
        self.s.line(list(xs), [self.H-v for v in ys], lw, color, **kw)

    def rect(self, x, ytop, w, h, fc='white', ec='none', lw=0, **kw):
        self.s.rect(x, self.H-ytop-h, w, h, fc=fc, ec=ec, lw=lw, **kw)

    def marker(self, x, y, m='o', ms=3.3, fc=INK, ec=INK, mew=0.7, zorder=3):
        self.ax.plot([x], [self.H-y], marker=m, ms=ms, mfc=fc, mec=ec, mew=mew, ls='none', zorder=zorder)

    def whisker(self, x0, x1, y, color=NULLW, lw=0.9, cap=1.4):
        self.line([x0, x1], [y, y], lw=lw, color=color)
        for v in (x0, x1):
            self.line([v, v], [y-cap, y+cap], lw=lw, color=color)

    def vwhisker(self, x, y0, y1, color=NULLW, lw=0.9, cap=2.0):
        self.line([x, x], [y0, y1], lw=lw, color=color)
        for v in (y0, y1):
            self.line([x-cap, x+cap], [v, v], lw=lw, color=color)

    def save(self, p, dpi=600):
        import matplotlib.pyplot as plt
        self.s.fig.savefig(p, transparent=True, dpi=dpi); plt.close(self.s.fig)


def clear(pdf, H, rects_top, center=False):
    rects = [(x0, H-y1, x1, H-y0) for x0, y0, x1, y1 in rects_top]
    return remove_regions(pdf, rects, center=center)


def render_png(pdf_path, png_path, dpi=150, boxes=None, H=None):
    import subprocess
    stem = str(png_path)[:-4]
    subprocess.run(['pdftoppm', '-r', str(dpi), '-png', '-singlefile', str(pdf_path), stem], check=True)
    if boxes:
        from PIL import Image, ImageDraw
        im = Image.open(png_path).convert('RGB'); d = ImageDraw.Draw(im); k = dpi/72
        for x0, y0, x1, y1 in boxes:
            d.rectangle([x0*k, y0*k, x1*k, y1*k], outline=(220, 30, 30), width=1)
        im.save(str(png_path)[:-4]+'.tiling.png')


def prune_xobjects(pdf_):
    """Drop /XObject resource entries that no content stream (page or nested form) still draws with Do."""
    seen = set()
    def walk(obj):
        if obj.objgen in seen and obj.objgen != (0, 0):
            return
        seen.add(obj.objgen)
        res = obj.get('/Resources')
        xo = res.get('/XObject') if res is not None else None
        try:
            ops = pikepdf.parse_content_stream(obj)
        except Exception:
            return
        used = {str(o[0][0]) for o in ops if str(o[1]) == 'Do'}
        if xo is None:
            return
        for k in list(xo.keys()):
            if k not in used:
                del xo[k]
            elif xo[k].get('/Subtype') == '/Form':
                walk(xo[k])
    for pg_ in pdf_.pages:
        walk(pg_.obj)


def shift_block(pdf, master, H, W, pre_clear, region, dx=0.0, dy=0.0, scale=1.0):
    """Move master content lying inside `region` (top-origin x0, y0, x1, y1) by (dx, dy), optionally scaled about
    the region's top-left corner. A copy of the master with the same `pre_clear` removals is cropped to the region
    (objects fully outside are deleted, so nothing rides along hidden) and placed as a clipped form; the originals
    inside the region are removed from `pdf`. Returns (removed_in_page, kept_in_copy)."""
    x0, y0, x1, y1 = region
    blk = pikepdf.open(master)
    if pre_clear:
        clear(blk, H, pre_clear)
    kept = clear(blk, H, [(-20, -20, x0-0.05, H+20), (x1+0.05, -20, W+20, H+20), (-20, -20, W+20, y0-0.05), (-20, y1+0.05, W+20, H+20)])
    prune_xobjects(blk)
    bp = blk.pages[0]; box = [x0, H-y1, x1, H-y0]; bp.MediaBox = box; bp.CropBox = box
    n = len(clear(pdf, H, [(x0-0.3, y0-0.3, x1+0.3, y1+0.3)]))
    w, h = (x1-x0)*scale, (y1-y0)*scale
    pdf.pages[0].add_overlay(bp, pikepdf.Rectangle(x0+dx, H-(y0+dy+h), x0+dx+w, H-(y0+dy)))
    return n


INFO = {}
def info(panel, mark, n, channels):
    """Information count bookkeeping: `n` plotted marks of one kind, each encoding `channels` data channels."""
    INFO.setdefault(panel, []).append(dict(mark=mark, n=int(n), channels=int(channels)))


def save_info(path):
    json.dump(INFO, open(path, 'w'), indent=1)
