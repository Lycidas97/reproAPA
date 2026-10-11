"""Shared set-up for the v11 supplementary figures (2026-10-10).
Same house style as the consolidated supplement: notes/2026-10-06-supplement-polyapipe/tools/supdraw.py (page 453.5 pt,
4.5 pt margin, 11.3 pt panel gap, size/line/marker tokens, named colours, layout checks, receipts).
New analyses of 2026-10-10 are read directly from their notes/ folders; every read is SHA-256 recorded in the receipt.
Printed S numbers come only from NUM below (v11 order, reviews/2026-10-10-divergence-writing/supplement-plan.md)."""
import sys, json, hashlib
from pathlib import Path
import pandas as pd
FIGDIR = Path(__file__).resolve().parent
ROOT = FIGDIR.parents[2]
sys.path.insert(0, str(ROOT/'notes/2026-10-06-supplement-polyapipe/tools'))
import supdraw                                           # noqa: E402
from supdraw import *                                    # noqa: E402,F401
OUT = FIGDIR.parent/'out'            # defined after the star import (supdraw has its own HERE/OUT)
# v11 numbering of the new figures (content code -> printed number)
NUM = {'strat': 10, 'naive': 12, 'cycle': 18, 'mech': 19, 'profile': 20}


def name(code):
    return f'SuppFig.S{NUM[code]}'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


class Reads:
    """Direct reads of notes/ tables (outside the r2 contract bundle), each SHA-256 recorded."""
    def __init__(self):
        self.reads = {}

    def df(self, rel, **kw):
        p = ROOT/rel
        self.reads[rel] = sha(p)
        kw.setdefault('sep', '\t')
        return pd.read_csv(p, **kw)

    def receipt(self):
        return dict(direct_reads=self.reads, statistics_recomputed_by_renderer=False, production_replaced=False,
                    note='new 2026-10-10 analyses read from their notes/ folders; not yet registered in a plot contract')


# ------------------------------------------------------------------ vector reuse of kept panels (as in the v11 main figures)
import subprocess, tempfile
sys.path.insert(0, str(ROOT/'notes/2026-10-06-main-figures-polyapipe/builders/tools-local'))
import pikepdf                                            # noqa: E402
from v5draw import clear, prune_xobjects                  # noqa: E402


def crop(src, H, region, extra_clear=(), W=PAGE_W):
    """Vector crop of region (x0, y0, x1, y1; top origin) of a one-page PDF of height H; extra_clear boxes are blanked."""
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


def new_page(H, W=PAGE_W):
    p = pikepdf.new(); p.add_blank_page(page_size=(W, H)); return p


def finish(dst, F, H, data, extra):
    """Overlay the newly drawn panels (F) on the page with the kept crops, then write PDF, PNG, receipt and values."""
    import json as _json
    chk = F.check()
    tmp = Path(tempfile.mkdtemp()); F.t.save(tmp/'overlay.pdf')
    ov = pikepdf.open(tmp/'overlay.pdf')
    dst.pages[0].add_overlay(ov.pages[0], pikepdf.Rectangle(0, 0, PAGE_W, H))
    OUT.mkdir(exist_ok=True)
    out = OUT/f'{F.name}.pdf'; dst.save(out)
    try:
        from pdfclean import finalize
        finalize(out, label=F.name, log=OUT/'sizes.tsv')
    except Exception as e:
        extra['pdfclean_error'] = repr(e)
    subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(out), str(OUT/F.name)], check=True, capture_output=True)
    rec = data.receipt()
    rec.update(figure=F.name, page_pt=[PAGE_W, H], panels_new=F.panels, checks=chk, pdf_sha256=sha(out),
               design_tokens=dict(LW=supdraw.LW, MS=supdraw.MS, MEW=supdraw.MEW, ALPHA=supdraw.ALPHA, used=supdraw.TOKENS_USED))
    rec.update(extra)
    (OUT/f'{F.name}.receipt.json').write_text(_json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    (OUT/f'{F.name}.values.json').write_text(_json.dumps(F.values, indent=1, ensure_ascii=False, default=str))
    bad = {k: len(chk[k]) for k in ('outside', 'text_overlaps', 'letter_zone', 'line_aspect_bad', 'colour_unnamed')}
    print(F.name, bad, 'slack', chk['slack']); print('  v', chk['internal_vertical_strip'], 'h', chk['internal_horizontal_strip'])
    for k in ('outside', 'text_overlaps', 'letter_zone'):
        if chk[k]: print(' ', k, chk[k][:6])
    return chk
