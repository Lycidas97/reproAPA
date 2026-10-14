"""Supplementary Figure S20 (v11): state-resolved APA profile classes and partial endpoint retention.
a–c former main Figure 2 a–c (Fig/Figure.2.pdf), vector crop, letters unchanged
d   former Supplementary Figure S7c (Fig/SuppFig.S7.pdf = consolidated G), vector crop, relettered c -> d
Former main Figure 2 d–f and former S6, S7 a,b,d,e, S8 and S11 are not carried into the v11 supplement (author decision
2026-10-10: one condensed profile figure supporting the Discussion statement on endpoint retention).
No data are redrawn; the receipt records source SHA-256 and regions."""
import json, subprocess, tempfile
from pathlib import Path
from common import *                                     # noqa
from common import supdraw, name, OUT, ROOT, crop, place, new_page, sha
import pikepdf

F2 = ROOT/'revision-backups/20261011-before-v11-apply/Fig/Figure.2.pdf'; H2 = 556.74
import os
R2 = (0.0, 0.0, PAGE_W, 322.0)                  # former Figure 2 rows a–c incl. their keys
# 2026-10-12 text cleanup: the formula block of former Figure 2b (L, D, T definitions) is removed (now in the caption), and
# d is cropped from a composed page holding only former S7 c (consolidated G, D07 C) rebuilt without its column title (V11_S20_D).
F2_CLEAR = [(121.0, 188.8, 215.0, 218.6),       # 'L = |R−S| + |E−R| + |P−E|', 'D = |P−S|', 'T = 1 − D/L'
            (356.0, 187.0, 453.0, 198.3),       # 2026-10-12 S20c: '|ΔDWUI| / L' and 'Median D/L' ran together (2.5 pt gap); redrawn below
            (412.0, 274.0, 417.8, 283.5)]       # 2026-10-12 S20c: dot-axis '0' touched the bar-axis '1' (2.1 pt gap); '.5' drawn at the dotted line instead
S20C_HDR_TOP = float(os.environ.get('S20C_HDR_TOP', '189.1'))
if os.environ.get('V11_S20_D'):
    S7 = Path(os.environ['V11_S20_D'])
    with pikepdf.open(S7) as _p: H7 = float(_p.pages[0].MediaBox[3])
    R7 = (0.0, 4.0, PAGE_W, H7-4.5); L7 = (4.0, 4.0, 13.5, 16.5)
else:
    S7 = ROOT/'revision-backups/20261011-before-v11-apply/Fig/SuppFig.S7.pdf'; H7 = 384.4
    R7 = (0.0, 125.8, PAGE_W, 273.6); L7 = (4.0, 126.0, 13.5, 138.5)
DY7 = 322.0+GAP-(R7[1]+0.5)
H = R7[3]+DY7+4.5
dst = new_page(H)
place(dst, H, crop(F2, H2, R2, extra_clear=F2_CLEAR), R2, 0.0, 0.0)
place(dst, H, crop(S7, H7, R7, extra_clear=[L7]), R7, 0.0, R7[1]+DY7)
F = supdraw.Fig(H, name('profile'))
F.t.text_top(4.5+0.2, R7[1]+0.5+0.4+DY7, 'd', 11.0, font='Bol')
# S20c column headers (2026-10-12): bar block 358.2–408.3 pt, dot column 414.9–444.4 pt in the former Figure 2 coordinates.
# The bar header loses its spaces and is centred on the bars; the dot-column header is right-aligned to the 4.5 pt margin.
F.text((358.2+408.3)/2, S20C_HDR_TOP, '|ΔDWUI|/L', 8.0, ha='center')
F.text(PAGE_W-4.5, S20C_HDR_TOP, 'Median D/L', 8.0, ha='right')
F.text((414.9+444.4)/2, 274.7, '.5', 8.0, ha='center')     # dot column 0–1 spans 414.9–444.4 pt; the dotted line sits at .5
tmp = Path(tempfile.mkdtemp()); F.t.save(tmp/'ov.pdf')
ov = pikepdf.open(tmp/'ov.pdf')
dst.pages[0].add_overlay(ov.pages[0], pikepdf.Rectangle(0, 0, PAGE_W, H))
OUT.mkdir(exist_ok=True)
out = OUT/f'{F.name}.pdf'; dst.save(out)
try:
    from pdfclean import finalize
    finalize(out, label=F.name, log=OUT/'sizes.tsv')
except Exception as e:
    print('pdfclean', e)
subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(out), str(OUT/F.name)], check=True, capture_output=True)
def rel(p): return os.path.relpath(Path(p).resolve(), ROOT.resolve())
rec = dict(figure=F.name, page_pt=[PAGE_W, H], pdf_sha256=sha(out), redrawn=False,
           kept=[dict(source=rel(F2), source_sha=sha(F2), region=R2, dy=0.0, letters='a, b, c unchanged', cleared=F2_CLEAR),
                 dict(source=rel(S7), source_sha=sha(S7), region=R7, dy=DY7, letters={'c': 'd'})],
           not_carried=['former Figure 2 d–f', 'former S6', 'former S7 a, b, d, e', 'former S8', 'former S11'])
(OUT/f'{F.name}.receipt.json').write_text(json.dumps(rec, indent=1))
print(F.name, 'H', round(H, 1), 'mm', round(H/72*25.4, 1))
