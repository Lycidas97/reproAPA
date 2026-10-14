"""Objective crowding scan.
H: numeric labels (ticks, counts, q) on one baseline closer than 4 pt (space ≈ 2 pt counts as merged).
V: vertically stacked text (x-overlap > 1 pt) with clearance < 1.0 pt between glyph boxes."""
import sys, re, pdfplumber
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
src = sys.argv[1] if len(sys.argv) > 1 else 'Fig'
NUM = re.compile(r'^[−+\-–]?[\d.,/%]+(pp)?$|^[−+]?\.\d+$')
for n in (1, 2, 3, 4):
    pg = pdfplumber.open(ROOT/src/f'Figure.{n}.pdf').pages[0]
    ws = [w for w in pg.extract_words(x_tolerance=1.0, extra_attrs=['size', 'upright']) if w['upright']]
    print(f'##### Figure {n}')
    for i, a in enumerate(ws):
        for b in ws[i+1:]:
            if abs(a['top']-b['top']) < 1.0:
                if NUM.match(a['text']) and NUM.match(b['text']):
                    l, r = (a, b) if a['x0'] < b['x0'] else (b, a)
                    g = r['x0']-l['x1']
                    if g < 4.0: print(f"  H {g:4.1f}pt y={a['top']:6.1f}  '{l['text']}' '{r['text']}'  x={l['x0']:.0f}")
            else:
                ox = min(a['x1'], b['x1'])-max(a['x0'], b['x0'])
                if ox > 1:
                    up, lo = (a, b) if a['top'] < b['top'] else (b, a)
                    g = lo['top']-up['bottom']
                    if -3 < g < 1.0: print(f"  V {g:4.1f}pt x={max(a['x0'],b['x0']):6.1f}  '{up['text']}'(y{up['top']:.0f}) / '{lo['text']}'(y{lo['top']:.0f})")
