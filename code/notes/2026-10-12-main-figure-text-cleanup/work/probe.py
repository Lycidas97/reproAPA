import sys, pdfplumber
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
src, x0, y0, x1, y1 = sys.argv[1], *map(float, sys.argv[2:6])
pg = pdfplumber.open(ROOT/src).pages[0]
c = pg.crop((x0, y0, x1, y1), strict=False)
for w in c.extract_words(extra_attrs=['size'], x_tolerance=1.5):
    print(f"W {w['text']!r:24} x {w['x0']:.1f}-{w['x1']:.1f} y {w['top']:.1f}-{w['bottom']:.1f} s{w['size']:.1f}")
for kind in ('rects', 'curves', 'lines'):
    for o in getattr(c, kind):
        if o['x1']-o['x0'] > 200 and o['bottom']-o['top'] > 100: continue
        print(kind[0].upper(), f"x {o['x0']:.1f}-{o['x1']:.1f} y {o['top']:.1f}-{o['bottom']:.1f} fill={o.get('non_stroking_color')} st={o.get('stroking_color')}")
