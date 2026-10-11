"""Figure 5f: remove the second placement of the panel-letter crop (/BKj… inside /t5wOQ…).

The master stacks the 5f header group twice: the full group (cm 0, +0.037) and a crop of the letter region
(cm 0, -7.403). The crop repeats the 'f' at the same spot (hidden under the first one) but brings its white
letter knockout 7.44 pt lower, where it cuts a 0.48 pt notch into the top of the Guo2019-CIP colour bar
(x 5.0-11.5 pt, top 398.16-398.64 pt). Dropping that Do call leaves the letter and removes the notch.
Checks: exactly one Do removed; render differs only inside the notch box; 'f' glyph still present; then pdfclean.finalize.
"""
import sys, subprocess, tempfile
from pathlib import Path
import numpy as np, pikepdf, pdfplumber
from PIL import Image
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'notes/2026-10-06-main-figures-polyapipe/builders/tools-local'))
import pdfclean
SRC = ROOT/'Fig/Figure.5.pdf'; OUT = HERE/'final/Figure.5.pdf'

def render(p, dpi=300):
    t = Path(tempfile.mkdtemp())/'r'
    subprocess.run(['pdftoppm', '-r', str(dpi), '-png', '-singlefile', str(p), str(t)], check=True, capture_output=True)
    return np.asarray(Image.open(str(t)+'.png').convert('RGB'), dtype=int)

pdf = pikepdf.open(SRC); xo = pdf.pages[0].Resources.XObject
host = [n for n in xo.keys() if n.startswith('/t5wOQ')]; assert len(host) == 1
form = xo[host[0]]; ops = pikepdf.parse_content_stream(form)
drop = [i for i, (o, op) in enumerate(ops) if str(op) == 'Do' and str(o[0]).startswith('/BKjCUTUAw')]
assert len(drop) == 1, drop
form.write(pikepdf.unparse_content_stream([x for i, x in enumerate(ops) if i not in drop]))
pdf.save(OUT)

a, b = render(SRC), render(OUT)
ys, xs = np.nonzero(np.abs(a-b).sum(2))
box = [round(v*72/300, 2) for v in (xs.min(), ys.min(), xs.max()+1, ys.max()+1)]
print('changed px @300dpi', len(xs), 'box pt', box)
assert box[0] >= 4.5 and box[2] <= 12.0 and box[1] >= 397.5 and box[3] <= 399.0, box
pg = pdfplumber.open(OUT).pages[0]
assert [c['text'] for c in pg.chars if c['x0'] < 20 and 380 < c['top'] < 400] == ['f']
assert not [r for r in pg.rects if r['non_stroking_color'] in [(1, 1, 1), (1.0, 1.0, 1.0)] and r['bottom'] > 398.16 and r['x0'] < 12 and r['top'] < 401.16]
print(pdfclean.finalize(OUT, label='Figure.5.pdf (f-letter knockout fix)'))
c = render(OUT)
assert int(np.abs(b-c).max()) == 0
print('ok')
