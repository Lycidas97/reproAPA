"""List painted objects and words in the moving area that are not fully inside exactly one region."""
import sys, pdfplumber
sys.path.insert(0, '.')
from pathlib import Path
import plan
ROOT = Path('.').resolve().parents[1]
pg = pdfplumber.open(ROOT/plan.SRC).pages[0]
R = plan.regions()
def inside(o, r): return o['x0'] >= r[0]-0.01 and o['x1'] <= r[2]+0.01 and o['top'] >= r[1]-0.01 and o['bottom'] <= r[3]+0.01
moving = lambda o: (o['top'] > 197 and o['x0'] >= 296) or o['top'] > 289
bad = []
for o in pg.rects+pg.lines+pg.curves+pg.extract_words(x_tolerance=1.5):
    if not moving(o): continue
    hits = [r[2] for r in R if inside(o, r[0])]
    if len(hits) != 1:
        bad.append((o.get('object_type', 'word'), round(o['x0'], 1), round(o['x1'], 1), round(o['top'], 2), round(o['bottom'], 2), o.get('text', ''), hits))
print('page height', pg.height, '| objects not in exactly one region:', len(bad))
for b in bad: print('  ', b)
