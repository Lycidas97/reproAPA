"""Every painted object and character of the source must reappear in the output at its planned displacement
(same type, size and x; top + dy of its region; dy = 0 outside the moving regions). Redrawn e lines are checked apart."""
import sys, collections, json
sys.path.insert(0, '.')
import pdfplumber, plan
from pathlib import Path
ROOT = Path('.').resolve().parents[1]
src = pdfplumber.open(ROOT/plan.SRC).pages[0]; out = pdfplumber.open('out/SuppFig.S19.pdf').pages[0]
R = plan.regions()
def dy_of(o):
    cx, cy = (o['x0']+o['x1'])/2, (o['top']+o['bottom'])/2
    for (x0, y0, x1, y1), d, _ in R:
        if x0 <= cx <= x1 and y0 <= cy <= y1: return d
    return 0.0
def key(o, d=0.0, t=None): return (t or o['object_type'], round(o['x0'], 1), round(o['x1'], 1), round(o['top']+d, 1), round(o['bottom']+d, 1), o.get('text', ''))
skip = lambda o: (abs(o['x0']-401.42) < 0.05 and o['bottom']-o['top'] > 30) or (abs(o['top']-264.6) < 0.05 and abs(o['x0']-371.8) < 0.05 and o['x1']-o['x0'] > 50)
exp = collections.Counter(); got = collections.Counter()
for kind in ('rects', 'lines', 'curves', 'chars'):
    for o in getattr(src, kind):
        if kind != 'chars' and skip(o): continue
        exp[key(o, dy_of(o))] += 1
    for o in getattr(out, kind): got[key(o)] += 1
missing = exp-got; extra = got-exp
print('objects expected', sum(exp.values()), '| missing', sum(missing.values()), '| extra', sum(extra.values()))
for k, v in list(missing.items())[:12]: print('  missing', k, v)
for k, v in list(extra.items())[:12]: print('  extra  ', k, v)
json.dump(dict(expected=sum(exp.values()), missing=[list(k)+[v] for k, v in missing.items()], extra=[list(k)+[v] for k, v in extra.items()]), open('out/verify.json', 'w'), indent=1, ensure_ascii=False)
