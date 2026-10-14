"""Pixel checks (144 dpi, so the 7.5 pt move is a whole number of pixels): changed regions per figure; for S19 the region below the move is compared after
shifting back by 7.5 pt, so only the intended changes remain."""
import subprocess, tempfile, json
from pathlib import Path
import numpy as np
from PIL import Image
HERE = Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]; DPI = 144; k = DPI/72   # 7.5 pt = 15 px exactly
def img(p):
    o = Path(tempfile.mkdtemp())/'i'
    subprocess.run(['pdftoppm', '-r', str(DPI), '-png', '-singlefile', str(p), str(o)], check=True, capture_output=True)
    return np.asarray(Image.open(str(o)+'.png').convert('L')).astype(int)
def regions(m, cell=12):
    ys, xs = np.nonzero(m); cells = sorted({(int(y/k//cell), int(x/k//cell)) for y, x in zip(ys, xs)}); boxes = []
    for r, c in cells:
        for b in boxes:
            if b[0]-1 <= c <= b[2]+1 and b[1]-1 <= r <= b[3]+1:
                b[0], b[1], b[2], b[3] = min(b[0], c), min(b[1], r), max(b[2], c), max(b[3], r); break
        else: boxes.append([c, r, c, r])
    return [[b[0]*cell, b[1]*cell, (b[2]+1)*cell, (b[3]+1)*cell] for b in boxes]
res = {}
for f in ('Figure.1', 'Figure.4'):
    A, B = img(ROOT/f'Fig/{f}.pdf'), img(HERE/f'out/{f}.pdf'); m = np.abs(A-B) > 40
    res[f] = regions(m); print(f, int(m.sum()), res[f])
A, B = img(ROOT/'Fig/SuppFig.S19.pdf'), img(HERE/'out/SuppFig.S19.pdf')
cut = int(round(289.9*k)); sh = int(round(7.5*k))
top = np.abs(A[:cut]-B[:cut]) > 40
bot = np.abs(A[cut:]-B[cut+sh:cut+sh+A.shape[0]-cut]) > 40
res['S19'] = dict(top=regions(top), below_after_unshift=int(bot.sum()))
print('S19 top', int(top.sum()), regions(top), '| below the cut, after undoing the 7.5 pt move:', int(bot.sum()), 'px')
json.dump(res, open(HERE/'out/verify.json', 'w'), indent=1)
