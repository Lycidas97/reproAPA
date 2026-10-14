"""Candidates vs production: (1) crowding/overlap scan at the edited spots, (2) 150-dpi pixel diff regions,
(3) page size unchanged, (4) word inventory difference."""
import re, json, subprocess, tempfile, collections
from pathlib import Path
import numpy as np, pdfplumber
from PIL import Image
HERE = Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]
NUM = re.compile(r'^[−+\-–]?[\d.,/%×]+(pp)?$|^[−+]?\.\d+$')
def scan(p):
    pg = pdfplumber.open(p).pages[0]
    up = [w for w in pg.extract_words(x_tolerance=1.0, extra_attrs=['size', 'upright']) if w['upright']]
    out = set()
    for i, a in enumerate(up):
        for b in up[i+1:]:
            if abs(a['top']-b['top']) >= 1.0:
                ox = min(a['x1'], b['x1'])-max(a['x0'], b['x0'])
                if ox > 1:
                    u, d = (a, b) if a['top'] < b['top'] else (b, a); g = d['top']-u['bottom']
                    if -3 < g < 1.0: out.add(f"V {g:.1f} '{u['text']}'/'{d['text']}'")
    return out, (float(pg.width), float(pg.height))
res = {}; k = 150/72
for p in sorted((HERE/'out').glob('SuppFig.S*.pdf'), key=lambda q: int(q.stem.split('S')[-1])):
    f = p.stem; src = ROOT/f'Fig/{f}.pdf'
    a, sa = scan(src); b, sb = scan(p)
    d = Path(tempfile.mkdtemp())
    for tag, q in [('a', src), ('b', p)]:
        subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(q), str(d/tag)], check=True, capture_output=True)
    A = np.asarray(Image.open(d/'a.png').convert('RGB')).astype(int); B = np.asarray(Image.open(d/'b.png').convert('RGB')).astype(int)
    m = np.abs(A-B).sum(2) > 30; ys, xs = np.nonzero(m)
    # cluster changed pixels into 12-pt cells
    cells = sorted({(int(y/k//12), int(x/k//12)) for y, x in zip(ys, xs)})
    boxes = []
    for r, c in cells:
        for bx in boxes:
            if bx[0]-1 <= c <= bx[2]+1 and bx[1]-1 <= r <= bx[3]+1:
                bx[0], bx[1], bx[2], bx[3] = min(bx[0], c), min(bx[1], r), max(bx[2], c), max(bx[3], r); break
        else: boxes.append([c, r, c, r])
    regions = [[b0*12, b1*12, (b2+1)*12, (b3+1)*12] for b0, b1, b2, b3 in boxes]
    res[f] = dict(size_same=sa == sb, changed_px=int(m.sum()), regions=regions,
                  crowd_fixed=sorted(a-b), crowd_new=sorted(b-a))
    print(f, 'size_same', sa == sb, 'px', int(m.sum()), 'regions(pt)', regions)
    if a-b: print('   resolved:', sorted(a-b)[:8])
    if b-a: print('   NEW:', sorted(b-a)[:8])
json.dump(res, open(HERE/'out/verify.json', 'w'), indent=1, ensure_ascii=False)
