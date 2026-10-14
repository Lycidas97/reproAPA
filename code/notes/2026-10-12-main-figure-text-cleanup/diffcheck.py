"""Pixel diff production vs candidate at 150 dpi (top-origin; a taller candidate is compared on the production
height and its extra bottom strip must be blank): connected changed regions in pt."""
import subprocess, tempfile, json
from pathlib import Path
import numpy as np
from PIL import Image
res = {}; k = 150/72
for n in (1, 2, 3, 4):
    d = Path(tempfile.mkdtemp())
    for tag, p in [('a', f'../../Fig/Figure.{n}.pdf'), ('b', f'out/Figure.{n}.pdf')]:
        subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', p, str(d/tag)], check=True, capture_output=True)
    A = np.asarray(Image.open(d/'a.png').convert('RGB')).astype(int); B = np.asarray(Image.open(d/'b.png').convert('RGB')).astype(int)
    extra = None
    if B.shape[0] > A.shape[0]:
        strip = B[A.shape[0]:]; extra = dict(rows=int(strip.shape[0]), nonwhite=int((strip.sum(2) < 750).sum())); B = B[:A.shape[0]]
    m = (np.abs(A-B).sum(2) > 30)
    g = int(10*k); H, W = m.shape; gh, gw = (H+g-1)//g, (W+g-1)//g
    cell = np.zeros((gh, gw), bool); ys, xs = np.nonzero(m); cell[ys//g, xs//g] = True
    seen = np.zeros_like(cell); boxes = []
    for i in range(gh):
        for j in range(gw):
            if cell[i, j] and not seen[i, j]:
                st = [(i, j)]; seen[i, j] = True; cs = []
                while st:
                    a, b = st.pop(); cs.append((a, b))
                    for da in (-1, 0, 1):
                        for db in (-1, 0, 1):
                            u, v = a+da, b+db
                            if 0 <= u < gh and 0 <= v < gw and cell[u, v] and not seen[u, v]:
                                seen[u, v] = True; st.append((u, v))
                r0 = min(a for a, _ in cs); r1 = max(a for a, _ in cs)+1; c0 = min(b for _, b in cs); c1 = max(b for _, b in cs)+1
                boxes.append([round(c0*g/k), round(r0*g/k), round(c1*g/k), round(r1*g/k)])
    res[n] = dict(changed_px=int(m.sum()), regions=boxes, extra_bottom=extra)
    print(n, int(m.sum()), boxes, extra or '')
json.dump(res, open('out/pixel-diff.json', 'w'), indent=1)
