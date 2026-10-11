"""Panel containment and density audit for the v5 candidates.

containment: every text word, rect, line, curve and image is assigned to the panel holding its centre;
             it is flagged when its box crosses that panel's rectangle (tolerance 0.2 pt) and listed as
             'tight' when it comes within 1.5 pt of an internal boundary.
density:     rendered at 72 dpi; per panel, ink fraction (non-near-white pixels) and empty fraction
             (share of 12 x 12 pt cells with no ink, after removing the panel-letter cell).
Usage: python3 audit.py [figure numbers...]
"""
import sys, json, subprocess, tempfile
from pathlib import Path
import numpy as np
import pdfplumber
from PIL import Image

HERE = Path(__file__).resolve().parent
PANELS = {   # shared boundaries 1.5 pt left of / 1.0 pt above the locked panel-letter boxes
    1: dict(A=(0, 0, 453.5, 320.1), B=(0, 320.1, 226.5, 464.8), C=(0, 464.8, 226.5, 599.8), D=(226.5, 320.1, 453.5, 599.8)),
    2: dict(A=(0, 0, 453.5, 231.7), B=(0, 231.7, 232.4, 372.5), C=(232.4, 231.7, 453.5, 372.5),
            D=(0, 372.5, 453.5, 491.6), E=(0, 491.6, 329.6, 613.3), F=(329.6, 491.6, 453.5, 613.3)),
    3: dict(A=(0, 0, 235.7, 146.1), B=(235.7, 0, 453.5, 146.1), C=(0, 146.1, 453.5, 436.9), D=(0, 436.9, 232.8, 640.8),
            E=(232.8, 436.9, 453.5, 640.8)),
    4: dict(A=(0, 0, 302.8, 172.0), B=(302.8, 0, 453.5, 172.0), C=(0, 172.0, 240.7, 330.0),
            E=(0, 330.0, 240.7, 593.7), D=(240.7, 172.0, 453.5, 593.7)),   # 2026-10-05 polyApipe relayout: C+E | D, no F
    5: dict(A=(0, 0, 333.7, 103.9), B=(0, 103.9, 333.7, 191.6), C=(333.7, 0, 453.5, 191.6), D=(0, 191.6, 239.0, 377.8),
            E=(239.0, 191.6, 453.5, 377.8), F=(0, 377.8, 453.5, 527.2)),
}


def panels_for(fig, pdf=None):
    """Panel rectangles for a figure; a retiled figure carries its moved rectangles in Figure.N.panels.json."""
    import json as _j
    p = HERE.parent/f'figure{fig}'/f'Figure.{fig}.panels.json'
    if pdf is not None and Path(pdf).name == f'Figure.{fig}.pdf' and p.exists():
        return {k: tuple(v) for k, v in _j.load(open(p)).items()}
    return PANELS[fig]


def objects(pg):
    out = []
    for w in pg.extract_words(keep_blank_chars=False, x_tolerance=1.2, extra_attrs=['size']):
        if w['bottom'] < 0 or w['top'] > pg.height:
            continue
        if w['size'] >= 10.0 and len(w['text']) == 1 and w['text'].isupper():
            continue                     # locked panel letters (their font box overhangs the anchor box)
        if w['size'] < 6.0:          # rotated labels split into per-glyph words; checked through chars below
            continue
        out.append(('text', w['text'], w['x0'], w['top'], w['x1'], w['bottom']))
    for c in pg.chars:
        if not c.get('upright', True):
            out.append(('vtext', c['text'], c['x0'], c['top'], c['x1'], c['bottom']))
    W, H = pg.width, pg.height
    for kind in ('rects', 'lines', 'curves', 'images'):
        for o in getattr(pg, kind):
            col = o.get('non_stroking_color')
            if isinstance(col, (int, float)):
                col = (col,)
            if kind in ('rects', 'curves') and not o.get('stroke') and col is not None and len(col) in (1, 3) and min(col) > 0.97:
                continue                  # white backgrounds / halos
            if (o['x1']-o['x0']) > 0.8*W and (o['bottom']-o['top']) > 0.4*H:
                continue                  # page-size clip or background paths
            if o['bottom'] < 0 or o['top'] > H or o['x1'] < 0 or o['x0'] > W:
                continue                  # fully off page (clipped away)
            out.append((kind[:-1], '', o['x0'], o['top'], o['x1'], o['bottom']))
    return out


def containment(fig, pdf):
    pg = pdfplumber.open(pdf).pages[0]; P = panels_for(fig, pdf)
    over, tight = [], []
    for kind, txt, x0, y0, x1, y1 in objects(pg):
        cx, cy = (x0+x1)/2, (y0+y1)/2
        home = [k for k, (a, b, c, d) in P.items() if a <= cx <= c and b <= cy <= d]
        if not home:
            over.append((kind, txt, 'no panel', [round(v, 1) for v in (x0, y0, x1, y1)])); continue
        k = home[0]; a, b, c, d = P[k]
        ex = max(a-x0, b-y0, x1-c, y1-d)
        if ex > 0.2:
            over.append((kind, txt, k, [round(v, 1) for v in (x0, y0, x1, y1)], round(ex, 1)))
        else:
            # distance to boundaries shared with another panel (page edges excluded)
            m = []
            if a > 0.5: m.append(x0-a)
            if b > 0.5: m.append(y0-b)
            if c < pg.width-0.5: m.append(c-x1)
            if d < pg.height-0.5: m.append(d-y1)
            if m and min(m) < 1.0:
                tight.append((kind, txt, k, [round(v, 1) for v in (x0, y0, x1, y1)], round(min(m), 1)))
    return over, tight


def density(fig, pdf):
    tmp = Path(tempfile.mkdtemp())/'r'
    subprocess.run(['pdftoppm', '-r', '72', '-png', '-singlefile', str(pdf), str(tmp)], check=True, capture_output=True)
    im = np.asarray(Image.open(str(tmp)+'.png').convert('L'), dtype=float)
    ink = im < 235
    rgb = Image.open(str(tmp)+'.png').convert('RGB')
    big = Path(tempfile.mkdtemp())/'b'
    subprocess.run(['pdftoppm', '-r', '110', '-png', '-singlefile', str(pdf), str(big)], check=True, capture_output=True)
    from PIL import ImageDraw
    bi = Image.open(str(big)+'.png').convert('RGB'); dr = ImageDraw.Draw(bi); kk = 110/72
    for k, (a, b, c, d) in panels_for(fig, pdf).items():
        dr.rectangle([a*kk, b*kk, c*kk-1, d*kk-1], outline=(220, 30, 30), width=1)
    bi.save(HERE.parent/f'figure{fig}'/f'Figure.{fig}.draft.panels.png')
    res = {}
    for k, (a, b, c, d) in panels_for(fig, pdf).items():
        sub = ink[int(b):int(d), int(a):int(c)]
        cells = []
        for yy in range(0, sub.shape[0]-6, 12):
            for xx in range(0, sub.shape[1]-6, 12):
                if yy < 18 and xx < 18:
                    continue
                cells.append(sub[yy:yy+12, xx:xx+12].any())
        res[k] = dict(ink=round(float(sub.mean()), 3), empty_cells=round(1-float(np.mean(cells)), 3))
    return res


if __name__ == '__main__':
    figs = [int(v) for v in sys.argv[1:]] or [1, 2, 3, 4, 5]
    report = {}
    for f in figs:
        pdf = HERE.parent/f'figure{f}'/f'Figure.{f}.draft.pdf'
        over, tight = containment(f, pdf); dens = density(f, pdf)
        report[f] = dict(overflow=over, tight=tight, density=dens)
        print(f'== Figure {f}')
        for k, v in dens.items():
            print(f'   {k}: ink {v["ink"]:.3f}  empty cells {v["empty_cells"]:.2f}')
        print(f'   overflow {len(over)}')
        for o in over:
            print('     ', o)
        print(f'   tight text {len(tight)}')
        for o in tight:
            print('     ', o)
    json.dump(report, open(HERE.parent/'audit.json', 'w'), indent=1)
