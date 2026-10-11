"""Page-level text collision check: word boxes that overlap (> 0.3 pt in both axes) in a candidate but not in its source."""
import sys, json, itertools
from pathlib import Path
import pdfplumber
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
def boxes(pdf):
    pg = pdfplumber.open(pdf).pages[0]
    ws = pg.extract_words(keep_blank_chars=False, use_text_flow=False)
    return [(w['text'], w['x0'], w['top'], w['x1'], w['bottom']) for w in ws]
def hits(B):
    out = []
    for a, b in itertools.combinations(B, 2):
        ox = min(a[3], b[3])-max(a[1], b[1]); oy = min(a[4], b[4])-max(a[2], b[2])
        if ox > 0.3 and oy > 0.3: out.append((a[0], b[0], round(a[1], 1), round(a[2], 1)))
    return out
res = {}
for cand, src in (('Figure.3', 'Figure.3'), ('Figure.4', 'Figure.5'), ('Figure.5', 'Figure.4')):
    hc = hits(boxes(HERE/'out'/f'{cand}.pdf')); hs = hits(boxes(ROOT/'Fig'/f'{src}.pdf'))
    pairs_s = {(a, b) for a, b, *_ in hs}
    new = [h for h in hc if (h[0], h[1]) not in pairs_s]
    res[cand] = dict(candidate=len(hc), source=len(hs), new=new)
    print(cand, 'overlaps cand/src', len(hc), len(hs), 'new:', new[:10])
(HERE/'out/overlap-check.json').write_text(json.dumps(res, indent=1, ensure_ascii=False))
