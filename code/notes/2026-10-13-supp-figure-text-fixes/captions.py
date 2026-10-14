r"""Text edits accompanying the 2026-10-13 supplementary figure fixes: the S3f legend now says 'primary fit', so the
caption and the supplementary results use the same term. Each replacement must match exactly once; --apply writes."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
CAP = ROOT/'supplement/supplementary-figures/main.tex'
RES = ROOT/'supplement/supplementary-results/expression-states-pas-measurement-and-the-global-apa-direction.tex'
EDITS = [
 (CAP, 'route pairs: the original fit, its centered null 95\\% range', 'route pairs: the primary fit, its centered null 95\\% range', 'S3f caption 1'),
 (CAP, 'The null range belongs to the original fit and', 'The null range belongs to the primary fit and', 'S3f caption 2'),
 (RES, 'route pairs with the original fit and', 'route pairs with the primary fit and', 'results 1'),
 (RES, 'The original-fit exceptions', 'The primary-fit exceptions', 'results 2'),
]
texts = {p: p.read_text() for p in {e[0] for e in EDITS}}
for p, old, new, note in EDITS:
    n = texts[p].count(old); print(f'[{note}] {p.name}: matches={n}'); assert n == 1
    texts[p] = texts[p].replace(old, new)
if '--apply' in sys.argv:
    for p, t in texts.items(): p.write_text(t)
    print('written')
