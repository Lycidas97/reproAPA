"""Structural check of the S13 candidate blocks (same method as verify_main.py)."""
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
src = (HERE/'verify_main.py').read_text()
g = {'__file__': str(HERE/'verify_main.py')}
exec(src[:src.index('rc = json.loads')], g)
objs = g['objs']
rc = json.loads((HERE/'out/receipt.json').read_text())['SuppFig.S13']
F = ROOT/'Fig'; O = HERE/'out/SuppFig.S13.pdf'; out = []
import pdfplumber
S = F/'SuppFig.S13.pdf'
lets = [(c['x0']-0.7, c['top']-0.7, c['x1']+0.7, c['bottom']+0.7) for c in pdfplumber.open(S).pages[0].chars if c['size'] >= 9.5 and 'Bol' in c['fontname'] and len(c['text']) == 1]
for k in rc['kept']:
    R = k['region']
    a = set(objs(S, R, 0, 0, lets)); bb = set(objs(O, R, 0, 0, lets))
    out.append(dict(block=f"{k['old']} -> {k['new']}", source=len(a), cand=len(bb), missing=len(a-bb), extra=len(bb-a)))
    print(out[-1])
(HERE/'out/verify-supp.json').write_text(json.dumps(out, indent=1))
