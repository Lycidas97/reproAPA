"""Figure 4 in-place recolouring (2026-10-11 caption review); no geometry or text changes.
c  Trip12 call-out (label, leader line, ring) from the P-state pink to ink: pink means the P state elsewhere.
e  decomposition terms from blues (blue means proximal elsewhere) to neutral slates: family position #6B7B85,
   family mass #C5CDD2; outside-family tan unchanged.
f  'net gain' from blue to the orange used for 'Gained' in 4b (bars, legend square, per-gene dots); net loss stays dark.
"""
import json
from fixlib import Job, INK
J = Job(4)
rec = J.run(recolors=[
    ([(345.0, 60.0, 405.0, 132.0)], {'#CC79A7': INK}, '4c Trip12 call-out'),
    ([(0.0, 336.0, 136.0, 425.0)], {'#2D6A96': '#6B7B85', '#9CC3DD': '#C5CDD2'}, '4e decomposition terms'),
    ([(140.0, 336.0, 453.6, 425.0)], {'#0072B2': '#D55E00'}, '4f net gain'),
])
print(json.dumps({k: rec[k] for k in ('edits', 'words_removed', 'words_added', 'sha256')}, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig4.json', 'w'), indent=1, ensure_ascii=False)
