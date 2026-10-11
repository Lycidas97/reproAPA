"""Figure 3 in-place edits (2026-10-11 caption review).
Panel letters d and e swap so that they follow the order of first citation in the text (PDLIM7 is discussed right
after 3c, the homologous-PAS correspondence afterwards): left-bottom PDLIM7 panel = d, right correspondence panel = e.
The PDLIM7 title now names the species of each column: 'PDLIM7 (human)' over the Xing2020 column and
'Pdlim7 (mouse)' centred over the two mouse columns.
"""
import json
from fixlib import Job, chars_of, textwidth, INK, MUTED
J = Job(3)
J.replace('d', (238, 168, 250, 183), 'e', note='correspondence panel, was d')
J.replace('e', (3, 295, 13, 310), 'd', note='PDLIM7 panel, was e')
title = chars_of(J.src, 'PDLIM7/Pdlim7', (18, 295, 80, 306))
J.remove((19.5, 296.0, 80.0, 305.5), '3e title (redrawn with species)')
B1 = (32.5, 92.2); B23 = (100.2, 227.6)


def draw(t):
    g1, s1 = 'PDLIM7', ' (human)'
    t.text(title['origin'], title['base'], g1, 8.0, font='BolIta', color=INK)
    t.text(title['origin']+textwidth(g1, 8.0, 'BolIta'), title['base'], s1, 7.5, color=MUTED)
    g2, s2 = 'Pdlim7', ' (mouse)'
    wtot = textwidth(g2, 8.0, 'BolIta')+textwidth(s2, 7.5)
    x = (B23[0]+B23[1])/2-wtot/2
    t.text(x, title['base'], g2, 8.0, font='BolIta', color=INK)
    t.text(x+textwidth(g2, 8.0, 'BolIta'), title['base'], s2, 7.5, color=MUTED)


rec = J.run(draw=draw)
print(json.dumps({k: rec[k] for k in ('removed_kinds', 'words_removed', 'words_added', 'sha256')}, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig3.json', 'w'), indent=1, ensure_ascii=False)
