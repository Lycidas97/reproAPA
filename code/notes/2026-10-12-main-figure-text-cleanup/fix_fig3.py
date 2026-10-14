"""Figure 3 (2026-10-12 main-figure text cleanup).
3d  the species suffixes ' (human)' and ' (mouse)' added on 2026-10-11 are removed: the gene symbols alone
    (PDLIM7 over the Xing2020 column, Pdlim7 centred over the two mouse columns) identify the example, and the
    caption names species and studies. 'Long' (transcript-row label) sat 0.2 pt above the y tick label '100'; it
    moves up 0.8 pt.
3e  x tick labels '+.2' and '+.4' stood 3.2 pt apart; '+.4' is dropped (tick mark stays), leaving −.2, 0, +.2."""
import json
from fixlib import Job, chars_of, textwidth, INK
J = Job(3)
g1 = chars_of(J.src, 'PDLIM7', (18, 295, 52, 306))
J.remove((19.5, 296.0, 195.0, 305.5), '3d gene titles with species suffixes (redrawn as gene symbols only)')
J.move_word('Long', (3, 322, 22, 331.5), -0.8, note='3d transcript-row label clear of y tick 100')
J.drop_word('+.4', (370, 305, 386, 315), note='3e x tick label (crowded with +.2)')
B23 = (100.2, 227.6)   # the two mouse columns (as in the 2026-10-11 script)


def draw(t):
    t.text(g1['origin'], g1['base'], 'PDLIM7', 8.0, font='BolIta', color=INK)
    t.text((B23[0]+B23[1])/2, g1['base'], 'Pdlim7', 8.0, font='BolIta', color=INK, ha='center')


rec = J.run(draw=draw)
print(json.dumps({k: rec[k] for k in ('removed_kinds', 'words_removed', 'words_added', 'sha256')}, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig3.json', 'w'), indent=1, ensure_ascii=False)
