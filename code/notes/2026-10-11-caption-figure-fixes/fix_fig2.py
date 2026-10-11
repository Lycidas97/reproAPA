"""Figure 2 in-place edits (2026-10-11 caption review).
a  count-column header 'Obs./mean' -> 'n/exp.' (same term as 2d and Figure 3).
b  upper axis gets its title 'Obs./null mean' on the tick-label line (left of '1'); the Human/Mouse key moves into the
   column headers of the lower block ('● Human', '○ Mouse'); the split-half rows get a muted subheading like the wells
   rows and become one row 'S and P halves' with 'recurrent/observed' values (no spaces around '/').
d  title '3′UTR' -> '3′ UTR'.
e  legend loses the q entries (every plotted test has q < .05; values are printed at right).
"""
import json
from fixlib import Job, chars_of, textwidth, INK, MUTED
J = Job(2)
# a
J.replace('Obs./mean', (130, 6, 175, 15), 'n/exp.', align='right')
# b: key row out, inline axis title in
J.remove((249.5, 56.0, 316.0, 64.5), '2b Human/Mouse key row (moved into lower-block column headers)')
tick1 = chars_of(J.src, '1', (335, 55, 343, 65))
hum = chars_of(J.src, 'Human', (330, 65, 360, 75)); mou = chars_of(J.src, 'Mouse', (398, 65, 426, 75))
# b: split-half rows
for w, box in [('S', (241, 110, 248, 120)), ('P', (241, 119.5, 248, 129))]:
    pass
J.remove((240.5, 110.8, 453.0, 128.8), '2b same-state halves rows (redrawn as subheading + one row)')
wells = chars_of(J.src, 'Wells,', (248, 84, 275, 93))
s_lab = chars_of(J.src, 'state', (294, 110, 313, 120))
v10 = chars_of(J.src, '10', (318, 110, 330, 120)); v45 = chars_of(J.src, '45', (385, 110, 397, 120))
p_row = chars_of(J.src, 'state', (294, 119.5, 313, 129))
# d
title = chars_of(J.src, 'Strict', (20, 272, 45, 282))
J.remove((21.5, 273.0, 165.0, 281.6), '2d title (redrawn with 3′ UTR)')
# e
J.remove((249.5, 457.5, 323.0, 466.5), '2e legend q entries')


def draw(t):
    size = 7.5
    t.text(tick1['x0']-6.0, tick1['base'], 'Obs./null mean', size, color=INK, ha='right')
    for c, filled in [(hum, True), (mou, False)]:
        t.marker(c['x0']-4.05, c['top']+3.8, 'o', 3.1, fc=INK if filled else 'white', ec=INK, mew=0.7)
    t.text(wells['origin'], s_lab['base'], 'Same-state halves, median recurrent/observed', size, color=MUTED)
    t.text(p_row['x1'], p_row['base'], 'S and P halves', size, color=INK, ha='right')
    t.text(v10['origin'], p_row['base'], '10/424, 31/404', size, color=INK)
    t.text(v45['origin'], p_row['base'], '45/515, 35/513', size, color=INK)
    t.text(title['origin'], title['base'], 'Strict tandem 3′ UTR, recurrent P–S genes', size, color=INK)


rec = J.run(draw=draw)
print(json.dumps({k: rec[k] for k in ('removed_kinds', 'words_removed', 'words_added', 'sha256')}, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig2.json', 'w'), indent=1, ensure_ascii=False)

# ---- stage 2: panel order follows first citation in the text.
# Old rows: R1 a|b (ink 4.5–127.75), R2 c (141.0–258.75), R3 d|e (270.25–465.75); page 470.4 pt.
# New rows: R1 unchanged, old R3 moved up to ink top 141.0 (same 13.25 pt visible gap), old R2 below it with the same
# 11.5 pt gap; page height unchanged. Letters then become c = strict tandem 3′ UTR, d = APA–GEX, e = composition.
import shutil
from fixlib import reorder_rows, OUT, Job as _Job
stage1 = OUT/'stage'; stage1.mkdir(exist_ok=True)
shutil.copy(OUT/'Figure.2.pdf', stage1/'Figure.2.edits.pdf')
W = 453.5433
placed = reorder_rows(stage1/'Figure.2.edits.pdf', stage1/'Figure.2.reordered.pdf', [
    ((0, 0, W, 134.5), 0.0),
    ((0, 264.0, W, 470.4), -129.25),
    ((0, 134.5, W, 263.0), 207.0),
])
K = _Job(2, src=stage1/'Figure.2.reordered.pdf', name='Figure.2')
K.replace('d', (3, 139, 14, 154), 'c', note='strict tandem 3′ UTR, was d')
K.replace('e', (234, 139, 242, 154), 'd', note='APA–GEX, was e')
K.replace('c', (3, 346, 14, 361), 'e', note='composition, was c')
rec2 = K.run()
rec2['row_moves'] = dict(placed=placed, note='R1 0–134.5 kept; R3 264.0–470.4 dy −129.25; R2 134.5–263.0 dy +207.0')
print(json.dumps({k: rec2[k] for k in ('removed_kinds', 'words_removed', 'words_added', 'sha256')}, indent=1, ensure_ascii=False))
json.dump(dict(stage1=rec, stage2=rec2), open('out/receipt.fig2.json', 'w'), indent=1, ensure_ascii=False)
