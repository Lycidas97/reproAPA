"""Figure 2 (2026-10-12 main-figure text cleanup).
a  the key row moves down 2.5 pt (it stood 0.7 pt under the x tick labels).
b  lower block: the subheading 'Same-state halves, median recurrent/observed' (added 2026-10-11) is removed and its
   row is renamed 'Same-state halves'; the subheading 'Wells, mouse P–M (DWUI q < .05)' becomes 'Well-level q < .05,
   mouse P–M' ('Wells' was an in-house term; the DWUI definition is in the caption). The six remaining lines are
   re-spaced at an even 11.36 pt pitch (was 9.1 pt for seven lines); the last line aligns with the moved 2a key.
c  title 'Strict tandem 3′ UTR, recurrent P–S genes' removed (the caption opens with this definition).
d  q column: eight values in 58 pt overlapped each other and the 'q' head. Each test now shows its q at the 0.02
   threshold (the value the text cites), centred on its row; the four q at 0.05 move to the caption.
e  the key row moves down 2 pt (it stood 0.7 pt under the x tick labels); the page gains 2 pt at the bottom
   (470.4 -> 472.4 pt) so that the bottom margin stays as before."""
import json
from fixlib import Job, chars_of, enlarge_bottom, INK, MUTED, OUT
J = Job(2)
# ---- b lower block
wells = chars_of(J.src, 'Wells,', (248, 84, 275, 93))
J.remove((250.6, 84.4, 366.0, 92.5), '2b wells subheading (redrawn shorter)')
J.remove((250.6, 111.2, 407.2, 119.4), '2b same-state halves subheading')
halves = chars_of(J.src, 'halves', (285, 119.5, 313, 129.5))
v10 = chars_of(J.src, '10/424,', (318, 119.5, 345, 129.5)); v45 = chars_of(J.src, '45/515,', (385, 119.5, 412, 129.5))
J.remove((261.0, 120.3, 437.2, 128.5), '2b same-state halves row (redrawn lower with new label)')
TOP = dict(boot=75.6, wells=84.7, all=93.3, rec=102.4, halves=120.6)
NEW = dict(boot=77.66, wells=89.02, all=100.38, rec=111.74, halves=123.1)
d = {k: round(NEW[k]-TOP[k], 2) for k in TOP}
# ---- c title
J.remove((22.3, 144.2, 165.1, 152.2), '2c title')
# ---- d q column: keep the 0.02 q (upper of each pair), centred on the row; drop the 0.05 q
# (text, top) of the q at 0.02 (upper of each pair, kept) and at 0.05 (lower, dropped); row centre = marker midpoint
ROWS = [(('.0002', 246.4), ('.0002', 254.2), 253.95),   # human identity
        (('.0008', 260.9), ('.022', 268.8), 268.4),     # human direction
        (('.0002', 279.4), ('.0002', 287.2), 286.95),   # mouse identity
        (('.0027', 293.9), ('.0008', 301.8), 301.4)]    # mouse direction
for (w02, t02), (w05, t05), ctr in ROWS:
    J.drop_word(w05, (429, t05+2, 451, t05+6), note='2d q at 0.05 (moved to caption)')
    J.move_word(w02, (429, t02+2, 451, t02+6), round(ctr-4.0-t02, 2), align='right', note='2d q at 0.02, centred on its row')

def draw(t):
    t.text(wells['origin'], wells['base']+d['wells'], 'Well-level q < .05, mouse P–M', wells['size'], color=wells['color'])
    t.text(halves['x1'], halves['base']+d['halves'], 'Same-state halves', halves['size'], color=halves['color'], ha='right')
    t.text(v10['origin'], v10['base']+d['halves'], '10/424, 31/404', v10['size'], color=v10['color'])
    t.text(v45['origin'], v45['base']+d['halves'], '45/515, 35/513', v45['size'], color=v45['color'])


# Shifts run bottom-up: each shift clears its source region on the page (including forms placed by earlier shifts),
# so a block moved down must not land in a region that is cleared later.
rec = J.run(draw=draw, shifts=[
    ((0.0, 458.25, 453.5, 470.4), 2.0, '2e key row'),
    ((0.0, 120.3, 233.0, 129.6), 2.5, '2a key row'),
    ((236.0, 101.75, 453.5, 110.8), d['rec'], '2b recurrent-members row'),
    ((236.0, 92.75, 453.5, 101.6), d['all'], '2b all-tested row'),
    ((236.0, 74.6, 453.5, 84.3), d['boot'], '2b bootstrap row'),
])
rec['row_tops'] = dict(old=TOP, new=NEW)
rec['page'] = dict(before=[0, 0, 453.5433, 470.4], after=enlarge_bottom(OUT/'Figure.2.pdf', 2.0, 'Figure.2'))
import hashlib, subprocess
rec['sha256'] = hashlib.sha256((OUT/'Figure.2.pdf').read_bytes()).hexdigest()
subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(OUT/'Figure.2.pdf'), str(OUT/'Figure.2')], check=True, capture_output=True)
print(json.dumps({k: rec[k] for k in ('removed_kinds', 'words_removed', 'words_added', 'sha256', 'page')}, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig2.json', 'w'), indent=1, ensure_ascii=False)
