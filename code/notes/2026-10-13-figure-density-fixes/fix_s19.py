"""S19e: the six rows (Hs/Mm proximal rank, rank difference, log span) sat at a 7.0 pt pitch with 7.0 pt labels
(zero clearance). The rows are spread to an 8.5 pt pitch: row i moves down 1.5*i pt with its label, point, q value
and band shading; the dotted zero line is redrawn to the new span; the axis, ticks and axis title of e move down
7.5 pt; everything from the f/g row down moves 7.5 pt; the page gains 7.5 pt at the bottom (504.9 -> 512.4 pt).
Plotted x positions are unchanged."""
import json
from fixlib import Job, ROOT, OUT, enlarge_bottom
J = Job(19, src=ROOT/'Fig/SuppFig.S19.pdf', name='SuppFig.S19')
TOP0, PITCH, STEP, X0, X1 = 221.6, 7.0, 1.5, 296.5, 453.54
DROP = STEP*5                                   # 7.5 pt
J.remove((401.0, 221.3, 401.85, 264.9), 'e dotted zero line (redrawn over the new span)')
shifts = [((0.0, 289.9, X1, 504.9), DROP, 'f–i rows'), ((X0, 264.45, X1, 289.85), DROP, 'e axis, ticks and title')]
for i in range(5, 0, -1):                       # bottom-up: a moved row must not land in a region cleared later
    top = TOP0+PITCH*i
    shifts.append(((X0, top-0.02, X1, top+PITCH+0.02), STEP*i, f'e row {i}'))


def draw(t):
    lw = 0.4; on, off = 0.4, 0.6
    t.line([401.42, 401.42], [TOP0, 264.6+DROP], lw=lw, color='#87969E', dashes=(on/lw, off/lw))
    # the e x-axis baseline is cleared with its band but does not survive the block copy (the ticks do); redraw it
    t.line([371.8, 431.0], [264.6+DROP, 264.6+DROP], lw=0.5, color='#87969E')


rec = J.run(draw=draw, shifts=shifts)
rec['page'] = dict(before=[0, 0, 453.54, 504.9], after=enlarge_bottom(OUT/'SuppFig.S19.pdf', DROP, 'SuppFig.S19'))
import hashlib, subprocess
rec['sha256'] = hashlib.sha256((OUT/'SuppFig.S19.pdf').read_bytes()).hexdigest()
subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(OUT/'SuppFig.S19.pdf'), str(OUT/'SuppFig.S19')], check=True, capture_output=True)
json.dump(rec, open('out/receipt.S19.json', 'w'), indent=1, ensure_ascii=False)
print('removed', rec['words_removed'], 'added', rec['words_added'], rec['page']['after'])
