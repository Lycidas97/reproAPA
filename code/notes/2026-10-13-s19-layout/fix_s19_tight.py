"""S19 tightened (2026-10-13), from the 10-12 production S19 (revision-backups/20261013-before-figure-density-fixes).
Supersedes notes/2026-10-13-figure-density-fixes/fix_s19.py, which spread e but made row 2 taller and left c with an
empty band below. Here (see plan.py):
- e moves up 3.2 pt (d-to-e gap back to the 11.3 pt panel gap) and its rows go from a 7.0 to an 8.0 pt pitch
  (1.0 pt between labels), so row 2 ends only 1.8 pt lower than before;
- g's four bars go from a 15 to a 12 pt pitch, its axis, title and key follow;
- h's four rows go from a 13.5 to an 11.5 pt pitch, its axis, title and key follow, and i's key moves with them so the
  two keys stay on one line;
- the page loses 13.2 pt (504.9 -> 491.7 pt).
x positions and all plotted values are unchanged; e's dotted zero line and x-axis baseline are redrawn to the new span."""
import json, hashlib, subprocess
from fixlib import Job, ROOT, OUT, enlarge_bottom
import plan
J = Job(19, src=ROOT/plan.SRC, name='SuppFig.S19')
J.exact += [((401.42, 221.6, 401.42, 264.6), 'e dotted zero line'), ((371.8, 264.6, 431.0, 264.6), 'e x-axis baseline')]
Y0 = 221.6+plan.U; YA = 264.6+plan.U+5*plan.EP


def draw(t):
    t.line([401.42, 401.42], [Y0, YA], lw=0.4, color='#87969E', dashes=(0.4/0.4, 0.6/0.4))
    t.line([371.8, 431.0], [YA, YA], lw=0.5, color='#87969E')


rec = J.run(draw=draw, shifts=plan.regions())
total = plan.S4-3*plan.HP
rec['page'] = dict(before=[0, 0, 453.54, 504.9], after=enlarge_bottom(OUT/'SuppFig.S19.pdf', total, 'SuppFig.S19'))
rec['sha256'] = hashlib.sha256((OUT/'SuppFig.S19.pdf').read_bytes()).hexdigest()
rec['plan'] = dict(U=plan.U, EP=plan.EP, S3=plan.S3, GP=plan.GP, S4=plan.S4, HP=plan.HP)
subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(OUT/'SuppFig.S19.pdf'), str(OUT/'SuppFig.S19')], check=True, capture_output=True)
json.dump(rec, open('out/receipt.S19.json', 'w'), indent=1, ensure_ascii=False)
print('exact', rec.get('exact_removed'), 'removed', rec['words_removed'], 'added', rec['words_added'], 'page', rec['page']['after'])
