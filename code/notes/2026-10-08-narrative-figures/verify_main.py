"""Structural check of kept blocks: every word / rect / line / curve of the source region must reappear in the candidate,
shifted by the block offset (tol 0.05 pt), and nothing else may sit in the moved region except the redrawn letters."""
import json, sys
from pathlib import Path
import pdfplumber
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
def objs(pdf, region, dx=0.0, dy=0.0, ignore=()):
    pg = pdfplumber.open(pdf).pages[0]; x0, y0, x1, y1 = region; out = []
    for kind, L in (('char', pg.chars), ('rect', pg.rects), ('line', pg.lines), ('curve', pg.curves)):
        for o in L:
            cx, cy = (o['x0']+o['x1'])/2, (o['top']+o['bottom'])/2
            if not (x0 <= cx <= x1 and y0 <= cy <= y1): continue
            if o['x1']-o['x0'] > 400 and o['bottom']-o['top'] > 400: continue   # page-sized clip/background rects
            if any(a <= cx <= c and b <= cy <= d for a, b, c, d in ignore): continue
            out.append((kind, o.get('text', ''), round(o['x0']+dx, 1), round(o['top']+dy, 1), round(o['x1']+dx, 1), round(o['bottom']+dy, 1)))
    return sorted(out)
rc = json.loads((HERE/'out/receipt.json').read_text())
F = ROOT/'Fig'; O = HERE/'out'
checks = []
def chk(name, src, cand, region, dx, dy, ign_src=(), ign_cand=()):
    a = objs(src, region, dx, dy, ign_src)
    b = objs(cand, (region[0]+dx, region[1]+dy, region[2]+dx, region[3]+dy), 0, 0, ign_cand)
    sa, sb = set(a), set(b)
    miss = sorted(sa-sb); extra = sorted(sb-sa)
    checks.append(dict(check=name, source_objects=len(a), candidate_objects=len(b), missing=len(miss), extra=len(extra), miss_sample=miss[:5], extra_sample=extra[:5]))
L3 = (0, 0, 0, 0)
# Figure 3
r3 = rc['Figure.3']; dy = [b for b in r3['blocks'] if 'dy' in b][0]['dy']
chk('F3 a,b', F/'Figure.3.pdf', O/'Figure.3.pdf', (0, 0, 453.6, 137.0), 0, 0)
cl = [(305.0, 155.0, 324.5, 164.6), (305.0, 270.0, 324.5, 279.6)]
for b in r3['blocks']:
    if 'to' in b and b['region'][1] > 140:
        R = b['region']; chk(f'F3 scatter {R[1]}', F/'Figure.3.pdf', O/'Figure.3.pdf', (R[0], R[1], R[2], R[3]), b['to'][0]-R[0], b['to'][1]-R[1], ign_src=cl)
chk('F3 d,e', F/'Figure.3.pdf', O/'Figure.3.pdf', (0, 414.0, 453.6, 618.4), 0, dy)
# Figure 4 (old 5)
r4 = rc['Figure.4']; B = {b['panel']: b for b in r4['blocks']}
src5 = F/'Figure.5.pdf'
chk('F4 a (old 5a)', src5, O/'Figure.4.pdf', (0, 0, 326.0, 105.0), 0, 0)
chk('F4 c (old 5c)', src5, O/'Figure.4.pdf', (326.0, 0, 453.6, 197.0), 0, 0)
lf = (3.0, 384.0, 9.0, 398.0); lb = (3.0, 109.0, 12.0, 123.0)
dF = B['old f -> d']['dy']; bx, bd = B['old b -> f']['dx'], B['old b -> f']['dy']
chk('F4 d (old 5f)', src5, O/'Figure.4.pdf', (0, 380.0, 453.6, 523.38), 0, dF, ign_src=[lf], ign_cand=[(lf[0], lf[1]+dF, lf[2], lf[3]+dF)])
chk('F4 f (old 5b)', src5, O/'Figure.4.pdf', (0, 105.0, 326.0, 197.0), bx, bd, ign_src=[lb], ign_cand=[(lb[0]+bx, lb[1]+bd, lb[2]+bx, lb[3]+bd)])
# Figure 5 (old 4)
r5 = rc['Figure.5']; dy5 = r5['dy']
lets = [(3.0, 169.0, 12.0, 183.0), (239.0, 169.0, 249.0, 183.0), (2.5, 296.0, 11.0, 310.0)]
chk('F5 a-c (old 4c-e)', F/'Figure.4.pdf', O/'Figure.5.pdf', (0, 164.0, 453.6, 481.12), 0, dy5, ign_src=lets, ign_cand=[(a, b+dy5, c, d+dy5) for a, b, c, d in lets])
out = dict(checks=checks, all_pass=all(c['missing'] == 0 and c['extra'] == 0 for c in checks))
(HERE/'out/verify.json').write_text(json.dumps(out, indent=1, ensure_ascii=False))
for c in checks: print(c['check'], c['source_objects'], c['candidate_objects'], 'missing', c['missing'], 'extra', c['extra'], c['miss_sample'][:2], c['extra_sample'][:2])
print('ALL PASS' if out['all_pass'] else 'CHECK')
