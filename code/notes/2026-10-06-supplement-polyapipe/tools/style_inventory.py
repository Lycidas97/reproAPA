"""Design-language inventory of the consolidated figures: every type size, line width, marker size/edge, alpha, dash and
colour actually drawn in the kept panels, classified against the house tokens (STYLE.zh.md rules 6, 7, 9).
Runs the source scripts through compose.run_source (same patches) and reads the recorded calls; nothing is rendered.
Usage: python3 tools/style_inventory.py [plan letters]   -> verification/style-consolidated.json
"""
import sys, json
from collections import Counter, defaultdict
from pathlib import Path
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import numpy as np
from matplotlib.colors import to_hex, to_rgba
import compose
from compose import SPEC, PRINT, PNUM, run_source
import supdraw as S

NAMED = {}
def add(name, c): NAMED.setdefault(to_hex(c).upper(), name)
for k, v in dict(INK=S.INK, MUTED=S.MUTED, AXC=S.AXC, NULLBAR=S.NULLBAR, NULLW=S.NULLW, MISSING=S.MISSING, BAND=S.BAND,
                 PROX=S.PROX, DIST=S.DIST, white='white').items(): add(k, v)
for k, v in S.SC.items(): add('state '+k, v)
for k, (n, c, m) in S.STUDY.items(): add('study '+n, c)
for k, (n, c) in S.DATASET.items(): add('dataset '+n, c)
for i, c in enumerate(S.SEQ_STOPS): add(f'seq stop {i}', c)
for i, c in enumerate(S.DIV): add(f'div stop {i}', c)
SEQ = [to_rgba(S.seq(t)) for t in np.linspace(0, 1, 501)]
DIVR = [to_rgba(S.DIV_CMAP(t)) for t in np.linspace(0, 1, 501)]


def classify(c):
    if c is None or (isinstance(c, str) and c == 'none'): return 'none'
    try: rgba = to_rgba(c)
    except Exception: return f'?{c}'
    h = to_hex(rgba).upper()
    if h in NAMED: return NAMED[h]
    d = lambda ramp: min(sum((a-b)**2 for a, b in zip(rgba[:3], r[:3])) for r in ramp)
    if d(SEQ) < 1e-5: return 'seq ramp (continuous)'
    if d(DIVR) < 1e-5: return 'diverging ramp (continuous)'
    return 'OTHER '+h


def inventory(cid):
    sp = SPEC[cid]
    if 'external' in sp: return None
    fids = list(dict.fromkeys(f for band in sp['bands'] for f, _ in band))
    keep = {f: set(L for band in sp['bands'] for g, Ls in band if g == f for L in Ls) for f in fids}
    inv = dict(size=Counter(), lw=Counter(), ms=Counter(), mew=Counter(), alpha=Counter(), dash=Counter(), colour=Counter(),
               raw=Counter())
    other = defaultdict(set)
    for f in fids:
        R, _ = run_source(f, sp.get('patches', {}).get(f, []))
        for L in keep[f]:
            for name, a, k in R.calls[L]:
                if name in ('text', 'vtext', 'vlabel'):
                    size = a[3]; col = a[5] if name == 'text' else (a[4] if name == 'vtext' else a[5])
                    inv['size'][size] += 1; c = classify(col); inv['colour'][('text', c)] += 1
                    if c.startswith('OTHER'): other[c].add(f'{f}{L} text')
                elif name == 'line':
                    lw, col, ls, alpha = a[2], a[3], a[4], a[5]
                    inv['lw'][lw] += 1; inv['alpha'][alpha] += 1; inv['dash'][ls if isinstance(ls, (str, type(None))) else 'RAW '+str(ls)] += 1
                    c = classify(col); inv['colour'][('line', c)] += 1
                    if c.startswith('OTHER'): other[c].add(f'{f}{L} line')
                    if k: inv['raw'][('line kw', tuple(sorted(k)))] += 1
                elif name == 'rect':
                    fc, ec, lw = a[4], a[5], a[6]
                    if lw: inv['lw'][lw] += 1
                    for role, c_ in (('fill', fc), ('edge', ec)):
                        c = classify(c_); inv['colour'][('rect '+role, c)] += 1
                        if c.startswith('OTHER'): other[c].add(f'{f}{L} rect {role}')
                    if k: inv['raw'][('rect kw', tuple(sorted(k.items())))] += 1
                elif name == 'marker':
                    ms, fc, ec, mew = a[3], a[4], a[5], a[6]
                    inv['ms'][ms] += 1; inv['mew'][mew if mew is not None else 'auto'] += 1
                    for role, c_ in (('fill', fc), ('edge', ec if ec is not None else fc)):
                        c = classify(c_); inv['colour'][('marker '+role, c)] += 1
                        if c.startswith('OTHER'): other[c].add(f'{f}{L} marker {role}')
                elif name == 'points':
                    ms, col, alpha, mew, fc = a[2], a[3], a[4], a[7], a[8]
                    inv['ms'][ms] += 1; inv['alpha'][alpha] += 1; inv['raw'][('points mew', mew)] += 1
                    c = classify(col if fc is None else fc); inv['colour'][('points', c)] += 1
                    if c.startswith('OTHER'): other[c].add(f'{f}{L} points')
                elif name == 'dot':
                    inv['raw'][('dot (area-encoded)', round(a[4], 2))] += 1
                    c = classify(a[3]); inv['colour'][('dot', c)] += 1
                    if c.startswith('OTHER'): other[c].add(f'{f}{L} dot')
                elif name == 'poly':
                    fc, ec, lw, alpha = a[1], a[2], a[3], a[4]
                    inv['raw'][('poly lw/alpha', lw, alpha)] += 1
                    c = classify(fc); inv['colour'][('poly', c)] += 1
                    if c.startswith('OTHER'): other[c].add(f'{f}{L} poly')
    js = {k: {str(kk): v for kk, v in c.items()} for k, c in inv.items()}
    js['other_colours'] = {c: sorted(w) for c, w in other.items()}
    return js


if __name__ == '__main__':
    out = {}
    for cid in (sys.argv[1:] or PRINT):
        r = inventory(cid)
        if r is None: continue
        out[f'S{PNUM[cid]} ({cid})'] = r
        print(f'S{PNUM[cid]} ({cid}) sizes {sorted(r["size"])} lw {sorted(r["lw"])} ms {sorted(r["ms"])} mew {sorted(r["mew"])} '
              f'alpha {sorted(r["alpha"])} dash {sorted(r["dash"])}')
        print('   raw', r['raw'])
        if r['other_colours']: print('   OTHER', {c: w[:4] for c, w in r['other_colours'].items()})
    p = HERE.parent/'verification'/('style-consolidated.json' if len(sys.argv) == 1 else 'style-consolidated-partial.json')
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False))
