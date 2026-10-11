"""Captions of the 21 consolidated figures (consolidated/captions/S<n>.en.txt) and their checks.

Unchanged figures (one source, all panels, same letters): the standalone caption (captions/S<old>.en.txt) with the new
number; cross-references to other supplementary figures are renumbered. Merged figures: hand-written captions already in
consolidated/captions/ (rewritten from the source captions for the new letters; text of removed panels dropped).
Checks (as tools/check_captions.py): the panel letters cited match the receipt; every number appears in values.json or
the receipt notes (within the printed rounding). Writes consolidated/captions/CAPTIONS.en.md and captions.json.
Usage: python3 tools/consolidated_captions.py [--build]   (--build regenerates the unchanged captions first)
"""
import sys, json, re, glob
from pathlib import Path
HERE = Path(__file__).resolve().parent; BASE = HERE.parent
sys.path.insert(0, str(HERE))
from compose import SPEC, PRINT, PNUM, OLD, cdir
CAP = BASE/'consolidated'/'captions'


def unchanged(cid):
    sp = SPEC[cid]
    if 'external' in sp: return sp['external']
    b = sp['bands']
    return b[0][0][0] if len(b) == 1 and len(b[0]) == 1 and not sp.get('patches') else None


# old printed number -> new printed number, only where the letters are unchanged
RENUM = {OLD[unchanged(c)]: PNUM[c] for c in PRINT if unchanged(c)}

if '--build' in sys.argv:
    for c in PRINT:
        f = unchanged(c)
        if not f: continue
        t = (BASE/'captions'/f'S{OLD[f]}.en.txt').read_text().strip()
        assert t.startswith(f'Figure S{OLD[f]}. ')
        body = t[len(f'Figure S{OLD[f]}. '):]
        def rn(m):
            n = int(m.group(1)); assert n in RENUM, (c, m.group(0)); return f'Figure S{RENUM[n]}'
        body = re.sub(r'Figure S(\d+)', rn, body)
        (CAP/f'S{PNUM[c]}.en.txt').write_text(f'Figure S{PNUM[c]}. {body}\n')


def nums(o, out):
    if isinstance(o, dict):
        for k, v in o.items(): nums(k, out); nums(v, out)
    elif isinstance(o, list):
        for v in o: nums(v, out)
    elif isinstance(o, bool): pass
    elif isinstance(o, (int, float)): out.add(o)
    elif isinstance(o, str):
        for t in re.findall(r'-?\d+(?:\.\d+)?', o): out.add(float(t))


md, js, report = [], {}, []
for c in PRINT:
    n = PNUM[c]
    cap = (CAP/f'S{n}.en.txt').read_text().strip()
    assert cap.startswith(f'Figure S{n}. '), n
    r = json.loads((cdir(c)/f'SuppFig.S{n}.receipt.json').read_text())
    pan = sorted(p.lower() for p in r['panels'])
    refs = set()
    body0 = re.sub(r'Figure S?\d+[a-z]?(?:, [a-z](?![a-z]))*', '', cap)          # cross-references are not panel citations
    for m in re.finditer(r'\(([a-z])(?:([,–-]) ?([a-z]))?\)', body0):
        a, sep, b = m.groups()
        if b and sep == '–': refs |= {chr(i) for i in range(ord(a), ord(b)+1)}
        else: refs |= {a} | ({b} if b else set())
    for m in re.finditer(r'Panels ([a-z])(?:–([a-z]))?', body0):
        a, b = m.groups(); refs |= {chr(i) for i in range(ord(a), ord(b or a)+1)}
    miss, extra = sorted(set(pan)-refs), sorted(refs-set(pan))
    vals = set(); nums(json.loads((cdir(c)/f'SuppFig.S{n}.values.json').read_text()), vals); nums(r.get('notes', []), vals)
    body = re.sub(r'Figure S\d+\.', '', cap); body = re.sub(r'Figure \d[a-z]?|Figure S\d+[a-z, ]*', '', body)
    un = []
    for t in re.findall(r'(?<![A-Za-z0-9])[−-]?\d[\d,]*(?:\.\d+)?|(?<![\d])\.\d+', body):
        s = t.replace(',', '').replace('−', '-')
        try: x = float(s)
        except ValueError: continue
        tol = 0.5*10**-(len(s.split('.')[1]) if '.' in s else 0)+1e-9
        if not any(abs(abs(v)-abs(x)) <= tol for v in vals if isinstance(v, (int, float))): un.append(t)
    title = cap.split('. ', 2)[1] if cap.count('. ') >= 2 else ''
    report.append(dict(S=n, plan=c, panels=pan, missing=miss, extra=extra, words=len(cap.split()), unmatched=list(dict.fromkeys(un))))
    md.append(f'## Figure S{n}\n\n{cap}\n'); js[f'S{n}'] = dict(plan_letter=c, caption=cap, words=len(cap.split()))
(CAP/'CAPTIONS.en.md').write_text('# Supplementary figure legends (consolidated, 21 figures)\n\n'+'\n'.join(md))
(CAP/'captions.json').write_text(json.dumps(js, indent=1, ensure_ascii=False))
(BASE/'verification'/'captions-consolidated.json').write_text(json.dumps(report, indent=1, ensure_ascii=False))
for x in report:
    print(f"S{x['S']:<3}{x['plan']} panels={''.join(x['panels'])} missing={x['missing']} extra={x['extra']} words={x['words']}"
          + (f"  unmatched: {' '.join(x['unmatched'])}" if x['unmatched'] else ''))
