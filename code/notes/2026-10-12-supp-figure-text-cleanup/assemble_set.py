"""Collect the rebuilt S1–S21 under their v11 names into <run>/out/set (candidate only; Fig/ is written by promote.py)."""
import shutil, hashlib, subprocess, json
from pathlib import Path
import pikepdf
RUN = Path(__file__).resolve().parent; O = RUN/'out'; SET = O/'set'; SET.mkdir(parents=True, exist_ok=True)
SRC = {1: 'compose/A/SuppFig.S1', 2: 'compose/B/SuppFig.S2', 3: 'compose/C/SuppFig.S3', 4: 'compose/D/SuppFig.S4', 5: 'compose/E/SuppFig.S5',
       6: 'compose/J/SuppFig.S6', 7: 'compose/V7/SuppFig.S7', 8: 'compose/V8/SuppFig.S8', 9: 'compose/M/SuppFig.S9', 10: 'v11/SuppFig.S10',
       11: 'compose/N/SuppFig.S11', 12: 'v11/SuppFig.S12', 13: 'compose/P/SuppFig.S13', 14: 'compose/T/SuppFig.S14', 15: 'compose/S/SuppFig.S15',
       16: 'compose/Q/SuppFig.S16', 17: 'compose/R/SuppFig.S17', 18: 'v11/SuppFig.S18', 19: 'v11/SuppFig.S19', 20: 'v11/SuppFig.S20',
       21: 'a01/current/SuppFig.S31'}
rows = []
for n, s in SRC.items():
    src = O/(s+'.pdf'); dst = SET/f'SuppFig.S{n}.pdf'
    shutil.copyfile(src, dst)
    for ext in ('receipt.json', 'values.json'):
        if (O/(s+'.'+ext)).exists(): shutil.copyfile(O/(s+'.'+ext), SET/f'SuppFig.S{n}.{ext}')
    with pikepdf.open(dst) as p: w, h = (float(v) for v in p.pages[0].MediaBox[2:])
    subprocess.run(['pdftoppm', '-r', '110', '-png', '-singlefile', str(dst), str(SET/f'SuppFig.S{n}')], check=True, capture_output=True)
    rows.append(dict(figure=f'S{n}', source=f'notes/2026-10-12-supp-figure-text-cleanup/out/{s}.pdf', sha256=hashlib.sha256(dst.read_bytes()).hexdigest(),
                     width_pt=round(w, 2), height_pt=round(h, 2), height_mm=round(h/72*25.4, 1)))
(SET/'manifest.json').write_text(json.dumps(rows, indent=1))
for r in rows: print(r['figure'], r['height_mm'], r['width_pt'], r['sha256'][:12])
subprocess.run(['qpdf', '--empty', '--pages']+[str(SET/f'SuppFig.S{n}.pdf') for n in range(1, 22)]+['--', str(SET/'Supplement-S1-S21.pdf')], check=True)
