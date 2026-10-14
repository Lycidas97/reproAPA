"""Copy the rebuilt S1–S21 from out/set into Fig/ and update Fig/figure-manifest.tsv (source_candidate, sha256).
Pre-change copies: revision-backups/20261012-before-supp-figure-text-cleanup/Fig/."""
import csv, hashlib, shutil, json
from pathlib import Path
RUN = Path(__file__).resolve().parent; ROOT = RUN.parents[1]
SET = RUN/'out'/'set'; FIG = ROOT/'Fig'; MAN = FIG/'figure-manifest.tsv'
rows = list(csv.DictReader(MAN.open(), delimiter='\t'))
log = []
for r in rows:
    if not r['figure'].startswith('Figure S'): continue
    n = r['figure'].split('S')[1]
    src = SET/f'SuppFig.S{n}.pdf'; dst = FIG/r['file']
    old = hashlib.sha256(dst.read_bytes()).hexdigest()
    shutil.copyfile(src, dst)
    new = hashlib.sha256(dst.read_bytes()).hexdigest()
    r['source_candidate'] = f'notes/2026-10-12-supp-figure-text-cleanup/out/set/SuppFig.S{n}.pdf'; r['sha256'] = new
    log.append(dict(figure=r['figure'], old_sha256=old, new_sha256=new))
with MAN.open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
(RUN/'promote.log.json').write_text(json.dumps(log, indent=1))
print(len(log), 'figures promoted')
