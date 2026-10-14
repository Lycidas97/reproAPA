"""Promote out/SuppFig.S19.pdf to Fig/ and update Fig/figure-manifest.tsv. Pre-change copy:
revision-backups/20261013-before-s19-tighten/ (the 10-13 03:00 version); the 10-12 version is in
revision-backups/20261013-before-figure-density-fixes/."""
import csv, hashlib, json, shutil
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
man = ROOT/'Fig/figure-manifest.tsv'; rows = list(csv.reader(open(man), delimiter='\t')); hdr = rows[0]
f = 'SuppFig.S19.pdf'; src, dst = HERE/'out'/f, ROOT/'Fig'/f
old = sha(dst); shutil.copyfile(src, dst); new = sha(dst); assert new == sha(src)
hit = [r for r in rows[1:] if r[hdr.index('file')] == f]; assert len(hit) == 1
hit[0][hdr.index('source_candidate')] = str(src.relative_to(ROOT)); hit[0][hdr.index('sha256')] = new
with open(man, 'w', newline='') as fh: csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(rows)
json.dump([dict(file='Fig/'+f, old_sha=old, new_sha=new)], open(HERE/'promote.log.json', 'w'), indent=1)
print(f, old[:12], '->', new[:12])
