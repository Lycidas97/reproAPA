"""Promote out/Figure.1.pdf, out/Figure.4.pdf, out/SuppFig.S19.pdf to Fig/ and update Fig/figure-manifest.tsv.
Pre-change copies: revision-backups/20261013-before-figure-density-fixes/."""
import csv, hashlib, json, shutil
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
man = ROOT/'Fig/figure-manifest.tsv'; rows = list(csv.reader(open(man), delimiter='\t')); hdr = rows[0]
i_file, i_src, i_sha = hdr.index('file'), hdr.index('source_candidate'), hdr.index('sha256')
log = []
for f in ('Figure.1.pdf', 'Figure.4.pdf', 'SuppFig.S19.pdf'):
    src, dst = HERE/'out'/f, ROOT/'Fig'/f
    old = sha(dst); shutil.copyfile(src, dst); new = sha(dst); assert new == sha(src)
    hit = [r for r in rows[1:] if r[i_file] == f]; assert len(hit) == 1
    hit[0][i_src] = str(src.relative_to(ROOT)); hit[0][i_sha] = new
    log.append(dict(file=f'Fig/{f}', old_sha=old, new_sha=new))
with open(man, 'w', newline='') as fh: csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(rows)
json.dump(log, open(HERE/'promote.log.json', 'w'), indent=1)
for e in log: print(e['file'], e['old_sha'][:12], '->', e['new_sha'][:12])
