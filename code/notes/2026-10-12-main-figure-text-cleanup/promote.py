"""Promote out/Figure.1–4.pdf to Fig/ and update Fig/figure-manifest.tsv (source_candidate, sha256).
Pre-change copies: revision-backups/20261012-before-main-figure-text-cleanup/."""
import csv, hashlib, json, shutil
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
man = ROOT/'Fig/figure-manifest.tsv'
rows = list(csv.reader(open(man), delimiter='\t'))
hdr = rows[0]; i_file, i_src, i_sha = hdr.index('file'), hdr.index('source_candidate'), hdr.index('sha256')
log = []
for n in (1, 2, 3, 4):
    src = HERE/f'out/Figure.{n}.pdf'; dst = ROOT/f'Fig/Figure.{n}.pdf'
    old = sha(dst); shutil.copyfile(src, dst); new = sha(dst); assert new == sha(src)
    for r in rows[1:]:
        if r[i_file] == f'Figure.{n}.pdf':
            r[i_src] = str(src.relative_to(ROOT)); r[i_sha] = new
    log.append(dict(file=f'Fig/Figure.{n}.pdf', old_sha=old, new_sha=new))
with open(man, 'w', newline='') as f:
    csv.writer(f, delimiter='\t', lineterminator='\n').writerows(rows)
json.dump(log, open(HERE/'promote.log.json', 'w'), indent=1)
for e in log: print(e['file'], e['old_sha'][:12], '->', e['new_sha'][:12])
