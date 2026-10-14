"""Stream-verify reproAPA-source-data-v11d.tar.gz against its own FILES.tsv, against v11c and against the current working
tree. Writes zenodo/reproapa-final-data-20261012-v11d/stream-verification.json."""
import csv, hashlib, io, json, tarfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
DST = ROOT/'zenodo/reproapa-final-data-20261012-v11d'; ARCH = DST/'files/reproAPA-source-data-v11d.tar.gz'
def v11c_dir(root):                               # v11c moves to zenodo/superseded/ once v11d is current
    for d in (root/'zenodo/reproapa-final-data-20261012-v11c', root/'zenodo/superseded/reproapa-final-data-20261012-v11c'):
        if d.is_dir(): return d
    raise SystemExit('v11c upload set not found')
OLDF = v11c_dir(ROOT)/'files/FILES.tsv'
R = 'reproAPA-v11d/'
seen = {}; links = {}; types = {'file': 0, 'link': 0}; ledger = None
with tarfile.open(ARCH, 'r:gz') as t:               # full read also checks the gzip CRC
    for m in t:
        rel = m.name[len(R):]
        assert m.name.startswith(R), m.name
        if m.isfile():
            b = t.extractfile(m).read(); types['file'] += 1
            seen[rel] = (len(b), hashlib.sha256(b).hexdigest())
            if rel == 'FILES.tsv': ledger = b
        elif m.islnk():
            tgt = m.linkname[len(R):]; assert tgt in seen, (rel, tgt)
            seen[rel] = seen[tgt]; links[rel] = tgt; types['link'] += 1
        else: raise SystemExit(m.name)
rows = {r['path']: r for r in csv.DictReader(io.StringIO(ledger.decode()), delimiter='\t')}
assert (DST/'files/FILES.tsv').read_bytes() == ledger
probs = []
for p, r in rows.items():
    if p not in seen: probs.append(('missing', p)); continue
    if seen[p] != (int(r['bytes']), r['sha256']): probs.append(('hash', p))
    if r['stored_as'] != links.get(p, p): probs.append(('stored_as', p))
extra = sorted(set(seen) - set(rows) - {'FILES.tsv'})
old = {r['path']: r for r in csv.DictReader(open(OLDF), delimiter='\t')}
refresh = json.load(open(DST/'refresh-v11d.json'))
changed_vs_old = sorted(p for p in old if p in rows and old[p]['sha256'] != rows[p]['sha256'])
# every project/ member must equal the current working-tree file
stale = [p for p, r in rows.items() if p.startswith('project/')
         and hashlib.sha256((ROOT/p[8:]).read_bytes()).hexdigest() != r['sha256']]
fig = {}
for f in sorted((ROOT/'Fig').glob('*.pdf')):
    k = 'project/Fig/'+f.name
    fig[f.name] = k in rows and rows[k]['sha256'] == hashlib.sha256(f.read_bytes()).hexdigest()
man = rows.get('project/Fig/figure-manifest.tsv')
wb = 'project/table/Supplementary_Tables_S1-S10.xlsx'
res = dict(archive=ARCH.name, members=len(seen), files=types['file'], links=types['link'], ledger_rows=len(rows),
           problems=probs, members_not_in_ledger=extra,
           changed_vs_v11c_equals_planned_updates=(set(changed_vs_old) == set(refresh['updated']) | {'README.md'}),
           changed_vs_v11c=len(changed_vs_old), removed_vs_v11c=sorted(set(old) - set(rows)), added_vs_v11c=sorted(set(rows) - set(old)),
           project_members_stale_vs_working_tree=stale,
           current_figures_hash_match=all(fig.values()) and len(fig) == 25, figures_checked=len(fig),
           figure_manifest_match=man is not None and man['sha256'] == hashlib.sha256((ROOT/'Fig/figure-manifest.tsv').read_bytes()).hexdigest(),
           workbook_unchanged=rows[wb]['sha256'] == old[wb]['sha256'],
           figures_unchanged=all(rows[p]['sha256'] == old[p]['sha256'] for p in old if p.startswith('project/Fig/')),
           code_snapshot_unchanged=all(rows[p]['sha256'] == old[p]['sha256'] for p in old if p.startswith('code/')))
res['status'] = 'PASS' if (not probs and not extra and not stale and res['changed_vs_v11c_equals_planned_updates']
                           and res['current_figures_hash_match'] and res['figure_manifest_match'] and res['workbook_unchanged']
                           and res['figures_unchanged'] and res['code_snapshot_unchanged'] and not res['removed_vs_v11c']) else 'FAIL'
(DST/'stream-verification.json').write_text(json.dumps(res, indent=1)+'\n')
print(json.dumps({k: v for k, v in res.items() if k not in ('problems', 'added_vs_v11c')}, indent=1)[:2500]); print('problems', probs[:10])
