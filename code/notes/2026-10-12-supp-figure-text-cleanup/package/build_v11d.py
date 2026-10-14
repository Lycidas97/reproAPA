"""Refresh the consolidated source-data archive v11c -> v11d (2026-10-12, round 6 of the supplementary-figure revision).
Same method as build_v11c.py: the v11c tar.gz is streamed member by member (no extraction); project files whose working-tree
bytes changed are replaced, the round-6 edit records, the independent revision check and the v11d packaging scripts are
added, and everything else (code/ snapshot, workbook, figures, frozen inputs) is copied unchanged. Identical payloads are
stored once as tar hard links. Writes the ledger, README, LICENSE and refresh-v11d.json into
zenodo/reproapa-final-data-20261012-v11d/.

Usage: python3 build_v11d.py      (uncompressed tar -> $V11D_TAR, default /tmp/reproAPA-source-data-v11d.tar)
       gzip -6 -n -c $V11D_TAR > zenodo/reproapa-final-data-20261012-v11d/files/reproAPA-source-data-v11d.tar.gz
       python3 verify_v11d.py"""
import csv, hashlib, io, json, os, tarfile, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
def v11c_dir(root):                               # v11c moves to zenodo/superseded/ once v11d is current
    for d in (root/'zenodo/reproapa-final-data-20261012-v11c', root/'zenodo/superseded/reproapa-final-data-20261012-v11c'):
        if d.is_dir(): return d
    raise SystemExit('v11c upload set not found')
V11C = v11c_dir(ROOT)
OLD = V11C/'files/reproAPA-source-data-v11c.tar.gz'
DST = ROOT/'zenodo/reproapa-final-data-20261012-v11d'; FILES = DST/'files'; FILES.mkdir(parents=True, exist_ok=True)
OLD_ROOT, NEW_ROOT = 'reproAPA-v11c', 'reproAPA-v11d'
TAR = Path(os.environ.get('V11D_TAR', '/tmp/reproAPA-source-data-v11d.tar'))
ORIGIN = 'project working tree 2026-10-12 (supplementary-figure revision, round 6)'
RUN = 'notes/2026-10-12-supp-figure-text-cleanup'
EXPECTED_UPDATES = {                                     # the round-6 changes; anything else changed is an error
    'project/notes/2026-10-10-supplement-v11/figures/profile.py',
    f'project/{RUN}/README.zh.md',
    f'project/{RUN}/out/set/SuppFig.S20.receipt.json',
    'project/reviews/2026-10-12-supp-figure-caption-audit/REPORT.zh.md',
    'project/supplement/supplementary-figures/main.tex'}


def sha(b): return hashlib.sha256(b).hexdigest()


# ---------------------------------------------------------------- old ledger and the update/add plan
old_rows = {}
with tarfile.open(OLD, 'r:gz') as t:
    for m in t:
        if m.name == f'{OLD_ROOT}/FILES.tsv':
            old_rows = {r['path']: r for r in csv.DictReader(io.TextIOWrapper(t.extractfile(m), 'utf-8'), delimiter='\t')}
            break
assert len(old_rows) == 2897, len(old_rows)
assert (V11C/'files/FILES.tsv').read_bytes().decode().count('\n') == 2898
update = {}
for p, r in old_rows.items():
    if not p.startswith('project/'): continue
    b = (ROOT/p[8:]).read_bytes()
    if sha(b) != r['sha256']: update[p] = b
assert set(update) == EXPECTED_UPDATES, sorted(set(update) ^ EXPECTED_UPDATES)
add_paths = sorted(str(p.relative_to(ROOT)) for p in (ROOT/RUN/'edits').glob('*') if p.suffix in ('.py', '.json')
                   and 'project/'+str(p.relative_to(ROOT)) not in old_rows)   # round-6 edit scripts and their receipts
add_paths += ['reviews/2026-10-12-supp-figure-revision-check/REPORT.zh.md',
              f'{RUN}/package/build_v11d.py', f'{RUN}/package/verify_v11d.py']
add = {}
for a in add_paths:
    key = 'project/'+a
    assert key not in old_rows, key
    add[key] = (ROOT/a).read_bytes()

README = (ROOT/RUN/'package/README.v11d.md').read_bytes()
LICENSE = (ROOT/'LICENSE').read_bytes()
assert sha(LICENSE) == old_rows['LICENSE']['sha256']

# ---------------------------------------------------------------- final ledger content (stored_as filled while writing)
final = {}
for p, r in old_rows.items():
    if p in update: final[p] = dict(path=p, bytes=len(update[p]), sha256=sha(update[p]), origin=ORIGIN, origin_member=p[8:])
    elif p == 'README.md': final[p] = dict(path=p, bytes=len(README), sha256=sha(README), origin='consolidated release', origin_member='README.md')
    else: final[p] = dict(path=p, bytes=int(r['bytes']), sha256=r['sha256'], origin=r['origin'], origin_member=r['origin_member'])
for p, b in add.items(): final[p] = dict(path=p, bytes=len(b), sha256=sha(b), origin=ORIGIN, origin_member=p[8:])

# ---------------------------------------------------------------- stream-write
stored = {}                       # sha -> member path (relative to root) holding the payload
old_payload = {}                  # old payload of replaced members that other members link to
linked_targets = {r['stored_as'] for r in old_rows.values() if r['stored_as'] != r['path']}
now = int(time.time())


def put(out, rel, data, src=None):
    h = sha(data)
    ti = tarfile.TarInfo(f'{NEW_ROOT}/{rel}')
    if src is not None: ti.mode, ti.mtime, ti.uid, ti.gid, ti.uname, ti.gname = src.mode, src.mtime, src.uid, src.gid, src.uname, src.gname
    else: ti.mode, ti.mtime = 0o644, now
    if h in stored:
        ti.type = tarfile.LNKTYPE; ti.linkname = f'{NEW_ROOT}/{stored[h]}'; ti.size = 0; out.addfile(ti)
    else:
        ti.size = len(data); out.addfile(ti, io.BytesIO(data)); stored[h] = rel
    final[rel]['stored_as'] = stored[h]


with tarfile.open(OLD, 'r:gz') as t, tarfile.open(TAR, 'w', format=tarfile.GNU_FORMAT) as out:
    for m in t:
        rel = m.name[len(OLD_ROOT)+1:]
        if rel == 'FILES.tsv': continue                    # written last, after stored_as is known
        if m.isfile():
            data = t.extractfile(m).read()
            if rel in update:
                if rel in linked_targets: old_payload[rel] = data
                put(out, rel, update[rel], m)
            elif rel == 'README.md': put(out, rel, README, m)
            else:
                assert sha(data) == old_rows[rel]['sha256'], rel
                put(out, rel, data, m)
        elif m.islnk():
            target = m.linkname[len(OLD_ROOT)+1:]
            if rel in update: put(out, rel, update[rel], m)
            else:
                h = old_rows[rel]['sha256']
                if h in stored:
                    ti = tarfile.TarInfo(f'{NEW_ROOT}/{rel}'); ti.mode, ti.mtime, ti.uid, ti.gid, ti.uname, ti.gname = m.mode, m.mtime, m.uid, m.gid, m.uname, m.gname
                    ti.type = tarfile.LNKTYPE; ti.linkname = f'{NEW_ROOT}/{stored[h]}'; out.addfile(ti); final[rel]['stored_as'] = stored[h]
                else:                                          # its payload holder was replaced: store the old payload here
                    data = old_payload[target]; assert sha(data) == h, rel
                    put(out, rel, data, m)
        else:
            raise SystemExit(f'unexpected member type {m.name}')
    for rel in sorted(add): put(out, rel, add[rel])
    cols = ['path', 'bytes', 'sha256', 'stored_as', 'origin', 'origin_member']
    buf = io.StringIO(); w = csv.DictWriter(buf, fieldnames=cols, delimiter='\t', lineterminator='\n'); w.writeheader()
    for p in sorted(final): w.writerow({c: final[p][c] for c in cols})
    ledger = buf.getvalue().encode()
    ti = tarfile.TarInfo(f'{NEW_ROOT}/FILES.tsv'); ti.size, ti.mode, ti.mtime = len(ledger), 0o644, now; out.addfile(ti, io.BytesIO(ledger))
(FILES/'FILES.tsv').write_bytes(ledger)
(FILES/'README.md').write_bytes(README); (FILES/'LICENSE').write_bytes(LICENSE)
json.dump(dict(base='v11c', base_archive_sha256=hashlib.sha256(OLD.read_bytes()).hexdigest(), updated=sorted(update), added=sorted(add),
               relinked_from_replaced_holders=sorted(r for r in final if final[r].get('stored_as') and final[r]['stored_as'] != r
                                                     and old_rows.get(r, {}).get('stored_as') in update),
               members=len(final), unique_payloads=len(stored)), open(DST/'refresh-v11d.json', 'w'), indent=1)
print('members', len(final), 'payloads', len(stored), 'updated', len(update), 'added', len(add))
