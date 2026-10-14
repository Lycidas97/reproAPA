"""Refresh the consolidated source-data archive v11b -> v11c (2026-10-12 supplementary-figure revision).
Streams the v11b tar.gz member by member (no extraction): project files whose working-tree bytes changed are replaced,
the 12 Oct revision entry points and records are added, everything else (code/ snapshot, workbook, frozen inputs) is copied
unchanged. Identical payloads are stored once as tar hard links, as in v11b. Writes into zenodo/reproapa-final-data-20261012-v11c/.
Usage: python3 build_v11c.py; then gzip -6 -n -c $V11C_TAR > zenodo/reproapa-final-data-20261012-v11c/files/reproAPA-source-data-v11c.tar.gz"""
import csv, hashlib, io, json, os, tarfile, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT/'zenodo/reproapa-final-data-20261011-v11b/files/reproAPA-source-data-v11b.tar.gz'
DST = ROOT/'zenodo/reproapa-final-data-20261012-v11c'; FILES = DST/'files'; FILES.mkdir(parents=True, exist_ok=True)
OLD_ROOT, NEW_ROOT = 'reproAPA-v11b', 'reproAPA-v11c'
ARCH = FILES/'reproAPA-source-data-v11c.tar.gz'
# the uncompressed tar is written to V11C_TAR (default /tmp) and gzip-compressed into ARCH by a separate step (gzip -6)
TAR = Path(os.environ.get('V11C_TAR', '/tmp/reproAPA-source-data-v11c.tar'))
ORIGIN = 'project working tree 2026-10-12 (supplementary-figure revision)'
RUN = 'notes/2026-10-12-supp-figure-text-cleanup'

def sha(b): return hashlib.sha256(b).hexdigest()

# ---------------------------------------------------------------- old ledger and the update/add plan
old_rows = {}
with tarfile.open(OLD, 'r:gz') as t:
    for m in t:
        if m.name == f'{OLD_ROOT}/FILES.tsv':
            old_rows = {r['path']: r for r in csv.DictReader(io.TextIOWrapper(t.extractfile(m), 'utf-8'), delimiter='\t')}
            break
assert len(old_rows) == 2688, len(old_rows)
update = {}
for p, r in old_rows.items():
    if not p.startswith('project/'): continue
    f = ROOT/p[8:]
    b = f.read_bytes()
    if sha(b) != r['sha256']: update[p] = b
add_paths = [f'{RUN}/{n}' for n in ('README.zh.md', 'apply_edits.py', 'assemble_set.py', 'build_compose.py', 'promote.py', 'run_a01.py',
                                     'promote.log.json', 'promote-audit.log.json', 'package/build_v11c.py')]
add_paths += sorted(str(p.relative_to(ROOT)) for p in (ROOT/RUN/'edits').glob('*') if p.suffix in ('.py', '.json'))
add_paths += sorted(str(p.relative_to(ROOT)) for p in (ROOT/RUN/'out/set').glob('SuppFig.S*') if p.suffix in ('.pdf', '.json'))
add_paths += [f'{RUN}/out/set/manifest.json']
for k, stem in (('X19', 'SuppFig.S19-fi'), ('X20', 'SuppFig.S20-d')):          # inputs of mech.py (S19 f-i) and profile.py (S20 d)
    add_paths += [f'{RUN}/out/compose/{k}/{stem}.{e}' for e in ('pdf', 'receipt.json', 'values.json')]
add_paths += ['reviews/2026-10-12-supp-figure-caption-audit/REPORT.zh.md']
add = {}
for a in add_paths:
    key = 'project/'+a
    assert key not in old_rows, key
    add[key] = (ROOT/a).read_bytes()

README = (ROOT/RUN/'package/README.v11c.md').read_bytes()
LICENSE = (ROOT/'LICENSE').read_bytes()

# ---------------------------------------------------------------- final ledger content (stored_as filled while writing)
final = {}
for p, r in old_rows.items():
    if p in update: final[p] = dict(path=p, bytes=len(update[p]), sha256=sha(update[p]), origin=ORIGIN, origin_member=p[8:])
    elif p == 'README.md': final[p] = dict(path=p, bytes=len(README), sha256=sha(README), origin='consolidated release', origin_member='README.md')
    else: final[p] = dict(path=p, bytes=int(r['bytes']), sha256=r['sha256'], origin=r['origin'], origin_member=r['origin_member'])
for p, b in add.items(): final[p] = dict(path=p, bytes=len(b), sha256=sha(b), origin=ORIGIN, origin_member=p[8:])

# ---------------------------------------------------------------- stream-write
stored = {}                       # sha -> member path (relative to root) holding the payload
old_payload = {}                  # old payload of updated paths that other members link to
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
        if rel == 'FILES.tsv': continue                    # written last (after stored_as is known), same as other top files
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
    # ledger last
    cols = ['path', 'bytes', 'sha256', 'stored_as', 'origin', 'origin_member']
    buf = io.StringIO(); w = csv.DictWriter(buf, fieldnames=cols, delimiter='\t', lineterminator='\n'); w.writeheader()
    for p in sorted(final): w.writerow({c: final[p][c] for c in cols})
    ledger = buf.getvalue().encode()
    ti = tarfile.TarInfo(f'{NEW_ROOT}/FILES.tsv'); ti.size, ti.mode, ti.mtime = len(ledger), 0o644, now; out.addfile(ti, io.BytesIO(ledger))
(FILES/'FILES.tsv').write_bytes(ledger)
(FILES/'README.md').write_bytes(README); (FILES/'LICENSE').write_bytes(LICENSE)
json.dump(dict(updated=sorted(update), added=sorted(add), relinked_from_replaced_holders=sorted(r for r in final if final[r].get('stored_as') and final[r]['stored_as'] != r and old_rows.get(r, {}).get('stored_as') in update),
               members=len(final), unique_payloads=len(stored)), open(DST/'refresh-20261012.json', 'w'), indent=1)
print('members', len(final), 'payloads', len(stored), 'updated', len(update), 'added', len(add))
