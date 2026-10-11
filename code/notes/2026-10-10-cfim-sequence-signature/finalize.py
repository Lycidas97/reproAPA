"""Freeze the reviewed, complete H2 delivery without changing results."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    checks=json.loads((HERE/'verification.json').read_text());assert checks['status']=='PASS' and checks['skipped_checks']==0
    tests=(HERE/'tests.log').read_text();assert 'Ran 7 tests' in tests and tests.rstrip().endswith('OK')
    protocol=json.loads((HERE/'protocol.json').read_text());assert digest(HERE/'PROTOCOL.zh.md')==protocol['protocol_sha256']
    audit=json.loads((HERE/'analysis_audit.json').read_text())
    for rel,h in audit['inputs'].items():assert digest(ROOT/rel)==h,rel
    files=[p for p in HERE.iterdir() if p.is_file() and p.name!='manifest.json' and '.partial.' not in p.name]
    manifest=dict(status='complete',scope='frozen internal CFIm sequence susceptibility; exploratory; no external knockdown; no manuscript edits',protocol_sha256=protocol['protocol_sha256'],tests={'passed':7,'skipped':0},verification=checks,interpretation='No primary multiplicity-adjusted support; human paired proximal identities fixed in primary strata. Sequence structure is not CFIm activity.',artifacts=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p),bytes=p.stat().st_size) for p in sorted(files)])
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'status':manifest['status'],'artifacts':len(files),'tests':manifest['tests']},indent=2))


if __name__=='__main__':main()
