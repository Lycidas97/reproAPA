"""2026-10-12 review fix (item 5): rewrite the 'source' fields of the existing S20 receipts as project-relative paths,
matching what profile.py now writes (edits/review_profile_receipt.py). PDFs are not touched; a rerun of profile.py to a
scratch folder gave identical pixels (200 dpi) and identical text to Fig/SuppFig.S20.pdf, only the PDF bytes differ.
Usage: python3 <this file>   (from the project root); log in review_receipt_paths.applied.json"""
import json, os, hashlib
from pathlib import Path
ROOT = Path.cwd(); RUN = 'notes/2026-10-12-supp-figure-text-cleanup/'
F2 = 'revision-backups/20261011-before-v11-apply/Fig/Figure.2.pdf'
def fix(src):
    if src == 'Fig/Figure.2.pdf': return F2
    if '/mnt/reproAPA/' in src: return os.path.normpath(src.split('/mnt/reproAPA/', 1)[1])
    return src
log = []
for f in ['out/s20c/SuppFig.S20.receipt.json', 'out/set/SuppFig.S20.receipt.json', 'out/v11/SuppFig.S20.receipt.json']:
    p = ROOT/RUN/f; raw = p.read_bytes(); d = json.loads(raw); ch = []
    for k in d['kept']:
        new = fix(k['source'])
        if new != k['source']:
            assert (ROOT/new).exists(), (f, new)
            same = hashlib.sha256((ROOT/new).read_bytes()).hexdigest() == k['source_sha']
            # out/v11 is the 09:27 intermediate, superseded by out/s20c; its d source (X20 of 09:23) was rebuilt at 09:58,
            # so the file at that path no longer has the recorded source_sha. The path is fixed, the SHA is left as recorded.
            assert same or f.startswith('out/v11/'), (f, new)
            ch.append([k['source'], new, 'source_sha matches' if same else 'file since rebuilt; recorded source_sha kept']); k['source'] = new
    p.write_text(json.dumps(d, indent=1))
    log.append(dict(file=RUN+f, sha_before=hashlib.sha256(raw).hexdigest(), sha_after=hashlib.sha256(p.read_bytes()).hexdigest(), changed=ch))
Path(__file__).with_suffix('.applied.json').write_text(json.dumps(log, indent=1))
for r in log: print(r['file'], len(r['changed']))
