"""Apply exact-once text replacements to figure scripts.
Usage: python3 apply_edits.py edits/<name>.py   (the file defines TARGET (path relative to project root) and EDITS = [(old, new), ...])
Each old string must occur exactly once in the current target; the result is logged to edits/<name>.applied.json.
Re-running after success is a no-op only if every new string is already present and no old string is."""
import sys, json, runpy, hashlib
from pathlib import Path
RUN = Path(__file__).resolve().parent; ROOT = RUN.parents[1]
spec = runpy.run_path(sys.argv[1]); tgt = ROOT/spec['TARGET']; s = tgt.read_text()
before = hashlib.sha256(s.encode()).hexdigest()
done = 0
for i, (old, new) in enumerate(spec['EDITS']):
    n = s.count(old)
    if n == 1:
        s = s.replace(old, new); done += 1
    elif n == 0 and new in s:
        print(f'edit {i}: already applied')
    else:
        raise SystemExit(f'edit {i}: old occurs {n} times in {tgt}:\n{old[:200]}')
tgt.write_text(s)
log = dict(target=spec['TARGET'], sha_before=before, sha_after=hashlib.sha256(s.encode()).hexdigest(), applied=done, total=len(spec['EDITS']))
Path(sys.argv[1]).with_suffix('.applied.json').write_text(json.dumps(log, indent=1))
print(log)
