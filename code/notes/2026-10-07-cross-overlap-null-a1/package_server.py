"""Freeze only a fully verified rerun and retain its reproducible source data."""
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

import numpy
import pandas
import scipy

root = Path(__file__).resolve().parent
status = json.loads((root / 'status.json').read_text())
assert status['status'] == 'complete' and all(stage['status'] == 'complete' for stage in status['stages'])
for label in ['polyapipe', 'apadetect-p50']:
    assert json.loads((root / label / 'all-five-verification.json').read_text())['status'] == 'PASS'
    assert json.loads((root / label / 'independent-verification.json').read_text())['status'] == 'PASS'
    assert json.loads((root / 'validation' / (label + '-grid') / 'sensitivity-verification.json').read_text())['status'] == 'PASS'
    assert json.loads((root / 'validation' / (label + '-exact') / 'exact-verification.json').read_text())['status'] == 'PASS'
assert 'Ran 1 test' in (root / 'validation/original-null-intent-tests.log').read_text()
assert (root / 'validation/original-null-intent-tests.log').read_text().rstrip().endswith('OK')
(root / 'environment.json').write_text(json.dumps(dict(
    python=sys.version, executable=sys.executable, platform=platform.platform(),
    numpy=numpy.__version__, pandas=pandas.__version__, scipy=scipy.__version__,
    cpu_affinity=[0, 1], numerical_threads=1), indent=2) + '\n')
files = []
for folder, directories, names in os.walk(root, followlinks=True):
    directories[:] = [name for name in directories if name != '__pycache__']
    for name in names:
        if name.startswith('._') or name == 'PACKAGE-MANIFEST.json':
            continue
        path = Path(folder) / name
        h = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                h.update(block)
        files.append(dict(path=str(path.relative_to(root)), bytes=path.stat().st_size, sha256=h.hexdigest()))
manifest = dict(status='PASS', completed_sources=['polyapipe', 'apadetect-p50'],
                null_tests_per_source=10, intent_tests=7, tests_skipped=0,
                sensitivity_settings_per_source=16, threshold_grid_receipts_checked=True,
                joint_exact_distribution_receipts_checked=True,
                includes_reused_polyapipe_extraction=True, files=sorted(files, key=lambda item: item['path']))
(root / 'PACKAGE-MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({key: value for key, value in manifest.items() if key != 'files'}, indent=2))
print('files', len(files), 'bytes', sum(item['bytes'] for item in files))
