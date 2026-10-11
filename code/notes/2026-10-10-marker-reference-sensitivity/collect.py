"""Collect immutable inputs from the configured server without changing it."""
from pathlib import Path
import hashlib
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    config_path = ROOT / 'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'
    config = json.loads(config_path.read_text())
    code = 'CONFIG = ' + repr(config) + '\n' + (HERE / 'extract_remote.py').read_text()
    partial = HERE / 'inputs.partial.npz'
    with partial.open('wb') as out, (HERE / 'extract.log').open('wb') as log:
        subprocess.run([
            'ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', '213docker',
            'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 '
            '/root/miniconda3/envs/scanpy/bin/python -'
        ], input=code.encode(), stdout=out, stderr=log, check=True)
    partial.rename(HERE / 'inputs.npz')
    (HERE / 'collection_manifest.json').write_text(json.dumps({
        'scope': 'all native MC; fixed marker union; all primary PAS counts; no APA effects',
        'config_sha256': hashlib.sha256(config_path.read_bytes()).hexdigest(),
        'extract_script_sha256': hashlib.sha256((HERE / 'extract_remote.py').read_bytes()).hexdigest(),
        'collector_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'inputs_sha256': hashlib.sha256((HERE / 'inputs.npz').read_bytes()).hexdigest()
    }, indent=2) + '\n')
    print((HERE / 'extract.log').read_text())


if __name__ == '__main__':
    main()
