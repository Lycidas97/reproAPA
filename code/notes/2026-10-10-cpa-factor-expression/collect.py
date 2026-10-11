"""Freeze read-only factor extraction, without changing upstream analysis files."""
from pathlib import Path
import hashlib
import json
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    config_path = ROOT / 'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'
    config = json.loads(config_path.read_text())
    protocol = json.loads((HERE / 'protocol.json').read_text())
    code = 'CONFIG = ' + repr(config) + '\nPROTOCOL = ' + repr(protocol) + '\n' + (HERE / 'extract_remote.py').read_text()
    partial = HERE / 'expression_inputs.partial.npz'
    command = 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/scanpy/bin/python -'
    with partial.open('wb') as output, (HERE / 'extract.log').open('wb') as log:
        subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', '213docker', command], input=code.encode(), stdout=output, stderr=log, check=True)
    partial.rename(HERE / 'expression_inputs.npz')
    print((HERE / 'extract.log').read_text())
    sha = hashlib.sha256(config_path.read_bytes()).hexdigest()
    (HERE / 'collection_provenance.json').write_text(json.dumps({'config_path': str(config_path), 'config_sha256': sha, 'remote_read_only': True, 'Liu_GSE150311': 'supplied by the parallel all-gene naive-expression task'}, indent=2) + '\n')


if __name__ == '__main__':
    main()
