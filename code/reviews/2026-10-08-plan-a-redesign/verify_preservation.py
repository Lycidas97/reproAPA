"""Verify the two frozen preceding deliverables without modifying them."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    checks = []
    for name in ['2026-10-08-plan-a-analysis', '2026-10-08-plan-a-statistical-strictness']:
        directory = ROOT / 'reviews' / name
        manifest_path = directory / 'completion-manifest.json'
        manifest = json.loads(manifest_path.read_text())
        files = manifest['files']
        for relative, record in files.items():
            path = directory / relative
            assert path.stat().st_size == record['bytes'], str(path)
            assert sha(path) == record['sha256'], str(path)
        sources = manifest.get('source_files', {})
        for relative, expected in sources.items():
            assert sha(ROOT / relative) == expected, relative
        checks.append(dict(deliverable=name, files_verified=len(files),
                           source_files_verified=len(sources),
                           manifest_sha256=sha(manifest_path), status='PASS'))
    receipt = dict(status='PASS', checks=checks,
                   total_original_files=sum(r['files_verified'] + r['source_files_verified'] for r in checks),
                   code_sha256=sha(Path(__file__)))
    (HERE / 'preservation-verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
