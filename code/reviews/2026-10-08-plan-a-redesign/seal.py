"""Seal verified redesign outputs and check their source receipts and links."""
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    receipts = [
        'background/results/verification.json',
        'background/cross-audit.json',
        'cross-species/verification.json',
        'cross-species/background-audit-receipt.json',
        'cycle/validation.json',
        'background/cycle-final-audit.json',
        'background/cycle-final-results-audit.json',
        'background/cycle-final-code-audit.json',
        'cross-species/cycle-numerical-audit.json',
        'correction-verification.json',
        'preservation-verification.json',
    ]
    for relative in receipts:
        item = json.loads((HERE / relative).read_text())
        assert item['status'].lower() in ['pass', 'passed'], relative
    hashes = []
    for relative in ['background/results/input-audit.json', 'cross-species/input-audit.json']:
        item = json.loads((HERE / relative).read_text())
        hashes.extend(item['input_sha256'].items())
    hashes.extend(json.loads((HERE / 'cycle/input_sha256.json').read_text()).items())
    hashes.extend(json.loads((HERE / 'cycle/execution-receipt.json').read_text())['sha256'].items())
    for path, expected in hashes:
        assert sha(Path(path)) == expected, path
    cycle_manifest = json.loads((HERE / 'cycle/artifact-manifest.json').read_text())
    for entry in cycle_manifest['files']:
        path = HERE / 'cycle' / entry['path']
        assert path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256'], str(path)
    cycle_complete = json.loads((HERE / 'cycle/completion-receipt.json').read_text())
    assert sha(HERE / 'cycle/artifact-manifest.json') == cycle_complete['manifest_sha256']
    assert sha(HERE / 'cycle/REPORT.zh.md') == cycle_complete['report_sha256']
    cross_manifest = json.loads((HERE / 'cross-species/deliverable-manifest.json').read_text())
    for relative, expected in cross_manifest['files_SHA256'].items():
        assert sha(HERE / 'cross-species' / relative) == expected, relative
    correction = json.loads((HERE / 'correction-verification.json').read_text())
    assert sha(HERE / 'integrate.py') == correction['code_sha256']
    for relative, expected in correction['input_sha256'].items():
        assert sha(HERE / relative) == expected, relative
    assert sha(HERE / 'PROTOCOL.md') == '948de02088f563425c95120b2749f881c28318692678b75da7cfbdc0597bce82'
    assert sha(HERE / 'CYCLE-FREEZE-ADDENDUM.md') == '6c1b1e29122e0e554f56dec9acaac0badbb3b84f12590b7a847889cc1b140743'
    report = (HERE / 'REPORT.zh.md').read_text()
    assert '尚未封存' not in report and '待模块' not in report and '待全模块' not in report
    links = []
    for path in HERE.rglob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            target = target.strip('<>')
            if not target.startswith(('https://', 'http://', '#', 'mailto:', 'app://', 'codex://')):
                local = re.sub(r':\d+$', '', target)
                file = Path(local) if local.startswith('/') else (path.parent / local).resolve()
                if file != HERE / 'completion-manifest.json':
                    assert file.is_file(), (path, target)
                links.append(file)
    tests = pd.read_csv(HERE / 'corrected-tests.tsv', sep='\t')
    families = pd.read_csv(HERE / 'test-families.tsv', sep='\t')
    assert len(tests) == 46 and tests.test_id.is_unique
    assert tests.loc[tests.raw_p.isna(), 'q_BH_family'].isna().all()
    files = {}
    for path in sorted(HERE.rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.name != 'completion-manifest.json':
            files[str(path.relative_to(HERE))] = dict(bytes=path.stat().st_size, sha256=sha(path))
    manifest = dict(status='complete', finished_utc=datetime.now(timezone.utc).isoformat(),
                    rerun_items=[1, 2, 3, 8], primary='polyApipe contract-r2',
                    alternative='APADETECT point filter plus p50',
                    post_initial_results_method_redesign=True,
                    manuscript_modified=False, full_PAS_recalling=False,
                    scientific_thresholds_frozen_before_module_inference=True,
                    implementation_repairs_documented=True,
                    planned_positions=46, estimable_p=int(tests.raw_p.notna().sum()),
                    unestimable_NA=int(tests.raw_p.isna().sum()),
                    families=families.to_dict('records'), checks=receipts,
                    receipt_input_hash_entries_verified=len(hashes),
                    module_manifest_files_verified={
                        'cycle': len(cycle_manifest['files']),
                        'cross_species': len(cross_manifest['files_SHA256'])},
                    local_report_links_verified=len(links),
                    ignored_generated_files=['__pycache__'], files=files)
    (HERE / 'completion-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    assert all(path.is_file() for path in links)
    print(json.dumps({k: v for k, v in manifest.items() if k != 'files'}, indent=2))
    print('Sealed files:', len(files))


if __name__ == '__main__':
    main()
