"""Verify sealed subanalyses and unchanged manuscript before delivery."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    cache = {}

    def digest(path):
        path = path.resolve()
        if path not in cache:
            with path.open('rb') as stream:
                cache[path] = hashlib.file_digest(stream, 'sha256').hexdigest()
        return cache[path]

    frozen = json.loads((HERE / 'freeze_manifest.json').read_text())
    assert digest(HERE / 'PROTOCOL.zh.md') == frozen['protocol_sha256']
    for name, expected in frozen['manuscript_unchanged_checks'].items():
        assert digest(ROOT / name) == expected, name
    verified = []
    for sub in ['sm-core', 'contrast-dose', 'wot-dose']:
        folder = HERE / sub
        manifest = json.loads((folder / 'manifest.json').read_text())
        for field in ['inputs', 'input_sha256', 'outputs', 'output_sha256', 'files']:
            for name, expected in manifest.get(field, {}).items():
                path = Path(name)
                if not path.is_absolute():
                    path = ROOT / path if name.startswith('notes/') else folder / path
                assert digest(path) == expected, str(path)
                verified.append(str(path.relative_to(ROOT)))
    files = ['input_validation.json', 'sm_raw_independent_validation.json',
             'sm_aggregate_independent_validation.json', 'contrast_independent_validation.json',
             'wot_independent_validation.json', 'sm-core/delivery_validation.json',
             'contrast-dose/validation.json', 'wot-dose/validation.json']
    for name in files:
        result = json.loads((HERE / name).read_text())
        # The WOT standalone checker uses its established lowercase status schema.
        expected_status = 'passed' if name == 'wot-dose/validation.json' else 'PASS'
        assert result['status'] == expected_status, name
        assert result.get('tests_skipped', 0) == 0, name
        if name == 'wot-dose/validation.json':
            assert len(result['checks']) == result['n_checks'] == 4957
            assert all(check['passed'] for check in result['checks'])
    validation = dict(status='PASS', protocol_sha256=frozen['protocol_sha256'],
                      root_independent_reconstruction='complete', validation_files=files,
                      subanalysis_hash_records_verified=len(verified), distinct_hashed_files=len(cache),
                      manuscript_source_hashes_unchanged=len(frozen['manuscript_unchanged_checks']),
                      v11_submitted_by_this_task=False, W2_started_by_this_task=False,
                      unavailable_comparisons_explicit=True, species_population_inference=False, tests_skipped=0)
    (HERE / 'validation.json').write_text(json.dumps(validation, ensure_ascii=False, indent=2) + '\n')
    outputs = {str(path.relative_to(HERE)): digest(path) for path in sorted(HERE.rglob('*'))
               if path.is_file() and '__pycache__' not in path.parts and path != HERE / 'manifest.json'}
    manifest = dict(status='complete_analysis_verified', protocol_sha256=frozen['protocol_sha256'],
                    verification_scope=validation, outputs=outputs,
                    manuscript_source_sha256=frozen['manuscript_unchanged_checks'])
    (HERE / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(validation, ensure_ascii=False, indent=2))
    print('root manifest files', len(outputs))


if __name__ == '__main__':
    main()
