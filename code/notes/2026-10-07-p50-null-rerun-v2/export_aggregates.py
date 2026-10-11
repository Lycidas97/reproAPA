"""Export a fixed whitelist of aggregate p50 results after every check passes."""
import hashlib
import json
from pathlib import Path
import shutil

root = Path(__file__).resolve().parent
status = json.loads((root / 'status.json').read_text())
assert status['status'] == 'complete'
assert all(stage['status'] == 'complete' and stage['exit_code'] == 0 for stage in status['stages'])
assert 'Ran 8 tests' in (root / 'intent-tests.log').read_text()
assert (root / 'intent-tests.log').read_text().rstrip().endswith('OK')
paths = ['apadetect-p50/' + name for name in [
    'all-five-null-comparison.tsv', 'all-five-verification.json', 'independent-verification.json',
    'mobility-diagnostics.tsv', 'original-nulls/overlap_null_results.tsv',
    'original-nulls/overlap_threshold_layer_sensitivity.tsv',
    'inference/expression-coverage-null-results.tsv', 'inference/balance.tsv',
    'inference/calibration-power.tsv', 'inference/bin-definitions.tsv',
    'inference/coverage-and-missingness.tsv', 'inference/original-model-replay.tsv',
    'ladder/cross_species_ladder.tsv', 'ladder/cross_species_power.tsv',
    'ladder/cross_species_balance.tsv']]
paths += ['eligibility-verification.tsv', 'verification.json', 'intent-tests.log',
          'validation/grid/sensitivity-verification.json',
          'validation/exact/exact-verification.json', 'validation/exact/joint-exact-results.tsv']
for relative in ['apadetect-p50/all-five-verification.json', 'apadetect-p50/independent-verification.json',
                 'verification.json', 'validation/grid/sensitivity-verification.json',
                 'validation/exact/exact-verification.json']:
    assert json.loads((root / relative).read_text())['status'] == 'PASS', relative
out = root / 'aggregate-handoff'
assert not out.exists()
entries = []
for relative in paths:
    path = root / relative
    if path.suffix == '.tsv':
        with path.open() as handle:
            header = handle.readline().rstrip().split('\t')
        assert not {'gene_id', 'human_gene_id', 'mouse_gene_id', 'metacell_id',
                    'metacell_ids', 'barcode', 'library_key'} & set(header), relative
    target = out / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
    entries.append(dict(path=relative, bytes=target.stat().st_size,
        sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
config = json.loads((root / 'apadetect-p50/config.json').read_text())
manifest = dict(status='PASS', source='apadetect-p50',
    identity={key: config[key] for key in ['run_id', 'input_fingerprint', 'reference_variant']},
    null_tests=10, sensitivity_settings=16, ladder_tests=10, ladder_simulation_rows=90,
    intent_tests=8, tests_skipped=0, stages=[stage['stage'] for stage in status['stages']],
    native_covariate_extraction_reused_from_p50=True,
    code_sha256={str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in sorted((root / 'code').rglob('*.py')) if not path.name.startswith('._')},
    contains_identifiers_or_contributor_data=False, contains_matrices_or_permutation_arrays=False,
    server_directory=str(root), prose_changes=False, files=entries)
(out / 'HANDOFF-MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({key: value for key, value in manifest.items() if key != 'files'}, indent=2))
print('files', len(entries), 'bytes', sum(item['bytes'] for item in entries))
