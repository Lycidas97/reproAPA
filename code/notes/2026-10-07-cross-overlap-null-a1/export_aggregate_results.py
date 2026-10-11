"""Export only aggregate statistics and verification reports, with no identifiers."""
import hashlib
import json
from pathlib import Path
import shutil

root = Path(__file__).resolve().parent
destination = root / 'aggregate-handoff'
assert not destination.exists()
assert json.loads((root / 'status.json').read_text())['status'] == 'complete'
paths = []
for source in ['polyapipe', 'apadetect-p50']:
    paths += [source + '/' + name for name in [
        'all-five-null-comparison.tsv', 'all-five-verification.json',
        'independent-verification.json', 'mobility-diagnostics.tsv',
        'original-nulls/overlap_null_results.tsv',
        'original-nulls/overlap_threshold_layer_sensitivity.tsv',
        'inference/expression-coverage-null-results.tsv',
        'inference/balance.tsv', 'inference/calibration-power.tsv',
        'inference/bin-definitions.tsv', 'inference/coverage-and-missingness.tsv',
        'inference/original-model-replay.tsv',
    ]]
    paths += ['validation/' + source + suffix for suffix in [
        '-grid/sensitivity-verification.json', '-exact/exact-verification.json',
        '-exact/joint-exact-results.tsv',
    ]]
entries = []
for relative in paths:
    path = root / relative
    if path.suffix == '.tsv':
        header = path.open().readline().rstrip().split('\t')
        assert not {'gene_id', 'human_gene_id', 'mouse_gene_id', 'metacell_id',
                    'metacell_ids', 'barcode', 'library_key'} & set(header)
    target = destination / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
    entries.append(dict(path=relative, bytes=target.stat().st_size,
                        sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
sources = {}
for source in ['polyapipe', 'apadetect-p50']:
    config = json.loads((root / source / 'config.json').read_text())
    sources[source] = {key: config[key] for key in ['run_id', 'input_fingerprint', 'reference_variant']}
manifest = dict(status='PASS', payload='aggregate_statistics_and_verification_reports',
                sources=sources, null_tests_per_source=10, intent_tests=7, tests_skipped=0,
                sensitivity_settings_per_source=16, contains_contributor_data=False,
                contains_cell_or_gene_identifiers=False,
                contains_count_matrices_or_permutation_arrays=False,
                complete_server_package_manifest_sha256=hashlib.sha256((root / 'PACKAGE-MANIFEST.json').read_bytes()).hexdigest(),
                files=entries)
(destination / 'HANDOFF-MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({key: value for key, value in manifest.items() if key != 'files'}, indent=2))
print('files', len(entries), 'bytes', sum(item['bytes'] for item in entries))
