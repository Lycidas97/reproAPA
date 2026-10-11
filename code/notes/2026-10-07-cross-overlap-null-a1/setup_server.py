"""Bind corrected scripts to the original, immutable source runs."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
code = root / 'code'
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
for label in ['polyapipe', 'apadetect-p50']:
    config = json.loads((root / (label + '.source.config.json')).read_text())
    original_config_sha256 = digest(root / (label + '.source.config.json'))
    names = [path.name for path in code.glob('*.py')]
    config['code_sha256'] = {key: value for key, value in config['code_sha256'].items()
                             if Path(key).name not in names}
    config['code_sha256'].update({str(path): digest(path) for path in code.glob('*.py')})
    for relative, name in [
        ('notes/2026-09-29-core-sensitivity/figure4_matched_overlap.py', 'figure4_matched_overlap.py'),
        ('notes/2026-10-02-null-model-controls/common.py', 'common.py'),
        ('notes/2026-10-02-null-model-controls/cross_species.py', 'cross_species.py'),
    ]:
        path = code / 'numerics' / name
        config['source_overrides'][relative] = dict(path=str(path), sha256=digest(path))
    config['a1_rerun'] = dict(original_config_sha256=original_config_sha256,
                              correction='Add actual measurable effect study count to each species null stratum',
                              original_analysis_root=config['analysis_root'])
    destination = root / label
    destination.mkdir(exist_ok=True)
    (destination / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
    if label == 'polyapipe' and not (destination / 'extraction').exists():
        (destination / 'extraction').symlink_to(
            '/root/nfsdata/apadetect/analysis_projects/reproapa_cross_overlap_expression_coverage_20261006_v1/extraction')

