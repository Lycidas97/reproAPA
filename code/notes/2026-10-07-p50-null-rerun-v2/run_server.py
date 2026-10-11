"""Repeat all five p50 nulls and masked ladder simulations in a fresh directory."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

root = Path(__file__).resolve().parent
prior = root.parent / 'reproapa_cross_overlap_null_a1_20261007_v1'
folder = root / 'apadetect-p50'
code = root / 'code'
python = sys.executable
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
assert not (root / 'status.json').exists(), 'Use a fresh run directory'
folder.mkdir()
environment = os.environ.copy()
environment.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
environment['PYTHONPATH'] = str(code)
available = os.sched_getaffinity(0)
os.sched_setaffinity(0, {min(available)})
status = dict(status='running', source='apadetect-p50', stages=[], prose_changes=False)


def save():
    (root / 'status.json').write_text(json.dumps(status, indent=2) + '\n')


def run(name, arguments):
    item = dict(stage=name, status='running', command=[python, '-u'] + list(map(str, arguments)))
    status['stages'].append(item)
    save()
    with (root / (name + '.log')).open('w') as log:
        result = subprocess.run(item['command'], env=environment, stdout=log, stderr=subprocess.STDOUT)
    item.update(status='complete' if result.returncode == 0 else 'failed', exit_code=result.returncode)
    save()
    if result.returncode:
        raise RuntimeError(name)


try:
    original_path = prior / 'apadetect-p50/config.json'
    original = json.loads(original_path.read_text())
    (folder / 'original-run.config.json').write_bytes(original_path.read_bytes())
    config = original.copy()
    names = {path.name for path in code.glob('*.py')}
    config['code_sha256'] = {key: value for key, value in original['code_sha256'].items() if Path(key).name not in names}
    config['code_sha256'].update({str(path): digest(path) for path in code.glob('*.py')})
    config['source_overrides'] = original['source_overrides'].copy()
    for relative, name in [
        ('notes/2026-09-29-core-sensitivity/figure4_matched_overlap.py', 'figure4_matched_overlap.py'),
        ('notes/2026-10-02-null-model-controls/common.py', 'common.py'),
        ('notes/2026-10-02-null-model-controls/cross_species.py', 'cross_species.py')]:
        path = code / 'numerics' / name
        config['source_overrides'][relative] = dict(path=str(path), sha256=digest(path))
    config['p50_null_rerun'] = dict(original_config_sha256=digest(original_path),
        correction='Actual measurable effect count in all five nulls, ladder and injection strata',
        immutable_extraction=str(prior / 'apadetect-p50/extraction'))
    config_path = folder / 'config.json'
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    (folder / 'source.config.json').write_bytes(config_path.read_bytes())
    (folder / 'extraction').symlink_to(prior / 'apadetect-p50/extraction')

    run('intent-tests', ['-m', 'unittest', 'test_cross_overlap_expression_coverage',
        'test_cross_overlap_null_eligibility', 'test_cross_controls_eligibility'])
    run('original-nulls', [code / 'downstream_cross.py', '--config', config_path,
        '--out', folder / 'original-nulls', '--node', 'cross_overlap', '--cached-overlap'])
    run('extended-nulls', [code / 'cross_overlap_expression_coverage.py', '--config', config_path,
        '--out', folder / 'inference', '--stage', 'infer', '--extraction', folder / 'extraction'])
    run('extended-verification', [code / 'verify_cross_overlap_expression_coverage.py',
        '--root', folder, '--config', config_path])
    run('all-five-verification', [code / 'verify_cross_overlap_null_eligibility.py',
        '--root', folder, '--config', config_path])
    run('grid-verification', [code / 'verify_cross_overlap_sensitivity.py', '--root', folder,
        '--config', config_path, '--out', root / 'validation/grid'])
    run('joint-exact', [code / 'verify_cross_overlap_joint_exact.py', '--root', folder,
        '--out', root / 'validation/exact'])

    (folder / 'analysis').mkdir()
    (folder / 'analysis/cross_overlap').symlink_to(folder / 'original-nulls')
    ladder_config = config.copy()
    ladder_config['analysis_root'] = str(folder / 'analysis')
    ladder_path = folder / 'ladder.config.json'
    ladder_path.write_text(json.dumps(ladder_config, indent=2) + '\n')
    run('masked-ladder', [code / 'downstream_cross.py', '--config', ladder_path,
        '--out', folder / 'ladder', '--node', 'cross_controls'])
    run('ladder-verification', [code / 'verify_ladder.py', '--root', root, '--source', 'apadetect-p50'])
    status['status'] = 'complete'
except Exception:
    status.update(status='failed', error=traceback.format_exc())
    raise
finally:
    save()
