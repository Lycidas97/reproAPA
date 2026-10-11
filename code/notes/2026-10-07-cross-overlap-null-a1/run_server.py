"""Complete polyApipe first, then rerun p50 from its own contributors."""
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

root = Path(__file__).resolve().parent
python = sys.executable
environment = os.environ.copy()
environment.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
os.sched_setaffinity(0, {0, 1})
status = dict(status='running', stages=[], source_order=['polyapipe', 'apadetect-p50'])


def save():
    (root / 'status.json').write_text(json.dumps(status, indent=2) + '\n')


def run(name, arguments):
    item = dict(stage=name, status='running', command=[python, '-u'] + arguments)
    status['stages'].append(item)
    save()
    with (root / (name + '.log')).open('w') as log:
        result = subprocess.run(item['command'], env=environment, stdout=log, stderr=subprocess.STDOUT)
    item.update(status='complete' if result.returncode == 0 else 'failed', exit_code=result.returncode)
    save()
    if result.returncode:
        raise RuntimeError(name)


try:
    run('intent-tests', [str(root / 'code/test_cross_overlap_expression_coverage.py')])
    for label in status['source_order']:
        folder = root / label
        config = str(folder / 'config.json')
        if label == 'apadetect-p50':
            run(label + '-extract', [str(root / 'code/cross_overlap_expression_coverage.py'),
                '--config', config, '--out', str(folder / 'extraction'), '--stage', 'extract'])
        run(label + '-original-nulls', [str(root / 'code/downstream_cross.py'),
            '--config', config, '--out', str(folder / 'original-nulls'), '--node', 'cross_overlap', '--cached-overlap'])
        run(label + '-extended-nulls', [str(root / 'code/cross_overlap_expression_coverage.py'),
            '--config', config, '--out', str(folder / 'inference'), '--stage', 'infer',
            '--extraction', str(folder / 'extraction')])
        run(label + '-extended-verification', [str(root / 'code/verify_cross_overlap_expression_coverage.py'),
            '--root', str(folder), '--config', config])
        run(label + '-all-five-verification', [str(root / 'code/verify_cross_overlap_null_eligibility.py'),
            '--root', str(folder), '--config', config])
    status['status'] = 'complete'
except Exception:
    status['status'] = 'failed'
    status['error'] = traceback.format_exc()
    raise
finally:
    save()
