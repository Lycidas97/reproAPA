"""Wait for audited annotation, then run all seven study-native core analyses.

This explicitly stops at core effects. It never marks the manuscript's 49-node
rerun complete or stages figures before the remaining adapters have executed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
from polyapipe_prepare import digest


def status(root, phase, **extra):
    record = dict(phase=phase, time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                  full_manuscript_rerun_complete=False, **extra)
    temp = root/'core_status.tmp.json'
    temp.write_text(json.dumps(record, indent=2)+'\n')
    temp.replace(root/'core_status.json')
    print(json.dumps(record), flush=True)


def main(args):
    root=args.run.resolve()
    code=Path(__file__).parent.resolve()
    source_hashes={str(p):digest(p) for p in sorted(code.glob('polyapipe_*.py'))}
    frozen=json.loads((root/'intake.native.json').read_text())
    assert len(frozen['studies']) == 7 and sum(x['libraries'] for x in frozen['studies']) == 161
    status(root,'waiting_for_structural_annotation', annotation_pid=args.annotation_pid,
           code_sha256=source_hashes)
    marker=root/'annotation.complete.json'
    while not marker.exists():
        try:
            os.kill(args.annotation_pid,0)
        except ProcessLookupError:
            raise RuntimeError('Annotation process exited without completion; inspect annotation.log')
        time.sleep(30)
    annotation=json.loads(marker.read_text())
    assert annotation['status']=='complete' and annotation['rows']==sum(x['atlas_pas'] for x in frozen['studies'])
    assert all(digest(p)==h for p,h in source_hashes.items()), 'Core code changed while annotation was running'
    folders=sorted((root/'studies').iterdir())
    assert len(folders)==7
    status(root,'running_core_effects', studies=[p.name for p in folders], workers=args.workers,
           settings=['supported50_p08','supported50_p09','supported50_p08_mc5','expanded_primary_p08'])

    def run(folder):
        with (folder/'effects.log').open('w') as log:
            subprocess.run([sys.executable,str(code/'polyapipe_effects.py'),'--study',str(folder)],
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        result=json.loads((folder/'effects/complete.json').read_text())
        assert result['status']=='core_effects_complete_not_full_rerun'
        return dict(study=folder.name, summaries=result['summaries'])

    completed=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run,p) for p in folders]
        for future in as_completed(futures):
            completed.append(future.result())
            status(root,'running_core_effects', completed=completed, total=7)
    status(root,'core_effects_complete_remaining_manuscript_analyses_pending', completed=completed,
           code_sha256=source_hashes, figures_replaced=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    p.add_argument('--annotation-pid',type=int,required=True);p.add_argument('--workers',type=int,default=2)
    args=p.parse_args()
    try:
        main(args)
    except Exception as error:
        status(args.run,'failed',error=repr(error))
        raise
