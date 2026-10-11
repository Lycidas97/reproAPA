"""Run explicitly bound statistical jobs after the native core producer finishes."""
import argparse
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from polyapipe_prepare import digest


def main(args):
    config_path=args.config.resolve();config=json.loads(config_path.read_text())
    root=Path(config['analysis_root']);root.mkdir(exist_ok=bool(args.phase_id),parents=True)
    status_name=('dispatch-'+args.phase_id+'-status') if args.phase_id else 'dispatch-status'
    bindings=config['jobs'];completed={};pending=set(bindings);active={}
    def status(phase,**extra):
        p=root/(status_name+'.tmp.json');p.write_text(json.dumps(dict(phase=phase,completed=list(completed),
            pending=sorted(pending),running=list(active.values()),full_manuscript_rerun_complete=False,**extra),indent=2)+'\n')
        p.replace(root/(status_name+'.json'))
    def core_ready():
        path=Path(config['native_run'])/'core_status.json'
        if not path.exists():
            os.kill(args.core_pid,0)
            return False
        core=json.loads(path.read_text())
        if core['phase']=='failed':raise RuntimeError(core)
        if core['phase']=='core_effects_complete_remaining_manuscript_analyses_pending':return True
        os.kill(args.core_pid,0)
        return False
    status('waiting_for_dependencies')
    def job(node):
        entry=bindings[node]
        for path,h in config['code_sha256'].items():assert digest(path)==h,f'Code changed: {path}'
        out=root/node;out.mkdir(exist_ok=False)
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
        with (out/'command.log').open('w') as log:
            command=[entry.get('python',sys.executable),entry['script'],'--config',str(config_path),'--out',str(out),*entry['arguments']]
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
        result=json.loads((out/'result.json').read_text())
        assert result['node']==node and result['status']=='complete' and all(result['checks'].values())
        for artifact in result['artifacts']:assert digest(out/artifact['path'])==artifact['sha256']
        return digest(out/'result.json')
    external=set().union(*(set(v['depends_on']) for v in bindings.values()))-set(bindings)
    def available_dependencies():
        for node in external-set(completed):
            path=root/node/'result.json'
            if not path.exists():continue
            receipt=json.loads(path.read_text())
            assert receipt['node']==node and receipt['status']=='complete' and all(receipt['checks'].values())
            assert receipt['run_id']==config['run_id'] and receipt['input_fingerprint']==config['input_fingerprint']
            for a in receipt['artifacts']:assert digest(path.parent/a['path'])==a['sha256']
            completed[node]=digest(path)
        for source in config.get('upstream_dispatch_status',[]):
            if Path(source).exists():
                upstream=json.loads(Path(source).read_text())
                if upstream['phase']=='failed':raise RuntimeError(dict(upstream_failure=source,detail=upstream))
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        while pending or active:
            available_dependencies()
            ready=core_ready()
            for node in list(bindings):
                if node not in pending or len(active)>=args.workers:continue
                if bindings[node].get('requires_core',True) and not ready:continue
                if not set(bindings[node]['depends_on'])<=set(completed):continue
                pending.remove(node);active[executor.submit(job,node)]=node
            status('running_bound_statistics')
            if not active:
                assert not ready or external-set(completed),'No runnable bound job; dependency cycle or missing binding'
                time.sleep(30);continue
            done,_=wait(active,timeout=30,return_when=FIRST_COMPLETED)
            for future in done:
                node=active.pop(future)
                try:completed[node]=future.result()
                except Exception as error:
                    status('failed',failed_node=node,error=repr(error));raise
    status('bound_statistics_complete',receipt_sha256=completed)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True)
    p.add_argument('--phase-id');p.add_argument('--core-pid',type=int,required=True);p.add_argument('--workers',type=int,default=2)
    args=p.parse_args()
    try:main(args)
    except Exception as error:
        config=json.loads(args.config.read_text());root=Path(config['analysis_root'])
        name='dispatch-'+args.phase_id+'-status.json' if args.phase_id else 'dispatch-status.json'
        target=root/name
        if target.exists():
            record=json.loads(target.read_text());record.update(phase='failed',error=repr(error))
            temporary=target.with_suffix('.tmp.json');temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(target)
        raise
