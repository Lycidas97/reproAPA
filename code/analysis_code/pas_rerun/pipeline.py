"""Caller-independent statistical DAG. Rendering is a separate consumer.

Existing --config/--out scientific adapters retain their estimators. Resume is
allowed only for matching inputs, code, dependency receipts and output hashes.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)


def topological(jobs):
    done, order = set(), []
    for n,j in jobs.items():
        require(set(j['depends_on']) <= set(jobs), f'Unbound dependency for {n}')
    while len(done) < len(jobs):
        ready = [n for n,j in jobs.items() if n not in done and set(j['depends_on']) <= done]
        require(bool(ready), 'Dependency cycle')
        order.extend(ready); done.update(ready)
    return order


def verified_artifacts(root, result):
    require(result.get('status') == 'complete' and bool(result.get('checks')) and
            all(v is True for v in result['checks'].values()), 'Invalid scientific receipt')
    require(bool(result.get('artifacts')), 'Missing scientific artifacts')
    for a in result['artifacts']:
        path = (root / a['path']).resolve()
        require(path.is_relative_to(root.resolve()), 'Artifact escapes node directory')
        require(digest(path) == a['sha256'], f'Changed node artifact: {path}')


def job_fingerprint(config, node, dependencies):
    payload = dict(run_id=config['run_id'], inputs=config['input_fingerprint'],
                   job=config['jobs'][node], config=config, dependencies=dependencies)
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def validate_cached(root, config, node, dependencies):
    marker = root / node / 'pipeline-receipt.json'
    if not marker.exists():
        return False
    record = read(marker)
    require(record['fingerprint'] == job_fingerprint(config,node,dependencies), f'Stale cache: {node}; use a new run directory')
    result_path = root / node / 'result.json'
    require(digest(result_path) == record['result_sha256'], f'Changed result receipt: {node}')
    result=read(result_path)
    require(result['run_id'] == config['run_id'] and result['input_fingerprint'] == config['input_fingerprint'],
            f'Cross-run receipt: {node}')
    verified_artifacts(root/node,result)
    return True


def validate_inputs(config):
    require(config.get('prepared') is True, 'Unbound plan: prepare with a complete native input before running')
    for path,h in config['code_sha256'].items():
        require(digest(path)==h, f'Code changed: {path}')
    for path,h in config['input_files'].items():
        require(digest(path)==h, f'Input changed: {path}')
    for key in [k for k in config.get('resources',{}) if k.endswith('_resume')]:
        require(key not in config.get('resources',{}), f'Old PAS cache forbidden: {key}')
    require(not any(k in config for k in ['inherited_previous_nodes','upstream_dispatch_status']), 'Inherited old-run state forbidden')


def cpu_slots(config,node):
    children=config.get('node_workers',{}).get(node,1)
    if node=='core':children=config.get('core_workers',2)
    if node=='count_inference':children=config.get('count_workers',4)
    require(isinstance(children,int) and children>=1, f'Invalid worker count: {node}')
    return children+1 if children>1 else 1


def execute(config_path, workers=2, selected=None, cpu_budget=None):
    config_path = Path(config_path).resolve(); config=read(config_path)
    order=topological(config['jobs']); validate_inputs(config)
    budget=cpu_budget or config.get('cpu_budget')
    require(budget is None or isinstance(budget,int) and budget>=1,'Invalid CPU budget')
    if budget is not None:
        require(all(cpu_slots(config,n)<=budget for n in order),'A node exceeds the CPU budget')
        if hasattr(os,'sched_getaffinity'):
            allowed=sorted(os.sched_getaffinity(0))
            os.sched_setaffinity(0,allowed[:min(budget,len(allowed))])
    if selected:
        wanted=set(selected)
        require(wanted<=set(order), 'Unknown selected node')
        while True:
            expanded=wanted|{d for n in wanted for d in config['jobs'][n]['depends_on']}
            if expanded==wanted:break
            wanted=expanded
        order=[n for n in order if n in wanted]
    root=Path(config['analysis_root']); root.mkdir(parents=True,exist_ok=True)
    lock_path=root/'pipeline.lock'
    # Cross-process exclusive lock; no duplicate scheduler on the same outputs.
    import fcntl
    with lock_path.open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        done={}; pending=set(order); active={}
        def status(phase,**extra):
            write(root/'pipeline-status.json',dict(phase=phase,run_id=config['run_id'],
                completed=list(done),pending=sorted(pending),running=list(active.values()),
                cpu_budget=budget,reserved_cpu_slots=sum(cpu_slots(config,n) for n in active.values()),
                full_manuscript_complete=False,**extra))
        def job(node,dependencies):
            if validate_cached(root,config,node,dependencies):return digest(root/node/'result.json')
            out=root/node
            require(not out.exists(), f'Incomplete output {node}; preserve it and use a fresh run or explicit repair')
            out.mkdir()
            j=config['jobs'][node]
            command=[j.get('python',sys.executable),j['script'],'--config',str(config_path),'--out',str(out),*j.get('arguments',[])]
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
            with (out/'command.log').open('w') as log:
                subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            result=read(out/'result.json')
            require(result['node']==node and result['run_id']==config['run_id'] and
                    result['input_fingerprint']==config['input_fingerprint'], f'Wrong result identity: {node}')
            verified_artifacts(out,result)
            write(out/'pipeline-receipt.json',dict(fingerprint=job_fingerprint(config,node,dependencies),
                result_sha256=digest(out/'result.json'),dependencies=dependencies,command=command))
            return digest(out/'result.json')
        status('running')
        failed_node=None
        try:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                while pending or active:
                    for node in order:
                        dependencies=config['jobs'][node]['depends_on']
                        used=sum(cpu_slots(config,n) for n in active.values())
                        if (node in pending and set(dependencies)<=set(done) and len(active)<workers and
                            (budget is None or used+cpu_slots(config,node)<=budget)):
                            pending.remove(node)
                            active[pool.submit(job,node,{d:done[d] for d in dependencies})]=node
                    require(bool(active), 'No runnable node')
                    status('running')
                    finished,_=wait(active,return_when=FIRST_COMPLETED)
                    for f in finished:
                        node=active.pop(f)
                        try:done[node]=f.result()
                        except Exception as error:
                            failed_node=node
                            # Publish the failure before executor shutdown waits for siblings.
                            status('failed',error=repr(error),failed_node=node,draining_nodes=list(active.values()))
                            raise
                    status('running')
            validate_inputs(config)
            status('statistics_complete' if selected is None else 'selected_nodes_complete', receipt_sha256=done,
                   rendering_performed=False)
        except Exception as error:
            status('failed',error=repr(error),failed_node=failed_node);raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['plan','run'])
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--cpu-budget',type=int,help='Total CPU slots, including nested study workers')
    p.add_argument('--nodes',help='Comma-separated targets; dependencies are included')
    a=p.parse_args()
    if a.command=='plan':
        c=read(a.config);print(json.dumps(dict(run_id=c['run_id'],order=topological(c['jobs'])),indent=2))
    else:
        require(a.workers>=1,'workers must be positive')
        execute(a.config,a.workers,a.nodes.split(',') if a.nodes else None,a.cpu_budget)
