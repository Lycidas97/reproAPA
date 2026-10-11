"""Launch explicit analysis entry points with durable stdout/stderr and PID receipts."""
import os
import subprocess
import sys
from common import *

scripts={'A01':'a01_quality.py','A02':'a02_graphs.py','A03':'a03_features.py',
         'A04':'a04_annotation.py','A05':'a05_apa.py','A06':'a06_doublets.py','A07':'a07_partitions.py'}
assert (RUN/'shared/ready.json').exists()
registry_path=RUN/'launch_registry.json'
registry=json.loads(registry_path.read_text()) if registry_path.exists() else {}
retry='--retry' in sys.argv
for job in [v for v in sys.argv[1:] if v!='--retry']:
    assert job in scripts, job
    history=[]
    if job in registry:
        assert retry, f'Already submitted: {job}'
        old=registry[job]
        assert json.loads((RUN/job/'status.json').read_text())['phase']=='FAILED'
        proc=Path('/proc')/str(old['pid'])/'cmdline'
        assert not proc.exists() or not proc.read_bytes(), 'Do not replace a live worker.'
        history=old.get('prior_attempts',[])+[{k:v for k,v in old.items() if k!='prior_attempts'}]
        os.replace(RUN/job/'run.log',RUN/job/f'run.attempt{len(history)}.log')
    script=RUN/'runtime'/scripts[job];assert script.is_file()
    runtime='/root/miniconda3/envs/single-cell-qc/bin/python' if job=='A07' else str(MOUSE/'.g2g-venv/bin/python')
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_NUM_THREADS='2',
               VECLIB_MAXIMUM_THREADS='2',CUDA_VISIBLE_DEVICES='',PYTHONUNBUFFERED='1')
    if job!='A07':
        env['LD_LIBRARY_PATH']='/root/miniconda3/envs/saturn-alignment/lib'+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    folder=RUN/job;folder.mkdir(exist_ok=True)
    command=['nice','-n','19','ionice','-c','3',runtime,str(script)]
    with (folder/'run.log').open('xb') as log:
        process=subprocess.Popen(command,cwd=RUN/'runtime',env=env,stdout=log,stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL,start_new_session=True)
    registry[job]=dict(pid=process.pid,command=command,script_sha256=sha(script),submitted_unix=time.time(),
        log=str(folder/'run.log'),status=str(folder/'status.json'),threads_per_job=2,prior_attempts=history)
    write_json(registry_path,registry)
    print(json.dumps({job:registry[job]}),flush=True)
