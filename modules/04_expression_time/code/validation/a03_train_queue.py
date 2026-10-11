"""Run the fixed four-arm scVI matrix only when a GPU is unoccupied."""
import os
import subprocess
import sys
import time
from common import *


def available():
    text=subprocess.check_output(['nvidia-smi','--query-gpu=index,memory.used,utilization.gpu',
        '--format=csv,noheader,nounits'],text=True)
    return [int(v[0]) for line in text.strip().splitlines()
            if len(v:=[s.strip() for s in line.split(',')])==3 and int(v[1])<1024 and int(v[2])<10]


def run():
    subprocess.run([sys.executable,str(RUN/'runtime/prepare_balanced_training_pilot.py')],check=True)
    arms=[(features,batch) for features in ['original','symmetric'] for batch in ['species','species_study']]
    for scope,seeds in [('pilot',[SEED]),('full',[SEED,SEED+1,SEED+2])]:
        for seed in seeds:
            for features,batch in arms:
                name=f'{scope}_{features}_{batch}_{seed}'
                target=RUN/'A03/models'/name
                if (target/'complete.json').exists():continue
                last=0
                while True:
                    candidates=available()
                    if candidates:
                        time.sleep(2)
                        confirmed=available()
                        common=[d for d in candidates if d in confirmed]
                        if common:device=common[0];break
                    if time.time()-last>300:
                        progress('A03','WAITING_GPU',next_model=name,completed_input_preparation=True,
                            free_device_rule='memory.used < 1024 MiB and utilization < 10% twice')
                        last=time.time()
                    time.sleep(30)
                env=os.environ.copy();env['CUDA_VISIBLE_DEVICES']=str(device)
                target.mkdir(parents=True,exist_ok=True)
                command=[sys.executable,str(RUN/'runtime/train_one.py'),'--scope',scope,'--features',features,
                         '--batch',batch,'--seed',str(seed),'--output',str(target)]
                progress('A03','RUNNING',stage='scVI_training',model=name,physical_gpu=device)
                with (target/'training.log').open('ab') as log:
                    subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    finish('A03','model_matrix',pilot_models=4,full_models=12,
        next_steps=['A01/A04 evaluation of all new representations; do not select by visual mixing alone'])


if __name__=='__main__': run()
