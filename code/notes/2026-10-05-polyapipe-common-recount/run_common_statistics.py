"""Own one fresh common-catalogue preparation, core and main statistics run."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def write(path, value):
    tmp=path.with_suffix('.tmp.json');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def main(root, handoff, intake_config, environment, cpu_budget=96, resume_producer_pid=None):
    root=Path(root).resolve(); release=Path(__file__).resolve().parent;native=root/'prepared-native';series=root/'run-series-v2-primary'
    sys.path.insert(0,str(release))
    from polyapipe_prepare import digest
    from prepare_pipeline import prepare
    from pipeline import execute
    h=json.loads(Path(handoff).read_text());inc=json.loads(Path(intake_config).read_text())
    assert digest(handoff)=='7d62ba8a81f7601f82e17461a92d9c26d1c20b3a61a093ab55e38729097f6c70'
    for name in ['real-library-pilot.json','full-framework-preflight.json','annotation-policy-pilot.json']:
        assert json.loads((root/name).read_text())['status']=='PASS', f'Unverified pilot: {name}'
    series.mkdir(exist_ok=resume_producer_pid is not None)
    available=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,set(available[:cpu_budget]))
    assert len(os.sched_getaffinity(0))<=cpu_budget
    cfgpath=series/'runs/polyapipe_common/pipeline.config.json';bundle=cfgpath.parent/'statistics/bundle'
    command=Path('/proc/self/cmdline').read_bytes().decode().split('\0')[:-1]
    controller=dict(pid=os.getpid(),started_at_utc=datetime.now(timezone.utc).isoformat(),command=command,
        start_ticks=Path('/proc/self/stat').read_text().split()[21])
    state=dict(schema_version=1,source_id='polyapipe_common',caller='polyApipe',reference_variant='v2_primary',
        count_input_fingerprint=h['count_input_fingerprint'],count_manifest_sha256=digest(handoff),
        prepared_native=str(native),preparation_receipt=None,source_series=str(series),
        source_runs=str(series/'source-runs.json'),source_config=None,source_bundle=None,
        main_controller=controller,statistics_started=False,core_started=False,cpu_budget=cpu_budget,
        library_aggregation_workers=32,annotation_workers=32,DAG_workers=32,supplement_auto_continuation=False,
        supplement_owner='task_5',original_apadetect_25_routes_modified=False,production_replaced=False)
    def status(phase,**extra):
        state.update(status=phase,current_stage=phase,checked_at_utc=datetime.now(timezone.utc).isoformat(),**extra)
        write(root/'main-handoff-receipt.json',state)
        write(series/'source-runs.json',dict(status=phase,routes=[dict(source_id='polyapipe_common',
            config=state['source_config'],bundle=state['source_bundle'],status=phase)]))
    def run(name,cmd):
        with (root/(name+'.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    try:
        status('VERIFYING_RELEASE')
        for name,sha in json.loads((release/'release-manifest.json').read_text())['files'].items():
            assert digest(release/name)==sha, f'Changed release: {name}'
        for item in [inc['reference'],inc['annotation_code']]:
            assert digest(item['path'])==item['sha256']
        status('PREPARING_COMMON_COUNTS')
        if resume_producer_pid is None:
            exe=release/'polyapipe_common_aggregate'
            run('compile',['g++','-O3','-std=c++17',str(release/'polyapipe_common_aggregate.cpp'),'-lz','-o',str(exe)])
            run('prepare',[sys.executable,'-u',str(release/'polyapipe_common_prepare.py'),'--handoff',str(handoff),
                '--frozen',inc['frozen_native_run'],'--out',str(native),'--executable',str(exe),'--workers','32'])
        else:
            previous=json.loads((root/'launch.json').read_text())
            assert not Path('/proc/'+str(previous['pid'])).exists(), 'Original controller must be stopped first'
            assert not (cfgpath.parent/'statistics').exists(), 'Cannot resume an active or partial scientific DAG'
            producer=Path('/proc/'+str(resume_producer_pid))
            if producer.exists():
                args=producer.joinpath('cmdline').read_bytes().decode().split('\0')
                assert 'polyapipe_common_prepare.py' in args[2] and str(native) in args
            while not (native/'intake.native.json').exists():
                assert producer.exists(), 'Preparation producer exited without a complete intake'
                time.sleep(15)
            prepared=json.loads((native/'intake.native.json').read_text())
            assert prepared['count_input_fingerprint']==h['count_input_fingerprint']
            assert sum(s['molecules'] for s in prepared['studies'])==h['umi_total']
            assert not (native/'annotation_chunks').exists(), 'An annotation process may already own this input'
            state['resumed_counts_without_reaggregation']=True
        status('ANNOTATING_COMMON_PAS')
        code=str(Path(inc['annotation_code']['path']).parents[2])
        run('annotation',[sys.executable,'-u',str(release/'polyapipe_common_annotate.py'),'--root',str(native),
            '--code',code,'--reference',inc['reference']['path'],'--bedtools','/root/nfsdata/conda/envs/apasim/bin/bedtools',
            '--equivalence','/root/nfsdata/apadetect/analysis_projects/reproapa_polyapipe_rerun_20261003_v1/batch_benchmark/benchmark.json',
            '--workers','32'])
        status('BINDING_MAIN_STATISTICS',preparation_receipt=str(native/'intake-complete.json'),
            preparation_receipt_sha256=digest(native/'intake-complete.json'))
        cfg=prepare(environment,native,cfgpath.parent,python=sys.executable,core_workers=7)
        cfg.update(cpu_budget=cpu_budget,count_workers=32,source_id='polyapipe_common',
            execution_revision='common_catalogue_single_route_parallel_v1',node_workers=dict(
                state_vectors=28,profiles=14,intermediate=21,measurement=14,metacell=7,pas_qc=7,coordinate_states=4,count_inference=32))
        # These fields are frozen before any scientific node or receipt exists.
        write(cfgpath,cfg)
        status('MAIN_CONTROLLER_VALIDATING_INPUTS',source_config=str(cfgpath),source_bundle=str(bundle),
            input_fingerprint=cfg['input_fingerprint'])
        with ThreadPoolExecutor(max_workers=1) as pool:
            future=pool.submit(execute,cfgpath,workers=32,cpu_budget=cpu_budget)
            while not future.done():
                started=(cfgpath.parent/'statistics/core/command.log').exists()
                status('CORE_AND_MAIN_STATISTICS_RUNNING' if started else 'MAIN_CONTROLLER_VALIDATING_INPUTS',
                    statistics_started=started,core_started=started)
                time.sleep(15)
            future.result()
        status('MAIN_STATISTICS_COMPLETE',full_main_statistics_complete=True,supplement_statistics_complete=False)
    except Exception as error:
        status('FAILED',error=repr(error),failed_stage=state.get('current_stage'));raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['root','handoff','intake-config','environment']:p.add_argument('--'+name,required=True)
    p.add_argument('--resume-producer-pid',type=int)
    a=p.parse_args();main(a.root,a.handoff,a.intake_config,a.environment,resume_producer_pid=a.resume_producer_pid)
