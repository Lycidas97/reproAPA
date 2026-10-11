"""Intent: each processing family progresses; existing calculations are not restarted."""
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch
import pandas as pd
from pipeline import write,read,digest
from apadetect_prepare import derive_route
from parallel_processing import derive_route_parallel
from run_processing_families import family_plan,run_group,wait_for_adopted,FAMILIES

ROOT=Path(__file__).resolve().parents[2]


class ProcessingFamilyTests(unittest.TestCase):
    def test_each_family_gets_its_own_slots(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);gate=threading.Barrier(5);seen=[];lock=threading.Lock()
            rows=[dict(route='base',family=FAMILIES[0])]+[dict(route=f'r{f}_{i}',family=f) for f in FAMILIES[1:] for i in range(2)]
            groups=family_plan(dict(analysis_routes=rows),{'base':{}})
            def work(row):
                with lock:seen.append(row['family'])
                gate.wait(timeout=5)
                return dict(source_id=row['route'],bundle=str(root/row['route']))
            with ThreadPoolExecutor(max_workers=3) as pool:
                futures=[pool.submit(run_group,group,min(2,len(group)),work,root/(f+'.json')) for f,group in groups.items()]
                results=[f.result() for f in futures]
            self.assertEqual({f:seen.count(f) for f in FAMILIES},{FAMILIES[0]:1,FAMILIES[1]:2,FAMILIES[2]:2})
            self.assertEqual(sum(len(r['results']) for r in results),5)
            self.assertTrue(all(not r['failures'] for r in results))

    def test_failed_route_does_not_block_other_family_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            def work(row):
                if row['route']=='bad':raise ValueError('invalid source fixture')
                return dict(bundle=str(root/row['route']))
            result=run_group([dict(route=r) for r in ['bad','good','next']],2,work,root/'status.json')
            self.assertEqual(len(result['results']),2);self.assertEqual(len(result['failures']),1)
            self.assertEqual(read(root/'status.json')['phase'],'partial_failure')

    def test_parallel_derivation_preserves_serial_eligibility_and_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'base';base.mkdir();rows=[];outputs={};view='reject_original_point'
            for i in range(3):
                sid='study'+str(i);d=base/'studies'/sid;d.mkdir(parents=True)
                f=pd.DataFrame(dict(atlas_pas_id=['p1','p2','p3'],coordinate_supported=[True,True,True],
                     quantifiable=[True,False,True],keep_for_apa_primary=[True,True,True],
                     keep_for_sensitivity=[True,True,True],keep_for_pas_quant=[True,True,True],
                     matched_source_count_50nt=[1,1,0],reject_original_point=[False,False,True]))
                f.to_csv(d/'features.tsv.gz',sep='\t',index=False)
                (d/'metacell_by_pas_counts.npz').write_bytes(b'immutable count fixture')
                rows.append(dict(study_id=sid,atlas_pas=3,primary_pas=3,supported50_primary_pas=2))
                outputs[sid]={p.name:digest(p) for p in d.iterdir()}
            write(base/'intake.native.json',dict(run_id='fixture__base',analysis_route='base',studies=rows,outputs=outputs))
            write(base/'intake-complete.json',dict(status='READY_FOR_CORE'))
            write(base/'annotation.complete.json',dict(support_reference_variant='v2_primary'))
            write(base/'frozen_native_membership_audit.json',dict(status='fixture'))
            (base/'frozen_native_membership_audit.tsv').write_text('fixture\n')
            before={str(p):digest(p) for p in base.rglob('*') if p.is_file()}
            derive_route(base,view,root/'serial');derive_route_parallel(base,view,root/'parallel',2)
            a=read(root/'serial/intake.native.json');b=read(root/'parallel/intake.native.json')
            self.assertEqual(a['studies'],b['studies']);self.assertEqual(a['run_id'],b['run_id'])
            for row in rows:
                rel=Path('studies')/row['study_id']
                pd.testing.assert_frame_equal(pd.read_csv(root/'serial'/rel/'features.tsv.gz',sep='\t'),pd.read_csv(root/'parallel'/rel/'features.tsv.gz',sep='\t'))
                self.assertEqual((root/'parallel'/rel/'metacell_by_pas_counts.npz').resolve(),(base/rel/'metacell_by_pas_counts.npz').resolve())
            self.assertTrue(all(digest(Path(p))==h for p,h in before.items()))

    def test_core_concurrency_is_configurable_without_changing_default(self):
        from prepare_pipeline import prepare
        env=ROOT/'notes/2026-10-05-apadetect-analysis-intake/pipeline.environment.json'
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            cfg=prepare(env,root/'input',root/'parallel',plan_only=True,core_workers=7)
            self.assertEqual(cfg['core_workers'],7)
            default=prepare(env,root/'input',root/'default',plan_only=True)
            self.assertEqual(default['core_workers'],2)
            with self.assertRaisesRegex(ValueError,'positive integer'):prepare(env,root/'input',root/'invalid',plan_only=True,core_workers=0)
            self.assertFalse((root/'invalid').exists())

    def test_sigint_drains_active_v9_work_without_dispatching_pending_route(self):
        # SIGINT to the controller PID (not its process group) unwinds the main
        # executor into shutdown(wait=True), leaving its active workers running.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);script=root/'controller.py'
            script.write_text('''from pathlib import Path
import sys,time
sys.path.insert(0,CODE)
import run_sources
root=Path(sys.argv[1])
def execute(config,workers):
 name=Path(config).stem;(root/(name+'.started')).touch()
 while not (root/(name+'.release')).exists():time.sleep(.02)
 (root/(name+'.finished')).touch()
run_sources.execute=execute
runs=[dict(source_id=n,config=str(root/(n+'.json')),bundle=str(root/n),identity={}) for n in ['a','b','pending']]
run_sources.execute_sources(runs,root,1,2)
'''.replace('CODE',repr(str(Path(__file__).parent))))
            with (root/'log').open('w') as log:
                proc=subprocess.Popen([sys.executable,str(script),str(root)],stdout=log,stderr=log,start_new_session=True)
                try:
                    until=time.monotonic()+10
                    while not all((root/(n+'.started')).exists() for n in ['a','b']):
                        self.assertIsNone(proc.poll());self.assertLess(time.monotonic(),until);time.sleep(.02)
                    os.kill(proc.pid,signal.SIGINT);time.sleep(.15)
                    (root/'a.release').touch();time.sleep(.15)
                    self.assertFalse((root/'pending.started').exists());self.assertIsNone(proc.poll())
                    (root/'b.release').touch();proc.wait(timeout=5)
                    self.assertNotEqual(proc.returncode,0)
                    self.assertTrue(all((root/(n+'.finished')).exists() for n in ['a','b']))
                    self.assertFalse((root/'pending.started').exists())
                finally:
                    if proc.poll() is None:proc.kill();proc.wait()

    def test_adopted_failure_is_reported_without_recomputation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);write(root/'config.json',dict(analysis_route='base',analysis_root=str(root/'stats')))
            write(root/'stats/pipeline-status.json',dict(phase='failed',error='fixture failed'))
            with patch('run_processing_families.execute') as compute:
                with self.assertRaisesRegex(ValueError,'Adopted route failed'):
                    wait_for_adopted('base',dict(config=str(root/'config.json'),bundle=str(root/'bundle'),source_stage=str(root/'stage.json')),{},root/'outstage.json',.01)
                compute.assert_not_called()

    def test_adoption_waits_for_final_DAG_acceptance_even_when_bundle_exists(self):
        from test_result_records import bundle
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);result,identity,_=bundle(root,'done',1.)
            cfg=dict(identity,analysis_root=str(root/'stats'))
            write(root/'config.json',cfg);write(root/'stats/pipeline-status.json',dict(phase='running',completed=[]))
            controller=dict(pid=1,start_ticks='fixture',argv=['fixture'])
            def finish(_):write(root/'stats/pipeline-status.json',dict(phase='statistics_complete'))
            with patch('run_processing_families.process_identity',return_value=controller),patch('run_processing_families.time.sleep',side_effect=finish) as wait:
                record=wait_for_adopted('base',dict(config=str(root/'config.json'),bundle=str(result),source_stage=str(root/'stage.json')),controller,root/'outstage.json',.01)
                wait.assert_called_once()
            self.assertTrue(record['adopted']);self.assertEqual(record['identity'],identity)


if __name__=='__main__':unittest.main()
