"""Parallelism preserves contributors, axes, draws and complete input hashes."""
import json
import os
from pathlib import Path
import tempfile
import time
import unittest

import numpy as np
import pandas as pd
from scipy import sparse
from downstream_runtime import Context
from downstream_parallel import ordered_map
from downstream_bundle import reconstruct_setting
from downstream_profiles import reconstruct_profile_setting
from downstream_intermediate import extract_study
from downstream_measurement import measurement_block
from polyapipe_effects import study
from polyapipe_prepare import save
from pipeline import digest,write,execute,cpu_slots
ROOT=Path(__file__).resolve().parents[2]


def barrier_task(ctx,item):
    (ctx.out/str(item)).write_text('ready');end=time.monotonic()+20
    while not all((ctx.out/str(i)).exists() for i in [0,1]):
        if time.monotonic()>end:raise RuntimeError('Tasks did not overlap')
        time.sleep(.01)
    ctx.register(ctx.config['test_input']);return os.getpid()


def fixture(root):
    for sid in ['GSE_A','GSE_B']:
        p=root/'native/studies'/sid;p.mkdir(parents=True)
        f=pd.DataFrame(dict(atlas_pas_id=['p0','p1','p2'],feature_index_1based=[1,2,3],
            primary_gene_id=['new_gene']*3,primary_gene_symbol=['New']*3,keep_for_apa_primary=[True]*3,
            chrom=['1']*3,strand=['+']*3,cut_site=[100,200,300],matched_source_count_50nt=[0,1,1],species=['human']*3))
        save(f,p/'features.tsv.gz')
        meta=pd.DataFrame(dict(species=['human']*8,study_id=[sid]*8,metacell_id=[f'm{i}' for i in range(8)],
            library_key=['target']*4+['source']*4,route=['R']*8,node=['P']*4+['S']*4,purity=[.85]*8,
            n_cells=[10]*8,matrix_row=np.arange(8)))
        save(meta,p/'metacells.tsv.gz');save(meta,p/'selected_metacells.tsv.gz')
        sparse.save_npz(p/'metacell_by_pas_counts.npz',sparse.csr_matrix(np.array([[0,6,54]]*3+[[0,1,1]]+[[0,54,6]]*3+[[0,1,1]])))
        save(pd.DataFrame([dict(species='human',study_id=sid,route='R',library_key_a='target',library_key_b='source',
            time_a='D16',time_b='D0',sample_unit_a='T',sample_unit_b='S',contrast='P_S',design='true_source')]),p/'endpoint_library_design.tsv')
        save(pd.DataFrame([dict(route='R',library_key=l,sample_unit=u,replicate_status='replicate',time=t)
            for l,u,t in [('target','T','D16'),('source','S','D0')]]),p/'sampling_design.tsv')
        study(p,[('supported50_p08',True,.8,3)])
        g=root/'gex'/sid;g.mkdir(parents=True);gm=meta.copy();gm['total_gex']=1000;save(gm,g/'metacells.tsv.gz')
        save(pd.DataFrame(dict(gene_id=['new_gene'],gene_symbol=['New'])),g/'genes.tsv.gz')
        sparse.save_npz(g/'metacell_by_gene_counts.npz',sparse.csr_matrix(np.array([10,10,10,999,20,20,20,999])[:,None]))
    cycle=root/'cycle.npz';np.savez(cycle,empty=np.array([]))
    return dict(project_root=str(ROOT),native_run=str(root/'native'),frozen_gex=str(root/'gex'),analysis_root=str(root/'analysis'),
        run_id='fixture',input_fingerprint='fixture',resources={'expression_cycle_inputs':dict(path=str(cycle),sha256=digest(cycle))})


def compare(a,b):
    if isinstance(a,pd.DataFrame):pd.testing.assert_frame_equal(a,b,check_exact=True)
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    else:assert a==b


class SingleRouteParallelTests(unittest.TestCase):
    def test_reconstruction_keeps_gene_qualification_and_GEX_contributors_exact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=fixture(root);tasks=[(s,'supported50_p08') for s in ['GSE_A','GSE_B']];values={}
            for workers in [1,2]:
                c=root/f'config{workers}.json';write(c,dict(config,node_workers=dict(state_vectors=workers,profiles=workers,intermediate=workers)))
                for node,fn in [('state_vectors',reconstruct_setting),('profiles',reconstruct_profile_setting),('intermediate',extract_study)]:
                    ctx=Context(c,node,root/f'w{workers}'/node);values[workers,node]=ordered_map(ctx,fn,tasks)
                    self.assertTrue(any(p.endswith('metacell_by_pas_counts.npz') for p in ctx.inputs))
            for node in ['state_vectors','profiles','intermediate']:compare(values[1,node],values[2,node])
            pair=values[2,'state_vectors'][0][0][0].iloc[0]
            self.assertEqual(pair.mean_GEX,15000);self.assertEqual(json.loads(pair.PAS_ids),['p1','p2']);self.assertAlmostEqual(pair.delta_DWUI,.4)

    def test_bootstrap_thinning_and_same_state_draws_are_exact(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=fixture(root);c=root/'vectors.json';write(c,config);ctx=Context(c,'state_vectors',root/'analysis/state_vectors')
            records=ordered_map(ctx,reconstruct_setting,[(s,'supported50_p08') for s in ['GSE_A','GSE_B']])
            for i,n in enumerate(['pair_effects','route_effects','PAS_axes','library_states']):ctx.table(pd.concat([r[i][0] for r in records]),n+'.primary.tsv.gz')
            for workers in [1,2]:
                c=root/f'm{workers}.json';write(c,dict(config,node_workers=dict(measurement=workers)))
                ctx=Context(c,'measurement',root/f'measurement{workers}')
                ordered_map(ctx,measurement_block,[(s,k) for s in ['GSE_A','GSE_B'] for k in ['native','same-state']])
            left=root/'measurement1';right=root/'measurement2';files=list(left.rglob('*.tsv.gz'));self.assertTrue(files)
            for f in files:pd.testing.assert_frame_equal(pd.read_csv(f,sep='\t'),pd.read_csv(right/f.relative_to(left),sep='\t'),check_exact=True)

    def test_workers_overlap_and_merge_input_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input.txt';source.write_text('frozen');c=root/'config.json'
            write(c,dict(project_root=str(ROOT),native_run=str(root),frozen_gex=str(root),analysis_root=str(root),
                run_id='fixture',input_fingerprint='fixture',test_input=str(source),node_workers=dict(test=2)))
            ctx=Context(c,'test',root/'out');self.assertEqual(len(set(ordered_map(ctx,barrier_task,[0,1]))),2)
            self.assertEqual(ctx.inputs[str(source.resolve())],digest(source))

    def test_core_import_rejects_changed_artifacts_and_keeps_source_immutable(self):
        from prepare_parallel_run import prepare
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);package=root/'old';package.mkdir();code=package/'code';code.mkdir();native=root/'native';native.mkdir()
            write(native/'intake.native.json',dict(run_id='r',analysis_route='base'))
            script=code/'core.py';script.write_text("""import argparse,json,pathlib,hashlib
p=argparse.ArgumentParser();p.add_argument('--config');p.add_argument('--out');a=p.parse_args();c=json.loads(pathlib.Path(a.config).read_text());o=pathlib.Path(a.out)
v=o/'value.tsv';v.write_text('value\\n1\\n');h=hashlib.sha256(v.read_bytes()).hexdigest()
(o/'result.json').write_text(json.dumps(dict(node='core',run_id=c['run_id'],input_fingerprint=c['input_fingerprint'],status='complete',checks=dict(current_input=True),artifacts=[dict(path=v.name,sha256=h)])))
""")
            contracts=package/'contracts';contracts.mkdir();(contracts/'table.tsv').write_text('frozen')
            cfg=dict(prepared=True,run_id='r',input_fingerprint='i',caller='APADETECT',reference_variant='v2_primary',analysis_route='base',
                prepared_native_source=str(native),analysis_root=str(package/'statistics'),contracts_root=str(contracts),
                jobs=dict(core=dict(script=str(script),depends_on=[],arguments=[])),code_sha256={str(script):digest(script)},input_files={},resources={})
            c=package/'pipeline.config.json';write(c,cfg);execute(c)
            before={str(p):digest(p) for p in package.rglob('*') if p.is_file()}
            new=prepare(c,root/'parallel',96)
            self.assertEqual(before,{str(p):digest(p) for p in package.rglob('*') if p.is_file()})
            audit=json.loads((new.parent/'core-import.json').read_text())
            self.assertEqual(audit['source_result_sha256'],audit['copied_result_sha256']);self.assertFalse(audit['downstream_results_imported'])
            (package/'statistics/core/value.tsv').write_text('corrupted')
            with self.assertRaisesRegex(ValueError,'Changed node artifact'):prepare(c,root/'bad',96)

    def test_scheduler_includes_nested_workers_in_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);script=root/'node.py'
            script.write_text("""import argparse,json,pathlib,hashlib,time
p=argparse.ArgumentParser();p.add_argument('--config');p.add_argument('--out');a=p.parse_args();c=json.loads(pathlib.Path(a.config).read_text());o=pathlib.Path(a.out)
start=time.time();time.sleep(.15);end=time.time();v=o/'times.json';v.write_text(json.dumps([start,end]));h=hashlib.sha256(v.read_bytes()).hexdigest()
(o/'result.json').write_text(json.dumps(dict(node=o.name,run_id=c['run_id'],input_fingerprint=c['input_fingerprint'],status='complete',checks=dict(intent=True),artifacts=[dict(path=v.name,sha256=h)])))
""")
            cfg=dict(prepared=True,run_id='r',input_fingerprint='i',analysis_root=str(root/'run'),
                jobs={n:dict(script=str(script),depends_on=[],arguments=[]) for n in ['a','b','c']},
                code_sha256={str(script):digest(script)},input_files={},resources={},node_workers=dict(a=2,b=2,c=1))
            self.assertEqual(cpu_slots(cfg,'a'),3)
            self.assertEqual(cpu_slots(dict(cfg,count_workers=32),'count_inference'),33)
            p=root/'config.json';write(p,cfg)
            affinity=os.sched_getaffinity(0) if hasattr(os,'sched_getaffinity') else None
            try:execute(p,workers=3,cpu_budget=4)
            finally:
                if affinity is not None:os.sched_setaffinity(0,affinity)
            periods={n:json.loads((root/'run'/n/'times.json').read_text()) for n in cfg['jobs']}
            self.assertLessEqual(periods['a'][1],periods['b'][0]);self.assertEqual(json.loads((root/'run/pipeline-status.json').read_text())['phase'],'statistics_complete')

if __name__=='__main__':unittest.main()
