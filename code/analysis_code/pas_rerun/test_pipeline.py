"""Intent tests: isolation, input invalidation, complete scope and route independence."""
import json
from pathlib import Path
import tempfile
import sys
import unittest
import pandas as pd
from pipeline import digest,write,execute,topological,validate_inputs
from apadetect_prepare import apply_view
from prepare_pipeline import make_jobs
from export_statistics_bundle import panel_contracts
from downstream_axis_extensions import weighted_rank_corr
from downstream_same_pas import covariance_terms,contrast_metrics
from downstream_candidate_screen import low_shift
import numpy as np
ROOT=Path(__file__).resolve().parents[2]

class PipelineTests(unittest.TestCase):
    def test_unbound_plan_cannot_be_run(self):
        with self.assertRaisesRegex(ValueError,'Unbound'):validate_inputs({'prepared':False})

    def test_no_old_PAS_cache_is_allowed(self):
        c=dict(prepared=True,code_sha256={},input_files={},resources={'intermediate_resume':{}})
        with self.assertRaisesRegex(ValueError,'Old PAS cache'):validate_inputs(c)

    def test_entire_manuscript_is_bound_to_executable_acyclic_graph(self):
        here=ROOT/'notes/2026-10-05-apadetect-analysis-intake';env=json.loads((here/'pipeline.environment.json').read_text())
        jobs=make_jobs(env,ROOT/'analysis_code/pas_rerun',sys.executable)
        self.assertEqual(topological(jobs)[-1],'bundle');self.assertNotIn('recount',jobs);self.assertNotIn('caller_comparison',jobs)
        panels,evidence=panel_contracts(here/'contracts',jobs);self.assertEqual((len(panels),len(evidence)),(164,35))
        self.assertTrue({'same_pas','naive_state','shape_ortholog','candidate_screen','geometry'}<=set(jobs))

    def test_parallel_protection_can_restore_sites_removed_in_ASNR_route(self):
        f=pd.DataFrame(dict(quantifiable=[True,True,False],coordinate_supported=[True]*3,
            keep_for_apa_primary=[True]*3,keep_for_sensitivity=[True]*3,keep_for_pas_quant=[True]*3,
            reject_original_point=[True,False,False],reject_original_point_terminal_p50=[False,False,False]))
        base=apply_view(f,'base');asnr=apply_view(base,'reject_original_point');protected=apply_view(base,'reject_original_point_terminal_p50')
        self.assertEqual(asnr.keep_for_apa_primary.tolist(),[False,True,False]);self.assertEqual(protected.keep_for_apa_primary.tolist(),[True,True,False])
        self.assertTrue(base.keep_for_apa_primary.iloc[0]);self.assertEqual(base.index.tolist(),protected.index.tolist())

    def test_direct_paired_metrics_and_covariance_identity(self):
        rng=np.random.default_rng(20);x=rng.normal(size=(30,2));y=x+rng.normal(size=(30,2));w=rng.integers(1,4,30)
        xx=np.column_stack([x,x[:,0]-x[:,1]]);yy=np.column_stack([y,y[:,0]-y[:,1]])
        z=contrast_metrics(xx,yy,w,'PMS');self.assertAlmostEqual(z[3],z[2]-z[0]);self.assertAlmostEqual(z[4],z[2]-z[1])
        c=covariance_terms(x,y,w);self.assertAlmostEqual(c[5:9].sum(),c[9])
        from scipy.stats import spearmanr
        self.assertAlmostEqual(weighted_rank_corr(x[:,0],y[:,0],w),spearmanr(np.repeat(x[:,0],w),np.repeat(y[:,0],w)).statistic)

    def test_new_candidate_identity_not_old_count_or_name(self):
        native=pd.DataFrame([dict(species='mouse',route=r,contrast='P_M',design='same_library',gene_id='new',gene_symbol='NEW',delta_DWUI=.001,MPRO=.2,n_target_libraries=3) for r in ['CIP','WOT']])
        regions=pd.DataFrame(columns=['gene_id','contrast','mask_a','mask_b']);well=pd.DataFrame(columns=['gene_id','metric','q_BH_joint'])
        q=low_shift(native,regions,well);self.assertEqual(q.gene_id.tolist(),['new']);self.assertEqual(q.strict_region_pairs.tolist(),[0])

    def test_verified_resume_and_modified_input_or_result_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);script=root/'node.py';source=root/'input.txt';source.write_text('first')
            script.write_text("""import argparse,json,hashlib,pathlib
p=argparse.ArgumentParser();p.add_argument('--config');p.add_argument('--out');a=p.parse_args();c=json.loads(pathlib.Path(a.config).read_text());o=pathlib.Path(a.out)
x=o/'value.tsv';x.write_text('value\\n1\\n');h=hashlib.sha256(x.read_bytes()).hexdigest()
(o/'result.json').write_text(json.dumps(dict(node=o.name,run_id=c['run_id'],input_fingerprint=c['input_fingerprint'],status='complete',checks=dict(intent=True),artifacts=[dict(path=x.name,sha256=h)])))
""")
            c=dict(prepared=True,run_id='r',input_fingerprint='i',analysis_root=str(root/'run'),jobs={n:dict(script=str(script),arguments=[],depends_on=[] if n=='a' else ['a']) for n in ['a','b']},code_sha256={str(script):digest(script)},input_files={str(source):digest(source)},resources={})
            cfg=root/'config.json';write(cfg,c);execute(cfg);receipt=root/'run/b/pipeline-receipt.json';before=receipt.stat().st_mtime_ns
            execute(cfg);self.assertEqual(before,receipt.stat().st_mtime_ns)
            source.write_text('changed')
            with self.assertRaisesRegex(ValueError,'Input changed'):execute(cfg)
            source.write_text('first');(root/'run/a/value.tsv').write_text('value\n999\n')
            with self.assertRaisesRegex(ValueError,'Changed node artifact'):execute(cfg)

    def test_renderer_consumes_recorded_values_and_rejects_changed_table(self):
        from render_from_bundle import render
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);bundle=root/'bundle';bundle.mkdir()
            table=bundle/'result.tsv';table.write_text('label\testimate\tlo\thi\nA\t-0.4\t-0.7\t0.1\n')
            write(bundle/'table-manifest.json',dict(tables=[dict(path='result.tsv',sha256=digest(table))]))
            write(bundle/'plot-contract.json',dict(panels=[dict(panel='F4D',tables=['result.tsv'])]))
            write(bundle/'bundle-complete.json',dict(status='statistics_bundle_complete',run_id='UNIT_TEST_ONLY',
                table_manifest_sha256=digest(bundle/'table-manifest.json'),plot_contract_sha256=digest(bundle/'plot-contract.json')))
            recipe=root/'recipe.json';write(recipe,dict(panels=[dict(panel='F4D',table='result.tsv',kind='forest',label='label',
                value='estimate',low='lo',high='hi',xlabel='Recorded estimate',title='Unit-test fixture')]))
            before=digest(table);render(bundle,recipe,root/'render')
            receipt=json.loads((root/'render/render-receipt.json').read_text())
            self.assertFalse(receipt['statistics_recomputed']);self.assertEqual(digest(table),before)
            self.assertEqual(receipt['tables']['result.tsv'],before)
            table.write_text(table.read_text().replace('-0.4','0.9'))
            with self.assertRaises(AssertionError):render(bundle,recipe,root/'invalid-render')

if __name__=='__main__':unittest.main()
