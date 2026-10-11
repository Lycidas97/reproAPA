import csv
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

from aggregate import aggregate_counts
from intake import digest, load_matrix, validate
from run import audit_baseline, binding_errors, execute, load_workflow, verify_result
from stage import check_package, stage, read_tsv


def write_table(path, rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


class RerunTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.x=sparse.csr_matrix(np.array([[3,2,0],[0,0,0],[1,0,4],[99,20,0]],dtype=np.int64))
        sparse.save_npz(self.root/'counts.npz',self.x)
        self.cells=[dict(cell_index_1based=i+1,original_barcode=b,cell_key='L:'+b) for i,b in enumerate('abcd')]
        self.members=[dict(cell_id='L:'+b,barcode=b,metacell_id=mc,library_key='L') for b,mc in [('c','m1'),('a','m1'),('b','m2')]]
        self.metacells=[dict(metacell_id='m2',n_cells=1),dict(metacell_id='m1',n_cells=2)]
        write_table(self.root/'cells.tsv',self.cells)
        write_table(self.root/'membership.tsv',self.members)
        write_table(self.root/'metacells.tsv',self.metacells)
        write_table(self.root/'features.tsv',[dict(atlas_pas_id='p'+str(i),feature_index_1based=i,chrom='chr1',cut_site=100*i,strand='-',primary_gene_id='g1',keep_for_apa_primary='True') for i in range(1,4)])
        write_table(self.root/'status.tsv',[dict(atlas_pas_id='p'+str(i),terminal_status='observed',projected_umi=n) for i,n in enumerate([103,22,4],1)])
        write_table(self.root/'scope.tsv',[dict(study='GSE_TEST',library_id='L',species='human',expression_atlas_cells=3)])
        (self.root/'resource.json').write_text('{}\n')
        self.m=dict(schema_version=1,preserve_expression_framework=True,pas_release_id='synthetic-test-only',caller='synthetic',quantification=dict(version='test',command='fixture',deduplication_unit='PAS-cell-UMI'),cut_coordinate_system='0-based-genomic-base',
                    studies=[dict(study_id='GSE_TEST',assembly='GRCh38',features=self.desc('features.tsv'))],
                    libraries=[dict(library_id='L',library_key='L',study_id='GSE_TEST',orientation='cell_by_pas',counts=self.desc('counts.npz'),cells=self.desc('cells.tsv'),membership=self.desc('membership.tsv'),metacells=self.desc('metacells.tsv'),feature_status=self.desc('status.tsv'))],
                    resources={k:[self.desc('resource.json')] for k in ['gex','state_annotations','embeddings','library_design','contrast_design','genome_fasta','gene_annotation','external_pas','orthologs','sequence_resources','recount_inputs','external_validation','legacy_rosters']})
        self.manifest=self.root/'intake.json';self.save()

    def tearDown(self):
        self.tmp.cleanup()

    def desc(self,name):
        p=self.root/name
        return dict(path=str(p),sha256=digest(p))

    def save(self):
        self.manifest.write_text(json.dumps(self.m))

    def validate(self):
        self.save();return validate(self.manifest,scope_path=self.root/'scope.tsv')

    def test_valid_input_keeps_zero_cells_and_extra_called_cells(self):
        result=self.validate()
        self.assertEqual(result['frozen_cells'],3)
        self.assertEqual(result['libraries'][0]['cells'],4)

    def test_aggregation_uses_barcode_not_old_pas_row_and_excludes_nonatlas_cells(self):
        got=aggregate_counts(self.x,self.cells,self.members,self.metacells)
        np.testing.assert_array_equal(got.toarray(),[[0,0,0],[4,2,4]])

    def test_matrix_market_orientation_matches_npz(self):
        path=self.root/'feature_by_cell.mtx';mmwrite(path,self.x.T)
        self.assertEqual((load_matrix(path,'pas_by_cell',4,3)!=self.x).nnz,0)
        with self.assertRaisesRegex(ValueError,'dimension mismatch'):
            load_matrix(path,'cell_by_pas',4,3)

    def test_noninteger_count_is_not_accepted_as_umi(self):
        sparse.save_npz(self.root/'fraction.npz',self.x.astype(float)/2)
        with self.assertRaisesRegex(ValueError,'Non-integer'):
            load_matrix(self.root/'fraction.npz','cell_by_pas',4,3)

    def test_missing_library_cannot_silently_shrink_scope(self):
        self.m['libraries']=[]
        with self.assertRaisesRegex(ValueError,'scope differs'):
            self.validate()

    def test_missing_zero_cell_is_not_an_expression_membership_change(self):
        write_table(self.root/'cells.tsv',[r for r in self.cells if r['original_barcode']!='b'])
        for i,r in enumerate([r for r in self.cells if r['original_barcode']!='b'],1):
            r['cell_index_1based']=i
        write_table(self.root/'cells.tsv',[r for r in self.cells if r['original_barcode']!='b'])
        self.m['libraries'][0]['cells']=self.desc('cells.tsv')
        with self.assertRaisesRegex(ValueError,'Frozen cells missing'):
            self.validate()

    def test_projection_zero_cannot_hide_counts(self):
        rows=[dict(atlas_pas_id='p'+str(i),terminal_status='not_observed',projected_umi=n) for i,n in enumerate([103,22,4],1)]
        write_table(self.root/'status.tsv',rows);self.m['libraries'][0]['feature_status']=self.desc('status.tsv')
        with self.assertRaisesRegex(ValueError,'absent PAS has nonzero'):
            self.validate()

    def test_changed_input_hash_fails(self):
        (self.root/'resource.json').write_text('{"changed":true}')
        with self.assertRaisesRegex(ValueError,'hash mismatch'):
            self.validate()

    def test_counts_must_match_status_totals(self):
        self.x[0,0]=4;sparse.save_npz(self.root/'counts.npz',self.x)
        self.m['libraries'][0]['counts']=self.desc('counts.npz')
        with self.assertRaisesRegex(ValueError,'count conservation'):
            self.validate()

    def adapter(self,fail=False,stale=False):
        path=self.root/'adapter.py'
        path.write_text('''import hashlib,json,os
from pathlib import Path
p=Path(os.environ['PAS_NODE_DIR'])
(p/'data.tsv').write_text('value\\n1\\n')
(p/'audit.json').write_text('{"fixture":true}')
r=dict(node=os.environ['PAS_NODE_ID'],run_id=os.environ['PAS_RUN_ID'],input_fingerprint=os.environ['PAS_INPUT_FINGERPRINT'],status='complete',checks={'fixture':True},artifacts=[{'path':name,'sha256':hashlib.sha256((p/name).read_bytes()).hexdigest(),'role':role} for name,role in [('data.tsv','source_data'),('audit.json','audit')]])
'''+("r['run_id']='OLD_RUN'\n" if stale else '')+"(p/'result.json').write_text(json.dumps(r))\n"+("raise SystemExit(9)\n" if fail else ''))
        return dict(argv=[sys.executable,str(path)],code_files=[self.desc('adapter.py')])

    def test_full_binding_is_required(self):
        w={'nodes':[dict(id='first',depends_on=[]),dict(id='second',depends_on=['first'])]}
        problems=binding_errors(w,{'nodes':{'first':self.adapter()}},self.root)
        self.assertTrue(any('second: no adapted command' in p for p in problems))

    def test_failed_node_stops_dependents_and_cannot_mark_complete(self):
        w={'nodes':[dict(id='first',depends_on=[]),dict(id='second',depends_on=['first'])]}
        b={'nodes':{'first':self.adapter(fail=True),'second':dict(argv=['SHOULD_NOT_RUN'])}}
        with self.assertRaises(Exception):
            execute(w,b,self.manifest,self.root/'run',scope_path=self.root/'scope.tsv')
        state=json.loads((self.root/'run/run.json').read_text())
        self.assertEqual(state['status'],'failed');self.assertFalse((self.root/'run/second').exists())

    def test_old_run_receipt_is_rejected(self):
        w={'nodes':[dict(id='first',depends_on=[])]};b={'nodes':{'first':self.adapter(stale=True)}}
        with self.assertRaisesRegex(ValueError,'different node/run/inputs'):
            execute(w,b,self.manifest,self.root/'run',scope_path=self.root/'scope.tsv')

    def test_execution_success_and_no_stale_directory_reuse(self):
        w={'nodes':[dict(id='first',depends_on=[]),dict(id='second',depends_on=['first'])]}
        a=self.adapter();b={'nodes':{'first':a,'second':a}}
        result=execute(w,b,self.manifest,self.root/'run',scope_path=self.root/'scope.tsv')
        self.assertEqual(result['status'],'complete');self.assertEqual(result['completed'],['first','second'])
        with self.assertRaisesRegex(ValueError,'already exists'):
            execute(w,b,self.manifest,self.root/'run',scope_path=self.root/'scope.tsv')

    def test_production_graph_covers_every_registered_analysis(self):
        w=load_workflow(Path(__file__).with_name('workflow.json'))
        ids={r['id'] for r in w['nodes']}
        self.assertEqual(len(ids),49)
        self.assertIn('caller_comparison',ids)
        self.assertEqual(set(next(r for r in w['nodes'] if r['id']=='render')['depends_on']),ids-{'render','tables','text','release'})

    def make_package(self):
        home=Path(__file__).parent;out=self.root/'delivery_run';out.mkdir()
        nodes=[n['id'] for n in load_workflow(home/'workflow.json')['nodes'] if n['id']!='release']
        (out/'run.json').write_text(json.dumps(dict(status='running',active_node='release',completed=nodes,run_id='fixture',input_fingerprint='fixture')))
        files=[]
        for node in nodes:
            folder=out/node;folder.mkdir()
            (folder/'data.tsv').write_text('fixture\n1\n');(folder/'audit.json').write_text('{}')
            artifacts=[dict(path=name,sha256=digest(folder/name),role=role) for name,role in [('data.tsv','source_data'),('audit.json','audit')]]
            def add(name,content,kind=None,target=None,**kw):
                path=folder/name;path.write_text(content)
                artifacts.append(dict(path=name,sha256=digest(path),role='source_data'))
                descriptor=dict(path=str(path),sha256=digest(path))
                if kind:
                    files.append(dict(**descriptor,kind=kind,target=target,**kw))
                return descriptor
            if node=='render':
                for row in read_tsv(home/'figures.tsv'):
                    if row['action']=='rebuild':
                        add(Path(row['target']).name,'%PDF-1.7\nSYNTHETIC, NOT A RENDERED FIGURE\n',kind='figure',target=row['target'])
            if node=='tables':
                for i in range(8,31):
                    add(f'S{i}.tsv','fixture\n1\n',kind='table',target=f'table/S{i}.tsv',table_ids=[i])
            if node=='text':
                rows=read_tsv(home/'text-targets.tsv')
                review=add('review.tsv','id\tbaseline_sha256\tstatus\tevidence\n'+''.join(f'{r["id"]}\t{r["sha256"]}\tverified_unchanged\tsynthetic test only\n' for r in rows))
            (folder/'result.json').write_text(json.dumps(dict(node=node,run_id='fixture',input_fingerprint='fixture',status='complete',checks={'synthetic_fixture':True},artifacts=artifacts)))
        package=self.root/'package.json';package.write_text(json.dumps(dict(run_id='fixture',input_fingerprint='fixture',files=files,text_review=review)))
        return package,out

    def test_package_refuses_missing_figure_even_when_run_nodes_completed(self):
        package,out=self.make_package();data=json.loads(package.read_text());data['files'].pop(0);package.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'figure replacements'):
            check_package(package,out)

    def test_package_refuses_incomplete_text_review(self):
        package,out=self.make_package();data=json.loads(package.read_text())
        review=Path(data['text_review']['path']);review.write_text(review.read_text().replace('verified_unchanged','pending',1))
        data['text_review']['sha256']=digest(review);package.write_text(json.dumps(data))
        # The changed artifact is rejected even before the unresolved claim is examined.
        with self.assertRaisesRegex(ValueError,'hash mismatch'):
            check_package(package,out)

    def test_staging_uses_actual_filename_not_display_number(self):
        package,out=self.make_package()
        # Baseline mutation is tested separately; this fixture exercises output naming.
        with patch('stage.audit_baseline',return_value=0):candidate=stage(package,out,self.root/'candidate')
        rows=read_tsv(candidate/'Fig/figure-manifest.tsv')
        row=next(r for r in rows if r['figure']=='Figure S6')
        self.assertEqual(row['file'],'SuppFig.S5.pdf')
        self.assertEqual(row['sha256'],digest(candidate/'Fig'/row['file']))
        self.assertEqual(json.loads((candidate/'PAS-RERUN.json').read_text())['status'],'staged_not_compiled')
        self.assertFalse((candidate/'table/Supplementary-Tables-S1-S15.xlsx').exists())

    def test_baseline_guard_rejects_changed_source_even_if_outputs_are_valid(self):
        source=self.root/'source.py';source.write_text('version one\n')
        write_table(self.root/'baseline-files.tsv',[dict(file='source.py',sha256=digest(source))])
        with patch('run.HERE',self.root),patch('run.ROOT',self.root):
            self.assertEqual(audit_baseline(),1)
            source.write_text('version two\n')
            with self.assertRaisesRegex(ValueError,'Baseline changed'):audit_baseline()


if __name__=='__main__':
    unittest.main()
