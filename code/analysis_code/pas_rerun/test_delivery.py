import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from downstream_runtime import Context
from downstream_tables import verified_sources,workbook_input
from polyapipe_prepare import digest


class DeliveryIntentTests(unittest.TestCase):
    def test_candidate_release_rejects_missing_pages_and_changed_pdf(self):
        import pymupdf
        from finalize_release import check_pages
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);folder=root/'output/pdf';folder.mkdir(parents=True)
            image=root/'page.png';image.write_bytes(b'reviewed-page-fixture')
            review=dict(documents={})
            for name in ['main','supplement']:
                document=pymupdf.open();document.new_page();document.new_page()
                pdf=folder/f'{name}.pdf';document.save(pdf);document.close()
                review['documents'][name]=dict(sha256=digest(pdf),pages=[dict(page=i,layout='pass',notes='Synthetic validation fixture',image=str(image),image_sha256=digest(image)) for i in [1,2]])
            check_pages(root,review,lambda p:p)
            page=review['documents']['main']['pages'].pop()
            with self.assertRaises(AssertionError):check_pages(root,review,lambda p:p)
            review['documents']['main']['pages'].append(page)
            with (folder/'main.pdf').open('ab') as f:f.write(b'\n% changed after review\n')
            with self.assertRaisesRegex(AssertionError,'PDF changed'):check_pages(root,review,lambda p:p)

    def test_local_state_exclusion_resolves_namespace_and_keeps_metacell_order(self):
        from downstream_runtime import local_state_labels
        frame=pd.DataFrame(dict(metacell_id=['a','b'],local_state_id=['WOT::common_to_2i::L00','C012']))
        self.assertEqual(local_state_labels(frame,['b','a']).tolist(),['C012','L00'])
        with self.assertRaises(AssertionError):local_state_labels(frame,['missing'])

    def test_named_seed_query_uses_the_seed_identity_not_a_mirna_number(self):
        from downstream_example_inclusion import named_sites
        f=pd.DataFrame(dict(category=['miRNA','miRNA','RBP'],factor=['AAGTGCT','TCACAAA','Pcbp1'],
            names=['mmu-miR-291a-3p;mmu-miR-302a-3p','mmu-miR-23a-3p','Pcbp1']))
        self.assertEqual(named_sites(f,'RAD23B').tolist(),[True,False,False])
        self.assertEqual(named_sites(f,'HMGN2').tolist(),[False,True,False])
        self.assertEqual(named_sites(f,'EIF4EBP2').tolist(),[False,False,True])

    def test_paired_inclusion_requires_support_in_the_common_libraries(self):
        from downstream_regulatory import observed_regions
        f=pd.DataFrame(dict(contrast=['MP','SP'],observed_both5=[True,True],
            paired_observed_both5=[False,True],effect20=[.2,.2],effect50=[.2,.2],
            separate_effect20=[.3,.3],separate_effect50=[.3,.3],strict_change=[True,True]))
        paired=observed_regions(f)
        self.assertEqual(paired.observed_both5.tolist(),[False,True])
        self.assertEqual(paired.strict_change.tolist(),[False,True])
        self.assertEqual(paired.effect20.tolist(),[.2,.2])
        separate=observed_regions(f,common=False)
        self.assertEqual(separate.strict_change.tolist(),[True,True])
        self.assertEqual(separate.effect20.tolist(),[.3,.3])
        self.assertEqual(f.strict_change.tolist(),[True,True])

    def test_endpoint_site_union_excludes_other_contrasts_and_duplicate_factors(self):
        from render_tables_tex import endpoint_site_union
        sites=pd.DataFrame([dict(species='mouse',category='miRNA',gene_id=gene,
            chrom='chr1',strand='+',start0=start,end0=start+7,contrast=contrast,factor=factor)
            for gene,start,contrast,factor in [('g1',10,'SP','a'),('g1',10,'MP','a'),
                ('g1',10,'SP','b'),('g1',20,'MP','c'),('g2',30,'RM','a')]])
        summary=endpoint_site_union(sites)
        self.assertEqual(summary.unique_physical_sites.tolist(),[2])
        self.assertEqual(summary.genes.tolist(),[1])
        self.assertEqual(len(sites),5)

    def test_unavailable_planned_null_tests_remain_inside_plot(self):
        import numpy as np
        import matplotlib.pyplot as plt
        from render_recurrence import null_forest,grouped_null
        f=pd.DataFrame(dict(observed=[2,3,4],null_mean=[np.nan,2,np.nan],
            null_p025=[np.nan,1,np.nan],null_p975=[np.nan,3,np.nan],
            q_metric_8=[np.nan,.2,np.nan],species=['human']*3,contrast=['P_S']*3,model=['a','b','c']))
        for draw in [lambda ax:null_forest(ax,f,['a','b','c']),
                     lambda ax:grouped_null(ax,f,'model','q_metric_8')]:
            fig,ax=plt.subplots();draw(ax);fig.canvas.draw()
            low,high=sorted(ax.get_ylim())
            self.assertLess(low,-.285);self.assertGreater(high,2)
            unavailable=sum(t.get_text()=='Unavailable' for t in ax.texts)+sum(line.get_marker()=='x' for line in ax.lines)
            self.assertEqual(unavailable,2)
            plt.close(fig)

    def test_frozen_embedding_resolves_WOT_alias_without_losing_cells(self):
        from render_common import join_frozen_scores
        cells=pd.DataFrame(dict(cell_id=['w','g'],study_id=['WOT_GSE115943','GSE103221'],metacell_id=['m','m']))
        scores=pd.DataFrame(dict(study_id=['GSE103221','GSE115943'],metacell_id=['m','m'],value=[.1,.2]))
        joined=join_frozen_scores(cells,scores)
        self.assertEqual(joined.cell_id.tolist(),['w','g'])
        self.assertEqual(joined.study_id.tolist(),cells.study_id.tolist())
        self.assertEqual(joined.value.tolist(),[.2,.1])
        with self.assertRaises(AssertionError):join_frozen_scores(cells,scores.iloc[:1])

    def test_control_archive_requires_every_statistic_including_range_sensitivity(self):
        from archive_reference_control import archive
        import zipfile
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'config.json';workflow=root/'workflow.json'
            config.write_text(json.dumps(dict(reference_variant='v3_ge62_comparison',analysis_root=str(root),run_id='control',input_fingerprint='control-inputs')))
            workflow.write_text(json.dumps(dict(nodes=[dict(id='recurrence')])) )
            for node in ['recurrence','figure_scores','cross_sharing_range10']:
                folder=root/node;folder.mkdir();items=[]
                for name,role in [('values.tsv','source_data'),('audit.json','audit')]:
                    path=folder/name;path.write_text('fixture\n');items.append(dict(path=name,sha256=digest(path),role=role))
                (folder/'result.json').write_text(json.dumps(dict(node=node,run_id='control',input_fingerprint='control-inputs',status='complete',checks=dict(verified=True),artifacts=items)))
            receipt=root/'cross_sharing_range10/result.json';saved=receipt.read_text();receipt.unlink()
            out=root/'control.zip'
            with self.assertRaises(FileNotFoundError):archive(config,workflow,out)
            self.assertFalse(out.exists())
            receipt.write_text(saved)
            c=json.loads(config.read_text());c['extra_control_nodes']=['shared_projection'];config.write_text(json.dumps(c))
            with self.assertRaises(FileNotFoundError):archive(config,workflow,out)
            self.assertFalse(out.exists())
            c['extra_control_nodes']=[];config.write_text(json.dumps(c));archive(config,workflow,out)
            with zipfile.ZipFile(out) as z:
                manifest=json.loads(z.read('MANIFEST.json'))
                self.assertEqual(set(manifest['completed_nodes']),{'recurrence','figure_scores','cross_sharing_range10'})
                self.assertEqual(manifest['run_id'],'control')

    def test_text_release_rejects_pending_claims_and_changed_frozen_prose(self):
        from finalize_text import check_review
        with tempfile.TemporaryDirectory() as tmp:
            baseline=Path(tmp)/'baseline';split=Path(tmp)/'split'
            baseline.mkdir();split.mkdir()
            (baseline/'main.tex').write_text('Frozen GEX result.\n')
            target=dict(id='T0001',file='main.tex',sha256='baseline',text='Frozen GEX result.')
            row=dict(id='T0001',file='main.tex',baseline_sha256='baseline',status='pending',evidence='Awaiting results')
            with self.assertRaises(AssertionError):check_review([row],[target],split,baseline)
            row.update(status='verified_unchanged',evidence='Frozen expression framework retained')
            check_review([row],[target],split,baseline)
            (split/'main.tex').write_text('Altered GEX result.\n')
            with self.assertRaises(AssertionError):check_review([row],[target],split,baseline)

    def test_small_probabilities_preserve_significant_digits_before_excel_export(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            context=SimpleNamespace(out=Path(tmp))
            values=[1.234567890123456e-12,3.234567890123456e-100,0.,None]
            path=workbook_input(context,pd.DataFrame({'p':values}),'S21')
            result=json.loads((Path(tmp)/path).read_text())
            self.assertEqual([row[0] for row in result['data']],values)

    def test_workbook_validation_cannot_accept_missing_as_zero_or_tiny_p_as_zero(self):
        from finalize_tables import same_value
        self.assertFalse(same_value(0.,None))
        self.assertFalse(same_value(0.,1e-100))
        self.assertFalse(same_value(1.000000000000001,1))
        self.assertTrue(same_value(None,None))
        self.assertTrue(same_value(0.,0.))

    def test_source_connection_uses_designated_studies_and_original_PM_direction(self):
        from render_recurrence import matched_flow
        roster=pd.DataFrame([dict(species=sp,gene_id=g,metric='WPP',threshold=.02,contrast='P_M',
            series_matching=True,direction='negative') for sp,g in [('human','h'),('mouse','m')]])
        vectors=pd.DataFrame([dict(species=sp,gene_id=g,setting='supported50_p08',study_id=sid,route=route,
            DWUI_PM=-.1,DWUI_PS=ps) for sp,g,sid,route,ps in [('human','h','Xing','GSE118258',-.1),
            ('mouse','m','Guo','CIP',.1),('mouse','m','WOT','WOT',.1)]])
        flow,detail=matched_flow(roster,vectors);d=detail.set_index('species')
        self.assertTrue(d.loc['human','complete']);self.assertTrue(d.loc['human','additional_PS'])
        self.assertTrue(d.loc['mouse','retained_PM']);self.assertFalse(d.loc['mouse','additional_PS'])
        self.assertEqual(d.loc['human','required_studies'],1)
        self.assertEqual(d.loc['mouse','required_studies'],2)

    def test_leaving_a_library_rechecks_joint_PAS_support_before_reporting_influence(self):
        from unittest.mock import patch
        import numpy as np
        from downstream_example_inclusion import fixed_source_influence
        rows=[dict(study_id='study',gene_id='gene',gene_symbol='Eif4ebp2',setting='supported50_p08',route='CIP',
            library_key=lib,state=state,time='D0' if state=='S' else 'D4',sample_unit=lib,value=value)
            for lib,state,value in [('source','S',.1),('P1','P',.2),('P2','P',.3)]]
        class Capture:
            analysis=Path('/fixture')
            def __init__(self):self.tables={}
            def read(self,path):return pd.DataFrame(rows)
            def table(self,frame,name,**kwargs):self.tables[name]=pd.DataFrame(frame)
        class FakeStudy:
            def __init__(self,*args):self.pairs=pd.DataFrame([dict(route='CIP',design='true_source',library_key_b='source')])
            def axis(self,gene):return pd.DataFrame(dict(position=[0.,1.]))
            def gene_counts(self,axis):return None
            def state(self,data,route,lib,state):
                return dict(endpoint_eligible=True,site_counts=np.array([20,4 if lib=='source' else 3]),site_mc=np.array([2,1]))
        ctx=Capture()
        with patch('downstream_example_inclusion.Study',FakeStudy):fixed_source_influence(ctx)
        f=ctx.tables['fixed_source_DWUI_influence.tsv'];f=f[(f.omission_kind=='library')&(f.state=='P')].set_index('omitted')
        self.assertEqual(f.loc['none','status'],'measured')
        self.assertEqual(f.loc['P2','n_joint_pas'],1)
        self.assertEqual(f.loc['P2','status'],'joint_PAS_support_failed')
        self.assertTrue(pd.isna(f.loc['P2','delta_DWUI']))

    def test_cis_descriptive_controls_retain_reused_gene_weights_and_unmatched_model_genes(self):
        from downstream_cis import descriptive_controls
        from downstream_cis_models import CIS
        class Capture:
            analysis=Path('/fixture')
            def __init__(self,matches):self.matches=matches;self.tables={}
            def read(self,path):return self.matches.copy()
            def table(self,frame,name):self.tables[name]=pd.DataFrame(frame)
        identity=dict(layer='supported50',contrast='P-S',discovery='human',tau=.02)
        data=pd.DataFrame([dict(identity,human_gene_id=g,Y=y,**{f:v for f in CIS})
            for g,y,v in [('case1',1,1),('case2',0,2),('control',1,3),('unmatched',0,4)]])
        matched=pd.DataFrame([dict(identity,human_gene_id=g,match_id=k,role=role,match_weight=w,validation_same_direction=y)
            for g,k,role,w,y in [('case1',1,'recurrent',1,True),('control',1,'control',.25,True),
                ('case2',2,'recurrent',1,False),('control',2,'control',.75,True)]])
        ctx=Capture(matched);descriptive_controls(ctx,data)
        counts=ctx.tables['descriptive_class_counts.tsv'].iloc[0]
        self.assertEqual(counts.model_n,4)
        self.assertEqual(counts['Matched background'],1)
        distributions=ctx.tables['feature_distribution_summary.tsv']
        control=distributions[(distributions.contrast=='P-S')&(distributions.discovery=='human')&
            (distributions['class']=='Matched background')&(distributions.feature==CIS[0])].iloc[0]
        self.assertEqual(control.n_observed,1);self.assertEqual(control.weight_sum,1.)
        self.assertEqual(control['median'],3.)

    def test_separate_figure_runs_preserve_prior_verified_outputs_and_exact_panel_arrays(self):
        from render_common import Renderer
        import numpy as np
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);config=root/'config.json'
            config.write_text(json.dumps(dict(analysis_root=str(root),run_id='fixture',input_fingerprint='fixture',resources={})))
            for number in [1,2]:
                r=Renderer(config,root/'render');r.start('Figure '+str(number),60,40);r.label('A',2,2)
                ax=r.axes((12,12,40,20));ax.plot([0,1],[np.nan,-.00003]);ax.scatter([0,1],[0,1])
                r.save('Fig/Figure.'+str(number)+'.pdf','Synthetic data: registry/precision test only')
            progress=json.loads((root/'render/render-progress.json').read_text())
            self.assertEqual(len(progress['figures']),2)
            self.assertFalse(any(a['visual_review_complete'] for a in progress['figures']))
            with np.load(root/'render/panel-data/Figure_1/A.npz') as arrays:
                self.assertTrue(np.isnan(arrays['axes0_line0_y'][0]))
                self.assertEqual(arrays['axes0_line0_y'][1],-.00003)

    def test_interval_pairing_uses_common_physical_libraries_and_preserves_singleton_uncertainty(self):
        from downstream_example_inclusion import feature_libraries,mean_summary
        def row(lib,counts):return dict(library=lib,counts=counts,total=sum(counts),coverage=1.)
        model=dict(profiles={'M':dict(libraries=[row('shared',[16,4]),row('M_only',[0,40])]),
            'P':dict(libraries=[row('shared',[4,16]),row('P_only',[40,0])])})
        libs,pairs=feature_libraries(model,[1,0],dict(shared='d1',M_only='d1',P_only='d2'),20)
        self.assertEqual(pairs.library_key.tolist(),['shared'])
        self.assertAlmostEqual(pairs.delta.iloc[0],-.6)
        summary=mean_summary(pairs,'delta')
        self.assertTrue(pd.isna(summary['leave_one_min']))
        _,deep=feature_libraries(model,[1,0],dict(shared='d1',M_only='d1',P_only='d2'),50)
        self.assertTrue(deep.empty)

    def test_matched_reference_flow_keeps_unmeasured_candidates_and_requires_two_studies(self):
        from render_recurrence import matched_flow
        members=pd.DataFrame(dict(metric=['WPP']*3,threshold=[.02]*3,series_matching=[True]*3,
            contrast=['P_M']*3,species=['mouse']*3,gene_id=['a','b','c'],direction=['negative']*3))
        rows=[]
        for gene,values in [('a',[-.1,-.1]),('b',[-.1]),('c',[-.1,.1])]:
            for i,value in enumerate(values):rows.append(dict(species='mouse',gene_id=gene,setting='supported50_p08',study_id=f'G{i}',route=['CIP','WOT'][i],DWUI_PM=value,DWUI_PS=-.2))
        flow,detail=matched_flow(members,pd.DataFrame(rows))
        self.assertEqual(flow[flow.species.eq('mouse')].n.tolist(),[3,2,1,1])
        self.assertEqual(len(detail),3)
        self.assertFalse(detail.set_index('gene_id').loc['b','complete'])

    def test_stale_or_changed_results_cannot_enter_new_supplement(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);config=root/'config.json'
            config.write_text(json.dumps(dict(project_root=str(root),native_run=str(root),frozen_gex=str(root),
                analysis_root=str(root/'analysis'),run_id='new',input_fingerprint='new-counts')))
            context=Context(config,'tables',root/'tables')
            folder=context.analysis/'recurrence';folder.mkdir(parents=True)
            data=folder/'summary.tsv';data.write_text('p\n0.01\n')
            receipt=dict(node='recurrence',status='complete',checks={'new_counts':True},run_id='old',input_fingerprint='new-counts',
                artifacts=[dict(path=data.name,sha256=digest(data))])
            result=folder/'result.json';result.write_text(json.dumps(receipt))
            with self.assertRaises(AssertionError):verified_sources(context,['recurrence'])
            receipt['run_id']='new';result.write_text(json.dumps(receipt))
            context=Context(config,'tables',root/'tables')
            artifacts,index=verified_sources(context,['recurrence'])
            self.assertTrue(index.native_polyapipe_evidence.all())
            data.write_text('p\n0.99\n')
            with self.assertRaises(AssertionError):verified_sources(context,['recurrence'])

    def test_legacy_audit_is_never_labelled_native_and_missing_numeric_is_not_zero(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);config=root/'config.json'
            config.write_text(json.dumps(dict(project_root=str(root),native_run=str(root),frozen_gex=str(root),
                analysis_root=str(root/'analysis'),run_id='new',input_fingerprint='new-counts')))
            context=Context(config,'tables',root/'tables');folder=context.analysis/'consensus';folder.mkdir(parents=True)
            data=folder/'historical.tsv';data.write_text('old_n\n10\n')
            (folder/'result.json').write_text(json.dumps(dict(node='consensus',status='complete',checks={'labelled':True},run_id='new',
                input_fingerprint='new-counts',artifacts=[dict(path=data.name,sha256=digest(data))])))
            _,index=verified_sources(context,['consensus'])
            self.assertFalse(index.native_polyapipe_evidence.any())
            self.assertTrue(index.evidence_scope.eq('historical_original_pipeline').all())
            path=workbook_input(context,pd.DataFrame(dict(gene=['A','B'],effect=[None,0.])),'S8')
            values=json.loads((context.out/path).read_text())
            self.assertIsNone(values['data'][0][1]);self.assertEqual(values['data'][1][1],0.)


if __name__=='__main__':unittest.main()
