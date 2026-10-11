"""Intent checks for new native external and cis analyses."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize,brentq
from scipy.stats import chi2

import test_statistics_extension as fixtures
from downstream_external_pas import usable_matches,endpoint_samples,ENDPOINT
from downstream_cis import select_pairs
from downstream_cis_models import SOURCE,folds,auc,design_status,models,REC,CIS
from downstream_liftover import lift_points


class RemainingIntentTests(unittest.TestCase):
    context=fixtures.ExtensionIntentTests.context

    def test_external_sequence_cannot_silently_use_wrong_assembly_or_interval(self):
        import json
        from downstream_external_structure import checked_mm9_sequence
        data=dict(genome='mm9',chrom='chr4',start=10,end=17,dna='aagtgct')
        self.assertEqual(checked_mm9_sequence(json.dumps(data),('chr4',10,17)),'AAGTGCT')
        for changes in [dict(genome='mm10'),dict(chrom='chr5'),dict(start=11),dict(end=18),dict(dna='AAGTGC')]:
            with self.assertRaises(AssertionError):checked_mm9_sequence(json.dumps(data|changes),('chr4',10,17))

    def test_candidate_state_tracking_keeps_unmeasured_settings_and_entrance_scope(self):
        from downstream_apa_gex import candidate_relations
        effects=pd.DataFrame([dict(species='mouse',axis_id='mouse_pair',setting=setting,gene_id=gene,
            study_id=sid,contrast='P_S',design='true_source',APA_delta=.05,GEX_delta=gex,GEX_all_pairs_measurable=True)
            for setting,gene in [('primary','g'),('strict','other')]
            for sid,gex in [('a',.5),('b',-.5)]])
        members=pd.DataFrame([dict(species='mouse',metric='WPP',threshold=.02,series_matching=True,
            contrast='P_M',gene_id='g',gene_symbol='G',direction='positive')])
        out={};ctx=SimpleNamespace(analysis=Path('/candidate-test'),
            read=lambda p:effects if 'matched_route' in str(p) else members,
            table=lambda f,name,columns:out.update(result=pd.DataFrame(f,columns=columns)))
        candidate_relations(ctx);f=out['result'];self.assertEqual(len(f),20)
        self.assertTrue(f[f.setting.eq('strict')].status.eq('fewer_than_two_studies').all())
        row=f[f.setting.eq('primary')&f.contrast.eq('P_S')].iloc[0]
        self.assertTrue(row.recurrent and row.GEX_opposite)
        self.assertFalse(row.source_interpretation_allowed)

    def test_external_endpoint_identity_excludes_unmanifested_perturbations(self):
        manifest=pd.DataFrame(dict(count_column=ENDPOINT,group=['MEF','MEF','iPSC','iPSC'],
            gsm=['gsm1','gsm2','gsm3','gsm4'],sra_run=['run1','run2','run3','run4']))
        selected=endpoint_samples(manifest,list(reversed(ENDPOINT))+['Nudt_D3_R1','Day6_SSEAI+'])
        self.assertEqual(selected.index.tolist(),ENDPOINT)
        with self.assertRaisesRegex(AssertionError,'absent from counts'):
            endpoint_samples(manifest,ENDPOINT[:-1])
        with self.assertRaisesRegex(AssertionError,'absent from manifest'):
            endpoint_samples(manifest.iloc[:-1],ENDPOINT)
        with self.assertRaisesRegex(AssertionError,'Duplicate'):
            endpoint_samples(pd.concat([manifest,manifest.iloc[:1]]),ENDPOINT)

    def test_external_nearest_ties_and_collisions_do_not_enter_common_axis(self):
        ext=pd.DataFrame([dict(COORDS='a',gene_id='g',mapped_chrom='chr1',orientation='+',source_strand='+',coordinate_base_sensitive=False,mapped_cut=100,alternate_cut=101),
            dict(COORDS='b',gene_id='g',mapped_chrom='chr1',orientation='+',source_strand='+',coordinate_base_sensitive=False,mapped_cut=101,alternate_cut=100),
            dict(COORDS='c',gene_id='g',mapped_chrom='chr1',orientation='+',source_strand='+',coordinate_base_sensitive=False,mapped_cut=150,alternate_cut=150)])
        ref=pd.DataFrame(dict(primary_gene_id=['g']*3,chrom=['chr1']*3,strand=['+']*3,cut_site=[100,140,160],atlas_pas_id=['x','y','z']))
        out=usable_matches(ext,ref)
        self.assertEqual(out.compatible50.tolist(),[True,True,False])
        self.assertFalse(out.usable.any())

    def test_liftover_follows_retained_base_and_rejects_multiple_destinations(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            ctx=SimpleNamespace(out=directory,resource=lambda key:directory/key,
                table=lambda frame,name:frame.to_csv(directory/name,sep='\t',index=False))
            f=pd.DataFrame(dict(id=['plus','minus','ambiguous'],chrom=['chr1']*3,strand=['+','-','+'],cut_site=[100,100,200]))
            def fake(argv,**kwargs):
                source=pd.read_csv(argv[-4],sep='\t',header=None)
                self.assertEqual(source[1].tolist(),[99,100,199])
                Path(argv[-2]).write_text('chr2\t1000\t1001\tplus\t0\t-\nchr2\t1000\t1001\tminus\t0\t+\nchr2\t2000\t2001\tambiguous\t0\t+\nchr3\t2000\t2001\tambiguous\t0\t+\n')
                Path(argv[-1]).write_text('')
            with patch('downstream_liftover.subprocess.run',fake):out=lift_points(ctx,f,'chain','fixture')
            self.assertEqual(out.mapped_cut.iloc[:2].tolist(),[1000,1001])
            self.assertEqual(out.mapping_status.tolist(),['unique','unique','ambiguous'])
            self.assertTrue(pd.isna(out.mapped_cut.iloc[2]))

    def test_weak_one_to_one_match_cannot_steal_a_reciprocal_pair(self):
        f=pd.DataFrame([dict(human_id='h1',mouse_id='m1',reciprocal25=True,reciprocal50=True,forward_distance=24,backward_distance=20),
            dict(human_id='h2',mouse_id='m1',reciprocal25=False,reciprocal50=False,forward_distance=1,backward_distance=90),
            dict(human_id='h1',mouse_id='m2',reciprocal25=False,reciprocal50=True,forward_distance=30,backward_distance=30)])
        out=select_pairs(f)
        self.assertEqual(out[['human_id','mouse_id']].to_numpy().tolist(),[['h1','m1']])

    def test_paired_cv_folds_keep_classes_and_auc_ties(self):
        y=np.r_[np.zeros(23),np.ones(17)];a=folds(y,42);b=folds(y,42)
        np.testing.assert_array_equal(a,b)
        for value in [0,1]:self.assertLessEqual(np.ptp(np.bincount(a[y==value],minlength=5)),1)
        self.assertEqual(auc(y,np.ones(len(y))),.5)
        self.assertEqual(auc(y,y),1.)

    def test_Firth_gradient_and_degenerate_design_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'cis');math=ctx.functions(SOURCE,['value_grad','fit','design','profile_ci','verify'],dict(expit=expit,minimize=minimize,brentq=brentq,chi2=chi2))
            self.assertTrue(math.verify()['complete_separation_finite'])
            f=pd.DataFrame(dict(a=[1.,2.,3.],b=[2.,4.,6.]))
            self.assertEqual(design_status(f,['a','b'],math)[0],'rank_deficient')
            f['b']=np.nan
            self.assertEqual(design_status(f,['a','b'],math)[0],'all_missing_predictor')

    def test_empty_cis_universe_keeps_all_planned_test_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'cis')
            frame=pd.DataFrame(columns=['layer','contrast','discovery','tau','human_gene_id','Y',*REC,*CIS])
            models(ctx,frame)
            out=pd.read_csv(ctx.out/'recomputed_models.tsv',sep='\t')
            self.assertEqual(len(out),32)
            self.assertTrue(out.p_boot.eq(1).all())
            self.assertTrue(out.q_boot_planned_32.eq(1).all())
            self.assertTrue(out.status.eq('fewer_than_three_genes').all())

    def test_source_endpoint_is_qualified_independently_of_P_M_pairs(self):
        from downstream_inference import source_states,full_study_axes
        from downstream_runtime import Study
        meta=pd.DataFrame(dict(route=['WOT']*6,library_key=['source']*3+['target']*3,node=['S']*3+['P']*3,purity=[1.]*6,n_cells=[10]*6,metacell_id=list('abcdef')))
        axis=pd.DataFrame(dict(position=[0.,1.],feature_index_1based=[1,2]))
        study=Study.__new__(Study);study.id='GSE115943';study.meta=meta
        study.genes={'measured':axis,'unqualified':axis};study.axis=lambda g:axis
        study.gene_counts=lambda a:np.tile([10,10],(6,1))
        axes=full_study_axes(study);out=source_states(study,axes)
        self.assertEqual(set(axes),{('GSE115943','measured'),('GSE115943','unqualified')})
        self.assertEqual(set(out.state),{'S','P'})
        self.assertEqual(len(out),4)

    def test_empty_cross_species_background_retains_56_planned_tests(self):
        from downstream_recurrence import cross_species_current,numerical
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'cross_sharing')
            events=pd.DataFrame(columns=['metric','threshold','contrast','unit','gene_id','event'])
            orth=pd.DataFrame(columns=['human_gene_id','mouse_gene_id'])
            cross_species_current(ctx,numerical(ctx),events,orth,100,42)
            out=pd.read_csv(ctx.out/'sharing_summary.tsv',sep='\t')
            self.assertEqual(len(out),56)
            self.assertTrue(out.q_joint_all_56.eq(1).all())
            self.assertEqual(int(out.q_joint_primary_24.notna().sum()),24)
            self.assertTrue(out.status.eq('no_common_measurable_orthologs').all())

    def test_empty_Range10_background_keeps_its_own_24_DWUI_tests(self):
        from downstream_recurrence import cross_species_current,numerical
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'cross_sharing_range10')
            events=pd.DataFrame(columns=['metric','threshold','contrast','unit','gene_id','event'])
            orth=pd.DataFrame(columns=['human_gene_id','mouse_gene_id'])
            cross_species_current(ctx,numerical(ctx),events,orth,100,42,model='Range10')
            out=pd.read_csv(ctx.out/'sharing_summary.tsv',sep='\t')
            self.assertEqual(len(out),24)
            self.assertTrue(out.q_DWUI_all_24.eq(1).all())
            self.assertEqual(int(out.q_primary_8.notna().sum()),8)
            self.assertNotIn('q_joint_primary_24',out)

    def test_new_trajectory_themes_use_new_names_and_no_forced_old_six(self):
        from downstream_profile_functional import choose_themes
        rows=[]
        for sp in ['human','mouse']:
            rows.append(dict(species=sp,shared_term='GO:NEW',database='GO',cluster='C1',term_id='GO:NEW',term_name='New function',
                background_genes=';'.join(map(str,range(20))),hits=6,fold_enrichment=2,p_value=.01,q_value=.1))
        chosen,audit=choose_themes(pd.DataFrame(rows))
        self.assertEqual(chosen.label.tolist(),['New function'])
        self.assertEqual(len(chosen),1)
        empty=pd.DataFrame(rows);empty['hits']=0
        self.assertEqual(len(choose_themes(empty)[0]),0)

    def test_full_trajectory_functional_background_is_not_the_selected_cohort(self):
        from downstream_profile_functional import functional_ranges
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'profile_gex');rows=[];members=[]
            for sp,routes in [('mouse',['WOT','CIP']),('human',['GSE118258'])]:
                for gene in range(40):
                    if gene<12:members.append(dict(species=sp,gene_id=str(gene),cluster='C'+str(gene//4+1)))
                    for route in routes:
                        for setting in ['supported50_p08','expanded_primary_p08']:
                            for j,state in enumerate('SREP'):
                                rows.append(dict(species=sp,route=route,gene_id=str(gene),setting=setting,depth=20,node=state,
                                    position_index=.1+.001*gene*j,gex_log2_cpm1=1+gene*.05+j*.01,mean_pas_umi=30+gene,
                                    qualified_metacells=10+gene%7,within_library_apa_range=.2))
            ann={db:{sp:{str(g):{'T'} if g<20 else {'U'} for g in range(40)} for sp in ['mouse','human']} for db in ['GO','KEGG']}
            names={db:{t:{'name':t+' new'} for t in ['T','U']} for db in ['GO','KEGG']}
            with patch('downstream_profile_functional.annotations',lambda ctx:(ann,names,None)):
                functional_ranges(ctx,pd.DataFrame(rows),pd.DataFrame(members))
            audit=pd.read_csv(ctx.out/'enrichment_background_audit.tsv',sep='\t')
            self.assertTrue(audit.evaluable_background.eq(40).all())
            self.assertTrue(audit.foreground.eq(4).all())
            rank=pd.read_csv(ctx.out/'rank_enrichment_all_tests.tsv.gz',sep='\t')
            self.assertEqual(len(rank),48)
            self.assertTrue(rank.n_background.eq(40).all())


    def test_common_caller_projection_conserves_molecules_and_excludes_ties(self):
        from downstream_caller_comparison import nearest,project_old
        from scipy import sparse
        old=pd.DataFrame(dict(primary_gene_id=['g']*5,chrom=['chr1']*5,strand=['-']*5,
            atlas_pas_id=list('abcde'),cut_site=[101,102,150,199,300]))
        new=pd.DataFrame(dict(primary_gene_id=['g']*2,chrom=['chr1']*2,strand=['-']*2,
            atlas_pas_id=['x','y'],cut_site=[100,200],position=[1.,0.]))
        counts=sparse.csr_matrix([[2,3,7,11,13]])
        shared,projected,mapping=project_old(old,counts,new)
        self.assertTrue(mapping.set_index('site_id').loc['c','tie'])
        np.testing.assert_array_equal(projected.toarray(),[[5,11]])
        np.testing.assert_array_equal(shared.position,[1.,0.])
        new['primary_gene_id']='other'
        self.assertTrue(nearest(old,new).distance.isna().all())

    def test_matched_caller_effect_requires_same_qualified_metacells(self):
        from downstream_caller_comparison import numerics
        from scipy import sparse
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'caller_comparison');math=numerics(ctx)
            meta=pd.DataFrame(dict(metacell_id=list('abcdef'),route=['r']*6,library_key=['L']*6,node=['P']*3+['M']*3,purity=[1.]*6,n_cells=[10]*6))
            sites=pd.DataFrame(dict(site_id=['x','y'],chrom=['chr1']*2,strand=['+']*2,cut_site=[100,200],primary_axis_position=[.2,.8]))
            tool=sparse.csr_matrix([[20,0],[20,0],[0,20],[0,20],[0,20],[20,0]])
            reference=sparse.csr_matrix([[20,0],[20,0],[0,0],[0,20],[0,20],[20,0]])
            design=pd.DataFrame([dict(species='mouse',study_id='x',route='r',contrast='P_M',design='same_library',library_key_a='L',library_key_b='L',time_a='D1',time_b='D1',sample_unit_a='L',sample_unit_b='L')])
            frozen=design.assign(gene_id='g',gene_symbol='g',weight=1.,delta_WPP=.2)
            pairs,routes=math.gene_effects('g','g',sites,tool,sites,tool,reference,meta,design,frozen,True)
            out=pd.DataFrame(routes)
            common=out[out.view.eq('matched_sites_common_metacells')]
            self.assertTrue(common.delta_primary_axis.isna().all())
            self.assertTrue(common.tool_qualified_pairs.eq(0).all())
            self.assertEqual(len(pairs),3)


    def test_intermediate_missingness_is_not_converted_to_zero_or_reselected_sensitivity(self):
        from downstream_process_classes import process_classes
        from polyapipe_prepare import digest
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'intermediate');rec=ctx.analysis/'recurrence';rec.mkdir()
            pd.DataFrame([dict(species='mouse',metric='WPP',threshold=.02,contrast='P_S',series_matching=True,gene_id='g',gene_symbol='G',direction='positive')]).to_csv(rec/'recurrence_members.tsv.gz',sep='\t',index=False)
            rows=[]
            for setting in ['shared_p08','shared_p09']:
                for sid,route in [('GSE103221','CIP'),('GSE115943','WOT')]:
                    for state,v in [('S',.1),('R',.2 if route=='CIP' else .05),('E',.3),('P',.5)]:
                        if setting=='shared_p09' and route=='WOT' and state=='R':continue
                        rows.append(dict(axis_id='GSE103221__GSE115943',species='mouse',study_id=sid,route=route,setting=setting,
                            gene_id='g',gene_symbol='G',library_key=sid+state,time='D0' if state=='S' else 'D10',sample_unit=sid+state,state=state,APA=v))
            f=pd.DataFrame(rows);native=f[f.setting.eq('shared_p08')].copy();native['setting']='supported50_p08';native.rename(columns={'APA':'value'}).to_csv(ctx.out/'library_state_values.tsv.gz',sep='\t',index=False)
            resource=Path(tmp)/'expression.tsv';pd.DataFrame(dict(species=['mouse'],study_id=['GSE103221'],node=['R'],metacell_id=['mc'],support_q95=[True],support_q90=[False])).to_csv(resource,sep='\t',index=False)
            ctx.config['resources']['state_expression_support']=dict(path=str(resource),sha256=digest(resource))
            process_classes(ctx,f)
            outcomes=pd.read_csv(ctx.out/'process/fixed_cohort_outcomes.tsv',sep='\t')
            row=outcomes[outcomes.variant.eq('shared_p09')&outcomes.scope.eq('requalified')&outcomes.state.eq('R')].iloc[0]
            self.assertEqual(row.fixed_cohort,1);self.assertEqual(row.unavailable,1)
            classes=pd.read_csv(ctx.out/'process/same_series_process_classes.tsv.gz',sep='\t')
            base=classes[classes.variant.eq('shared_p08')&classes.scope.eq('requalified')&classes.state.eq('R')]
            self.assertEqual(base.pair_class.tolist(),['direction_difference'])


if __name__=='__main__':unittest.main()
