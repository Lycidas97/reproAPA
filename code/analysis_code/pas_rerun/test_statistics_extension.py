"""Scientific boundary tests for the remaining statistical adapters."""
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_runtime import Context
from downstream_continuous import effects_frame,metrics,prepare,add_bins,null_record
from downstream_global import backgrounds,fits_full_primary
from downstream_shared_states import clusters,conditional,compare,state_groups
from downstream_apa_gex import permutations
from downstream_consensus import consensus,joint_draw
from downstream_examples import shared_regions,query_features
from downstream_external_structure import retention,mapped_interval,numerics
from downstream_pas_qc import sequence,region_class
from downstream_coordinate_states import rebuild
from downstream_regulatory import F5

ROOT=Path(__file__).resolve().parents[2]


class ExtensionIntentTests(unittest.TestCase):
    def test_cached_metadata_groups_preserve_exclusions_and_missing_support(self):
        meta=pd.DataFrame(dict(purity=[.95,.85,.95,.95,.95,.95],node=['R','R','R','S','E','P'],
            route=['WOT']*5+['reference'],library_key=['a','a','b','b','b','c'],
            metacell_id=list('abcdef'),n_cells=[10]*6))
        support=pd.DataFrame(dict(support_q95=[True,False,True],support_q90=[True,True,True]),index=list('abc'))
        groups=state_groups(meta,'GSE115943',np.array(['L00','L01','L01','L01','L01','L01']),support)
        rows=lambda setting:sorted(int(i) for _,ix in groups[setting] for i in ix)
        self.assertEqual(rows('shared_p08'),[0,1,2,3,4])
        self.assertEqual(rows('shared_p09'),[0,2,3,4])
        self.assertEqual(rows('exclude_L00'),[1,2,3,4])
        self.assertEqual(rows('support_q95'),[0,2])
        self.assertEqual(rows('support_q90'),[0,1,2])
        # Grouping must not pre-apply a gene-specific depth gate or its library size:
        # different PAS axes can leave different eligible contributors in each group.
        self.assertEqual(groups['shared_p08'][0][1].tolist(),[0,1])

    def context(self,tmp,node):
        root=Path(tmp);p=root/'config.json'
        p.write_text(json.dumps(dict(project_root=str(ROOT),native_run=str(root/'native'),frozen_gex=str(root/'gex'),
            analysis_root=str(root/'analysis'),run_id='fixture',input_fingerprint='fixture',resources={})))
        return Context(p,node,root/'analysis'/node)

    def test_consensus_preserves_mode_dependence_and_missing_study_mask(self):
        matrix=np.array([[[1,1,1,1,1,1],[1,1,1,1,1,1]],
                         [[1,1,1,-1,1,1],[1,1,1,-1,1,1]],
                         [[9,9,9,9,9,9],[1,1,1,1,1,1]]],np.int8)
        self.assertEqual(consensus(matrix)[0].tolist(),[True,False,False])
        rng=np.random.default_rng(42)
        for _ in range(30):
            shuffled=joint_draw(matrix,{0:[np.array([0,1])],1:[np.array([0,1,2])]},rng)
            np.testing.assert_array_equal(shuffled[2,0],matrix[2,0])
            for unit in [0,1]:
                self.assertEqual(sorted(map(tuple,shuffled[:,unit,:])),sorted(map(tuple,matrix[:,unit,:])))

    def test_native_boundary_uncertainty_cannot_be_reported_as_exact_site_loss(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'site_scan')
            core=ctx.functions(F5+'full-cohort-v2/sources/execution-scripts/core.py',['masks','change'])
            pas=[dict(member_cut_site_min=90,member_cut_site_max=110),
                 dict(member_cut_site_min=200,member_cut_site_max=200)]
            mask=core.masks(95,100,'+',pas)
            # The shifted proximal cut may retain this site, so zero remains possible.
            np.testing.assert_allclose(core.change([.8,.2],[.2,.8],mask),[-.6,0])
            exact=core.masks(95,100,'+',[dict(member_cut_site_min=90,member_cut_site_max=90),pas[1]])
            np.testing.assert_allclose(core.change([.8,.2],[.2,.8],exact),[-.6,-.6])

    def test_state_exclusion_requalifies_libraries_and_keeps_outside_family_denominator(self):
        f=pd.DataFrame(dict(atlas_pas_id=['a','b','c','outside'],feature_index_1based=[1,2,3,4],
            primary_gene_id=['g']*4,matched_source_count_50nt=[1]*4,
            keep_for_apa_primary=[True]*4,transcript_model_conflict=[False]*4,
            multi_gene_ambiguous=[False]*4,internal_priming_flag=[False]*4))
        study=SimpleNamespace(f=f,counts=sparse.csc_matrix(np.tile([10,20,30,40],(6,1))))
        raw=dict(gene_id='g',pas=[dict(atlas_pas_id=p) for p in 'abc'])
        annotation=dict(gene_id='g',pas=[dict(atlas_pas_id=p) for p in ['A','B']])
        mapping=pd.DataFrame(dict(native_pas_id=list('abc'),annotation_id=['A','A','B']))
        meta=pd.DataFrame(dict(node=['S']*3+['E']*3,library_key=['source']*3+['middle']*3,n_cells=[10]*6))
        strict=lambda frame:pd.Series(True,index=frame.index)
        native,annotated=rebuild(study,raw,annotation,mapping,meta,strict)
        np.testing.assert_array_equal(annotated['profiles']['S']['counts'],[90,90])
        self.assertEqual(annotated['profiles']['S']['coverage'],.6)
        self.assertEqual(native['profiles']['S']['total'],annotated['profiles']['S']['total'])
        _,excluded=rebuild(study,raw,annotation,mapping,meta.drop(index=5),strict)
        self.assertNotIn('E',excluded['profiles'])
        np.testing.assert_array_equal(excluded['profiles']['S']['counts'],annotated['profiles']['S']['counts'])
        from downstream_coordinate_states import prepare_family
        prepared=prepare_family(study,raw,annotation,mapping,strict)
        _,reused=rebuild(study,raw,annotation,mapping,meta.drop(index=5),strict,prepared)
        self.assertEqual(set(reused['profiles']),set(excluded['profiles']))
        np.testing.assert_array_equal(reused['profiles']['S']['counts'],excluded['profiles']['S']['counts'])
        self.assertEqual(reused['profiles']['S']['coverage'],excluded['profiles']['S']['coverage'])

    def test_missing_intermediate_cannot_be_imputed_as_late_onset(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'coordinate_states')
            onset=ctx.functions(F5+'site-recurrence-v1/supplement-compressed-v1/prepare.py',['first_observed'])
            self.assertEqual(onset.first_observed({20:np.array([.8,np.nan,.5,.2]),50:np.array([.8,.7,.5,.2])},.1),'incomplete')
            self.assertEqual(onset.first_observed({20:np.array([.8,.75,.5,.2]),50:np.array([.8,.75,.5,.2])},.1),'E')

    def test_new_shared_interval_is_intersection_not_historical_feature(self):
        rows=[]
        for r in [(100,250,180,300,190,197),(100,250,180,300,210,217)]:
            rows.append(dict(species='mouse',contrast='SP',gene_id='new',gene_symbol_a='Rad23b',chrom='chr4',strand='+',
                region_id_a='a',region_id_b='b',region_start0_a=r[0],region_end0_a=r[1],region_start0_b=r[2],region_end0_b=r[3],
                start0=r[4],end0=r[5],category='miRNA',factor='seed',recurrent_strict=True))
        sites=pd.DataFrame(rows);regions=shared_regions(sites);self.assertEqual(len(regions),1)
        self.assertEqual(regions[['shared_start0','shared_end0']].iloc[0].tolist(),[180,250])
        features=query_features(regions,sites);self.assertEqual(len(features),3)
        self.assertFalse(features.end0.sub(features.start0).eq(415).any())

    def test_exonic_overlap_cannot_substitute_for_terminal_inclusion(self):
        self.assertEqual(retention([(0,50),(100,200)],'+',20,30),('partial_or_spliced',10))
        self.assertEqual(retention([(0,50),(100,200)],'+',150,170),('included',20))
        self.assertEqual(retention([(0,50),(100,200)],'-',20,30),('included',10))
        self.assertEqual(retention([(0,50),(100,200)],'+',60,80),('excluded',0))

    def test_liftover_ambiguity_is_explicitly_unmeasurable(self):
        with tempfile.TemporaryDirectory() as tmp:
            math=numerics(self.context(tmp,'longread'))
            forward={'chr1':[(0,100,'chrA',200,'+',1000,'chain')]};reverse={'chrA':[(200,300,'chr1',0,'+',1000,'back')]}
            self.assertEqual(mapped_interval(math,forward,reverse,'chr1',10,20)['start0'],210)
            forward['chr1'].append((0,100,'chrA',200,'+',1000,'duplicate'))
            self.assertEqual(mapped_interval(math,forward,reverse,'chr1',10,20)['status'],'unmappable_or_nonreciprocal')

    def test_sequence_centers_on_retained_base_on_both_strands(self):
        class Fasta:
            references=['chr1']
            def get_reference_length(self,chrom):return 400
            def fetch(self,chrom,start,end):return ('AACGT'*80)[start:end]
        plus=sequence(Fasta(),'1',151,'+');minus=sequence(Fasta(),'chr1',150,'-')
        expected=('AACGT'*80)[50:251]
        self.assertEqual(plus,expected)
        self.assertEqual(minus,expected.translate(str.maketrans('ACGT','TGCA'))[::-1])
        model={'transcript':[(0,300)],'exon':[(0,300)],'CDS':[(50,200)]}
        self.assertEqual(region_class(model,210,'+'),'three_prime_UTR')
        self.assertEqual(region_class(model,40,'-'),'three_prime_UTR')

    def test_unmeasured_setting_cannot_shrink_fixed_background_requirement(self):
        q=pd.DataFrame([dict(species='mouse',study_id='s',route='r',contrast='P_S',design='true_source',setting=s,gene_id='g',delta=.2)
                        for s in ['supported50_p08','supported50_p09','expanded_primary_p08']])
        orth=pd.DataFrame(dict(human_gene_id=['h'],mouse_gene_id=['g']))
        summary,_=backgrounds(q,orth,[('mouse','s','r','P_S','true_source')])
        self.assertTrue(summary[summary.background.eq('fixed_settings_within_route')].n.eq(0).all())
        self.assertEqual(len(summary),16)

    def test_tandem_requires_full_primary_and_retained_minus_base(self):
        axis=pd.DataFrame(dict(chrom=['1','1'],strand=['-','-'],cut_site=[99,199]))
        utr=pd.DataFrame(dict(chrom=['1'],strand=['-'],utr_start=[100],utr_end=[200],transcript_id=['t']))
        self.assertEqual(len(fits_full_primary(axis,utr)),1)
        axis.loc[2]=['1','-',250]
        self.assertEqual(len(fits_full_primary(axis,utr)),0)

    def test_matched_states_use_intersection_and_measurement_gate(self):
        rows=[]
        for code in ['W','H1']:
            for gene in ['all','no_R','failed_P']:
                for node,value in [('S',.6),('R',.4),('P',.2)]:
                    if gene=='no_R' and node=='R':continue
                    rows.append(dict(layer='supported50',code=code,human_gene_id=gene,node=node,position_index=value,
                        measurement_pass=not(gene=='failed_P' and node=='P')))
        nodes=pd.DataFrame(rows);paired=effects_frame(nodes,'supported50','W','H1',['R','P'])
        self.assertEqual(paired.index.tolist(),['all'])
        np.testing.assert_allclose(paired.iloc[0].to_numpy(),[-.2,-.4,-.2,-.4])
        self.assertEqual(null_record(np.nan,np.zeros(9))['p_greater'],1.)
        self.assertTrue(np.isnan(metrics(prepare([],[]))).all())

    def test_constant_discovery_magnitude_keeps_one_stratum(self):
        x=pd.DataFrame(dict(disc_meta=[.2]*3,disc_direction=[1]*3,disc_eligible=[2]*3,val_eligible=[2]*3))
        self.assertEqual(add_bins(x).effect_bin.nunique(),1)

    def test_common_pas_does_not_bridge_chain_and_sums_many_to_one(self):
        def axis(s,values):return pd.DataFrame(dict(study_id=[s]*len(values),cut_site=values,atlas_pas_id=[s+str(i) for i in range(len(values))],chrom=['1']*len(values),strand=['+']*len(values)))
        groups=clusters([axis('a',[0,20,90]),axis('b',[40,120])],['a','b'])
        self.assertEqual([len(g) for g in groups],[3,2])
        self.assertTrue(all(max(r['cut_site'] for r in g)-min(r['cut_site'] for r in g)<=50 for g in groups))

    def test_conditional_collinearity_is_not_reported_as_a_fitted_effect(self):
        x=np.column_stack([np.ones(6),[0,0,0,1,1,1],[0,0,0,1,1,1]])
        self.assertTrue(np.isnan(conditional(x,np.arange(6))['adjusted']))

    def test_common_library_triangle_preserves_exact_contributors(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'intermediate');rows=[]
            for lib,node,apa,gex,counts in [('src','S',.7,3,[30,70]),('target','P',.2,4,[80,20]),('target','M',.4,5,[60,40])]:
                rows.append(dict(axis_id='a__b',species='mouse',study_id='a',route='route',setting='shared_p08',gene_id='new',gene_symbol='New',
                    library_key=lib,state=node,time='D0' if lib=='src' else 'D10',sample_unit=lib,replicate_status='unknown',APA=apa,GEX=gex,
                    GEX_all=gex,GEX_count=100,GEX_CPM=10,n_mc=3,n_cells=30,pas_counts=100,site_counts=json.dumps(counts),site_mc='[3, 3]'))
            compare(ctx,pd.DataFrame(rows))
            tri=pd.read_csv(ctx.out/'shared/common_window_triangle.tsv.gz',sep='\t')
            self.assertEqual(len(tri),1)
            self.assertLess(abs(tri.APA_delta_triangle_error.iloc[0]),1e-12)
            self.assertEqual(tri.GEX_delta_triangle_error.iloc[0],0)

    def test_empty_APA_GEX_null_retains_both_shuffle_directions(self):
        columns=['species','guan_series_pattern']
        columns += [field+'_'+suffix for suffix in ['1','2'] for field in ['alo','ahi','glo','ghi','meas','GEX','depth','PAS','reference_WPP','reference_WPP_min','reference_WPP_max']]
        q=pd.DataFrame(columns=columns)
        for field in columns:
            if field not in ['species','guan_series_pattern']:q[field]=q[field].astype(float)
        rows,_=permutations(q,'human',.02,'direction_specific',5,draws=9)
        self.assertEqual(len(rows),8)
        self.assertTrue(all(r['p_enrich']==1 and r['observed']==0 for r in rows))


if __name__=='__main__':unittest.main()
