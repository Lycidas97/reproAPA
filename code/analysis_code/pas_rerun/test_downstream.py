import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_runtime import Context,bh
from downstream_bundle import run
from downstream_recurrence import numerical,events_from_routes,recurrence_null
from polyapipe_prepare import save
from polyapipe_effects import study

ROOT=Path(__file__).resolve().parents[2]


class DownstreamIntentTests(unittest.TestCase):
    def test_cached_groups_never_cache_gene_depth_or_conflate_profile_cell_denominators(self):
        from downstream_runtime import Study
        study=Study.__new__(Study)
        study.meta=pd.DataFrame(dict(route=['r']*6,library_key=['l']*6,node=['S']*6,
            purity=[.85,.95,.95,.95,.95,.95],n_cells=[4,4,4,4,4,20],metacell_id=list('abcdef')))
        first=np.array([[10,10]]*5+[[1,1]])
        a=study.state(first,'r','l','S')
        self.assertEqual(a['n_cells'],20)
        self.assertFalse(a['endpoint_eligible'])
        self.assertTrue(a['profile_eligible'])
        second=first.copy();second[:5]=0;second[-1]=[30,30]
        b=study.state(second,'r','l','S')
        self.assertEqual(b['rows'].tolist(),[5])
        self.assertFalse(b['profile_eligible'])
        c=study.state(first,'r','l','S','supported50_p09',exclude=['b'])
        self.assertEqual(c['rows'].tolist(),[2,3,4])
        self.assertEqual(c['group_cells'],32)
        d=study.state(first,'r','l','S',depth=40)
        self.assertEqual(d['n_mc'],0)

    def context(self,tmp,node):
        p=Path(tmp);config=p/'config.json'
        config.write_text(json.dumps(dict(project_root=str(ROOT),native_run=str(p/'native'),
            frozen_gex=str(p/'gex'),analysis_root=str(p/'analysis'),run_id='fixture',input_fingerprint='fixture',resources={})))
        return Context(config,node,p/'analysis'/node)

    def test_BH_preserves_order_and_uses_new_family_size(self):
        np.testing.assert_allclose(bh([.04,.001,.03]),[.04,.003,.04])
        np.testing.assert_allclose(bh([.04,.001,.03,1]),[.04*4/3,.004,.04*4/3,1])
        self.assertEqual(len(bh([])),0)

    def test_native_vector_bridge_requalifies_expression_contributors(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'state_vectors');p=ctx.native/'studies/GSE_TEST';p.mkdir(parents=True)
            f=pd.DataFrame(dict(atlas_pas_id=['p0','p1','p2'],feature_index_1based=[1,2,3],
                primary_gene_id=['new_gene']*3,primary_gene_symbol=['New']*3,keep_for_apa_primary=[True]*3,
                chrom=['1']*3,strand=['+']*3,cut_site=[100,200,300],matched_source_count_50nt=[0,1,1],species=['human']*3))
            save(f,p/'features.tsv.gz')
            meta=pd.DataFrame(dict(species=['human']*8,study_id=['GSE_TEST']*8,metacell_id=[f'm{i}' for i in range(8)],
                library_key=['target']*4+['source']*4,route=['R']*8,node=['P']*4+['S']*4,purity=[.85]*8,
                n_cells=[10]*8,matrix_row=np.arange(8)))
            save(meta,p/'metacells.tsv.gz');save(meta,p/'selected_metacells.tsv.gz')
            data=np.array([[0,6,54]]*3+[[0,1,1]]+[[0,54,6]]*3+[[0,1,1]])
            sparse.save_npz(p/'metacell_by_pas_counts.npz',sparse.csr_matrix(data))
            design=pd.DataFrame([dict(species='human',study_id='GSE_TEST',route='R',library_key_a='target',
                library_key_b='source',time_a='D16',time_b='D0',sample_unit_a='T',sample_unit_b='S',contrast='P_S',design='true_source')])
            save(design,p/'endpoint_library_design.tsv');save(design,p/'sampling_design.tsv')
            study(p,[('supported50_p08',True,.8,3)])
            g=ctx.gex/'GSE_TEST';g.mkdir(parents=True)
            gm=meta.copy();gm['total_gex']=1000;save(gm,g/'metacells.tsv.gz')
            save(pd.DataFrame(dict(gene_id=['new_gene'],gene_symbol=['New'])),g/'genes.tsv.gz')
            sparse.save_npz(g/'metacell_by_gene_counts.npz',sparse.csr_matrix(np.array([10,10,10,999,20,20,20,999])[:,None]))
            with patch('downstream_bundle.STUDIES',['GSE_TEST']),patch('downstream_bundle.SETTINGS',{'supported50_p08':(True,.8,3,20)}):run(ctx)
            effect=pd.read_csv(ctx.out/'route_effects.primary.tsv.gz',sep='\t').iloc[0]
            self.assertAlmostEqual(effect.mean_GEX,15000) # PAS-disqualified high-GEX MCs must not leak in.
            self.assertAlmostEqual(effect.delta_DWUI,.4) # Unsupported proximal primary site still defines native axis.
            self.assertEqual(json.loads(effect.PAS_ids),['p1','p2'])
            self.assertTrue((ctx.out/'result.json').exists())

    def test_full_vector_sign_test_cancels_opposing_PAS_shifts(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'pool_inference')
            math=ctx.functions('notes/2026-09-30-multimetric-statistics/pool_inference/analyze.py',['sign_statistics'])
            d=np.array([[.3,-.3],[-.3,.3],[.2,-.2],[-.2,.2]])
            effects,p,mean,null,extreme=math.sign_statistics(d,np.ones(4)/4,np.array([[0,.5,0,1],[1,1,1,0]]))
            self.assertEqual(effects[-1],0.)
            self.assertEqual(p[-1],1.) # Mean of per-unit MPRO would falsely call this large.

    def test_study_vote_rejects_disagreeing_series_and_BH_keeps_planned_positions(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'recurrence');math=numerical(ctx)
            rows=[]
            for study_id,route,d in [('G1','a',.08),('G1','b',-.08),('G2','c',.08)]:
                rows.append(dict(species='human',contrast='P_S',study_id=study_id,route=route,gene_id='new',gene_symbol='New',
                    mean_GEX=100,APA_coverage=30,PAS_number=2,delta_WPP=d,delta_RWUI=d,MPRO=abs(d)*2,
                    reference_WPP=.5,reference_RWUI=.6,rwui_min_weight=.5,rwui_max_weight=1.,
                    reference_concentration=.5,reference_min_usage=.5,reference_max_usage=.5))
            e=events_from_routes(pd.DataFrame(rows),math)
            self.assertTrue(e[e.unit.eq('G1')&e.metric.eq('WPP')].event.eq(0).all())
            self.assertEqual(e[e.unit.eq('G1')].series_pattern.iloc[0],'a;b')
            s,m,*_=recurrence_null(e,math,n=10)
            self.assertEqual(len(s),24)
            self.assertTrue(s.groupby('metric').size().eq(8).all())
            self.assertTrue(s[s.status.ne('tested')].empirical_p.eq(1).all())
            self.assertFalse((m.metric.eq('WPP')).any())


    def test_common_axis_does_not_chain_across_fifty(self):
        from downstream_measurement import common_axes
        frame=pd.DataFrame(dict(study_id=['GSE103221','GSE115943']*3,
            primary_gene_id=['g']*6,atlas_pas_id=list('abcdef'),strand=['+']*6,
            chrom=['1']*6,cut_site=[0,40,45,85,90,130]))
        result=common_axes(frame)
        self.assertTrue(result.cluster_span.le(50).all())
        self.assertEqual(result.cluster_id.nunique(),2)
        self.assertTrue(result.groupby('cluster_id').study_id.nunique().eq(2).all())

    def test_common_axis_uses_cross_study_clusters_and_conserves_selected_counts(self):
        from downstream_measurement import collapse_common_axis
        axis=pd.DataFrame(dict(study_id=['s']*4,atlas_pas_id=list('abcd'),primary_gene_id=['g']*4,
            cluster_id=['native1','native2','native3','excluded'],position=[0,.1,.9,1]))
        mapping=axis.iloc[:3,:3].copy()
        mapping['cluster_id']=['common1','common1','common2']
        mapping['common_position']=[1,1,0]  # Minus-strand order must survive aggregation.
        mapping['common_cut_site']=[20,20,200]
        regions,counts=collapse_common_axis(axis,sparse.csr_matrix([[2,3,7,100],[0,5,1,200]]),mapping)
        self.assertEqual(regions.atlas_pas_id.tolist(),['common1','common2'])
        self.assertEqual(regions['rank'].tolist(),[2.,1.])
        np.testing.assert_array_equal(counts.toarray(),[[5,7],[5,1]])

    def test_intermediate_inference_excludes_unassigned_cells_without_relabeling(self):
        from downstream_intermediate import qualified_state_groups
        meta=pd.DataFrame(dict(route=['WOT']*9,library_key=['u']*3+['r']*3+['p']*3,
            node=['U']*3+['R']*3+['P']*3,purity=[1.]*6+[.85]*3,n_cells=[40]*9))
        self.assertEqual(set(qualified_state_groups(meta,.8)),{('WOT','r','R'),('WOT','p','P')})
        self.assertEqual(set(qualified_state_groups(meta,.9)),{('WOT','r','R')})
        self.assertEqual(meta.node.eq('U').sum(),3)

    def test_expression_symbol_resolution_keeps_gene_axis_and_existing_names(self):
        from downstream_runtime import resolve_gene_symbols
        genes=pd.DataFrame(dict(gene_id=['ENSMUSG1','ENSMUSG2.3','ENSMUSG3','ENSMUSG4'],
            gene_symbol=['ENSMUSG1','ENSMUSG2.3','OldName',None]),index=[8,2,5,9])
        symbols=resolve_gene_symbols(genes,{'ENSMUSG1':'Pcbp1','ENSMUSG2':'Mcm5','ENSMUSG3':'NewName'})
        self.assertEqual(symbols.index.tolist(),[8,2,5,9])
        self.assertEqual(symbols.tolist(),['Pcbp1','Mcm5','OldName',''])
        self.assertEqual(genes.gene_symbol.iloc[0],'ENSMUSG1')

    def test_exact_functional_tail_preserves_stratum_targets(self):
        from collections import defaultdict
        from functools import lru_cache
        from scipy.special import logsumexp
        from scipy.stats import hypergeom
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'functional')
            math=ctx.functions('notes/2026-10-02-null-model-controls/functional.py',
                ['component','distribution','tail','rejection_probability'],
                dict(lru_cache=lru_cache,hypergeom=hypergeom,logsumexp=logsumexp))
            # A fixed all-success stratum contributes an offset; only the second is mobile.
            dist=math.distribution([(3,3,2),(4,2,2)])
            self.assertAlmostEqual(math.tail(dist,4),1/6)
            self.assertEqual(math.tail(dist,2),1.)
            self.assertAlmostEqual(math.rejection_probability(dist,dist,.2),1/6)
            tilted=math.distribution([(3,3,2),(4,2,2)],odds=4)
            self.assertGreater(math.rejection_probability(dist,tilted,.2),1/6)

    def test_species_balance_and_candidate_selection_use_new_genes(self):
        from downstream_profiles import feature,fit
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'profiles')
            math=ctx.functions('notes/2026-10-02-reviewer-followup/trajectory_noise.py',['selection','choose'])
            genes=np.array(['new1','new2','new3']);raw=np.array([[.8,.6,.4,.2],[.2,.4,.6,.8],[.2,.6,.8,.3]])
            data={'WOT':(genes,raw),'CIP':(genes,raw),'GSE118258':(genes,raw),
                '178325_0618':(genes,raw[:,[0,1,3]]),'GSE150311':(genes,raw[:,[0,2,3]])}
            selected=math.choose(data)
            self.assertEqual(set(selected['human']),set(genes))
            x,species,ids=feature(data,selected)
            model,labels,centers,weights=fit(x,species)
            self.assertAlmostEqual(weights[species=='human'].sum(),.5)
            self.assertAlmostEqual(weights[species=='mouse'].sum(),.5)
            self.assertTrue(np.diff(centers[:,-1]-centers[:,0]).min()>=0)

    def test_depth_sensitivity_reports_constant_candidates_without_assigning_a_class(self):
        from downstream_profiles import sensitivity_members,feature
        genes=np.array(['retained','constant','missing'])
        raw=np.array([[.8,.6,.4,.2],[.5,.5,.5,.5]])
        data={r:(genes[:2],raw.copy()) for r in ['CIP','WOT','GSE118258']}
        data['WOT'][1][1]=[.8,.6,.4,.2]
        members={s:genes.copy() for s in ['human','mouse']}
        chosen,audit=sensitivity_members(data,members)
        for sp in members:np.testing.assert_array_equal(chosen[sp],['retained'])
        self.assertEqual(len(audit),9)  # Every original candidate remains in each required route's denominator.
        self.assertEqual(audit[audit.gene_id.eq('missing')].status.unique().tolist(),['missing_complete_profile'])
        self.assertEqual(audit[audit.gene_id.eq('constant')&audit.route.eq('CIP')].status.iloc[0],'constant_or_nonfinite_profile')
        self.assertEqual(audit[audit.gene_id.eq('constant')&audit.route.eq('WOT')].status.iloc[0],'classifiable')
        self.assertTrue(np.isfinite(feature(data,chosen)[0]).all())

    def test_missing_series_cannot_reduce_fixed_majority_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'dwui_pool')
            math=ctx.functions('notes/2026-09-14-integration-reliability/evaluate_reaggregated_apa_effects.py',['study_votes','majority'])
            vote,state=math.study_votes(np.array([[.02,.02],[np.nan,-.02]]),.01)
            np.testing.assert_array_equal(vote,[0,0])
            direction,pos,neg,need=math.majority(np.array([[1],[1],[0],[0]]))
            self.assertEqual(need,3)
            self.assertEqual(direction[0],0)

    def test_regulatory_minus_boundary_and_catalog_many_to_one(self):
        from downstream_regulatory import retained_base,numerical
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'intervals');structure,migration,paired=numerical(ctx)
            # A minus PAS at interbase 100 retains [100,101), not [99,100).
            cut=retained_base(100,'-')
            self.assertEqual(cut,101)
            self.assertEqual(migration.mask_for(100,101,'-',np.array([cut])).tolist(),[1])
            self.assertEqual(migration.mask_for(99,100,'-',np.array([cut])).tolist(),[0])
            raw=np.array([3,7,11]);mapping=np.array([[1,0],[1,0],[0,1]])
            np.testing.assert_array_equal(raw@mapping,[10,11])
            # Mapping counts does not renormalize or lose molecules.
            self.assertEqual((raw@mapping).sum(),raw.sum())

    def test_site_factor_identity_does_not_duplicate_physical_site(self):
        from downstream_sites import site_inventory
        d=pd.DataFrame([dict(species='mouse',route='CIP',gene_id='g',gene_symbol='G',category=cat,chrom='chr1',strand='-',
            start0=100,end0=110,contrast='SP',mask_json='[[1,1],[0,0]]',exact=True,names='RBP1',feature_groups='motif',
            eq20_lo=-.2,eq50_lo=-.3) for cat in ['RBP','RBP_transfer']])
        regions=pd.DataFrame([dict(species='mouse',route='CIP',gene_id='g',contrast='SP',start0=100,end0=120,
            mask='[1,0]',observed_both5=True,effect20=-.2,effect50=-.3,region_id='r')])
        result=site_inventory(d,regions)
        self.assertEqual(len(result),1)
        self.assertEqual(result.provenance.iloc[0],'RBP;RBP_transfer')
        self.assertTrue(result.strict_changed.iloc[0])
        result=site_inventory(d.iloc[:0],regions.iloc[:0])
        self.assertEqual(len(result),0)
        self.assertEqual(result.strict_changed.dtype,bool)

    def test_cross_species_thresholds_keep_missing_mouse_study_as_failure(self):
        from downstream_cross import proximal,qcodes
        f=pd.DataFrame(dict(W=[-.1,-.1],C=[-.1,np.nan]))
        self.assertEqual(proximal(f,['W','C']).tolist(),[True,False])
        h=pd.DataFrame(dict(H1=[-.1,-.1],H2=[-.1,-.1],H3=[np.nan,.03],H4=[np.nan,np.nan]))
        self.assertEqual(proximal(h,['H1','H2','H3','H4']).tolist(),[True,False])
        self.assertEqual(qcodes(pd.Series([1.,1.,np.nan])).tolist(),[0,0,-1])

    def test_background_control_choice_never_uses_motif_density(self):
        from downstream_motifs import helpers
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'motif_background');num,read,choose=helpers(ctx)
            c=np.array([[100,.5,.01,.1],[200,.5,.01,.1],[300,.5,.1,.1]])
            # Identical matching scores use position, not downstream site counts.
            self.assertEqual(choose.choose(c,.05,.2)[0],100)
            self.assertIsNone(choose.choose(c,.001,.2))
            r=num.paired(np.array([1.,2.]),('small',))
            self.assertEqual(r['p'],1.)
            self.assertEqual(r['status'],'fewer_than_10_genes')

    def test_no_changed_intervals_keeps_all_motif_test_positions(self):
        from downstream_motifs import primary,controls,retention
        from downstream_regulatory import ROUTES
        with tempfile.TemporaryDirectory() as tmp:
            ctx=self.context(tmp,'motif_background')
            region=ctx.analysis/'intervals/all_coordinate_regions.tsv.gz';region.parent.mkdir()
            save(pd.DataFrame(columns=['species','route','gene_id','contrast','start0','end0','effect20','effect50','observed_both5','paired_observed_both5','strict_change']),region)
            for sp,routes in ROUTES.items():
                for sid,route in routes:
                    p=ctx.analysis/f'intervals/annotation/{sp}/{route}/models.json';p.parent.mkdir(parents=True);p.write_text('{}')
                p=ctx.analysis/f'site_scan/{sp}/query.fa';p.parent.mkdir(parents=True);p.write_text('')
                p=ctx.analysis/f'site_scan/full-cohort-v2/data/{sp}/site_effects.tsv.gz';p.parent.mkdir(parents=True)
                save(pd.DataFrame(columns=['gene_id','category','chrom','strand','start0','end0']),p)
            primary(ctx)
            result=pd.read_csv(ctx.out/'density_summary.tsv',sep='\t')
            self.assertEqual(len(result),8);self.assertTrue(result.p.eq(1).all())
            ctx=Context(ctx.config_path,'motif_controls',ctx.analysis/'motif_controls');controls(ctx)
            grid=pd.read_csv(ctx.out/'grid_summary.tsv',sep='\t')
            self.assertEqual(len(grid),240)
            ladder=pd.read_csv(ctx.out/'background_ladder.tsv',sep='\t')
            self.assertEqual(len(ladder[ladder.cohort.eq('six_background_common')]),48)
            self.assertTrue(ladder.p.eq(1).all())
            ctx=Context(ctx.config_path,'weighted_retention',ctx.analysis/'weighted_retention');retention(ctx)
            self.assertTrue(pd.read_csv(ctx.out/'retention_summary.tsv',sep='\t').genes.eq(0).all())


if __name__=='__main__':unittest.main()
