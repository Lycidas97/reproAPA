"""Protect v2-only test families, estimands, empty cohorts and exact delivery gates."""
import json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from downstream_direction import correction_families,family_registry
from downstream_display import rank_bins,distributions
from downstream_mirna import grouped_difference,influence,position_tests
from export_statistics_bundle import require_numerical_products


class CompletionTests(unittest.TestCase):
    def test_reference_gate_accepts_both_caller_receipts_but_rejects_v3_conflicts(self):
        from prepare_pipeline import primary_reference
        primary_reference({'reference_variant':'v2_primary'}, {})
        primary_reference({}, {'support_reference_variant':'v2_primary'})
        for manifest,annotation in [({},{}),({'reference_variant':'v3_ge62_comparison'},{}),
            ({'reference_variant':'v2_primary'},{'support_reference_variant':'v3_ge62_comparison'})]:
            with self.assertRaises(ValueError):primary_reference(manifest,annotation)

    def test_v2_family_recomputed_without_v3_or_q160(self):
        a=pd.DataFrame([dict(version='v2',metric='DWUI',threshold=.02,pair=str(i),scope='pair',
            test='distal_minus_proximal_excess_rate' if i<24 else 'log_distal_vs_proximal_relative_enrichment',p_primary=(i+1)/100) for i in range(48)])
        b=a.iloc[:32].copy();a,b=correction_families(a,b)
        np.testing.assert_allclose(np.r_[a.q_complete_80,b.q_complete_80],false_discovery_control(np.r_[a.p_primary,b.p_primary]))
        self.assertEqual(family_registry(a,b).test_id.nunique(),80)
        self.assertFalse(any('160' in c for c in a))
        a.loc[0,'version']='v3'
        with self.assertRaises(AssertionError):correction_families(a,b)

    def test_rank_bin_SE_uses_current_membership_and_source_order_ties(self):
        f=pd.DataFrame([dict(design='PMS',tolerance=25,mass_min=.9,pair=pair,human_gene_id='g'+str(i),
            **{f'rank_{co}_{sp}':0. if sp=='h' else i for co in ['PM','PS','MS'] for sp in ['h','m']}) for pair in ['H1_C','H1_W'] for i in range(16)])
        bins,members,summary=rank_bins(f)
        q=bins[bins.pair.eq('H1_C')&bins.contrast.eq('PM')]
        np.testing.assert_allclose(q.y_mean,np.arange(8)*2+.5);np.testing.assert_allclose(q.y_se,.5)
        self.assertTrue(summary.rho.isna().all())
        self.assertEqual(members[members.pair.eq('H1_C')&members.contrast.eq('PM')].bin.tolist(),np.repeat(np.arange(8),2).tolist())
        bins,_,summary=rank_bins(f.iloc[:3]);self.assertTrue(bins.y_mean.isna().all())
        self.assertTrue(summary.status.eq('fewer_genes_than_planned_bins').all())

    def test_weighted_quantiles_and_histograms_keep_original_weights(self):
        f=pd.DataFrame(dict(group=['a']*3,value=[0.,.5,1.],weight=[1.,2.,1.]))
        points,hist,summary=distributions(f,['group'],'value','weight',np.linspace(0,1,13))
        np.testing.assert_allclose(points.cumulative,[.25,.75,1.]);self.assertEqual(hist.weight.sum(),4.)
        self.assertAlmostEqual(summary.q50.iloc[0],.25)
        f.value=1.;_,h,s=distributions(f,['group'],'value');self.assertEqual(h.weight.sum(),3);self.assertEqual(s.q50.iloc[0],1.)

    def test_empty_external_group_is_unknown_not_zero_and_influence_is_dynamic(self):
        f=pd.DataFrame(dict(gene_id=['a','a','b','b'],group=['somatic','PSC']*2,difference=[8,3,4,6]))
        np.testing.assert_allclose(grouped_difference(f,'somatic','PSC',4,3),[1,-1])
        self.assertEqual(len(grouped_difference(f,'somatic','PSC',0,3)),0)
        a=influence(pd.Series([2.,4.,9.],index=['a','b','c']));np.testing.assert_allclose(a.leave_one_out_mean,[6.5,5.5,3.])
        self.assertTrue(influence(pd.Series([7.],index=['a'])).leave_one_out_mean.isna().all())

    def test_no_eligible_windows_keeps_all_position_test_slots_without_zero_effect(self):
        with tempfile.TemporaryDirectory() as temp:
            tables={};ctx=SimpleNamespace(out=Path(temp),table=lambda f,n,**kw:tables.update({n:f}))
            state=dict(sets={sp:dict(barrier={'A'},ESCC={'B'},somatic={'A'},PSC={'B'}) for sp in ['human','mouse']},
                windows=pd.DataFrame(columns=['species','contrast','gene_id']),contexts={},hits={})
            q=position_tests(ctx,None,state)
            self.assertEqual(len(q),16);self.assertTrue(q.p.eq(1).all());self.assertTrue(q.effect.isna().all());self.assertTrue(q.status.eq('no_genes').all())

    def test_other_tables_cannot_mask_missing_required_statistic(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp);(p/'required-statistics-tables.tsv').write_text('path\trequired_columns\nexternal_mirna/primary_tests.tsv\tp;q_BH8\n')
            rows=[dict(path='tables/external_mirna/seed_universe.tsv',columns=['p','q_BH8'])]
            with self.assertRaisesRegex(ValueError,'Missing required'):require_numerical_products(p,rows)
            rows[0]['path']='tables/external_mirna/primary_tests.tsv';rows[0]['columns']=['p']
            with self.assertRaisesRegex(ValueError,'Missing statistical columns'):require_numerical_products(p,rows)
            rows[0]['columns'].append('q_BH8');self.assertEqual(len(require_numerical_products(p,rows)),1)


if __name__=='__main__':unittest.main()
