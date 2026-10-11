"""Omitted studies and Monte Carlo precision must alter the correct universe."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from downstream_inference import mc_stability
from downstream_recurrence import numerical,omit_study,leave_one_study_nulls,cross_strata,THRESHOLDS
from downstream_runtime import Context,bh
from downstream_cohorts import complete_triplets


class InferenceCompletionTests(unittest.TestCase):
    def make_tests(self,p,exact,n=999):
        return pd.DataFrame(dict(gene_id=[str(i) for i in range(len(p))],gene_symbol='g',
            metric='DWUI',p=p,q_BH_joint=bh(p),exact_sign_enumeration=exact,sign_patterns=n))

    def test_MC_bounds_keep_primary_p_and_exact_tail_probabilities(self):
        original=self.make_tests([.001,.05,1.],[False,True,False]);before=original.copy(deep=True)
        out=mc_stability(original)
        pd.testing.assert_frame_equal(original,before)
        self.assertEqual(out.p_lower_simultaneous95.iloc[0],0)
        self.assertEqual(out.p_upper_simultaneous95.iloc[2],1)
        self.assertEqual(out.p_lower_simultaneous95.iloc[1],.05)
        self.assertEqual(out.p_upper_simultaneous95.iloc[1],.05)
        self.assertTrue(out.stable_joint_q05_reject.iloc[0])
        self.assertTrue(out.q_BH_joint_lower.le(out.q_BH_joint_upper).all())

    def test_larger_MC_family_requires_wider_simultaneous_bounds(self):
        one=mc_stability(self.make_tests([.01],[False]))
        two=mc_stability(self.make_tests([.01,.01],[False,False]))
        self.assertLess(two.p_lower_simultaneous95.iloc[0],one.p_lower_simultaneous95.iloc[0])
        self.assertGreater(two.p_upper_simultaneous95.iloc[0],one.p_upper_simultaneous95.iloc[0])
        exact=mc_stability(self.make_tests([.01,.2],[True,True]))
        np.testing.assert_array_equal(exact.q_BH_joint,exact.q_BH_joint_lower)
        np.testing.assert_array_equal(exact.q_BH_joint,exact.q_BH_joint_upper)

    def events(self):
        rows=[]
        for metric,cuts in THRESHOLDS.items():
            for cut in cuts:
                for unit,genes in [('A',['g','h','only_A']),('B',['g','h']),('C',['g'])]:
                    for gene in genes:
                        rows.append(dict(species='human',contrast='P_M',metric=metric,threshold=cut,
                            unit=unit,gene_id=gene,gene_symbol=gene,event=int(gene=='g'),
                            mean_GEX=1.,APA_coverage=10.,PAS_number=2,baseline=.5,
                            capacity_unsigned=1.,capacity_positive=.5,capacity_negative=.5,
                            n_routes=1,series_pattern='r'))
        out=pd.DataFrame(rows)
        out['eligible_units']=out.groupby(['metric','threshold','gene_id']).unit.transform('nunique')
        return out

    def test_two_qualified_endpoint_pairs_do_not_imply_a_qualified_MS_pair(self):
        base=dict(species='mouse',study_id='A',route='r',setting='supported50_p08',gene_id='g',gene_symbol='G',
            library_key_a='target',PAS_ids='["a", "b"]',n_mc_b=3,target_usage='[0.5, 0.5]',
            reference_usage='[1.0, 0.0]',molecules_b='[50, 0]',site_mc_b='[3, 0]')
        source=pd.DataFrame([base]);mesenchymal=pd.DataFrame([base])
        self.assertEqual(len(complete_triplets(source,mesenchymal)),0)
        source.loc[0,'molecules_b']='[45, 5]';source.loc[0,'site_mc_b']='[3, 1]'
        mesenchymal.loc[0,'molecules_b']='[45, 5]';mesenchymal.loc[0,'site_mc_b']='[3, 1]'
        out=complete_triplets(source,mesenchymal)
        self.assertEqual(len(out),1);self.assertEqual(out.MS_joint_PAS.iloc[0],2)

    def test_omission_rebuilds_eligibility_without_imputing_missing_genes(self):
        events=self.events();out=omit_study(events,'human','A')
        self.assertEqual(set(out.gene_id),{'g','h'})
        self.assertEqual(set(out[out.gene_id.eq('g')].eligible_units),{2})
        self.assertEqual(set(out[out.gene_id.eq('h')].eligible_units),{1})
        self.assertEqual(len(out),18)
        self.assertEqual(set(events[events.gene_id.eq('g')].eligible_units),{3})

    def test_LOSO_rebuilds_nulls_and_retains_planned_unavailable_tests(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'config.json'
            config.write_text(json.dumps(dict(project_root=str(Path(__file__).resolve().parents[2]),
                native_run=tmp,frozen_gex=tmp,analysis_root=tmp,run_id='fixture',input_fingerprint='fixture')))
            math=numerical(Context(config,'fixture',root/'out'))
            summary,strata,diagnostics,archive=leave_one_study_nulls(self.events(),math,n=40)
            self.assertEqual(len(summary),72)
            unavailable=summary[summary.contrast.eq('P_S')]
            self.assertTrue(unavailable.status.eq('omitted_study_not_evaluable').all())
            self.assertTrue(unavailable.empirical_p.eq(1).all())
            self.assertEqual(len(archive),36)
            self.assertTrue(all(len(draw)==40 for draw in archive.values()))
            st=pd.concat(strata)
            self.assertFalse(st.unit.eq(st.omitted_study).any())
            self.assertFalse(st[st.omitted_study.eq('A')].gene_id.eq('only_A').any())
            self.assertTrue(st[st.gene_id.eq('g')].eligible_units.eq(2).all())
            self.assertTrue(all(len(frame)>0 for frame in diagnostics))

    def test_range_sensitivity_separates_same_median_different_series_ranges(self):
        from types import SimpleNamespace
        frame=pd.DataFrame(dict(baseline_min=[.1,.2,1.],baseline_max=[.9,.8,1.],baseline=[.5,.5,1.],
            eligible_units=2,n_routes=2,mean_GEX_bin=0,APA_coverage_bin=0,PAS_number_bin=0,
            can_negative=True,can_positive=False,baseline_bin=0,stratum='original'))
        math=SimpleNamespace(stratify=lambda f,*a,**kw:f.copy())
        out=cross_strata(math,frame,'WPP',.02,'Range10')
        self.assertNotEqual(out.stratum.iloc[0],out.stratum.iloc[1])
        self.assertEqual(out.baseline_max_bin.iloc[2],9)
        self.assertEqual(out.baseline_min_bin.iloc[2],9)
        pd.testing.assert_frame_equal(cross_strata(math,frame,'WPP',.02,'Q5'),frame)


if __name__=='__main__':unittest.main()
