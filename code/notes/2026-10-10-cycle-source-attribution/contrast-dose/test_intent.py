"""Intent checks for contributor matching and date/unit weighting."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
spec=importlib.util.spec_from_file_location('analysis',Path(__file__).with_name('analyze.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Intent(unittest.TestCase):
    def test_source_unit_and_capture_counts_do_not_reweight_target_dates(self):
        q=pd.DataFrame([['D1','A','D0','S1'],['D1','A','D0','S1'],['D1','A','D0','S2'],['D2','B','D0','S1']],columns=['time_a','sample_unit_a','time_b','sample_unit_b'])
        np.testing.assert_allclose(m.weights(q),[.125,.125,.25,.5],rtol=0,atol=1e-14)
    def test_exact_gene_MC_differs_from_cohort_and_is_MC_equal(self):
        s=pd.DataFrame(dict(study_id=['study']*4,metacell_id=['p1','p2','s1','s2'],route=['route']*4,library_key=['P','P','S','S'],node=['P','P','S','S'],Cycle=[1.,5.,2.,8.],Adjusted_Cycle=[.2,2.,.4,3.]))
        p=pd.DataFrame(dict(study_id=['study'],route=['route'],contrast=['P_S'],library_key_a=['P'],library_key_b=['S'],metacell_ids_a=['["p1"]'],metacell_ids_b=['["s2"]']))
        out,_,_=m.cycle_pairs(p,s)
        self.assertEqual(out.exact_Cycle_delta.iloc[0],-7.)
        self.assertEqual(out.cohort_Cycle_delta.iloc[0],-2.)
        self.assertEqual(out.exact_Adjusted_Cycle_delta.iloc[0],-2.8)
if __name__=='__main__':unittest.main()
