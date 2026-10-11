import itertools,unittest
import numpy as np
import pandas as pd
from downstream_replicate_diagnostics import confidence_intervals,six_well_indices,exact_vector_test,run

class ReplicateDiagnosticsTests(unittest.TestCase):
    def test_degenerate_interval_is_unavailable_without_changing_test(self):
        f=pd.DataFrame(dict(SE=[0,.1,np.nan],df=[np.nan,8,np.nan],CI95_low=[0,-.2,np.nan],CI95_high=[0,.2,np.nan],p=[1,.2,.3],effect=[0,.1,.2]))
        z=confidence_intervals(f)
        self.assertTrue(z.loc[0,['CI95_low','CI95_high']].isna().all())
        self.assertEqual(z.CI_status.tolist(),['degenerate_observed_between_unit_variance','estimable','not_provided'])
        pd.testing.assert_frame_equal(z[['p','effect']],f[['p','effect']])
        self.assertEqual(z.loc[1,'CI95_low'],-.2)
    def test_six_units_preserve_equal_date_estimand_and_full_vector_null(self):
        q=pd.DataFrame(dict(time_a=['D1','D1','D2','D2','D3','D3','D4','D4'],sample_unit_a=[str(i) for i in range(8)]))
        ix,w,dates=six_well_indices(q,dict(D1=2,D2=0,D3=3,D4=1),{str(i):7-i for i in range(8)})
        self.assertEqual(len(set(q.sample_unit_a.iloc[ix])),6)
        self.assertEqual(q.time_a.iloc[ix].value_counts().sort_values().tolist(),[1,1,2,2])
        np.testing.assert_allclose(pd.Series(w).groupby(q.time_a.iloc[ix].to_numpy()).sum(),.25)
        d=np.array([[.1,-.1],[-.2,.2],[.15,-.15],[.05,-.05],[-.3,.3],[.2,-.2]])
        linear=np.array([[0,.5,0,1],[1,1,1,0]])
        effect,p=exact_vector_test(d,w,linear)
        values=[]
        for pattern in itertools.product([-1,1],repeat=6):
            v=sum(ww*ss*x for ww,ss,x in zip(w,pattern,d));values.append([*abs(v@linear),v.max()-v.min()])
        expected=(np.array(values)>=np.r_[abs(effect[:4]),effect[4]]-1e-12).mean(0)
        np.testing.assert_array_equal(p,expected)
        self.assertTrue((p>=.03125).all())
    def test_current_cohort_empty_is_explicit_not_old_fallback(self):
        class Context:
            def __init__(self):
                from pathlib import Path
                self.analysis=Path('/unused');self.tables={}
            def read(self,path):
                if 'gene_metric' in str(path):return pd.DataFrame(columns=['gene_id','gene_symbol','metric','SE','df','CI95_low','CI95_high','effect','p','q_BH_joint'])
                if 'qualified_well' in str(path):return pd.DataFrame(columns=['gene_id'])
                return pd.DataFrame(columns=['study_id','primary_gene_id','atlas_pas_id'])
            def table(self,f,name,columns=None):self.tables[name]=f if isinstance(f,pd.DataFrame) else pd.DataFrame(f,columns=columns)
        c=Context();run(c)
        s=c.tables['well_six_unit_downsample_summary.tsv']
        self.assertTrue(s.status.eq('no_eligible_current_well_genes').all())
        self.assertTrue(s.genes.eq(0).all())

if __name__=='__main__':unittest.main()
