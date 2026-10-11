import unittest
import numpy as np
import pandas as pd
from analyze import association,aggregate_contexts,library_states
class IntentTests(unittest.TestCase):
 def test_fixed_gene_panel_prevents_composition_artifact(self):
  # A distal gene appears only in a low-cycle MC. Varying rosters invent shortening.
  observed=np.array([[.2,.8],[.2,np.nan]])
  self.assertGreater(np.nanmean(observed[0]),np.nanmean(observed[1]))
  panel=np.isfinite(observed).all(0)
  self.assertTrue(np.array_equal(observed[:,panel].mean(1),[.2,.2]))
 def test_covariates_remove_read_depth_confounding(self):
  # Both scores covary with depth; independent orthogonal residuals contain no relation.
  n=100;t=np.arange(n)*2*np.pi/n;depth=np.linspace(-1,1,n)
  c=depth[:,None];base=np.column_stack([np.ones(n),depth])
  noise1=np.sin(t)-base@np.linalg.lstsq(base,np.sin(t),rcond=None)[0]
  noise2=np.cos(t)-base@np.linalg.lstsq(base,np.cos(t),rcond=None)[0]
  noise2-=noise1*(noise1@noise2)/(noise1@noise1)
  r=association(10*depth+noise1,10*depth+noise2,c)
  self.assertGreater(r['correlation_unadjusted'],.95)
  self.assertLess(abs(r['partial_correlation']),1e-10)
 def test_library_weight_not_MC_weight(self):
  rows=[]
  for lib,corr,states in [('large',.4,['S','M','P']),('small',-.4,['P'])]:
   for state in states:
    rows.append(dict(species='mouse',study_id='fixture',route='route',score='Cycle',APA_variant='mean_DWUI',state=state,library_key=lib,time='D1',status='available',partial_correlation=corr,beta_per_score_SD=corr,beta_per_score_unit=corr))
  route,library=aggregate_contexts(pd.DataFrame(rows));self.assertTrue(library.status.eq('descriptive').all());r=route[route.scope.eq('all_SMP')].iloc[0]
  self.assertAlmostEqual(r.partial_correlation,0)
  self.assertEqual(r.n_libraries,2)
 def test_available_library_estimates_are_not_labelled_unavailable(self):
  rows=[dict(species='mouse',study_id='fixture',route='route',score='Cycle',APA_variant='mean_DWUI',state='P',library_key='well1',time='D1',status='available',partial_correlation=.2,beta_per_score_SD=.1,beta_per_score_unit=.05)]
  route,library=aggregate_contexts(pd.DataFrame(rows))
  self.assertEqual(library.iloc[0].status,'descriptive')
  self.assertEqual(route[route.scope.eq('P')].iloc[0].status,'descriptive')
 def test_rank_is_effective_covariate_count(self):
  t=np.linspace(-1,1,30);r=association(np.sin(3*t),np.cos(2*t),np.column_stack([t,t,np.ones(30)]))
  self.assertEqual(r['covariate_rank'],2)
  self.assertEqual(r['residual_df'],27)
 def test_constant_residual_not_false_null(self):
  t=np.linspace(-1,1,30);r=association(t,t*t,t[:,None])
  self.assertEqual(r['status'],'constant_residual_score_or_APA')
  self.assertTrue(np.isnan(r['partial_correlation']))
 def test_gene_denominator_has_all_genes(self):
  # Cycle markers alone are not a valid library denominator.
  counts=np.array([10,10]);all_gene_totals=np.array([100,1000])
  score=np.log1p(counts/all_gene_totals*1e4)
  self.assertGreater(score[0],score[1])
if __name__=='__main__':unittest.main()
