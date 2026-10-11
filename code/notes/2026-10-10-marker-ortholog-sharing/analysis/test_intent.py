"""Intent fixtures for unchanged voting and prespecified correction opportunities."""
import importlib.util,unittest
from pathlib import Path
import numpy as np,pandas as pd
p=Path(__file__).resolve().parent/'analyze.py';s=importlib.util.spec_from_file_location('new_marker_production',p);M=importlib.util.module_from_spec(s);s.loader.exec_module(M)

class Intent(unittest.TestCase):
 def test_mean_above_threshold_is_not_all_series_event(self):
  f=pd.DataFrame({'effect_min':[.019],'effect_max':[.03],'capacity_negative':[.5],'capacity_positive':[.5]})
  self.assertEqual(M.S.threshold_events(f,.02).event.iloc[0],0)
 def test_threshold_exact_not_capacity_tolerance(self):
  f=pd.DataFrame({'effect_min':[-.0200000001],'effect_max':[-.0199999999],'capacity_negative':[.0199999999],'capacity_positive':[.5]})
  self.assertEqual(M.S.threshold_events(f,.02).event.iloc[0],0)
  self.assertTrue(M.S.threshold_events(f,.02).can_negative.iloc[0])
 def test_two_study_recurrence_allows_opposing_third_study(self):
  f=pd.DataFrame({'gene_id':['g']*3,'gene_symbol':['G']*3,'unit':['a','b','c'],'event':[-1,-1,1]})
  r=M.A.observed_members(f);self.assertEqual(r.direction.iloc[0],'proximal');self.assertEqual(r.distal_support.iloc[0],1)
 def test_all_prespecified_slots_and_PS_once(self):
  rows=[]
  definitions=['original_matched','human_Mstar','both_Mstar'];thresholds=[.01,.015,.02,.025,.03,.035,.04,.045,.05]
  for d in definitions:
   for c in ['P_S','P_M']:
    for t in thresholds:
     for m in ['measurement_PAS_capacity','full_assay']:
      for h in ['proximal','distal']:
       for mouse in ['proximal','distal']:
        special=d=='human_Mstar' and c=='P_M' and t==.02 and m=='full_assay' and h==mouse=='proximal'
        rows.append(dict(definition=d,contrast=c,threshold=t,model=m,human_direction=h,mouse_direction=mouse,tested=True,p_for_correction=.001 if special else 1.))
  out=M.corrections(pd.DataFrame(rows));hit=out[out.p_for_correction.eq(.001)].iloc[0]
  self.assertAlmostEqual(hit.q8,.008);self.assertAlmostEqual(hit.q72,.072);self.assertAlmostEqual(hit.q144,.144);self.assertAlmostEqual(hit.q288,.288)
  self.assertTrue(out[out.contrast.eq('P_S')].groupby(['threshold','model','human_direction','mouse_direction']).q288.nunique().eq(1).all())

if __name__=='__main__':unittest.main(verbosity=2)
