"""Meaningful numerical examples for balance, depth selection and small-n inference."""
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('cycle_under_test',HERE/'run_cycle.py')
cycle=importlib.util.module_from_spec(spec);spec.loader.exec_module(cycle)

def main():
    x=np.array([0.,1.,2.,3.]);y=np.array([0.,1.,1.,3.]);w=np.full(4,.25)
    assert abs(cycle.weighted_ks(x,w,y,w)-ks_2samp(x,y).statistic)<1e-12
    a=dict(rows=np.arange(4),cells=np.full(4,10),bins=np.arange(4),Cycle_S=x,Cycle_G2M=np.array([0.,1.,0.,1.]))
    fitted,info=cycle.optimize({'P':a,'M':a});assert info['objective']<1e-12 and cycle.balance({'P':a,'M':a},fitted)[0]
    b=dict(a,Cycle_S=x+10);fitted,info=cycle.optimize({'P':a,'M':b});assert fitted is None and info['solver_status']=='no_common_mean_support'
    p=dict(a,Cycle_S=np.zeros(4),Cycle_G2M=np.ones(4));m=dict(p,Cycle_S=np.ones(4))
    fitted,info=cycle.optimize({'P':p,'M':m});assert fitted is None and 'inconsistent_full_equalities' in info['message']
    def gene_state(rows,score):
        c=np.full((3,2),10);return dict(rows=np.array(rows),cells=np.full(3,10),bins=np.arange(3),Cycle_S=np.array(score,float),Cycle_G2M=np.array([0.,1.,0.]),c=c,p=np.full((3,2),.5))
    # Different depth-qualified slices of initially equal states must fail re-balance.
    states={'P':gene_state([1,2,3],[1,2,3]),'M':gene_state([0,1,2],[0,1,2])}
    context=dict(passed=True,weights={'P':w,'M':w},sd={'PM_Cycle_S':x.std(),'PM_Cycle_G2M':.5})
    usage,status,audit,_,_=cycle.gene_adjust(states,context)
    assert usage is None and status=='gene_failed' and 'SMD_PM_Cycle_S' in json.loads(audit['failures'])
    assert cycle.inference([1.,2.],'toy')['p']==.5
    zero=cycle.inference([0.,2.],'toy');assert zero['p']==1 and zero['min_attainable_p']==1 and zero['generic_signflip_lower_bound']==.5
    result=dict(status='passed',examples=['independent_weighted_KS','uniform_QP_optimum','true_common_mean_nonoverlap','QR_constant_state_inconsistency',
        'depth_selection_requires_new_balance_gate','exact_signflip_two_libraries','zeros_change_actual_attainable_p'],not_claiming_real_gene_gates_passed=True)
    (HERE/'numeric-tests.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
