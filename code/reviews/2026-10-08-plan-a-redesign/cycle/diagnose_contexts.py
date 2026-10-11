"""Describe frozen gate failures and independently check raw two-marker support."""
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
CALLERS=['polyapipe','apadetect_p50'];SCORES=['Cycle_S','Cycle_G2M']

def archive():
    result={}
    with np.load(ROOT/'notes/2026-09-26-state-comparability-audit/expression_inputs.npz',allow_pickle=False) as p:
        for sid in ['GSE103221','GSE115943','GSE118258']:
            m=pd.DataFrame(json.loads(str(p[sid+'_meta_json'])));sp=m.species.iloc[0]
            symbols=p[sp+'_marker_symbols'].tolist();programs=json.loads(str(p[sp+'_program_json']))
            for score in SCORES:m[score]=p[sid+'_marker_logcp10k'][:,[symbols.index(g.upper()) for g in programs[score].split()]].mean(1)
            result[sid]=m
    with np.load(ROOT/'reviews/2026-10-08-plan-a-analysis/cycle/additional_cycle_scores.npz',allow_pickle=False) as p:
        for sid in ['GSE150311','GSE221739']:
            m=pd.DataFrame(json.loads(str(p[sid+'_meta_json'])))
            for score in SCORES:m[score]=p[sid+'_'+score]
            result[sid]=m
    return result

def main():
    meta=archive();summaries=[];support=[];failing=[];states=[]
    for caller in CALLERS:
        solvers=pd.read_csv(HERE/caller/'context_solvers.tsv',sep='\t').fillna({'source_library_key':''})
        pairs=pd.read_csv(HERE/caller/'context_pair_balance.tsv',sep='\t')
        state=pd.read_csv(HERE/caller/'context_state_balance.tsv',sep='\t')
        for key,q in solvers.groupby(['study_id','scope','method']):
            sid,scope,method=key;pr=pairs[pairs.study_id.eq(sid)&pairs.scope.eq(scope)&pairs.method.eq(method)];sr=state[state.study_id.eq(sid)&state.scope.eq(scope)&state.method.eq(method)]
            def maximum(suffix):
                vals=pr[[c for c in pr if c.endswith(suffix)]].to_numpy().ravel();vals=vals[~np.isnan(vals)]
                return float(np.max(vals)) if len(vals) else np.nan
            flags=[z for x in q.failures for z in json.loads(x)]
            summaries.append(dict(caller=caller,study_id=sid,scope=scope,method=method,n_contexts=len(q),n_balanced=int(q.status.eq('balanced').sum()),
                n_no_common_mean_support=int(q.solver_status.eq('no_common_mean_support').sum()),n_mean_QP_optimal=int(q.solver_status.eq('optimal').sum()),
                n_fixed_bins_fitted=int(q.solver_status.eq('fixed_bins').sum()),n_no_common_bin=int(q.solver_status.eq('no_common_2x2_bin').sum()),
                n_numerical_solver_failure=int(q.solver_status.isin(['linprog_failure','SLSQP_failure','zero_trim_constraint_failure']).sum()),
                n_KS_failure=sum(any(z.startswith('KS_') for z in json.loads(x)) for x in q.failures),
                n_variance_failure=sum(any(z.startswith('variance_ratio_') for z in json.loads(x)) for x in q.failures),
                n_ESS_failure=sum(any(z.startswith('insufficient_ESS_') for z in json.loads(x)) for x in q.failures),
                n_SMD_failure=sum(any(z.startswith('SMD_') for z in json.loads(x)) for x in q.failures),
                max_weighted_KS=maximum('_KS'),max_variance_ratio=maximum('_variance_ratio'),max_SMD=maximum('_SMD'),min_ESS=float(sr.ESS.min()) if len(sr) else np.nan,
                failures_counted_nonexclusive=True))
        for r in solvers.itertuples():
            for reason in json.loads(r.failures):failing.append(dict(caller=caller,study_id=r.study_id,route=r.route,scope=r.scope,method=r.method,context_id=r.context_id,reason=reason))
            if r.method!='continuous':continue
            m=meta[r.study_id];m=m[m.purity.ge(.8)&m.route.eq(r.route)];by={n:m[m.library_key.eq(lib)&m.node.eq(n)] for n,lib in [('P',r.library_key),('M',r.library_key)]+([('S',r.source_library_key)] if r.scope=='PMS' else [])}
            sizes=[len(q) for q in by.values()];off=np.cumsum([0]+sizes);names=list(by);A=[];b=[]
            for j,n in enumerate(names):
                row=np.zeros(off[-1]);row[off[j]:off[j+1]]=1;A.append(row);b.append(1.)
                q=by[n]
                for score in SCORES:states.append(dict(caller=caller,context_id=r.context_id,study_id=r.study_id,scope=r.scope,library_key=r.library_key,source_library_key=r.source_library_key,state=n,
                    score=score,n_mc=len(q),n_cells=int(q.n_cells.sum()),min=float(q[score].min()),max=float(q[score].max()),mean=float(q[score].mean()),variance=float(q[score].var(ddof=0))))
            intervals={score:min(q[score].max() for q in by.values())-max(q[score].min() for q in by.values()) for score in SCORES}
            for score in SCORES:
                for j in range(1,len(names)):
                    row=np.zeros(off[-1]);row[:off[1]]=by[names[0]][score];row[off[j]:off[j+1]]=-by[names[j]][score];A.append(row);b.append(0.)
            result=linprog(np.zeros(off[-1]),A_eq=np.asarray(A),b_eq=np.asarray(b),bounds=(0,None),method='highs')
            observed=r.solver_status!='no_common_mean_support'
            assert bool(result.success)==observed,(r.context_id,result.message,r.solver_status)
            support.append(dict(caller=caller,study_id=r.study_id,scope=r.scope,context_id=r.context_id,original_solver_status=r.solver_status,
                independent_raw_scale_LP_status=int(result.status),independent_mean_feasible=bool(result.success),
                S_marginal_interval_overlap=float(intervals['Cycle_S']),G2M_marginal_interval_overlap=float(intervals['Cycle_G2M']),
                simple_one_marker_nonoverlap=bool(min(intervals.values())<0),verified=True))
    for data,name in [(summaries,'context_failure_summary.tsv'),(failing,'context_failure_reasons.tsv'),(support,'independent_mean_support.tsv'),(states,'context_raw_marker_ranges.tsv')]:
        pd.DataFrame(data).to_csv(HERE/name,sep='\t',index=False,na_rep='NA')
    (HERE/'context_diagnostics.json').write_text(json.dumps(dict(status='passed',independent_raw_scale_support_checks=len(support),
        mean_feasible=sum(r['independent_mean_feasible'] for r in support),mean_infeasible=sum(not r['independent_mean_feasible'] for r in support),
        any_accepted_context=any(r['n_balanced']>0 for r in summaries),no_new_APA_or_p=True,no_new_adjustment_model=True),indent=2)+'\n')
    print(pd.DataFrame(summaries).to_string(index=False))

if __name__=='__main__':main()
