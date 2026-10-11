"""Read-only numerical audit; synthetic examples, frozen metadata, saved outputs."""
import hashlib
import importlib.util
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent
OTHER=HERE.parent/'cycle'
ROOT=HERE.parents[2]
READS={}


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def register(path):
    READS[str(path)]=sha(path);return path


def read(path):
    return pd.read_csv(register(path),sep='\t',float_precision='round_trip',low_memory=False)


def constraints(parts,scaled=False):
    sizes=[len(q) for q in parts];off=np.cumsum([0]+sizes);A=[];b=[]
    for j in range(len(parts)):
        row=np.zeros(off[-1]);row[off[j]:off[j+1]]=1;A.append(row);b.append(1.)
    for score in ['Cycle_S','Cycle_G2M']:
        xx=np.concatenate([q[score] for q in parts]);center=xx.mean() if scaled else 0.;sd=xx.std() if scaled else 1.
        if sd==0:continue
        for j in range(1,len(parts)):
            row=np.zeros(off[-1]);row[:off[1]]=(parts[0][score]-center)/sd
            row[off[j]:off[j+1]]=-(parts[j][score]-center)/sd;A.append(row);b.append(0.)
    return np.array(A),np.array(b),np.repeat(sizes,sizes)


def toys(module):
    ks_cases=sign_cases=0
    arrays=list(itertools.product([0.,1.,2.],repeat=2));weights=[np.array([k/4,1-k/4]) for k in range(5)]
    for x0,y0 in itertools.product(arrays,repeat=2):
        x,y=np.array(x0),np.array(y0);grid=np.unique(np.r_[x,y])
        for wx,wy in itertools.product(weights,repeat=2):
            expected=max(abs(wx[x<=t].sum()-wy[y<=t].sum()) for t in grid)
            assert abs(module.weighted_ks(x,wx,y,wy)-expected)<1e-14;ks_cases+=1
    cases=[v for n in [2,3,4] for v in itertools.product([-1,0,1],repeat=n)]
    cases += [(1,2,3,4,5),(-1,0,2,0,-3),(1,2,3,4,5,6),(0,0,0,0,0,0)]
    for values in cases:
        n=len(values);sums=[sum(a*b for a,b in zip(sign,values)) for sign in itertools.product([-1,1],repeat=n)]
        expected=sum(abs(z)>=abs(sum(values)) for z in sums)/len(sums)
        minimum=sum(abs(z)==max(map(abs,sums)) for z in sums)/len(sums)
        result=module.inference(values,'independent_synthetic_audit',str(values))
        assert abs(result['p']-expected)<1e-14 and abs(result['min_attainable_p']-minimum)<1e-14;sign_cases+=1
    for values in [[],[0],[1]]:
        assert np.isnan(module.inference(values,'toy_na')['p'])
    qp_cases=[]
    for pp,mm in [([0,1,2],[1,2,3]),([0,0,0],[0,1,2]),([0,1,2],[3,4,5]),([0,0,0],[1,1,1])]:
        raw={name:dict(rows=np.arange(3),cells=np.full(3,10),Cycle_S=np.array(x,float),Cycle_G2M=2*np.array(x,float)+1)
             for name,x in [('P',pp),('M',mm)]}
        parts=[pd.DataFrame({s:raw[name][s] for s in ['Cycle_S','Cycle_G2M']}) for name in raw]
        A,b,d=constraints(parts);u=1/d;best=np.inf
        for mask in range(1,2**len(u)):
            free=np.flatnonzero([(mask>>i)&1 for i in range(len(u))]);af=A[:,free]
            adjustment=np.linalg.pinv((af/d[free])@af.T,rcond=1e-12)@(b-af@u[free])
            v=np.zeros(len(u));v[free]=u[free]+af.T@adjustment/d[free]
            if (v>=-1e-10).all() and np.max(np.abs(A@v-b))<1e-9:
                best=min(best,float(np.sum(d*(v-u)**2)))
        got,info=module.optimize(raw)
        if np.isfinite(best):
            assert got is not None;w=np.concatenate(list(got.values()))
            assert abs(np.sum(d*(w-u)**2)-best)<1e-8 and np.max(np.abs(A@w-b))<1e-8
        else:assert got is None and info['solver_status']=='no_common_mean_support'
        qp_cases.append(dict(P=pp,M=mm,exact_active_face_minimum=best if np.isfinite(best) else None,solver_status=info['solver_status']))
    return dict(weighted_KS_cases=ks_cases,exact_signflip_cases=sign_cases,NA_cases=3,QP_active_face_cases=qp_cases)


def main():
    for path in [OTHER/'PROTOCOL.md',OTHER.parent/'CYCLE-FREEZE-ADDENDUM.md',OTHER/'run_cycle.py',OTHER/'execution-receipt.json',OTHER/'IMPLEMENTATION-REPAIRS.json']:register(path)
    execution=json.loads((OTHER/'execution-receipt.json').read_text())
    assert all(sha(p)==h for p,h in execution['sha256'].items())
    spec=importlib.util.spec_from_file_location('audited_cycle',OTHER/'run_cycle.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    toy_results=toys(module);archive=module.cycle_archive();READS.update(module.INPUTS)
    original=register(ROOT/'reviews/2026-10-08-plan-a-analysis/cycle/pas_subsets.npz');additional=register(OTHER/'additional_pas_subsets.npz')
    support_rows=[];fit_rows=[];effect_rows=[];gate_rows=[]
    with np.load(original,allow_pickle=False) as old,np.load(additional,allow_pickle=False) as new:
        for caller in ['polyapipe','apadetect_p50']:
            solvers=read(OTHER/caller/'context_solvers.tsv');weights=read(OTHER/caller/'context_weights.tsv.gz')
            pairs=read(OTHER/caller/'context_pair_balance.tsv');states=read(OTHER/caller/'context_state_balance.tsv')
            metas={}
            for sid in module.STUDIES:
                pack=old if sid in module.STUDIES[:3] else new
                meta=pd.DataFrame(json.loads(str(pack[caller+'__'+sid+'__meta_json'])))
                meta=meta.merge(archive[sid][['metacell_id','Cycle_S','Cycle_G2M']],on='metacell_id',validate='one_to_one')
                metas[sid]=meta[meta.purity>=.8]
            for row in solvers.itertuples():
                f=metas[row.study_id];source='' if pd.isna(row.source_library_key) else row.source_library_key
                definitions=[('P',row.library_key),('M',row.library_key)]+([('S',source)] if row.scope=='PMS' else [])
                parts=[f[(f.route==row.route)&(f.library_key==lib)&(f.node==name)] for name,lib in definitions]
                assert sum(map(len,parts))==row.n_metadata
                if row.method=='continuous':
                    A,b,d=constraints(parts);rank=np.linalg.matrix_rank(A,tol=1e-10);aug=np.linalg.matrix_rank(np.column_stack([A,b]),tol=1e-10)
                    lp=linprog(np.zeros(len(d)),A_eq=A,b_eq=b,bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-9})
                    is_feasible=lp.success
                    assert is_feasible==(row.solver_status=='optimal')
                    certificates=[s for s in ['Cycle_S','Cycle_G2M'] if max(q[s].min() for q in parts)>min(q[s].max() for q in parts)]
                    if not is_feasible:assert lp.status==2 and certificates
                    support_rows.append(dict(caller=caller,context_id=row.context_id,scope=row.scope,full_LP_feasible=is_feasible,
                        rank_A=rank,rank_Ab=aug,full_b_consistent=rank==aug,interval_separation_scores=';'.join(certificates)))
                q=weights[weights.context_id==row.context_id]
                if len(q)==0:continue
                ordered=[q[q.state==name] for name,_ in definitions]
                for part,saved in zip(parts,ordered):assert set(part.metacell_id)==set(saved.metacell_id)
                ww={name:z.weight.to_numpy().copy() for (name,_),z in zip(definitions,ordered)}
                stats={};fail=[]
                for name,z in zip([v[0] for v in definitions],ordered):
                    w=ww[name];w[w<=1e-8]=0;w/=w.sum();positive=w>0
                    ess=1/(w@w);nc=int(z.n_cells.to_numpy()[positive].sum())
                    if positive.sum()<3:fail.append('insufficient_metacells_'+name)
                    if nc<30:fail.append('insufficient_cells_'+name)
                    if ess<3-1e-10:fail.append('insufficient_ESS_'+name)
                    stats[name]={}
                    for score in ['Cycle_S','Cycle_G2M']:
                        x=z[score].to_numpy();mu=w@x;var=w@((x-mu)**2);stats[name][score]=(mu,var,x.var())
                recorded=pairs[pairs.context_id==row.context_id].iloc[0]
                max_ks_error=0.
                for (a,_),(b,_) in itertools.combinations(definitions,2):
                    qa=q[q.state==a];qb=q[q.state==b]
                    for score in ['Cycle_S','Cycle_G2M']:
                        x=qa[score].to_numpy();y=qb[score].to_numpy();wa,wb=ww[a],ww[b]
                        ks=max(abs(wa[x<=t].sum()-wb[y<=t].sum()) for t in np.unique(np.r_[x,y]))
                        mx,vx,rx=stats[a][score];my,vy,ry=stats[b][score];sd=np.sqrt((rx+ry)/2)
                        smd=abs(mx-my)/sd if sd>0 else 0. if abs(mx-my)<=1e-12 else np.inf
                        ratio=max(vx,vy)/min(vx,vy) if min(vx,vy)>1e-14 else 1. if max(vx,vy)<=1e-14 and abs(mx-my)<=1e-12 else np.inf
                        key=a+b+'_'+score;max_ks_error=max(max_ks_error,abs(ks-recorded[key+'_KS']))
                        assert np.isclose(ks,recorded[key+'_KS'],atol=1e-10)
                        if not np.isfinite(smd) or smd>.1+1e-10:fail.append('SMD_'+key)
                        if not np.isfinite(ratio) or ratio>2+1e-10:fail.append('variance_ratio_'+key)
                        if ks>.2+1e-10:fail.append('KS_'+key)
                assert sorted(set(fail))==json.loads(row.failures)
                gate_rows.append(dict(caller=caller,context_id=row.context_id,status=row.status,reconstructed_failures=json.dumps(sorted(set(fail))),max_KS_error=max_ks_error))
                if row.method=='continuous':
                    A,b,d=constraints(ordered,scaled=True);w=np.concatenate([z.weight for z in ordered]);u=1/d;g=2*d*(w-u);positive=w>1e-8
                    lam=np.linalg.lstsq(A[:,positive].T,-g[positive],rcond=1e-10)[0];dw=np.maximum(0,u-A.T@lam/(2*d))
                    obj=float(np.sum(d*(w-u)**2));dual=float(np.sum(d*(dw-u)**2)+lam@(A@dw-b));res=g+A.T@lam
                    eq=float(np.max(np.abs(A@w-b)));gap=obj-dual
                    assert w.min()>=0 and eq<=1e-8 and abs(gap)<1e-7
                    fit_rows.append(dict(caller=caller,context_id=row.context_id,objective=obj,objective_minus_dual=gap,
                        max_full_equality_error=eq,positive_stationarity_max=float(np.abs(res[positive]).max()),
                        active_negative_violation=float(max(0,-res[~positive].min())) if (~positive).any() else 0.))
            pm=read(OTHER/caller/'pm_gene_library_effects.tsv.gz');dual=read(OTHER/caller/'dual_triplet_effects.tsv.gz')
            assert not pm.status.eq('measurable').any() and not dual.status.eq('measurable').any()
            for f in [pm,dual]:
                assert f.filter(regex='^(raw_context_common|raw_support_common|adjusted)_').isna().all().all()
            error=float((dual.raw_all_PS-dual.raw_all_PM-dual.raw_all_MS).abs().max());assert error<1e-10
            for name in ['pm_candidate_library_changes.tsv','dual_candidate_gene_library_margins.tsv','dual_library_margin_changes.tsv']:
                empty=read(OTHER/caller/name);assert len(empty)==0 and len(empty.columns)>0
            effect_rows.append(dict(caller=caller,PM_raw_rows=len(pm),PMS_raw_rows=len(dual),measurable_rows=0,raw_all_PMS_additivity_error=error,
                raw_context_common_additivity_observable_rows=0,raw_support_common_additivity_observable_rows=0,adjusted_additivity_observable_rows=0))
    tests=read(OTHER/'raw-tests.tsv');assert tests.groupby('family').size().to_dict()=={'A2-cycle-alignment-redesign':16,'A8-cycle-PM-redesign':10}
    assert tests.test_id.is_unique and len(tests)==26 and tests.n_libraries.eq(0).all() and tests.p.isna().all()
    assert tests[['mean','ci_low','ci_high']].isna().all().all()
    for name,rows in [('cycle-support-audit.tsv',support_rows),('cycle-QP-audit.tsv',fit_rows),('cycle-gate-audit.tsv',gate_rows),('cycle-effects-audit.tsv',effect_rows)]:
        pd.DataFrame(rows).to_csv(HERE/name,sep='\t',index=False,float_format='%.17g')
    assert all(sha(path)==fingerprint for path,fingerprint in READS.items())
    result=dict(status='PASS',utc=datetime.now(timezone.utc).isoformat(),audited_final_cycle_code_sha256=sha(OTHER/'run_cycle.py'),
        toy_tests=toy_results,continuous_contexts=len(support_rows),full_LP_feasible=sum(r['full_LP_feasible'] for r in support_rows),
        full_LP_infeasible=sum(not r['full_LP_feasible'] for r in support_rows),infeasible_interval_certificates_complete=True,
        QP_fits_checked=len(fit_rows),max_dual_gap=max(r['objective_minus_dual'] for r in fit_rows),
        max_equality_error=max(r['max_full_equality_error'] for r in fit_rows),distribution_gates_checked=len(gate_rows),
        planned_positions=26,valid_p_values=0,NA_positions=26,adjusted_additivity_observable_rows=0,
        scope_limit='No new models, APA effects, PAS rerun, extraction or scientific p; missing adjusted effects explicitly unavailable.',
        other_module_SHA_unchanged=True,read_sha256=READS,audit_code_sha256=sha(__file__))
    (HERE/'cycle-numerical-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
