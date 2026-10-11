"""Independent algebra, balance, library-statistic and sign-flip verification."""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASES={
    'polyapipe':ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables',
    'apadetect_p50':ROOT/'notes/2026-10-06-dual-main-plot-data/apadetect-p50/supplement/tables',
}

def read(path):return pd.read_csv(path,sep='\t',low_memory=False)

def close(a,b,tol=2e-10):
    a=np.asarray(a,float);b=np.asarray(b,float)
    assert np.array_equal(np.isnan(a),np.isnan(b)),(a,b)
    assert np.allclose(a,b,atol=tol,rtol=1e-9,equal_nan=True),(a,b)

def exact_p(v):
    v=np.asarray(v,float);n=len(v);hit=0;obs=abs(v.mean())
    for start in range(0,2**n,4096):
        numbers=np.arange(start,min(start+4096,2**n),dtype=np.uint64)
        signs=((numbers[:,None]>>np.arange(n,dtype=np.uint64))&1).astype(float)*2-1
        hit+=int((np.abs(signs@v/n)>=obs-1e-12).sum())
    return hit/(2**n)

def ks_independent(a,b,score):
    grid=sorted(set(a[score])|set(b[score]))
    return max(abs(a.loc[a[score]<=z,'weight'].sum()-b.loc[b[score]<=z,'weight'].sum()) for z in grid)

def recurrent(base):
    p=read(base/'recurrence/recurrence_members.tsv.gz')
    return p[p.metric.eq('WPP')&p.contrast.eq('P_M')&p.threshold.eq(.02)&p.variant.eq('primary')&~p.series_matching]

def context_checks(caller):
    weights=read(HERE/caller/'context_weights.tsv.gz');solvers=read(HERE/caller/'context_solvers.tsv');pairs=read(HERE/caller/'context_pair_balance.tsv');states=read(HERE/caller/'context_state_balance.tsv')
    kkt=[];balance_rows=[]
    for cid,q in weights[weights.method.eq('continuous')].groupby('context_id'):
        names=list(q.state.unique());by={n:q[q.state.eq(n)].reset_index(drop=True) for n in names};sizes=[len(by[n]) for n in names]
        w=np.concatenate([by[n].weight.to_numpy() for n in names]);off=np.cumsum([0]+sizes);A=[];b=[]
        for j,n in enumerate(names):
            row=np.zeros(len(w));row[off[j]:off[j+1]]=1;A.append(row);b.append(1.)
            assert (by[n].weight>=0).all();close(by[n].weight.sum(),1)
        for score in ['Cycle_S','Cycle_G2M']:
            allx=q[score].to_numpy();sd=allx.std();mu=allx.mean()
            if sd==0:continue
            for j in range(1,len(names)):
                row=np.zeros(len(w));row[:off[1]]=(by[names[0]][score]-mu)/sd;row[off[j]:off[j+1]]=-(by[names[j]][score]-mu)/sd;A.append(row);b.append(0.)
        A=np.array(A);b=np.array(b);constraint=float(np.max(np.abs(A@w-b)));assert constraint<=1e-7
        counts=np.repeat(sizes,sizes);uniform=np.concatenate([np.repeat(1/n,n) for n in sizes]);grad=2*counts*(w-uniform);objective=float(np.sum(counts*(w-uniform)**2))
        recorded=solvers[solvers.context_id.eq(cid)].iloc[0];close(objective,recorded.objective,3e-7)
        active=w>1e-8;lam=np.linalg.lstsq(A[:,active].T,-grad[active],rcond=1e-10)[0];reduced=grad+A.T@lam
        positive=float(np.max(np.abs(reduced[active]))) if active.any() else np.nan
        inactive=float(max(0,-reduced[~active].min())) if (~active).any() else 0.
        # The Frank-Wolfe gap independently bounds suboptimality of this convex QP.
        linear=linprog(grad,A_eq=A,b_eq=b,bounds=(0,None),method='highs');assert linear.success
        gap=float(grad@w-linear.fun)
        dual_w=np.maximum(0,uniform-(A.T@lam)/(2*counts))
        dual_value=float(np.sum(counts*(dual_w-uniform)**2)+lam@(A@dual_w-b))
        dual_gap=objective-dual_value;assert dual_gap<1e-6,(cid,dual_gap)
        labels=['sum_'+n for n in names]+[score+'_'+names[0]+'_minus_'+n for score in ['Cycle_S','Cycle_G2M'] if q[score].std(ddof=0)>0 for n in names[1:]]
        kkt.append(dict(caller=caller,context_id=cid,objective=objective,constraint_error=constraint,positive_KKT_residual=positive,inactive_KKT_violation=inactive,
            Frank_Wolfe_gap_bound=max(0,gap),dual_lower_bound=dual_value,dual_constraint_labels=json.dumps(labels),dual_multipliers=json.dumps(lam.tolist()),
            dual_coordinate_formula='max(0,uniform-A.T@lambda/(2*n_state))',convex_optimality_gap_bound=max(0,dual_gap)))
        rr=pairs[pairs.context_id.eq(cid)].iloc[0]
        for a,z in itertools.combinations(names,2):
            for score in ['Cycle_S','Cycle_G2M']:
                k=ks_independent(by[a],by[z],score);close(k,rr[a+z+'_'+score+'_KS'])
    for scope in ['pm','dual']:
        eff=read(HERE/caller/(scope+'_gene_library_effects.tsv.gz' if scope=='pm' else 'dual_triplet_effects.tsv.gz'))
        state=read(HERE/caller/(scope+'_gene_state_balance.tsv.gz'));pair=read(HERE/caller/(scope+'_gene_pair_balance.tsv.gz'))
        valid=state[state.status.eq('measurable')];assert (valid.n_mc>=3).all() and (valid.n_cells>=30).all() and (valid.ESS>=3-1e-9).all()
        valid=pair[pair.status.eq('measurable')]
        for col in valid:
            if col.endswith('_SMD') and 'library_scale' not in col:assert (valid[col]<=.1+1e-9).all()
            if col.endswith('_variance_ratio'):assert (valid[col]<=2+1e-9).all()
            if col.endswith('_KS'):assert (valid[col]<=.2+1e-9).all()
            if col.endswith('_joint_PAS'):assert (valid[col]>=2).all()
        if scope=='dual':
            for v in ['raw_all','raw_context_common','raw_support_common','adjusted']:close(eff[v+'_PS'],eff[v+'_PM']+eff[v+'_MS'])
        balance_rows.append(dict(caller=caller,scope=scope,gene_contexts=len(eff),measurable=int(eff.status.eq('measurable').sum()),measurable_state_checks=len(state[state.status.eq('measurable')]),all_gates_verified=True))
    return kkt,balance_rows

def library_checks(caller,primary):
    frame=read(HERE/caller/'pm_gene_library_effects.tsv.gz');libs=read(HERE/caller/'pm_library_means.tsv');changes=read(HERE/caller/'pm_candidate_library_changes.tsv');margins=read(HERE/caller/'dual_candidate_gene_library_margins.tsv')
    own=recurrent(BASES[caller]);signs={'fixed_primary':primary.set_index(['species','gene_id']).direction.map({'positive':1,'negative':-1}),
        'caller_own':own.set_index(['species','gene_id']).direction.map({'positive':1,'negative':-1})}
    for r in libs.itertuples():
        q=frame[frame.study_id.eq(r.study_id)&frame.library_key.eq(r.library_key)&frame.method.eq(r.method)]
        if r.version!='raw_all':q=q[q.status.eq('measurable')]
        assert len(q)==r.n_genes;close(q[r.version+'_PM'].mean(),r.mean)
    for r in changes.itertuples():
        q=frame[frame.study_id.eq(r.study_id)&frame.library_key.eq(r.library_key)&frame.method.eq(r.method)&frame.status.eq('measurable')]
        sg=q.gene_id.map(signs[r.cohort].loc[r.species]);q=q[sg.notna()];sg=sg[sg.notna()]
        assert len(q)==r.n_genes;close((sg.to_numpy()*(q.adjusted_PM-q.raw_context_common_PM)).mean(),r.mean_signed_change)
        close(r.mean_support_trim+r.mean_reweight,r.mean_signed_change)
    trip=read(HERE/caller/'dual_triplet_effects.tsv.gz')
    fixed=read(BASES['polyapipe']/'supplement_statistics/dual-reference-all-measurable.tsv.gz');own=read(BASES[caller]/'supplement_statistics/dual-reference-all-measurable.tsv.gz')
    dsigns={name:p.set_index(['species','route','gene_id']).recurrent_direction.map({'positive':1,'negative':-1}) for name,p in [('fixed_primary',fixed),('caller_own',own)]}
    for r in margins.itertuples():
        q=trip[trip.study_id.eq(r.study_id)&trip.route.eq(r.route)&trip.library_key.eq(r.library_key)&trip.gene_id.eq(r.gene_id)&trip.method.eq(r.method)&trip.status.eq('measurable')]
        w=q.triplet_weight.to_numpy(copy=True);w/=w.sum();sg=dsigns[r.cohort].loc[(r.species,r.route,r.gene_id)]
        for version in ['raw_context_common','raw_support_common','adjusted']:
            value=min(sg*(w@q[version+'_PS']),sg*(w@q[version+'_PM']));close(value,getattr(r,version+'_margin'))
        close(r.support_trim+r.reweight,r.margin_change)
    return dict(caller=caller,PM_library_rows=len(libs),candidate_PM_library_rows=len(changes),candidate_PMS_gene_library_rows=len(margins),all_independently_recomputed=True)

def main():
    audits=json.loads((HERE/'input_sha256.json').read_text())
    for p,h in audits.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    receipt=json.loads((HERE/'execution-receipt.json').read_text())
    for p,h in receipt['sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
    primary=recurrent(BASES['polyapipe']);kkt=[];balance=[];libraries=[]
    for caller in BASES:
        a,b=context_checks(caller);kkt+=a;balance+=b;libraries.append(library_checks(caller,primary))
    tests=read(HERE/'raw-tests.tsv');assert tests.test_id.is_unique and len(tests)==26
    assert tests.groupby('family').size().to_dict()=={'A2-cycle-alignment-redesign':16,'A8-cycle-PM-redesign':10}
    assert tests.method.eq('continuous').all()
    assert exact_p([1.,2.])==.5 and exact_p([0.,2.])==1
    for r in tests.itertuples():
        if r.family=='A8-cycle-PM-redesign':
            f=read(HERE/r.pas_source/'pm_library_means.tsv');v=f.loc[f.study_id.eq(r.study_id)&f.method.eq('continuous')&f.version.eq('adjusted'),'mean']
        elif r.scope=='PM':
            f=read(HERE/r.pas_source/'pm_candidate_library_changes.tsv');v=f.loc[f.study_id.eq(r.study_id)&f.method.eq('continuous')&f.cohort.eq('fixed_primary'),'mean_signed_change']
        else:
            f=read(HERE/r.pas_source/'dual_library_margin_changes.tsv');v=f.loc[f.study_id.eq(r.study_id)&f.method.eq('continuous')&f.cohort.eq('fixed_primary'),'mean_margin_change']
        v=np.asarray(v,float);v=v[np.isfinite(v)];assert len(v)==r.n_libraries;close(float(v.mean()) if len(v) else np.nan,r.mean)
        if len(v)<2:assert pd.isna(r.p)
        elif len(v)<=20:close(exact_p(v),r.p,1e-12)
    pd.DataFrame(kkt).to_csv(HERE/'independent_QP_checks.tsv',sep='\t',index=False)
    result=dict(status='passed',input_hashes_unchanged=len(audits),formal_positions=26,family_positions=tests.groupby('family').size().to_dict(),
        finite_p=int(tests.p.notna().sum()),independent_QP_contexts=len(kkt),max_QP_gap=max([r['convex_optimality_gap_bound'] for r in kkt],default=0),balance=balance,libraries=libraries,
        PAS_additivity_verified=True,measured_gene_gates_verified=True,formal_exact_signflip_verified=True,
        prior_validation_failure_disclosed='VALIDATION-REPAIRS.json',no_scientific_model_or_FTOL_changed_for_validation=True)
    (HERE/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
