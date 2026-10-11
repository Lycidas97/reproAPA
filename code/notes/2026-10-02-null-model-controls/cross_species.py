"""Fixed ortholog cohort, identical independent-species shuffle mechanism."""
import numpy as np
import pandas as pd
from scipy.stats import nchypergeom_fisher, hypergeom
from common import *

MODELS = {
    'L0_measurement':['effect_n'],
    'L1_capacity':['effect_n','prox_capacity'],
    'L2_baseline':['effect_n','prox_capacity','ref_q5'],
    'L2_PAS_first':['effect_n','prox_capacity','pas_q3'],
    'L3_full':['effect_n','prox_capacity','ref_q5','pas_q3'],
}

def independent(x,y,hg,mg,n,random):
    result=[]
    for start in range(0,n,500):
        size=min(500,n-start)
        h=binary_draw(x,hg,size,random);m=binary_draw(y,mg,size,random)
        result.append((h&m).sum(1))
    return np.concatenate(result)

def injected(x,y,hg,mg,n,odds,random):
    # Aggregate overlap draws are sufficient; every species/stratum event count
    # remains fixed. Human labels vary independently before conditional tilt.
    out=[]
    for start in range(0,n,500):
        size=min(500,n-start);h=binary_draw(x,hg,size,random)
        overlap=np.zeros(size,dtype=int)
        for ii in mg:
            N=len(ii); draws=int(y[ii].sum()); K=h[:,ii].sum(1)
            for k in np.unique(K):
                use=K==k;lo=max(0,draws+int(k)-N);hi=min(draws,int(k))
                if lo==hi: overlap[use]+=lo
                else:
                    dist=hypergeom if odds==1 else nchypergeom_fisher
                    args=(N,int(k),draws) if odds==1 else (N,int(k),draws,odds)
                    overlap[use]+=dist.rvs(*args,size=int(use.sum()),random_state=random)
        out.append(overlap)
    return np.concatenate(out)

def main():
    f=read(ROOT/'notes/2026-09-29-core-sensitivity/figure4-matched-overlap/ortholog_matched_features.tsv.gz')
    old=read(ROOT/'notes/2026-09-29-core-sensitivity/figure4-matched-overlap/overlap_null_results.tsv')
    assert len(f)==16092 and not f.human_gene_id.duplicated().any() and not f.mouse_gene_id.duplicated().any()
    summary=[];balance=[];strata=[];power=[];archives={};replay=[]
    for contrast,expected_observed in [('PS',18),('PM',10)]:
        x=f[f'h_{contrast}'].to_numpy(bool); y=f[f'm_{contrast}'].to_numpy(bool)
        observed=int((x&y).sum());assert observed==expected_observed
        cells={};nulls={}
        for model,cols in MODELS.items():
            hs=groups(f,[f'h_{contrast}_{c}' for c in cols]);ms=groups(f,[f'm_{contrast}_{c}' for c in cols]);cells[model]=(hs,ms)
            hp,hd=mobility(x,hs);mp,md=mobility(y,ms)
            mu=float(hp@mp)
            null=independent(x,y,hs,ms,20000,rng('cross',contrast,model));nulls[model]=null
            archives[f'{contrast}|{model}']=null
            row=dict(contrast=contrast,model=model,n_orthologs=len(f),human_events=int(x.sum()),mouse_events=int(y.sum()),
                observed=observed,expected=mu,mc_mean=null.mean(),null_sd=null.std(ddof=1),excess=observed-mu,
                p=probability(null,observed),null_lo=np.quantile(null,.025),null_hi=np.quantile(null,.975),
                **{'human_'+k:v for k,v in hd.items()},**{'mouse_'+k:v for k,v in md.items()})
            summary.append(row)
            for side,lab,probs,ss in [('h',x,hp,hs),('m',y,mp,ms)]:
                for j,ii in enumerate(ss):
                    strata.append(dict(contrast=contrast,model=model,side=side,stratum=j,n=len(ii),events=int(lab[ii].sum())))
                for col in ['effect_n','prox_capacity','ref_median','pas_q3']:
                    vals=f[f'{side}_{contrast}_{col}'].to_numpy(float);valid=np.isfinite(vals)
                    selected=lab&valid; weighted=probs*valid
                    balance.append(dict(contrast=contrast,model=model,side=side,covariate=col,
                        selected_mean=vals[selected].mean(),null_selected_mean=np.nansum(probs*vals)/weighted.sum(),
                        selected_missing=int((lab&~valid).sum()),expected_missing=float(probs[~valid].sum())))
                required=2 if contrast=='PS' else 1
                row[f'{side}_expected_events_below_capacity']=float(probs[f[f'{side}_{contrast}_prox_capacity'].lt(required)].sum())
            print('CROSS',contrast,model,observed,mu,row['p'],flush=True)
        # Same complete model, covariates and shuffling as the earlier run.
        before=old[old.contrast.eq(contrast)&old.model.eq('independent_q5_reference_capacity')].iloc[0]
        after=[r for r in summary if r['contrast']==contrast and r['model']=='L3_full'][0]
        assert abs(before.analytic_null_mean-after['expected'])<1e-10
        replay.append(dict(contrast=contrast,old_expected=before.analytic_null_mean,new_expected=after['expected'],observed=observed))
        for odds in [1,2,4,8]:
            sim=injected(x,y,*cells['L3_full'],3000,odds,rng('cross_inject',contrast,odds))
            archives[f'{contrast}|injection_OR{odds}']=sim
            for model,reference in nulls.items():
                p=probability(reference,sim)
                for alpha in [.05,.025]:
                    power.append(dict(contrast=contrast,generator='L3_full',tested_model=model,odds_ratio=odds,alpha=alpha,
                        experiment='conditional_injection',simulated_mean=sim.mean(),simulated_sd=sim.std(ddof=1),**rate(p<=alpha)))
        # Separate null calibration for each model; not the reference draws.
        for model,(hs,ms) in cells.items():
            sim=independent(x,y,hs,ms,3000,rng('cross_self_calibration',contrast,model))
            power.append(dict(contrast=contrast,generator=model,tested_model=model,odds_ratio=1,alpha=.05,
                experiment='self_calibration',simulated_mean=sim.mean(),simulated_sd=sim.std(ddof=1),**rate(probability(nulls[model],sim)<=.05)))
    s=pd.DataFrame(summary);s['q_two_contrasts']=s.groupby('model').p.transform(bh);s['q_all_10']=bh(s.p)
    save(s,'cross_species_ladder.tsv');save(pd.DataFrame(balance),'cross_species_balance.tsv')
    save(pd.DataFrame(strata),'cross_species_strata.tsv.gz');save(pd.DataFrame(power),'cross_species_power.tsv')
    save(pd.DataFrame(replay),'cross_species_replay.tsv');np.savez_compressed(HERE/'cross_species_draws.npz',**archives)
    manifest('cross_species_manifest.json',status='PASS',reference_draws=20000,injection_draws=3000)

if __name__=='__main__':main()
