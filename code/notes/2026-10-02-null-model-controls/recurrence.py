"""Nested controls and conditional label-injection power for DWUI recurrence."""
import importlib.util
import numpy as np
import pandas as pd
from common import *

spec=importlib.util.spec_from_file_location('frozen_null',ROOT/'notes/2026-09-30-multimetric-statistics/nulls/analyze.py')
N=importlib.util.module_from_spec(spec);spec.loader.exec_module(N)
BASE=['eligible_units','series_pattern','can_negative','can_positive']
MEAS=['mean_GEX_bin','APA_coverage_bin']
MODELS={'L0_feasible_series':BASE,'L1_measurement':BASE+MEAS,
    'L2_baseline':BASE+MEAS+['baseline_bin'],
    'L2_PAS_first':BASE+MEAS+['PAS_number_bin'],
    'L3_full':BASE+MEAS+['baseline_bin','PAS_number_bin']}

def simulate(units,ng,n,random,strength=0):
    result=[]
    for start in range(0,n,250):
        size=min(250,n-start)
        plus=np.zeros((size,ng),np.int8);minus=np.zeros_like(plus)
        if strength:
            latent_negative=random.normal(0,strength,(size,ng))
            latent_positive=random.normal(0,strength,(size,ng))
        for f,gi,cells in units:
            labels=f.event.to_numpy(np.int8)
            draw=np.broadcast_to(labels,(size,len(f))).copy()
            for ii in cells:
                lab=labels[ii];neg=int((lab==-1).sum());pos=int((lab==1).sum())
                if max(neg,pos,len(ii)-neg-pos)==len(ii):continue
                if not strength:
                    draw[:,ii]=random.permuted(np.broadcast_to(lab,(size,len(ii))),axis=1)
                else:
                    d=np.zeros((size,len(ii)),np.int8)
                    if neg:
                        scores=latent_negative[:,gi[ii]]+random.gumbel(size=(size,len(ii)))
                        ix=np.argpartition(scores,-neg,axis=1)[:,-neg:]
                        np.put_along_axis(d,ix,-1,axis=1)
                    if pos:
                        scores=latent_positive[:,gi[ii]]+random.gumbel(size=(size,len(ii)))
                        scores[d==-1]=-np.inf
                        ix=np.argpartition(scores,-pos,axis=1)[:,-pos:]
                        np.put_along_axis(d,ix,1,axis=1)
                    assert ((d==-1).sum(1)==neg).all() and ((d==1).sum(1)==pos).all()
                    draw[:,ii]=d
            assert not ((draw==-1)&~f.can_negative.to_numpy()).any()
            assert not ((draw==1)&~f.can_positive.to_numpy()).any()
            plus[:,gi]+=draw==1;minus[:,gi]+=draw==-1
        result.append((np.maximum(plus,minus)>=2).sum(1))
    return np.concatenate(result)

def main():
    source=read(ROOT/'notes/2026-10-02-reviewer-followup/nulls/native_study_events.tsv.gz')
    old=read(ROOT/'notes/2026-10-02-reviewer-followup/nulls/recurrence_summary.tsv')
    summaries=[];diagnostics=[];power=[];replay=[];frozen=[];archives={}
    for key,f in source[source.metric.eq('WPP')].groupby(['species','contrast','threshold']):
        context=dict(zip(['species','contrast','threshold'],key));genes=sorted(f.gene_id.unique());index={g:i for i,g in enumerate(genes)}
        observed_labels=f.pivot(index='gene_id',columns='unit',values='event')
        observed=int(np.maximum(observed_labels.eq(1).sum(1),observed_labels.eq(-1).sum(1)).ge(2).sum())
        prepared=[]
        for unit,sub in f.groupby('unit'):
            sf=N.stratify(sub,'WPP',key[2]);sf['unit']=unit
            # Preserve the exact old strata and verify the full partition.
            sf['legacy_stratum']=sf.stratum+'|'+sf.series_pattern
            new=sf[MODELS['L3_full']].astype(str).agg('|'.join,axis=1)
            assert pd.crosstab(sf.legacy_stratum,new).gt(0).sum(1).eq(1).all()
            assert pd.crosstab(sf.legacy_stratum,new).gt(0).sum(0).eq(1).all()
            prepared.append((sf,np.array([index[g] for g in sf.gene_id])))
            frozen.append(sf.assign(**context))
        nulls={};units_by_model={}
        for model,cols in MODELS.items():
            units=[];frozen_labels=0;changed=0.;n_strata=0
            for sf,gi in prepared:
                cells=groups(sf,cols);units.append((sf,gi,cells));n_strata+=len(cells)
                for j,ii in enumerate(cells):
                    freq=sf.event.iloc[ii].value_counts().reindex([-1,0,1],fill_value=0).to_numpy()
                    moving=1-(freq/len(ii))@(freq/len(ii));changed+=len(ii)*moving
                    frozen_labels+=len(ii)*(np.count_nonzero(freq)==1)
                    diagnostics.append(dict(**context,model=model,unit=sf.unit.iloc[0],stratum=j,n=len(ii),
                        negative=freq[0],zero=freq[1],positive=freq[2],expected_changed=len(ii)*moving))
            units_by_model[model]=units
            null=simulate(units,len(genes),10000,rng('recurrence',*key,model));nulls[model]=null
            archives['|'.join(map(str,(*key,model)))]=null
            row=dict(**context,model=model,n_genes=len(genes),n_study_genes=len(f),n_studies=f.unit.nunique(),observed=observed,
                expected=null.mean(),null_sd=null.std(ddof=1),excess=observed-null.mean(),p=probability(null,observed),
                null_lo=np.quantile(null,.025),null_hi=np.quantile(null,.975),strata=n_strata,
                frozen_label_fraction=frozen_labels/len(f),expected_changed_fraction=changed/len(f))
            summaries.append(row)
            print('RECURRENCE',key,model,observed,round(null.mean(),3),row['p'],flush=True)
        before=old[old.variant.eq('native')&old.metric.eq('WPP')&old.series_matching&old.species.eq(key[0])&old.contrast.eq(key[1])&old.threshold.eq(key[2])].iloc[0]
        full=nulls['L3_full'];se=np.sqrt((full.var(ddof=1)+before.null_sd**2)/10000)
        assert observed==before.observed and abs(full.mean()-before.null_mean)<6*se
        replay.append(dict(**context,observed=observed,old_mean=before.null_mean,new_mean=full.mean(),difference_se=(full.mean()-before.null_mean)/se))
        for strength in [0,.5,1,2]:
            sim=simulate(units_by_model['L3_full'],len(genes),1000,rng('recurrence_inject',*key,strength),strength)
            archives['|'.join(map(str,(*key,'inject',strength)))]=sim
            for model,reference in nulls.items():
                p=probability(reference,sim)
                for alpha in [.05,.05/8]:
                    power.append(dict(**context,generator='L3_full',tested_model=model,shared_logweight_sd=strength,alpha=alpha,
                        simulated_mean=sim.mean(),excess_over_full_null=sim.mean()-full.mean(),
                        **rate(p<=alpha)))
    s=pd.DataFrame(summaries);s['q_eight']=s.groupby('model').p.transform(bh);s['q_all_40']=bh(s.p)
    save(s,'recurrence_ladder.tsv');save(pd.DataFrame(diagnostics),'recurrence_strata.tsv.gz')
    save(pd.concat(frozen),'recurrence_frozen_covariates.tsv.gz');save(pd.DataFrame(power),'recurrence_power.tsv')
    save(pd.DataFrame(replay),'recurrence_replay.tsv');np.savez_compressed(HERE/'recurrence_draws.npz',**archives)
    manifest('recurrence_manifest.json',status='PASS',reference_draws=10000,injection_draws=1000)

if __name__=='__main__':main()
