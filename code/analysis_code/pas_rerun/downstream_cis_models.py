"""Firth models with planned families, paired folds, and explicit unestimable fits."""
from itertools import product
import json

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize,brentq
from scipy.stats import chi2,rankdata

from downstream_runtime import bh,seed

SOURCE='github/reproapa-analysis/modules/09_cross_species/code/cis_supplement/scripts/recompute_models.py'
BASE=['abs_effect','disc_eligible','val_eligible','disc_pas']
REC=BASE+['disc_score']
CIS=['driver_conserved','frac_conserved_discovery_pas','A_g','phastcons_driver_mean','motif_rank_driver_mean']
FAMILIES={'M0':BASE,'M1':REC,'PAS':REC+CIS[:2],'Architecture':REC+['A_g'],'Sequence':REC+CIS[3:],'All cis':REC+CIS}
KEYS=['layer','contrast','discovery','tau']


def folds(y,random_seed,n=5):
    """One deterministic class-stratified partition shared by every model."""
    y=np.asarray(y);assert set(y)<={0,1} and min(np.bincount(y.astype(int),minlength=2))>=n
    rng=np.random.default_rng(random_seed);labels=np.empty(len(y),int)
    for value in [0,1]:
        ix=rng.permutation(np.flatnonzero(y==value));labels[ix]=np.arange(len(ix))%n
    return labels


def auc(y,p):
    y=np.asarray(y);p=np.asarray(p);n1=int(y.sum());n0=len(y)-n1
    if not n0 or not n1:return np.nan
    return float((rankdata(p)[y==1].sum()-n1*(n1+1)/2)/(n1*n0))


def score(y,p):
    p=np.clip(p,1e-12,1-1e-12)
    return dict(auc=auc(y,p),logloss=float(-np.mean(y*np.log(p)+(1-y)*np.log1p(-p))))


def design_status(frame,columns,math):
    if len(frame)<3:return 'fewer_than_three_genes',None
    if not np.isfinite(frame[columns].to_numpy(float)).any(axis=0).all():return 'all_missing_predictor',None
    X,_,names=math.design(frame,frame,columns)
    if np.linalg.matrix_rank(X)<X.shape[1]:return 'rank_deficient',None
    return 'estimable',(X,names)


def models(ctx,data):
    math=ctx.functions(SOURCE,['value_grad','fit','design','profile_ci','verify'],dict(expit=expit,minimize=minimize,brentq=brentq,chi2=chi2))
    qa=math.verify();results=[];coef=[];cv=[];foldrows=[];oof=[];nulls=[];intervals=[];maximum_score=0.
    for key in product(['supported50','expanded_primary'],['P-S','P-M'],['human','mouse'],[.01,.02,.03,.05]):
        info=dict(zip(KEYS,key));g=data
        for column,value in info.items():g=g[g[column].eq(value)]
        g=g.sort_values('human_gene_id').reset_index(drop=True);assert g.human_gene_id.is_unique
        y=g.Y.to_numpy(float);primary=key[0]=='supported50' and key[3]==.02
        status,designed=design_status(g,REC+CIS,math)
        row=dict(info,n=len(g),positive=int(y.sum()),status=status,p_boot=1.,bootstrap_draws=0,primary=primary)
        if designed is not None:
            X,names=designed;fixed={names.index(c)+1:0 for c in CIS if c in names}
            row.update(retained_terms=';'.join(names),dropped_constants=';'.join(c for c in REC+CIS if c not in names),df=len(fixed))
            if not fixed:row['status']='no_variable_cis_predictor'
            else:
                b,obj,err=math.fit(X,y);bn,on,errn=math.fit(X,y,fixed,b);maximum_score=max(maximum_score,err,errn)
                assert np.isfinite(obj) and obj<1e99 and np.isfinite(on) and on<1e99
                stat=max(0.,2*(on-obj));rng=np.random.default_rng(seed('cis_bootstrap',*key));null=[]
                for draw in range(4999):
                    yy=rng.binomial(1,expit(X@bn));bb,oo,ss=math.fit(X,yy,start=bn);nb,no,ns=math.fit(X,yy,fixed,start=bn)
                    maximum_score=max(maximum_score,ss,ns);null.append(max(0.,2*(no-oo)))
                row.update(status='estimated',lrt=stat,p_asymptotic=chi2.sf(stat,len(fixed)),p_boot=(1+int((np.array(null)>=stat).sum()))/5000,bootstrap_draws=4999)
                nulls.extend(dict(info,draw=i,lrt_null=v) for i,v in enumerate(null))
                for term in CIS:
                    r=dict(info,term=term,status='constant',beta=np.nan,ci_low=np.nan,ci_high=np.nan)
                    if term in names:
                        j=names.index(term)+1;lo,hi=math.profile_ci(X,y,b,obj,j) if primary else [np.nan,np.nan]
                        r.update(status='estimated',beta=b[j],ci_low=lo,ci_high=hi)
                    coef.append(r)
        results.append(row)
        families=FAMILIES if primary else {n:FAMILIES[n] for n in ['M1','All cis']}
        predictions={name:np.full((10,len(g)),np.nan) for name in families}
        cv_allowed=len(y)>0 and min(np.bincount(y.astype(int),minlength=2))>=5
        for repeat in range(10):
            labels=folds(y,2200+repeat) if cv_allowed else None
            if labels is not None:foldrows.extend(dict(info,repeat=repeat,human_gene_id=gid,fold=int(fold)) for gid,fold in zip(g.human_gene_id,labels))
            for family,columns in families.items():
                reason='fewer_than_five_genes_in_one_outcome_class' if labels is None else 'estimated'
                if labels is not None:
                    for fold in range(5):
                        train=np.flatnonzero(labels!=fold);test=np.flatnonzero(labels==fold)
                        fit_status,_=design_status(g.iloc[train],columns,math)
                        if fit_status!='estimable':reason='training_'+fit_status;break
                        xx,xt,names=math.design(g.iloc[train],g.iloc[test],columns)
                        bt,obj,err=math.fit(xx,y[train]);assert obj<1e99;maximum_score=max(maximum_score,err)
                        predictions[family][repeat,test]=expit(xt@bt)
                pred=predictions[family][repeat]
                if reason!='estimated':pred[:]=np.nan
                rec=dict(info,repeat=repeat,model=family,status=reason,auc=np.nan,logloss=np.nan)
                if reason=='estimated':rec.update(score(y,pred))
                cv.append(rec)
                oof.extend(dict(info,repeat=repeat,model=family,human_gene_id=gid,Y=yy,prediction=pp) for gid,yy,pp in zip(g.human_gene_id,y,pred))
        if primary and len(g):
            reference=predictions['M1'];rng=np.random.default_rng(seed('cis_oof_bootstrap',*key))
            for family,pred in predictions.items():
                ready=np.isfinite(pred).all() and np.isfinite(reference).all()
                observed={m:np.mean([score(y,p)[m]-score(y,b)[m] for p,b in zip(pred,reference)]) if ready else np.nan for m in ['auc','logloss']}
                sampled=[]
                if ready:
                    # Resample genes jointly across repeats and paired models; folds are fixed.
                    for draw in range(2000):
                        ix=rng.integers(len(g),size=len(g));yy=y[ix]
                        sampled.append([np.mean([score(yy,p[ix])[m]-score(yy,b[ix])[m] for p,b in zip(pred,reference)]) for m in ['auc','logloss']])
                for i,metric in enumerate(['auc','logloss']):
                    values=np.array(sampled)[:,i] if sampled else np.array([]);finite=values[np.isfinite(values)]
                    bounds=np.quantile(finite,[.025,.975]) if len(finite) else [np.nan,np.nan]
                    intervals.append(dict(info,model=family,reference='M1',metric=metric,delta=observed[metric],lo=bounds[0],hi=bounds[1],
                        status='conditional_fixed_folds_gene_bootstrap' if ready else 'not_estimable',draws=2000 if ready else 0,finite_draws=len(finite)))
        ctx.table(results,'model_progress.tsv')
        print('CIS',json.dumps(row,default=str),flush=True)
    result=pd.DataFrame(results);assert len(result)==32
    result['q_boot_planned_32']=bh(result.p_boot);result['q_boot_primary_four']=np.nan
    result.loc[result.primary,'q_boot_primary_four']=bh(result.loc[result.primary,'p_boot'])
    for frame,name,columns in [(result,'recomputed_models.tsv',None),(coef,'recomputed_coefficients.tsv',KEYS+['term','beta','status']),
        (cv,'recomputed_cv.tsv',KEYS+['model','repeat','status','auc','logloss']),
        (foldrows,'paired_cv_folds.tsv.gz',KEYS+['repeat','human_gene_id','fold']),
        (oof,'out_of_fold_predictions.tsv.gz',KEYS+['repeat','model','human_gene_id','Y','prediction']),
        (nulls,'conditional_bootstrap_nulls.tsv.gz',KEYS+['draw','lrt_null']),
        (intervals,'paired_cv_uncertainty.tsv',KEYS+['model','metric','delta','lo','hi','status'])]:ctx.table(frame,name,columns=columns)
    (ctx.out/'model_numerical_qa.json').write_text(json.dumps(dict(qa,max_fitted_score_error=maximum_score,planned_tests=32,
        bootstrap='conditional parametric null restricted in the full model penalty space',
        cross_validation='10 deterministic class-stratified 5-fold partitions; identical folds across all models; training-only imputation and scaling',
        uncertainty='2000 gene bootstrap draws with paired predictions and folds held fixed; not a biological replicate interval'),indent=2)+'\n')
