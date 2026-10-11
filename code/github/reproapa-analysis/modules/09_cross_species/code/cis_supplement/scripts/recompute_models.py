"""Corrected Firth fits, same-space nested tests and paired out-of-fold evaluation.
Read only frozen snapshots; all outputs are confined to this project.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize,brentq
from scipy.stats import chi2,false_discovery_control
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score,log_loss
P=Path(__file__).resolve().parents[1]
D=P/'sources/upstream/data/figure4-cis-architecture-v27/annot'
BASE=['abs_effect','disc_eligible','val_eligible','disc_pas']
REC=BASE+['disc_score']
CIS=['driver_conserved','frac_conserved_discovery_pas','A_g','phastcons_driver_mean','motif_rank_driver_mean']
FAMILIES={'M0':BASE,'M1':REC,'PAS':REC+CIS[:2],'Architecture':REC+['A_g'],'Sequence':REC+CIS[3:],'All cis':REC+CIS}
KEYS=['layer','contrast','discovery','tau']

def value_grad(X,y,b):
    eta=X@b;mu=expit(eta);w=mu*(1-mu)
    I=(X*w[:,None]).T@X
    sign,ld=np.linalg.slogdet(I)
    if sign<=0:return 1e100,np.zeros(len(b))
    inv=np.linalg.inv(I)
    h=w*np.einsum('ij,ij->i',X@inv,X)
    val=np.sum(np.logaddexp(0,eta)-y*eta)-.5*ld
    grad=-X.T@(y-mu+h*(.5-mu))
    return val,grad

def fit(X,y,fixed=None,start=None):
    fixed={} if fixed is None else fixed
    b=np.zeros(X.shape[1]) if start is None else np.array(start).copy()
    for j,v in fixed.items():b[j]=v
    free=np.array([j for j in range(len(b)) if j not in fixed])
    def fun(v):
        z=b.copy();z[free]=v
        f,g=value_grad(X,y,z);return f,g[free]
    opt=minimize(fun,b[free],method='BFGS',jac=True,options={'gtol':1e-7,'maxiter':250})
    b[free]=opt.x
    err=np.max(np.abs(fun(opt.x)[1]))
    if not np.isfinite(opt.fun) or err>1e-4:
        raise RuntimeError(f'Fit did not converge: score={err}; {opt.message}')
    return b,float(opt.fun),float(err)

def design(train,test,cols):
    a=train[cols].to_numpy(float);b=test[cols].to_numpy(float)
    med=np.nanmedian(a,axis=0)
    if not np.isfinite(med).all():raise ValueError('All-missing predictor')
    a=np.where(np.isfinite(a),a,med);b=np.where(np.isfinite(b),b,med)
    mean=a.mean(axis=0);sd=a.std(axis=0,ddof=1);keep=sd>1e-10
    # Drop constants consistently within a fit; inference records retained names.
    names=[c for c,k in zip(cols,keep) if k]
    return np.column_stack([np.ones(len(a)),((a-mean)/np.where(keep,sd,1))[:,keep]]),np.column_stack([np.ones(len(b)),((b-mean)/np.where(keep,sd,1))[:,keep]]),names

def profile_ci(X,y,b,opt,j):
    target=chi2.ppf(.95,1)
    def f(v):return 2*(fit(X,y,{j:v},b)[1]-opt)-target
    out=[]
    for s in [-1,1]:
        width=.2
        while f(b[j]+s*width)<0:
            width*=1.7
            if width>50:raise RuntimeError('CI failed to bracket')
        lo,hi=sorted([b[j],b[j]+s*width]);out.append(brentq(f,lo,hi,xtol=1e-5))
    return out

def verify():
    from scipy.optimize._numdiff import approx_derivative
    rng=np.random.default_rng(9);X=np.column_stack([np.ones(50),rng.normal(size=(50,2))]);y=rng.binomial(1,.3,50);b=np.array([-.5,.3,-.2])
    fd=approx_derivative(lambda z:value_grad(X,y,z)[0],b).ravel();analytic=value_grad(X,y,b)[1]
    assert np.max(abs(fd-analytic))<1e-6
    X0=np.ones((20,1));y0=np.r_[np.ones(3),np.zeros(17)]
    b0,_,_=fit(X0,y0);assert abs(expit(b0[0])-3.5/21)<1e-7
    sep=np.column_stack([np.ones(20),np.r_[np.zeros(10),np.ones(10)]])
    bs,_,_=fit(sep,np.r_[np.zeros(10),np.ones(10)]);assert np.isfinite(bs).all()
    bf,vf,_=fit(X,y)
    ind=minimize(lambda z:value_grad(X,y,z)[0],np.zeros(3),method='BFGS')
    assert np.max(abs(bf-ind.x))<1e-4
    return {'gradient_max_error':float(np.max(abs(fd-analytic))),'intercept_only_exact_solution':True,'complete_separation_finite':True,'independent_finite_difference_optimizer_agrees':True}

def run():
    started=time.time();qa=verify();mi=pd.read_csv(D/'model_input_with_seq_annot.tsv.gz',sep='\t')
    for c in CIS:
        if c+'_missing' in mi:mi.loc[mi[c+'_missing']==1,c]=np.nan
    results=[];coefs=[];cv=[];nulls=[];maxscore=0
    B=4999;REPEATS=10
    for config_idx,(key,g) in enumerate(mi.groupby(KEYS,sort=True)):
        info=dict(zip(KEYS,key));y=g.Y.to_numpy(float)
        X,_,names=design(g,g,REC+CIS)
        b,obj,score=fit(X,y);maxscore=max(maxscore,score)
        fixed={names.index(c)+1:0 for c in CIS if c in names}
        bn,on,score=fit(X,y,fixed,b);maxscore=max(maxscore,score)
        stat=max(0,2*(on-obj));rng=np.random.default_rng(922000+config_idx)
        null=[]
        for k in range(B):
            yy=rng.binomial(1,expit(X@bn))
            bb,oo,ss=fit(X,yy,start=bn);nb,no,ns=fit(X,yy,fixed,start=bn)
            null.append(max(0,2*(no-oo)));maxscore=max(maxscore,ss,ns)
        p=(1+sum(v>=stat for v in null))/(B+1)
        results.append({**info,'n':len(g),'positive':int(y.sum()),'lrt':stat,'df':len(fixed),'p_asymptotic':chi2.sf(stat,len(fixed)),'p_boot':p,'bootstrap_draws':B})
        nulls.extend({**info,'draw':i,'lrt_null':v} for i,v in enumerate(null))
        primary=key[0]=='supported50' and key[3]==.02
        for term in CIS:
            j=names.index(term)+1
            lo,hi=profile_ci(X,y,b,obj,j) if primary else [np.nan,np.nan]
            coefs.append({**info,'term':term,'beta':b[j],'ci_low':lo,'ci_high':hi})
        families=FAMILIES if primary else {n:FAMILIES[n] for n in ['M1','All cis']}
        for repeat in range(REPEATS):
            folds=list(StratifiedKFold(n_splits=5,shuffle=True,random_state=2200+repeat).split(g,y))
            for family,cols in families.items():
                pred=np.empty(len(g))
                for train,test in folds:
                    xx,xt,_=design(g.iloc[train],g.iloc[test],cols)
                    bt,_,ss=fit(xx,y[train]);maxscore=max(maxscore,ss);pred[test]=expit(xt@bt)
                cv.append({**info,'repeat':repeat,'model':family,'auc':roc_auc_score(y,pred),'logloss':log_loss(y,pred)})
        print(json.dumps({**info,'p_boot':p,'elapsed_s':round(time.time()-started)},ensure_ascii=False),flush=True)
    r=pd.DataFrame(results);r['q_boot_32']=false_discovery_control(r.p_boot)
    r.to_csv(P/'data/recomputed_models.tsv',sep='\t',index=False)
    pd.DataFrame(coefs).to_csv(P/'data/recomputed_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(cv).to_csv(P/'data/recomputed_cv.tsv',sep='\t',index=False)
    pd.DataFrame(nulls).to_csv(P/'data/bootstrap_nulls.tsv.gz',sep='\t',index=False)
    qa.update({'configs':len(r),'B':B,'cv_repeats':REPEATS,'cv_folds':5,'max_fitted_score_error':maxscore,'elapsed_seconds':time.time()-started,'bootstrap':'conditional parametric null under restricted M1 in full M2 penalty space','multiple_testing':'BH across all 32 combined-cis tests'})
    (P/'analysis/recompute_qa.json').write_text(json.dumps(qa,indent=2))
    print(r[(r.layer=='supported50')&(r.tau==.02)].to_string(index=False),flush=True)
if __name__=='__main__':run()
