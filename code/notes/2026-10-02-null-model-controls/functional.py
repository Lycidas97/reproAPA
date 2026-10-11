"""Exact nested functional nulls using the already verified frozen strata."""
import sys, importlib.util
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.special import logsumexp
from common import *

PREV=ROOT/'notes/2026-10-02-null-result-diagnostics'
spec=importlib.util.spec_from_file_location('previous_exact',PREV/'functional_exact.py')
F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)

@lru_cache(maxsize=200000)
def component(N,K,n,odds=1):
    lo=max(0,n+K-N);hi=min(n,K);x=np.arange(lo,hi+1)
    logs=hypergeom.logpmf(x,N,K,n)+x*np.log(odds)
    values=np.exp(logs-logsumexp(logs))
    return lo,values

def distribution(params,odds=1):
    offset=0;values=np.array([1.])
    for N,K,n in params:
        if not n or not K:continue
        lo,p=component(int(N),int(K),int(n),odds)
        offset+=lo;values=np.convolve(values,p)
    assert abs(values.sum()-1)<1e-8
    return offset,values/values.sum()

def tail(dist,hits):
    lo,p=dist
    return min(1.,float(p[max(0,hits-lo):].sum()))

def rejection_probability(test,generated,alpha):
    tlo,tp=test;glo,gp=generated
    tails=np.r_[np.cumsum(tp[::-1])[::-1],0.]
    positions=np.clip(np.arange(glo,glo+len(gp))-tlo,0,len(tp))
    return float(gp[tails[positions]<=alpha].sum())

def setup(bg,comparison):
    split=bg.stratum.str.split('|',expand=True)
    if comparison=='G2_vs_G1':
        assert split.shape[1]==5
        fields={'L0_eligibility':[0],'L1_measurement':[0,2,3],
            'L2_baseline':[0,2,3,4],'L2_PAS_first':[0,2,3,1],'L3_full':[0,1,2,3,4]}
    else:
        assert split.shape[1]==2
        # G0 has no APA-architecture covariates. Keep its expression-matched
        # test unchanged in later stages so every stage retains all 28 families.
        fields={'L0_eligibility':[0],'L1_measurement':[0,1],
            'L2_baseline':[0,1],'L2_PAS_first':[0,1],'L3_full':[0,1]}
    result={}
    for model,cols in fields.items():
        code,_=pd.factorize(split[cols].agg('|'.join,axis=1),sort=True)
        N=np.bincount(code);n=np.bincount(code,weights=bg.target.astype(int)).astype(int)
        result[model]=(code,N,n)
    return result

def main():
    frozen=read(PREV/'functional_frozen_strata.tsv.gz')
    terms,aliases,ancestors=F.B.ontology()
    annotations={}
    for db in ['GO','KEGG']:
        annotations[db]={sp:(F.B.load_go_genes(sp,terms,aliases,ancestors) if db=='GO' else F.K.load_kegg(sp)[0]) for sp in ['mouse','human']}
    for p in [F.G/'scripts/run_go.py',F.G/'scripts/run_go_nested.py',F.G/'scripts/run_kegg_nested.py']+list((F.G/'sources').glob('*'))+list((F.G/'sources_kegg').glob('*')):
        if p.is_file():INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    power=[];mobility=[];replay=[]
    for db in ['GO','KEGG']:
        old=read(PREV/f'{db.lower()}_exact.tsv.gz');out=[]
        for aid,previous in old.groupby('analysis_id',sort=False):
            species,contrast,comparison,target=aid.split('|')
            bg=frozen[frozen.database.eq(db)&frozen.analysis_id.eq(aid)].reset_index(drop=True)
            assert len(bg)==previous.N.iloc[0] and bg.target.sum()==previous.n.iloc[0]
            models=setup(bg,comparison);full=list(models)[-1];index={g:i for i,g in enumerate(bg.gene_id)}
            family=F.N.term_family(set(bg.gene_id),annotations[db][species]);assert set(family)==set(previous.term_id)
            ordered=previous.sort_values(['K','term_id']).reset_index(drop=True)
            selected=set(ordered.iloc[np.round(np.linspace(0,len(ordered)-1,5)).astype(int)].term_id)
            parts=[];max_error=0.
            for model,(codes,N,n) in models.items():
                mobility.append(dict(database=db,analysis_id=aid,model=model,n_genes=len(bg),target=int(bg.target.sum()),strata=len(N),
                    frozen_gene_fraction=float(N[(n==0)|(n==N)].sum()/N.sum()),
                    expected_changed_fraction=float((2*n*(1-n/N)).sum()/N.sum())))
            for row in previous.itertuples():
                ii=np.array([index[g] for g in family[row.term_id]],int)
                assert len(ii)==row.K and int(bg.target.iloc[ii].sum())==row.hits
                tests={};parameters={}
                for model,(codes,N,n) in models.items():
                    K=np.bincount(codes[ii],minlength=len(N));params=list(zip(N,K,n));parameters[model]=params
                    dist=distribution(params);tests[model]=dist
                    p=tail(dist,row.hits);mu=float((n*K/N).sum())
                    sd=np.sqrt(sum(nn*(kk/nnn)*(1-kk/nnn)*(nnn-nn)/(nnn-1) for nnn,kk,nn in params if nnn>1))
                    parts.append(dict(database=db,analysis_id=aid,term_id=row.term_id,term_name=row.term_name,model=model,
                        N=row.N,n=row.n,K=row.K,hits=row.hits,expected=mu,null_sd=sd,excess=row.hits-mu,p=p))
                    if model==full:max_error=max(max_error,abs(p-row.p_exact))
                if row.term_id in selected:
                    for odds in [1,2,4,8]:
                        gen=distribution(parameters[full],odds)
                        values=np.arange(gen[0],gen[0]+len(gen[1]));mean=float(values@gen[1])
                        for model,test in tests.items():
                            for label,alpha in [('nominal',.05),('bonferroni_terms',.05/len(previous))]:
                                power.append(dict(database=db,analysis_id=aid,term_id=row.term_id,term_name=row.term_name,K=row.K,n=row.n,N=row.N,
                                    generator=full,tested_model=model,odds_ratio=odds,threshold=label,alpha=alpha,
                                    generated_mean=mean,rejection_probability=rejection_probability(test,gen,alpha)))
            assert max_error<1e-9,(db,aid,max_error)
            replay.append(dict(database=db,analysis_id=aid,terms=len(previous),max_abs_full_p_error=max_error))
            out.extend(parts)
            print('FUNCTIONAL',db,aid,'terms',len(previous),'models',len(models),'replay',max_error,flush=True)
        result=pd.DataFrame(out)
        result['q_comparison']=result.groupby(['analysis_id','model']).p.transform(bh)
        # Comparable full models have different names for G0 and G1 backgrounds.
        result['q_database_model']=result.groupby('model').p.transform(bh)
        result['q_all_new_positions']=bh(result.p)
        save(result,f'{db.lower()}_ladder.tsv.gz')
    save(pd.DataFrame(power),'functional_power.tsv.gz');save(pd.DataFrame(mobility),'functional_mobility.tsv')
    save(pd.DataFrame(replay),'functional_replay.tsv')
    manifest('functional_manifest.json',status='PASS',power='exact conditional Fisher odds ratios; term-size quantile panel',pandas=pd.__version__)

if __name__=='__main__':main()
