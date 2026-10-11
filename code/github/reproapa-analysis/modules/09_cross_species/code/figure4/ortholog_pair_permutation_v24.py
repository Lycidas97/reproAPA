from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr,hypergeom
P=Path(__file__).resolve().parents[1];D=P/'data/figure4-ortholog-null-v24';A=P/'analysis/figure4-ortholog-null-v24';D.mkdir(exist_ok=True);A.mkdir(exist_ok=True)
N=P/'data/figure4-program-v9/all_mapped_nodes.tsv.gz';O=P/'data/figure4-program-v9/ortholog_universe.tsv';nd=pd.read_csv(N,sep='\t');ort=pd.read_csv(O,sep='\t');
codes=['W','C','H1','H2','H3','H4']; mice=['W','C']; hums=['H1','H2','H3','H4'];tau=.02;B=5000;rng=np.random.default_rng(240919)
n=nd[nd.layer.eq('supported50')&nd.measurement_pass&nd.code.isin(codes)].copy();
# Effects indexed by native gene, then vectorized recurrence.
def mat(contrast):
 a,b=contrast; z=n.pivot_table(index='gene_id',columns=['code','node'],values='position_index',aggfunc='first'); out=pd.DataFrame(index=z.index)
 for c in codes: out[c]=z[(c,b)]-z[(c,a)] if (c,b) in z and (c,a) in z else np.nan
 return out
sp={'PS':mat('SP'),'PM':mat('MP')}
pm=pd.read_csv(P/'data/figure4-program-v9/paired_pm_study_effects.tsv.gz',sep='\t'); pm=pm[pm.layer.eq('supported50')&pm.quality_pass]
pm=pm[(pm.n_sample_units>=2)&pm.code.isin(mice)|((pm.n_sample_units>=1)&pm.code.isin(hums))]
pmw=pm.pivot_table(index='gene_id',columns='code',values='mean_delta',aggfunc='first'); sp['PM']=pmw.reindex(columns=codes)
def hit(x,cs):
 q=x[cs]
 return (q.notna().sum(axis=1).ge(2)&q.le(-tau).sum(axis=1).ge(2)&q.ge(tau).sum(axis=1).eq(0)) if len(cs)>2 else q.le(-tau).all(axis=1)
rows=[]
for r in ort.itertuples(index=False): rows.append({'human_gene_id':r.human_gene_id,'mouse_gene_id':r.mouse_gene_id,'human_symbol':r.human_gene_symbol})
f=pd.DataFrame(rows)
for c in ['PS','PM']:
 f['h_'+c]=f.human_gene_id.map(hit(sp[c],hums)).fillna(False).astype(bool).to_numpy();f['m_'+c]=f.mouse_gene_id.map(hit(sp[c],mice)).fillna(False).astype(bool).to_numpy()
# features for matched null
for p,cs in [('h',hums),('m',mice)]:
 for metric,col in [('n',None),('pas','supported_pas_in_route')]:
  if metric=='n': vals=n[n.code.isin(cs)].groupby('gene_id').code.nunique()
  else: vals=n[n.code.isin(cs)].groupby('gene_id')[col].median()
  ids=f.human_gene_id if p=='h' else f.mouse_gene_id; f[p+'_'+metric]=ids.map(vals).fillna(-1).round().astype(int).to_numpy()
def perm(x,y,groups=None):
 x=np.asarray(x,bool);y=np.asarray(y,bool);groups=[np.arange(len(x))] if groups is None else groups; null=np.zeros(B,dtype=int)
 for g in groups:
  if len(g)==1: null += int(x[g][0] and y[g][0])
  else: null += hypergeom.rvs(len(g),int(y[g].sum()),int(x[g].sum()),size=B,random_state=rng)
 obs=int((x&y).sum());return dict(observed=obs,expected=null.mean(),lo=np.quantile(null,.025),hi=np.quantile(null,.975),p=(1+(null>=obs).sum())/(B+1),n=len(x),n_strata=len(groups))
out=[]
for c in ['PS','PM']:
 out.append(dict(contrast=c,null='unmatched',**perm(f['h_'+c],f['m_'+c])))
 z=f.copy();z['stratum']=z[['h_n','m_n','h_pas','m_pas']].astype(str).agg('|'.join,axis=1);groups=[g.to_numpy() for _,g in z.groupby('stratum').groups.items()];out.append(dict(contrast=c,null='matched_meas_pas',**perm(f['h_'+c],f['m_'+c],groups)))
pd.DataFrame(out).to_csv(D/'overlap_permutation_results.tsv',sep='\t',index=False)
# Continuous, all measured orthologs, no endpoint threshold selection.
cont=[]
for c in ['PS','PM']:
 x=f.human_gene_id.map(sp[c][hums].mean(axis=1));y=f.mouse_gene_id.map(sp[c][mice].mean(axis=1));q=pd.DataFrame({'x':x,'y':y}).dropna();x=q.x.to_numpy();y=q.y.to_numpy();rs=spearmanr(x,y).statistic;rp=pearsonr(x,y).statistic;ss=np.mean(x*y>0);null=np.empty((B,3))
 for i in range(B):
  yp=rng.permutation(y);null[i]=[spearmanr(x,yp).statistic,pearsonr(x,yp).statistic,np.mean(x*yp>0)]
 for j,m,o in [(0,'spearman',rs),(1,'pearson',rp),(2,'same_sign',ss)]:cont.append(dict(contrast=c,metric=m,n=len(x),observed=o,null_mean=null[:,j].mean(),null_lo=np.quantile(null[:,j],.025),null_hi=np.quantile(null[:,j],.975),p=(1+(null[:,j]>=o).sum())/(B+1)))
pd.DataFrame(cont).to_csv(D/'continuous_concordance_permutation.tsv',sep='\t',index=False);f.to_csv(D/'ortholog_pair_features.tsv',sep='\t',index=False)
(A/'manifest.json').write_text(json.dumps({'frozen':'supported50 tau=.02; mouse W/C proximal; human >=2 proximal and no distal','universe':len(f),'permutations':B,'matched_strata':'joint measured-study and median PAS count bins','input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [N,O]}},indent=2)+'\n')
print(pd.DataFrame(out).to_string(index=False));print(pd.DataFrame(cont).to_string(index=False))
