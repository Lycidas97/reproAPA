"""Predefined complete-aUTR motif tests; fixed labels, fixed subpopulation, exact strata."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from prepare_events import HERE,ROOT,save,sha
MOTIFS={'ESCC_7mer_m8':'AGCACTT','ESCC_7mer_A1':'GCACTTA','ESCC_8mer':'AGCACTTA','U4':'TTTT'}
def occurrence_offsets(sequence,word):return [i for i in range(len(sequence)-len(word)+1) if sequence[i:i+len(word)]==word]
def bh(values):
 a=np.asarray(values,float);order=np.argsort(a);q=np.minimum.accumulate((a[order]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];out=np.zeros(len(a));out[order]=np.minimum(q,1);return out

def main():
 protocol=json.loads((HERE/'protocol.json').read_text());assert sha(HERE/'PROTOCOL.zh.md')==protocol['protocol_sha256']
 r=pd.read_csv(ROOT/'notes/2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz',sep='\t');r=r[r.contrast.eq('P_M')][['human_gene_id','mouse_gene_id','human_gene_symbol','mouse_gene_symbol','foreground_002']].copy()
 inv=pd.read_csv(HERE/'gene_interval_inventory.tsv.gz',sep='\t');seq=pd.read_csv(HERE/'full_aUTR_sequences.tsv',sep='\t');assert seq.complete_ACGT.all();i=inv.set_index(['species','gene_id']);features=[];hits=[]
 for row in seq.to_dict('records'):
  s=row['sequence'];feat={k:row[k] for k in ['species','gene_id','chrom','strand','start0','end0']};feat.update(aUTR_length=len(s),GC=sum(s.count(x) for x in 'GC')/len(s))
  for key,word in MOTIFS.items():
   offsets=occurrence_offsets(s,word);feat[key+'_count']=len(offsets)
   for off in offsets:
    lo=int(row['start0'])+off if row['strand']=='+' else int(row['end0'])-off-len(word)
    hits.append(dict(species=row['species'],gene_id=row['gene_id'],motif=key,word=word,aUTR_offset0=off,chrom=row['chrom'],strand=row['strand'],genomic_start0=lo,genomic_end0=lo+len(word)))
  feat['ESCC_ge7']=feat['ESCC_7mer_m8_count']+feat['ESCC_7mer_A1_count']>0;feat['ESCC_8mer']=feat['ESCC_8mer_count']>0
  feat['U4_density']=feat['U4_count']/max(1,len(s)-3);features.append(feat)
 f=pd.DataFrame(features);save(f,'species_motif_features.tsv.gz');save(pd.DataFrame(hits),'motif_hits.tsv.gz')
 for sp in ('human','mouse'):
  g=inv[inv.species.eq(sp)].drop(columns=['species','foreground_002','gene_symbol']);g=g.rename(columns={c:sp+'_'+c for c in g.columns if c!='gene_id'}).rename(columns={'gene_id':sp+'_gene_id'});r=r.merge(g,on=sp+'_gene_id',how='left',validate='one_to_one')
 eligible=r.human_status.eq('eligible')&r.mouse_status.eq('eligible');r['paired_eligible']=eligible;save(r,'paired_coverage.tsv.gz');b=r[eligible].copy().reset_index(drop=True)
 for sp in ('human','mouse'):
  ff=f[f.species.eq(sp)][['gene_id','ESCC_ge7','ESCC_8mer','U4_density']].rename(columns={c:sp+'_'+c for c in ['ESCC_ge7','ESCC_8mer','U4_density']}).rename(columns={'gene_id':sp+'_gene_id'});b=b.merge(ff,on=sp+'_gene_id',validate='one_to_one')
 b['both_ESCC_ge7']=b.human_ESCC_ge7&b.mouse_ESCC_ge7;b['both_ESCC_8mer']=b.human_ESCC_8mer&b.mouse_ESCC_8mer;b['both_U4_min_density']=np.minimum(b.human_U4_density,b.mouse_U4_density)
 endpoints=['human_ESCC_ge7','mouse_ESCC_ge7','both_ESCC_ge7','human_ESCC_8mer','mouse_ESCC_8mer','both_ESCC_8mer','human_U4_density','mouse_U4_density','both_U4_min_density']
 b['loglength_mean']=(np.log(b.human_aUTR_length)+np.log(b.mouse_aUTR_length))/2;b['GC_mean']=(b.human_GC+b.mouse_GC)/2;b['baseline_mean']=(b.human_baseline+b.mouse_baseline)/2;b['PAS_number_mean']=(b.human_PAS_number+b.mouse_PAS_number)/2
 cuts={c:float(b[c].median()) for c in ['loglength_mean','GC_mean','baseline_mean','PAS_number_mean']}
 for c,cut in cuts.items():b[c+'_bin']=b[c].gt(cut).astype(int)
 b['paired_measurement_pattern']=b.human_study_pattern+'|'+b.mouse_study_pattern
 Y=b[endpoints].astype(float).to_numpy();fg=b.foreground_002.to_numpy(dtype=bool);obs=Y[fg].mean(axis=0);counts=Y[fg].sum(axis=0);nfg=int(sum(fg));rng=np.random.default_rng(protocol['seed']);N=protocol['n_permutations'];results=[];strata=[];balances=[];draws={};modelinfo=[]
 covars=['loglength_mean','GC_mean','baseline_mean','PAS_number_mean']
 for sp in ('human','mouse'):
  for c in ('aUTR_length','GC','baseline','PAS_number','mean_GEX','APA_coverage'):
   name=sp+'_'+c
   if c in ('aUTR_length','mean_GEX','APA_coverage'):
    name+='__log1p';b[name]=np.log1p(b[sp+'_'+c])
   covars.append(name)
 for model in ('coarse','full'):
  columns=['paired_measurement_pattern','loglength_mean_bin','GC_mean_bin']+(['baseline_mean_bin','PAS_number_mean_bin'] if model=='full' else [])
  b['stratum_'+model]=b[columns].astype(str).agg('|'.join,axis=1)
  null=np.zeros((N,len(endpoints)));weights=np.zeros(len(b));exchangeable=0;nconditional=0
  for s,g in b.groupby('stratum_'+model,sort=True):
   inds=g.index.to_numpy();k=int(sum(fg[inds]));m=len(inds);weights[inds]=k/m
   strata.append(dict(model=model,stratum=s,n_background=m,n_foreground=k,exchangeable=(0<k<m),human_gene_ids=';'.join(g.human_gene_id)))
   if not k:continue
   if k==m:
    null+=Y[inds].sum(axis=0);continue
   exchangeable+=k;nconditional+=1
   for start in range(0,N,10000):
    size=min(10000,N-start);u=rng.random((size,m));choose=np.argpartition(u,k-1,axis=1)[:,:k];null[start:start+size]+=Y[inds][choose].sum(axis=1)
  null/=nfg;draws[model]=null;expected_feature_weight_sum=weights.sum();assert abs(expected_feature_weight_sum-nfg)<1e-10
  for c in covars:
   values=b[c].to_numpy(float);sd=np.nanstd(values,ddof=0);actual=np.nanmean(values[fg]);valid=np.isfinite(values);expected=float(sum(weights[valid]*values[valid])/sum(weights[valid]));smd=float((actual-expected)/sd) if sd else 0.
   balances.append(dict(model=model,covariate=c,foreground_mean=actual,null_expected_mean=expected,background_SD=sd,SMD=smd))
  maxsmd=max(abs(z['SMD']) for z in balances if z['model']==model);qualified=nfg>=5 and len(b)>=20 and exchangeable>=5;info=dict(model=model,n_background=len(b),n_foreground=nfg,exchangeable_foreground=exchangeable,exchangeable_strata=nconditional,max_abs_SMD=maxsmd,formal_claim_eligible=qualified,balanced=maxsmd<=.20);modelinfo.append(info)
  for j,name in enumerate(endpoints):
   expected=float(null[:,j].mean());p=(1+np.count_nonzero(null[:,j]>=obs[j]-1e-15))/(N+1)
   results.append(info|dict(endpoint=name,statistic_unit='foreground_fraction' if 'ESCC' in name else 'mean_overlapping_4mer_density',observed=float(obs[j]),observed_foreground_sum=float(counts[j]),expected=expected,expected_foreground_sum=expected*nfg,observed_minus_expected=float(obs[j]-expected),enrichment_ratio=float(obs[j]/expected) if expected else np.nan,p_upper=float(p),null_sd=float(null[:,j].std())))
 res=pd.DataFrame(results);assert len(res)==protocol['planned_family_size'];res['q18']=bh(res.p_upper);res['positive_evidence_gate']=res.formal_claim_eligible&res.balanced&res.q18.le(.05)
 save(res,'enrichment_tests.tsv');save(pd.DataFrame(strata),'stratum_inventory.tsv');save(pd.DataFrame(balances),'covariate_balance.tsv');save(b,'paired_motif_background.tsv.gz');np.savez_compressed(HERE/'null_draws.npz',**draws)
 for c in ['loglength_mean','GC_mean','baseline_mean','PAS_number_mean']:assert c in cuts
 (HERE/'analysis_audit.json').write_text(json.dumps(dict(protocol_sha256=protocol['protocol_sha256'],paired_bg=int(len(b)),paired_fg=nfg,binary_cuts=cuts,models=modelinfo,n_permutations=N,seed=protocol['seed'],planned_family_size=len(res),negative_noneligible_not_imputed=True),indent=2)+'\n')
 print('paired',len(b),'foreground',nfg);print(pd.DataFrame(modelinfo).to_string(index=False));print(res[['model','endpoint','observed_foreground_sum','expected','p_upper','q18','positive_evidence_gate']].to_string(index=False));print(b[b.foreground_002][['human_gene_symbol','human_ESCC_ge7','mouse_ESCC_ge7','human_ESCC_8mer','mouse_ESCC_8mer']].to_string(index=False))
if __name__=='__main__':main()
