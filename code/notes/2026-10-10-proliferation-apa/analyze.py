"""Fixed library/state association; MCs are measurement units, not replicates."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
KEY=['species','study_id','route']; SCORES=['Cycle_S','Cycle_G2M','Cycle','Adjusted_S','Adjusted_G2M','Adjusted_Cycle']
COV=['panel_GEX','panel_PAS','GEX_per_cell','PAS_per_cell','MC_cells']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(f,name):
 pd.DataFrame(f).to_csv(HERE/name,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
def seed(name):return 20261010+int(hashlib.sha256(str(name).encode()).hexdigest()[:12],16)
def residual(x,c):
 c=np.column_stack([np.ones(len(x)),c]);sd=np.std(c[:,1:],axis=0);c[:,1:]=np.divide(c[:,1:]-c[:,1:].mean(0),sd,out=np.zeros_like(c[:,1:]),where=sd>1e-12)
 rank=np.linalg.matrix_rank(c);return x-c@np.linalg.lstsq(c,x,rcond=None)[0],rank
def association(x,y,c):
 x=np.asarray(x,float);y=np.asarray(y,float);rx,rank=residual(x,c);ry,_=residual(y,c)
 n=len(x);rawsd=x.std();sx=rx.std();sy=ry.std();rec=dict(n=n,covariate_rank=rank,residual_df=n-rank-1,status='available',correlation_unadjusted=np.nan,partial_correlation=np.nan,beta_per_score_unit=np.nan,beta_per_score_SD=np.nan,score_residual_fraction=(sx/rawsd)**2 if rawsd>1e-12 else np.nan,APA_residual_fraction=(sy/y.std())**2 if y.std()>1e-12 else np.nan)
 if rawsd>1e-12 and y.std()>1e-12:rec['correlation_unadjusted']=np.corrcoef(x,y)[0,1]
 if n-rank-1<5:rec['status']='insufficient_residual_df';return rec
 if sx<=1e-12 or sy<=1e-12:rec['status']='constant_residual_score_or_APA';return rec
 beta=rx@ry/(rx@rx);rec.update(partial_correlation=(rx@ry)/np.sqrt((rx@rx)*(ry@ry)),beta_per_score_unit=beta,beta_per_score_SD=beta*rawsd)
 return rec

def library_states(meta,design):
 rows=[]
 for key,q in meta.groupby(KEY+['library_key','node'],sort=True):
  assert q.experimental_time.nunique(dropna=False)==1, key
  rec=dict(zip(KEY+['library_key','state'],key));rec.update(n_metacells=len(q),n_cells=int(q.n_cells.sum()),qualified=len(q)>=3 and q.n_cells.sum()>=30,time=str(q.experimental_time.iloc[0]))
  for score in SCORES:rec[score]=q[score].mean()
  d=design[design.route.eq(rec['route'])&design.library_key.eq(rec['library_key'])&design.state.eq(rec['state'])]
  rec['documented_source']=rec['state']=='S' and d.time.astype(str).isin(['D0','source','Stage_I_day_0']).any()
  assert d.sample_unit.nunique()<=1 and d.replicate_status.nunique()<=1
  rec['sample_unit']=';'.join(d.sample_unit.astype(str).unique());rec['replicate_status']=';'.join(d.replicate_status.astype(str).unique())
  rows.append(rec)
 return pd.DataFrame(rows)
def score_contrasts(libs):
 out=[];pair_rows=[]
 for key,q in libs.groupby(KEY):
  valid=q[q.qualified]
  for a,b in [('P','M'),('P','S'),('M','S')]:
   paired=valid[valid.state.eq(a)].merge(valid[valid.state.eq(b)],on='library_key',suffixes=('_target','_reference'),validate='one_to_one')
   for score in SCORES:
    rec=dict(zip(KEY,key),contrast=a+'_'+b,score=score,n_pairs=len(paired),delta=np.nan,time_equal_delta=np.nan,status='absent_paired_states')
    if len(paired):
     delta=paired[score+'_target']-paired[score+'_reference'];rec.update(delta=delta.mean(),time_equal_delta=delta.groupby(paired.time_target).mean().mean(),status='descriptive_paired',n_positive=int((delta>0).sum()),n_negative=int((delta<0).sum()))
     pair_rows.extend(dict(zip(KEY,key),contrast=a+'_'+b,score=score,library_key=z.library_key,time=z.time_target,sample_unit=z.sample_unit_target,replicate_status=z.replicate_status_target,delta=float(v)) for z,v in zip(paired.itertuples(),delta))
     if len(paired)>=2:
      rng=np.random.default_rng(seed((key,a,b,score)));boot=delta.to_numpy()[rng.integers(len(delta),size=(10000,len(delta)))].mean(1)
      rec.update(library_boot_p025=np.quantile(boot,.025),library_boot_p975=np.quantile(boot,.975))
    out.append(rec)
  # Unpaired documented source and observed S kept separate.
  for a in ['P','M']:
   target=valid[valid.state.eq(a)]
   for source_name,reference in [('documented_source',valid[valid.documented_source]),('observed_S',valid[valid.state.eq('S')])]:
    for score in SCORES:
     rec=dict(zip(KEY,key),contrast=a+'_S_'+source_name,score=score,n_pairs=0,n_target_libraries=len(target),n_reference_libraries=len(reference),delta=np.nan,time_equal_delta=np.nan,status='absent_state')
     if len(target) and len(reference):rec.update(delta=target[score].mean()-reference[score].mean(),time_equal_delta=target.groupby('time')[score].mean().mean()-reference.groupby('time')[score].mean().mean(),status='descriptive_unpaired')
     out.append(rec)
 return pd.DataFrame(out),pd.DataFrame(pair_rows)

def aggregate_contexts(f):
 librows=[];out=[]
 for key,q in f.groupby(KEY+['score','APA_variant']):
  for scope in ['all_SMP','P','M','S']:
   z=q[q.status.eq('available')]
   if scope!='all_SMP':z=z[z.state.eq(scope)]
   rec=dict(zip(KEY+['score','APA_variant'],key),scope=scope,n_contexts=len(z),n_libraries=z.library_key.nunique(),status='unavailable')
   lib=z.groupby('library_key')[['partial_correlation','beta_per_score_SD','beta_per_score_unit']].mean()
   if len(lib):
    librows.extend(dict(rec,status='descriptive',library_key=k,n_states=int((z.library_key==k).sum()),**v.to_dict()) for k,v in lib.iterrows())
    rec.update(status='descriptive',partial_correlation=lib.partial_correlation.mean(),beta_per_score_SD=lib.beta_per_score_SD.mean(),n_negative_libraries=int((lib.partial_correlation<0).sum()),n_positive_libraries=int((lib.partial_correlation>0).sum()),min_library_correlation=lib.partial_correlation.min(),max_library_correlation=lib.partial_correlation.max())
    time=z.groupby(['time','library_key']).partial_correlation.mean().groupby('time').mean();rec['time_equal_correlation']=time.mean()
    if len(lib)>=2:
     rng=np.random.default_rng(seed((key,scope)));v=lib.partial_correlation.to_numpy();boot=v[rng.integers(len(v),size=(10000,len(v)))].mean(1)
     rec.update(library_boot_p025=np.quantile(boot,.025),library_boot_p975=np.quantile(boot,.975),leave_one_library_min=min((v.sum()-v)/(len(v)-1)),leave_one_library_max=max((v.sum()-v)/(len(v)-1)))
   out.append(rec)
 return pd.DataFrame(out),pd.DataFrame(librows)

def main():
 x=np.load(HERE/'inputs.npz',allow_pickle=False);assert str(x['protocol_sha256'])==sha(HERE/'PROTOCOL.zh.md')
 panels=pd.DataFrame(json.loads(str(x['panels_json'])));mapping=pd.DataFrame(json.loads(str(x['mapping_json'])));controls=pd.DataFrame(json.loads(str(x['controls_json'])))
 metas=[];design=[]
 for k in x.files:
  if k.endswith('_meta_json'):
   sid=k[:-10];m=pd.DataFrame(json.loads(str(x[k])));m.study_id=sid;metas.append(m)
   design.append(pd.DataFrame(json.loads(str(x[sid+'_design_json']))))
 meta=pd.concat(metas,ignore_index=True);design=pd.concat(design,ignore_index=True)
 libs=library_states(meta,design);contrast,pairs=score_contrasts(libs)
 state_rows=[]
 for key,q in libs[libs.qualified].groupby(KEY+['state']):
  r=dict(zip(KEY+['state'],key),n_libraries=len(q),n_times=q.time.nunique(),n_metacells=int(q.n_metacells.sum()),n_cells=int(q.n_cells.sum()))
  for score in SCORES:r[score]=q[score].mean();r['time_equal_'+score]=q.groupby('time')[score].mean().mean()
  state_rows.append(r)
 rows=[];mcrows=[];roster=[]
 for z in panels[panels.status.eq('available')].itertuples():
  cid=z.context_id;m=pd.DataFrame(json.loads(str(x[z.study_id+'_meta_json']))).iloc[x[cid+'_rows']].copy().reset_index(drop=True)
  m['study_id']=z.study_id
  assert m.experimental_time.nunique(dropna=False)==1,cid
  dw=x[cid+'_DWUI'];depth=x[cid+'_depth'];gex=x[cid+'_gex_counts'];total=m.total_gex.to_numpy();cells=m.n_cells.to_numpy()
  assert np.isfinite(dw).all() and np.all(depth>=20) and dw.shape==depth.shape==gex.shape
  panelgex=np.log1p(gex/total[:,None]*1e4).mean(1);panelpas=np.log1p(depth).mean(1)
  c=np.column_stack([panelgex,panelpas,np.log1p(total/cells),np.log1p(m.total_pas.to_numpy()/cells),np.log(cells)])
  outcomes={'mean_DWUI':dw.mean(1),'median_DWUI':np.median(dw,axis=1),'thin20_mean_DWUI':x[cid+'_thin20_mean'].mean(0)}
  for score in SCORES:
   for variant,y in outcomes.items():
    rec=dict(species=z.species,study_id=z.study_id,route=z.route,library_key=z.library_key,state=z.state,time=str(m.experimental_time.iloc[0]),context_id=cid,n_genes=z.n_genes,score=score,APA_variant=variant,**association(m[score].to_numpy(),y,c));rows.append(rec)
  rr=m.copy();rr['context_id']=cid;rr['panel_genes']=z.n_genes
  for j,name in enumerate(COV):rr[name]=c[:,j]
  for name,y in outcomes.items():rr[name]=y
  mcrows.append(rr)
  roster.extend(dict(context_id=cid,gene_id=g) for g in x[cid+'_gene_ids'])
 f=pd.DataFrame(rows)
 unit=libs[KEY+['library_key','state','sample_unit','replicate_status']]
 f=f.merge(unit,on=KEY+['library_key','state'],validate='many_to_one')
 summary,libassoc=aggregate_contexts(f)
 library_unit=unit.drop(columns='state').drop_duplicates()
 assert not library_unit.duplicated(KEY+['library_key']).any()
 libassoc=libassoc.merge(library_unit,on=KEY+['library_key'],validate='many_to_one')
 save(meta,'metacell_cycle_scores.tsv.gz');save(mapping,'marker_mapping.tsv');save(controls,'expression_matched_controls.tsv');save(panels,'panel_eligibility.tsv');save(roster,'fixed_panel_genes.tsv.gz')
 save(libs,'library_state_scores.tsv');save(state_rows,'route_state_scores.tsv');save(contrast,'state_score_contrasts.tsv');save(pairs,'paired_library_score_differences.tsv');save(f,'context_associations.tsv');save(summary,'route_associations.tsv');save(libassoc,'library_associations.tsv');save(pd.concat(mcrows,ignore_index=True),'metacell_apa_cycle.tsv.gz')
 (HERE/'remote_input_sha256.json').write_text(str(x['provenance_json'])+'\n')
 print(panels.groupby(['study_id','status']).size().to_string()); print(contrast[contrast.contrast.eq('P_M')&contrast.score.isin(['Cycle','Adjusted_Cycle'])][KEY+['score','n_pairs','delta','time_equal_delta']].to_string(index=False))
 print(summary[summary.score.eq('Cycle')&summary.APA_variant.eq('mean_DWUI')][KEY+['scope','n_contexts','n_libraries','partial_correlation','n_negative_libraries','library_boot_p025','library_boot_p975']].to_string(index=False))
if __name__=='__main__':main()
