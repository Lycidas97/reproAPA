from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import t
P=Path(__file__).resolve().parent;M=P.parents[1];O=M/'outputs'/P.name;R=P/'raw';O.mkdir(exist_ok=True)
keys=['species','study_id','route','setting','gene_id','gene_symbol'];values=['APA','GEX','GEX_all']
def save(d,n):d.to_csv(O/n,sep='\t',index=False,float_format='%.10g')
def bh(p):
 a=np.asarray(p);o=np.argsort(a);v=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];r=np.empty(len(a));r[o]=np.minimum(v,1);return r

def pair(a,b,on,contrast,design):
 z=a.merge(b,on=on,suffixes=('_a','_b'))
 if z.empty:return z
 def qualify(r):
  ca=np.asarray(json.loads(r.site_counts_a));cb=np.asarray(json.loads(r.site_counts_b));na=np.asarray(json.loads(r.site_mc_a));nb=np.asarray(json.loads(r.site_mc_b));assert len(ca)==len(cb)
  return int(((ca+cb>=10)&(na+nb>=2)).sum())
 z['n_joint_PAS']=[qualify(r) for r in z.itertuples()];z=z[z.n_joint_PAS>=2].copy()
 for v in values:z[v+'_delta']=z[v+'_a']-z[v+'_b']
 z['contrast']=contrast;z['design']=design
 if 'GEX_count_a' in z:
  z['GEX_measurable']=(z.GEX_count_a+z.GEX_count_b>=20)&(z[['GEX_CPM_a','GEX_CPM_b']].max(axis=1)>=1)
 else:z['GEX_measurable']=False
 for c in ['library_key','time','sample_unit','replicate_status']:
  if c in on:z[c+'_a']=z[c];z[c+'_b']=z[c]
 retain=keys+['contrast','design','n_joint_PAS','GEX_measurable']+[v+'_delta' for v in values]+[c+s for c in ['library_key','time','sample_unit','replicate_status','n_mc','n_cells','pas_counts'] for s in ['_a','_b']]
 return z[retain]
def collapse_window(w):
 out=[]
 for kk,gg in w.groupby(keys):
  rr=dict(zip(keys,kk))
  for vv in values:rr[vv]=gg.groupby(['time','sample_unit'])[vv].mean().groupby('time').mean().mean()
  for vv in ['site_counts','site_mc']:rr[vv]=json.dumps(np.asarray([json.loads(x) for x in gg[vv]]).sum(axis=0).tolist())
  for vv in ['n_mc','n_cells','pas_counts','GEX_count']:
   if vv in gg:rr[vv]=gg[vv].sum()
  if 'GEX_CPM' in gg:rr['GEX_CPM']=gg.GEX_CPM.mean()
  for vv in ['library_key','time']:rr[vv]=';'.join(sorted(gg[vv].unique()))
  rr['sample_unit']='intermediate_window';rr['replicate_status']='aggregated_window_not_replicate';out.append(rr)
 return pd.DataFrame(out)
if __name__ == '__main__':
 allpairs=[];designs=[];fullbg=[]
 for folder in sorted(R.glob('gse*')):
  d=pd.read_csv(folder/'apa_gex_states.tsv.gz',sep='\t');d=d[~d.route.str.contains('reference')].copy();dd=pd.read_csv(folder/'sampling_design.tsv',sep='\t');dd['source_type']=np.where(dd.state.eq('S')&dd.time.isin(['D0','source','Stage_I_day_0']),'true_source',np.where(dd.state.eq('S'),'induced_or_mixed_S','not_source'));dd['target_type']=np.where(dd.state.eq('P'),np.where(dd.route.str.contains('reference'),'reference_PSC','experimental_P'),'not_target');designs.append(dd)
  for a,b in [('P','M'),('P','S'),('R','S'),('E','S'),('M','S'),('E','R'),('P','E'),('P','R'),('E','M'),('R','M')]:
   z=pair(d[d.state.eq(a)],d[d.state.eq(b)],keys+['library_key','time','sample_unit','replicate_status'],a+'_'+b,'same_library')
   if len(z):allpairs.append(z)
  for route,w in d.groupby('route'):
   source=w[w.state.eq('S')&w.time.isin(['D0','source','Stage_I_day_0'])]
   typ='true_source'
   if source.empty and route=='GSE221739':source=w[w.state.eq('S')&w.time.eq('D3')];typ='induced_D3_reference'
   if source.empty:continue
   for a in ['R','E','P','M']:
    z=pair(w[w.state.eq(a)],source,keys,a+'_S',typ)
    if len(z):allpairs.append(z)
   for b in ['E','R']:
    inter=w[w.state.eq(b)]
    if inter.empty or w[w.state.eq('P')].empty:continue
    if set(inter.library_key)&set(w[w.state.eq('P')].library_key):continue
    z=pair(w[w.state.eq('P')],collapse_window(inter),keys,'P_'+b,'target_vs_intermediate_window')
    if len(z):allpairs.append(z)
  print('PAIRS',folder.name,sum(map(len,allpairs)),flush=True)
  # Full GEX universe independent of APA eligibility, for endpoint comparisons.
  g=pd.read_csv(folder/'full_gex_states.tsv.gz',sep='\t');g=g[~g.route.str.contains('reference')]
  for route,w in g.groupby('route'):
   for con in ['P_S','P_M']:
    aa=w[w.state.eq('P')];bb=w[w.state.eq(con[-1])]
    if con=='P_S':bb=bb[bb.time.isin(['D0','source','Stage_I_day_0'])]
    if con=='P_M':
     common=set(aa.library_key)&set(bb.library_key);aa=aa[aa.library_key.isin(common)];bb=bb[bb.library_key.isin(common)]
    if aa.empty or bb.empty:continue
    x=aa.groupby('gene_id').agg(count_a=('count','sum'),cpm_a=('CPM','max'));y=bb.groupby('gene_id').agg(count_b=('count','sum'),cpm_b=('CPM','max'));q=x.join(y,how='outer').fillna(0);q=q[(q.count_a+q.count_b>=20)&(q[['cpm_a','cpm_b']].max(axis=1)>=1)].reset_index();q['route']=route;q['contrast']=con;q['species']=w.species.iloc[0];q['study_id']=w.study_id.iloc[0];fullbg.append(q)
 sampling=pd.concat(designs,ignore_index=True);provenance=pd.read_csv(M/'analysis/figure3-state-statistics-preflight-2026-09-13/library_sampling_units.tsv',sep='\t')[['library_key','unit_basis','source_url']];save(sampling.merge(provenance,on='library_key',how='left',validate='many_to_one'),'sampling_design.tsv');save(pd.concat(fullbg,ignore_index=True),'GEX_endpoint_backgrounds.tsv.gz')
 z=pd.concat(allpairs,ignore_index=True);save(z,'matched_library_pair_effects.tsv.gz');del allpairs
 # Equal weight source units then source dates, then target units and target dates. Keep designs separate.
 k=keys+['contrast','design'];v=[x+'_delta' for x in values]
 u=z.groupby(k+['time_a','sample_unit_a','time_b','sample_unit_b'],dropna=False)[v].mean().reset_index();u=u.groupby(k+['time_a','sample_unit_a','time_b'])[v].mean().groupby(k+['time_a','sample_unit_a'])[v].mean().reset_index();date=u.groupby(k+['time_a'])[v].mean().reset_index();save(date,'matched_date_effects.tsv.gz');save(u,'matched_target_unit_effects.tsv.gz')
 a=date.groupby(k).agg(APA_delta=('APA_delta','mean'),GEX_delta=('GEX_delta','mean'),GEX_all_delta=('GEX_all_delta','mean'),n_dates=('time_a','size'),date_APA_min=('APA_delta','min'),date_APA_max=('APA_delta','max'),date_GEX_min=('GEX_delta','min'),date_GEX_max=('GEX_delta','max')).reset_index();a['cell_selection_effect']=a.GEX_delta-a.GEX_all_delta;save(a,'matched_study_effects.tsv.gz')
 # Native-axis candidate screen. At least two independent studies, preserve route contradictions.
 e=a[((a.contrast=='P_S')&(a.design=='true_source'))|((a.contrast=='P_M')&(a.design=='same_library'))].copy();entries=[]
 for kk,w in e.groupby(['species','setting','gene_id','gene_symbol','contrast']):
  by=w.groupby('study_id').agg(delta=('APA_delta','mean'),lo=('APA_delta','min'),hi=('APA_delta','max'))
  pos=by.index[(by.lo>=.02)].tolist();neg=by.index[(by.hi<=-.02)].tolist();conf=by.index[(by.lo<=-.02)&(by.hi>=.02)].tolist();direction=1 if len(pos)>=2 else (-1 if len(neg)>=2 else 0)
  entries.append(dict(zip(['species','setting','gene_id','gene_symbol','contrast'],kk))|dict(n_studies_evaluable=len(by),positive_studies=';'.join(pos),negative_studies=';'.join(neg),within_study_conflicts=';'.join(conf),direction=direction,candidate=bool(direction),n_opposite_studies=len(neg) if direction==1 else len(pos) if direction==-1 else 0))
 c=pd.DataFrame(entries);save(c,'native_candidate_support.tsv');primary=c[c.setting.eq('supported50_p08')&c.candidate];union=primary.groupby(['species','gene_id','gene_symbol']).agg(entrances=('contrast',lambda q:';'.join(sorted(q)))).reset_index();save(union,'candidate_union.tsv');save(e.merge(union[['species','gene_id']],on=['species','gene_id']),'candidate_endpoint_effects.tsv.gz')
 # Native state relations; shared-axis version will supersede these.
 rels=[]
 for r in primary.itertuples():
  w=a[(a.species==r.species)&(a.gene_id==r.gene_id)&a.setting.eq('supported50_p08')]
  for rr in ['R_S','E_S','P_E','P_R','E_M','R_M','M_S']:
   q=w[w.contrast.eq(rr)]
   if rr.endswith('_S'):q=q[q.design.eq('true_source')]
   else:q=q[q.design.eq('same_library')]
   by=q.groupby('study_id').APA_delta.mean();rels.append(dict(species=r.species,gene_id=r.gene_id,gene_symbol=r.gene_symbol,entrance=r.contrast,relation=rr,reference_direction=r.direction,evaluable_studies=len(by),same_direction_studies=int((by*r.direction>=.02).sum()),opposite_studies=int((by*r.direction<=-.02).sum()),source_interpretation_allowed=r.contrast=='P_S'))
 save(pd.DataFrame(rels),'native_state_relation_evidence.tsv')
 # Formal fixed-date biological-replicate inference. No timepoints counted as replicates.
 tests=[];main=z[z.setting.eq('supported50_p08')&z.replicate_status_a.eq('biological_independent')&z.replicate_status_b.eq('biological_independent')&z.GEX_measurable]
 # Paired state differences within the same biological library.
 for kk,w in main[main.design.eq('same_library')].groupby(k+['time_a']):
  ww=w.groupby('sample_unit_a')[v].mean();n=len(ww)
  if n<2:continue
  row=dict(zip(k+['time_a'],kk))|dict(n_a=n,n_b=n,df=n-1,inference='paired_biological_units')
  for metric in ['APA','GEX']:
   val=ww[metric+'_delta'];est=val.mean();se=val.std(ddof=1)/np.sqrt(n);row[metric+'_estimate']=est;row[metric+'_se']=se
  tests.append(row)
 # Source-target comparisons within fixed target dates; source and target wells independent.
 # Recover source/target levels from state records rather than treating Cartesian pair differences as replicates.
 valid_source_pair_units={kk:set(zip(ww.sample_unit_a,ww.sample_unit_b)) for kk,ww in main[main.design.eq('true_source')].groupby(['route','gene_id','contrast','time_a'])}
 for folder in sorted(R.glob('gse*')):
  d=pd.read_csv(folder/'apa_gex_states.tsv.gz',sep='\t');d=d[d.setting.eq('supported50_p08')&d.replicate_status.eq('biological_independent')]
  for kk,w in d.groupby(keys):
   ss=w[w.state.eq('S')&w.time.eq('D0')].groupby('sample_unit')[values].mean()
   if len(ss)<2:continue
   for (st,tm),dd in w[w.state.isin(['R','E','P','M'])].groupby(['state','time']):
    aa=dd.groupby('sample_unit')[values].mean()
    if len(aa)<2:continue
    # Require gene/date pair eligibility already established, including joint multi-PAS.
    valid_pairs=valid_source_pair_units.get((kk[2],kk[4],st+'_S',tm),set())
    if not all((ia,ib) in valid_pairs for ia in aa.index for ib in ss.index):continue
    row=dict(zip(keys,kk))|dict(contrast=st+'_S',design='true_source',time_a=tm,n_a=len(aa),n_b=len(ss),inference='Welch_independent_biological_units')
    for metric in ['APA','GEX']:
     va=aa[metric].var()/len(aa);vb=ss[metric].var()/len(ss);se=np.sqrt(va+vb);df=(va+vb)**2/(va**2/(len(aa)-1)+vb**2/(len(ss)-1)) if se>0 else np.nan;row[metric+'_estimate']=aa[metric].mean()-ss[metric].mean();row[metric+'_se']=se;row[metric+'_df']=df
    tests.append(row)
 test=pd.DataFrame(tests)
 for metric in ['APA','GEX']:
  df=test.get(metric+'_df',pd.Series(np.nan,index=test.index)).fillna(test['df']);se=test[metric+'_se'];est=test[metric+'_estimate'];ok=(se>0)&df.notna();test[metric+'_df']=df;test[metric+'_CI90_low']=est-t.ppf(.95,df)*se;test[metric+'_CI90_high']=est+t.ppf(.95,df)*se;test[metric+'_p']=np.where(ok,2*t.sf(abs(est/se),df),np.nan)
 for margin in [.15,.25,.5]:
  se=test.GEX_se;df=test.GEX_df;est=test.GEX_estimate;test['GEX_TOST_p_'+str(margin)]=np.where(se>0,np.maximum(t.sf((est+margin)/se,df),t.cdf((est-margin)/se,df)),np.nan)
 for col in ['APA_p','GEX_p']+['GEX_TOST_p_'+str(m) for m in [.15,.25,.5]]:
  test[col.replace('_p','_q')]=test.groupby(['species','contrast','design'])[col].transform(lambda x:pd.Series(bh(x.dropna()),index=x.dropna().index).reindex(x.index))
 save(test,'fixed_date_biological_inference.tsv.gz');save(test.groupby(['species','route','contrast','design']).agg(gene_date_tests=('gene_id','size'),genes=('gene_id','nunique'),n_units_min=('n_a','min')).reset_index(),'equivalence_feasibility.tsv')
 summary=dict(native_candidate_counts=union.groupby('species').size().to_dict(),native_entrances=primary.groupby(['species','contrast']).size().to_dict(),library_pairs=len(z),study_effects=len(a),formal_tests=len(test),GEX_equiv_q05=int((test['GEX_TOST_q_0.25']<.05).sum()))
 (O/'summary.json').write_text(json.dumps({k:str(v) if isinstance(v,dict) and any(isinstance(x,tuple) for x in v) else v for k,v in summary.items()},indent=2));print(summary)
