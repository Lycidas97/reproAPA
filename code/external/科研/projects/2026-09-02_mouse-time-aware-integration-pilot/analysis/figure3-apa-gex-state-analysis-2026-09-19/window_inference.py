"""Uncertainty for the date-balanced effect actually shown in the proposed figure.
Biological wells nested in fixed dates provide variance; dates are not replicates.
Shared S0 variance enters source-target contrasts once, without dilution by dates.
"""
import argparse
import numpy as np,pandas as pd
from scipy.stats import t
import analyze as core
args=argparse.ArgumentParser();args.add_argument('--shared',action='store_true');shared=args.parse_args().shared
O=core.O;keys=(['axis_id'] if shared else [])+core.keys;setting='shared_p08' if shared else 'supported50_p08';prefix='shared_' if shared else ''
z=pd.read_csv(O/('shared_matched_library_effects.tsv.gz' if shared else 'matched_library_pair_effects.tsv.gz'),sep='\t');z=z[z.setting.eq(setting)&z.replicate_status_a.eq('biological_independent')&z.replicate_status_b.eq('biological_independent')&z.GEX_measurable];rows=[]
def component_test(means,variances,ns,weights):
 means=np.asarray(means);vv=np.asarray(variances)/np.asarray(ns)*np.square(weights);est=np.dot(means,weights);variance=vv.sum();denom=np.sum(vv**2/(np.asarray(ns)-1));df=variance**2/denom if denom>0 else np.nan
 return est,np.sqrt(variance),df
for kk,w in z[z.design.eq('same_library')].groupby(keys+['contrast','design']):
 u=w.groupby(['time_a','sample_unit_a'])[['APA_delta','GEX_delta']].mean().reset_index();assert u.groupby('sample_unit_a').time_a.nunique().max()==1
 stats=u.groupby('time_a')[['APA_delta','GEX_delta']].agg(['mean','var','count']);stats=stats[stats[('APA_delta','count')]>=2];nd=len(stats)
 if nd<2:continue
 row=dict(zip(keys+['contrast','design'],kk))|dict(n_dates=nd,n_bio_units=int(stats[('APA_delta','count')].sum()),source_units=0,dates=';'.join(stats.index),inference='independent biological wells nested in fixed dates; equal-date target mean',source_covariance='not_applicable')
 for metric in ['APA','GEX']:
  ss=stats[metric+'_delta'];est,se,df=component_test(ss['mean'],ss['var'],ss['count'],np.repeat(1/nd,nd));row[metric+'_estimate']=est;row[metric+'_se']=se;row[metric+'_df']=df;row[metric+'_date_min']=ss['mean'].min();row[metric+'_date_max']=ss['mean'].max();loo=(ss['mean'].sum()-ss['mean'])/(nd-1);row[metric+'_leave_date_min']=loo.min();row[metric+'_leave_date_max']=loo.max()
 rows.append(row)
validkeys=keys+['contrast','time_a'];valid={kk:set(zip(ww.sample_unit_a,ww.sample_unit_b)) for kk,ww in z[z.design.eq('true_source')].groupby(validkeys)}
if shared:s=pd.read_csv(O/'state_effect_master.tsv.gz',sep='\t');sources=[s]
else:sources=[pd.read_csv(f,sep='\t') for f in core.R.glob('gse*/apa_gex_states.tsv.gz') if '115943' in str(f)]
for s in sources:
 s=s[s.setting.eq(setting)&s.replicate_status.eq('biological_independent')]
 for kk,w in s.groupby(keys):
  ss=w[w.state.eq('S')&w.time.eq('D0')].groupby('sample_unit')[['APA','GEX']].mean();ns=len(ss)
  if ns<2:continue
  for state,ww in w[w.state.isin(['R','E','P','M'])].groupby('state'):
   ds=[]
   for tm,wwd in ww.groupby('time'):
    au=wwd.groupby('sample_unit')[['APA','GEX']].mean();pairs=valid.get(kk+(state+'_S',tm),set())
    if len(au)<2 or not all((ia,ib) in pairs for ia in au.index for ib in ss.index):continue
    ds.append((tm,au))
   nd=len(ds)
   if nd<2:continue
   row=dict(zip(keys,kk))|dict(contrast=state+'_S',design='true_source',n_dates=nd,n_bio_units=ns+sum(len(au) for _,au in ds),source_units=ns,dates=';'.join(tm for tm,_ in ds),inference='independent source wells and target wells nested in fixed dates',source_covariance='single S0 variance component retained across all target dates')
   for metric in ['APA','GEX']:
    means=[ss[metric].mean()]+[au[metric].mean() for _,au in ds];variances=[ss[metric].var()]+[au[metric].var() for _,au in ds];counts=[ns]+[len(au) for _,au in ds];weights=[-1]+[1/nd]*nd;est,se,df=component_test(means,variances,counts,weights);row[metric+'_estimate']=est;row[metric+'_se']=se;row[metric+'_df']=df;deltas=np.asarray(means[1:])-means[0];row[metric+'_date_min']=deltas.min();row[metric+'_date_max']=deltas.max();loo=(deltas.sum()-deltas)/(nd-1);row[metric+'_leave_date_min']=loo.min();row[metric+'_leave_date_max']=loo.max()
   rows.append(row)
r=pd.DataFrame(rows)
for metric in ['APA','GEX']:
 est=r[metric+'_estimate'];se=r[metric+'_se'];df=r[metric+'_df'];r[metric+'_CI95_low']=est-t.ppf(.975,df)*se;r[metric+'_CI95_high']=est+t.ppf(.975,df)*se;r[metric+'_CI90_low']=est-t.ppf(.95,df)*se;r[metric+'_CI90_high']=est+t.ppf(.95,df)*se;r[metric+'_p']=np.where(se>0,2*t.sf(abs(est/se),df),np.nan)
for margin in [.15,.25,.5]:r['GEX_TOST_p_'+str(margin)]=np.where(r.GEX_se>0,np.maximum(t.sf((r.GEX_estimate+margin)/r.GEX_se,r.GEX_df),t.cdf((r.GEX_estimate-margin)/r.GEX_se,r.GEX_df)),np.nan)
for col in ['APA_p','GEX_p']+['GEX_TOST_p_'+str(m) for m in [.15,.25,.5]]:r[col.replace('_p','_q')]=r.groupby(['species','contrast','design'])[col].transform(lambda x:pd.Series(core.bh(x.dropna()),index=x.dropna().index).reindex(x.index))
r['scope']='replicate-qualified fixed-date window; not all dates necessarily observable for every gene';core.save(r,prefix+'window_biological_inference.tsv.gz');print(r.groupby(['contrast','design']).agg(genes=('gene_id','size'),APA_q05=('APA_q',lambda x:(x<.05).sum()),GEX_equiv_q05=('GEX_TOST_q_0.25',lambda x:(x<.05).sum())).to_string())
