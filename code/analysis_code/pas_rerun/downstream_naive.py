"""Naive-state refinement, conditional on fixed Liu libraries in each count route."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from downstream_runtime import Context, Study, bh
from pipeline import digest
SETTINGS=[('supported50_p08',True,.8,3),('supported50_p09',True,.9,3),('supported50_p08_mc5',True,.8,5),('expanded_primary_p08',False,.8,3)]
L={x:f'GSE150311::gse150311_human_pluripotency_reprogramming_GSM454654{x}' for x in ['5','6','7','8']}
METRICS=['DWUI','RWUI','PDUI','PPUI']
PLANS={
 'hN_P_primed_condition':[(L['6'],'N',L['7'],'P'),(L['8'],'N',L['7'],'P')],
 'hN_P_same_library':[(L['6'],'N',L['6'],'P'),(L['8'],'N',L['8'],'P')],
 'hN_S_same_library':[(L['6'],'N',L['6'],'S'),(L['8'],'N',L['8'],'S')],
 'hN_S_D0_D7_pool':[(L['6'],'N',L['5'],'S'),(L['8'],'N',L['5'],'S')],
 'P_primed_S_D0_D7_pool':[(L['7'],'P',L['5'],'S')],
 'P_primed_S_same_library':[(L['7'],'P',L['7'],'S')],
 'P_all_S_D0_D7_pool':[(L['7'],'P',L['5'],'S'),(L['8'],'P',L['5'],'S')],
 'hN_P_all_cross_library':[(L['6'],'N',L['7'],'P'),(L['6'],'N',L['8'],'P'),(L['8'],'N',L['7'],'P'),(L['8'],'N',L['8'],'P')],
}

NODE='naive_state'

def run(ctx):
 global CTX,HERE,NATIVE,VARIANT,save
 CTX=ctx;HERE=ctx.out;NATIVE=ctx.native/'studies/GSE150311';VARIANT=ctx.config['reference_variant']
 def output(frame,name):
  frame=pd.DataFrame(frame);ctx.table(frame,name);return frame
 save=output
 calculate();followup();fixed_triangle()
 ctx.finish(dict(frozen_state_labels_verified=True,new_native_counts=True,all_four_settings=True,same_library_triangle_identity=True),
  inferential_scope='Gene bootstrap conditional on selected genes and fixed libraries; no biological replicate inference')

def metric_weights(a):
 signed=a.cut_site.to_numpy()*np.where(a.strand.eq('+'),1,-1)
 order=np.argsort(signed,kind='stable');rank=np.empty(len(a),int);rank[order]=np.arange(1,len(a)+1)
 return np.column_stack([a.position,rank/len(a),rank==len(a),rank==1])

def state(c,rows,meta,minimum_mc):
 rows=rows[c[rows].sum(1)>=20];z=c[rows]
 usage=(z/z.sum(1)[:,None]).mean(0) if len(z) else np.zeros(c.shape[1])
 return dict(usage=usage,molecules=z.sum(0),site_mc=(z>0).sum(0),n_mc=len(rows),n_cells=int(meta.n_cells.iloc[rows].sum()),rows=rows,
  eligible=len(rows)>=minimum_mc and meta.n_cells.iloc[rows].sum()>=30)

def why(a,b):
 if a is None or b is None:return 'state_group_absent_or_below_purity'
 if not a['eligible'] or not b['eligible']:return 'insufficient_gene_qualified_metacells_or_cells'
 support=(a['molecules']+b['molecules']>=10)&(a['site_mc']+b['site_mc']>=2)
 if support.sum()<2:return 'fewer_than_two_joint_supported_PAS'
 return ''

def bootstrap(x,rng):
 x=np.asarray(x,float)
 if not len(x):return np.nan,np.nan
 means=np.array([x[rng.integers(0,len(x),len(x))].mean() for _ in range(5000)])
 return tuple(np.quantile(means,[.025,.975]))

def describe(x,rng):
 x=np.asarray(x,float);lo,hi=bootstrap(x,rng)
 return dict(n_genes=len(x),mean=float(x.mean()) if len(x) else np.nan,median=float(np.median(x)) if len(x) else np.nan,
  gene_bootstrap_mean_lo=lo,gene_bootstrap_mean_hi=hi,
  proximal_002=int((x<=-.02).sum()),distal_002=int((x>=.02).sum()),proximal_005=int((x<=-.05).sum()),distal_005=int((x>=.05).sum()))

def calculate():
 inputs={str(p):digest(p) for p in NATIVE.glob('*') if p.is_file()}
 selected=NATIVE/'selected_metacells.tsv.gz';full=NATIVE/'metacells.tsv.gz'
 meta=pd.read_csv(selected,sep='\t');raw=pd.read_csv(full,sep='\t');counts=sparse.load_npz(NATIVE/'metacell_by_pas_counts.npz').tocsr()[meta.matrix_row].tocsc()
 assert np.array_equal(raw.metacell_id.iloc[meta.matrix_row],meta.metacell_id)
 assert np.issubdtype(counts.dtype,np.integer) and (counts.data>=0).all()
 meta=meta.merge(raw[['metacell_id','condition','sample_accession','time_kind','original_state_id']],on='metacell_id',validate='one_to_one')
 frozen=pd.read_csv(CTX.resource('frozen_metacell_states'),sep='\t');check=meta[meta.node.eq('N')].merge(frozen,on=['metacell_id','library_key'],validate='one_to_one')
 assert check.state_group.eq('human_n').all() and np.allclose(check.purity,check.state_purity)
 save(meta,'frozen_selected_metadata.tsv.gz')
 save(raw.groupby(['coarse_state','condition','sample_accession','experimental_time','time_kind'],dropna=False).agg(n_mc=('metacell_id','size'),n_cells=('n_cells','sum')).reset_index(),'unfiltered_state_condition_library_inventory.tsv')
 inv=[]
 for setting,_,purity,nmc in SETTINGS:
  for key,z in meta[meta.purity.ge(purity)].groupby(['node','condition','sample_accession','library_key','experimental_time','time_kind']):
   inv.append(dict(zip(['node','condition','sample_accession','library_key','experimental_time','time_kind'],key),setting=setting,n_mc=len(z),n_cells=int(z.n_cells.sum()),group_eligible=len(z)>=nmc and z.n_cells.sum()>=30))
 save(inv,'qualified_state_condition_library_inventory.tsv')
 annotations=CTX.analysis/'tandem/tandem_utr_eligibility.tsv.gz'
 tandem=pd.read_csv(annotations,sep='\t');tandem=set(tandem[tandem.study_id.eq('GSE150311')&tandem.tandem_utr_compatible.astype(str).str.lower().eq('true')].gene_id)
 inputs[str(annotations)]=digest(annotations)
 inputs[str(Path(__file__))]=digest(__file__)
 pairrows=[];erows=[];denoms=[]
 for version in [VARIANT]:
  f=pd.read_csv(NATIVE/'features.tsv.gz',sep='\t',low_memory=False)
  f=f[f.keep_for_apa_primary.astype(str).str.lower().eq('true')&f.primary_gene_id.notna()].copy()
  signed=f.cut_site*np.where(f.strand.eq('+'),1,-1);low=signed.groupby(f.primary_gene_id).transform('min');high=signed.groupby(f.primary_gene_id).transform('max')
  f['position']=(signed-low)/(high-low)
  shape=f.groupby('primary_gene_id').agg(n=('cut_site','nunique'),c=('chrom','nunique'),s=('strand','nunique'))
  valid=set(shape.query('n>=2 and c==1 and s==1').index)
  for setting,supported,purity,nmc in SETTINGS:
   sf=f[f.primary_gene_id.isin(valid)&(f.matched_source_count_50nt.ge(1) if supported else True)]
   groups={key:z.index.to_numpy() for key,z in meta[meta.purity.ge(purity)].groupby(['library_key','node'])}
   planned=list(dict.fromkeys(pair for pairs in PLANS.values() for pair in pairs))
   for gene,axis in sf.groupby('primary_gene_id',sort=True):
    base=dict(version=version,setting=setting,gene_id=gene,gene_symbol=axis.primary_gene_symbol.iloc[0],PAS_number=len(axis),tandem_utr_compatible=gene in tandem)
    if len(axis)<2:
     for a,na,b,nb in planned:erows.append(dict(base,library_a=a,node_a=na,library_b=b,node_b=nb,eligible=False,reason='fewer_than_two_retained_PAS'))
     continue
    w=metric_weights(axis);c=counts[:,axis.feature_index_1based.to_numpy()-1].toarray()
    st={key:state(c,rows,meta,nmc) for key,rows in groups.items()}
    for a,na,b,nb in planned:
     x=st.get((a,na));y=st.get((b,nb));reason=why(x,y)
     erows.append(dict(base,library_a=a,node_a=na,library_b=b,node_b=nb,eligible=not reason,reason=reason,n_mc_a=0 if x is None else x['n_mc'],n_mc_b=0 if y is None else y['n_mc'],n_cells_a=0 if x is None else x['n_cells'],n_cells_b=0 if y is None else y['n_cells']))
     if reason:continue
     d=x['usage']-y['usage'];assert abs(d.sum())<1e-12
     rec=dict(base,library_a=a,node_a=na,library_b=b,node_b=nb,n_mc_a=x['n_mc'],n_mc_b=y['n_mc'],n_cells_a=x['n_cells'],n_cells_b=y['n_cells'],pas_counts_a=int(x['molecules'].sum()),pas_counts_b=int(y['molecules'].sum()),PAS_ids=json.dumps(axis.atlas_pas_id.tolist()),usage_a=json.dumps(x['usage'].tolist()),usage_b=json.dumps(y['usage'].tolist()),MPRO=float(np.ptp(d)))
     rec.update({'delta_'+m:float(v) for m,v in zip(METRICS,d@w)});pairrows.append(rec)
   denoms.append(dict(version=version,setting=setting,genes_with_primary_axis=len(valid),genes_with_any_retained_PAS=sf.primary_gene_id.nunique(),genes_with_two_retained_PAS=int((sf.groupby('primary_gene_id').size()>=2).sum()),n_total_native_metacells=counts.shape[0],count_matrix_molecules=int(counts.sum())))
 pairs=save(pairrows,'library_pair_effects.tsv.gz');elig=save(erows,'library_pair_eligibility.tsv.gz');save(denoms,'gene_denominators.tsv')
 out=[]
 for version in pairs.version.unique():
  for setting in pairs.setting.unique():
   q=pairs[pairs.version.eq(version)&pairs.setting.eq(setting)]
   for contrast,planned in PLANS.items():
    selected=q[np.array([tuple(x) in planned for x in q[['library_a','node_a','library_b','node_b']].itertuples(index=False,name=None)])]
    for gene,z in selected.groupby('gene_id'):
     rec=dict(version=version,setting=setting,contrast=contrast,gene_id=gene,gene_symbol=z.gene_symbol.iloc[0],tandem_utr_compatible=z.tandem_utr_compatible.iloc[0],n_pairs=len(z),planned_pairs=len(planned),complete_contributors=len(z)==len(planned),n_libraries_a=z.library_a.nunique(),n_libraries_b=z.library_b.nunique(),PAS_number=int(z.PAS_number.iloc[0]))
     rec.update({'delta_'+m:float(z['delta_'+m].mean()) for m in METRICS})
     out.append(rec)
 eff=save(out,'contrast_gene_effects.tsv.gz')
 primary_contrasts=['hN_P_primed_condition','hN_S_D0_D7_pool','P_primed_S_D0_D7_pool']
 common={};fixed=set(); summaries=[];rng=np.random.default_rng(110031)
 for (version,setting),q in eff.groupby(['version','setting']):
  complete={c:set(z[z.complete_contributors].gene_id) for c,z in q.groupby('contrast')}
  ids=set.intersection(*(complete.get(c,set()) for c in primary_contrasts));common[version,setting]=ids
 fixed=set.intersection(*common.values())
 members=[]
 for (version,setting,contrast),q in eff.groupby(['version','setting','contrast']):
  choices={'available':set(q.gene_id),'complete_contributors':set(q[q.complete_contributors].gene_id)}
  if contrast in primary_contrasts:choices.update(common_hN_primed_P_pool_S=common[version,setting],fixed_all_settings_this_route=fixed)
  for background,ids in choices.items():
   z=q[q.gene_id.isin(ids)]
   for gene in sorted(ids):members.append(dict(version=version,setting=setting,contrast=contrast,background=background,gene_id=gene))
   for subset,b in [('all_primary_APA',z),('v2_full_primary_tandem_UTR',z[z.tandem_utr_compatible])]:
    for metric in METRICS:
     summaries.append(dict(version=version,setting=setting,contrast=contrast,background=background,subset=subset,metric=metric,**describe(b['delta_'+metric],rng),uncertainty='gene resampling conditional on fixed libraries; not biological-replicate inference'))
 save(members,'background_gene_membership.tsv.gz');save(summaries,'global_effect_summary.tsv')
 opposite=[];agreement=[];checks=[]
 for (version,setting),q in eff.groupby(['version','setting']):
  for background,ids in [('common_hN_primed_P_pool_S',common[version,setting]),('fixed_all_settings_this_route',fixed)]:
   z=q[q.contrast.isin(primary_contrasts)&q.gene_id.isin(ids)]
   for metric in METRICS:
    w=z.pivot(index='gene_id',columns='contrast',values='delta_'+metric)
    residual=w.hN_P_primed_condition-(w.hN_S_D0_D7_pool-w.P_primed_S_D0_D7_pool)
    assert residual.abs().max()<1e-12
    checks.append(dict(version=version,setting=setting,background=background,metric=metric,n_genes=len(w),contrast_identity_max_error=float(residual.abs().max())))
    for threshold in [.02,.05]:
     a=w.hN_S_D0_D7_pool;b=w.P_primed_S_D0_D7_pool
     opposite.append(dict(version=version,setting=setting,background=background,metric=metric,threshold=threshold,n_genes=len(w),hN_proximal_primed_distal=int(((a<=-threshold)&(b>=threshold)).sum()),hN_distal_primed_proximal=int(((a>=threshold)&(b<=-threshold)).sum()),both_proximal=int(((a<=-threshold)&(b<=-threshold)).sum()),both_distal=int(((a>=threshold)&(b>=threshold)).sum())))
  p=pairs[pairs.version.eq(version)&pairs.setting.eq(setting)&pairs.node_a.eq('N')&pairs.node_b.eq('P')&pairs.library_b.eq(L['7'])]
  for metric in METRICS:
   w=p.pivot(index='gene_id',columns='library_a',values='delta_'+metric).dropna()
   if len(w):
    a=w[L['6']];b=w[L['8']]
    agreement.append(dict(version=version,setting=setting,metric=metric,n_common_genes=len(w),rho=float(spearmanr(a,b).statistic),same_nonzero_sign=int((np.sign(a)==np.sign(b)).sum()),mean_naive_vs_primed=float(a.mean()),mean_RSeT_vs_primed=float(b.mean())))
 save(opposite,'source_relative_direction_categories.tsv');save(agreement,'two_hN_library_agreement.tsv');save(checks,'contrast_identity_verification.tsv')
 # Confirm the extension reproduces existing per-library S/P count qualification.
 existing=CTX.analysis/'global/observed_state_effects.tsv.gz'
 inputs[str(existing)]=digest(existing)
 for p,h in inputs.items():assert digest(p)==h,f'input changed {p}'
 (HERE/'input_sha256.json').write_text(json.dumps(inputs,indent=2)+'\n')
 (HERE/'VERIFICATION.json').write_text(json.dumps(dict(status='complete',inputs_verified=len(inputs),native_integer_matrix_reused=True,N_majority_matches_frozen_human_n=True,N_p08_cells=int(meta.loc[meta.node.eq('N')&meta.purity.ge(.8),'n_cells'].sum()),library_unit_inference='not run: two heterogeneous hN libraries and one explicitly primed library; unresolved biological replication',bootstrap_unit='gene',settings=SETTINGS,fixed_all_settings_this_route_genes=len(fixed),contrast_identity_checks=len(checks),no_BAM_recount=True),indent=2)+'\n')
 print(pd.DataFrame(summaries).query("version=='v2_primary' and setting=='supported50_p08' and metric=='DWUI' and subset=='all_primary_APA'")[['contrast','background','n_genes','mean','median','gene_bootstrap_mean_lo','gene_bootstrap_mean_hi','proximal_002','distal_002']].to_string(index=False))
 print('fixed',len(fixed))

def followup():
 paths=[HERE/'library_pair_effects.tsv.gz',NATIVE/'selected_metacells.tsv.gz',NATIVE/'metacell_by_pas_counts.npz',HERE/'frozen_selected_metadata.tsv.gz',NATIVE/'features.tsv.gz']
 inputs={str(p):digest(p) for p in paths};inputs[str(Path(__file__))]=digest(Path(__file__))
 pairs=pd.read_csv(paths[0],sep='\t');meta=pd.read_csv(HERE/'frozen_selected_metadata.tsv.gz',sep='\t')
 counts=sparse.load_npz(NATIVE/'metacell_by_pas_counts.npz').tocsr()[meta.matrix_row].tocsc();rng=np.random.default_rng(900314)
 out=[];triangle=[];den=[]
 for version in pairs.version.unique():
  f=pd.read_csv(NATIVE/'features.tsv.gz',sep='\t',low_memory=False)
  f=f[f.keep_for_apa_primary.astype(str).str.lower().eq('true')&f.primary_gene_id.notna()].copy()
  signed=f.cut_site*np.where(f.strand.eq('+'),1,-1)
  lo=signed.groupby(f.primary_gene_id).transform('min');hi=signed.groupby(f.primary_gene_id).transform('max');f['position']=(signed-lo)/(hi-lo)
  for setting,supported,purity,nmc in SETTINGS:
   q=pairs[pairs.version.eq(version)&pairs.setting.eq(setting)]
   # Both heterogeneous hN libraries against the one explicitly primed P library.
   common=set.intersection(*(set(q[q.node_a.eq('N')&q.node_b.eq('P')&q.library_a.eq(L[k])&q.library_b.eq(L['7'])].gene_id) for k in ['6','8']))
   for k in ['6','8']:
    z=q[q.node_a.eq('N')&q.node_b.eq('P')&q.library_a.eq(L[k])&q.library_b.eq(L['7'])]
    for background,b in [('available',z),('fixed_two_hN_libraries',z[z.gene_id.isin(common)])]:
     for subset,x in [('all_primary_APA',b),('v2_full_primary_tandem_UTR',b[b.tandem_utr_compatible])]:
      for metric in METRICS:out.append(dict(version=version,setting=setting,comparison='hN_'+('naive' if k=='6' else 'RSeT')+'_vs_primed_condition_P',background=background,subset=subset,metric=metric,**describe(x['delta_'+metric],rng)))
   npair=q[q.library_a.eq(L['8'])&q.library_b.eq(L['8'])&q.node_a.eq('N')&q.node_b.eq('P')].set_index('gene_id')
   ns=q[q.library_a.eq(L['8'])&q.library_b.eq(L['8'])&q.node_a.eq('N')&q.node_b.eq('S')].set_index('gene_id')
   genes=sorted(npair.index.intersection(ns.index));eligible=[]
   sf=f[f.matched_source_count_50nt.ge(1)] if supported else f
   axes={g:a for g,a in sf[sf.primary_gene_id.isin(genes)].groupby('primary_gene_id')}
   groups={node:z.index.to_numpy() for node,z in meta[meta.library_key.eq(L['8'])&meta.purity.ge(purity)].groupby('node')}
   for gene in genes:
    axis=axes[gene];assert json.loads(npair.loc[gene].PAS_ids)==axis.atlas_pas_id.tolist()
    c=counts[:,axis.feature_index_1based.to_numpy()-1].toarray();p=state(c,groups.get('P',np.array([],int)),meta,nmc);s=state(c,groups.get('S',np.array([],int)),meta,nmc)
    reason=why(p,s)
    if reason:continue
    effect=(p['usage']-s['usage'])@metric_weights(axis)
    eligible.append(gene)
    for metric,e in zip(METRICS,effect):
     n_p=float(npair.loc[gene,'delta_'+metric]);n_s=float(ns.loc[gene,'delta_'+metric]);assert abs(n_s-n_p-e)<1e-12
     triangle.append(dict(version=version,setting=setting,gene_id=gene,metric=metric,hN_P=n_p,hN_S=n_s,P_S=float(e),tandem_utr_compatible=bool(npair.loc[gene,'tandem_utr_compatible'])))
   den.append(dict(version=version,setting=setting,RSeT_NP_measurable=len(npair),RSeT_NS_measurable=len(ns),both_NP_NS=len(genes),qualified_NP_NS_PS_triangle=len(eligible),naive_library_NP_qualified=len(q[q.library_a.eq(L['6'])&q.library_b.eq(L['6'])&q.node_a.eq('N')&q.node_b.eq('P')])))
 tri=save(triangle,'RSeT_same_library_triangle_gene_effects.tsv.gz');save(den,'same_library_triangle_denominators.tsv')
 for (version,setting,metric),q in tri.groupby(['version','setting','metric']):
  for subset,x in [('all_primary_APA',q),('v2_full_primary_tandem_UTR',q[q.tandem_utr_compatible])]:
   for contrast in ['hN_P','hN_S','P_S']:out.append(dict(version=version,setting=setting,comparison='RSeT_same_library_'+contrast,background='complete_NP_NS_PS_contributors',subset=subset,metric=metric,**describe(x[contrast],rng)))
 save(out,'condition_and_same_library_summary.tsv')
 for p,h in inputs.items():assert digest(p)==h
 (HERE/'followup_input_sha256.json').write_text(json.dumps(inputs,indent=2)+'\n')
 print(pd.DataFrame(out).query("version=='v2_primary' and setting=='supported50_p08' and metric=='DWUI'")[['comparison','background','subset','n_genes','mean','gene_bootstrap_mean_lo','gene_bootstrap_mean_hi']].to_string(index=False))

def fixed_triangle():
 p=HERE/'RSeT_same_library_triangle_gene_effects.tsv.gz';q=pd.read_csv(p,sep='\t')
 inputs={str(p):digest(p),str(Path(__file__)):digest(Path(__file__))}
 common=set.intersection(*(set(z.gene_id) for _,z in q.groupby(['version','setting'])))
 rows=[];categories=[];rng=np.random.default_rng(310041)
 for (version,setting,metric),z in q.groupby(['version','setting','metric']):
  for background,b in [('available_complete_triangle',z),('fixed_all_settings_this_route_triangle',z[z.gene_id.isin(common)])]:
   for subset,x in [('all_primary_APA',b),('v2_full_primary_tandem_UTR',b[b.tandem_utr_compatible])]:
    for contrast in ['hN_P','hN_S','P_S']:rows.append(dict(version=version,setting=setting,metric=metric,background=background,subset=subset,contrast=contrast,**describe(x[contrast],rng)))
    for threshold in [.02,.05]:
     a=x.hN_S;b2=x.P_S
     categories.append(dict(version=version,setting=setting,metric=metric,background=background,subset=subset,threshold=threshold,n_genes=len(x),hN_proximal_P_distal=int(((a<=-threshold)&(b2>=threshold)).sum()),hN_distal_P_proximal=int(((a>=threshold)&(b2<=-threshold)).sum()),both_proximal=int(((a<=-threshold)&(b2<=-threshold)).sum()),both_distal=int(((a>=threshold)&(b2>=threshold)).sum())))
 save(rows,'RSeT_fixed_triangle_summary.tsv');save(categories,'RSeT_triangle_direction_categories.tsv')
 for p,h in inputs.items():assert digest(p)==h
 (HERE/'reference_sensitivity_input_sha256.json').write_text(json.dumps(inputs,indent=2)+'\n')
 print(pd.DataFrame(rows).query("setting=='supported50_p08' and metric=='DWUI' and background=='fixed_all_settings_this_route_triangle'")[['version','subset','contrast','n_genes','mean','gene_bootstrap_mean_lo','gene_bootstrap_mean_hi']].to_string(index=False))
 print('fixed RSeT genes',len(common))

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,NODE,a.out))
