"""Frozen-taxonomy diagnostics and explicit panel eligibility for Figure 2."""
from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import analyze_figure2_dynamic_architecture_v1 as base
P=base.P
O=P/'outputs/figure2-core-supplement-v1'
O.mkdir(exist_ok=True)
D=base.OUT
members=pd.read_csv(base.MEMBERS,sep='\t')
source=pd.read_csv(base.APA,sep='\t')
matched=pd.read_csv(base.MATCHED,sep='\t')
names=pd.read_csv(P/'outputs/figure2-dynamic-architecture-layout-v28/dataset_label_mapping.tsv',sep='\t')
names.to_csv(O/'dataset_label_mapping.tsv',sep='\t',index=False)
def save(d,n):d.to_csv(O/n,sep='\t',index=False)
# Confirm the old centered sensitivity corresponds to the current k=3 membership.
old=P/'outputs/figure2-supplement-layout-v5'
for n in ['resolution','centroid']:
 shutil.copy2(old/f'S2_{n}_source_k2-k6.tsv',O/f'S2_{n}.tsv')
cent=pd.read_csv(O/'S2_centroid.tsv',sep='\t')
expected=members.groupby(['cluster','species']).size().unstack()
assert np.array_equal(cent[cent.k.eq(3)].set_index('cluster')[['mouse','human']],expected[['mouse','human']])
z=np.load(base.FIT,allow_pickle=True);x=z['features'];sp=z['species'].astype(str);ids=z['gene_ids'].astype(str);w=z['weights']
ks=range(2,7)
centers={k:np.array([np.average(x[z[f'labels_k{k}']==j],axis=0,weights=w[z[f'labels_k{k}']==j]) for j in range(k)]) for k in ks}
def record(k,kind,rep,xx,labels,ix=None):
 ix=np.arange(len(x)) if ix is None else ix
 ww=np.array([.5/np.sum(sp[ix]==s) for s in sp[ix]])
 cc=np.array([np.average(xx[labels==j],axis=0,weights=ww[labels==j]) for j in range(k)])
 a,b=linear_sum_assignment(((centers[k][:,None]-cc[None])**2).sum(2))
 return dict(k=k,perturbation=kind,replicate=rep,ari=adjusted_rand_score(z[f'labels_k{k}'][ix],labels),centroid_rmse=float(np.sqrt(np.mean((centers[k][a]-cc[b])**2))))
records=[]
for k in ks:
 for rep,labels in enumerate(z[f'subsamples_k{k}']):
  ix=np.flatnonzero(labels>=0);records.append(record(k,'Gene subsampling',rep,x[ix],labels[ix],ix))
libs={};weights={}
for s in ['mouse','human']:
 a=np.load(P/f'data/apa-srep-common-v1/{s}/library_inputs.npz')
 assert np.array_equal(a['gene_ids'].astype(str),ids[sp==s])
 libs[s]=([a[f'values_{i}'] for i in range(len(base.ROUTES[s]))],[a[f'indices_{i}'] for i in range(len(base.ROUTES[s]))])
 weights[s]=np.load(P/f'data/apa-srep-common-v1/{s}/validation.npz')['weights']
for rep in range(32):
 blocks=[]
 for s in ['mouse','human']:
  arrays,indices=libs[s];views=[]
  for a,ix in zip(arrays,indices):
   ww=weights[s][rep,ix,None,None];finite=np.isfinite(a);den=(finite*ww).sum(0)
   v=np.divide((np.where(finite,a,0)*ww).sum(0),den,out=np.full(den.shape,np.nan),where=den>0)
   views.append(base.unit_centered(v))
  blocks.append(np.mean(views,axis=0))
 xx=np.concatenate(blocks)
 for k in ks:records.append(record(k,'Library reweighting',rep,xx,z[f'library_labels_k{k}'][rep]))
blocks=[]
for s in ['mouse','human']:
 views=[]
 for r in base.ROUTES[s]:
  a=source.query("layer=='expanded_primary'").loc[lambda a:a.species.eq(s)&a.route.eq(r)].pivot(index='gene_id',columns='node',values='position_index').reindex(index=ids[sp==s],columns=list('SREP')).to_numpy()
  views.append(base.unit_centered(a))
 blocks.append(np.mean(views,axis=0))
expanded=np.concatenate(blocks)
for k in ks:
 records.append(record(k,'Expanded PAS set',0,expanded,z[f'expanded_labels_k{k}']))
 model=KMeans(k,random_state=base.SEED+100,n_init=50,max_iter=500,tol=1e-7).fit(x,sample_weight=w)
 records.append(record(k,'Alternative seed',0,x,model.labels_))
stability=pd.DataFrame(records);save(stability,'S2_stability.tsv')
resolution=pd.read_csv(O/'S2_resolution.tsv',sep='\t')
for kind,col in [('Gene subsampling','subsample_median_ari'),('Library reweighting','library_median_ari'),('Expanded PAS set','expanded_ari'),('Alternative seed','different_seed_ari')]:
 assert np.allclose(stability[stability.perturbation.eq(kind)].groupby('k').ari.median(),resolution[col],atol=1e-10),(kind,col)
for file in ['leave_one_route_out_summary.tsv','leave_one_route_out_cluster_metrics.tsv']:
 shutil.copy2(D/file,O/file)
# Full distributions and threshold grids: fixed cohort and class labels.
metrics=pd.read_csv(D/'trajectory_gene_route_metrics.tsv.gz',sep='\t')
expanded_metrics=base.trajectory_architecture(members,source.assign(layer=source.layer.replace({'supported50':'unused','expanded_primary':'supported50'})))[0]
metrics['layer']='supported50';expanded_metrics['layer']='expanded_primary'
allmetrics=pd.concat([metrics,expanded_metrics],ignore_index=True)
allmetrics['SREP_amplitude']=allmetrics[[f'delta_{s}' for s in 'SREP']].max(axis=1)-allmetrics[[f'delta_{s}' for s in 'SREP']].min(axis=1)
save(allmetrics,'S3_gene_route_metrics.tsv.gz')
rows=[]
for (layer,spc,r,c),a in allmetrics.groupby(['layer','species','route','cluster']):
 for metric in ['path_length_L','SREP_amplitude']:
  for cutoff in np.round(np.arange(0,.201,.02),2):
   v=a[a[metric].ge(cutoff)]
   rows.append(dict(layer=layer,species=spc,route=r,cluster=c,threshold_metric=metric,cutoff=cutoff,n=len(v),n_initial=len(a),retained_fraction=len(v)/len(a),median_T=v.transient_fraction_T.median(),fraction_T_ge_0_5=v.transient_fraction_T.ge(.5).mean(),fraction_T_ge_0_75=v.transient_fraction_T.ge(.75).mean()))
save(pd.DataFrame(rows),'S3_threshold_sensitivity.tsv')
# Recompute Pearson and rank correspondence on exactly matched states.
genes=pd.read_csv(D/'APA_GEX_gene_route_correlations.tsv.gz',sep='\t')
m=matched.query("layer=='supported50'").drop(columns=['cluster','row_order'],errors='ignore').merge(members[['species','gene_id','cluster']],on=['species','gene_id'],validate='many_to_one')
wide=m[m.node.isin(list('SREP'))].pivot(index=['species','route','cluster','gene_id'],columns='node',values=['position_index','gex_log2_cpm1'])
rr=[]
for key,g in wide.iterrows():
 a=g['position_index'].reindex(list('SREP'));e=g['gex_log2_cpm1'].reindex(list('SREP'));ok=a.notna()&e.notna()
 if not ok.S or ok.sum()<3:continue
 aa=a[ok]-a.S;ee=e[ok]-e.S
 rho=spearmanr(aa,ee).statistic if np.ptp(aa)>1e-12 and np.ptp(ee)>1e-12 else np.nan
 rr.append(dict(zip(['species','route','cluster','gene_id'],key),spearman_rho=rho))
genes=genes.merge(pd.DataFrame(rr),on=['species','route','cluster','gene_id'],validate='one_to_one');save(genes,'S4_gene_correspondence.tsv.gz')
trans=pd.read_csv(D/'APA_GEX_transition_concordance.tsv.gz',sep='\t')
tr=trans.groupby(['species','route','cluster','transition']).agg(n=('gene_id','size'),n_genes=('gene_id','nunique'),same=('concordance',lambda s:s.eq('same').sum()),zero=('concordance',lambda s:s.eq('zero').sum())).reset_index();tr['same_fraction']=tr.same/tr.n
save(tr,'S4_transition_diagnostics.tsv')
counts=genes.groupby(['species','route','cluster','n_states']).agg(n=('r','count'),median_Pearson=('r','median'),median_Spearman=('spearman_rho','median')).reset_index();save(counts,'S4_state_count_summary.tsv')
# All 17 groups audited against their own APA and paired expression availability.
heat=pd.read_csv(P/'outputs/figure2-k3-srelative-v2/panel_AB_heat_values.tsv.gz',sep='\t')
elig=[]
for row in names.itertuples():
 r=row.route;spc=row.species
 a=source[source.layer.eq('supported50')&source.route.eq(r)&source.species.eq(spc)].merge(members[['species','gene_id']],on=['species','gene_id'])
 paired=m[m.route.eq(r)&m.species.eq(spc)]
 valid=paired[paired.position_index.notna()&paired.gex_log2_cpm1.notna()]
 h=heat[heat.route.eq(r)&heat.species.eq(spc)]
 g=genes[genes.route.eq(r)&genes.species.eq(spc)]
 t=trans[trans.route.eq(r)&trans.species.eq(spc)]
 obs=''.join(s for s in 'SREPM' if a[a.node.eq(s)].position_index.notna().any())
 joint=''.join(s for s in 'SREPM' if valid.node.eq(s).any())
 nA=h.loc[h.display_status.eq('estimable'),'gene_id'].nunique()
 nBC=metrics[metrics.route.eq(r)&metrics.species.eq(spc)].gene_id.nunique()
 nD=g.r.notna().sum();nE=t.gene_id.nunique()
 reason='eligible' if nD else ('S unavailable' if 'S' not in joint else 'fewer than 3 jointly observed SREP states' if len(set(joint)&set('SREP'))<3 else 'no nonconstant eligible gene profiles')
 elig.append(dict(species=spc,route=r,label=row.figure_label,APA_states=obs,paired_APA_GEX_states=joint,A_genes=nA,BC_genes=nBC,D_genes=int(nD),E_genes=nE,E_transitions=len(t),D_status=reason,E_status='eligible' if len(t) else 'no jointly observed adjacent SREP transition',GEX_source_present=bool(len(paired))))
elig=pd.DataFrame(elig);save(elig,'Table_S_state_analysis_group_eligibility.tsv')
# Outer join preserves classes/routes eligible for only one endpoint or GEX statistic.
summ=pd.read_csv(D/'trajectory_route_cluster_summary.tsv',sep='\t')
th=pd.read_csv(D/'endpoint_information_loss_threshold_summary.tsv',sep='\t')
summ=summ.merge(th[['species','route','cluster','fraction_T_ge_0_5','fraction_T_ge_0_75']],on=['species','route','cluster'])
gsum=genes.groupby(['species','route','cluster']).agg(n_corr=('r','count'),median_Pearson_r=('r','median'),median_Spearman_rho=('spearman_rho','median')).reset_index()
tsum=trans.groupby(['species','route','cluster']).agg(n_transitions=('concordance','size'),n_same=('concordance',lambda a:a.eq('same').sum())).reset_index();tsum['same_fraction']=tsum.n_same/tsum.n_transitions
summ=summ.merge(gsum,on=['species','route','cluster'],how='outer').merge(tsum,on=['species','route','cluster'],how='outer');save(summ,'Table_S_route_cluster_numerical_summary.tsv')
# Existing GO audit retained as tables/notes, with no new main biological claim.
go=P/'outputs/figure2-k3-cross-species-go-robustness-v1'
(O/'GO_audit').mkdir(exist_ok=True)
for f in go.iterdir():
 if f.suffix in ['.tsv','.gz','.md','.json']:shutil.copy2(f,O/'GO_audit'/f.name)
qa=dict(status='PASS',taxonomy='existing mean-centered unit-norm joint k=3; S-relative display only',n_metrics=len(metrics),n_expanded_metrics=len(expanded_metrics),n_correlations=int(genes.r.notna().sum()),n_transitions=len(trans),eligibility=elig.to_dict('records'),stability_medians_match_frozen=True,LOO_scope='conditional on frozen candidate cohort; label matching to full taxonomy; not a de novo selection validation')
(O/'analysis_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
print(elig.to_string(index=False))
print('Prepared',O)
