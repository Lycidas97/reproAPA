"""Build/freeze actual contributor covariates and effects before threshold labeling."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p,**k:pd.read_csv(p,sep='\t',float_precision='round_trip',**k)
def save(f,n):f.to_csv(HERE/n,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)
def main():
 freeze=json.loads((HERE/'effect_vector_freeze.json').read_text())
 for p,h in freeze['input_sha256'].items():assert sha(ROOT/p)==h,p
 for n,h in freeze['output_sha256'].items():assert sha(HERE/n)==h,n
 man=json.loads((HERE/'GEX_source_manifest.json').read_text());assert man['status']=='FROZEN_INTEGER_GROUP_AGGREGATION_PASS'
 assert sha(HERE/'gex_requests.tsv.gz')==man['request_sha256'] and sha(HERE/'group_GEX_aggregates.tsv.gz')==man['output_sha256']
 groups=read(HERE/'qualified_groups.tsv.gz');gex=read(HERE/'group_GEX_aggregates.tsv.gz')
 assert groups.group_id.is_unique and gex.group_id.is_unique
 groups=groups.merge(gex,on=['species','study_id','group_id','gene_id','n_mc','n_cells','qualified_MC_ids'],validate='one_to_one')
 assert len(groups)==man['group_rows']
 np.testing.assert_allclose(groups.GEX_CPM,groups.gene_count_sum/groups.all_gene_GEX_sum*1e6,atol=1e-10,rtol=1e-15)
 groups['mean_retained_pas_umi']=groups.qualified_MC_pas_depths.map(lambda v:float(np.mean(json.loads(v))))
 assert groups[['GEX_CPM','mean_retained_pas_umi']].notna().all().all()
 lookup=groups.set_index('group_id')
 pairs=read(HERE/'matched_pair_effects.tsv.gz');routes=read(HERE/'matched_route_effects.tsv.gz')
 for name,source in [('GEX_CPM','GEX_CPM'),('APA_coverage','mean_retained_pas_umi')]:
  for side,group_col in [('a','target_group'),('b','reference_group')]:pairs[name+'_'+side]=lookup[source].reindex(pairs[group_col]).to_numpy()
 pairs['mean_GEX']=(pairs.GEX_CPM_a+pairs.GEX_CPM_b)/2;pairs['APA_coverage']=(pairs.APA_coverage_a+pairs.APA_coverage_b)/2
 pairs['baseline']=np.array([np.array(json.loads(r.reference_usage))@np.array(json.loads(r.position)) for r in pairs.itertuples()])
 pairs['PAS_number']=pairs.PAS_ids.map(lambda x:len(json.loads(x)))
 for name in ['mean_GEX','APA_coverage','baseline','PAS_number']:assert pairs[name].notna().all() and np.isfinite(pairs[name]).all()
 cov=[];route_index=routes.set_index(['study_id','route','gene_id','contrast'])
 for key,part in pairs.groupby(['study_id','route','gene_id','contrast'],sort=True):
  r=route_index.loc[key]
  part=part.sort_values('library_key');ww=part.weight.to_numpy();assert abs(ww.sum()-1)<1e-12
  assert ww.tolist()==json.loads(r.weights) and part.library_key.tolist()==json.loads(r.contributor_libraries)
  row=r.to_dict();row.update(dict(zip(['study_id','route','gene_id','contrast'],key)))
  for n in ['mean_GEX','APA_coverage','baseline']:row[n]=float(ww@part[n])
  assert part.PAS_number.nunique()==1;row['PAS_number']=int(part.PAS_number.iloc[0]);row['gene_symbol']=part.gene_symbol.iloc[0]
  pos=np.array(json.loads(part.position.iloc[0]));row['position']=part.position.iloc[0];row['retained_position_min']=float(pos.min());row['retained_position_max']=float(pos.max())
  assert abs(row['baseline']-np.array(json.loads(r.reference_usage))@pos)<1e-12
  cov.append(row)
 route=pd.DataFrame(cov)
 records=[]
 for key,part in route.groupby(['species','contrast','study_id','gene_id'],sort=True):
  assert part.route.is_unique
  row=dict(zip(['species','source_contrast','unit','gene_id'],key));row.update(contrast='P_M',gene_symbol=part.gene_symbol.iloc[0],n_routes=len(part),series_pattern=';'.join(sorted(part.route)),metric='WPP',effect_mean=float(part.delta_DWUI.mean()),effect_min=float(part.delta_DWUI.min()),effect_max=float(part.delta_DWUI.max()),baseline=float(part.baseline.median()),baseline_min=float(part.baseline.min()),baseline_max=float(part.baseline.max()),capacity_negative=float(part.baseline.min()),capacity_positive=float((1-part.baseline).min()),capacity_unsigned=np.nan)
  for n in ['mean_GEX','APA_coverage','PAS_number']:row[n]=float(part[n].median())
  row['supported_catalogue_PAS_number']=row['PAS_number'];records.append(row)
 studies=pd.DataFrame(records)
 assert studies[['effect_min','effect_max','baseline','capacity_negative','capacity_positive','mean_GEX','APA_coverage','PAS_number']].notna().all().all()
 assert studies.groupby(['species','source_contrast']).unit.nunique().max()<=3
 # Qualification, not direction, must reproduce all 728 pair identities for each source contrast.
 bg=read(HERE/'fixed_background.tsv')
 for sp in ['human','mouse']:
  for contrast in ['D00','D01']:
   p=studies[studies.species.eq(sp)&studies.source_contrast.eq(contrast)]
   measured=p.groupby('gene_id').unit.nunique();assert measured.ge(2).all()
   assert set(measured.index)==set(bg[sp+'_gene_id'])
 save(groups,'qualified_group_covariates.tsv.gz');save(pairs,'pair_covariates.tsv.gz');save(route,'route_covariates.tsv.gz');save(studies,'study_covariates.tsv.gz')
 inputs=['effect_vector_freeze.json','GEX_source_manifest.json','group_GEX_aggregates.tsv.gz','qualified_groups.tsv.gz','matched_pair_effects.tsv.gz','matched_route_effects.tsv.gz','fixed_background.tsv']
 outputs=['qualified_group_covariates.tsv.gz','pair_covariates.tsv.gz','route_covariates.tsv.gz','study_covariates.tsv.gz']
 result={'status':'ACTUAL_COVARIATES_AND_EFFECTS_FROZEN_BEFORE_LABELS','background':len(bg),'group_rows':len(groups),'pair_rows':len(pairs),'route_rows':len(route),'study_rows':len(studies),'qualified_MC_GEX_actual_integer_source':True,'input_sha256':{n:sha(HERE/n) for n in inputs},'output_sha256':{n:sha(HERE/n) for n in outputs},'new_direction_labels_computed':False}
 (HERE/'covariate_freeze.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','background','group_rows','pair_rows','route_rows','study_rows']}),flush=True)
if __name__=='__main__':main()
