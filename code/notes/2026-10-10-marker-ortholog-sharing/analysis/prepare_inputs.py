"""Freeze matched effect vectors and exact qualified-MC requests before labels/null."""
import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
MARKER=ROOT/'notes/2026-10-10-marker-reference-sensitivity';OLD=ROOT/'notes/2026-10-10-ortholog-threshold-scan'
Q=HERE.parent/'qualification/qualification_members.tsv';KEY=['study_id','route','gene_id','library_key','setting','version','contrast']
SP={'GSE103221':'mouse','GSE115943':'mouse','GSE118258':'human','GSE150311':'human','GSE221739':'human'}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p,**k:pd.read_csv(p,sep='\t',float_precision='round_trip',**k)
def save(f,name):f.to_csv(HERE/name,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
def main():
 protocol=json.loads((HERE.parent/'protocol.json').read_text());pf=json.loads((HERE.parent/'protocol_freeze_sha256.json').read_text())
 for n,h in pf.items():assert sha(HERE.parent/n)==h,n
 sources=[HERE.parent/'PROTOCOL.zh.md',HERE.parent/'protocol.json',Q,MARKER/'unique_support_bridge_route_effects.tsv.gz',MARKER/'unique_support_pair_weights.tsv.gz',MARKER/'pair_effects.tsv.gz',MARKER/'gene_group_profiles.tsv.gz',MARKER/'primary_axis_features.tsv.gz',OLD/'fixed_ortholog_background.tsv',MARKER/'remote_input_sha256.json',ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json']
 hashes={str(p.relative_to(ROOT)):sha(p) for p in sources}
 q=read(Q);q=q[q.version.eq(protocol['marker_version'])&q.setting.eq(protocol['setting'])&q.old_background.eq('P_M')]
 assert len(q)==751
 old=read(OLD/'fixed_ortholog_background.tsv');old=old[old.contrast.eq('P_M')].drop(columns='contrast')
 bg=q[q.both_species_ge2].merge(old,on=['human_gene_id','mouse_gene_id'],validate='one_to_one').sort_values('human_gene_id').reset_index(drop=True)
 assert len(bg)==protocol['background_n']==728
 save(bg,'fixed_background.tsv')
 gids={sp:set(bg[sp+'_gene_id']) for sp in ['human','mouse']}
 routes=read(sources[3]);routes=routes[routes.version.eq(protocol['marker_version'])&routes.setting.eq(protocol['setting'])&routes.background.eq(protocol['bridge_background'])&routes.contrast.isin(['D00','D01'])].copy()
 routes['species']=routes.study_id.map(SP);routes=routes[[r.gene_id in gids.get(r.species,set()) for r in routes.itertuples()]].copy()
 assert routes.species.notna().all() and not routes.duplicated(['study_id','route','gene_id','contrast']).any()
 weights=read(sources[4]);weights=weights[weights.version.eq(protocol['marker_version'])&weights.setting.eq(protocol['setting'])&weights.background.eq(protocol['bridge_background'])&weights.contrast.isin(['D00','D01'])].copy()
 weights=weights[weights.study_id.isin(SP)]
 weights=weights[[r.gene_id in gids[SP[r.study_id]] for r in weights.itertuples()]].copy()
 assert not weights.duplicated(KEY).any()
 pairs=read(sources[5]);pairs=pairs.merge(weights[KEY+['weight']],on=KEY,validate='one_to_one').copy()
 assert len(pairs)==len(weights) and pairs.status.eq('descriptive').all(),pairs.status.unique()
 for _,r in routes.iterrows():
  pp=pairs[pairs.study_id.eq(r.study_id)&pairs.route.eq(r.route)&pairs.gene_id.eq(r.gene_id)&pairs.contrast.eq(r.contrast)].sort_values('library_key')
  assert len(pp)==r.n_libraries
  assert pp.library_key.tolist()==json.loads(r.contributor_libraries)
  np.testing.assert_allclose(pp.weight,np.array(json.loads(r.weights)),atol=0,rtol=0)
  a=pp.weight.to_numpy()@np.stack(pp.target_usage.map(json.loads));b=pp.weight.to_numpy()@np.stack(pp.reference_usage.map(json.loads))
  np.testing.assert_allclose(a,json.loads(r.target_usage),atol=1e-12,rtol=0);np.testing.assert_allclose(b,json.loads(r.reference_usage),atol=1e-12,rtol=0)
  pos=np.array(json.loads(pp.position.iloc[0]));assert abs((a-b)@pos-r.delta_DWUI)<1e-12
  assert all(v==pp.position.iloc[0] for v in pp.position) and all(v==r.PAS_ids for v in pp.PAS_ids)
 # Same P and contributor library IDs hold before computing any recurrence classes.
 z=pairs.set_index(KEY[:-1]+['contrast']).sort_index()
 for _,f in pairs.groupby(KEY[:-1]):
  assert set(f.contrast)=={'D00','D01'}
  assert f.target_group.nunique()==1 and f.position.nunique()==1
  aa,bb=f.sort_values('contrast').target_usage.map(json.loads)
  np.testing.assert_allclose(aa,bb,atol=1e-12,rtol=0)
 group_ids=set(pairs.target_group)|set(pairs.reference_group)
 cols=['species','study_id','route','library_key','setting','gene_id','group_id','cohort','status','n_mc','n_cells','qualified_MC_ids','qualified_MC_pas_depths','PAS_ids','position','pooled_molecules','usage']
 collected=[]
 for chunk in read(sources[6],usecols=cols,chunksize=100000):
  selected=chunk[chunk.group_id.isin(group_ids)]
  if len(selected):collected.append(selected)
 groups=pd.concat(collected,ignore_index=True).sort_values('group_id')
 assert groups.group_id.is_unique and set(groups.group_id)==group_ids and groups.status.eq('available').all()
 for r in groups.itertuples():
  ids=json.loads(r.qualified_MC_ids);depths=json.loads(r.qualified_MC_pas_depths)
  assert len(ids)==len(set(ids))==r.n_mc and len(depths)==r.n_mc and min(depths)>=20
  assert sum(depths)==sum(json.loads(r.pooled_molecules))
 save(groups,'qualified_groups.tsv.gz');save(pairs,'matched_pair_effects.tsv.gz');save(routes,'matched_route_effects.tsv.gz')
 request=groups[['species','study_id','group_id','gene_id','n_mc','n_cells','qualified_MC_ids']].copy();save(request,'gex_requests.tsv.gz')
 cfg=json.loads(sources[-1].read_text());remote=json.loads((MARKER/'remote_input_sha256.json').read_text())
 (HERE/'server_source_config.json').write_text(json.dumps({'frozen_gex':cfg['frozen_gex'],'expected_hashes':remote},indent=2)+'\n')
 result={'status':'MATCHED_VECTOR_AND_REQUEST_FREEZE_BEFORE_LABELS','background':len(bg),'pair_rows':len(pairs),'route_rows':len(routes),'group_rows':len(groups),'input_sha256':hashes,'output_sha256':{n:sha(HERE/n) for n in ['fixed_background.tsv','matched_pair_effects.tsv.gz','matched_route_effects.tsv.gz','qualified_groups.tsv.gz','gex_requests.tsv.gz','server_source_config.json']},'new_labels_computed':False}
 (HERE/'effect_vector_freeze.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','background','pair_rows','route_rows','group_rows']}),flush=True)
if __name__=='__main__':main()
