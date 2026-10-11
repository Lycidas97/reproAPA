"""Independent raw-subset arithmetic and frozen-metadata checks."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];B=ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
S=['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6'];G=['CDK1','TOP2A','NUSAP1','CDC20','CDC45','UBE2C','BIRC5','TPX2']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 x=np.load(HERE/'inputs.npz',allow_pickle=False);old=np.load(ROOT/'notes/2026-09-26-state-comparability-audit/expression_inputs.npz',allow_pickle=False)
 assert sha(HERE/'PROTOCOL.zh.md')==str(x['protocol_sha256'])
 checks=[];panels=pd.read_csv(HERE/'panel_eligibility.tsv',sep='\t');f=pd.read_csv(HERE/'context_associations.tsv',sep='\t');out=pd.read_csv(HERE/'metacell_apa_cycle.tsv.gz',sep='\t');roster=pd.read_csv(HERE/'fixed_panel_genes.tsv.gz',sep='\t')
 mapping=pd.read_csv(HERE/'marker_mapping.tsv',sep='\t')
 assert mapping.status.eq('unique').all() and len(mapping)==108
 orth=pd.read_csv(B/'cross_overlap/ortholog_matched_features.tsv.gz',sep='\t',usecols=['human_gene_id','mouse_gene_id'])
 oset=set(map(tuple,orth[['human_gene_id','mouse_gene_id']].to_numpy()))
 marker_pairs=[]
 for marker in S+G:
  h=mapping[mapping.species.eq('human')&mapping.marker.eq(marker)].gene_id.unique();m=mapping[mapping.species.eq('mouse')&mapping.marker.eq(marker)].gene_id.unique()
  assert len(h)==len(m)==1 and (h[0],m[0]) in oset
  marker_pairs.append(dict(marker=marker,human_gene_id=h[0],mouse_gene_id=m[0],pinned_one_to_one=True))
 pd.DataFrame(marker_pairs).to_csv(HERE/'cycle_marker_orthologs.tsv',sep='\t',index=False)
 checks.append(dict(check='18_cycle_markers_pinned_one_to_one',n=18,max_abs_error=0,passed=True))
 oldchecks=[]
 for k in x.files:
  if not k.endswith('_meta_json'):continue
  sid=k[:-10];m=pd.DataFrame(json.loads(str(x[k])));c=x[sid+'_marker_counts'];tot=x[sid+'_total_gex']
  assert (m.groupby(['route','library_key','node']).experimental_time.nunique(dropna=False)==1).all()
  raw=np.log1p(c/tot[:,None]*1e4);adj=raw-np.log1p(x[sid+'_control_counts']/tot[:,None]*1e4).reshape(len(m),18,20).mean(2)
  for name,v in [('Cycle_S',raw[:,:10].mean(1)),('Cycle_G2M',raw[:,10:].mean(1)),('Adjusted_S',adj[:,:10].mean(1)),('Adjusted_G2M',adj[:,10:].mean(1))]:
   assert np.allclose(m[name],v,rtol=0,atol=1e-12)
   checks.append(dict(check=sid+'_'+name+'_from_integer_GEX',n=len(m),max_abs_error=float(np.max(np.abs(m[name]-v))),passed=True))
  coverage=pd.read_csv(B/'metacell'/('WOT_'+sid if sid=='GSE115943' else sid)/'metacell_pas_coverage.tsv.gz',sep='\t').set_index('metacell_id').reindex(m.metacell_id)
  assert np.array_equal(coverage.total_pas_umi,m.total_pas)
  assert np.array_equal(coverage.n_cells,m.n_cells)
  checks.append(dict(check=sid+'_frozen_PAS_total_cells_identity',n=len(m),max_abs_error=0,passed=True))
  if sid+'_meta_json' in old.files:
   om=pd.DataFrame(json.loads(str(old[sid+'_meta_json'])));ix=pd.Index(om.metacell_id).get_indexer(m.metacell_id);assert np.all(ix>=0)
   sy=old[m.species.iloc[0]+'_marker_symbols'];jx=pd.Index(sy).get_indexer(S+G);assert np.all(jx>=0)
   z=old[sid+'_marker_logcp10k'][ix][:,jx];e=np.abs(z-raw)
   oldchecks.append(dict(study_id=sid,n_marker_values=e.size,max_abs_error=float(e.max()),match=bool(np.allclose(z,raw,atol=1e-12,rtol=0))))
   # The frozen original marker subset should be exactly recoverable.
   assert np.allclose(z,raw,atol=1e-12,rtol=0)
 axes=[]
 for chunk in pd.read_csv(B/'state_vectors/PAS_axes.primary.tsv.gz',sep='\t',usecols=['study_id','setting','primary_gene_id','feature_index_1based','position'],chunksize=100000):
  axes.append(chunk[chunk.setting.eq('supported50_p08')])
 axes=pd.concat(axes,ignore_index=True)
 axis_groups={(sid,g):z.sort_values('feature_index_1based').position.to_numpy() for (sid,g),z in axes.groupby(['study_id','primary_gene_id'])}
 raw_axes=pd.read_csv(HERE/'raw_primary_axis_features.tsv.gz',sep='\t')
 raw_groups={}
 axis_recovery=[]
 for (sid,g),z in raw_axes.groupby(['study_id','primary_gene_id']):
  assert z.chrom.nunique()==1 and z.strand.nunique()==1 and z.cut_site.nunique()>=2
  sign=1 if z.strand.iloc[0]=='+' else -1;coord=z.cut_site.to_numpy()*sign
  z=z.copy();z['recovered_position']=(coord-coord.min())/(coord.max()-coord.min())
  z=z[z.matched_source_count_50nt.ge(1)].sort_values('feature_index_1based')
  assert len(z)>=2
  raw_groups[sid,g]=z.recovered_position.to_numpy()
  covered=(sid,g) in axis_groups
  if covered:assert np.allclose(axis_groups[sid,g],raw_groups[sid,g],atol=1e-12,rtol=0)
  axis_recovery.append(dict(study_id=sid,gene_id=g,n_supported_PAS=len(z),n_full_primary_PAS=len(coord),current_endpoint_axis_table_covered=covered,raw_primary_coordinate_reconstruction_passed=True))
 pd.DataFrame(axis_recovery).to_csv(HERE/'primary_axis_recovery.tsv',sep='\t',index=False)
 for r in panels[panels.status.eq('available')].itertuples():
  cid=r.context_id;dw=x[cid+'_DWUI'];dep=x[cid+'_depth'];counts=x[cid+'_PAS_counts'];ax=x[cid+'_PAS_positions'];start=x[cid+'_PAS_starts'];genes=x[cid+'_gene_ids']
  assert dw.shape==(r.n_metacells,r.n_genes) and np.all(dep>=20)
  assert np.array_equal(genes,roster[roster.context_id.eq(cid)].gene_id.to_numpy())
  error=0.
  for j,g in enumerate(genes):
   c=counts[:,start[j]:start[j+1]];position=ax[start[j]:start[j+1]]
   assert np.allclose(position,raw_groups[r.study_id,g],atol=1e-12,rtol=0),(r.study_id,g)
   assert np.array_equal(c.sum(1),dep[:,j]) and ((c.sum(0)>=10)&((c>0).sum(0)>=2)).sum()>=2
   error=max(error,float(np.max(np.abs((c@position)/c.sum(1)-dw[:,j]))))
  assert error<1e-12
  z=out[out.context_id.eq(cid)].reset_index(drop=True);m=pd.DataFrame(json.loads(str(x[r.study_id+'_meta_json']))).iloc[x[cid+'_rows']].reset_index(drop=True)
  assert np.array_equal(z.metacell_id,m.metacell_id)
  assert np.allclose(z.mean_DWUI,dw.mean(1),atol=1e-12,rtol=0)
  assert np.allclose(z.median_DWUI,np.median(dw,axis=1),atol=1e-12,rtol=0)
  assert np.allclose(z.thin20_mean_DWUI,x[cid+'_thin20_mean'].mean(0),atol=1e-12,rtol=0)
  assert np.allclose(z.panel_GEX,np.log1p(x[cid+'_gex_counts']/m.total_gex.to_numpy()[:,None]*1e4).mean(1),atol=1e-12,rtol=0)
  checks.append(dict(check=cid+'_fixed_roster_raw_DWUI_depth',n=dw.size,max_abs_error=error,passed=True))
 # Library equal means must not weight prolific MC libraries more heavily.
 library=pd.read_csv(HERE/'library_associations.tsv',sep='\t');route=pd.read_csv(HERE/'route_associations.tsv',sep='\t')
 for r in route[route.status.eq('descriptive')].itertuples():
  z=library[library.study_id.eq(r.study_id)&library.route.eq(r.route)&library.score.eq(r.score)&library.APA_variant.eq(r.APA_variant)&library.scope.eq(r.scope)]
  assert z.library_key.is_unique and len(z)==r.n_libraries
  assert abs(z.partial_correlation.mean()-r.partial_correlation)<1e-12
 checks.append(dict(check='valid_library_rows_are_descriptive',n=len(library),max_abs_error=0,passed=True))
 checks.append(dict(check='all_route_library_equal_means',n=int(route.status.eq('descriptive').sum()),max_abs_error=0,passed=True))
 # Any failed context stays visible and never acquires an inferential zero.
 assert library.status.eq('descriptive').all()
 assert f.loc[f.status.ne('available'),'partial_correlation'].isna().all()
 checks.append(dict(check='unavailable_never_zero_effect',n=int(f.status.ne('available').sum()),max_abs_error=0,passed=True))
 pd.DataFrame(checks).to_csv(HERE/'verification.tsv',sep='\t',index=False);pd.DataFrame(oldchecks).to_csv(HERE/'frozen_marker_recovery.tsv',sep='\t',index=False)
 audit=dict(status='PASS',checks=len(checks),raw_or_frozen_values_checked=int(sum(z['n'] for z in checks)+sum(z['n_marker_values'] for z in oldchecks)),frozen_marker_studies=len(oldchecks),frozen_marker_max_error=max(z['max_abs_error'] for z in oldchecks),failed_contexts_explicit=True,all_library_state_dates_unique=True,raw_primary_axis_genes_checked=len(raw_groups),current_endpoint_axis_genes_missing=sum(not z['current_endpoint_axis_table_covered'] for z in axis_recovery),biological_replication_claim=False,protocol_sha256=sha(HERE/'PROTOCOL.zh.md'),tests_skipped=0)
 (HERE/'validation.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
