"""Read-only raw GEX/PAS extraction under a fixed H1 protocol."""
from pathlib import Path
import gzip, hashlib, io, json, re, sys
import numpy as np
import pandas as pd
from scipy import sparse

S=['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6']
G=['CDK1','TOP2A','NUSAP1','CDC20','CDC45','UBE2C','BIRC5','TPX2']
STUDIES={'GSE103221':('mouse',['CIP']),'GSE115943':('mouse',['WOT']), 'GSE118258':('human',['GSE118258']), 'GSE150311':('human',['GSE150311']), 'GSE178325':('human',['178325_0618','178325_0809']), 'GSE221739':('human',['GSE221739'])}
arrays={}; provenance={}; mapping=[]; panels=[]; controls=[]
def sha(p):
 with Path(p).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def read(p):
 provenance[str(p)]=sha(p);return pd.read_csv(p,sep='\t',low_memory=False)
annotation={}
for species in ['human','mouse']:
 p=Path(CONFIG['resources'][species+'_annotation_gtf']['path']);provenance[str(p)]=sha(p)
 assert provenance[str(p)]==CONFIG['resources'][species+'_annotation_gtf']['sha256']
 d={}
 with gzip.open(p,'rt') as h:
  for line in h:
   if line.startswith('#'):continue
   v=line.rstrip().split('\t')
   if v[2]!='gene':continue
   a=dict(re.findall(r'(\w+) "([^"]*)"',v[8]))
   if 'gene_name' in a:d[a['gene_id'].split('.')[0]]=a['gene_name']
 annotation[species]=d
for sid,(species,routes) in STUDIES.items():
 folder='WOT_GSE115943' if sid=='GSE115943' else sid
 root=Path(CONFIG['native_run'])/'studies'/folder
 meta=read(root/'selected_metacells.tsv.gz')
 meta=meta[meta.route.isin(routes)&meta.node.isin(['S','M','P'])&meta.purity.ge(.8)].reset_index(drop=True)
 full=read(root/'metacells.tsv.gz')
 assert np.array_equal(full.metacell_id.iloc[meta.matrix_row],meta.metacell_id)
 design=read(root/'sampling_design.tsv')
 gm=read(Path(CONFIG['frozen_gex'])/folder/'metacells.tsv.gz'); genes=read(Path(CONFIG['frozen_gex'])/folder/'genes.tsv.gz')
 ids=genes.gene_id.astype(str).str.split('.').str[0]
 assert ids.is_unique
 sy=genes.gene_symbol.fillna('').astype(str).copy();replace=ids.map(annotation[species]);bad=sy.eq('')|sy.str.split('.').str[0].eq(ids)
 sy.loc[bad&replace.notna()]=replace[bad&replace.notna()]
 ix=pd.Index(gm.metacell_id).get_indexer(meta.metacell_id);assert (ix>=0).all()
 xp=Path(CONFIG['frozen_gex'])/folder/'metacell_by_gene_counts.npz';provenance[str(xp)]=sha(xp)
 x=sparse.load_npz(xp).tocsr()[ix];tot=gm.total_gex.to_numpy()[ix]
 assert (x.data>=0).all() and np.issubdtype(x.dtype,np.integer) and np.array_equal(np.asarray(x.sum(1)).ravel(),tot)
 assert np.array_equal(gm.n_cells.iloc[ix].to_numpy(),meta.n_cells.to_numpy())
 lg=x.astype(float);lg.data=np.log1p(lg.data/np.repeat(tot,np.diff(lg.indptr))*1e4)
 pos=[]
 for marker in S+G:
  jj=np.flatnonzero(sy.str.upper().eq(marker));status='unique' if len(jj)==1 else 'missing' if len(jj)==0 else 'ambiguous'
  mapping.append(dict(study_id=sid,species=species,marker=marker,status=status,gene_id=';'.join(ids.iloc[jj]),gene_symbol=';'.join(sy.iloc[jj])))
  assert len(jj)==1,(sid,marker,status)
  pos.append(jj[0])
 marker_counts=x[:,pos].toarray();raw=np.log1p(marker_counts/tot[:,None]*1e4)
 means=np.asarray(lg.mean(0)).ravel();valid=~sy.str.upper().isin(S+G).to_numpy();adj=[];cids=[]
 for k,j in enumerate(pos):
  cand=np.flatnonzero(valid);order=np.lexsort((ids.iloc[cand].to_numpy(),np.abs(means[cand]-means[j])));jj=cand[order[:20]]
  adj.append(raw[:,k]-np.asarray(lg[:,jj].mean(1)).ravel());cids.extend(jj.tolist())
  controls.extend(dict(study_id=sid,marker=(S+G)[k],control_gene_id=ids.iloc[a],control_gene_symbol=sy.iloc[a],global_mean_logCP10K=means[a],marker_global_mean=means[j]) for a in jj)
 adj=np.column_stack(adj)
 meta['Cycle_S']=raw[:,:10].mean(1);meta['Cycle_G2M']=raw[:,10:].mean(1);meta['Cycle']=(meta.Cycle_S+meta.Cycle_G2M)/2
 meta['Adjusted_S']=adj[:,:10].mean(1);meta['Adjusted_G2M']=adj[:,10:].mean(1);meta['Adjusted_Cycle']=(meta.Adjusted_S+meta.Adjusted_G2M)/2
 meta['total_gex']=tot
 arrays[sid+'_meta_json']=np.array(meta.to_json(orient='records',double_precision=15)); arrays[sid+'_design_json']=np.array(design.to_json(orient='records'))
 arrays[sid+'_marker_counts']=marker_counts;arrays[sid+'_total_gex']=tot
 arrays[sid+'_control_counts']=x[:,cids].toarray();arrays[sid+'_control_gene_ids']=ids.iloc[cids].to_numpy(str)
 f=read(root/'features.tsv.gz');f=f[f.keep_for_apa_primary.astype(str).str.lower().eq('true')&f.primary_gene_id.notna()].copy()
 f['position']=np.nan
 shape=f.groupby('primary_gene_id').agg(n=('cut_site','nunique'),c=('chrom','nunique'),s=('strand','nunique'))
 for gene,part in f.groupby('primary_gene_id'):
  a=part.cut_site.to_numpy()*np.where(part.strand.eq('+'),1,-1)
  f.loc[part.index,'position']=(a-a.min())/(a.max()-a.min()) if np.ptp(a)>0 else np.nan
 f=f[f.primary_gene_id.isin(shape.query('n>=2 and c==1 and s==1').index)&f.matched_source_count_50nt.ge(1)].copy()
 groups=[(gene,z) for gene,z in f.groupby('primary_gene_id',sort=True) if len(z)>=2]
 pp=root/'metacell_by_pas_counts.npz';provenance[str(pp)]=sha(pp)
 pas=sparse.load_npz(pp).tocsr()[meta.matrix_row.to_numpy()]
 assert np.issubdtype(pas.dtype,np.integer) and (pas.data>=0).all()
 meta['total_pas']=np.asarray(pas.sum(1)).ravel()
 # Preserve all supported counts for independently replaying fixed panels.
 allcols=np.concatenate([z.feature_index_1based.to_numpy()-1 for _,z in groups]); ap=pas[:,allcols].toarray()
 starts=np.r_[0,np.cumsum([len(z) for _,z in groups])]
 gene_ids=np.array([g for g,z in groups]); positions=np.concatenate([z.position.to_numpy() for g,z in groups])
 gexidx=pd.Index(ids).get_indexer(gene_ids)
 for key,q in meta.groupby(['route','library_key','node'],sort=True):
  route,lib,node=key;rr=q.index.to_numpy();row=dict(species=species,study_id=sid,route=route,library_key=lib,state=node,n_metacells=len(rr),n_cells=int(q.n_cells.sum()))
  if len(rr)<12 or q.n_cells.sum()<30:
   panels.append(dict(row,status='insufficient_MC_for_within_state',n_genes=0));continue
  selected=[]
  for j,(g,z) in enumerate(groups):
   c=ap[rr,starts[j]:starts[j+1]]
   joint=(c.sum(0)>=10)&((c>0).sum(0)>=2)
   if np.all(c.sum(1)>=20) and joint.sum()>=2:selected.append(j)
  if len(selected)<20:
   panels.append(dict(row,status='fewer_than_20_common_genes',n_genes=len(selected)));continue
  if np.any(gexidx[selected]<0):
   panels.append(dict(row,status='panel_GEX_unmapped',n_genes=len(selected)));continue
  cid=f'{sid}_{len(panels):03d}';depth=[];dw=[];thin=[];counts=[];axis=[];pstarts=[0]
  for j in selected:
   c=ap[rr,starts[j]:starts[j+1]];a=positions[starts[j]:starts[j+1]];d=c.sum(1)
   depth.append(d);dw.append(c@a/d);counts.append(c);axis.extend(a);pstarts.append(pstarts[-1]+len(a))
   rng=np.random.default_rng(20261010+int(hashlib.sha256((cid+gene_ids[j]).encode()).hexdigest()[:12],16))
   draws=[]
   for run in range(20):
    left=np.full(len(rr),20);population=d.copy();chosen=[]
    for col in range(c.shape[1]-1):
     value=rng.hypergeometric(c[:,col],population-c[:,col],left)
     chosen.append(value);left-=value;population-=c[:,col]
    chosen.append(left);v=np.column_stack(chosen)
    assert np.all(v.sum(1)==20) and np.all(v<=c)
    draws.append(v@a/20)
   thin.append(np.stack(draws))
  arrays[cid+'_rows']=rr;arrays[cid+'_gene_ids']=gene_ids[selected];arrays[cid+'_DWUI']=np.column_stack(dw);arrays[cid+'_depth']=np.column_stack(depth)
  arrays[cid+'_gex_counts']=x[rr][:,gexidx[selected]].toarray();arrays[cid+'_PAS_counts']=np.column_stack(counts);arrays[cid+'_PAS_positions']=np.array(axis);arrays[cid+'_PAS_starts']=np.array(pstarts)
  arrays[cid+'_thin20_mean']=np.stack(thin,axis=2).mean(2)
  panels.append(dict(row,status='available',n_genes=len(selected),context_id=cid,gene_list_sha256=hashlib.sha256('|'.join(gene_ids[selected]).encode()).hexdigest()))
 arrays[sid+'_meta_json']=np.array(meta.to_json(orient='records',double_precision=15))
 print(f'{sid}: {len(meta)} MC, available panels {sum(z.get("study_id")==sid and z["status"]=="available" for z in panels)}',file=sys.stderr,flush=True)
arrays['mapping_json']=np.array(json.dumps(mapping));arrays['panels_json']=np.array(json.dumps(panels));arrays['controls_json']=np.array(json.dumps(controls));arrays['provenance_json']=np.array(json.dumps(provenance,sort_keys=True));arrays['protocol_sha256']=np.array(PROTOCOL_SHA)
b=io.BytesIO();np.savez_compressed(b,**arrays);sys.stdout.buffer.write(b.getvalue())
