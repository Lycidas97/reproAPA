"""Read only frozen all-gene counts on server; return exact group aggregates only."""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
TASK=Path(sys.argv[1]);cfg=json.loads((TASK/'server_source_config.json').read_text())
requests=pd.read_csv(TASK/'gex_requests.tsv.gz',sep='\t');assert requests.group_id.is_unique
expected=cfg['expected_hashes'];base=Path(cfg['frozen_gex'])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
rows=[];sources={};summary=[]
for sid,part in requests.groupby('study_id',sort=True):
 start=time.time();folder='WOT_GSE115943' if sid=='GSE115943' else sid;root=base/folder
 paths=[root/n for n in ['metacells.tsv.gz','genes.tsv.gz','metacell_by_gene_counts.npz']]
 for path in paths:
  assert str(path) in expected,('source_hash_unfrozen',str(path))
  actual=sha(path);assert actual==expected[str(path)],('source_changed',str(path))
  sources[str(path)]={'sha256':actual,'bytes':path.stat().st_size}
 meta=pd.read_csv(paths[0],sep='\t');genes=pd.read_csv(paths[1],sep='\t')
 ids=genes.gene_id.astype(str).str.split('.').str[0];assert ids.is_unique and meta.metacell_id.is_unique
 matrix=sparse.load_npz(paths[2]).tocsr();assert matrix.shape==(len(meta),len(genes))
 assert np.issubdtype(matrix.dtype,np.integer) and (matrix.data>=0).all()
 total=meta.total_gex.to_numpy();assert np.array_equal(np.asarray(matrix.sum(1)).ravel(),total)
 mi=pd.Index(meta.metacell_id);gi=pd.Index(ids)
 cache={}
 for gene,gg in part.groupby('gene_id',sort=True):
  col=gi.get_loc(gene);counts=matrix[:,col].toarray().ravel()
  for r in gg.itertuples():
   mids=json.loads(r.qualified_MC_ids);ix=mi.get_indexer(mids)
   assert len(ix)==r.n_mc and (ix>=0).all() and len(set(mids))==len(mids)
   assert int(meta.n_cells.iloc[ix].sum())==r.n_cells
   num=int(counts[ix].sum());den=int(total[ix].sum());assert den>0
   rows.append(dict(species=r.species,study_id=sid,group_id=r.group_id,gene_id=gene,n_mc=r.n_mc,n_cells=r.n_cells,qualified_MC_ids=r.qualified_MC_ids,GEX_matrix_rows=json.dumps(ix.tolist(),separators=(',',':')),GEX_gene_column=int(col),gene_count_sum=num,all_gene_GEX_sum=den,GEX_CPM=num/den*1e6,source_count_path=str(paths[2]),source_count_sha256=sources[str(paths[2])]['sha256'],source_metacell_path=str(paths[0]),source_metacell_sha256=sources[str(paths[0])]['sha256'],source_gene_path=str(paths[1]),source_gene_sha256=sources[str(paths[1])]['sha256']))
 del matrix
 summary.append({'study_id':sid,'groups':len(part),'genes':part.gene_id.nunique(),'seconds':time.time()-start});print(json.dumps(summary[-1]),flush=True)
out=pd.DataFrame(rows).sort_values('group_id');assert len(out)==len(requests) and set(out.group_id)==set(requests.group_id)
out.to_csv(TASK/'group_GEX_aggregates.tsv.gz',sep='\t',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
(TASK/'GEX_source_manifest.json').write_text(json.dumps({'status':'FROZEN_INTEGER_GROUP_AGGREGATION_PASS','sources':sources,'group_rows':len(out),'study_summary':summary,'request_sha256':sha(TASK/'gex_requests.tsv.gz'),'output_sha256':sha(TASK/'group_GEX_aggregates.tsv.gz')},indent=2)+'\n')
print('GEX_AGGREGATION_COMPLETE',len(out),flush=True)
