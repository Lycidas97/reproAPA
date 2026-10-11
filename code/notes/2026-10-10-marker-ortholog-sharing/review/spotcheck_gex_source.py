"""Fifteen fixed-order source-matrix checks; remote reads only, scalar returns."""
from pathlib import Path
import base64,hashlib,json,subprocess

import pandas as pd

HERE=Path(__file__).resolve().parent
ANALYSIS=HERE.parent/'analysis'
requests=pd.read_csv(ANALYSIS/'gex_requests.tsv.gz',sep='\t')
picked=[]
for sid,part in requests.groupby('study_id',sort=True):
    q=part.sort_values('group_id').reset_index(drop=True)
    assert len(q)>3
    for i in [0,len(q)//2,len(q)-1]:picked.append(q.iloc[i].to_dict())
assert len(picked)==15
payload={'configuration':json.loads((ANALYSIS/'server_source_config.json').read_text()),'requests':picked}
(HERE/'fixed15_GEX_requests.json').write_text(json.dumps(payload,indent=2)+'\n')
encoded=base64.b64encode(json.dumps(payload).encode()).decode()
remote=r'''
from pathlib import Path
import base64,hashlib,json
import numpy as np,pandas as pd
from scipy import sparse
task=json.loads(base64.b64decode('PAYLOAD'))
cfg=task['configuration'];req=pd.DataFrame(task['requests'])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()
out={'status':'PASS','remote_writes':False,'selection':'group_id sorting: first, middle, last per study','source_matrix_downloads':False,'results':[],'source_sha256':{}}
for sid,g in req.groupby('study_id',sort=True):
 folder=Path(cfg['frozen_gex'])/('WOT_GSE115943' if sid=='GSE115943' else sid)
 paths=[folder/n for n in ['metacells.tsv.gz','genes.tsv.gz','metacell_by_gene_counts.npz']]
 for p in paths:
  h=sha(p);assert h==cfg['expected_hashes'][str(p)];out['source_sha256'][str(p)]=h
 meta=pd.read_csv(paths[0],sep='\t');genes=pd.read_csv(paths[1],sep='\t')
 matrix=sparse.load_npz(paths[2]).tocsr()
 assert matrix.shape==(len(meta),len(genes)) and np.issubdtype(matrix.dtype,np.integer)
 geneids=genes.gene_id.astype(str).str.split('.').str[0]
 assert geneids.is_unique and meta.metacell_id.is_unique
 for r in g.itertuples():
  ids=json.loads(r.qualified_MC_ids)
  indices=pd.Index(meta.metacell_id).get_indexer(ids)
  assert (indices>=0).all() and len(indices)==r.n_mc==len(set(ids))
  col=int(pd.Index(geneids).get_loc(r.gene_id))
  block=matrix[indices]
  gene_count_sum=int(block[:,col].sum())
  all_gene_GEX_sum=int(block.sum())
  cells=int(meta.n_cells.iloc[indices].sum())
  assert cells==r.n_cells and all_gene_GEX_sum>0
  out['results'].append({'study_id':sid,'group_id':r.group_id,'gene_id':r.gene_id,
     'n_mc':r.n_mc,'n_cells':cells,'gene_column':col,
     'gene_count_sum':gene_count_sum,'all_gene_GEX_sum':all_gene_GEX_sum,
     'GEX_CPM':gene_count_sum/all_gene_GEX_sum*1e6})
 del matrix
print('SAFE_GEX_SPOTCHECK='+json.dumps(out),flush=True)
'''.replace('PAYLOAD',encoded)
command="bash -ic 'PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/apasim_analysis/bin/python -'"
r=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker',command],input=remote.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
safe=[x.split('SAFE_GEX_SPOTCHECK=',1)[1] for x in r.stdout.decode(errors='replace').splitlines() if x.startswith('SAFE_GEX_SPOTCHECK=')]
if r.returncode!=0 or not safe:
    (HERE/'source_GEX_spotcheck_status.json').write_text(json.dumps({'status':'FAILED','exit_code':r.returncode,'stderr_present':bool(r.stderr),'raw_startup_text_saved':False},indent=2)+'\n')
    raise RuntimeError('Source GEX spotcheck failed; no shell initialization values persisted')
out=json.loads(safe[-1]);assert len(out['results'])==15
(HERE/'source_GEX_spotcheck.json').write_text(json.dumps(out,indent=2)+'\n')
print('SOURCE_GEX_SPOTCHECK_PASS',len(out['results']),flush=True)
