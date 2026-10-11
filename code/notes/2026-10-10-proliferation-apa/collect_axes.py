"""Small read-only raw-coordinate snapshot for independent whole-primary axes QA."""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def main():
 x=np.load(HERE/'inputs.npz',allow_pickle=False);panels=pd.DataFrame(json.loads(str(x['panels_json'])))
 request={sid:sorted(set(str(g) for cid in q.context_id for g in x[cid+'_gene_ids'])) for sid,q in panels[panels.status.eq('available')].groupby('study_id')}
 config=json.loads((ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json').read_text());expected=json.loads(str(x['provenance_json']))
 remote=r'''
from pathlib import Path
import gzip,hashlib,io,json,sys
import pandas as pd
pieces=[]
for sid,genes in REQUEST.items():
 folder='WOT_GSE115943' if sid=='GSE115943' else sid
 path=Path(CONFIG['native_run'])/'studies'/folder/'features.tsv.gz'
 with path.open('rb') as h:value=hashlib.file_digest(h,'sha256').hexdigest()
 assert value==EXPECTED[str(path)]
 cols=['primary_gene_id','primary_gene_symbol','keep_for_apa_primary','matched_source_count_50nt','cut_site','strand','chrom','feature_index_1based','atlas_pas_id']
 f=pd.read_csv(path,sep='\t',usecols=cols,low_memory=False)
 f=f[f.keep_for_apa_primary.astype(str).str.lower().eq('true')&f.primary_gene_id.isin(genes)].copy();f['study_id']=sid
 assert set(f.primary_gene_id)==set(genes)
 pieces.append(f)
 print(sid,len(f),value,file=sys.stderr,flush=True)
b=io.BytesIO()
with gzip.GzipFile(fileobj=b,mode='wb',mtime=0) as h:h.write(pd.concat(pieces,ignore_index=True).to_csv(sep='\t',index=False).encode())
sys.stdout.buffer.write(b.getvalue())
'''
 code='CONFIG = '+repr(config)+'\nEXPECTED = '+repr(expected)+'\nREQUEST = '+repr(request)+'\n'+remote
 with (HERE/'raw_primary_axis_features.partial.tsv.gz').open('wb') as out,(HERE/'axes_extract.log').open('wb') as log:
  subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker','OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/scanpy/bin/python -'],input=code.encode(),stdout=out,stderr=log,check=True)
 (HERE/'raw_primary_axis_features.partial.tsv.gz').rename(HERE/'raw_primary_axis_features.tsv.gz');print((HERE/'axes_extract.log').read_text())
if __name__=='__main__':main()
