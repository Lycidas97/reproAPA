"""Extract frozen PAS subsets over read-only SSH; no files are written remotely."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
CALLERS = {
    'polyapipe': ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables',
    'apadetect_p50': ROOT/'notes/2026-10-06-dual-main-plot-data/apadetect-p50/supplement/tables',
}
STUDIES = ['GSE103221','GSE115943','GSE118258']

REMOTE = '''
import hashlib,io,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
CONFIG = __CONFIG__
arrays={};provenance={}
for rec in CONFIG:
    def checked(path, expected):
        h=hashlib.sha256()
        with open(path,'rb') as f:
            for chunk in iter(lambda:f.read(2**20),b''): h.update(chunk)
        value=h.hexdigest();assert value==expected,(path,value,expected)
        provenance[path]=value
    checked(rec['counts_path'],rec['counts_sha256'])
    checked(rec['meta_path'],rec['meta_sha256'])
    checked(rec['full_meta_path'],rec['full_meta_sha256'])
    meta=pd.read_csv(rec['meta_path'],sep='\\t')
    full=pd.read_csv(rec['full_meta_path'],sep='\\t')
    assert np.array_equal(full.metacell_id.iloc[meta.matrix_row],meta.metacell_id)
    matrix=sparse.load_npz(rec['counts_path']).tocsr()[meta.matrix_row.to_numpy()][:,rec['columns']].tocsr()
    assert np.issubdtype(matrix.data.dtype,np.integer) and (matrix.data>=0).all()
    key=rec['caller']+'__'+rec['study_id']+'__'
    arrays[key+'meta_json']=np.array(meta.to_json(orient='records'))
    arrays[key+'data']=matrix.data;arrays[key+'indices']=matrix.indices
    arrays[key+'indptr']=matrix.indptr;arrays[key+'shape']=np.array(matrix.shape)
    arrays[key+'columns']=np.array(rec['columns'])
    print(rec['caller'],rec['study_id'],matrix.shape,matrix.nnz,file=sys.stderr,flush=True)
arrays['provenance_json']=np.array(json.dumps(provenance))
np.savez_compressed(sys.stdout.buffer,**arrays)
'''

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--collect',action='store_true');args=parser.parse_args()
    records=[];local={}
    for caller,base in CALLERS.items():
        audit_path=base/'measurement/audit.json';axes_path=base/'state_vectors/PAS_axes.primary.tsv.gz'
        audit=json.loads(audit_path.read_text());axes=pd.read_csv(axes_path,sep='\t',usecols=['study_id','feature_index_1based'])
        for p in [audit_path,axes_path]:local[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
        for sid in STUDIES:
            folder='WOT_GSE115943' if sid=='GSE115943' else sid
            paths={Path(p).name:(p,h) for p,h in audit['inputs'].items() if f'/studies/{folder}/' in p}
            rec=dict(caller=caller,study_id=sid,columns=(axes[axes.study_id.eq(sid)].feature_index_1based-1).tolist())
            for stem,name in [('counts','metacell_by_pas_counts.npz'),('meta','selected_metacells.tsv.gz'),('full_meta','metacells.tsv.gz')]:
                rec[stem+'_path'],rec[stem+'_sha256']=paths[name]
            assert len(rec['columns'])==len(set(rec['columns']))
            records.append(rec)
    (HERE/'extraction_manifest.json').write_text(json.dumps(dict(remote=records,local=local),indent=2)+'\n')
    script=REMOTE.replace('__CONFIG__',repr(records));(HERE/'extract_remote.py').write_text(script)
    if args.collect:
        with (HERE/'pas_subsets.npz').open('wb') as out,(HERE/'collection.log').open('w') as log:
            subprocess.run(['ssh','-o','BatchMode=yes','213docker','python3 -'],input=script.encode(),stdout=out,stderr=log,check=True)
        print('Collected', (HERE/'pas_subsets.npz').stat().st_size,'bytes')

if __name__=='__main__':main()
