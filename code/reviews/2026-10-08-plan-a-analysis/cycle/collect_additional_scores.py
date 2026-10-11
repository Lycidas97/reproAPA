"""Read-only remote marker extraction for the other three frozen studies."""
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CONFIG = ROOT / 'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'
ARCHIVE = ROOT / 'notes/2026-09-26-state-comparability-audit/expression_inputs.npz'

REMOTE = r'''
import gzip,hashlib,json,re,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
params = __PARAMS__
arrays,provenance={},{}
def register(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for x in iter(lambda:f.read(2**20),b''):h.update(x)
    provenance[str(path)]=h.hexdigest()
    return path
for sid in ['GSE150311','GSE221739','GSE171127']:
    base=Path(params['frozen_gex'])/sid
    meta_path=Path(params['native_run'])/'studies'/sid/'selected_metacells.tsv.gz'
    genes=pd.read_csv(register(base/'genes.tsv.gz'),sep='\t')
    gmeta=pd.read_csv(register(base/'metacells.tsv.gz'),sep='\t')
    meta=pd.read_csv(register(meta_path),sep='\t')
    assert meta.metacell_id.is_unique and gmeta.metacell_id.is_unique
    ix=pd.Index(gmeta.metacell_id).get_indexer(meta.metacell_id)
    assert (ix>=0).all()
    matrix=sparse.load_npz(register(base/'metacell_by_gene_counts.npz')).tocsr()
    assert matrix.shape==(len(gmeta),len(genes))
    total=gmeta.total_gex.to_numpy()[ix]
    assert (total>0).all()
    sym=genes.gene_symbol.fillna('').str.upper()
    if sid=='GSE171127':
        # Follow Study.gex_symbols: ID placeholders use the pinned GEX GTF.
        annotation=params['mouse_annotation_gtf']
        gtf=register(Path(annotation['path']))
        assert provenance[str(gtf)]==annotation['sha256']
        names={}
        with gzip.open(gtf,'rt') as handle:
            for line in handle:
                if line.startswith('#'):continue
                fields=line.rstrip().split('\t')
                if fields[2]!='gene':continue
                attr=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
                if 'gene_name' not in attr:continue
                gid=attr['gene_id'].split('.')[0]
                assert gid not in names or names[gid]==attr['gene_name']
                names[gid]=attr['gene_name']
        mapped=genes.gene_id.str.split('.').str[0].map(names)
        placeholder=genes.gene_symbol.eq(genes.gene_id)
        sym.loc[placeholder]=mapped.loc[placeholder].str.upper()
    marker_details={}
    for program,markers in params['programs'].items():
        ii=[]
        for name in markers.split():
            found=np.flatnonzero(sym.eq(name))
            assert len(found)==1,(sid,name,'marker missing or duplicated',len(found))
            ii.append(int(found[0]))
        counts=matrix[ix][:,ii].toarray()
        scores=np.log1p(counts/total[:,None]*1e4).mean(1)
        assert np.isfinite(scores).all()
        arrays[sid+'_'+program]=scores
        marker_details[program]=dict(markers=markers,gene_ids=genes.gene_id.iloc[ii].tolist())
    arrays[sid+'_meta_json']=np.array(meta.to_json(orient='records'))
    print(sid,len(meta),json.dumps(marker_details),file=sys.stderr,flush=True)
    del matrix
arrays['provenance_json']=np.array(json.dumps(provenance))
np.savez_compressed(sys.stdout.buffer,**arrays)
'''


def main():
    config = json.loads(CONFIG.read_text())
    with np.load(ARCHIVE, allow_pickle=False) as a:
        programs = json.loads(str(a['human_program_json']))
    params = {key: config[key] for key in ['frozen_gex', 'native_run']}
    params['programs'] = {key: programs[key] for key in ['Cycle_S', 'Cycle_G2M']}
    params['mouse_annotation_gtf'] = config['resources']['mouse_annotation_gtf']
    script = REMOTE.replace('__PARAMS__', repr(params))
    (HERE/'additional_scores_remote.py').write_text(script)
    partial = HERE/'additional_cycle_scores.partial.npz'
    with partial.open('wb') as output, (HERE/'additional_scores_collection.log').open('w') as log:
        subprocess.run(['ssh', '-o', 'BatchMode=yes', '213docker', 'python3 -'],
                       input=script.encode(), stdout=output, stderr=log, check=True)
    with np.load(partial, allow_pickle=False) as check:
        assert all(sid+'_'+key in check.files for sid in ['GSE150311','GSE221739','GSE171127']
                   for key in ['meta_json','Cycle_S','Cycle_G2M'])
    partial.replace(HERE/'additional_cycle_scores.npz')
    provenance = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [CONFIG, ARCHIVE, Path(__file__)]}
    provenance['params'] = params
    (HERE/'additional_scores_extraction_manifest.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (HERE/'additional_cycle_scores.complete.json').write_text(json.dumps(dict(
        status='complete', studies=['GSE150311','GSE221739','GSE171127'],
        sha256=hashlib.sha256((HERE/'additional_cycle_scores.npz').read_bytes()).hexdigest()), indent=2)+'\n')
    print('Additional score archive:', (HERE/'additional_cycle_scores.npz').stat().st_size, 'bytes')


if __name__ == '__main__':
    main()
