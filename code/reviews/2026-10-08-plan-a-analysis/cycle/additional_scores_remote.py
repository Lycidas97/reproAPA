
import gzip,hashlib,json,re,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
params = {'frozen_gex': '/root/nfsdata/apadetect/analysis_projects/reproapa_polyapipe_rerun_20261003_v1/frozen_gex', 'native_run': '/root/nfsdata/apadetect/analysis_projects/reproapa_polyapipe_common_statistics_20261005/run-series-v2-primary/runs/polyapipe_common/statistics/core/native', 'programs': {'Cycle_S': 'MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6', 'Cycle_G2M': 'CDK1 TOP2A NUSAP1 CDC20 CDC45 UBE2C BIRC5 TPX2'}, 'mouse_annotation_gtf': {'path': '/root/nfsdata/apadetect/data/genome/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz', 'sha256': 'dfc17ee4c85f3d7665f48fe758f9a1109c3cc4c38a4185943e9e9678fae45fe6'}}
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
