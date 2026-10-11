"""Export unchanged all-gene MC counts on the new run's fixed MC registry."""
import argparse
import json
import os
from pathlib import Path

os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from polyapipe_prepare import digest, save

BASE = Path('/root/nfsdata/apadetect/analysis_projects')
MOUSE = BASE/'mouse_time_aware_integration_pilot_v1/results_mouse_metacell_gex_completion_v1'
HUMAN = BASE/'human_time_aware_integration_v1/results_metacell_v1/release/human_metacell_all_gene_counts.h5ad'


def run(native, out):
    out.mkdir(parents=True, exist_ok=False)
    inputs, summaries = {}, []
    human = ad.read_h5ad(HUMAN, backed='r')
    inputs[str(HUMAN)] = digest(HUMAN)
    try:
        for study in sorted((native/'studies').iterdir()):
            mp=study/'metacells.tsv.gz'
            inputs[str(mp)]=digest(mp)
            meta=pd.read_csv(mp,sep='\t')
            # Species is fixed independently of PAS annotation eligibility.
            species=pd.read_csv(study/'coordinate_features.tsv.gz',sep='\t',nrows=1).species.iloc[0]
            blocks, registry, row_order, genes = [], None, [], None
            for library, ids in meta.groupby('library_key',sort=False).groups.items():
                path=MOUSE/library/'all_gene_counts.h5ad' if species=='mouse' else HUMAN
                a=ad.read_h5ad(path,backed='r') if species=='mouse' else human
                inputs[str(path)]=inputs.get(str(path)) or digest(path)
                gi=pd.Index(a.var.gene_ids.astype(str) if 'gene_ids' in a.var else a.var_names.astype(str))
                normalized=gi.str.replace(r'\.[0-9]+$','',regex=True)
                if registry is None:
                    registry=pd.Index(sorted(set(normalized)))
                    symbol=next((c for c in ['gene_symbols','gene_symbol','symbol','gene_name'] if c in a.var),None)
                    genes=pd.DataFrame(dict(gene_id=normalized,gene_symbol=a.var[symbol].astype(str).to_numpy() if symbol else normalized))
                    genes=genes.groupby('gene_id',sort=True).gene_symbol.agg(lambda x:';'.join(sorted(set(x)))).reindex(registry).reset_index()
                    genes.columns=['gene_id','gene_symbol']
                assert set(registry)==set(normalized), (study.name,library,'GEX feature universes differ')
                rows=a.obs_names.get_indexer(meta.loc[ids,'metacell_id'])
                assert (rows>=0).all() and np.array_equal(a.obs.iloc[rows].n_cells,meta.loc[ids,'n_cells'])
                x=a.X[rows].tocsr()
                assert (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
                x=x.astype(np.int64)
                collapse=sparse.csr_matrix((np.ones(len(gi),dtype=np.int64),(np.arange(len(gi)),registry.get_indexer(normalized))),shape=(len(gi),len(registry)))
                mapped=x@collapse
                assert np.array_equal(np.asarray(x.sum(1)).ravel(),np.asarray(mapped.sum(1)).ravel())
                blocks.append(mapped);row_order.extend(ids)
                if species=='mouse':a.file.close()
            counts=sparse.vstack(blocks,format='csr')[np.argsort(row_order)]
            totals=np.asarray(counts.sum(1)).ravel()
            assert len(totals)==len(meta) and (totals>0).all()
            dest=out/study.name;dest.mkdir()
            sparse.save_npz(dest/'metacell_by_gene_counts.npz',counts)
            save(genes,dest/'genes.tsv.gz')
            save(meta.assign(total_gex=totals),dest/'metacells.tsv.gz')
            summaries.append(dict(study_id=study.name,metacells=len(meta),genes=len(genes),molecules=int(totals.sum())))
            print(json.dumps(summaries[-1]),flush=True)
    finally:
        human.file.close()
    assert all(digest(p)==h for p,h in inputs.items())
    (out/'complete.json').write_text(json.dumps(dict(status='PASS',inputs=inputs,studies=summaries,
        checks=dict(fixed_MC_identity=True,integer_count_conservation=True,all_gene_denominator=True)),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native-run',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.native_run,a.out)
