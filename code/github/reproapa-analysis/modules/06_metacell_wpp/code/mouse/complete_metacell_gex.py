import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def aggregate(matrix, cells):
    if cells.barcode.duplicated().any() or cells.cell_id.duplicated().any():
        raise ValueError('Duplicate cell identity')
    if not pd.Index(cells.barcode).isin(matrix.obs_names).all():
        raise ValueError('Missing exact GEX barcodes')
    x = matrix[cells.barcode].X.tocsr()
    if not np.isfinite(x.data).all() or (x.data < 0).any() or not np.equal(x.data, np.floor(x.data)).all():
        raise ValueError('GEX input is not nonnegative integer counts')
    x = x.astype(np.int64)
    ids = pd.Index(sorted(cells.metacell_id.unique()), name='metacell_id')
    codes = ids.get_indexer(cells.metacell_id)
    g = sparse.csr_matrix((np.ones(len(cells), dtype=np.int64), (codes, np.arange(len(cells)))), shape=(len(ids),len(cells)))
    y = g @ x
    if not np.array_equal(np.asarray(y.sum(0)), np.asarray(x.sum(0))):
        raise ValueError('Gene count conservation failed')
    obs = pd.DataFrame({'n_cells':np.asarray(g.sum(1)).ravel()}, index=ids)
    obs['gex_umi'] = np.asarray(y.sum(1)).ravel()
    obs['detected_genes'] = np.asarray((y > 0).sum(1)).ravel()
    return ad.AnnData(X=y, obs=obs, var=matrix.var.copy())


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--assignments', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    m = pd.read_csv(a.manifest, sep='\t').set_index('library_key')
    record = m.loc[a.library]
    cells = pd.read_csv(a.assignments, sep='\t')
    if len(cells) != int(record.n_cells) or set(cells.library_key) != {a.library}:
        raise ValueError('Assignment/manifest identity mismatch')
    path = Path(record.matrix_path)
    matrix = sc.read_10x_h5(path) if path.suffix == '.h5' else sc.read_10x_mtx(path)
    matrix.var_names_make_unique()
    result = aggregate(matrix, cells)
    result.obs['library_key'] = a.library
    a.output.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(a.output/'all_gene_counts.h5ad', compression='lzf')
    result.obs.to_csv(a.output/'gex_coverage.tsv', sep='\t')
    rows = []
    for field in ['local_state_id','trajectory_phase','fate_axis','biology_label','confidence']:
        frame = cells[['metacell_id',field]].fillna({field:'unresolved'})
        grouped = frame.groupby(['metacell_id',field], observed=True).size().rename('n_cells').reset_index().rename(columns={field:'label'})
        grouped['fraction'] = grouped.n_cells / grouped.metacell_id.map(result.obs.n_cells)
        grouped['field'] = field
        rows.append(grouped)
    composition = pd.concat(rows, ignore_index=True)
    if not np.allclose(composition.groupby(['metacell_id','field']).fraction.sum(), 1):
        raise ValueError('Composition fractions do not sum to one')
    composition.to_csv(a.output/'state_composition.tsv.gz', sep='\t', index=False)
    summary = dict(library_key=a.library, cells=len(cells), metacells=result.n_obs,
                   genes=result.n_vars, umi=int(result.X.sum()), exact_barcode_match=True,
                   gene_counts_conserved=True, composition_closed=True)
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (a.output/'complete.ok').write_text('GEX_COMPOSITION_COMPLETE\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
