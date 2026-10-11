"""Frozen mouse aggregation helper, reused for human PAS counts."""

import anndata as ad
import numpy as np
import pandas as pd
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

