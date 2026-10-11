#!/usr/bin/env python3
"""Materialize and independently check the frozen raw-count reference."""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import yaml

from prepare_inputs import decode


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    config = yaml.safe_load(Path(p.parse_args().config).read_text())
    root = Path(config['output_root'])
    prep = root / 'prepared'
    assert json.loads((prep / 'preflight.json').read_text())['status'] == 'PASS'
    target = root / 'reference'
    target.mkdir(parents=True, exist_ok=True)
    if (target / 'validation.json').exists():
        raise FileExistsError('Reference already validated; inspect before rebuilding')
    cells = pd.read_csv(prep / 'reference_cells.tsv.gz', sep='\t', dtype=str).fillna('')
    libraries = pd.read_csv(prep / 'libraries.tsv', sep='\t', dtype=str).set_index('library_key')
    features = pd.read_csv(prep / 'shared_features.tsv', sep='\t', dtype=str).fillna('')
    chunks, parts, totals, checks = [], [], [], []
    for key, group in cells.groupby('library_key', sort=True):
        path = libraries.loc[key, 'input_path']
        with h5py.File(path, 'r') as handle:
            raw = handle['matrix']
            barcodes = decode(raw['barcodes'][:])
            genes = decode(raw['features/id'][:])
            row_index = pd.Index(barcodes).get_indexer(group.barcode)
            col_index = pd.Index(genes).get_indexer(features.gene_id)
            assert (row_index >= 0).all() and (col_index >= 0).all()
            matrix = sparse.csc_matrix((raw['data'][:], raw['indices'][:], raw['indptr'][:]),
                                       shape=tuple(raw['shape'][:])).T.tocsr()
            chosen = matrix[row_index][:, col_index].tocsr()
            assert np.isfinite(chosen.data).all() and (chosen.data >= 0).all()
            assert np.equal(chosen.data, np.floor(chosen.data)).all()
            full_totals = np.asarray(matrix[row_index].sum(axis=1)).ravel().astype(np.int64)
            shared_totals = np.asarray(chosen.sum(axis=1)).ravel().astype(np.int64)
            assert (shared_totals > 0).all()
            chunks.append(chosen.astype(np.int32))
            parts.append(group)
            totals.extend(shared_totals.tolist())
            checks.append({'library_key': key, 'cells': len(group), 'shared_gene_umis': int(shared_totals.sum()),
                           'all_input_feature_umis': int(full_totals.sum()), 'shared_umi_fraction': float(shared_totals.sum() / full_totals.sum())})
        print(f'{key}: {len(group)} cells extracted', flush=True)
    obs = pd.concat(parts, ignore_index=True).set_index('cell_id', drop=False)
    matrix = sparse.vstack(chunks, format='csr')
    a = ad.AnnData(X=matrix, obs=obs, var=features.set_index('gene_id', drop=False))
    a.uns['expression_scale'] = 'raw_integer_counts'
    a.uns['feature_contract'] = 'exact_ensembl_id_intersection_without_zero_imputation'
    a.uns['used_for_model_training'] = False
    a.write_h5ad(target / 'human_reference_counts.h5ad', compression='gzip')
    # Reopen the saved artifact: ID order, shape and per-cell count sums must close.
    reread = ad.read_h5ad(target / 'human_reference_counts.h5ad')
    assert reread.obs_names.equals(a.obs_names) and reread.var_names.equals(a.var_names)
    assert reread.shape == (len(cells), len(features))
    assert set(reread.obs_names) == set(cells.cell_id)
    assert np.array_equal(np.asarray(reread.X.sum(axis=1)).ravel(), np.asarray(totals))
    assert reread.obs.study_id.nunique() == 4 and reread.obs.library_key.nunique() == 59
    pd.DataFrame(checks).to_csv(target / 'library_count_audit.tsv', sep='\t', index=False)
    report = {'status': 'PASS', 'reference_cells': a.n_obs, 'shared_genes': a.n_vars,
              'libraries': 59, 'all_saved_cell_ids_match': True, 'all_saved_gene_ids_match': True,
              'per_cell_counts_conserved': True, 'raw_integer_counts': True,
              'minimum_library_shared_umi_fraction': min(x['shared_umi_fraction'] for x in checks),
              'reference_h5ad_sha256': hashlib.sha256((target / 'human_reference_counts.h5ad').read_bytes()).hexdigest(),
              'joint_model_trained': False}
    (target / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
