#!/usr/bin/env python3
"""Project one library's retained cells through the fitted, unchanged encoder."""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import scvi
import torch
import yaml

from prepare_inputs import decode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--library-key', required=True)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    root = Path(config['output_root'])
    prep, model_root = root / 'prepared', root / 'model'
    target = root / 'projection/libraries'
    target.mkdir(parents=True, exist_ok=True)
    output = target / (args.library_key + '.npz')
    report_path = target / (args.library_key + '.json')
    if report_path.exists():
        raise FileExistsError('Refusing to overwrite a validated library projection')
    torch.set_num_threads(args.threads)
    manifest = pd.read_csv(prep / 'libraries.tsv', sep='\t', dtype=str).set_index('library_key')
    library = manifest.loc[args.library_key]
    columns = ['cell_id', 'barcode', 'study_id', 'library_key', 'library_id', 'qc_tier']
    parts = [pd.read_csv(prep / (tier + '_cells.tsv.gz'), sep='\t', dtype=str, usecols=columns)
             for tier in ['primary', 'sensitivity']]
    cells = pd.concat(parts, ignore_index=True)
    cells = cells[cells.library_key.eq(args.library_key)].sort_values('cell_id').set_index('cell_id', drop=False)
    assert len(cells) and cells.index.is_unique
    source = Path(library.input_path)
    provenance = json.loads((prep / 'source_provenance.json').read_text())
    frozen_source = next(row for row in provenance['raw_matrix_stats'] if row['path'] == str(source))
    assert source.stat().st_size == frozen_source['bytes'] and source.stat().st_mtime_ns == frozen_source['mtime_ns']
    reference = ad.read_h5ad(model_root / 'human_scvi_reference.h5ad')
    assert reference.n_vars == config['model']['n_hvg']
    with h5py.File(source, 'r') as handle:
        raw = handle['matrix']
        rows = pd.Index(decode(raw['barcodes'][:])).get_indexer(cells.barcode)
        cols = pd.Index(decode(raw['features/id'][:])).get_indexer(reference.var_names)
        assert (rows >= 0).all() and (cols >= 0).all()
        matrix = sparse.csc_matrix((raw['data'][:], raw['indices'][:], raw['indptr'][:]),
                                    shape=tuple(raw['shape'][:])).T.tocsr()
        counts = matrix[rows][:, cols].tocsr()
    del matrix
    assert np.isfinite(counts.data).all() and (counts.data >= 0).all()
    assert np.equal(counts.data, np.floor(counts.data)).all()
    assert (np.asarray(counts.sum(axis=1)).ravel() > 0).all()
    query = ad.AnnData(X=counts, obs=cells, var=reference.var.copy())
    query.layers['counts'] = query.X
    assert query.var_names.equals(reference.var_names)
    reference_rows = reference.obs_names.get_indexer(query.obs_names)
    overlap = reference_rows >= 0
    assert (query.layers['counts'][overlap] != reference.layers['counts'][reference_rows[overlap]]).nnz == 0
    model = scvi.model.SCVI.load(model_root / 'scvi_model', adata=reference, accelerator='cpu')
    latent = model.get_latent_representation(query, batch_size=2048).astype(np.float32)
    assert latent.shape == (len(cells), config['model']['n_latent']) and np.isfinite(latent).all()
    expected = reference.obsm['X_scVI'][reference_rows[overlap]]
    max_abs = float(np.abs(latent[overlap] - expected).max())
    assert np.allclose(latent[overlap], expected, atol=0.01, rtol=0.001), max_abs
    latent[overlap] = expected
    np.savez_compressed(output, cell_id=query.obs_names.to_numpy(dtype=str), latent=latent)
    with np.load(output, allow_pickle=False) as saved:
        assert np.array_equal(saved['cell_id'], query.obs_names.to_numpy(dtype=str))
        assert np.array_equal(saved['latent'], latent)
    report = {'status': 'PASS', 'library_key': args.library_key, 'cells': len(cells),
              'qc_tiers': cells.qc_tier.value_counts().to_dict(), 'latent_dimensions': latent.shape[1],
              'reference_overlap': int(overlap.sum()), 'reference_raw_hvg_counts_exactly_equal': True,
              'reference_latent_max_abs_before_replacement': max_abs,
              'reference_comparison_atol': 0.01, 'reference_comparison_rtol': 0.001,
              'saved_reference_latent_exactly_preserved': True,
              'frozen_encoder': True, 'query_finetuned': False,
              'raw_source_matches_preparation_stat': True,
              'model_sha256': hashlib.sha256((model_root / 'scvi_model/model.pt').read_bytes()).hexdigest(),
              'projection_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
