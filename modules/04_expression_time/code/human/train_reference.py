#!/usr/bin/env python3
"""Fit the mouse workflow's scVI settings on the frozen human raw-count reference."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scvi
import torch
import yaml


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    config = yaml.safe_load(Path(p.parse_args().config).read_text())
    root = Path(config['output_root'])
    spec = config['model']
    assert json.loads((root / 'reference/validation.json').read_text())['status'] == 'PASS'
    target = root / 'model'
    target.mkdir(parents=True, exist_ok=True)
    if (target / 'scvi_model').exists():
        raise FileExistsError('Refusing to overwrite a fitted model')
    assert torch.cuda.is_available(), 'A CUDA GPU is required by this model configuration'
    scvi.settings.seed = config['seed']
    torch.set_float32_matmul_precision('high')
    a = ad.read_h5ad(root / 'reference/human_reference_counts.h5ad')
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=spec['n_hvg'], flavor=spec['hvg_flavor'],
                                layer='counts', batch_key=spec['hvg_batch_key'], subset=True)
    scvi.model.SCVI.setup_anndata(a, layer='counts', batch_key=spec['batch_key'])
    model = scvi.model.SCVI(a, n_latent=spec['n_latent'], gene_likelihood=spec['gene_likelihood'])
    model.train(max_epochs=spec['max_epochs'], early_stopping=spec['early_stopping'],
                accelerator='gpu', devices=1, batch_size=spec['batch_size'])
    a.obsm['X_scVI'] = model.get_latent_representation()
    assert np.isfinite(a.obsm['X_scVI']).all()
    model.save(target / 'scvi_model', overwrite=False)
    a.write_h5ad(target / 'human_scvi_reference.h5ad', compression='gzip')
    a.var[['gene_id', 'gene_symbol']].to_csv(target / 'hvg_features.tsv', sep='\t', index=False)
    report = {'status': 'MODEL_FIT_REQUIRES_BIOLOGICAL_REVIEW', 'cells': a.n_obs, 'genes': a.n_vars,
              'latent_dimensions': spec['n_latent'], 'batch_key': spec['batch_key'],
              'time_used_in_scvi': False, 'clusters_used_in_scvi': False, 'sensitivity_used_to_fit': False,
              'full_primary_projection_complete': False, 'joint_layout_complete': False}
    (target / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
