#!/usr/bin/env python3
"""Check saved expression outputs against frozen upstream identities and QC scope."""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    config = yaml.safe_load(Path(parser.parse_args().config).read_text())
    root = Path(config['output_root'])
    target = root / 'expression_diagnostics'
    summary = json.loads((target / 'summary.json').read_text())
    primary = ad.read_h5ad(target / 'human_primary_scvi.h5ad')
    sensitivity = ad.read_h5ad(target / 'human_sensitivity_scvi.h5ad')
    checks = {}

    def check(name, value):
        checks[name] = bool(value)
        if not value:
            raise ValueError(name)

    check('full_primary_count', primary.n_obs == 215114)
    check('full_sensitivity_count', sensitivity.n_obs == 42878)
    ids = primary.obs_names.append(sensitivity.obs_names)
    check('all_retained_cells_unique', ids.is_unique)
    technical = pd.read_csv(root / 'prepared/technical_exclude_cells.tsv.gz', sep='\t', usecols=['cell_id'])
    check('technical_cells_not_in_expression_objects', not set(ids).intersection(technical.cell_id))
    for tier, obj in [('primary', primary), ('sensitivity', sensitivity)]:
        frozen = pd.read_csv(root / 'prepared' / (tier + '_cells.tsv.gz'), sep='\t', dtype=str).fillna('').set_index('cell_id', drop=False)
        check(tier + '_all_frozen_cell_ids_match', obj.obs_names.is_unique and set(obj.obs_names) == set(frozen.index))
        frame = frozen.loc[obj.obs_names]
        check(tier + '_all_frozen_metadata_fields_unchanged', all(obj.obs[column].astype(str).eq(frame[column]).all() for column in frozen.columns))
        check(tier + '_finite_twenty_dimensional_latent', obj.obsm['X_scVI'].shape == (obj.n_obs, 20) and np.isfinite(obj.obsm['X_scVI']).all())
        check(tier + '_latent_not_mislabeled_as_gene_expression', obj.X is None and obj.n_vars == 0)
    for study, spec in config['studies'].items():
        path = Path(config['clustering_root']) / spec['clustering_release'] / 'frozen/cell_clusters.tsv.gz'
        original = pd.read_csv(path, sep='\t', dtype=str, usecols=['cell_id', 'cluster_id']).set_index('cell_id')
        saved = primary.obs[primary.obs.study_id.eq(study)]
        check(study + '_original_release_labels_unchanged', set(saved.index) == set(original.index) and
              saved.cluster_id.astype(str).eq(original.loc[saved.index, 'cluster_id']).all())
    check('all_59_libraries_and_62_states_present', primary.obs.library_key.nunique() == 59 and primary.obs.original_state_id.nunique() == 62)
    graph = primary.obsp['connectivities'].tocsr()
    check('primary_only_graph_finite_and_nonempty', graph.shape == (215114, 215114) and graph.nnz > 0 and np.isfinite(graph.data).all())
    distances = primary.obsp['distances'].tocsr()
    rows = np.repeat(np.arange(primary.n_obs), np.diff(distances.indptr))
    nonself = rows != distances.indices
    study = primary.obs.study_id.astype(str).to_numpy()
    observed_cross = np.mean(study[rows[nonself]] != study[distances.indices[nonself]])
    neighbor_counts = pd.read_csv(target / 'study_neighbor_counts.tsv', sep='\t', index_col=0)
    check('neighbor_statistics_exclude_all_self_edges', summary['neighbor_metrics_exclude_self'] and
          summary['nonself_directed_neighbor_edges'] == int(nonself.sum()) == int(neighbor_counts.to_numpy().sum()))
    check('cross_study_neighbor_rate_recomputed_from_saved_graph', np.isclose(observed_cross, summary['cross_study_neighbor_fraction'], atol=1e-12, rtol=0))
    check('primary_umap_finite', primary.obsm['X_umap'].shape == (215114, 2) and np.isfinite(primary.obsm['X_umap']).all())
    check('sensitivity_not_fitted_to_primary_graph_or_layout', len(sensitivity.obsp) == 0 and 'X_umap' not in sensitivity.obsm)
    check('time_and_original_labels_not_graph_inputs', not summary['time_used_to_fit_graph'] and not summary['clusters_used_to_fit_graph'])
    model_hash = hashlib.sha256((root / 'model/scvi_model/model.pt').read_bytes()).hexdigest()
    check('both_outputs_use_same_frozen_model', primary.uns['model_sha256'] == sensitivity.uns['model_sha256'] == model_hash)
    check('no_query_finetuning', not primary.uns['query_finetuned'] and not sensitivity.uns['query_finetuned'])
    artifacts = {}
    for name in ['human_primary_scvi.h5ad', 'human_sensitivity_scvi.h5ad', 'primary_coordinates.tsv.gz']:
        path = target / name
        artifacts[name] = {'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    report = {'status': 'TECHNICAL_PASS', 'phase': 'FULL_EXPRESSION_INTEGRATION_FIRST_PASS', 'checks': checks,
              'primary_cells': primary.n_obs, 'sensitivity_cells': sensitivity.n_obs, 'model_sha256': model_hash,
              'artifacts': artifacts, 'human_ot_solved': False, 'biological_review': 'RECORDED_SEPARATELY'}
    (target / 'independent_acceptance.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'checks_passed': sum(checks.values()), 'checks_total': len(checks)}, indent=2))


if __name__ == '__main__':
    main()
