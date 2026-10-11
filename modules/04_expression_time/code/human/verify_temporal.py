"""Independent input, directed-support and saved-layout audit, plus matched diagnostics."""
import argparse
import hashlib
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
import yaml


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    args = p.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    root, source = Path(cfg['output_root']), Path(cfg['input_root'])
    out = root / 'validation'
    out.mkdir(exist_ok=True)
    checks = []

    def check(name, value, detail=''):
        checks.append(dict(check=name, passed=bool(value), detail=str(detail)))
        assert value, f'{name}: {detail}'

    a = ad.read_h5ad(source / 'expression_diagnostics/human_primary_scvi.h5ad')
    saved = ad.read_h5ad(root / 'layout/human_time_constrained.h5ad')
    primary = pd.read_csv(source / 'prepared/primary_cells.tsv.gz', sep='\t', dtype=str).fillna('')
    sensitivity = pd.read_csv(source / 'prepared/sensitivity_cells.tsv.gz', sep='\t', usecols=['cell_id'])
    check('all_primary_once', saved.n_obs == cfg['expected_cells'] and saved.obs_names.is_unique and saved.obs_names.equals(a.obs_names))
    check('exact_registered_order', list(saved.obs_names) == primary.cell_id.tolist())
    check('no_sensitivity', not set(saved.obs_names).intersection(sensitivity.cell_id))
    check('all_original_metadata_exact', all(saved.obs[c].astype(str).tolist() == a.obs[c].astype(str).tolist() for c in a.obs.columns))
    check('62_states_59_libraries', saved.obs.original_state_id.nunique() == 62 and saved.obs.library_key.nunique() == 59)
    check('latent_unchanged', np.array_equal(saved.obsm['X_scVI'], a.obsm['X_scVI']))
    check('expression_initialization_preserved', np.array_equal(saved.obsm['X_umap_expression'], a.obsm['X_umap']))
    check('model_unchanged', saved.uns['model_sha256'] == a.uns['model_sha256'])
    node_order = pd.read_csv(root / 'expression/node_order.tsv.gz', sep='\t')
    check('expression_kernel_node_order_exact', node_order.cell_id.tolist() == list(a.obs_names))
    expression_report = json.loads((root / 'expression/summary.json').read_text())
    check('expression_kernel_input_hash', expression_report['input_sha256'] == hashlib.sha256(
          (source / 'expression_diagnostics/human_primary_scvi.h5ad').read_bytes()).hexdigest())
    check('sorted_conditions_retained', saved.obs.loc[saved.obs.library_key.astype(str).str.contains('D16_'), 'condition'].nunique() == 2)
    plan_path = source / 'prepared/temporal_edge_plan.tsv'
    plan = pd.read_csv(plan_path, sep='\t').fillna('')
    check('34_registered_pairs', len(plan) == cfg['expected_pairs'])
    lib_codes, libraries = pd.factorize(a.obs.library_key, sort=True)
    allowed = np.zeros((len(libraries), len(libraries)), dtype=bool)
    time_cells, source_cells = set(), set()
    for i, row in plan.iterrows():
        directory = root / 'pairs' / f'{i:02d}'
        report = json.loads((directory / 'summary.json').read_text())
        path = directory / 'coupling.npz'
        check(f'pair_{i:02d}_converged_and_provenance', report['status'] == 'PASS' and report['converged']
              and report['plan_sha256'] == hashlib.sha256(plan_path.read_bytes()).hexdigest()
              and report['coupling_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
              and report['config_sha256'] == hashlib.sha256(Path(args.config).read_bytes()).hexdigest())
        with np.load(path) as b:
            s, t = [a.obs_names.get_indexer(b[key]) for key in ['source', 'target']]
            check(f'pair_{i:02d}_all_registered_cells', all(np.array_equal(ix, np.flatnonzero(a.obs.library_key.isin(json.loads(row[key]))))
                  for ix, key in [(s, 'source_libraries'), (t, 'target_libraries')]))
            transport = b['transport']
            check(f'pair_{i:02d}_valid_mass', transport.shape == (len(s), len(t)) and np.isfinite(transport).all()
                  and (transport >= 0).all() and (transport.sum(0) > 0).all() and (transport.sum(1) > 0).all()
                  and np.allclose(b['source_mass'], 1. / len(s)) and np.allclose(b['target_mass'], 1. / len(t)))
            check(f'pair_{i:02d}_within_study', a.obs.iloc[np.r_[s, t]].study_id.nunique() == 1)
            if row.relation == 'actual_day':
                d0 = pd.to_numeric(a.obs.iloc[s].actual_day).unique()
                d1 = pd.to_numeric(a.obs.iloc[t].actual_day).unique()
                check(f'pair_{i:02d}_actual_day_forward', len(d0) == len(d1) == 1 and d0[0] < d1[0]
                      and np.isclose(d1[0] - d0[0], float(row.elapsed_days)))
            else:
                check(f'pair_{i:02d}_stage_not_global_day', row.elapsed_days == '' and report['solver_keys'] == [0., 1.]
                      and pd.to_numeric(a.obs.iloc[np.r_[s, t]].actual_day, errors='coerce').isna().all())
                # Stage/day order is lexicographic only within this recorded route.
                sf, tf = a.obs.iloc[s], a.obs.iloc[t]
                s_stage = pd.to_numeric(sf.stage_order, errors='coerce').dropna().unique()
                t_stage = pd.to_numeric(tf.stage_order, errors='coerce').dropna().unique()
                if len(s_stage) and len(t_stage):
                    correct = len(s_stage) == len(t_stage) == 1 and s_stage[0] <= t_stage[0]
                    if correct and s_stage[0] == t_stage[0]:
                        sd = pd.to_numeric(sf.within_stage_day, errors='coerce').dropna().unique()
                        td = pd.to_numeric(tf.within_stage_day, errors='coerce').dropna().unique()
                        correct = len(sd) == len(td) == 1 and sd[0] < td[0]
                    check(f'pair_{i:02d}_stage_order_forward', correct)
            allowed[np.ix_(np.unique(lib_codes[s]), np.unique(lib_codes[t]))] = True
            time_cells.update(s); time_cells.update(t); source_cells.update(s)
        print(f'Audited OT {i + 1}/34', flush=True)
    graphs = {name: sparse.load_npz(root / folder / file).tocsr() for name, folder, file in [
        ('local', 'expression', 'local_kernel.npz'), ('cross', 'expression', 'cross_kernel.npz'),
        ('forward', 'layout', 'forward_temporal_kernel.npz'), ('temporal', 'layout', 'temporal_kernel.npz')]}
    study = a.obs.study_id.astype(str).to_numpy()
    for name, graph in graphs.items():
        graph.eliminate_zeros()
        rr = np.repeat(np.arange(a.n_obs), np.diff(graph.indptr)); cc = graph.indices
        totals = np.asarray(graph.sum(1)).ravel()
        check(name + '_valid_kernel', graph.shape == (a.n_obs, a.n_obs) and np.isfinite(graph.data).all()
              and (graph.data >= 0).all() and np.allclose(totals[totals > 0], 1) and not (rr == cc).any())
        if name == 'local':
            check('local_same_study', (study[rr] == study[cc]).all())
        if name == 'cross':
            check('cross_different_study', (study[rr] != study[cc]).all())
        if name in ['forward', 'temporal']:
            permit = allowed if name == 'forward' else allowed | allowed.T
            check(name + '_only_registered_library_relations', permit[lib_codes[rr], lib_codes[cc]].all())
            check(name + '_no_mixed_time_reference', not (study[rr] == 'GSE150311').any())
            expected = source_cells if name == 'forward' else time_cells
            check(name + '_full_eligible_coverage', set(np.flatnonzero(totals > 0)) == expected)
    shared = np.flatnonzero(a.obs.sample_accession.astype(str).eq('GSM5387808'))
    check('shared_start_once', len(shared) == 3663 and a.obs_names[shared].is_unique)
    for _, row in plan[plan.route_id.str.startswith('gse178325_0330')].iterrows():
        targets = np.flatnonzero(a.obs.library_key.isin(json.loads(row.target_libraries)))
        check('shared_start_route_weight_' + row.route_id, np.allclose(np.asarray(graphs['forward'][shared][:, targets].sum(1)).ravel(), .25))
    flows = pd.read_csv(root / 'layout/state_condition_flows.tsv.gz', sep='\t')
    check('flows_all_34_pairs', set(flows.pair) == set(range(34)))
    check('flows_capture_mass_conserved', np.allclose(flows.groupby('pair').mass_fraction.sum(), 1))
    check('flows_source_conditionals', np.allclose(flows.groupby(['pair', 'source_state_condition']).source_fraction.sum(), 1))
    check('flows_target_conditionals', np.allclose(flows.groupby(['pair', 'target_state_condition']).target_fraction.sum(), 1))
    check('d16_flow_strata_separate', flows.loc[flows.pair.eq(3), 'target_state_condition'].str.contains('TRA-1-60_positive').any()
          and flows.loc[flows.pair.eq(3), 'target_state_condition'].str.contains('TRA-1-60_negative').any())
    temporal = graphs['temporal']
    edge = temporal.tocoo()
    metrics, study_metrics, time_metrics = [], [], []
    studies = sorted(set(study))
    state = a.obs.original_state_id.astype(str).to_numpy()
    for name, xy in [('expression', a.obsm['X_umap']), ('control', saved.obsm['X_umap_control']), ('time', saved.obsm['X_umap'])]:
        check(name + '_finite_coordinates', xy.shape == (a.n_obs, 2) and np.isfinite(xy).all())
        d, j = cKDTree(xy).query(xy, k=31, workers=8)
        keep = j != np.arange(a.n_obs)[:, None]
        pick = np.argsort(~keep, axis=1, kind='stable')[:, :30]
        j = np.take_along_axis(j, pick, axis=1); d = np.take_along_axis(d, pick, axis=1)
        check(name + '_no_self_in_metrics', not (j == np.arange(a.n_obs)[:, None]).any())
        cross = (study[:, None] != study[j]).mean(1)
        same_study = study[:, None] == study[j]
        purity = ((state[:, None] == state[j]) & same_study).sum(1) / np.maximum(same_study.sum(1), 1)
        radius = np.maximum(d[:, -1], 1e-6)
        stretch = np.linalg.norm(xy[edge.row] - xy[edge.col], axis=1) / np.sqrt(radius[edge.row] * radius[edge.col])
        temporal_distance = float(np.average(stretch, weights=edge.data))
        metrics.append(dict(layout=name, knn_2d_cross_study_fraction=float(cross.mean()),
                            mean_same_state_fraction_within_study=float(purity.mean()),
                            mean_time_edge_distance_over_local_radius=temporal_distance,
                            temporal_mass_within_2_local_radii=float(np.average(stretch <= 2, weights=edge.data))))
        for s in studies:
            mask = study == s
            study_metrics.append(dict(layout=name, study_id=s, cells=int(mask.sum()),
                                      knn_2d_cross_study_fraction=float(cross[mask].mean()),
                                      same_state_fraction_within_study=float(purity[mask].mean())))
        for s in ['GSE118258', 'GSE221739']:
            mask = study == s
            ix = np.flatnonzero(mask)
            _, near = cKDTree(xy[ix]).query(xy[ix], k=31, workers=8)
            keep = near != np.arange(len(ix))[:, None]
            near = np.take_along_axis(near, np.argsort(~keep, axis=1, kind='stable')[:, :30], axis=1)
            days = pd.to_numeric(a.obs.iloc[ix].actual_day).to_numpy()
            time_metrics.append(dict(layout=name, study_id=s, mean_neighbor_day_gap=float(np.abs(days[:, None] - days[near]).mean())))
        if name != 'expression':
            graph = sparse.load_npz(root / 'layout' / f'{name}_graph.npz')
            check(name + '_symmetric_graph', (graph - graph.T).nnz == 0 and not graph.diagonal().any())
            weights = cfg['time_weights'] if name == 'time' else cfg['control_weights']
            expected = weights[0] * graphs['local'] + weights[1] * temporal + weights[2] * graphs['cross']
            totals = np.asarray(expected.sum(1)).ravel()
            expected = sparse.diags(1. / np.maximum(totals, 1e-30)) @ expected
            expected = ((expected + expected.T) * .5).tocsr()
            expected.data /= expected.data.max()
            difference = (graph - expected).data
            check(name + '_declared_kernel_weights', not len(difference) or np.max(np.abs(difference)) < 1e-6)
            if name == 'time':
                check('h5ad_graph_exact', (graph != saved.obsp['connectivities']).nnz == 0)
            count, labels = connected_components(graph, directed=False)
            metrics[-1]['graph_components'] = int(count)
            metrics[-1]['largest_component_fraction'] = float(np.bincount(labels).max() / a.n_obs)
    pd.DataFrame(checks).to_csv(out / 'checks.tsv', sep='\t', index=False)
    pd.DataFrame(metrics).to_csv(out / 'layout_comparison.tsv', sep='\t', index=False)
    pd.DataFrame(study_metrics).to_csv(out / 'study_layout_comparison.tsv', sep='\t', index=False)
    pd.DataFrame(time_metrics).to_csv(out / 'actual_day_neighborhoods.tsv', sep='\t', index=False)
    result = dict(status='PASS', checks=len(checks), passed=sum(c['passed'] for c in checks), pairs=len(plan),
                  cells=a.n_obs, time_eligible_cells=len(time_cells), expression_only_cells=a.n_obs - len(time_cells),
                  layouts=metrics, interpretation='Technical acceptance and descriptive diagnostics; no independent lineage truth',
                  h5ad_sha256=hashlib.sha256((root / 'layout/human_time_constrained.h5ad').read_bytes()).hexdigest())
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
