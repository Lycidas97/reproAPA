"""Full-cell conditional OT sampling and one matched time/no-time layout pair."""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
from temporal_graph import log, normalize


def sample(matrix, rng, draws):
    import numpy as np
    from scipy import sparse
    chosen = np.empty((len(matrix), draws), np.int32)
    for start in range(0, len(matrix), 256):
        block = matrix[start:start + 256]
        assert np.isfinite(block).all() and (block >= 0).all() and (block.sum(1) > 0).all()
        cdf = np.cumsum(block / block.sum(1, keepdims=True), axis=1)
        cdf[:, -1] = 1.
        for i in range(len(block)):
            chosen[start + i] = np.searchsorted(cdf[i], rng.random(draws))
    result = sparse.csr_matrix((np.full(chosen.size, 1. / draws, np.float32),
        (np.repeat(np.arange(len(matrix)), draws), chosen.ravel())), shape=matrix.shape)
    assert np.allclose(np.asarray(result.sum(1)).ravel(), 1.)
    return result.tocoo()


def state_flows(t, src_obs, tgt_obs, pair):
    import numpy as np
    import pandas as pd
    from scipy import sparse
    groups = []
    for obs in [src_obs, tgt_obs]:
        labels = obs.original_state_id.astype(str) + '|' + obs.condition.astype(str)
        codes, names = pd.factorize(labels, sort=True)
        groups.append((codes, names, np.bincount(codes)))
    s, sn, sc = groups[0]; d, dn, dc = groups[1]
    membership = sparse.csr_matrix((np.ones(len(s)), (s, np.arange(len(s)))), shape=(len(sn), len(s)))
    by_src = membership @ t
    flow = np.stack([by_src[:, d == j].sum(1) for j in range(len(dn))], axis=1)
    assert np.isclose(flow.sum(), t.sum())
    rows = []
    for i in range(len(sn)):
        for j in range(len(dn)):
            rows.append(dict(pair=pair, source_state_condition=sn[i], target_state_condition=dn[j],
                             source_cells=int(sc[i]), target_cells=int(dc[j]),
                             mass_fraction=float(flow[i, j] / flow.sum()),
                             source_fraction=float(flow[i, j] / flow[i].sum()),
                             target_fraction=float(flow[i, j] / flow[:, j].sum())))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    p.add_argument('--threads', type=int, default=32)
    args = p.parse_args()
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:args.threads])
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components
    from umap.umap_ import find_ab_params, simplicial_set_embedding
    import yaml

    cfg = yaml.safe_load(Path(args.config).read_text())
    root = Path(cfg['output_root'])
    out = root / 'layout'
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'summary.json').exists():
        raise FileExistsError('Refusing to overwrite completed temporal layouts')
    a = ad.read_h5ad(Path(cfg['input_root']) / 'expression_diagnostics/human_primary_scvi.h5ad')
    ids, n = a.obs_names, a.n_obs
    assert n == cfg['expected_cells'] and ids.is_unique and a.obs.original_state_id.nunique() == cfg['expected_states']
    plan = pd.read_csv(Path(cfg['input_root']) / 'prepared/temporal_edge_plan.tsv', sep='\t').fillna('')
    local, cross = [sparse.load_npz(root / 'expression' / (name + '_kernel.npz')) for name in ['local', 'cross']]
    scale = np.load(root / 'expression/scaling.npz')
    z = (a.obsm['X_scVI'] - scale['mean']) / scale['scale']
    rng = np.random.default_rng(cfg['seed'])
    rr, cc, vv, fr, fc, fv, summaries, flows = [], [], [], [], [], [], [], []
    for i, row in plan.iterrows():
        directory = root / 'pairs' / f'{i:02d}'
        report = json.loads((directory / 'summary.json').read_text())
        assert report['converged'] and report['status'] == 'PASS'
        assert report['route_id'] == row.route_id and report['time_from'] == row.time_from and report['time_to'] == row.time_to
        bundle = np.load(directory / 'coupling.npz')
        src, tgt = ids.get_indexer(bundle['source']), ids.get_indexer(bundle['target'])
        assert (src >= 0).all() and (tgt >= 0).all()
        for ix, key in [(src, 'source_libraries'), (tgt, 'target_libraries')]:
            expected = np.flatnonzero(a.obs.library_key.isin(json.loads(row[key])))
            assert np.array_equal(ix, expected), 'Full-cell pair membership differs from frozen plan'
        t = bundle['transport']
        assert t.shape == (len(src), len(tgt))
        flows.extend(state_flows(t, a.obs.iloc[src], a.obs.iloc[tgt], i))
        for direction, (aa, bb, matrix) in enumerate([(src, tgt, t), (tgt, src, t.T)]):
            part = sample(matrix, rng, cfg['draws'])
            r, c = aa[part.row], bb[part.col]
            rr.append(r); cc.append(c); vv.append(part.data)
            if direction == 0:
                fr.append(r); fc.append(c); fv.append(part.data)
        summaries.append(report)
        bundle.close()
        del t
        log(f'Time kernel: {i + 1}/{len(plan)} full OT pairs sampled')
    temporal = normalize(sparse.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(n, n))).astype(np.float32)
    forward = normalize(sparse.csr_matrix((np.concatenate(fv), (np.concatenate(fr), np.concatenate(fc))), shape=(n, n))).astype(np.float32)
    sparse.save_npz(out / 'temporal_kernel.npz', temporal)
    sparse.save_npz(out / 'forward_temporal_kernel.npz', forward)
    pd.DataFrame(summaries).to_csv(out / 'pair_metrics.tsv', sep='\t', index=False)
    pd.DataFrame(flows).merge(plan.reset_index(names='pair')[['pair', 'route_id', 'time_from', 'time_to', 'relation']], on='pair').to_csv(out / 'state_condition_flows.tsv.gz', sep='\t', index=False)
    del rr, cc, vv, fr, fc, fv, forward
    gc.collect()
    initial = a.obsm['X_umap'].copy()
    output = ad.AnnData(obs=a.obs.copy())
    output.obsm['X_scVI'] = a.obsm['X_scVI'].copy()
    output.obsm['X_umap_expression'] = initial
    results = []
    for name, weights in [('time', cfg['time_weights']), ('control', cfg['control_weights'])]:
        mixture = normalize(weights[0] * local + weights[1] * temporal + weights[2] * cross)
        graph = ((mixture + mixture.T) * .5).tocsr()
        graph.data /= graph.data.max()
        graph.eliminate_zeros()
        assert np.isfinite(graph.data).all() and (graph - graph.T).nnz == 0 and not graph.diagonal().any()
        sparse.save_npz(out / f'{name}_graph.npz', graph)
        components = int(connected_components(graph, directed=False)[0])
        log(f'Optimizing {name}: {n} free cells, {graph.nnz // 2} edges, {components} components')
        aa, bb = find_ab_params(cfg['spread'], cfg['min_dist'])
        xy, _ = simplicial_set_embedding(z, graph.copy(), n_components=2, initial_alpha=1., a=aa, b=bb,
            gamma=1., negative_sample_rate=5, n_epochs=cfg['epochs'], init=initial.copy(),
            random_state=np.random.RandomState(cfg['seed']), metric='euclidean', metric_kwds={},
            densmap=False, densmap_kwds={}, output_dens=False, parallel=True, verbose=True)
        assert xy.shape == (n, 2) and np.isfinite(xy).all()
        output.obsm['X_umap' if name == 'time' else 'X_umap_control'] = xy
        np.save(out / f'{name}_coordinates.npy', xy)
        coords = a.obs.copy()
        coords[['UMAP1', 'UMAP2']] = xy
        coords.to_csv(out / f'{name}_coordinates.tsv.gz', sep='\t', index=False)
        if name == 'time':
            output.obsp['connectivities'] = graph
        results.append(dict(name=name, weights=weights, components=components, undirected_edges=graph.nnz // 2,
                            coordinates_sha256=hashlib.sha256((out / f'{name}_coordinates.npy').read_bytes()).hexdigest()))
        log(f'{name} layout saved')
    output.uns['model_sha256'] = a.uns['model_sha256']
    output.uns['time_constraints'] = 'Forward OT on registered within-study adjacent times; symmetrized for soft layout constraints'
    output.uns['temporal_config_sha256'] = hashlib.sha256(Path(args.config).read_bytes()).hexdigest()
    output.write_h5ad(out / 'human_time_constrained.h5ad', compression='gzip')
    result = dict(status='COMPUTATION_COMPLETE_VALIDATION_PENDING', cells=n, pairs=len(summaries),
                  original_states=cfg['expected_states'], layouts=results, all_cells_free=True,
                  sensitivity_used=False, labels_used_to_fit=False, initialization='same_prior_expression_coordinates',
                  parallel_umap_bitwise_reproducible=False, sampling_draws=cfg['draws'],
                  time_connected_cells=int((np.asarray(temporal.sum(1)).ravel() > 0).sum()),
                  precise_common_time_or_lineage_validated=False)
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    log('Both layouts complete; independent validation and visual review required')


if __name__ == '__main__':
    main()
