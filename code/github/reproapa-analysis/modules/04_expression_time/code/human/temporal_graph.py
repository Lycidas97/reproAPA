"""Study-local expression kernels and exact cross-study MNN, without time labels."""
import argparse
import gc
import hashlib
import itertools
import json
import os
from pathlib import Path
import time


def log(message):
    print(time.strftime('%H:%M:%S'), message, flush=True)


def normalize(matrix):
    import numpy as np
    from scipy import sparse
    totals = np.asarray(matrix.sum(1)).ravel()
    return (sparse.diags(1. / np.maximum(totals, 1e-30)) @ matrix).tocsr()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    p.add_argument('--threads', type=int, default=16)
    args = p.parse_args()
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:args.threads])
    import anndata as ad
    import numpy as np
    import pandas as pd
    from pynndescent import NNDescent
    from scipy import sparse
    from scipy.spatial import cKDTree
    from sklearn.neighbors import NearestNeighbors
    import yaml

    cfg = yaml.safe_load(Path(args.config).read_text())
    out = Path(cfg['output_root']) / 'expression'
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'summary.json').exists():
        raise FileExistsError('Refusing to overwrite completed expression kernels')
    source = Path(cfg['input_root']) / 'expression_diagnostics/human_primary_scvi.h5ad'
    a = ad.read_h5ad(source)
    obs = a.obs
    assert a.n_obs == cfg['expected_cells'] and obs.cell_id.is_unique
    z = np.asarray(a.obsm['X_scVI'], dtype=np.float32)
    mean, scale = z.mean(0), np.maximum(z.std(0), 1e-6)
    z = (z - mean) / scale
    np.savez(out / 'scaling.npz', mean=mean, scale=scale)
    n, k = len(z), cfg['local_neighbors']
    radius = np.zeros(n, dtype=np.float32)
    rows, cols, weights = [], [], []
    for group in sorted(obs[cfg['local_group']].unique()):
        ix = np.flatnonzero(obs[cfg['local_group']].eq(group))
        model = NNDescent(z[ix], n_neighbors=40, random_state=cfg['seed'], n_jobs=args.threads, low_memory=True)
        model.prepare()
        j, d = model.neighbor_graph
        keep = j != np.arange(len(ix))[:, None]
        chosen = np.argsort(~keep, axis=1, kind='stable')[:, :k]
        j, d = np.take_along_axis(j, chosen, axis=1), np.take_along_axis(d, chosen, axis=1)
        assert (j >= 0).all() and np.isfinite(d).all() and not (j == np.arange(len(ix))[:, None]).any()
        radius[ix] = d.max(1)
        rows.append(np.repeat(ix, k)); cols.append(ix[j.ravel()])
        weights.append(np.exp(-(d / max(float(np.median(d)), 1e-6)) ** 2).ravel())
        log(f'Local kernel {group}: {len(ix)} cells')
        del model
        gc.collect()
    local = sparse.csr_matrix((np.concatenate(weights), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    local = normalize(local.maximum(local.T)).astype(np.float32)
    sparse.save_npz(out / 'local_kernel.npz', local)
    np.save(out / 'local_radius.npy', radius)
    indices = {s: np.flatnonzero(obs.study_id.eq(s)) for s in sorted(obs.study_id.unique())}
    trees = {s: cKDTree(z[ix]) for s, ix in indices.items()}
    rows, cols, weights, checks = [], [], [], []
    rng = np.random.default_rng(cfg['seed'])
    k = cfg['cross_neighbors']
    for s1, s2 in itertools.combinations(indices, 2):
        aa, bb = indices[s1], indices[s2]
        d, j = trees[s2].query(z[aa], k=k, eps=0, workers=args.threads)
        bd, back = trees[s1].query(z[bb], k=k, eps=0, workers=args.threads)
        recalls = []
        for src, dst, candidates in [(aa, bb, j), (bb, aa, back)]:
            test = rng.choice(len(src), min(100, len(src)), replace=False)
            exact = NearestNeighbors(n_neighbors=k, algorithm='brute', n_jobs=args.threads).fit(z[dst]).kneighbors(z[src[test]], return_distance=False)
            recalls.append(float(np.mean([len(set(x) & set(y)) / k for x, y in zip(exact, candidates[test])])))
        assert min(recalls) >= .99, (s1, s2, recalls)
        forward = sparse.csr_matrix((np.ones(j.size, np.int8), (np.repeat(np.arange(len(aa)), k), j.ravel())), shape=(len(aa), len(bb)))
        reverse = sparse.csr_matrix((np.ones(back.size, np.int8), (back.ravel(), np.repeat(np.arange(len(bb)), k))), shape=forward.shape)
        mutual = forward.multiply(reverse).tocoo()
        r, c = aa[mutual.row], bb[mutual.col]
        distance = np.linalg.norm(z[r] - z[c], axis=1)
        ratio = distance / np.maximum(np.sqrt(radius[r] * radius[c]), 1e-6)
        valid = ratio <= cfg['cross_radius_gate']
        rows.append(r[valid]); cols.append(c[valid]); weights.append(np.exp(-ratio[valid] ** 2))
        checks.append(dict(study_a=s1, study_b=s2, mutual_candidates=len(r), retained=int(valid.sum()),
                           recall_a_to_b=recalls[0], recall_b_to_a=recalls[1]))
        log(str(checks[-1]))
    cross = sparse.csr_matrix((np.concatenate(weights), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    cross = normalize(cross.maximum(cross.T)).astype(np.float32)
    sparse.save_npz(out / 'cross_kernel.npz', cross)
    pd.DataFrame(checks).to_csv(out / 'cross_checks.tsv', sep='\t', index=False)
    obs[['cell_id', 'study_id', 'library_key', 'original_state_id']].to_csv(out / 'node_order.tsv.gz', sep='\t', index=False)
    result = dict(status='PASS', cells=n, local_group=cfg['local_group'], local_edges=local.nnz,
                  cross_edges=cross.nnz, time_or_cluster_labels_used=False, exact_cross_checks=checks,
                  input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  config_sha256=hashlib.sha256(Path(args.config).read_bytes()).hexdigest())
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
