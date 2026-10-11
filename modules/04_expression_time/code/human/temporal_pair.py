"""Full-cell, route-scaled adjacent-time OT; mirrors the accepted mouse solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True)
    p.add_argument('--pair', type=int, required=True)
    p.add_argument('--threads', type=int, default=8)
    args = p.parse_args()
    cpus = sorted(os.sched_getaffinity(0))
    start = (args.pair % 4) * args.threads
    os.sched_setaffinity(0, cpus[start:start + args.threads])
    import anndata as ad
    import numpy as np
    import pandas as pd
    import yaml
    import moscot
    from moscot.problems.time import TemporalProblem

    cfg = yaml.safe_load(Path(args.config).read_text())
    root = Path(cfg['input_root'])
    out = Path(cfg['output_root']) / 'pairs' / f'{args.pair:02d}'
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'summary.json').exists():
        raise FileExistsError('Refusing to overwrite a completed OT pair')
    plan_path = root / 'prepared/temporal_edge_plan.tsv'
    plan = pd.read_csv(plan_path, sep='\t').fillna('')
    assert len(plan) == cfg['expected_pairs']
    row = plan.iloc[args.pair]
    a = ad.read_h5ad(root / 'expression_diagnostics/human_primary_scvi.h5ad')
    assert a.n_obs == cfg['expected_cells'] and a.obs_names.is_unique
    assert a.obs.qc_tier.eq('primary_pass').all()
    route = plan[plan.route_id.eq(row.route_id)]
    libraries = set()
    for values in route[['source_libraries', 'target_libraries']].to_numpy().ravel():
        libraries.update(json.loads(values))
    route_mask = a.obs.library_key.isin(libraries).to_numpy()
    z = np.asarray(a.obsm['X_scVI'], dtype=float)
    mean = z[route_mask].mean(0)
    scale = np.maximum(z[route_mask].std(0), 1e-6)
    src = np.flatnonzero(a.obs.library_key.isin(json.loads(row.source_libraries)))
    tgt = np.flatnonzero(a.obs.library_key.isin(json.loads(row.target_libraries)))
    assert len(src) and len(tgt) and not np.intersect1d(src, tgt).size
    assert a.obs.iloc[np.r_[src, tgt]].study_id.nunique() == 1
    assert not a.obs.iloc[src].study_id.eq('GSE150311').any()
    if row.relation == 'actual_day':
        times = [float(row.time_from[1:]), float(row.time_to[1:])]
        assert times[0] < times[1] and np.isclose(times[1] - times[0], float(row.elapsed_days))
        for ix, value in zip([src, tgt], times):
            assert pd.to_numeric(a.obs.iloc[ix].actual_day).eq(value).all()
    else:
        assert row.relation == 'ordered_stage_without_elapsed_days' and row.elapsed_days == ''
        assert pd.to_numeric(a.obs.iloc[np.r_[src, tgt]].actual_day, errors='coerce').isna().all()
        times = [0., 1.]  # Pair-order solver keys only; not elapsed days.
    obs = a.obs.iloc[np.r_[src, tgt]].copy()
    obs['solver_order'] = np.r_[np.full(len(src), times[0]), np.full(len(tgt), times[1])]
    obs['mass'] = np.r_[np.full(len(src), 1. / len(src)), np.full(len(tgt), 1. / len(tgt))]
    latent = (z[np.r_[src, tgt]] - mean) / scale
    pair = ad.AnnData(X=latent, obs=obs)
    pair.obsm['X_scVI'] = latent
    print(f'PAIR {args.pair}: {row.route_id} {row.time_from} -> {row.time_to}: {len(src)} x {len(tgt)}', flush=True)
    tp = TemporalProblem(pair).prepare(time_key='solver_order', joint_attr='X_scVI', a='mass', b='mass')
    tp.solve(epsilon=cfg['epsilon'], tau_a=cfg['tau'], tau_b=cfg['tau'], scale_cost='mean',
             max_iterations=cfg['max_iterations'], threshold=cfg['threshold'], device='cpu')
    solution = tp[times[0], times[1]].solution
    t = np.asarray(solution.transport_matrix, dtype=float)
    assert t.shape == (len(src), len(tgt)) and np.isfinite(t).all() and (t >= 0).all()
    assert (t.sum(1) > 0).all() and (t.sum(0) > 0).all()
    conditional = t / t.sum(1, keepdims=True)
    assert np.allclose(conditional.sum(1), 1)
    np.savez_compressed(out / 'coupling.npz', transport=t, source=a.obs_names[src].to_numpy(str),
                        target=a.obs_names[tgt].to_numpy(str), route_mean=mean, route_scale=scale,
                        source_mass=obs.mass.to_numpy()[:len(src)], target_mass=obs.mass.to_numpy()[len(src):])
    entropy = -(conditional * np.log(np.maximum(conditional, 1e-300))).sum(1) / np.log(max(2, len(tgt)))
    summary = dict(pair=args.pair, route_id=row.route_id, time_from=row.time_from, time_to=row.time_to,
                   relation=row.relation, elapsed_days=row.elapsed_days, solver_keys=times,
                   source_cells=len(src), target_cells=len(tgt), converged=bool(solution.converged),
                   median_entropy=float(np.median(entropy)),
                   median_top5_probability=float(np.median(np.partition(conditional, -5, axis=1)[:, -5:].sum(1))),
                   transport_mass=float(t.sum()), mass_mode=cfg['mass_mode'],
                   moscot_version=moscot.__version__, epsilon=cfg['epsilon'], tau=cfg['tau'],
                   iterations_limit=cfg['max_iterations'], threshold=cfg['threshold'],
                   plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest(),
                   config_sha256=hashlib.sha256(Path(args.config).read_bytes()).hexdigest(),
                   coupling_sha256=hashlib.sha256((out / 'coupling.npz').read_bytes()).hexdigest(),
                   expression_model_sha256=a.uns['model_sha256'])
    (out / 'diagnostics.json').write_text(json.dumps(summary, indent=2) + '\n')
    assert summary['converged'], 'Nonconverged OT retained for diagnosis, not accepted'
    summary['status'] = 'PASS'
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'PAIR {args.pair} PASS', flush=True)


if __name__ == '__main__':
    main()
