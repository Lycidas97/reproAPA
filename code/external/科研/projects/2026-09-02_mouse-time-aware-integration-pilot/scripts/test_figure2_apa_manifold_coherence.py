"""Test C1-C3 WPP coherence on the frozen within-study expression graph."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


PROJECT = Path(__file__).resolve().parents[1]
WORK = PROJECT.parents[1]
INTEGRATION = WORK / "projects/2026-09-10_human-mouse-time-aware-umap/data/results_integrated_species_scvi_v2"
STATE_ROOT = WORK / "projects/2026-09-10_human-mouse-time-aware-umap/data/results_raw_cell_species_integrated_v3_faiss2/state_and_trajectory_display_v13"
APA = PROJECT / "outputs/figure2-k3-group-umap-v1/k3_group_wpp_metacells.tsv.gz"
OUT = PROJECT / "outputs/figure2-apa-manifold-coherence-v1"
CLASSES = ("C1", "C2", "C3")


def bh(p):
    order = np.argsort(p)
    n = len(p)
    result = np.empty(n)
    result[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(result, 1)


def state_assignment():
    states = pd.read_csv(STATE_ROOT / "raw_cell_stateflow_groups.tsv.gz", sep="\t",
                         usecols=["cell_id", "state_group", "library_key"])
    roots = [PROJECT / "data/metacell-pas-v7-v1",
             WORK / "projects/2026-09-09_human-time-aware-integration/data/results_metacell_pas_v2"]
    mapping = pd.concat([pd.read_csv(root / "all_cell_to_pas_and_metacell.tsv.gz", sep="\t",
                                     usecols=["cell_id", "metacell_id"]) for root in roots], ignore_index=True)
    cells = mapping.merge(states, on="cell_id", validate="one_to_one")
    assert len(cells) == 551201 and cells.cell_id.is_unique
    assert cells.groupby("metacell_id").library_key.nunique().eq(1).all()
    counts = pd.crosstab(cells.metacell_id, cells.state_group)
    return pd.DataFrame({"metacell_id": counts.index,
                         "state_group": counts.idxmax(axis=1).to_numpy(),
                         "state_purity": (counts.max(axis=1) / counts.sum(axis=1)).to_numpy(),
                         "library_key": cells.groupby("metacell_id").library_key.first().reindex(counts.index).to_numpy()})


def pas_support():
    roots = [PROJECT / "data/metacell-pas-v7-v1/studies",
             WORK / "projects/2026-09-09_human-time-aware-integration/data/results_metacell_pas_v2/studies"]
    frames = [pd.read_csv(study / "metacells.tsv.gz", sep="\t",
                          usecols=["metacell_id", "primary_pas_umi"])
              for root in roots for study in sorted(root.iterdir())]
    result = pd.concat(frames, ignore_index=True)
    assert len(result) == 11104 and result.metacell_id.is_unique
    return result


def permuted_moran(graph, values, blocks, permutations, seed):
    n = len(values)
    centered = values - values.mean()
    denom = float(centered @ centered)
    assert n > 20 and denom > 0
    weight_sum = float(graph.sum())
    factor = n / (weight_sum * denom)
    observed = float(factor * (centered @ (graph @ centered)))
    rng = np.random.default_rng(seed)
    block_indices = [np.flatnonzero(blocks == b) for b in pd.unique(blocks)] if blocks is not None else None
    null = np.empty(permutations)
    for start in range(0, permutations, 64):
        count = min(64, permutations - start)
        shuffled = np.broadcast_to(centered[:, None], (n, count)).copy()
        if block_indices is None:
            shuffled = rng.permuted(shuffled, axis=0)
        else:
            for idx in block_indices:
                if len(idx) > 1:
                    shuffled[idx] = rng.permuted(shuffled[idx], axis=0)
        null[start:start + count] = factor * np.einsum("ij,ij->j", shuffled, graph @ shuffled)
    return dict(moran_i=observed, null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                excess_i=float(observed - null.mean()),
                permutation_z=float((observed - null.mean()) / null.std(ddof=1)) if null.std(ddof=1) else np.nan,
                empirical_p=float((1 + (null >= observed - 1e-12).sum()) / (1 + permutations)),
                null_q95=float(np.quantile(null, .95)), null_q99=float(np.quantile(null, .99)))


def main(permutations):
    OUT.mkdir(parents=True, exist_ok=True)
    nodes = pd.read_csv(INTEGRATION / "common_metacell_nodes.tsv.gz", sep="\t",
                        usecols=["metacell_id", "study_id", "species_key", "library_key"], low_memory=False)
    apa = pd.read_csv(APA, sep="\t")
    assert len(nodes) == len(apa) == 11104
    assert nodes.metacell_id.is_unique and apa.metacell_id.is_unique
    labels = state_assignment()
    assert len(labels) == len(nodes) and labels.metacell_id.is_unique
    joined = nodes.merge(apa, on="metacell_id", validate="one_to_one").merge(
        labels, on="metacell_id", validate="one_to_one")
    joined = joined.merge(pas_support(), on="metacell_id", validate="one_to_one")
    assert len(joined) == len(nodes)
    assert joined.metacell_id.eq(nodes.metacell_id).all()
    assert joined.library_key_x.eq(joined.library_key_y).all()
    assert joined.species_key.eq(joined.species).all()
    graph = sparse.load_npz(INTEGRATION / "local_expression_kernel.npz").tocsr()
    assert graph.shape == (len(nodes), len(nodes))
    graph = graph.maximum(graph.T)
    graph.data[:] = 1.0
    graph.setdiag(0)
    graph.eliminate_zeros()
    assert (graph - graph.T).nnz == 0
    assert np.all(joined.study_id.to_numpy()[graph.nonzero()[0]] == joined.study_id.to_numpy()[graph.nonzero()[1]])
    records = []
    for si, study in enumerate(pd.unique(joined.study_id)):
        study_pure = joined.study_id.eq(study) & joined.state_purity.ge(.8)
        print(study, 'state-pure', int(study_pure.sum()), flush=True)
        for ci, klass in enumerate(CLASSES):
            take = study_pure & joined[f"{klass}_centered"].notna()
            original = np.flatnonzero(take)
            local = graph[original][:, original].tocsr()
            connected = np.diff(local.indptr) > 0
            isolated = int((~connected).sum())
            original = original[connected]
            local = local[connected][:, connected].tocsr()
            frame = joined.iloc[original].copy()
            assert local.nnz and np.diff(local.indptr).min() > 0
            state = frame.state_group.to_numpy()
            library_state = (frame.library_key_x + "|" + frame.state_group).to_numpy()
            value = frame[f"{klass}_centered"].to_numpy(float)
            assert np.isfinite(value).all()
            residual = value - pd.Series(value).groupby(state).transform('mean').to_numpy()
            design = pd.get_dummies(frame.state_group, dtype=float).to_numpy()
            design = np.column_stack([design, np.log1p(frame[f'{klass}_n'].to_numpy(float)),
                                      np.log1p(frame.primary_pas_umi.to_numpy(float))])
            coverage_residual = value - design @ np.linalg.lstsq(design, value, rcond=None)[0]
            for ti, (test, input_value, blocks) in enumerate([
                    ('total_unrestricted', value, None),
                    ('state_residual_within_state', residual, state),
                    ('state_residual_within_library_state', residual, library_state),
                    ('state_and_coverage_residual_within_library_state', coverage_residual, library_state)]):
                result = permuted_moran(local, input_value, blocks, permutations,
                                        seed=260924 + 1000 * si + 100 * ci + ti)
                records.append(dict(study_id=study, species=frame.species_key.iloc[0],
                                    cluster=klass, test=test, n_metacells=len(frame),
                                    n_state_pure_before_wpp_filter=int(study_pure.sum()),
                                    n_missing_wpp=int((study_pure & joined[f'{klass}_centered'].isna()).sum()),
                                    n_isolated_after_wpp_filter=isolated,
                                    n_undirected_edges=local.nnz // 2,
                                    n_exchangeability_blocks=len(pd.unique(blocks)) if blocks is not None else 1,
                                    permutable_metacells=int(sum(len(v) for _, v in pd.Series(
                                        blocks).groupby(blocks) if len(v) > 1)) if blocks is not None else len(frame),
                                    permutations=permutations, **result))
    table = pd.DataFrame(records)
    core = table.test.isin(['total_unrestricted', 'state_residual_within_state'])
    table.loc[core, 'q_bh'] = bh(table.loc[core, 'empirical_p'].to_numpy())
    for test in ['state_residual_within_library_state',
                 'state_and_coverage_residual_within_library_state']:
        take = table.test.eq(test)
        table.loc[take, 'q_bh'] = bh(table.loc[take, 'empirical_p'].to_numpy())
    table.to_csv(OUT / 'moran_results.tsv', sep='\t', index=False, float_format='%.9g')
    qa = dict(status='PASS', analysis_unit='metacell', graph='frozen within-study 20-neighbor expression-only scVI kernel; binary symmetrized',
              graph_source=str(INTEGRATION / 'local_expression_kernel.npz'),
              map_latent_relation='metacell expression latent predates current raw-cell UMAP; this is an expression-manifold pilot, not an exact raw-cell-layout graph test',
              n_metacells_total=len(nodes), n_state_pure=int(joined.state_purity.ge(.8).sum()),
              n_studies=table.study_id.nunique(), n_clusters=3, n_tests=len(table),
              permutations=permutations, primary_bh_family=int(core.sum()),
              library_state_bh_family=21, coverage_adjusted_bh_family=21,
              coverage_covariates=['log1p(evaluable genes in group)', 'log1p(primary PAS UMI)'],
              source_wpp=str(APA), source_state=str(STATE_ROOT / 'raw_cell_stateflow_groups.tsv.gz'))
    (OUT / 'qa.json').write_text(json.dumps(qa, indent=2) + '\n')
    print(table[['study_id','cluster','test','moran_i','excess_i','empirical_p','q_bh']].to_string(index=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--permutations', type=int, default=9999)
    args = parser.parse_args()
    main(args.permutations)
