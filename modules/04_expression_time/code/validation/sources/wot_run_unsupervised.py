#!/usr/bin/env python3
"""Build and freeze a program-blind WOT GEX topology before overlays."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import platform
import time
import warnings
from pathlib import Path

import anndata as ad
import igraph as ig
import leidenalg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import umap
from pynndescent import NNDescent
from scipy import sparse
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.spatial.distance import cdist
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score


STUDY = "WOT"
STUDY_ID = "GSE115943_GSE122662"
EXPECTED_TEMPORAL = 226_692
EXPECTED_ROUTES = {
    "common_early": 92_214,
    "2i_late": 66_417,
    "serum_late": 68_061,
}
SCOPES = {
    "common_trunk": ("common_early",),
    "common_to_2i": ("common_early", "2i_late"),
    "common_to_serum": ("common_early", "serum_late"),
}
N_FEATURES = 33_696
N_HVG = 2_500
N_PCS = 30
PRIMARY_K = 30
SENSITIVITY_K = 50
SEED = 218
INDEPENDENT_SEED = 219


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gex(path: Path) -> ad.AnnData:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        data = sc.read_10x_h5(path, gex_only=True)
    data.var["gene_symbol"] = data.var_names.astype(str)
    gene_ids = data.var["gene_ids"].astype(str)
    if gene_ids.duplicated().any() or not data.obs_names.is_unique:
        raise RuntimeError(f"Non-unique GEX identity in {path}")
    data.var_names = gene_ids
    return data


def neighbor_arrays(
    x: np.ndarray,
    cache: Path,
    seed: int,
    token: str,
) -> tuple[np.ndarray, np.ndarray]:
    if cache.is_file():
        stored = np.load(cache)
        cached_token = str(stored["cache_token"].item())
        if cached_token != token:
            raise RuntimeError(f"Stale neighbor cache token: {cache}")
        indices = stored["indices"]
        distances = stored["distances"]
        if indices.shape != (len(x), SENSITIVITY_K):
            raise RuntimeError(f"Stale neighbor cache shape: {cache}")
        print(f"Loaded {cache.name} {indices.shape}", flush=True)
        return indices, distances
    started = time.time()
    index = NNDescent(
        x,
        n_neighbors=SENSITIVITY_K + 1,
        metric="euclidean",
        random_state=seed,
        n_jobs=-1,
        low_memory=True,
        verbose=True,
    )
    raw_indices, raw_distances = index.neighbor_graph
    rows = np.arange(len(x))[:, None]
    work = raw_distances.copy()
    work[raw_indices == rows] = np.inf
    order = np.argsort(work, axis=1)
    indices = np.take_along_axis(raw_indices, order, axis=1)[:, :SENSITIVITY_K]
    distances = np.take_along_axis(work, order, axis=1)[:, :SENSITIVITY_K]
    if (indices < 0).any() or not np.isfinite(distances).all():
        raise RuntimeError(f"Invalid PyNNDescent result: {cache.name}")
    indices = indices.astype(np.int32, copy=False)
    distances = distances.astype(np.float32, copy=False)
    np.savez(
        cache,
        indices=indices,
        distances=distances,
        cache_token=np.asarray(token),
    )
    print(
        f"Built {cache.name} {indices.shape} in {(time.time()-started)/60:.1f} min",
        flush=True,
    )
    del index, raw_indices, raw_distances, work, order
    gc.collect()
    return indices, distances


def connectivity(
    indices: np.ndarray, distances: np.ndarray, k: int
) -> sparse.csr_matrix:
    n = len(indices)
    rows = np.repeat(np.arange(n, dtype=np.int32), k)
    cols = indices[:, :k].reshape(-1)
    local_scale = np.maximum(np.median(distances[:, :k], axis=1), 1.0e-6)
    weights = np.exp(
        -distances[:, :k].reshape(-1) / np.repeat(local_scale, k)
    ).astype(np.float32)
    graph = sparse.csr_matrix((weights, (rows, cols)), shape=(n, n))
    graph = graph.maximum(graph.T).tocsr()
    graph.setdiag(0.0)
    graph.eliminate_zeros()
    return graph


def graph_stats(
    scope: str,
    variant: str,
    indices: np.ndarray,
    distances: np.ndarray,
    k: int,
    graph: sparse.csr_matrix,
) -> dict[str, object]:
    components, component_labels = connected_components(graph, directed=False)
    counts = np.bincount(component_labels)
    degree = np.diff(graph.indptr)
    directed = len(indices) * k
    undirected = graph.nnz // 2
    mutual_pairs = max(directed - undirected, 0)
    return {
        "study": STUDY,
        "scope": scope,
        "record_type": "graph_variant",
        "variant": variant,
        "n_cells": len(indices),
        "k": k,
        "components": int(components),
        "giant_component_fraction": float(counts.max() / len(indices)),
        "median_degree": float(np.median(degree)),
        "p05_degree": float(np.quantile(degree, 0.05)),
        "p95_degree": float(np.quantile(degree, 0.95)),
        "mutual_directed_neighbor_fraction": float(2 * mutual_pairs / directed),
        "median_neighbor_distance": float(np.median(distances[:, :k])),
        "p95_neighbor_distance": float(np.quantile(distances[:, :k], 0.95)),
    }


def leiden_partition(graph: sparse.csr_matrix, seed: int) -> tuple[np.ndarray, float]:
    upper = sparse.triu(graph, k=1).tocoo()
    edge_list = list(zip(upper.row.tolist(), upper.col.tolist()))
    community_graph = ig.Graph(n=graph.shape[0], edges=edge_list, directed=False)
    community_graph.es["weight"] = upper.data.astype(float).tolist()
    partition = leidenalg.find_partition(
        community_graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=0.45,
        seed=seed,
        n_iterations=-1,
    )
    labels = np.asarray(partition.membership, dtype=np.int32)
    modularity = float(partition.modularity)
    del upper, edge_list, community_graph, partition
    gc.collect()
    return labels, modularity


def sampled_neighbor_jaccard(
    left: np.ndarray, right: np.ndarray, seed: int, n_sample: int = 5_000
) -> float:
    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(left), size=min(n_sample, len(left)), replace=False)
    values = []
    for index in chosen:
        a = set(left[index, :PRIMARY_K].tolist())
        b = set(right[index, :PRIMARY_K].tolist())
        values.append(len(a & b) / len(a | b))
    return float(np.median(values))


def local_intrinsic_dimension(
    x: np.ndarray, indices: np.ndarray, seed: int, n_sample: int = 1_500
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(x), size=min(n_sample, len(x)), replace=False)
    dimensions = []
    for index in chosen:
        local = x[indices[index, :PRIMARY_K]]
        local = local - local.mean(axis=0, keepdims=True)
        eigenvalues = np.linalg.eigvalsh(local.T @ local)
        eigenvalues = np.maximum(eigenvalues[::-1], 0.0)
        total = eigenvalues.sum()
        dimensions.append(
            int(np.searchsorted(np.cumsum(eigenvalues) / total, 0.90) + 1)
            if total > 0
            else 0
        )
    return float(np.median(dimensions)), float(np.quantile(dimensions, 0.90))


def skeleton_metrics(x: np.ndarray, seed: int) -> dict[str, float]:
    results: dict[str, float] = {}
    for n_clusters in (20, 30, 40):
        model = MiniBatchKMeans(
            n_clusters=n_clusters,
            random_state=seed,
            batch_size=4096,
            n_init=5,
        ).fit(x)
        distances = cdist(model.cluster_centers_, model.cluster_centers_)
        np.fill_diagonal(distances, 0.0)
        tree = minimum_spanning_tree(distances).tocoo()
        degree = np.bincount(
            np.concatenate([tree.row, tree.col]), minlength=n_clusters
        )
        results[f"kmeans{n_clusters}_mst_max_degree"] = float(degree.max())
        results[f"kmeans{n_clusters}_mst_branch_nodes"] = float((degree >= 3).sum())
        results[f"kmeans{n_clusters}_mst_leaf_fraction"] = float(
            (degree == 1).mean()
        )
    return results


def cluster_summary(
    scope: str,
    cells: pd.DataFrame,
    positions: np.ndarray,
    x: np.ndarray,
    embedding: np.ndarray,
    labels: np.ndarray,
    indices: np.ndarray,
) -> pd.DataFrame:
    same_neighbor = (labels[indices[:, :PRIMARY_K]] == labels[:, None]).mean(axis=1)
    frame = pd.DataFrame(
        {
            "cluster": labels,
            "PC1": x[:, 0],
            "PC2": x[:, 1],
            "embedding_x": embedding[:, 0],
            "embedding_y": embedding[:, 1],
            "same_cluster_neighbor_fraction": same_neighbor,
        }
    )
    summary = (
        frame.groupby("cluster", sort=True)
        .agg(
            n_cells=("cluster", "size"),
            PC1_centroid=("PC1", "mean"),
            PC2_centroid=("PC2", "mean"),
            embedding_x_centroid=("embedding_x", "mean"),
            embedding_y_centroid=("embedding_y", "mean"),
            median_same_cluster_neighbor_fraction=(
                "same_cluster_neighbor_fraction",
                "median",
            ),
        )
        .reset_index()
    )
    summary.insert(0, "scope", scope)
    summary["cluster"] = summary["cluster"].map(lambda value: f"L{value:02d}")
    summary["fraction_of_scope"] = summary["n_cells"] / len(positions)
    return summary


def save_unlabeled_figure(
    coordinate_frames: list[pd.DataFrame], output_dir: Path
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    rng = np.random.default_rng(SEED)
    for axis, frame in zip(axes, coordinate_frames):
        if len(frame) > 60_000:
            chosen = rng.choice(len(frame), size=60_000, replace=False)
            shown = frame.iloc[chosen]
        else:
            shown = frame
        cluster_codes = pd.Categorical(shown["unsupervised_cluster"]).codes
        axis.scatter(
            shown["embedding_x"],
            shown["embedding_y"],
            c=cluster_codes,
            cmap="turbo",
            s=0.35,
            alpha=0.42,
            linewidths=0,
            rasterized=True,
        )
        axis.set_title(f"{frame['scope'].iloc[0]} (label-free Leiden colors)")
        axis.set_xlabel("UMAP1")
        axis.set_ylabel("UMAP2")
        axis.set_xticks([])
        axis.set_yticks([])
    fig.suptitle("UST01 Phase A: program-blind route-local GEX topology")
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / "figures" / f"unlabeled_topology.{suffix}", dpi=220)
    plt.close(fig)


def save_sensitivity_figure(metrics: pd.DataFrame, output_dir: Path) -> None:
    scope_rows = metrics.loc[metrics["record_type"].eq("scope_stability")].copy()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    fields = [
        ("neighbor_jaccard_seed218_vs_seed219", "median 30-NN Jaccard", (0, 1)),
        ("leiden_ARI_seed218_vs_seed219", "Leiden ARI", (0, 1)),
        ("median_local_dimension_90pct", "local dimension (90% variance)", None),
    ]
    for axis, (field, title, limits) in zip(axes, fields):
        axis.bar(scope_rows["scope"], scope_rows[field], color="#4C78A8")
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=25)
        if limits is not None:
            axis.set_ylim(*limits)
        axis.grid(axis="y", alpha=0.25)
    fig.suptitle("UST01 Phase A: scope and parameter stability")
    for suffix in ("png", "pdf"):
        fig.savefig(
            output_dir / "figures" / f"scope_parameter_sensitivity.{suffix}",
            dpi=220,
        )
    plt.close(fig)


def software_table() -> pd.DataFrame:
    packages = [
        "anndata",
        "igraph",
        "leidenalg",
        "matplotlib",
        "numpy",
        "pandas",
        "pynndescent",
        "scanpy",
        "scikit-learn",
        "scipy",
        "umap-learn",
    ]
    rows = [{"software": "python", "version": platform.python_version()}]
    for package in packages:
        rows.append(
            {
                "software": package,
                "version": importlib.metadata.version(package),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    intermediate_dir = output_dir / "intermediate"
    figures_dir.mkdir(exist_ok=True)
    intermediate_dir.mkdir(exist_ok=True)
    started = time.time()

    issue_manifest = pd.read_csv(output_dir / "issue_inputs/UPLOAD_MANIFEST.tsv", sep="\t")
    for row in issue_manifest.itertuples(index=False):
        path = output_dir / "issue_inputs" / row.filename
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256(path) != row.sha256:
            raise RuntimeError(f"Locked issue input mismatch: {row.filename}")

    upstream = Path("analyses/DOC-218-wot-full-cell-expression-pseudotime-v4").resolve()
    cell_path = upstream / "cell_full_graph_pseudotime_v4.parquet"
    input_audit_path = upstream / "input_audit.tsv"
    allowed_columns = [
        "cell_id",
        "barcode",
        "study_id",
        "library_id",
        "sample_id",
        "route_id",
        "graph_in_sample",
    ]
    cells = pd.read_parquet(cell_path, columns=allowed_columns)
    cells = cells.loc[cells["graph_in_sample"]].drop(columns="graph_in_sample").copy()
    if len(cells) != EXPECTED_TEMPORAL or not cells["study_id"].eq(STUDY_ID).all():
        raise RuntimeError("Unexpected accepted temporal cell universe")
    if cells["cell_id"].duplicated().any() or cells.duplicated(["library_id", "barcode"]).any():
        raise RuntimeError("Duplicate accepted cell identity")
    observed_routes = cells["route_id"].value_counts().to_dict()
    if observed_routes != EXPECTED_ROUTES:
        raise RuntimeError(f"Unexpected temporal route counts: {observed_routes}")

    input_audit = pd.read_csv(input_audit_path, sep="\t")
    input_audit = input_audit.loc[input_audit["record_type"].eq("library_input")].copy()
    if input_audit["library_id"].duplicated().any():
        raise RuntimeError("Duplicate library in upstream input audit")
    input_audit = input_audit.sort_values("library_id").reset_index(drop=True)
    cells_by_library = {
        key: value.copy() for key, value in cells.groupby("library_id", sort=False)
    }

    feature_ids: np.ndarray | None = None
    feature_table: pd.DataFrame | None = None
    count_parts: list[sparse.csr_matrix] = []
    cell_parts: list[pd.DataFrame] = []
    provenance_rows: list[dict[str, object]] = []
    for row in input_audit.itertuples(index=False):
        raw_path = Path(row.filtered_h5)
        raw_hash = sha256(raw_path)
        selected = cells_by_library.get(row.library_id)
        n_selected = 0 if selected is None else len(selected)
        provenance_rows.append(
            {
                "record_type": "raw_GEX_library",
                "path": str(raw_path),
                "sha256": raw_hash,
                "bytes": raw_path.stat().st_size,
                "library_id": row.library_id,
                "n_selected_cells": n_selected,
                "columns_or_layer_read": "10x Gene Expression only",
                "use": "program-blind Phase A GEX",
                "access": "read_only",
            }
        )
        if selected is None:
            continue
        data = read_gex(raw_path)
        observed_ids = data.var_names.to_numpy(dtype=str)
        if feature_ids is None:
            feature_ids = observed_ids
            feature_table = data.var.copy()
        elif not np.array_equal(feature_ids, observed_ids):
            raise RuntimeError(f"Feature order mismatch: {row.library_id}")
        missing = sorted(set(selected["barcode"]) - set(data.obs_names))
        if missing:
            raise RuntimeError(f"Missing selected barcodes in {row.library_id}: {len(missing)}")
        count_parts.append(data[selected["barcode"].to_numpy(dtype=str), :].X.tocsr())
        cell_parts.append(selected.reset_index(drop=True))
        del data
    assert feature_ids is not None and feature_table is not None
    cells = pd.concat(cell_parts, ignore_index=True)
    counts = sparse.vstack(count_parts, format="csr").astype(np.float32)
    del count_parts, cell_parts
    gc.collect()
    if counts.shape != (EXPECTED_TEMPORAL, N_FEATURES):
        raise RuntimeError(f"Unexpected assembled raw GEX shape: {counts.shape}")

    totals = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise RuntimeError("Accepted temporal cell has zero GEX counts")
    normalized = counts.multiply((1.0e4 / totals)[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    del counts, totals
    gc.collect()

    hvg = ad.AnnData(X=normalized, var=feature_table.copy())
    sc.pp.highly_variable_genes(hvg, n_top_genes=5_000, flavor="seurat")
    symbols = hvg.var["gene_symbol"].astype(str)
    technical_excluded = symbols.str.lower().str.startswith("mt-") | symbols.str.match(
        r"^Rp[ls]", case=False
    )
    candidates = hvg.var.loc[~technical_excluded].copy()
    candidates = candidates.loc[np.isfinite(candidates["dispersions_norm"])]
    candidates = candidates.sort_values(
        ["dispersions_norm", "dispersions"], ascending=False
    )
    selected_ids = candidates.head(N_HVG).index.astype(str).tolist()
    if len(selected_ids) != N_HVG:
        raise RuntimeError("Could not select 2,500 program-blind HVGs")
    selected_indices = hvg.var_names.get_indexer(selected_ids)
    if (selected_indices < 0).any():
        raise RuntimeError("HVG identity failure")
    feature_manifest = pd.DataFrame(
        {
            "gene_id": hvg.var_names.astype(str),
            "gene_symbol": symbols.to_numpy(dtype=str),
            "dispersion_normalized": hvg.var["dispersions_norm"].to_numpy(),
            "mitochondrial_or_ribosomal_excluded": technical_excluded.to_numpy(),
        }
    )
    feature_manifest["primary_HVG"] = feature_manifest["gene_id"].isin(selected_ids)
    feature_manifest["selection_role"] = np.where(
        feature_manifest["primary_HVG"], "program_blind_HVG", "not_primary_HVG"
    )
    feature_manifest.to_csv(
        output_dir / "unsupervised_feature_manifest.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    del hvg
    gc.collect()

    dense = normalized[:, selected_indices].toarray().astype(np.float32, copy=False)
    del normalized
    gc.collect()
    scale_mean = dense.mean(axis=0, dtype=np.float64).astype(np.float32)
    scale_std = dense.std(axis=0, dtype=np.float64).astype(np.float32)
    if np.any(scale_std <= 0):
        raise RuntimeError("Zero-variance program-blind HVG")
    dense -= scale_mean[None, :]
    dense /= scale_std[None, :]
    np.clip(dense, -10.0, 10.0, out=dense)
    pca_model = PCA(
        n_components=N_PCS,
        svd_solver="randomized",
        random_state=SEED,
        iterated_power=4,
    )
    pca = pca_model.fit_transform(dense).astype(np.float32, copy=False)
    del dense
    gc.collect()
    pca_path = intermediate_dir / "program_blind_temporal_pca30.npy"
    np.save(pca_path, pca)
    pca_token = sha256(pca_path) + ":" + sha256(cell_path)
    pd.DataFrame(
        {
            "PC": np.arange(1, N_PCS + 1),
            "explained_variance": pca_model.explained_variance_,
            "explained_variance_ratio": pca_model.explained_variance_ratio_,
        }
    ).to_csv(
        output_dir / "pca_variance.tsv", sep="\t", index=False, lineterminator="\n"
    )

    coordinate_frames: list[pd.DataFrame] = []
    cluster_tables: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
    for scope_index, (scope, routes) in enumerate(SCOPES.items()):
        positions = np.flatnonzero(cells["route_id"].isin(routes).to_numpy())
        x = np.ascontiguousarray(pca[positions])
        primary_cache = intermediate_dir / f"knn_{scope}_program_blind_seed218_k50.npz"
        independent_cache = intermediate_dir / f"knn_{scope}_program_blind_seed219_k50.npz"
        primary_indices, primary_distances = neighbor_arrays(
            x, primary_cache, SEED, pca_token + f":{scope}:PCA30"
        )
        independent_indices, independent_distances = neighbor_arrays(
            x, independent_cache, INDEPENDENT_SEED, pca_token + f":{scope}:PCA30"
        )
        pca20_cache = intermediate_dir / f"knn_{scope}_program_blind_pca20_seed218_k50.npz"
        pca20_indices, pca20_distances = neighbor_arrays(
            np.ascontiguousarray(x[:, :20]),
            pca20_cache,
            SEED,
            pca_token + f":{scope}:PCA20",
        )

        primary_graph = connectivity(primary_indices, primary_distances, PRIMARY_K)
        k50_graph = connectivity(primary_indices, primary_distances, SENSITIVITY_K)
        independent_graph = connectivity(
            independent_indices, independent_distances, PRIMARY_K
        )
        pca20_graph = connectivity(pca20_indices, pca20_distances, PRIMARY_K)
        metric_rows.extend(
            [
                graph_stats(scope, "primary_PCA30_seed218_k30", primary_indices, primary_distances, PRIMARY_K, primary_graph),
                graph_stats(scope, "neighbor_k50", primary_indices, primary_distances, SENSITIVITY_K, k50_graph),
                graph_stats(scope, "independent_seed219_k30", independent_indices, independent_distances, PRIMARY_K, independent_graph),
                graph_stats(scope, "PCA20_seed218_k30", pca20_indices, pca20_distances, PRIMARY_K, pca20_graph),
            ]
        )

        labels, modularity = leiden_partition(primary_graph, SEED)
        independent_labels, _ = leiden_partition(independent_graph, INDEPENDENT_SEED)
        ari = float(adjusted_rand_score(labels, independent_labels))
        neighbor_jaccard = sampled_neighbor_jaccard(
            primary_indices, independent_indices, SEED + scope_index
        )
        pca20_jaccard = sampled_neighbor_jaccard(
            primary_indices, pca20_indices, SEED + 10 + scope_index
        )
        median_dimension, p90_dimension = local_intrinsic_dimension(
            x, primary_indices, SEED + scope_index
        )
        variances = np.var(x, axis=0, dtype=np.float64)
        participation_ratio = float(variances.sum() ** 2 / np.square(variances).sum())
        skeleton = skeleton_metrics(x, SEED)
        stability_row: dict[str, object] = {
            "study": STUDY,
            "scope": scope,
            "record_type": "scope_stability",
            "variant": "cross_variant",
            "n_cells": len(x),
            "k": PRIMARY_K,
            "n_leiden_clusters": int(labels.max() + 1),
            "leiden_modularity": modularity,
            "leiden_ARI_seed218_vs_seed219": ari,
            "neighbor_jaccard_seed218_vs_seed219": neighbor_jaccard,
            "neighbor_jaccard_PCA30_vs_PCA20": pca20_jaccard,
            "median_local_dimension_90pct": median_dimension,
            "p90_local_dimension_90pct": p90_dimension,
            "global_PC_participation_ratio": participation_ratio,
            "PC1_fraction_first5_variance": float(variances[0] / variances[:5].sum()),
            "PC1_PC2_fraction_first5_variance": float(variances[:2].sum() / variances[:5].sum()),
            **skeleton,
        }
        metric_rows.append(stability_row)

        mapper = umap.UMAP(
            n_neighbors=PRIMARY_K,
            n_components=2,
            min_dist=0.30,
            spread=1.0,
            metric="euclidean",
            random_state=SEED,
            n_jobs=1,
            low_memory=True,
            init="spectral",
            precomputed_knn=(
                primary_indices[:, :PRIMARY_K],
                primary_distances[:, :PRIMARY_K],
                None,
            ),
        )
        embedding = mapper.fit_transform(x).astype(np.float32, copy=False)
        frame = pd.DataFrame(
            {
                "cell_id": cells.iloc[positions]["cell_id"].to_numpy(),
                "scope": scope,
                "route_id": cells.iloc[positions]["route_id"].to_numpy(),
                "embedding_x": embedding[:, 0],
                "embedding_y": embedding[:, 1],
                "PC1": x[:, 0],
                "PC2": x[:, 1],
                "unsupervised_cluster": [f"L{value:02d}" for value in labels],
            }
        )
        coordinate_frames.append(frame)
        cluster_tables.append(
            cluster_summary(
                scope, cells, positions, x, embedding, labels, primary_indices
            )
        )
        del (
            x,
            primary_indices,
            primary_distances,
            independent_indices,
            independent_distances,
            pca20_indices,
            pca20_distances,
            primary_graph,
            k50_graph,
            independent_graph,
            pca20_graph,
            labels,
            independent_labels,
            embedding,
            mapper,
        )
        gc.collect()
        print(f"Completed Phase A scope {scope}", flush=True)

    coordinates = pd.concat(coordinate_frames, ignore_index=True)
    coordinates.to_csv(
        output_dir / "unsupervised_cell_coordinates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
        lineterminator="\n",
    )
    clusters = pd.concat(cluster_tables, ignore_index=True)
    clusters.to_csv(
        output_dir / "unsupervised_cluster_summary.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(
        output_dir / "topology_metrics.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    save_unlabeled_figure(coordinate_frames, output_dir)
    save_sensitivity_figure(metrics, output_dir)

    stability = metrics.loc[metrics["record_type"].eq("scope_stability")]
    stable_graph = bool(
        (stability["leiden_ARI_seed218_vs_seed219"] >= 0.75).all()
        and (stability["neighbor_jaccard_seed218_vs_seed219"] >= 0.40).all()
    )
    low_dimensional = bool((stability["median_local_dimension_90pct"] <= 3).all())
    topology_class = "unresolved"
    pseudotime_eligible = False
    metric_verdict = (
        "stable low-dimensional route-local graph; manual unlabeled review required"
        if stable_graph and low_dimensional
        else "no sufficiently stable root-free one-dimensional/branched scaffold across all route scopes"
    )
    decision = pd.DataFrame(
        [
            {
                "study": STUDY,
                "scope": "route-local full-cell graph",
                "topology_class": topology_class,
                "pseudotime_eligible": pseudotime_eligible,
                "orientation_source": "none_no_root_or_pseudotime_in_phase_A",
                "visual_verdict": "pending_manual_unlabeled_review",
                "metric_verdict": metric_verdict,
                "conflict_flag": False,
                "limitations": "UMAP is descriptive; route metadata defines scope; no biological replicate-level inference; program-blind graph rebuilt because upstream HVGs were program-informed",
            }
        ]
    )
    decision.to_csv(
        output_dir / "topology_decision.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    provenance_rows.extend(
        [
            {
                "record_type": "cell_scope_table",
                "path": str(cell_path),
                "sha256": sha256(cell_path),
                "bytes": cell_path.stat().st_size,
                "library_id": "all",
                "n_selected_cells": EXPECTED_TEMPORAL,
                "columns_or_layer_read": ",".join(allowed_columns),
                "use": "accepted temporal cell identity plus route scope only",
                "access": "read_only_column_projection",
            },
            {
                "record_type": "raw_input_locator",
                "path": str(input_audit_path),
                "sha256": sha256(input_audit_path),
                "bytes": input_audit_path.stat().st_size,
                "library_id": "all",
                "n_selected_cells": EXPECTED_TEMPORAL,
                "columns_or_layer_read": "library_id,filtered_h5",
                "use": "locate immutable 10x GEX inputs",
                "access": "read_only",
            },
        ]
    )
    pd.DataFrame(provenance_rows).to_csv(
        output_dir / "input_provenance.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    software_table().to_csv(
        output_dir / "environment_versions.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    artifact_paths = [
        output_dir / "analysis_plan.md",
        output_dir / "upstream_independence_audit.tsv",
        output_dir / "input_provenance.tsv",
        output_dir / "unsupervised_feature_manifest.tsv",
        output_dir / "pca_variance.tsv",
        output_dir / "unsupervised_cell_coordinates.tsv.gz",
        output_dir / "unsupervised_cluster_summary.tsv",
        output_dir / "topology_metrics.tsv",
        output_dir / "topology_decision.tsv",
        output_dir / "environment_versions.tsv",
        figures_dir / "unlabeled_topology.png",
        figures_dir / "unlabeled_topology.pdf",
        figures_dir / "scope_parameter_sensitivity.png",
        figures_dir / "scope_parameter_sensitivity.pdf",
        pca_path,
        *sorted(intermediate_dir.glob("knn_*.npz")),
        Path(__file__).resolve(),
    ]
    artifacts = [
        {
            "path": str(path.relative_to(output_dir)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in artifact_paths
    ]
    forbidden_fields = [
        "stage/day or observed time",
        "author/local-state labels",
        "marker or expression-program definitions/scores",
        "transferred reference/fate labels",
        "P0/P12 or local-state eligibility masks",
        "PAS/APA features, matrices, summaries, or results",
    ]
    freeze = {
        "contract": "UST01 unsupervised-first topology audit v1",
        "created_utc": pd.Timestamp.now(tz="UTC").isoformat(),
        "study": STUDY,
        "cell_universe": {
            "n_temporal_cells": EXPECTED_TEMPORAL,
            "route_counts": EXPECTED_ROUTES,
            "scope_counts": {
                scope: int(coordinates["scope"].eq(scope).sum()) for scope in SCOPES
            },
        },
        "expression": {
            "layer": "10x Gene Expression raw counts only",
            "normalization": "library_size_10000_log1p",
            "HVG": "top_2500_seurat_dispersion_program_blind_mt_ribo_technical_exclusion",
            "PCA": "scaled_clipped10_randomized_PCA30",
        },
        "parameters": {
            "primary_k": PRIMARY_K,
            "sensitivity_k": SENSITIVITY_K,
            "primary_seed": SEED,
            "independent_seed": INDEPENDENT_SEED,
            "scope_routes": SCOPES,
            "clustering": "Leiden_RBConfiguration_resolution_0.45",
            "embedding": "UMAP_on_frozen_primary_neighbors_min_dist_0.30",
            "root": None,
            "pseudotime": None,
        },
        "phase_A_allowed_inputs": [
            "raw GEX counts and gene identity",
            "accepted primary-pass temporal cell identity",
            "library/sample identity for alignment only",
            "route_id for preregistered graph scope only",
        ],
        "phase_A_forbidden_fields_not_read": forbidden_fields,
        "upstream_reuse_decision": "rebuild; preferred v3 absent and v4 HVG set program-informed",
        "topology_decision": decision.to_dict("records")[0],
        "artifacts": artifacts,
        "elapsed_minutes": (time.time() - started) / 60,
    }
    freeze_path = output_dir / "pre_overlay_freeze.json"
    freeze_path.write_text(
        json.dumps(freeze, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_hash = sha256(freeze_path)
    (output_dir / "pre_overlay_freeze.sha256").write_text(
        f"{freeze_hash}  pre_overlay_freeze.json\n", encoding="utf-8"
    )
    print(
        f"Phase A frozen: {freeze_hash}; elapsed={(time.time()-started)/60:.1f} min",
        flush=True,
    )


if __name__ == "__main__":
    main()
