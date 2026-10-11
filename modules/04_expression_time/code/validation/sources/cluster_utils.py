#!/usr/bin/env python3
from __future__ import annotations

import gc
import hashlib
import json
import math
import os
from pathlib import Path

import anndata as ad
import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
import scanpy as sc
from pynndescent import NNDescent
from scipy import sparse
from sklearn.decomposition import PCA


def ordered_sha256(values) -> str:
    return hashlib.sha256("".join(f"{value}\n" for value in values).encode("utf-8")).hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def write_marker(path: Path, **payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text("".join(f"{key}={value}\n" for key, value in payload.items()), encoding="utf-8")
    os.replace(temporary, path)


def canonicalize(labels) -> np.ndarray:
    labels = np.asarray(labels)
    first = sorted((int(np.flatnonzero(labels == value)[0]), value) for value in pd.unique(labels))
    mapping = {value: index for index, (_, value) in enumerate(first)}
    return np.asarray([mapping[value] for value in labels], dtype=np.int16)


def select_hvg(matrix: sparse.csr_matrix, features: pd.DataFrame, n_hvg: int, method: str):
    if n_hvg >= matrix.shape[1]:
        return np.arange(matrix.shape[1], dtype=int), {"requested_count": n_hvg, "selected_count": matrix.shape[1], "method": method}
    gene_ids = features["gene_id"].astype(str).to_numpy()
    if method == "scanpy_seurat_on_sha256_gene_discovery_split":
        source_index = np.flatnonzero([hashlib.sha256(value.encode()).digest()[0] % 5 != 0 for value in gene_ids])
    elif method in {"scanpy_seurat", "seurat_dispersion_top5000_then_rank_excluding_mt_ribo"}:
        source_index = np.arange(len(gene_ids))
    else:
        raise ValueError(f"unsupported hvg_flavor={method}")
    obj = ad.AnnData(X=matrix[:, source_index], var=features.iloc[source_index].set_index("gene_id")[["gene_symbol"]].copy())
    top = min(5000 if method.startswith("seurat_dispersion") else n_hvg, obj.n_vars)
    sc.pp.highly_variable_genes(obj, n_top_genes=top, flavor="seurat")
    if method.startswith("seurat_dispersion"):
        symbols = obj.var["gene_symbol"].astype(str)
        technical = symbols.str.lower().str.startswith("mt-") | symbols.str.match(r"^Rp[ls]", case=False)
        candidates = obj.var.loc[~technical & np.isfinite(obj.var["dispersions_norm"])].copy()
    else:
        finite = np.isfinite(obj.var["dispersions_norm"].to_numpy(float)) & np.isfinite(obj.var["dispersions"].to_numpy(float))
        candidates = obj.var.loc[finite].copy()
    candidates["gene_id_utf8"] = candidates.index.astype(str)
    if len(candidates) < n_hvg:
        raise ValueError(f"only {len(candidates)} finite HVG candidates for requested {n_hvg}")
    chosen = candidates.sort_values(
        ["dispersions_norm", "dispersions", "gene_id_utf8"],
        ascending=[False, False, True], kind="mergesort",
    ).index[:n_hvg]
    local = pd.Index(obj.var.index.astype(str)).get_indexer(chosen.astype(str))
    selected = np.sort(source_index[local].astype(int), kind="stable")
    return selected, {
        "requested_count": n_hvg,
        "selected_count": len(selected),
        "method": method,
        "selected_feature_order_sha256": ordered_sha256(features.iloc[selected]["gene_id"].astype(str)),
    }


def compute_pca(matrix, features, pipe):
    selected, audit = select_hvg(matrix, features, int(pipe["hvg_n_top_genes"]), str(pipe["hvg_flavor"]))
    dense = matrix[:, selected].toarray().astype(np.float32, copy=False)
    mean = dense.mean(axis=0, dtype=np.float64).astype(np.float32)
    std = dense.std(axis=0, dtype=np.float64).astype(np.float32)
    if np.any(std <= 0):
        raise ValueError("zero-variance selected feature")
    dense -= mean
    dense /= std
    np.clip(dense, -float(pipe.get("scale_max_value", 10)), float(pipe.get("scale_max_value", 10)), out=dense)
    solver_name = str(pipe["pca_solver"])
    solver = "randomized" if "randomized" in solver_name else "arpack"
    kwargs = {"iterated_power": 4} if solver_name == "sklearn_randomized_iterated_power_4" else {}
    n_components = int(pipe["pca_n_components_fitted"])
    model = PCA(n_components=n_components, svd_solver=solver, random_state=int(pipe["pca_seed"]), **kwargs)
    coordinates = model.fit_transform(dense).astype(np.float32, copy=False)
    del dense
    gc.collect()
    model_state = {
        "selected_feature_indices": selected.astype(np.int64),
        "feature_mean": mean,
        "feature_scale": std,
        "components": model.components_.astype(np.float32),
        "scale_max_value": np.asarray(float(pipe.get("scale_max_value", 10)), dtype=np.float32),
    }
    return coordinates, selected, model.explained_variance_ratio_.astype(np.float32), audit, model_state


def graph_from_coordinates(coords, k: int, seed: int, threads: int, backend: str):
    if len(coords) < 3:
        raise ValueError("at least three cells are required")
    k = min(int(k), len(coords) - 1)
    if backend == "pynndescent.NNDescent":
        width_requested = min(max(k + 1, 51), len(coords) - 1)
        model = NNDescent(coords, n_neighbors=width_requested, metric="euclidean", random_state=seed, n_jobs=threads, low_memory=True)
        indices, distances = model.neighbor_graph
        rows = np.arange(len(coords))[:, None]
        work = distances.copy()
        work[indices == rows] = np.inf
        order = np.argsort(work, axis=1)
        width = min(k, indices.shape[1] - 1)
        indices = np.take_along_axis(indices, order, axis=1)[:, :width]
        distances = np.take_along_axis(work, order, axis=1)[:, :width]
        row = np.repeat(np.arange(len(coords)), width)
        col = indices.reshape(-1)
        scale = np.maximum(np.median(distances, axis=1), 1e-6)
        weight = np.exp(-distances.reshape(-1) / np.repeat(scale, width)).astype(np.float32)
        directed = sparse.csr_matrix((weight, (row, col)), shape=(len(coords), len(coords)))
        graph = directed.maximum(directed.T).tocsr()
    elif backend == "scanpy.pp.neighbors_method_umap_transformer_auto":
        obj = ad.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
        obj.obsm["X_pca"] = coords
        sc.pp.neighbors(obj, n_neighbors=k, n_pcs=coords.shape[1], random_state=seed, method="umap", use_rep="X_pca")
        graph = obj.obsp["connectivities"].tocsr().astype(np.float32)
    else:
        raise ValueError(f"unsupported neighbor_backend={backend}")
    graph.setdiag(0)
    graph.eliminate_zeros()
    graph.sort_indices()
    return graph


def leiden(graph, resolution: float, seed: int, iterations: int, implementation: str):
    if implementation == "scanpy.tl.leiden_flavor_igraph":
        obj = ad.AnnData(X=np.zeros((graph.shape[0], 1), dtype=np.float32))
        sc.tl.leiden(obj, adjacency=graph, resolution=resolution, random_state=seed,
                     key_added="cluster", flavor="igraph", n_iterations=iterations, directed=False)
        return canonicalize(obj.obs["cluster"].astype(str).to_numpy())
    if implementation != "leidenalg.find_partition":
        raise ValueError(f"unsupported Leiden implementation={implementation}")
    upper = sparse.triu(graph, k=1).tocoo()
    model = ig.Graph(n=graph.shape[0], edges=list(zip(upper.row.tolist(), upper.col.tolist())), directed=False)
    model.es["weight"] = upper.data.astype(float).tolist()
    labels = leidenalg.find_partition(
        model, leidenalg.RBConfigurationVertexPartition, weights="weight",
        resolution_parameter=resolution, seed=seed, n_iterations=iterations,
    ).membership
    return canonicalize(labels)


def stratified_indices(strata, target_n: int, seed: int) -> np.ndarray:
    strata = pd.Series(strata, dtype=str)
    if strata.isna().any() or strata.isin({"", "NA", "nan"}).any():
        raise ValueError("sampling strata contain missing values")
    if target_n >= len(strata):
        return np.arange(len(strata), dtype=int)
    groups = [(name, group.index.to_numpy(dtype=int)) for name, group in strata.groupby(strata, sort=True)]
    exact = np.asarray([target_n * len(index) / len(strata) for _, index in groups])
    allocation = np.floor(exact).astype(int)
    allocation = np.maximum(allocation, 1)
    while allocation.sum() > target_n:
        candidates = np.flatnonzero(allocation > 1)
        allocation[candidates[np.argmax(allocation[candidates] - exact[candidates])]] -= 1
    remainder = target_n - allocation.sum()
    if remainder:
        order = np.argsort(-(exact - np.floor(exact)), kind="stable")
        for position in order[:remainder]:
            allocation[position] += 1
    rng = np.random.default_rng(seed)
    selected = []
    for (_, index), size in zip(groups, allocation):
        selected.extend(rng.choice(index, size=int(size), replace=False).tolist())
    return np.asarray(sorted(selected), dtype=int)
