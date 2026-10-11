#!/usr/bin/env python3
"""Build an all-raw-cell human/mouse time-constrained integrated map.

Every computational node and every plotted point is a primary single cell.
The pre-existing human and mouse all-cell time graphs are retained as frozen
within-species constraints. A new species-balanced scVI latent is used only to
make expression-only human--mouse MNN links. No metacell coordinate, edge, or
back-projection enters this workflow.
"""

import argparse
import gc
import json
import os
from pathlib import Path

import anndata as ad
import faiss
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
import scvi
from scipy import sparse
from scipy.linalg import orthogonal_procrustes
from scipy.sparse.csgraph import connected_components
import torch
from umap.umap_ import find_ab_params, simplicial_set_embedding


def row_normalize(matrix):
    matrix = sparse.csr_matrix(matrix, dtype=np.float32)
    total = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(total <= 0) or not np.isfinite(total).all():
        raise ValueError("A layout graph row has no finite positive mass")
    return (sparse.diags(1.0 / total) @ matrix).tocsr()


def decode(values):
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values], dtype=str)


def read_one_to_one_features(project, cfg):
    registry = pd.read_csv(project / cfg["orthologs"], sep="\t", dtype=str).fillna("")
    registry = registry.loc[registry.documented_one_to_one.str.lower().eq("true")].copy()
    if registry.mouse_gene_id.duplicated().any() or registry.human_gene_id.duplicated().any():
        raise ValueError("One-to-one registry is not unique")
    human = ad.read_h5ad(cfg["human"]["reference_counts"], backed="r")
    first_mouse = ad.read_h5ad(Path(cfg["mouse"]["library_hvg_root"]) / "GSM3195742.h5ad", backed="r")
    human_genes = set(human.var_names.astype(str))
    mouse_genes = set(first_mouse.var.gene_ids.astype(str))
    features = registry.loc[
        registry.human_gene_id.isin(human_genes) & registry.mouse_gene_id.isin(mouse_genes)
    ].copy().reset_index(drop=True)
    if len(features) < 1000:
        raise ValueError(f"Insufficient raw-cell shared feature set: {len(features)}")
    return features


def load_human_raw_cells(cfg, features):
    """Read exact primary cells and counts from the source 10x H5 files."""
    layout = ad.read_h5ad(cfg["human"]["time_h5ad"], backed="r")
    obs = layout.obs.copy().reset_index(drop=True)
    if len(obs) != cfg["expected"]["human_cells"] or not obs.cell_id.is_unique:
        raise ValueError("Human frozen time-layout cell registry is invalid")
    manifest = pd.read_csv(cfg["human"]["libraries"], sep="\t", dtype=str)
    manifest = manifest.set_index("library_key")
    parts, ids = [], []
    for library in sorted(obs.library_key.unique()):
        desired = obs.loc[obs.library_key.eq(library), ["cell_id", "barcode"]]
        if library not in manifest.index:
            raise ValueError(f"Human source manifest lacks {library}")
        source = Path(manifest.at[library, "input_path"])
        with h5py.File(source, "r") as handle:
            matrix = handle["matrix"]
            genes = decode(matrix["features/id"][:])
            barcodes = decode(matrix["barcodes"][:])
            feature_index = pd.Index(genes).get_indexer(features.human_gene_id)
            cell_index = pd.Index(barcodes).get_indexer(desired.barcode.astype(str))
            if (feature_index < 0).any() or (cell_index < 0).any():
                raise ValueError(f"Human raw source does not close feature/cell membership: {library}")
            full = sparse.csc_matrix(
                (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                shape=tuple(matrix["shape"][:]), dtype=np.float32,
            ).T.tocsr()
            parts.append(full[cell_index][:, feature_index].tocsr())
            ids.append(desired.cell_id.astype(str).to_numpy())
    counts = sparse.vstack(parts, format="csr")
    ordered_ids = np.concatenate(ids)
    position = pd.Index(ordered_ids).get_indexer(obs.cell_id.astype(str))
    if (position < 0).any() or len(np.unique(position)) != len(position):
        raise ValueError("Human raw source cell order cannot be reconciled to frozen time order")
    counts = counts[position]
    if counts.shape != (len(obs), len(features)) or not np.isfinite(counts.data).all():
        raise ValueError("Human raw count matrix is invalid")
    obs["species_key"] = "human"
    obs["route"] = obs.study_id.astype(str)
    return obs, counts, np.asarray(layout.obsm["X_umap"], dtype=np.float32)


def load_mouse_raw_cells(cfg, features):
    """Read frozen full-cell raw-HVG counts; they are not metacell aggregates."""
    obs = pd.read_csv(cfg["mouse"]["coordinates"], sep="\t", dtype={"cell_id": str, "library_key": str})
    if len(obs) != cfg["expected"]["mouse_cells"] or not obs.cell_id.is_unique:
        raise ValueError("Mouse frozen time-layout cell registry is invalid")
    root = Path(cfg["mouse"]["library_hvg_root"])
    parts, ids, reference = [], [], None
    for library in sorted(obs.library_key.unique()):
        current = ad.read_h5ad(root / f"{library}.h5ad")
        genes = current.var.gene_ids.astype(str).to_numpy()
        if reference is None:
            reference = genes
            feature_index = pd.Index(genes).get_indexer(features.mouse_gene_id)
            if (feature_index < 0).any():
                raise ValueError("Mouse raw-HVG source lacks requested one-to-one features")
        elif not np.array_equal(genes, reference):
            raise ValueError(f"Mouse raw-HVG feature order differs in {library}")
        expected = obs.loc[obs.library_key.eq(library), "cell_id"].astype(str).to_numpy()
        if set(current.obs.cell_id.astype(str)) != set(expected):
            raise ValueError(f"Mouse raw-HVG membership differs from frozen time layout: {library}")
        parts.append(sparse.csr_matrix(current.X[:, feature_index], dtype=np.float32))
        ids.append(current.obs.cell_id.astype(str).to_numpy())
    counts = sparse.vstack(parts, format="csr")
    ordered_ids = np.concatenate(ids)
    position = pd.Index(ordered_ids).get_indexer(obs.cell_id.astype(str))
    if (position < 0).any() or len(np.unique(position)) != len(position):
        raise ValueError("Mouse raw-HVG cell order cannot be reconciled to frozen time order")
    counts = counts[position]
    if counts.shape != (len(obs), len(features)) or not np.isfinite(counts.data).all():
        raise ValueError("Mouse raw count matrix is invalid")
    obs["species_key"] = "mouse"
    obs["route"] = obs.dataset.astype(str)
    return obs, counts, obs[["UMAP1", "UMAP2"]].to_numpy(np.float32)


def fit_species_scvi(counts, nodes, features, spec, seed, output):
    if counts.data.size and not np.allclose(counts.data, np.rint(counts.data)):
        raise ValueError("Raw-cell scVI input has non-integer counts")
    var = features[["human_gene_id", "mouse_gene_id", "human_gene_symbol", "mouse_gene_symbol"]].copy()
    var.index = pd.Index([f"ortholog_{i:05d}" for i in range(len(var))], name="ortholog_feature")
    data = ad.AnnData(X=counts, obs=nodes[["species_key"]].copy(), var=var)
    data.obs_names = nodes.cell_id.astype(str)
    scvi.settings.seed = seed
    torch.set_float32_matmul_precision("high")
    scvi.model.SCVI.setup_anndata(data, batch_key=spec["batch_key"])
    model = scvi.model.SCVI(data, n_latent=spec["n_latent"], gene_likelihood="nb")
    model.train(max_epochs=spec["max_epochs"], early_stopping=spec["early_stopping"],
                accelerator="gpu", devices=1, batch_size=spec["batch_size"])
    latent = model.get_latent_representation(batch_size=spec["batch_size"]).astype(np.float32)
    if latent.shape != (len(nodes), spec["n_latent"]) or not np.isfinite(latent).all():
        raise ValueError("Raw-cell species-balanced scVI latent is invalid")
    model.save(output / "species_balanced_raw_cell_scvi_model", overwrite=False)
    np.save(output / "species_balanced_raw_cell_scvi_latent.npy", latent)
    var.reset_index().to_csv(output / "one_to_one_raw_cell_features.tsv", sep="\t", index=False)
    return latent


def faiss_ivf_neighbors(database, queries, k, spec, seed):
    """Deterministic IVF-Flat search, scaled for the full raw-cell atlas."""
    if database.dtype != np.float32 or queries.dtype != np.float32:
        raise TypeError("FAISS inputs must be float32")
    nlist = min(int(spec["nlist"]), max(1, len(database) // 40))
    if nlist < 2:
        raise ValueError("FAISS database is too small for raw-cell neighbor search")
    faiss.omp_set_num_threads(32)
    quantizer = faiss.IndexFlatL2(database.shape[1])
    index = faiss.IndexIVFFlat(quantizer, database.shape[1], nlist, faiss.METRIC_L2)
    index.cp.seed = int(seed)
    rng = np.random.default_rng(seed)
    train_rows = rng.choice(len(database), min(len(database), max(100000, nlist * 40)), replace=False)
    index.train(database[train_rows])
    index.add(database)
    index.nprobe = min(int(spec["nprobe"]), nlist)
    distances = np.empty((len(queries), k), dtype=np.float32)
    neighbors = np.empty((len(queries), k), dtype=np.int64)
    for start in range(0, len(queries), int(spec["search_batch_size"])):
        stop = min(start + int(spec["search_batch_size"]), len(queries))
        d, found = index.search(queries[start:stop], k)
        distances[start:stop], neighbors[start:stop] = d, found
    if (neighbors < 0).any() or not np.isfinite(distances).all():
        raise ValueError("FAISS IVF search did not return complete finite neighbor lists")
    return np.sqrt(np.maximum(distances, 0.0)), neighbors


def sampled_exact_neighbor_recall(source, target, approximate, k, sample_size, seed):
    """Audit an approximate search without attempting an intractable all-pairs exact run."""
    rng = np.random.default_rng(seed)
    selected = rng.choice(len(source), min(len(source), int(sample_size)), replace=False)
    exact = faiss.IndexFlatL2(target.shape[1])
    exact.add(target)
    _, truth = exact.search(source[selected], k)
    recalls = [np.intersect1d(approximate[i], truth[j], assume_unique=False).size / k
               for j, i in enumerate(selected)]
    return float(np.mean(recalls))


def build_cross_species_mnn(latent, n_human, spec, seed):
    """Create expression-only reciprocal FAISS-neighbor links for all raw cells."""
    human, mouse = latent[:n_human], latent[n_human:]
    k = spec["neighbors"]
    radius_k = spec["radius_neighbors"] + 1
    human_radius, _ = faiss_ivf_neighbors(human, human, radius_k, spec, seed + 1)
    mouse_radius, _ = faiss_ivf_neighbors(mouse, mouse, radius_k, spec, seed + 2)
    distance, h_to_m = faiss_ivf_neighbors(mouse, human, k, spec, seed + 3)
    _, m_to_h = faiss_ivf_neighbors(human, mouse, k, spec, seed + 4)
    human_radius = human_radius[:, -1]
    mouse_radius = mouse_radius[:, -1]
    forward = sparse.csr_matrix((np.ones(h_to_m.size, np.int8),
                                 (np.repeat(np.arange(len(human)), k), h_to_m.ravel())),
                                shape=(len(human), len(mouse)))
    backward = sparse.csr_matrix((np.ones(m_to_h.size, np.int8),
                                  (m_to_h.ravel(), np.repeat(np.arange(len(mouse)), k))),
                                 shape=forward.shape)
    mutual = forward.multiply(backward).tocoo()
    actual = np.linalg.norm(human[mutual.row] - mouse[mutual.col], axis=1)
    ratio = actual / np.maximum(np.sqrt(human_radius[mutual.row] * mouse_radius[mutual.col]), 1e-8)
    keep = ratio <= spec["radius_gate"]
    rows, cols, values = mutual.row[keep], mutual.col[keep], np.exp(-ratio[keep] ** 2).astype(np.float32)
    n = len(latent)
    cross = sparse.csr_matrix((values, (rows, n_human + cols)), shape=(n, n), dtype=np.float32)
    cross = (cross + cross.T).tocsr()
    cross.sum_duplicates(); cross.eliminate_zeros()
    human_recall = sampled_exact_neighbor_recall(human, mouse, h_to_m, k, spec["exact_audit_cells"], seed + 5)
    mouse_recall = sampled_exact_neighbor_recall(mouse, human, m_to_h, k, spec["exact_audit_cells"], seed + 6)
    if min(human_recall, mouse_recall) < spec["minimum_exact_recall"]:
        raise ValueError(f"FAISS sampled exact neighbor recall is too low: human={human_recall:.3f}, mouse={mouse_recall:.3f}")
    check = {
        "method": "FAISS IVF-Flat approximate KNN; reciprocal links only",
        "nlist": int(spec["nlist"]), "nprobe": int(spec["nprobe"]),
        "sampled_exact_neighbor_recall_human_to_mouse": human_recall,
        "sampled_exact_neighbor_recall_mouse_to_human": mouse_recall,
        "human_to_mouse_neighbor_candidates": int(forward.nnz),
        "mouse_to_human_neighbor_candidates": int(backward.nnz),
        "mutual_candidates": int(mutual.nnz),
        "retained_undirected_edges": int(len(rows)),
        "radius_gate": float(spec["radius_gate"]),
    }
    return cross, check


def aligned_initialization(human_xy, mouse_xy, cross, n_human, seed):
    """Align the old all-cell layouts only as an initialization, never as output."""
    def standardize(x):
        center = np.median(x, axis=0)
        scale = np.maximum(np.sqrt(np.mean((x - center) ** 2)), 1e-6)
        return (x - center) / scale
    human = standardize(human_xy)
    mouse = standardize(mouse_xy)
    part = cross[:n_human, n_human:].tocoo()
    if part.nnz < 100:
        raise ValueError("Insufficient cross-species MNN anchors for layout initialization")
    rng = np.random.default_rng(seed)
    take = rng.choice(part.nnz, min(100000, part.nnz), replace=False)
    rotation, _ = orthogonal_procrustes(mouse[part.col[take]], human[part.row[take]])
    initial = np.vstack([human, mouse @ rotation]).astype(np.float32)
    if not np.isfinite(initial).all():
        raise ValueError("Joint UMAP initialization is non-finite")
    return initial


def route_colors(routes):
    colors = plt.get_cmap("tab20").colors
    return {route: colors[i % len(colors)] for i, route in enumerate(sorted(routes))}


def median_point(nodes, mask):
    point = nodes.loc[mask, ["UMAP1", "UMAP2"]].to_numpy()
    if not len(point):
        return None
    return np.median(point, axis=0), int(len(point))


def direct_human_arrows(nodes, cfg):
    edges = pd.read_csv(cfg["human"]["display_edges"], sep="\t", dtype=str)
    plan = pd.read_csv(cfg["human"]["temporal_plan"], sep="\t", dtype=str).fillna("")
    rows = []
    human = nodes.species_key.eq("human")
    for edge in edges.itertuples(index=False):
        candidates = plan.loc[(plan.route_id.eq(edge.route_id)) & (plan.time_from.eq(edge.source_time)) &
                              (plan.time_to.eq(edge.target_time))]
        if len(candidates) != 1:
            raise ValueError(f"Human direct display edge does not resolve uniquely: {edge}")
        pair = candidates.iloc[0]
        source_libraries = json.loads(pair.source_libraries)
        target_libraries = json.loads(pair.target_libraries)
        source = human & nodes.library_key.isin(source_libraries) & nodes.experimental_time.eq(edge.source_time) & nodes.condition.eq(edge.source_condition)
        target = human & nodes.library_key.isin(target_libraries) & nodes.experimental_time.eq(edge.target_time) & nodes.condition.eq(edge.target_condition)
        a, b = median_point(nodes, source), median_point(nodes, target)
        if a is None or b is None:
            raise ValueError(f"Human direct display endpoint has no raw cells: {edge}")
        rows.append({"species": "human", "route": edge.route_id, "source": f"{edge.route_id}|{edge.source_time}|{edge.source_condition}",
                     "target": f"{edge.route_id}|{edge.target_time}|{edge.target_condition}", "source_x": a[0][0], "source_y": a[0][1],
                     "target_x": b[0][0], "target_y": b[0][1], "source_cells": a[1], "target_cells": b[1],
                     "evidence": "frozen raw-cell human OT time-pair display"})
    return pd.DataFrame(rows)


def direct_mouse_arrows(nodes, cfg):
    edges = pd.read_csv(cfg["mouse"]["display_edges"], sep="\t")
    rows = []
    mouse = nodes.species_key.eq("mouse")
    for edge in edges.itertuples(index=False):
        source_route, source_state = edge.source.split("|", 1)
        target_route, target_state = edge.target.split("|", 1)
        # ``route`` is the frozen OT display grouping.  Most groups stay within
        # one dataset, but gse171127_late_pooled legitimately connects its D6
        # and D12 source datasets, so endpoint dataset IDs must be respected.
        source = mouse & nodes.dataset.eq(source_route) & nodes.local_state_id.eq(source_state) & np.isclose(nodes.day.astype(float), float(edge.day_from))
        target = mouse & nodes.dataset.eq(target_route) & nodes.local_state_id.eq(target_state) & np.isclose(nodes.day.astype(float), float(edge.day_to))
        a, b = median_point(nodes, source), median_point(nodes, target)
        if a is None or b is None:
            raise ValueError(f"Mouse direct display endpoint has no raw cells: {edge}")
        rows.append({"species": "mouse", "route": edge.route, "source": edge.source, "target": edge.target,
                     "source_dataset": source_route, "target_dataset": target_route,
                     "source_x": a[0][0], "source_y": a[0][1], "target_x": b[0][0], "target_y": b[0][1],
                     "source_cells": a[1], "target_cells": b[1], "mass": float(edge.mass),
                     "evidence": "frozen raw-cell mouse OT trajectory display"})
    return pd.DataFrame(rows)


def inspect_frozen_raw_cell_graph(path, expected_nodes, label):
    """Verify that a frozen constraint graph has exactly one node per raw cell."""
    graph = sparse.load_npz(path).tocsr()
    if graph.shape != (expected_nodes, expected_nodes):
        raise ValueError(f"{label} frozen time graph shape is {graph.shape}, not {expected_nodes} raw cells")
    if not graph.nnz or not np.isfinite(graph.data).all() or (graph - graph.T).nnz:
        raise ValueError(f"{label} frozen raw-cell time graph is invalid")
    return {"shape": list(graph.shape), "nnz": int(graph.nnz), "symmetric": True}


def render_raw_cell_map(nodes, arrows, output):
    output.mkdir(parents=True, exist_ok=True)
    xy = nodes[["UMAP1", "UMAP2"]].to_numpy()
    species = nodes.species_key.to_numpy()
    palette = route_colors(arrows.route.unique())
    fig, ax = plt.subplots(figsize=(16, 12))
    fig.subplots_adjust(left=.03, right=.73, bottom=.06, top=.93)
    rng = np.random.default_rng(20260910)
    order = rng.permutation(len(nodes))
    for key, color, marker in (("human", "#B7476A", "o"), ("mouse", "#23768A", "s")):
        chosen = order[species[order] == key]
        ax.scatter(*xy[chosen].T, s=.22, color=color, marker=marker, alpha=.18, linewidths=0,
                   rasterized=True, zorder=1)
    # A white under-stroke isolates each direct raw-cell OT backbone from the dense point cloud.
    for edge in arrows.itertuples(index=False):
        start, stop = (edge.source_x, edge.source_y), (edge.target_x, edge.target_y)
        for color, width, alpha, zorder in (("white", 4.8, .98, 4), (palette[edge.route], 2.35, 1., 5)):
            ax.add_patch(FancyArrowPatch(start, stop, arrowstyle="-|>", mutation_scale=14,
                                         connectionstyle="arc3,rad=0.08", linewidth=width, color=color,
                                         alpha=alpha, shrinkA=2, shrinkB=2, zorder=zorder))
    ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_title("All-cell human–mouse coembedding with frozen within-study OT trajectories", fontsize=16)
    route_handles = [Line2D([], [], color=palette[route], lw=2.4, label=route) for route in sorted(palette)]
    legend = ax.legend(handles=route_handles, title="Direct OT trajectory route", frameon=False,
                       loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=7.5, title_fontsize=8.5)
    ax.add_artist(legend)
    ax.legend(handles=[Line2D([], [], marker="o", color="#B7476A", linestyle="", label="Human raw cell"),
                       Line2D([], [], marker="s", color="#23768A", linestyle="", label="Mouse raw cell")],
              frameon=False, loc="lower left", bbox_to_anchor=(1.01, 0), fontsize=8.5)
    fig.text(.38, .015,
             "All 551,201 points are raw primary cells. Thick coloured paths are frozen raw-cell OT summaries; "
             "they remain within species and never define a common clock or cross-species lineage.",
             ha="center", fontsize=8.4, color="#4D5660")
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output / f"human_mouse_raw_cell_integrated_trajectories.{suffix}", dpi=250)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--reuse-latent", type=Path,
                        help="A completed raw-cell latent with the same frozen cell ordering; skips retraining.")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    output = args.output
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output}")
    if not torch.cuda.is_available():
        raise RuntimeError("A CUDA GPU is required for the raw-cell species-balanced scVI model")
    output.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:32])
    project = args.config.resolve().parent
    features = read_one_to_one_features(project, cfg)
    human_obs, human_counts, human_initial = load_human_raw_cells(cfg, features)
    mouse_obs, mouse_counts, mouse_initial = load_mouse_raw_cells(cfg, features)
    nodes = pd.concat([human_obs, mouse_obs], ignore_index=True, sort=False)
    if not nodes.cell_id.is_unique or len(nodes) != sum(cfg["expected"].values()):
        raise ValueError("Raw-cell node registry does not close")
    preflight = {
        "status": "PASS",
        "analysis_unit": "raw primary cell",
        "human_cells": int(len(human_obs)), "mouse_cells": int(len(mouse_obs)),
        "one_to_one_ortholog_features": int(len(features)),
        "human_raw_counts_shape": list(human_counts.shape), "mouse_raw_counts_shape": list(mouse_counts.shape),
        "metacells_used": False,
    }
    if args.preflight_only:
        preflight["human_frozen_raw_cell_time_graph"] = inspect_frozen_raw_cell_graph(
            cfg["human"]["time_graph"], len(human_obs), "Human"
        )
        preflight["mouse_frozen_raw_cell_time_graph"] = inspect_frozen_raw_cell_graph(
            cfg["mouse"]["time_graph"], len(mouse_obs), "Mouse"
        )
    (output / "preflight.json").write_text(json.dumps(preflight, indent=2) + "\n")
    if args.preflight_only:
        print(json.dumps(preflight, indent=2))
        return

    if args.reuse_latent:
        latent = np.load(args.reuse_latent).astype(np.float32, copy=False)
        if latent.shape != (len(nodes), cfg["scvi"]["n_latent"]) or not np.isfinite(latent).all():
            raise ValueError("Reused raw-cell latent does not match the current frozen cell registry")
        np.save(output / "species_balanced_raw_cell_scvi_latent.npy", latent)
        features.to_csv(output / "one_to_one_raw_cell_features.tsv", sep="\t", index=False)
        (output / "latent_provenance.json").write_text(json.dumps({
            "source": str(args.reuse_latent), "analysis_unit": "raw primary cell",
            "cell_order": "human frozen order followed by mouse frozen order"
        }, indent=2) + "\n")
    else:
        counts = sparse.vstack([human_counts, mouse_counts], format="csr")
        latent = fit_species_scvi(counts, nodes, features, cfg["scvi"], cfg["seed"], output)
        del counts
    del human_counts, mouse_counts
    gc.collect()
    cross, cross_audit = build_cross_species_mnn(latent, len(human_obs), cfg["cross_species_mnn"], cfg["seed"])
    human_time_graph = sparse.load_npz(cfg["human"]["time_graph"]).tocsr()
    mouse_time_graph = sparse.load_npz(cfg["mouse"]["time_graph"]).tocsr()
    if human_time_graph.shape != (len(human_obs), len(human_obs)):
        raise ValueError("Human frozen time graph does not match the raw-cell registry")
    if mouse_time_graph.shape != (len(mouse_obs), len(mouse_obs)):
        raise ValueError("Mouse frozen time graph does not match the raw-cell registry")
    within = sparse.block_diag((human_time_graph, mouse_time_graph), format="csr")
    if within.shape != cross.shape or not np.isfinite(within.data).all() or (within - within.T).nnz:
        raise ValueError("Frozen all-cell within-species time graph is invalid")
    mix = row_normalize(cfg["layout"]["within_time_graph_weight"] * within +
                        cfg["layout"]["cross_species_mnn_weight"] * cross)
    graph = ((mix + mix.T) * .5).tocsr()
    graph.data /= graph.data.max(); graph.eliminate_zeros()
    if not np.isfinite(graph.data).all() or (graph - graph.T).nnz:
        raise ValueError("All-cell integrated graph is invalid")
    initial = aligned_initialization(human_initial, mouse_initial, cross, len(human_obs), cfg["seed"])
    a, b = find_ab_params(cfg["layout"]["spread"], cfg["layout"]["min_dist"])
    xy, _ = simplicial_set_embedding(
        latent, graph.copy(), n_components=2, initial_alpha=1., a=a, b=b, gamma=1., negative_sample_rate=5,
        n_epochs=cfg["layout"]["epochs"], init=initial, random_state=np.random.RandomState(cfg["seed"]),
        metric="euclidean", metric_kwds={}, densmap=False, densmap_kwds={}, output_dens=False, parallel=True, verbose=True,
    )
    if xy.shape != (len(nodes), 2) or not np.isfinite(xy).all():
        raise ValueError("All-cell integrated UMAP is invalid")
    nodes[["UMAP1", "UMAP2"]] = xy
    arrows = pd.concat([direct_human_arrows(nodes, cfg), direct_mouse_arrows(nodes, cfg)], ignore_index=True)
    if not len(arrows):
        raise ValueError("No direct raw-cell trajectories were available for display")
    nodes.to_csv(output / "raw_cell_nodes.tsv.gz", sep="\t", index=False)
    arrows.to_csv(output / "raw_cell_trajectory_arrows.tsv", sep="\t", index=False)
    features.to_csv(output / "one_to_one_raw_cell_features.tsv.gz", sep="\t", index=False)
    sparse.save_npz(output / "cross_species_expression_mnn.npz", cross)
    sparse.save_npz(output / "integrated_time_graph.npz", graph)
    render_raw_cell_map(nodes, arrows, output / "figures")
    summary = {
        "status": "COMPUTATION_COMPLETE_VALIDATION_PENDING",
        **preflight,
        "scvi": {**cfg["scvi"], "labels_time_apa_used": False},
        "within_time_graph": "frozen full-cell human and mouse time graphs; no metacells",
        "within_time_graph_weight": cfg["layout"]["within_time_graph_weight"],
        "cross_species_mnn_weight": cfg["layout"]["cross_species_mnn_weight"],
        "cross_species_mnn": cross_audit,
        "integrated_graph_components": int(connected_components(graph, directed=False)[0]),
        "direct_raw_cell_trajectory_arrows": int(len(arrows)),
        "direct_raw_cell_trajectory_routes": int(arrows.route.nunique()),
        "cross_species_ot_or_common_clock_created": False,
        "cross_species_lineage_or_identity_validated": False,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
