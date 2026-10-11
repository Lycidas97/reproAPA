#!/usr/bin/env python3
"""Build a metacell-level human/mouse coembedding with frozen within-species OT.

The graph has three distinct layers: local expression edges within each study,
expression-only MNN edges between studies (including human--mouse), and the
already accepted within-species sampled OT kernels aggregated to metacells.
No state, time, or APA label is used to build the expression layers.
"""

import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from umap.umap_ import find_ab_params, simplicial_set_embedding


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_rows(matrix):
    total = np.asarray(matrix.sum(axis=1)).ravel()
    return (sparse.diags(1.0 / np.maximum(total, 1e-30)) @ matrix).tocsr()


def log_normalize(matrix):
    matrix = sparse.csr_matrix(matrix, dtype=np.float32)
    total = np.asarray(matrix.sum(axis=1)).ravel()
    scaled = sparse.diags(1e4 / np.maximum(total, 1.0)) @ matrix
    scaled.data = np.log1p(scaled.data)
    return scaled.astype(np.float32)


def read_orthologs(path):
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    eligible = table.documented_one_to_one.astype(str).str.lower().eq("true")
    table = table.loc[eligible & table.mouse_gene_id.ne("") & table.human_gene_id.ne("")].copy()
    if table.mouse_gene_id.duplicated().any() or table.human_gene_id.duplicated().any():
        raise ValueError("One-to-one ortholog registry is not unique after eligibility filtering")
    return table.sort_values(["mouse_gene_id", "human_gene_id"]).reset_index(drop=True)


def restrict_to_input_orthologs(orthologs, human_path, mouse_root):
    """Keep only registry pairs represented in both frozen all-gene matrices."""
    human = ad.read_h5ad(human_path, backed="r")
    human_ids = set(human.var_names.astype(str))
    example_mouse = next(Path(mouse_root).glob("*/all_gene_counts.h5ad"), None)
    if example_mouse is None:
        raise FileNotFoundError(f"No mouse metacell all-gene H5AD under {mouse_root}")
    mouse = ad.read_h5ad(example_mouse, backed="r")
    mouse_ids = set(mouse.var.gene_ids.astype(str))
    retained = orthologs.loc[
        orthologs.human_gene_id.isin(human_ids) & orthologs.mouse_gene_id.isin(mouse_ids)
    ].copy().reset_index(drop=True)
    if len(retained) < 1000:
        raise ValueError(f"Only {len(retained)} one-to-one orthologs are present in both inputs")
    return retained


def build_human_matrix(path, orthologs):
    data = ad.read_h5ad(path)
    if data.n_obs != 4332:
        raise ValueError(f"Unexpected human metacell count: {data.n_obs}")
    index = pd.Index(data.var_names.astype(str))
    pos = index.get_indexer(orthologs.human_gene_id)
    if (pos < 0).any():
        raise ValueError("Input-intersected orthologs unexpectedly missing from human matrix")
    obs = data.obs.copy().reset_index(names="metacell_id")
    obs["species_key"] = "human"
    obs["dataset"] = obs.dataset_id.astype(str)
    obs["route"] = obs.dataset.astype(str)
    obs["coarse_state"] = obs.coarse_state.astype(str)
    obs["experimental_time"] = obs.experimental_time.astype(str)
    counts = sparse.csr_matrix(data.X[:, pos], dtype=np.float32)
    if counts.shape != (len(obs), len(orthologs)):
        raise ValueError("Human ortholog matrix shape mismatch")
    return obs, counts


def build_mouse_matrix(registry_path, gex_root, orthologs, mouse_state_map):
    registry = pd.read_csv(registry_path, sep="\t", dtype={"metacell_id": str}).set_index("metacell_id")
    if len(registry) != 6772 or not registry.index.is_unique:
        raise ValueError("Frozen mouse metacell registry is incomplete")
    state_map = pd.read_csv(mouse_state_map, sep="\t", dtype=str)
    state_map = state_map[["dataset", "local_state_id", "coarse_state"]].drop_duplicates()
    if state_map.duplicated(["dataset", "local_state_id"]).any():
        raise ValueError("Mouse local-to-coarse mapping is ambiguous")
    parts, metadata = [], []
    reference_ids = None
    for library in sorted(registry.library_key.unique()):
        path = Path(gex_root) / library / "all_gene_counts.h5ad"
        if not path.exists():
            raise FileNotFoundError(path)
        data = ad.read_h5ad(path)
        if not set(data.obs_names.astype(str)).issubset(registry.index):
            raise ValueError(f"Mouse GEX rows do not belong to frozen registry: {library}")
        ids = data.var.gene_ids.astype(str).to_numpy()
        if reference_ids is None:
            reference_ids = ids
            positions = pd.Index(ids).get_indexer(orthologs.mouse_gene_id)
            if (positions < 0).any():
                raise ValueError("Input-intersected orthologs unexpectedly missing from mouse matrix")
        elif not np.array_equal(ids, reference_ids):
            raise ValueError(f"Mouse feature order differs in {library}")
        rows = registry.loc[data.obs_names.astype(str)].copy()
        if not rows.library_key.eq(library).all():
            raise ValueError(f"Mouse GEX library identity mismatch: {library}")
        parts.append(sparse.csr_matrix(data.X[:, positions], dtype=np.float32))
        metadata.append(rows)
    counts = sparse.vstack(parts, format="csr")
    obs = pd.concat(metadata).reset_index(names="metacell_id")
    if len(obs) != 6772 or not obs.metacell_id.is_unique:
        raise ValueError("Mouse GEX aggregation did not close frozen metacells")
    obs = obs.merge(state_map, on=["dataset", "local_state_id"], how="left", validate="many_to_one")
    if obs.coarse_state.isna().any():
        missing = obs.loc[obs.coarse_state.isna(), ["dataset", "local_state_id"]].drop_duplicates()
        raise ValueError(f"Mouse coarse state missing: {missing.to_dict('records')[:5]}")
    obs["species_key"] = "mouse"
    obs["route"] = obs.dataset.astype(str)
    obs["experimental_time"] = obs.experimental_time.astype(str)
    if counts.shape != (len(obs), len(orthologs)):
        raise ValueError("Mouse ortholog matrix shape mismatch")
    return obs, counts


def map_human_cells(time_h5ad, cell_to_metacell, metacell_index):
    source = ad.read_h5ad(time_h5ad, backed="r")
    cells = source.obs.cell_id.astype(str)
    mapping = pd.read_csv(cell_to_metacell, sep="\t", usecols=["cell_id", "metacell_id"], dtype=str)
    mapping = mapping.set_index("cell_id")
    if not mapping.index.is_unique or len(cells) != 215114:
        raise ValueError("Human frozen cell-to-metacell mapping is invalid")
    assigned = mapping.reindex(cells).metacell_id
    if assigned.isna().any():
        raise ValueError("Human time-layout cells are absent from the metacell mapping")
    code = metacell_index.get_indexer(assigned)
    if (code < 0).any():
        raise ValueError("Human mapping refers to non-existent metacells")
    return code


def map_mouse_cells(coordinate_path, cell_to_metacell_root, metacell_index):
    cells = pd.read_csv(coordinate_path, sep="\t", usecols=["cell_id", "library_key"], dtype=str)
    if len(cells) != 336087 or not cells.cell_id.is_unique:
        raise ValueError("Mouse time-layout coordinate registry is invalid")
    mappings = []
    for library in sorted(cells.library_key.unique()):
        path = Path(cell_to_metacell_root) / library / "cell_to_metacell.tsv.gz"
        if not path.exists() and library == "GSM3195672":
            path = Path("/root/nfsdata/apadetect/analysis_projects/mouse_time_aware_integration_pilot_v1/results_metacell_pilot_gsm3195672_v1/target_50/cell_to_metacell.tsv.gz")
        if not path.exists():
            raise FileNotFoundError(path)
        part = pd.read_csv(path, sep="\t", usecols=["cell_id", "metacell_id", "library_key"], dtype=str)
        if not part.library_key.eq(library).all():
            raise ValueError(f"Mouse membership source mismatch: {library}")
        mappings.append(part[["cell_id", "metacell_id"]])
    mapping = pd.concat(mappings, ignore_index=True).set_index("cell_id")
    if not mapping.index.is_unique:
        raise ValueError("Mouse cell-to-metacell mapping has duplicate cell IDs")
    assigned = mapping.reindex(cells.cell_id).metacell_id
    if assigned.isna().any():
        raise ValueError("Mouse time-layout cells are absent from the metacell mapping")
    code = metacell_index.get_indexer(assigned)
    if (code < 0).any():
        raise ValueError("Mouse mapping refers to non-existent metacells")
    return code


def aggregate_cell_kernel(path, cell_to_metacell, n_metacells):
    matrix = sparse.load_npz(path).tocoo()
    if matrix.shape[0] != len(cell_to_metacell) or matrix.shape[1] != len(cell_to_metacell):
        raise ValueError(f"Cell kernel shape does not match frozen cell order: {path}")
    result = sparse.csr_matrix(
        (matrix.data.astype(np.float32), (cell_to_metacell[matrix.row], cell_to_metacell[matrix.col])),
        shape=(n_metacells, n_metacells), dtype=np.float32)
    result.sum_duplicates()
    result.eliminate_zeros()
    return result


def local_kernel(z, groups, neighbors):
    n = len(z)
    rows, cols, values = [], [], []
    radius = np.zeros(n, dtype=np.float32)
    for group in sorted(pd.unique(groups)):
        index = np.flatnonzero(groups == group)
        k = min(neighbors, len(index) - 1)
        if k < 1:
            raise ValueError(f"Local graph group has fewer than two metacells: {group}")
        distances, found = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", n_jobs=-1).fit(z[index]).kneighbors(z[index])
        found, distances = found[:, 1:], distances[:, 1:]
        radius[index] = distances[:, -1]
        scale = max(float(np.median(distances)), 1e-6)
        rows.append(np.repeat(index, k))
        cols.append(index[found.ravel()])
        values.append(np.exp(-(distances / scale) ** 2).ravel())
    kernel = sparse.csr_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    return normalize_rows(kernel.maximum(kernel.T)).astype(np.float32), radius


def mnn_kernel(z, studies, species, radius, neighbors, gate, seed):
    n = len(z)
    indices = {study: np.flatnonzero(studies == study) for study in sorted(pd.unique(studies))}
    trees = {study: cKDTree(z[index]) for study, index in indices.items()}
    rows, cols, values, checks = [], [], [], []
    rng = np.random.default_rng(seed)
    for left, right in itertools.combinations(indices, 2):
        a, b = indices[left], indices[right]
        k = min(neighbors, len(a), len(b))
        distance, forward = trees[right].query(z[a], k=k, workers=-1)
        _, backward = trees[left].query(z[b], k=k, workers=-1)
        if k == 1:
            distance, forward = distance[:, None], forward[:, None]
            backward = backward[:, None]
        recalls = []
        for source, target, approximate in ((a, b, forward), (b, a, backward)):
            probe = rng.choice(len(source), min(100, len(source)), replace=False)
            exact = NearestNeighbors(n_neighbors=k, algorithm="brute", n_jobs=-1).fit(z[target]).kneighbors(z[source[probe]], return_distance=False)
            recalls.append(float(np.mean([len(set(x) & set(y)) / k for x, y in zip(exact, approximate[probe])])) )
        if min(recalls) < 0.99:
            raise AssertionError(f"Exact MNN check failed for {left}/{right}: {recalls}")
        fwd = sparse.csr_matrix((np.ones(forward.size, np.int8), (np.repeat(np.arange(len(a)), k), forward.ravel())), shape=(len(a), len(b)))
        back = sparse.csr_matrix((np.ones(backward.size, np.int8), (backward.ravel(), np.repeat(np.arange(len(b)), k))), shape=fwd.shape)
        mutual = fwd.multiply(back).tocoo()
        row, col = a[mutual.row], b[mutual.col]
        actual = np.linalg.norm(z[row] - z[col], axis=1)
        ratio = actual / np.maximum(np.sqrt(radius[row] * radius[col]), 1e-6)
        valid = ratio <= gate
        rows.append(row[valid])
        cols.append(col[valid])
        values.append(np.exp(-ratio[valid] ** 2))
        checks.append({
            "study_a": left, "study_b": right,
            "species_a": species[a[0]], "species_b": species[b[0]],
            "mutual_candidates": int(len(row)), "retained": int(valid.sum()),
            "recall_a_to_b": recalls[0], "recall_b_to_a": recalls[1]
        })
    kernel = sparse.csr_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n), dtype=np.float32)
    return normalize_rows(kernel.maximum(kernel.T)).astype(np.float32), pd.DataFrame(checks)


def semantic_group(species, state):
    shared = {
        "S": "Source-like", "R": "Exit/remodelling", "E": "Epithelial/mixed",
        "M": "ECM/mesenchymal retained"
    }
    if state in shared:
        return shared[state]
    if (species == "human" and state == "P") or (species == "mouse" and state == "T"):
        return "Pluripotency programme"
    return "Species-specific / unresolved"


def arrows_from_forward_kernel(kernel, nodes, species, minimum=0.15):
    index = np.flatnonzero(nodes.species_key.eq(species))
    part = kernel[index][:, index].tocsr()
    key = (nodes.loc[index, "route"].astype(str) + " | " + nodes.loc[index, "experimental_time"].astype(str)).to_numpy()
    code, names = pd.factorize(key, sort=True)
    membership = sparse.csr_matrix((np.ones(len(index)), (code, np.arange(len(index)))), shape=(len(names), len(index)))
    flow = (membership @ part @ membership.T).tocsr()
    positions = pd.DataFrame({"group": names, "x": 0.0, "y": 0.0, "n_metacells": 0})
    xy = nodes.loc[index, ["UMAP1", "UMAP2"]].to_numpy()
    for group, name in enumerate(names):
        point = xy[code == group]
        positions.loc[group, ["x", "y", "n_metacells"]] = [np.median(point[:, 0]), np.median(point[:, 1]), len(point)]
    rows = []
    for source in range(flow.shape[0]):
        start, stop = flow.indptr[source], flow.indptr[source + 1]
        targets, values = flow.indices[start:stop], flow.data[start:stop]
        if not len(targets) or values.sum() <= 0:
            continue
        candidates = [(target, value / values.sum()) for target, value in zip(targets, values) if target != source]
        if not candidates:
            continue
        target, fraction = max(candidates, key=lambda item: item[1])
        if fraction >= minimum:
            rows.append({"species": species, "source_group": names[source], "target_group": names[target], "source_fraction": fraction})
    arrows = pd.DataFrame(rows)
    if len(arrows):
        arrows = arrows.merge(positions.add_prefix("source_").rename(columns={"source_group": "source_group"}), on="source_group", how="left")
        arrows = arrows.merge(positions.add_prefix("target_").rename(columns={"target_group": "target_group"}), on="target_group", how="left")
    return arrows


def mouse_forward_kernel(path, metacell_index, n):
    table = pd.read_csv(path, sep="\t")
    row = metacell_index.get_indexer(table.source_metacell.astype(str))
    col = metacell_index.get_indexer(table.target_metacell.astype(str))
    if (row < 0).any() or (col < 0).any():
        raise ValueError("Frozen mouse top-3 OT file refers to unknown metacells")
    result = sparse.csr_matrix((table.conditional.to_numpy(np.float32), (row, col)), shape=(n, n))
    result.sum_duplicates()
    return result


def fit_layout(z, local, temporal, cross, weights, graph_cfg, seed):
    mixture = normalize_rows(weights[0] * local + weights[1] * temporal + weights[2] * cross)
    graph = ((mixture + mixture.T) * 0.5).tocsr()
    graph.data /= graph.data.max()
    graph.eliminate_zeros()
    if not np.isfinite(graph.data).all() or (graph - graph.T).nnz:
        raise ValueError("Common graph is not finite and symmetric")
    a, b = find_ab_params(graph_cfg["spread"], graph_cfg["min_dist"])
    xy, _ = simplicial_set_embedding(
        z, graph.copy(), n_components=2, initial_alpha=1.0, a=a, b=b, gamma=1.0,
        negative_sample_rate=5, n_epochs=graph_cfg["umap_epochs"], init="spectral",
        random_state=np.random.RandomState(seed), metric="euclidean", metric_kwds={},
        densmap=False, densmap_kwds={}, output_dens=False, parallel=True, verbose=False)
    if not np.isfinite(xy).all() or xy.shape != (len(z), 2):
        raise ValueError("UMAP did not return finite common coordinates")
    return xy, graph


def render(nodes, arrows, cross, out):
    colors = {
        "Source-like": "#C58B24", "Exit/remodelling": "#4C78A8",
        "Epithelial/mixed": "#67B6D0", "Pluripotency programme": "#208D79",
        "ECM/mesenchymal retained": "#826443", "Species-specific / unresolved": "#929AA4"
    }
    fig, axes = plt.subplots(1, 2, figsize=(18, 8.4), constrained_layout=True)
    rng = np.random.default_rng(20260910)
    xy = nodes[["UMAP1", "UMAP2"]].to_numpy()
    cross_coo = cross.tocoo()
    cross_mask = (nodes.species_key.to_numpy()[cross_coo.row] != nodes.species_key.to_numpy()[cross_coo.col]) & (cross_coo.row < cross_coo.col)
    cross_rows, cross_cols = cross_coo.row[cross_mask], cross_coo.col[cross_mask]
    if len(cross_rows):
        take = rng.choice(len(cross_rows), min(2000, len(cross_rows)), replace=False)
        segments = np.stack([xy[cross_rows[take]], xy[cross_cols[take]]], axis=1)
        for ax in axes:
            ax.add_collection(LineCollection(segments, colors="#6E5B9A", linewidths=0.20, alpha=0.055, rasterized=True))
    species_style = {"human": ("#B7476A", "o", "Human"), "mouse": ("#23768A", "s", "Mouse")}
    for species, (color, marker, label) in species_style.items():
        selected = nodes.species_key.eq(species).to_numpy()
        axes[0].scatter(xy[selected, 0], xy[selected, 1], s=7, c=color, marker=marker, alpha=.70,
                        linewidths=0, rasterized=True, label=f"{label} ({selected.sum():,} metacells)")
    for group, color in colors.items():
        selected = nodes.semantic_group.eq(group).to_numpy()
        if not selected.any():
            continue
        for species, (_, marker, _) in species_style.items():
            visible = selected & nodes.species_key.eq(species).to_numpy()
            if visible.any():
                axes[1].scatter(xy[visible, 0], xy[visible, 1], s=7, c=color, marker=marker, alpha=.72,
                                linewidths=0, rasterized=True)
    arrow_colors = {"human": "#8F274F", "mouse": "#14546B"}
    for ax in axes:
        for row in arrows.itertuples(index=False):
            ax.annotate("", xy=(row.target_x, row.target_y), xytext=(row.source_x, row.source_y),
                        arrowprops={"arrowstyle": "->", "color": arrow_colors[row.species], "lw": .55,
                                    "alpha": .65, "shrinkA": 3, "shrinkB": 3})
        ax.set_aspect("equal")
        ax.set_axis_off()
    axes[0].set_title("Species and study-internal time trajectories", fontsize=14)
    axes[1].set_title("Broad state semantics; circle = human, square = mouse", fontsize=14)
    axes[0].legend(frameon=False, loc="best", fontsize=10)
    state_handles = [Line2D([], [], marker="o", color=color, linestyle="", markersize=5, label=label) for label, color in colors.items()]
    axes[1].legend(handles=state_handles, frameon=False, loc="best", fontsize=8)
    fig.suptitle("Human–mouse reprogramming coembedding with frozen within-species time constraints", fontsize=17)
    fig.text(.5, .018,
             "Points are frozen study-local metacells. Purple links are a sparse display of expression-only cross-species MNN. "
             "Arrows summarize within-species OT only; they do not define a common clock or cross-species lineage.",
             ha="center", fontsize=8.7, color="#4D5660")
    out.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(out / f"human_mouse_time_aware_metacell_umap.{suffix}", dpi=250)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    out = args.output
    if (out / "summary.json").exists():
        raise FileExistsError(f"Refusing to overwrite completed output: {out}")
    out.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:32])

    project = args.config.resolve().parent
    registered_orthologs = read_orthologs(project / cfg["orthologs"])
    orthologs = restrict_to_input_orthologs(registered_orthologs, cfg["human"]["metacell_gex"], cfg["mouse"]["metacell_gex_root"])
    human_obs, human_counts = build_human_matrix(cfg["human"]["metacell_gex"], orthologs)
    mouse_obs, mouse_counts = build_mouse_matrix(cfg["mouse"]["metacell_registry"], cfg["mouse"]["metacell_gex_root"], orthologs, project / cfg["mouse_state_map"])
    nodes = pd.concat([human_obs, mouse_obs], ignore_index=True, sort=False)
    if len(nodes) != cfg["expected"]["human_metacells"] + cfg["expected"]["mouse_metacells"] or not nodes.metacell_id.is_unique:
        raise ValueError("Combined metacell registry does not close")
    nodes["semantic_group"] = [semantic_group(species, state) for species, state in zip(nodes.species_key, nodes.coarse_state)]
    nodes["study_group"] = nodes.species_key.astype(str) + "::" + nodes.route.astype(str)
    metacell_index = pd.Index(nodes.metacell_id.astype(str))

    expression = sparse.vstack([log_normalize(human_counts), log_normalize(mouse_counts)], format="csr")
    svd = TruncatedSVD(n_components=cfg["graph"]["pca_components"], random_state=cfg["seed"])
    latent = svd.fit_transform(expression).astype(np.float32)
    # Per-species centring/scaling removes global measurement scale, but not labels or time information.
    for species in ("human", "mouse"):
        index = np.flatnonzero(nodes.species_key.eq(species))
        mean, scale = latent[index].mean(axis=0), np.maximum(latent[index].std(axis=0), 1e-6)
        latent[index] = (latent[index] - mean) / scale
    np.save(out / "ortholog_expression_pca.npy", latent)
    orthologs.to_csv(out / "one_to_one_orthologs_used.tsv.gz", sep="\t", index=False)

    local, radius = local_kernel(latent, nodes.study_group.to_numpy(), cfg["graph"]["local_neighbors"])
    cross, mnn_checks = mnn_kernel(latent, nodes.study_group.to_numpy(), nodes.species_key.to_numpy(), radius,
                                   cfg["graph"]["cross_neighbors"], cfg["graph"]["cross_radius_gate"], cfg["seed"])
    human_codes = map_human_cells(cfg["human"]["time_h5ad"], cfg["human"]["cell_to_metacell"], metacell_index)
    mouse_codes = map_mouse_cells(cfg["mouse"]["cell_coordinates"], cfg["mouse"]["cell_to_metacell_root"], metacell_index)
    human_time = aggregate_cell_kernel(cfg["human"]["time_kernel"], human_codes, len(nodes))
    mouse_time = aggregate_cell_kernel(cfg["mouse"]["time_kernel"], mouse_codes, len(nodes))
    temporal = normalize_rows(human_time + mouse_time).astype(np.float32)
    sparse.save_npz(out / "local_expression_kernel.npz", local)
    sparse.save_npz(out / "cross_study_expression_mnn_kernel.npz", cross)
    sparse.save_npz(out / "within_species_time_kernel.npz", temporal)

    time_xy, time_graph = fit_layout(latent, local, temporal, cross, cfg["graph"]["time_weights"], cfg["graph"], cfg["seed"])
    control_xy, control_graph = fit_layout(latent, local, temporal, cross, cfg["graph"]["control_weights"], cfg["graph"], cfg["seed"])
    nodes[["UMAP1", "UMAP2"]] = time_xy
    nodes[["UMAP1_control", "UMAP2_control"]] = control_xy
    sparse.save_npz(out / "time_constrained_graph.npz", time_graph)
    sparse.save_npz(out / "no_time_control_graph.npz", control_graph)

    human_forward = aggregate_cell_kernel(cfg["human"]["forward_time_kernel"], human_codes, len(nodes))
    mouse_forward = mouse_forward_kernel(project / "sources/local/mouse-metacell-ot-top3-edges.tsv.gz", metacell_index, len(nodes))
    arrows = pd.concat([
        arrows_from_forward_kernel(human_forward, nodes, "human"),
        arrows_from_forward_kernel(mouse_forward, nodes, "mouse")
    ], ignore_index=True)
    arrows.to_csv(out / "within_species_trajectory_arrows.tsv", sep="\t", index=False)
    nodes.to_csv(out / "common_metacell_nodes.tsv.gz", sep="\t", index=False)
    mnn_checks.to_csv(out / "cross_study_mnn_checks.tsv", sep="\t", index=False)
    render(nodes, arrows, cross, out / "figures")

    components = int(connected_components(time_graph, directed=False)[0])
    cross_coo = cross.tocoo()
    cross_species_edges = int(((nodes.species_key.to_numpy()[cross_coo.row] != nodes.species_key.to_numpy()[cross_coo.col]) & (cross_coo.row < cross_coo.col)).sum())
    result = {
        "status": "COMPUTATION_COMPLETE_VALIDATION_PENDING",
        "analysis_unit": "frozen within-library metacell",
        "human_metacells": int(nodes.species_key.eq("human").sum()),
        "mouse_metacells": int(nodes.species_key.eq("mouse").sum()),
        "human_primary_cells_represented": cfg["expected"]["human_cells"],
        "mouse_primary_cells_represented": cfg["expected"]["mouse_cells"],
        "one_to_one_orthologs": int(len(orthologs)),
        "one_to_one_orthologs_in_registry": int(len(registered_orthologs)),
        "expression_representation": "library-size-normalized log ortholog counts; TruncatedSVD; per-species PCA centring/scaling",
        "labels_used_to_build_expression_layers": False,
        "time_used_to_build_expression_layers": False,
        "apa_used_to_build_expression_layers": False,
        "time_kernel": "frozen accepted within-species sampled OT kernels aggregated to frozen metacells and symmetrized only for layout",
        "time_weights": cfg["graph"]["time_weights"],
        "control_weights": cfg["graph"]["control_weights"],
        "time_graph_components": components,
        "time_graph_edges": int(time_graph.nnz // 2),
        "cross_species_mnn_undirected_edges": cross_species_edges,
        "display_arrows": int(len(arrows)),
        "cross_species_lineage_or_common_clock_validated": False,
        "input_sha256": {
            "orthologs": sha256(project / cfg["orthologs"]),
            "human_metacell_gex": sha256(cfg["human"]["metacell_gex"]),
            "mouse_metacell_registry": sha256(cfg["mouse"]["metacell_registry"]),
            "human_time_kernel": sha256(cfg["human"]["time_kernel"]),
            "mouse_time_kernel": sha256(cfg["mouse"]["time_kernel"])
        }
    }
    (out / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
