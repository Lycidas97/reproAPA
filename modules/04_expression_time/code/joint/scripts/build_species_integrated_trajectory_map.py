#!/usr/bin/env python3
"""Fit a species-balanced metacell coembedding and overlay frozen study trajectories.

The species covariate is used only by scVI to make a display-oriented shared
expression latent. Study, state, time, and APA are never given to scVI. The
two frozen within-species OT kernels are then added back solely as layout
constraints; no cross-species OT or arrows are created.
"""

import argparse
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
import scanpy as sc
import scvi
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import torch

from build_human_mouse_time_aware_umap import (
    aggregate_cell_kernel,
    arrows_from_forward_kernel,
    build_human_matrix,
    build_mouse_matrix,
    fit_layout,
    local_kernel,
    map_human_cells,
    map_mouse_cells,
    mnn_kernel,
    mouse_forward_kernel,
    normalize_rows,
    read_orthologs,
    restrict_to_input_orthologs,
    semantic_group,
)


def fit_species_balanced_scvi(counts, nodes, orthologs, spec, seed, output):
    """Learn a shared latent using raw one-to-one ortholog counts.

    `species_key` is deliberately the sole batch covariate. It is a display
    alignment choice, not a claim that the species effect is technical noise.
    """
    counts = sparse.csr_matrix(counts, dtype=np.float32)
    if counts.data.size and not np.allclose(counts.data, np.rint(counts.data)):
        raise ValueError("scVI input must retain integer raw counts")
    var = orthologs[["human_gene_id", "mouse_gene_id", "human_gene_symbol", "mouse_gene_symbol"]].copy()
    var.index = pd.Index([f"ortholog_{i:05d}" for i in range(len(var))], name="ortholog_feature")
    data = ad.AnnData(X=counts, obs=nodes[["species_key"]].copy(), var=var)
    data.layers["counts"] = data.X.copy()
    sc.pp.highly_variable_genes(
        data, n_top_genes=spec["n_hvg"], flavor="seurat_v3", layer="counts",
        batch_key=spec["batch_key"], subset=True,
    )
    if data.n_vars != spec["n_hvg"]:
        raise ValueError(f"Unexpected HVG count: {data.n_vars}")
    scvi.settings.seed = seed
    torch.set_float32_matmul_precision("high")
    scvi.model.SCVI.setup_anndata(data, layer="counts", batch_key=spec["batch_key"])
    model = scvi.model.SCVI(data, n_latent=spec["n_latent"], gene_likelihood="nb")
    model.train(
        max_epochs=spec["max_epochs"], early_stopping=spec["early_stopping"],
        accelerator="gpu", devices=1, batch_size=spec["batch_size"],
    )
    latent = model.get_latent_representation().astype(np.float32)
    if latent.shape != (len(nodes), spec["n_latent"]) or not np.isfinite(latent).all():
        raise ValueError("Species-balanced scVI did not return a finite full latent")
    data.var.reset_index().to_csv(output / "scvi_hvg_orthologs.tsv", sep="\t", index=False)
    model.save(output / "species_balanced_scvi_model", overwrite=False)
    np.save(output / "species_balanced_scvi_latent.npy", latent)
    return latent


def route_palette(routes):
    base = plt.get_cmap("tab20").colors
    return {route: base[i] for i, route in enumerate(sorted(routes))}


def draw_arrows(ax, arrows, colors):
    for row in arrows.itertuples(index=False):
        ax.annotate(
            "", xy=(row.target_x, row.target_y), xytext=(row.source_x, row.source_y),
            arrowprops={"arrowstyle": "->", "color": colors[row.route], "lw": 0.75,
                        "alpha": 0.80, "shrinkA": 3, "shrinkB": 3},
            zorder=4,
        )


def render(nodes, arrows, cross, out):
    """Render one shared coordinate system with study-resolved trajectory overlays."""
    out.mkdir(parents=True, exist_ok=True)
    xy = nodes[["UMAP1", "UMAP2"]].to_numpy()
    species = nodes.species_key.to_numpy()
    route_colors = route_palette(nodes.route.astype(str).unique())
    semantic_colors = {
        "Source-like": "#C58B24", "Exit/remodelling": "#4C78A8",
        "Epithelial/mixed": "#67B6D0", "Pluripotency programme": "#208D79",
        "ECM/mesenchymal retained": "#826443", "Species-specific / unresolved": "#929AA4",
    }
    fig, axes = plt.subplots(1, 3, figsize=(25.5, 8.8), constrained_layout=True)
    rng = np.random.default_rng(20260910)
    edge = cross.tocoo()
    mask = (species[edge.row] != species[edge.col]) & (edge.row < edge.col)
    rows, cols = edge.row[mask], edge.col[mask]
    if len(rows):
        chosen = rng.choice(len(rows), min(3000, len(rows)), replace=False)
        segments = np.stack([xy[rows[chosen]], xy[cols[chosen]]], axis=1)
        for ax in axes:
            ax.add_collection(LineCollection(segments, colors="#756297", linewidths=.18,
                                              alpha=.035, rasterized=True, zorder=1))

    style = {"human": ("#B7476A", "o", "Human"), "mouse": ("#23768A", "s", "Mouse")}
    for key, (color, marker, label) in style.items():
        selected = species == key
        axes[0].scatter(*xy[selected].T, s=8, c=color, marker=marker, alpha=.62, linewidths=0,
                        rasterized=True, label=f"{label} ({selected.sum():,} metacells)", zorder=2)
    for route in sorted(route_colors):
        selected = nodes.route.astype(str).eq(route).to_numpy()
        axes[1].scatter(*xy[selected].T, s=7, c=[route_colors[route]], alpha=.58, linewidths=0,
                        rasterized=True, label=route, zorder=2)
    for group, color in semantic_colors.items():
        selected = nodes.semantic_group.eq(group).to_numpy()
        for key, (_, marker, _) in style.items():
            visible = selected & (species == key)
            if visible.any():
                axes[2].scatter(*xy[visible].T, s=7, c=color, marker=marker, alpha=.65, linewidths=0,
                                rasterized=True, zorder=2)

    for ax in axes:
        draw_arrows(ax, arrows, route_colors)
        ax.set_aspect("equal")
        ax.set_axis_off()
    axes[0].set_title("Shared map: species", fontsize=14)
    axes[1].set_title("Study-specific OT trajectories", fontsize=14)
    axes[2].set_title("Broad state annotation (not used for alignment)", fontsize=14)
    axes[0].legend(frameon=False, loc="best", fontsize=10)
    route_handles = [Line2D([], [], color=route_colors[x], lw=2, label=x) for x in sorted(route_colors)]
    axes[1].legend(handles=route_handles, frameon=False, loc="best", fontsize=6.6, ncol=1)
    state_handles = [Line2D([], [], marker="o", color=color, linestyle="", markersize=5, label=name)
                     for name, color in semantic_colors.items()]
    axes[2].legend(handles=state_handles, frameon=False, loc="best", fontsize=7.5)
    fig.suptitle("Species-balanced human–mouse reprogramming map with frozen study-internal trajectories", fontsize=18)
    fig.text(
        .5, .018,
        "Each point is a frozen library-local metacell. scVI aligns species using raw one-to-one ortholog counts; "
        "study, time, state and APA were withheld from scVI. Coloured arrows summarize frozen OT within each study only; "
        "they do not create a common clock, cross-species OT, or lineage claim.",
        ha="center", fontsize=8.6, color="#4D5660",
    )
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(out / f"human_mouse_species_integrated_trajectories.{suffix}", dpi=250)
    plt.close(fig)


def render_study_trajectory_map(nodes, arrows, out):
    """Render the requested single shared map with trajectory overlays."""
    out.mkdir(parents=True, exist_ok=True)
    xy = nodes[["UMAP1", "UMAP2"]].to_numpy()
    route_colors = route_palette(nodes.route.astype(str).unique())
    fig, ax = plt.subplots(figsize=(16, 11.5))
    fig.subplots_adjust(left=.04, right=.73, bottom=.07, top=.93)
    for route in sorted(route_colors):
        for species, marker in (("human", "o"), ("mouse", "s")):
            selected = nodes.route.astype(str).eq(route).to_numpy() & nodes.species_key.eq(species).to_numpy()
            if selected.any():
                ax.scatter(*xy[selected].T, s=8, c=[route_colors[route]], marker=marker,
                           alpha=.64, linewidths=0, rasterized=True, zorder=2)
    draw_arrows(ax, arrows, route_colors)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("Human–mouse species-balanced coembedding with study-specific OT trajectories", fontsize=16)
    route_handles = [Line2D([], [], color=route_colors[x], lw=2, label=x) for x in sorted(route_colors)]
    legend1 = ax.legend(handles=route_handles, title="Input study / route", frameon=False,
                        loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=7.4, title_fontsize=8.5)
    ax.add_artist(legend1)
    ax.legend(handles=[Line2D([], [], marker="o", color="#444", linestyle="", label="Human metacell"),
                       Line2D([], [], marker="s", color="#444", linestyle="", label="Mouse metacell")],
              frameon=False, loc="lower left", bbox_to_anchor=(1.01, 0), fontsize=8)
    fig.text(.385, .017,
             "Species-balanced scVI from one-to-one ortholog counts. Arrows: frozen study-internal OT only "
             "(available for 7 of 11 routes); no time link crosses species.",
             ha="center", fontsize=8.2, color="#4D5660")
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(out / f"human_mouse_species_integrated_study_trajectory_map.{suffix}", dpi=250)
    plt.close(fig)


def cross_species_fraction(xy, species, neighbors=30):
    tree = cKDTree(xy)
    idx = tree.query(xy, k=neighbors + 1, workers=-1)[1][:, 1:]
    return float(np.mean(species[idx] != species[:, None]))


def relative_temporal_distance(xy, temporal):
    """Weighted temporal-edge length normalized by each node's 30-NN scale."""
    distances = cKDTree(xy).query(xy, k=31, workers=-1)[0][:, 1:]
    scale = np.maximum(np.median(distances, axis=1), 1e-8)
    upper = sparse.triu(temporal.maximum(temporal.T), k=1).tocoo()
    relative = np.linalg.norm(xy[upper.row] - xy[upper.col], axis=1) / np.sqrt(scale[upper.row] * scale[upper.col])
    return float(np.average(relative, weights=upper.data))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    output = args.output
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if not torch.cuda.is_available():
        raise RuntimeError("GPU is required for the species-balanced scVI display model")
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:32])

    project = args.config.resolve().parent
    registered = read_orthologs(project / cfg["orthologs"])
    orthologs = restrict_to_input_orthologs(registered, cfg["human"]["metacell_gex"], cfg["mouse"]["metacell_gex_root"])
    human_obs, human_counts = build_human_matrix(cfg["human"]["metacell_gex"], orthologs)
    mouse_obs, mouse_counts = build_mouse_matrix(
        cfg["mouse"]["metacell_registry"], cfg["mouse"]["metacell_gex_root"], orthologs,
        project / cfg["mouse_state_map"],
    )
    nodes = pd.concat([human_obs, mouse_obs], ignore_index=True, sort=False)
    if len(nodes) != cfg["expected"]["human_metacells"] + cfg["expected"]["mouse_metacells"] or not nodes.metacell_id.is_unique:
        raise ValueError("Frozen metacell registry does not close")
    nodes["semantic_group"] = [semantic_group(a, b) for a, b in zip(nodes.species_key, nodes.coarse_state)]
    nodes["study_group"] = nodes.species_key.astype(str) + "::" + nodes.route.astype(str)
    metacell_index = pd.Index(nodes.metacell_id.astype(str))
    raw_counts = sparse.vstack([human_counts, mouse_counts], format="csr")
    latent = fit_species_balanced_scvi(raw_counts, nodes, orthologs, cfg["scvi"], cfg["seed"], output)
    orthologs.to_csv(output / "one_to_one_orthologs_used.tsv.gz", sep="\t", index=False)

    local, radius = local_kernel(latent, nodes.study_group.to_numpy(), cfg["graph"]["local_neighbors"])
    cross, mnn_checks = mnn_kernel(
        latent, nodes.study_group.to_numpy(), nodes.species_key.to_numpy(), radius,
        cfg["graph"]["cross_neighbors"], cfg["graph"]["cross_radius_gate"], cfg["seed"],
    )
    human_codes = map_human_cells(cfg["human"]["time_h5ad"], cfg["human"]["cell_to_metacell"], metacell_index)
    mouse_codes = map_mouse_cells(cfg["mouse"]["cell_coordinates"], cfg["mouse"]["cell_to_metacell_root"], metacell_index)
    temporal = normalize_rows(
        aggregate_cell_kernel(cfg["human"]["time_kernel"], human_codes, len(nodes)) +
        aggregate_cell_kernel(cfg["mouse"]["time_kernel"], mouse_codes, len(nodes))
    ).astype(np.float32)
    for name, matrix in (("local_expression_kernel.npz", local), ("cross_study_expression_mnn_kernel.npz", cross),
                         ("within_species_time_kernel.npz", temporal)):
        sparse.save_npz(output / name, matrix)

    xy, graph = fit_layout(latent, local, temporal, cross, cfg["graph"]["time_weights"], cfg["graph"], cfg["seed"])
    control_xy, control_graph = fit_layout(
        latent, local, temporal, cross, cfg["graph"]["control_weights"], cfg["graph"], cfg["seed"],
    )
    nodes[["UMAP1", "UMAP2"]] = xy
    nodes[["UMAP1_control", "UMAP2_control"]] = control_xy
    sparse.save_npz(output / "time_constrained_graph.npz", graph)
    sparse.save_npz(output / "no_time_control_graph.npz", control_graph)

    human_forward = aggregate_cell_kernel(cfg["human"]["forward_time_kernel"], human_codes, len(nodes))
    mouse_forward = mouse_forward_kernel(project / "sources/local/mouse-metacell-ot-top3-edges.tsv.gz", metacell_index, len(nodes))
    arrows = pd.concat([
        arrows_from_forward_kernel(human_forward, nodes, "human"),
        arrows_from_forward_kernel(mouse_forward, nodes, "mouse"),
    ], ignore_index=True)
    arrows["route"] = arrows.source_group.str.split(" | ", regex=False).str[0]
    if not arrows.route.isin(nodes.route.astype(str)).all():
        raise ValueError("A displayed OT arrow cannot be assigned to a source study")
    nodes.to_csv(output / "common_metacell_nodes.tsv.gz", sep="\t", index=False)
    arrows.to_csv(output / "within_species_study_trajectory_arrows.tsv", sep="\t", index=False)
    mnn_checks.to_csv(output / "cross_study_mnn_checks.tsv", sep="\t", index=False)
    render(nodes, arrows, cross, output / "figures")
    render_study_trajectory_map(nodes, arrows, output / "figures")

    species = nodes.species_key.to_numpy()
    coo = cross.tocoo()
    cross_mask = (species[coo.row] != species[coo.col]) & (coo.row < coo.col)
    broad_match = nodes.semantic_group.to_numpy()[coo.row[cross_mask]] == nodes.semantic_group.to_numpy()[coo.col[cross_mask]]
    time_relative = relative_temporal_distance(xy, temporal)
    control_relative = relative_temporal_distance(control_xy, temporal)
    result = {
        "status": "DISPLAY_CANDIDATE_COMPLETE_VALIDATION_PENDING",
        "analysis_unit": "frozen within-library metacell",
        "human_metacells": int((species == "human").sum()),
        "mouse_metacells": int((species == "mouse").sum()),
        "one_to_one_orthologs": int(len(orthologs)),
        "scvi": {**cfg["scvi"], "labels_time_apa_used": False},
        "time_weights": cfg["graph"]["time_weights"],
        "control_weights": cfg["graph"]["control_weights"],
        "time_graph_components": int(connected_components(graph, directed=False)[0]),
        "time_graph_undirected_edges": int(graph.nnz // 2),
        "cross_species_mnn_undirected_edges": int(cross_mask.sum()),
        "cross_species_30nn_fraction": cross_species_fraction(xy, species),
        "control_cross_species_30nn_fraction": cross_species_fraction(control_xy, species),
        "time_edge_relative_distance": time_relative,
        "control_time_edge_relative_distance": control_relative,
        "time_edge_relative_distance_change_vs_control_pct": float(100 * (time_relative / control_relative - 1)),
        "cross_species_mnn_broad_semantic_match_fraction": float(broad_match.mean()) if len(broad_match) else None,
        "study_specific_trajectory_arrows": int(len(arrows)),
        "cross_species_ot_or_common_clock_created": False,
        "cross_species_lineage_or_identity_validated": False,
    }
    (output / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
