#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import SEACells


def mode_and_purity(values: pd.Series) -> tuple[str, float]:
    counts = values.fillna("missing").astype(str).value_counts()
    return str(counts.index[0]), float(counts.iloc[0] / counts.sum())


def summarize_assignments(cells: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    latent_columns = [column for column in cells if column.startswith("scvi_latent_")]
    for metacell_id, group in cells.groupby("metacell_id", sort=True, observed=True):
        local_state, local_state_purity = mode_and_purity(group.local_state_id)
        phase, phase_purity = mode_and_purity(group.trajectory_phase)
        fate, fate_purity = mode_and_purity(group.fate_axis)
        pseudotime = pd.to_numeric(group.route_pseudotime, errors="coerce").dropna()
        row: dict[str, object] = {
            "metacell_id": metacell_id,
            "n_cells": len(group),
            "local_state_id": local_state,
            "local_state_purity": local_state_purity,
            "trajectory_phase": phase,
            "trajectory_phase_purity": phase_purity,
            "fate_axis": fate,
            "fate_axis_purity": fate_purity,
            "route_pseudotime_q25": pseudotime.quantile(0.25) if len(pseudotime) else np.nan,
            "route_pseudotime_median": pseudotime.median() if len(pseudotime) else np.nan,
            "route_pseudotime_q75": pseudotime.quantile(0.75) if len(pseudotime) else np.nan,
            "route_pseudotime_iqr": (
                pseudotime.quantile(0.75) - pseudotime.quantile(0.25)
                if len(pseudotime)
                else np.nan
            ),
            "umap_1": group.umap_1.mean(),
            "umap_2": group.umap_2.mean(),
        }
        row.update({column: group[column].mean() for column in latent_columns})
        rows.append(row)
    return pd.DataFrame(rows).sort_values("metacell_id")


def aggregate_counts(adata: ad.AnnData, cells: pd.DataFrame, metadata: pd.DataFrame) -> ad.AnnData:
    labels = pd.Categorical(cells.metacell_id, categories=metadata.metacell_id)
    grouping = sparse.csr_matrix(
        (np.ones(adata.n_obs), (labels.codes, np.arange(adata.n_obs))),
        shape=(len(labels.categories), adata.n_obs),
    )
    counts = grouping @ adata.layers["counts"]
    result = ad.AnnData(X=counts.tocsr(), obs=metadata.set_index("metacell_id"), var=adata.var.copy())
    result.layers["counts"] = result.X.copy()
    return result


def plot_audit(cells: pd.DataFrame, metadata: pd.DataFrame, output: Path) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    states = sorted(cells.local_state_id.astype(str).unique())
    state_codes = pd.Categorical(cells.local_state_id.astype(str), categories=states).codes
    axes[0, 0].scatter(cells.umap_1, cells.umap_2, c=state_codes, s=2, cmap="tab20", alpha=0.45)
    axes[0, 0].set_title("Cells: frozen local state")

    axes[0, 1].scatter(cells.umap_1, cells.umap_2, c="#d9d9d9", s=1, alpha=0.25)
    points = axes[0, 1].scatter(
        metadata.umap_1,
        metadata.umap_2,
        c=metadata.route_pseudotime_median,
        s=20 + 0.4 * metadata.n_cells,
        cmap="viridis",
        edgecolor="black",
        linewidth=0.25,
    )
    figure.colorbar(points, ax=axes[0, 1], label="Median route pseudotime")
    axes[0, 1].set_title("Metacell centroids")

    points = axes[0, 2].scatter(
        cells.umap_1,
        cells.umap_2,
        c=cells.route_pseudotime,
        s=2,
        cmap="viridis",
        alpha=0.55,
    )
    figure.colorbar(points, ax=axes[0, 2], label="Route pseudotime")
    axes[0, 2].set_title("Cells: time projection")

    axes[1, 0].hist(metadata.n_cells, bins=20, color="#4477aa", edgecolor="white")
    axes[1, 0].set_xlabel("Cells per metacell")
    axes[1, 0].set_ylabel("Metacells")
    axes[1, 0].set_title("Aggregation size")

    axes[1, 1].hist(metadata.local_state_purity, bins=np.linspace(0, 1, 21), alpha=0.7, label="local state")
    axes[1, 1].hist(metadata.fate_axis_purity, bins=np.linspace(0, 1, 21), alpha=0.7, label="fate")
    axes[1, 1].set_xlabel("Dominant-label fraction")
    axes[1, 1].set_title("Biological label purity")
    axes[1, 1].legend(frameon=False)

    axes[1, 2].scatter(
        metadata.route_pseudotime_median,
        metadata.route_pseudotime_iqr,
        c=metadata.local_state_purity,
        cmap="magma",
        s=25,
        edgecolor="black",
        linewidth=0.2,
    )
    axes[1, 2].set_xlabel("Median route pseudotime")
    axes[1, 2].set_ylabel("Within-metacell pseudotime IQR")
    axes[1, 2].set_title("Temporal compactness")

    for axis in axes.flat[:3]:
        axis.set_xticks([])
        axis.set_yticks([])
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--latent", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--library-key", required=True)
    parser.add_argument("--target-cells", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    if adata.obs.library_key.astype(str).nunique() != 1 or str(adata.obs.library_key.iloc[0]) != args.library_key:
        raise ValueError("input AnnData does not contain exactly the requested library")
    if "counts" not in adata.layers:
        raise ValueError("raw GEX counts layer is required")
    if adata.obs_names.duplicated().any():
        raise ValueError("duplicate cell IDs in GEX input")

    annotations = pd.read_csv(args.annotation, sep="\t")
    annotations = annotations[annotations.cell_id.isin(adata.obs_names)].copy()
    annotation_fields = [
        "cell_id",
        "route_pseudotime",
        "local_state_id",
        "trajectory_phase",
        "fate_axis",
        "biology_label",
        "confidence",
    ]
    annotations = annotations[annotation_fields].drop_duplicates("cell_id")
    if len(annotations) != adata.n_obs:
        raise ValueError(f"annotation join is incomplete: {len(annotations)} of {adata.n_obs} cells")
    latent = pd.read_csv(args.latent, sep="\t")
    latent_columns = [column for column in latent if column.startswith("scvi_latent_")]
    latent = latent[["cell_id", *latent_columns]].drop_duplicates("cell_id")
    if len(latent) != adata.n_obs or set(latent.cell_id) != set(adata.obs_names):
        raise ValueError("scVI latent join does not exactly match the GEX cells")

    work = adata.copy()
    work.X = work.layers["counts"].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    n_components = min(30, work.n_obs - 1, work.n_vars - 1)
    sc.pp.pca(work, n_comps=n_components, random_state=args.seed)
    sc.pp.neighbors(work, n_neighbors=30, use_rep="X_pca", random_state=args.seed)
    sc.tl.umap(work, random_state=args.seed)

    n_metacells = max(2, math.ceil(work.n_obs / args.target_cells))
    model = SEACells.core.SEACells(
        work,
        build_kernel_on="X_pca",
        n_SEACells=n_metacells,
        use_gpu=False,
        verbose=True,
        n_waypoint_eigs=min(10, n_components - 1),
        n_neighbors=15,
        convergence_epsilon=1e-3,
        use_sparse=True,
    )
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=100)
    assignments = model.get_hard_assignments().rename(columns={"SEACell": "seacell_label"})
    if not assignments.index.equals(work.obs_names):
        assignments = assignments.reindex(work.obs_names)
    if assignments.seacell_label.isna().any():
        raise ValueError("SEACells left cells unassigned")

    label_order = sorted(assignments.seacell_label.astype(str).unique())
    label_map = {
        label: f"{args.library_key}.t{args.target_cells:03d}.MC{index:04d}"
        for index, label in enumerate(label_order, start=1)
    }
    cells = work.obs.copy()
    cells["metacell_id"] = assignments.seacell_label.astype(str).map(label_map).to_numpy()
    cells["umap_1"] = work.obsm["X_umap"][:, 0]
    cells["umap_2"] = work.obsm["X_umap"][:, 1]
    cells = cells.reset_index(drop=True).merge(annotations, on="cell_id", how="left", validate="one_to_one")
    cells = cells.merge(latent, on="cell_id", how="left", validate="one_to_one")
    if not cells.cell_id.astype(str).equals(pd.Series(adata.obs_names.astype(str))):
        raise ValueError("cell order changed while attaching annotation and latent metadata")

    metadata = summarize_assignments(cells)
    metacell_counts = aggregate_counts(adata, cells, metadata)
    metacell_counts.write_h5ad(args.output_dir / "metacell_gex_counts.h5ad", compression="lzf")
    cells.to_csv(args.output_dir / "cell_to_metacell.tsv.gz", sep="\t", index=False)
    metadata.to_csv(args.output_dir / "metacell_metadata.tsv", sep="\t", index=False)
    plot_audit(cells, metadata, args.output_dir / "gex_visual_audit.png")
    model.plot_convergence(save_as=args.output_dir / "seacells_convergence.png", show=False)

    summary = {
        "status": "GEX_METACELL_TECHNICAL_PASS_VISUAL_REVIEW_REQUIRED",
        "library_key": args.library_key,
        "target_cells_per_metacell": args.target_cells,
        "n_cells": adata.n_obs,
        "requested_metacells": n_metacells,
        "observed_metacells": len(metadata),
        "assigned_cells": len(cells),
        "assignment_fraction": float(cells.metacell_id.notna().mean()),
        "median_cells_per_metacell": float(metadata.n_cells.median()),
        "minimum_cells_per_metacell": int(metadata.n_cells.min()),
        "maximum_cells_per_metacell": int(metadata.n_cells.max()),
        "median_local_state_purity": float(metadata.local_state_purity.median()),
        "median_phase_purity": float(metadata.trajectory_phase_purity.median()),
        "median_fate_purity": float(metadata.fate_axis_purity.median()),
        "median_pseudotime_iqr": float(metadata.route_pseudotime_iqr.median()),
        "clustering_changed": False,
        "time_used_to_construct_metacells": False,
        "pas_used_to_construct_metacells": False,
    }
    (args.output_dir / "gex_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.output_dir / "gex.ok").write_text(summary["status"] + "\n")


if __name__ == "__main__":
    main()
