#!/usr/bin/env python3
"""Compute post-clustering cell-cycle and stress overlays without changing the graph."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import normalized_mutual_info_score

from cluster_utils import write_json, write_marker


SIGNATURES = {
    "s_phase": [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
        "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
        "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
        "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
        "CHAF1B", "BRIP1", "E2F8",
    ],
    "g2m_phase": [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
        "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
        "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1",
        "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2",
        "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3",
        "GAS2L3", "CBX5", "CENPA",
    ],
    "immediate_early_stress": [
        "FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "DUSP1", "DUSP2", "EGR1", "EGR2", "IER2",
        "IER3", "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "PPP1R15A", "GADD45A", "DDIT3", "NR4A1",
    ],
    "interferon_response": [
        "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI27", "IFI44", "IFI44L", "MX1", "MX2", "OAS1",
        "OAS2", "OAS3", "OASL", "RSAD2", "BST2", "IRF7", "STAT1", "DDX58", "IFIH1",
    ],
}


def signature_score(matrix: sparse.csr_matrix, indices: np.ndarray) -> np.ndarray:
    subset = matrix[:, indices].astype(np.float64)
    mean = np.asarray(subset.mean(axis=0)).ravel()
    mean_sq = np.asarray(subset.power(2).mean(axis=0)).ravel()
    std = np.sqrt(np.maximum(mean_sq - mean * mean, 0.0))
    keep = std > 1e-8
    if not np.any(keep):
        return np.zeros(matrix.shape[0], dtype=np.float32)
    weights = 1.0 / std[keep]
    centered = np.asarray(subset[:, keep].dot(weights)).ravel() - float(np.sum(mean[keep] * weights))
    return (centered / int(keep.sum())).astype(np.float32)


def add_program_review_flags(
    cluster_scores: pd.DataFrame,
    maximum_cell_cycle_polarization: float,
    maximum_positive_program_mean: float,
) -> pd.DataFrame:
    """Flag program-polarized clusters for review without changing their labels."""
    result = cluster_scores.copy()
    for column in [
        "phase_fraction_G1_like", "phase_fraction_S_like", "phase_fraction_G2M_like"
    ]:
        if column not in result:
            result[column] = 0.0
        result[column] = result[column].fillna(0.0).astype(float)
    result["phase_fraction_cycling_like"] = (
        result["phase_fraction_S_like"] + result["phase_fraction_G2M_like"]
    )
    result["cell_cycle_polarization"] = result[
        ["phase_fraction_G1_like", "phase_fraction_cycling_like"]
    ].max(axis=1)
    result["cell_cycle_polarized_candidate"] = (
        result["cell_cycle_polarization"] >= maximum_cell_cycle_polarization
    )
    result["positive_stress_or_interferon_mean"] = result[
        ["immediate_early_stress_mean", "interferon_response_mean"]
    ].max(axis=1)
    result["stress_or_interferon_high_candidate"] = (
        result["positive_stress_or_interferon_mean"] >= maximum_positive_program_mean
    )
    reasons = []
    for row in result.itertuples(index=False):
        row_reasons = []
        if row.cell_cycle_polarized_candidate:
            row_reasons.append("cell_cycle_polarized")
        if row.stress_or_interferon_high_candidate:
            row_reasons.append("stress_or_interferon_high")
        reasons.append(";".join(row_reasons) if row_reasons else "NONE")
    result["program_review_reasons"] = reasons
    result["requires_program_review"] = result["program_review_reasons"].ne("NONE")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--umap", required=True)
    parser.add_argument("--maximum-cell-cycle-polarization", type=float, required=True)
    parser.add_argument("--maximum-positive-program-mean", type=float, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    matrix = sparse.load_npz(args.matrix).tocsr()
    features = pd.read_csv(args.features, sep="\t", dtype=str)
    clusters = pd.read_csv(args.clusters, sep="\t", dtype=str)
    umap = pd.read_csv(args.umap, sep="\t", dtype={"cell_id": str})
    if matrix.shape[0] != len(clusters) or not clusters["cell_id"].equals(umap["cell_id"]):
        raise ValueError("program audit inputs have inconsistent cell order")
    symbols = features["gene_symbol"].astype(str).str.upper()
    symbol_to_index = pd.Series(np.arange(len(symbols)), index=symbols).groupby(level=0).first()
    scores = clusters[["cell_id", "cluster_id"]].copy()
    coverage_rows = []
    for name, genes in SIGNATURES.items():
        found = [gene for gene in genes if gene in symbol_to_index.index]
        if len(found) < 3:
            raise ValueError(f"signature {name} has insufficient genes in the expression feature universe")
        indices = symbol_to_index.loc[found].to_numpy(dtype=int)
        scores[name] = signature_score(matrix, indices)
        coverage_rows.append({
            "program": name,
            "requested_genes": len(genes),
            "found_genes": len(found),
            "coverage_fraction": len(found) / len(genes),
            "found_gene_symbols": ";".join(found),
        })
    scores["cell_cycle_phase"] = np.where(
        (scores["s_phase"] <= 0) & (scores["g2m_phase"] <= 0),
        "G1_like",
        np.where(scores["s_phase"] >= scores["g2m_phase"], "S_like", "G2M_like"),
    )
    scores.to_csv(
        target / "cell_program_scores.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(coverage_rows).to_csv(target / "signature_coverage.tsv", sep="\t", index=False)
    program_columns = list(SIGNATURES)
    cluster_scores = scores.groupby("cluster_id", sort=True)[program_columns].agg(["mean", "median"])
    cluster_scores.columns = [f"{program}_{stat}" for program, stat in cluster_scores.columns]
    phase = pd.crosstab(scores["cluster_id"], scores["cell_cycle_phase"], normalize="index")
    cluster_scores = cluster_scores.join(phase.add_prefix("phase_fraction_")).reset_index()
    cluster_scores = add_program_review_flags(
        cluster_scores,
        maximum_cell_cycle_polarization=args.maximum_cell_cycle_polarization,
        maximum_positive_program_mean=args.maximum_positive_program_mean,
    )
    cluster_scores.to_csv(target / "cluster_program_scores.tsv", sep="\t", index=False)
    nmi = float(normalized_mutual_info_score(scores["cluster_id"], scores["cell_cycle_phase"]))
    review_clusters = cluster_scores.loc[
        cluster_scores["requires_program_review"], "cluster_id"
    ].astype(str).tolist()
    status = "REVIEW" if review_clusters else "PASS"
    write_json(target / "summary.json", {
        "n_cells": len(scores),
        "n_clusters": int(scores["cluster_id"].nunique()),
        "cluster_cell_cycle_phase_normalized_mutual_information": nmi,
        "used_for_pca_neighbors_or_clustering": False,
        "score_definition": "mean per-gene z score on normalized log1p expression",
        "interpretation": "orthogonal overlay; not a reprogramming-state label",
        "maximum_cell_cycle_polarization_without_review": args.maximum_cell_cycle_polarization,
        "maximum_positive_program_mean_without_review": args.maximum_positive_program_mean,
        "n_program_review_clusters": len(review_clusters),
        "program_review_clusters": review_clusters,
        "review_semantics": "manual interpretation required; clusters and cells are retained unchanged",
        "status": status,
    })

    plot = umap[["umap_1", "umap_2"]].to_numpy(float)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    for ax, program in zip(axes.ravel(), program_columns):
        values = scores[program].to_numpy(float)
        low, high = np.quantile(values, [0.01, 0.99])
        image = ax.scatter(plot[:, 0], plot[:, 1], c=np.clip(values, low, high), s=0.8,
                           linewidths=0, cmap="coolwarm", rasterized=True)
        ax.set(title=program.replace("_", " "), xticks=[], yticks=[])
        fig.colorbar(image, ax=ax, shrink=0.75)
    fig.savefig(target / "program_overlay.png", dpi=220)
    plt.close(fig)
    write_marker(target / ".ok", status=status, n_cells=len(scores), cell_cycle_nmi=nmi,
                 n_program_review_clusters=len(review_clusters))


if __name__ == "__main__":
    main()
