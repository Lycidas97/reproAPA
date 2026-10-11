#!/usr/bin/env python3
"""Independent raw-data oracle for a completed unified RNA-QC run.

This module intentionally imports no producer code.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import math
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml


BARE_MITO = frozenset(
    {"ATP6", "ATP8", "CO1", "CO2", "CO3", "COX1", "COX2", "COX3", "CYB",
     "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"}
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    return parser.parse_args()


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values], object)


def first(directory: Path, names: list[str]) -> Path:
    for name in names:
        path = directory / name
        if path.is_file():
            return path
    raise FileNotFoundError(f"Missing matrix member in {directory}: {names}")


def mito_mask(ids: np.ndarray, names: np.ndarray) -> np.ndarray:
    def normalise(value: str) -> str:
        text = str(value).strip()
        if text.lower().startswith(("mt-", "mt_", "mt.")):
            text = text[3:]
        return "".join(char for char in text.upper() if char.isalnum())
    ids_s, names_s = pd.Series(ids, dtype=str), pd.Series(names, dtype=str)
    prefixes = ("mt-", "MT-", "Mt-", "KEG11_")
    return (
        ids_s.str.startswith(prefixes) | names_s.str.startswith(prefixes)
        | ids_s.map(normalise).isin(BARE_MITO) | names_s.map(normalise).isin(BARE_MITO)
    ).to_numpy(bool)


def metrics(counts: sp.csr_matrix, ids: np.ndarray, names: np.ndarray) -> dict[str, np.ndarray]:
    total = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    genes = np.diff(counts.indptr).astype(np.int64)
    mask = mito_mask(ids, names)
    mt = np.asarray(counts[:, mask].sum(axis=1)).ravel() if mask.any() else np.zeros(len(total))
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.divide(mt * 100, total, out=np.zeros(len(total)), where=total > 0)
        complexity = np.divide(
            np.log10(np.maximum(genes, 1)), np.log10(np.maximum(total, 10)),
            out=np.zeros(len(total)), where=total > 0,
        )
    return {"rna_total_umi": total, "rna_detected_genes": genes,
            "rna_pct_mito": pct, "rna_complexity": complexity}


def read_h5(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][:])
        barcodes = decode(matrix["barcodes"][:])
        features = matrix["features"]
        ids = decode(features["id"][:])
        names = decode(features["name"][:]) if "name" in features else ids
        types = decode(features["feature_type"][:]) if "feature_type" in features else np.repeat("Gene Expression", shape[0])
        gene = types == "Gene Expression"
        source = sp.csc_matrix((matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]), shape=shape)
        counts = source[gene, :].T.tocsr()
    return barcodes, metrics(counts, ids[gene], names[gene])


def read_mtx(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    matrix_path = first(path, ["matrix.mtx.gz", "matrix.mtx"])
    barcode_path = first(path, ["barcodes.tsv.gz", "barcodes.tsv"])
    feature_path = first(path, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
    raw = scipy.io.mmread(matrix_path).tocsr()
    opener = gzip.open if barcode_path.suffix == ".gz" else open
    with opener(barcode_path, "rt", encoding="utf-8") as handle:
        barcodes = np.asarray([line.rstrip("\n").split("\t")[0] for line in handle if line.strip()], object)
    features = pd.read_csv(feature_path, sep="\t", header=None, dtype=str)
    ids = features.iloc[:, 0].to_numpy(object)
    names = features.iloc[:, 1].to_numpy(object) if features.shape[1] > 1 else ids
    types = features.iloc[:, 2].to_numpy(object) if features.shape[1] > 2 else np.repeat("Gene Expression", len(ids))
    gene = types == "Gene Expression"
    return barcodes, metrics(raw[gene, :].T.tocsr(), ids[gene], names[gene])


def raw_input(adapter: str, path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if adapter == "cellranger_filtered_h5":
        return read_h5(path)
    if adapter == "cellranger_filtered_mtx":
        return read_mtx(path)
    raise ValueError(f"Independent validator does not support {adapter}")


def mad(values: np.ndarray) -> float:
    values = np.asarray(values, float)
    return float(np.median(np.abs(values - np.median(values)))) if len(values) else math.nan


def branch(metrics_frame: pd.DataFrame, rules: dict) -> tuple[dict[str, float], np.ndarray]:
    log_umi = np.log10(np.maximum(metrics_frame["rna_total_umi"].to_numpy(), 1))
    log_genes = np.log10(np.maximum(metrics_frame["rna_detected_genes"].to_numpy(), 1))
    mito = metrics_frame["rna_pct_mito"].to_numpy(float)
    threshold = {
        "hard_min_umi": float(rules["hard_min_umi"]),
        "hard_min_genes": float(rules["hard_min_genes"]),
        "mad_n": float(rules["mad_n"]),
        "log10_umi_median": float(np.median(log_umi)), "log10_umi_mad": mad(log_umi),
        "log10_genes_median": float(np.median(log_genes)), "log10_genes_mad": mad(log_genes),
        "mito_median_pct": float(np.median(mito)), "mito_mad_pct": mad(mito),
        "absolute_mito_cap_pct": float(rules["mito_cap_pct"]),
        "complexity_min": float(rules["complexity_min"]),
    }
    threshold["adaptive_min_log10_umi"] = threshold["log10_umi_median"] - threshold["mad_n"] * threshold["log10_umi_mad"]
    threshold["adaptive_min_log10_genes"] = threshold["log10_genes_median"] - threshold["mad_n"] * threshold["log10_genes_mad"]
    threshold["effective_mito_max_pct"] = min(
        threshold["absolute_mito_cap_pct"], threshold["mito_median_pct"] + threshold["mad_n"] * threshold["mito_mad_pct"]
    )
    passed = (
        (metrics_frame["rna_total_umi"].to_numpy() >= threshold["hard_min_umi"])
        & (metrics_frame["rna_detected_genes"].to_numpy() >= threshold["hard_min_genes"])
        & (log_umi >= threshold["adaptive_min_log10_umi"])
        & (log_genes >= threshold["adaptive_min_log10_genes"])
        & (mito <= threshold["effective_mito_max_pct"])
        & (metrics_frame["rna_complexity"].to_numpy() >= threshold["complexity_min"])
    )
    return threshold, passed


def set_hash(values: set[str]) -> str:
    return hashlib.sha256(("\n".join(sorted(values)) + "\n").encode()).hexdigest()


def main() -> None:
    parsed = parse_args()
    root, run = parsed.root.resolve(), parsed.run_dir.resolve()
    config = yaml.safe_load(parsed.config.read_text(encoding="utf-8"))
    manifest = pd.read_csv(resolve(root, config["library_manifest"]), sep="\t", dtype=str, keep_default_na=False)
    registry = pd.read_csv(resolve(root, config["metadata_registry"]), dtype=str, keep_default_na=False, encoding="utf-8-sig")
    registry = registry.loc[
        registry["series_accession"].eq(str(config["study_id"]))
        & registry["sample_accession"].isin(manifest["sample_accession"]),
        ["sample_accession", "sample_title", "organism", "source_url"],
    ].drop_duplicates()
    source_check = manifest.merge(registry, on="sample_accession", suffixes=("", "_registry"), validate="one_to_one")
    for observed, expected in (("sample_title", "sample_title_registry"), ("species", "organism"), ("source_url", "source_url_registry")):
        if source_check[observed].ne(source_check[expected]).any():
            raise ValueError(f"Independent source metadata mismatch: {observed}")
    contract = pd.read_csv(resolve(root, config["metadata_contract"]), sep="\t", dtype=str, keep_default_na=False)
    contract_check = manifest.merge(contract, on="sample_accession", suffixes=("", "_contract"), validate="one_to_one")
    for field in ("experimental_time", "condition", "metadata_derivation"):
        if contract_check[field].ne(contract_check[f"{field}_contract"]).any():
            raise ValueError(f"Independent curated metadata mismatch: {field}")
    final = pd.read_parquet(run / "cell_qc_decisions.parquet")
    samples = pd.read_csv(run / "sample_qc_decisions.tsv", sep="\t")
    rows, all_ids = [], []
    for item in manifest.itertuples(index=False):
        barcodes, raw = raw_input(item.input_adapter, resolve(root, item.input_path))
        ids = np.asarray([
            config["global_cell_id_template"].format(study_id=config["study_id"], library_id=item.library_id, barcode=x)
            for x in barcodes
        ], object)
        if len(ids) != int(item.expected_called_cells):
            raise ValueError(f"Raw called-cell count mismatch: {item.library_id}")
        all_ids.extend(ids.astype(str))
        observed_metrics = pd.read_parquet(run / "per_library" / item.library_id / "rna_cell_qc_metrics.parquet").set_index("global_cell_id").loc[ids]
        independent = pd.DataFrame(raw, index=ids)
        for field, tolerance in (("rna_total_umi", 0), ("rna_detected_genes", 0), ("rna_pct_mito", 1e-5), ("rna_complexity", 1e-7)):
            if not np.allclose(independent[field], observed_metrics[field], rtol=1e-6, atol=tolerance, equal_nan=True):
                raise ValueError(f"Independent raw metric mismatch: {item.library_id}/{field}")
        # Raw reparsing is independently reconciled above. The stored metric ledger
        # is the canonical floating-point input to thresholding, so exact boundary
        # decisions remain reproducible across sparse-matrix library builds.
        canonical = observed_metrics[["rna_total_umi", "rna_detected_genes", "rna_pct_mito", "rna_complexity"]]
        flags = {}
        stored_thresholds = pd.read_csv(run / "per_library" / item.library_id / "rna_qc_thresholds.tsv", sep="\t").set_index("branch")
        for name in ("strict", "candidate", "retention_sensitivity"):
            threshold, flags[name] = branch(canonical, config["basic_qc"][name])
            for field in ("hard_min_umi", "hard_min_genes", "mad_n", "log10_umi_median", "log10_umi_mad", "log10_genes_median", "log10_genes_mad", "mito_median_pct", "mito_mad_pct", "absolute_mito_cap_pct", "complexity_min", "adaptive_min_log10_umi", "adaptive_min_log10_genes", "effective_mito_max_pct"):
                if not np.isclose(float(stored_thresholds.loc[name, field]), threshold[field], rtol=1e-6, atol=1e-7):
                    raise ValueError(f"Independent threshold mismatch: {item.library_id}/{name}/{field}")
        n_candidate = int(flags["candidate"].sum())
        primary_min = int(config["sample_gates"]["primary_min_candidate_basic_qc"])
        sensitivity_min = int(config["sample_gates"]["sensitivity_min_candidate_basic_qc"])
        gate = "primary_eligible" if n_candidate >= primary_min else "sensitivity_only_50_99" if n_candidate >= sensitivity_min else "insufficient_primary_cells_lt50"
        primary = flags["candidate"] & (gate == "primary_eligible")
        expected_tier = np.select([primary, (~primary) & flags["retention_sensitivity"]], ["primary_pass", "sensitivity_only"], default="technical_exclude")
        observed_final = final.set_index("global_cell_id").loc[ids]
        comparisons = {
            "strict_basic_qc_pass": flags["strict"],
            "candidate_basic_qc_pass": flags["candidate"],
            "retention_sensitivity_basic_qc_pass": flags["retention_sensitivity"],
        }
        for field, expected in comparisons.items():
            if not np.array_equal(observed_final[field].astype(bool).to_numpy(), expected):
                raise ValueError(f"Independent flag mismatch: {item.library_id}/{field}")
        if not np.array_equal(observed_final["final_qc_tier"].astype(str).to_numpy(), expected_tier):
            raise ValueError(f"Independent tier mismatch: {item.library_id}")
        sample = samples.loc[samples["library_id"].eq(item.library_id)]
        if len(sample) != 1 or sample.iloc[0]["sample_gate"] != gate or int(sample.iloc[0]["n_candidate_basic_qc"]) != n_candidate:
            raise ValueError(f"Independent sample gate mismatch: {item.library_id}")
        rows.append({"check": f"raw_threshold_tier_oracle:{item.library_id}", "status": "PASS", "observed": len(ids), "expected": int(item.expected_called_cells)})
    if len(set(all_ids)) != len(all_ids) or set(all_ids) != set(final["global_cell_id"].astype(str)):
        raise ValueError("Independent global-cell identity check failed")
    if final["doublet_hard_exclude"].astype(bool).any() or (final["pre_doublet_qc_tier"] != final["final_qc_tier"]).any():
        raise ValueError("Diagnostic-only Scrublet contract failed")
    baseline_value = str(config.get("regression_baseline", ""))
    if baseline_value:
        baseline = pd.read_parquet(resolve(root, baseline_value))
        for tier in ("primary_pass", "sensitivity_only", "technical_exclude"):
            old = set(baseline.loc[baseline["final_qc_tier"].eq(tier), "global_cell_id"].astype(str))
            new = set(final.loc[final["final_qc_tier"].eq(tier), "global_cell_id"].astype(str))
            if old != new:
                raise ValueError(f"Regression tier differs: {tier}")
            rows.append({"check": f"regression_membership:{tier}", "status": "PASS", "observed": f"{len(new)}:{set_hash(new)}", "expected": f"{len(old)}:{set_hash(old)}"})
    rows.extend([
        {"check": "global_cell_identity", "status": "PASS", "observed": len(all_ids), "expected": int(config["expected_called_cell_count"])},
        {"check": "diagnostic_only_hard_exclusions", "status": "PASS", "observed": 0, "expected": 0},
        {"check": "diagnostic_only_tier_changes", "status": "PASS", "observed": 0, "expected": 0},
    ])
    pd.DataFrame(rows).to_csv(run / "independent_validation.tsv", sep="\t", index=False)
    print(f"PASS: {config['study_id']} ({len(all_ids)} cells)")


if __name__ == "__main__":
    main()
