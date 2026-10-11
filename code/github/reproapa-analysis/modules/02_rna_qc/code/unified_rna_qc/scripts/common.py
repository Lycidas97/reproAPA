#!/usr/bin/env python3
"""Configuration-driven primitives for unified RNA QC v3.2."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scrublet as scr
import yaml


BARE_MITO = frozenset(
    {"ATP6", "ATP8", "CO1", "CO2", "CO3", "COX1", "COX2", "COX3", "CYB",
     "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2"}
)
REQUIRED_MANIFEST_COLUMNS = {
    "study_id", "species", "dataset_id", "library_id", "sample_id",
    "sample_accession", "sample_title", "experimental_time", "condition",
    "author_state", "metadata_derivation", "source_url",
    "input_adapter", "input_path", "expected_called_cells",
    "basic_qc_evidence_adapter", "prior_basic_qc_path",
    "diagnostic_evidence_adapter", "prior_diagnostic_path",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def hash_path(path: Path) -> tuple[int, str]:
    if path.is_file():
        return path.stat().st_size, sha256_file(path)
    if not path.is_dir():
        raise FileNotFoundError(path)
    members = []
    total = 0
    for child in sorted(p for p in path.rglob("*") if p.is_file()):
        total += child.stat().st_size
        members.append(f"{child.relative_to(path)}\t{child.stat().st_size}\t{sha256_file(child)}")
    return total, sha256_text("\n".join(members) + "\n")


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError("Task configuration must be a mapping")
    return value


def resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def read_manifest(root: Path, config: dict[str, Any]) -> pd.DataFrame:
    path = resolve(root, str(config["library_manifest"]))
    if not path.is_file():
        raise FileNotFoundError(f"Missing library manifest: {path}")
    manifest = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = REQUIRED_MANIFEST_COLUMNS - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest columns missing: {sorted(missing)}")
    if len(manifest) != int(config["expected_library_count"]):
        raise ValueError("Configured expected_library_count does not match manifest")
    if set(manifest["study_id"]) != {str(config["study_id"])}:
        raise ValueError("Manifest study_id does not match configuration")
    if set(manifest["species"]) != {str(config["species"])}:
        raise ValueError("Manifest species does not match configuration")
    pattern = re.compile(str(config["library_id_pattern"]))
    bad = [value for value in manifest["library_id"] if pattern.fullmatch(value) is None]
    if bad:
        raise ValueError(f"Library identifiers violate library_id_pattern: {bad[:3]}")
    if manifest["library_id"].duplicated().any():
        raise ValueError("Duplicate library_id in manifest")
    registry_path = resolve(root, str(config["metadata_registry"]))
    registry = pd.read_csv(registry_path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    required_registry = {"series_accession", "sample_accession", "sample_title", "organism", "source_url"}
    if required_registry - set(registry.columns):
        raise ValueError(f"Metadata registry columns missing: {sorted(required_registry - set(registry.columns))}")
    registry = registry.loc[
        registry["series_accession"].eq(str(config["study_id"]))
        & registry["sample_accession"].isin(manifest["sample_accession"]),
        list(required_registry),
    ].drop_duplicates()
    if registry["sample_accession"].duplicated().any() or set(registry["sample_accession"]) != set(manifest["sample_accession"]):
        raise ValueError("Metadata registry does not provide one authoritative row per manifest sample")
    checked = manifest.merge(registry, on="sample_accession", suffixes=("", "_registry"), validate="one_to_one")
    comparisons = {
        "sample_title": "sample_title_registry",
        "species": "organism",
        "source_url": "source_url_registry",
    }
    for observed, expected in comparisons.items():
        bad = checked.loc[checked[observed].ne(checked[expected]), "sample_accession"].tolist()
        if bad:
            raise ValueError(f"Manifest {observed} disagrees with metadata registry: {bad[:3]}")
    if manifest["metadata_derivation"].str.strip().eq("").any():
        raise ValueError("Every sample requires a non-empty metadata_derivation")
    contract = pd.read_csv(resolve(root, str(config["metadata_contract"])), sep="\t", dtype=str, keep_default_na=False)
    contract_columns = {"sample_accession", "experimental_time", "condition", "metadata_derivation"}
    if contract_columns - set(contract.columns):
        raise ValueError(f"Metadata contract columns missing: {sorted(contract_columns - set(contract.columns))}")
    if contract["sample_accession"].duplicated().any() or set(contract["sample_accession"]) != set(manifest["sample_accession"]):
        raise ValueError("Metadata contract does not provide one row per manifest sample")
    curated = manifest.merge(contract[list(contract_columns)], on="sample_accession", suffixes=("", "_contract"), validate="one_to_one")
    for field in ("experimental_time", "condition", "metadata_derivation"):
        bad = curated.loc[curated[field].ne(curated[f"{field}_contract"]), "sample_accession"].tolist()
        if bad:
            raise ValueError(f"Manifest {field} disagrees with reviewed metadata contract: {bad[:3]}")
    manifest["expected_called_cells"] = manifest["expected_called_cells"].astype(int)
    if int(manifest["expected_called_cells"].sum()) != int(config["expected_called_cell_count"]):
        raise ValueError("Configured expected_called_cell_count does not match manifest")
    allowed_raw = set(config["adapters"]["raw_input"])
    allowed_basic = set(config["adapters"]["basic_qc_evidence"])
    allowed_diag = set(config["adapters"]["diagnostic_evidence"])
    for column, allowed in [
        ("input_adapter", allowed_raw),
        ("basic_qc_evidence_adapter", allowed_basic),
        ("diagnostic_evidence_adapter", allowed_diag),
    ]:
        invalid = sorted(set(manifest[column]) - allowed)
        if invalid:
            raise ValueError(f"Unsupported {column}: {invalid}")
    return manifest


def global_cell_ids(config: dict[str, Any], library_id: str, barcodes: np.ndarray) -> np.ndarray:
    template = str(config["global_cell_id_template"])
    try:
        values = np.asarray(
            [template.format(study_id=config["study_id"], library_id=library_id, barcode=b) for b in barcodes],
            dtype=object,
        )
    except KeyError as exc:
        raise ValueError(f"Unknown global-cell template key: {exc}") from exc
    if len(np.unique(values)) != len(values):
        raise ValueError("Duplicate global cell IDs")
    return values


def _decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values], dtype=object)


def _mito_mask(ids: np.ndarray, names: np.ndarray) -> np.ndarray:
    def normalise(value: str) -> str:
        text = str(value).strip()
        if text.lower().startswith(("mt-", "mt_", "mt.")):
            text = text[3:]
        return "".join(c for c in text.upper() if c.isalnum())
    ids_s = pd.Series(ids, dtype=str)
    names_s = pd.Series(names, dtype=str)
    prefixes = ("mt-", "MT-", "Mt-", "KEG11_")
    return (
        ids_s.str.startswith(prefixes) | names_s.str.startswith(prefixes)
        | ids_s.map(normalise).isin(BARE_MITO) | names_s.map(normalise).isin(BARE_MITO)
    ).to_numpy(bool)


def _metrics_from_matrix(counts: sp.csr_matrix, ids: np.ndarray, names: np.ndarray) -> dict[str, Any]:
    if counts.data.dtype.kind not in "iu" or (len(counts.data) and counts.data.min() < 0):
        raise ValueError("Raw count matrix must contain non-negative integers")
    total = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    detected = np.diff(counts.indptr).astype(np.int64)
    mt = _mito_mask(ids, names)
    mt_counts = np.asarray(counts[:, mt].sum(axis=1)).ravel() if mt.any() else np.zeros(len(total))
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.divide(mt_counts * 100.0, total, out=np.zeros(len(total)), where=total > 0)
        complexity = np.divide(
            np.log10(np.maximum(detected, 1)), np.log10(np.maximum(total, 10)),
            out=np.zeros(len(total)), where=total > 0,
        )
    return {
        "counts": counts,
        "rna_total_umi": total,
        "rna_detected_genes": detected,
        "rna_pct_mito": pct,
        "rna_complexity": complexity,
        "n_empty_called_cells": int((total == 0).sum()),
    }


def read_filtered_h5(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with h5py.File(path, "r") as handle:
        if "matrix" not in handle:
            raise ValueError("Missing /matrix in filtered H5")
        matrix = handle["matrix"]
        required = {"data", "indices", "indptr", "shape", "barcodes", "features"}
        if required - set(matrix.keys()):
            raise ValueError("Malformed Cell Ranger filtered H5")
        shape = tuple(int(x) for x in matrix["shape"][:])
        barcodes = _decode(matrix["barcodes"][:])
        features = matrix["features"]
        ids = _decode(features["id"][:])
        names = _decode(features["name"][:]) if "name" in features else ids
        types = _decode(features["feature_type"][:]) if "feature_type" in features else np.repeat("Gene Expression", shape[0])
        gene = types == "Gene Expression"
        source = sp.csc_matrix((matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]), shape=shape)
        counts = source[gene, :].T.tocsr()
    result = {
        "barcodes": barcodes,
        "counts_dtype": str(counts.dtype),
        "integer_counts": True,
        "n_features_total": shape[0],
        "n_gene_expression_features": int(gene.sum()),
        "n_called_cells": shape[1],
    }
    result.update(_metrics_from_matrix(counts, ids[gene], names[gene]))
    return result


def _first_existing(directory: Path, names: list[str]) -> Path:
    for name in names:
        path = directory / name
        if path.is_file():
            return path
    raise FileNotFoundError(f"Missing required adapter member in {directory}: {names}")


def read_filtered_mtx(path: Path) -> dict[str, Any]:
    if not path.is_dir():
        raise FileNotFoundError(path)
    matrix_path = _first_existing(path, ["matrix.mtx.gz", "matrix.mtx"])
    barcode_path = _first_existing(path, ["barcodes.tsv.gz", "barcodes.tsv"])
    feature_path = _first_existing(path, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
    raw = scipy.io.mmread(matrix_path).tocsr()
    if raw.data.dtype.kind not in "iu":
        rounded = np.rint(raw.data)
        if not np.array_equal(raw.data, rounded):
            raise ValueError("Matrix Market input is not integer-valued")
        raw.data = rounded.astype(np.int64)
    opener = gzip.open if barcode_path.suffix == ".gz" else open
    with opener(barcode_path, "rt", encoding="utf-8") as handle:
        barcodes = np.asarray([line.rstrip("\n").split("\t")[0] for line in handle if line.strip()], dtype=object)
    features = pd.read_csv(feature_path, sep="\t", header=None, dtype=str)
    ids = features.iloc[:, 0].to_numpy(object)
    names = features.iloc[:, 1].to_numpy(object) if features.shape[1] > 1 else ids
    types = features.iloc[:, 2].to_numpy(object) if features.shape[1] > 2 else np.repeat("Gene Expression", len(ids))
    gene = types == "Gene Expression"
    if raw.shape != (len(ids), len(barcodes)):
        raise ValueError("Matrix Market dimensions disagree with features/barcodes")
    if len(np.unique(barcodes)) != len(barcodes):
        raise ValueError("Duplicate barcodes")
    counts = raw[gene, :].T.tocsr()
    result = {
        "barcodes": barcodes,
        "counts_dtype": str(counts.dtype),
        "integer_counts": True,
        "n_features_total": len(ids),
        "n_gene_expression_features": int(gene.sum()),
        "n_called_cells": len(barcodes),
    }
    result.update(_metrics_from_matrix(counts, ids[gene], names[gene]))
    return result


def read_input(adapter: str, path: Path) -> dict[str, Any]:
    if adapter == "cellranger_filtered_h5":
        return read_filtered_h5(path)
    if adapter == "cellranger_filtered_mtx":
        return read_filtered_mtx(path)
    raise ValueError(f"Unsupported raw input adapter: {adapter}")


def mad(values: np.ndarray) -> float:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    return float(np.median(np.abs(values - np.median(values)))) if len(values) else math.nan


def compute_thresholds(metrics: pd.DataFrame, branch: str, rules: dict[str, Any]) -> dict[str, Any]:
    rule = rules[branch]
    log_umi = np.log10(np.maximum(metrics["rna_total_umi"].to_numpy(), 1))
    log_genes = np.log10(np.maximum(metrics["rna_detected_genes"].to_numpy(), 1))
    mito = metrics["rna_pct_mito"].to_numpy(float)
    result = {
        "branch": branch,
        "hard_min_umi": float(rule["hard_min_umi"]),
        "hard_min_genes": float(rule["hard_min_genes"]),
        "mad_n": float(rule["mad_n"]),
        "log10_umi_median": float(np.median(log_umi)),
        "log10_umi_mad": mad(log_umi),
        "log10_genes_median": float(np.median(log_genes)),
        "log10_genes_mad": mad(log_genes),
        "mito_median_pct": float(np.median(mito)),
        "mito_mad_pct": mad(mito),
        "absolute_mito_cap_pct": float(rule["mito_cap_pct"]),
        "complexity_min": float(rule["complexity_min"]),
    }
    result["adaptive_min_log10_umi"] = result["log10_umi_median"] - result["mad_n"] * result["log10_umi_mad"]
    result["adaptive_min_umi"] = float(10 ** result["adaptive_min_log10_umi"])
    result["adaptive_min_log10_genes"] = result["log10_genes_median"] - result["mad_n"] * result["log10_genes_mad"]
    result["adaptive_min_genes"] = float(10 ** result["adaptive_min_log10_genes"])
    result["adaptive_mito_max_pct"] = result["mito_median_pct"] + result["mad_n"] * result["mito_mad_pct"]
    result["effective_mito_max_pct"] = min(result["absolute_mito_cap_pct"], result["adaptive_mito_max_pct"])
    return result


def apply_flags(metrics: pd.DataFrame, threshold: dict[str, Any]) -> pd.DataFrame:
    log_umi = np.log10(np.maximum(metrics["rna_total_umi"].to_numpy(), 1))
    log_genes = np.log10(np.maximum(metrics["rna_detected_genes"].to_numpy(), 1))
    gates = {
        "pass_hard_min_umi": metrics["rna_total_umi"].to_numpy() >= threshold["hard_min_umi"],
        "pass_hard_min_genes": metrics["rna_detected_genes"].to_numpy() >= threshold["hard_min_genes"],
        "pass_adaptive_min_umi": log_umi >= threshold["adaptive_min_log10_umi"],
        "pass_adaptive_min_genes": log_genes >= threshold["adaptive_min_log10_genes"],
        "pass_mito": metrics["rna_pct_mito"].to_numpy() <= threshold["effective_mito_max_pct"],
        "pass_complexity": metrics["rna_complexity"].to_numpy() >= threshold["complexity_min"],
    }
    result = pd.DataFrame(gates)
    result["rna_qc_pass_pre_doublet"] = result[list(gates)].all(axis=1)
    labels = dict(zip(gates, ["low_umi_hard_floor", "low_genes_hard_floor", "low_umi_adaptive", "low_genes_adaptive", "high_mito", "low_complexity"], strict=True))
    result["rna_qc_fail_reason"] = [";".join(labels[c] for c, v in zip(gates, row, strict=True) if not v) for row in result[list(gates)].itertuples(index=False, name=None)]
    return result


def assign_tiers(metrics: pd.DataFrame, flags: dict[str, pd.DataFrame], config: dict[str, Any]) -> tuple[pd.DataFrame, str]:
    candidate = flags["candidate"]
    retention = flags["retention_sensitivity"]
    n_candidate = int(candidate["rna_qc_pass_pre_doublet"].sum())
    primary_min = int(config["sample_gates"]["primary_min_candidate_basic_qc"])
    sensitivity_min = int(config["sample_gates"]["sensitivity_min_candidate_basic_qc"])
    gate = "primary_eligible" if n_candidate >= primary_min else "sensitivity_only_50_99" if n_candidate >= sensitivity_min else "insufficient_primary_cells_lt50"
    primary = candidate["rna_qc_pass_pre_doublet"].to_numpy() & (gate == "primary_eligible")
    broad = retention["rna_qc_pass_pre_doublet"].to_numpy()
    tier = np.select([primary, (~primary) & broad], ["primary_pass", "sensitivity_only"], default="technical_exclude")
    identity = [
        "global_cell_id", "barcode", "study_id", "species", "dataset_id",
        "library_id", "sample_id", "sample_accession", "sample_title",
        "experimental_time", "condition", "author_state",
    ]
    out = metrics[identity].copy()
    for branch in ("strict", "candidate", "retention_sensitivity"):
        prefix = branch.replace("candidate", "candidate_basic_qc").replace("retention_sensitivity", "retention_sensitivity_basic_qc").replace("strict", "strict_basic_qc")
        out[f"{prefix}_pass"] = flags[branch]["rna_qc_pass_pre_doublet"].to_numpy()
        out[f"{prefix}_fail_reason"] = flags[branch]["rna_qc_fail_reason"].to_numpy()
    out["sample_gate"] = gate
    out["pre_doublet_qc_tier"] = tier
    out["pre_doublet_qc_reason"] = np.where(tier == "primary_pass", "candidate_basic_qc_pass_sample_gate_pass", np.where(tier == "sensitivity_only", "retention_sensitivity_basic_qc_pass", "retention_sensitivity_basic_qc_fail"))
    return out, gate


def run_scrublet(counts: sp.csr_matrix, settings: dict[str, Any], seed: int) -> dict[str, Any]:
    n_pcs = max(1, min(int(settings["n_prin_comps"]), counts.shape[0] - 1, counts.shape[1] - 1))
    model = scr.Scrublet(counts, expected_doublet_rate=float(settings["expected_doublet_rate"]), sim_doublet_ratio=float(settings["sim_doublet_ratio"]), random_state=int(seed))
    scores, predicted = model.scrub_doublets(
        synthetic_doublet_umi_subsampling=float(settings["synthetic_doublet_umi_subsampling"]),
        min_counts=int(settings["min_counts"]), min_cells=int(settings["min_cells"]),
        min_gene_variability_pctl=float(settings["min_gene_variability_pctl"]),
        n_prin_comps=n_pcs, verbose=False,
    )
    threshold = float(model.threshold_) if getattr(model, "threshold_", None) is not None else math.nan
    return {"seed": int(seed), "scores": np.asarray(scores), "predicted": np.asarray(predicted, bool), "simulated": np.asarray(model.doublet_scores_sim_), "threshold": threshold, "n_prin_comps": n_pcs}
