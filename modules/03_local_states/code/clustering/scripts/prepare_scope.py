#!/usr/bin/env python3
"""Prepare a frozen scope with explicit source, strata, and transform contracts."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse


FORBIDDEN_FIELDS = {
    "time", "stage", "marker", "program_score", "author_label", "route_outcome",
    "pas", "apa", "cross_study_state", "cross_study_mapping",
}


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values])


def ordered_sha256(values) -> str:
    return hashlib.sha256("".join(f"{value}\n" for value in values).encode("utf-8")).hexdigest()


def set_sha256(values) -> str:
    return ordered_sha256(sorted((str(value) for value in values), key=lambda value: value.encode("utf-8")))


def uses_legacy_gse171_frozen_contract(config: dict, scope_cfg: dict) -> bool:
    return config["study_id"] == "GSE171127" and bool(scope_cfg.get("frozen_execution_object"))


def infer_stratum(study_id: str, scope: str, cell_id: str) -> str:
    if study_id == "GSE118258":
        fields = str(cell_id).split("::")
        if len(fields) < 3 or not fields[1]:
            raise ValueError(f"cannot recover GSE118258 day/library from global_cell_id={cell_id}")
        return f"day_{fields[1]}"
    fields = str(cell_id).split("|")
    if len(fields) >= 3 and fields[1]:
        return fields[1]
    raise ValueError(f"no non-NA stratum rule for {study_id}/{scope}")


def load_membership(config: dict, scope: str, scope_cfg: dict) -> pd.DataFrame:
    mapping = config["membership_columns"]
    columns = [mapping["cell_id"]]
    if mapping.get("state_id"):
        columns.append(mapping["state_id"])
    for key in ["scope_id", "library_id", "barcode"]:
        if mapping.get(key):
            columns.append(mapping[key])
    frame = pd.read_csv(config["membership_path"], sep="\t", dtype=str, usecols=list(dict.fromkeys(columns)))
    if mapping.get("scope_id"):
        values = scope_cfg.get("membership_scope_values", [scope])
        frame = frame.loc[frame[mapping["scope_id"]].isin(values)].copy()
    rename = {mapping["cell_id"]: "cell_id"}
    if mapping.get("state_id"):
        rename[mapping["state_id"]] = "source_state_id"
    if mapping.get("scope_id"):
        rename[mapping["scope_id"]] = "source_scope_id"
    if mapping.get("library_id"):
        rename[mapping["library_id"]] = "library_id"
    if mapping.get("barcode"):
        rename[mapping["barcode"]] = "barcode"
    frame = frame.rename(columns=rename)
    if scope_cfg.get("deduplicate_membership_by_cell_id", False):
        frame = frame.drop_duplicates("cell_id", keep="first").copy()
    if "source_scope_id" not in frame:
        frame["source_scope_id"] = scope
    if "source_state_id" not in frame:
        frame["source_state_id"] = "NOT_USED"
    if "library_id" not in frame:
        frame["library_id"] = [infer_stratum(config["study_id"], scope, value) for value in frame["cell_id"]]
    if frame["library_id"].isna().any() or frame["library_id"].astype(str).str.strip().isin({"", "NA", "nan"}).any():
        raise ValueError("bootstrap stratum contains NA")
    if "barcode" not in frame:
        frame["barcode"] = frame["cell_id"].str.rsplit("|", n=1).str[-1].str.rsplit("::", n=1).str[-1]
    frame["scope_id"] = scope
    frame["study_id"] = config["study_id"]
    if frame["cell_id"].duplicated().any() or frame["cell_id"].isna().any():
        raise ValueError("invalid frozen membership")
    output_columns = ["study_id", "scope_id", "source_scope_id", "cell_id", "source_state_id", "library_id", "barcode"]
    metadata_path = config.get("library_metadata_path")
    if metadata_path:
        metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
        safe_fields = [
            field for field in ["species", "sample_accession", "sample_title", "experimental_time", "condition"]
            if field in metadata
        ]
        library_metadata = metadata[["library_id", *safe_fields]].drop_duplicates()
        if library_metadata["library_id"].duplicated().any():
            raise ValueError("library metadata is not one-to-one by library_id")
        frame = frame.merge(library_metadata, on="library_id", how="left", validate="many_to_one")
        if frame[safe_fields].isna().any().any():
            raise ValueError("membership does not close against library metadata")
        output_columns.extend(safe_fields)
    return frame[output_columns]


def normalize_log1p(raw: sparse.csr_matrix, target_sum: float) -> sparse.csr_matrix:
    totals = np.asarray(raw.sum(axis=1)).ravel().astype(np.float64)
    if np.any(totals <= 0):
        raise ValueError("frozen scope contains a zero-count GEX cell")
    matrix = sparse.diags((target_sum / totals).astype(np.float32)).dot(raw.astype(np.float32)).tocsr()
    matrix.data = np.log1p(matrix.data).astype(np.float32, copy=False)
    return matrix


def load_10x(path: str, membership: pd.DataFrame) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        feature_type = decode(group["features/feature_type"][:])
        keep = feature_type == "Gene Expression"
        raw = sparse.csc_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=tuple(int(v) for v in group["shape"][:])).T.tocsr()
        barcodes = decode(group["barcodes"][:])
        gene_ids = decode(group["features/id"][:])[keep]
        symbols = decode(group["features/name"][:])[keep]
    selected = pd.Index(barcodes).get_indexer(membership["barcode"])
    if np.any(selected < 0) or len(np.unique(selected)) != len(selected):
        raise ValueError("membership does not map one-to-one to 10x barcodes")
    return raw[selected][:, keep].tocsr(), pd.DataFrame({"gene_id": gene_ids, "gene_symbol": symbols})


def load_h5ad(path: str, membership: pd.DataFrame) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    obj = ad.read_h5ad(path)
    selected = pd.Index(obj.obs_names.astype(str)).get_indexer(membership["cell_id"])
    if np.any(selected < 0) or len(np.unique(selected)) != len(selected):
        raise ValueError("membership does not map one-to-one to H5AD obs_names")
    matrix = sparse.csr_matrix(obj.X[selected])
    symbols = obj.var.get("gene_symbol", pd.Series(obj.var_names.astype(str), index=obj.var_names)).astype(str).to_numpy()
    gene_ids = obj.var.get("gene_id", pd.Series(obj.var_names.astype(str), index=obj.var_names)).astype(str).to_numpy()
    return matrix, pd.DataFrame({"gene_id": gene_ids, "gene_symbol": symbols})


def load_wot(config: dict, membership: pd.DataFrame) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    cells = pd.read_parquet(config["expression_cells_path"], columns=["cell_index", "cell_id"])
    joined = membership[["cell_id"]].merge(cells, on="cell_id", how="left", validate="one_to_one")
    if joined["cell_index"].isna().any():
        raise ValueError("WOT membership does not close against cached cells")
    matrix = sparse.load_npz(config["expression_path"]).tocsr()[joined["cell_index"].to_numpy(dtype=int)]
    features = pd.read_parquet(config["expression_features_path"])[["gene_id", "gene_symbol"]].reset_index(drop=True)
    return matrix, features


def align_and_stack_sources(chunks: list[tuple[np.ndarray, sparse.csr_matrix, pd.DataFrame]], n_cells: int):
    if not chunks:
        raise ValueError("no expression source matched the study-wide membership")
    for _, _, features in chunks:
        if features["gene_id"].duplicated().any():
            raise ValueError("expression source contains duplicate gene IDs")
    common = set(chunks[0][2]["gene_id"].astype(str))
    for _, _, features in chunks[1:]:
        common &= set(features["gene_id"].astype(str))
    first_features = chunks[0][2]
    keep_ids = [gene_id for gene_id in first_features["gene_id"].astype(str) if gene_id in common]
    if not keep_ids:
        raise ValueError("expression sources have no common genes")
    matrices = []
    positions = []
    for row_positions, matrix, features in chunks:
        indexer = pd.Index(features["gene_id"].astype(str)).get_indexer(keep_ids)
        if np.any(indexer < 0):
            raise ValueError("common-gene alignment failed")
        matrices.append(matrix[:, indexer].tocsr())
        positions.append(row_positions)
    stacked = sparse.vstack(matrices, format="csr")
    stacked_positions = np.concatenate(positions)
    if len(stacked_positions) != n_cells or len(np.unique(stacked_positions)) != n_cells:
        raise ValueError("expression sources do not cover each membership cell exactly once")
    reorder = np.argsort(stacked_positions, kind="stable")
    if not np.array_equal(stacked_positions[reorder], np.arange(n_cells)):
        raise ValueError("expression source row contract is incomplete")
    features = first_features.set_index(first_features["gene_id"].astype(str)).loc[keep_ids].reset_index(drop=True)
    return stacked[reorder].tocsr(), features[["gene_id", "gene_symbol"]]


def load_multiple_sources(scope_cfg: dict, membership: pd.DataFrame):
    chunks = []
    assigned = np.zeros(len(membership), dtype=np.int8)
    for source in scope_cfg["expression_sources"]:
        selector = source.get("selector", "source_scope_id")
        value = str(source["value"])
        positions = np.flatnonzero(membership[selector].astype(str).eq(value).to_numpy())
        if not len(positions):
            raise ValueError(f"expression source {selector}={value} matched no cells")
        assigned[positions] += 1
        selected = membership.iloc[positions]
        if source["mode"] == "10x_h5":
            matrix, features = load_10x(source["path"], selected)
        elif source["mode"] == "h5ad":
            matrix, features = load_h5ad(source["path"], selected)
        else:
            raise ValueError(f"unsupported expression source mode: {source['mode']}")
        chunks.append((positions, matrix, features))
    if np.any(assigned != 1):
        raise ValueError("each study-wide membership cell must match exactly one expression source")
    return align_and_stack_sources(chunks, len(membership))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    if args.scope not in config["scopes"] or config["study_id"] == "GSE225142":
        raise ValueError("illegal scope/study")
    if FORBIDDEN_FIELDS & set(config.get("partition_selection_fields", [])):
        raise ValueError("forbidden biological field requested")
    scope_cfg = config["scopes"][args.scope]
    pipe = scope_cfg["pipeline"]
    membership = load_membership(config, args.scope, scope_cfg)
    source_membership_order_sha256 = ordered_sha256(membership["cell_id"])
    membership_set_sha256 = set_sha256(membership["cell_id"])
    if len(membership) != int(scope_cfg["n_cells"]):
        raise ValueError("scope membership count drift")
    if source_membership_order_sha256 != scope_cfg["source_membership_order_sha256"]:
        raise ValueError("source membership order hash drift")
    if membership_set_sha256 != scope_cfg["cell_set_sha256"]:
        raise ValueError("scope membership cell-set hash drift")
    frozen_object_sha256 = "NOT_APPLICABLE"
    legacy_gse171_contract = uses_legacy_gse171_frozen_contract(config, scope_cfg)
    if legacy_gse171_contract:
        frozen_path = Path(scope_cfg["frozen_execution_object"])
        with np.load(frozen_path, allow_pickle=True) as frozen:
            execution_order = frozen["global_cell_id"].astype(str)
        if set(execution_order) != set(membership["cell_id"]):
            raise ValueError("frozen executable cell set and accepted membership differ")
        execution_contract = config.get("gse171_execution_contract", "frozen_count_object_order")
        if execution_contract == "frozen_count_object_order":
            membership = membership.set_index("cell_id").loc[execution_order].reset_index()
        elif execution_contract == "historical_source_membership_order_compatibility_only":
            if config.get("execution_mode") != "builder_smoke":
                raise ValueError("historical GSE171 order is restricted to builder_smoke compatibility evidence")
        else:
            raise ValueError(f"unknown GSE171 execution contract: {execution_contract}")
        digest = hashlib.sha256()
        with frozen_path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        frozen_object_sha256 = digest.hexdigest()
    execution_cell_order_sha256 = ordered_sha256(membership["cell_id"])
    mode = config["expression_mode"]
    if scope_cfg.get("expression_sources"):
        raw, features = load_multiple_sources(scope_cfg, membership)
    elif mode == "scope_10x_h5":
        raw, features = load_10x(scope_cfg["expression_path"], membership)
    elif mode in {"scope_h5ad", "shared_h5ad"}:
        raw, features = load_h5ad(scope_cfg["expression_path"], membership)
    elif mode == "wot_log1p_cache":
        raw, features = load_wot(config, membership)
    else:
        raise ValueError(f"unsupported expression mode: {mode}")
    observed_filter = 0
    if int(pipe.get("filter_genes_min_cells", 0)):
        minimum = int(pipe["filter_genes_min_cells"])
        keep = np.asarray((raw > 0).sum(axis=0)).ravel() >= minimum
        observed_filter = minimum
        raw, features = raw[:, keep].tocsr(), features.loc[keep].reset_index(drop=True)
    already_log1p = bool(config.get("expression_already_log1p", False))
    if already_log1p:
        matrix = raw
    elif legacy_gse171_contract:
        # Historical SRG05 replay uses Scanpy on the frozen sparse count matrix.
        # Calling the same APIs preserves their sparse dtype/rounding behavior.
        normalized = ad.AnnData(X=raw.copy())
        sc.pp.normalize_total(normalized, target_sum=float(pipe["normalization_target_sum"]))
        sc.pp.log1p(normalized)
        matrix = sparse.csr_matrix(normalized.X).tocsr()
    else:
        matrix = normalize_log1p(raw, float(pipe["normalization_target_sum"]))
    if matrix.shape != (len(membership), len(features)):
        raise ValueError("prepared matrix dimension mismatch")
    target = Path(args.output_root) / "scopes" / args.scope / "prepared"
    target.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(target / "log1p_gex.npz", matrix.tocsr(), compressed=True)
    membership.to_csv(target / "cells.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    features.to_csv(target / "features.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    pd.DataFrame([
        {"component": "filter_genes_min_cells", "expected": str(pipe.get("filter_genes_min_cells", 0)), "observed": str(observed_filter), "status": "PASS"},
        {"component": "normalization", "expected": str(pipe["normalization_method"]), "observed": "precomputed_log1p" if already_log1p else "normalize_total", "status": "PASS"},
        {"component": "transform", "expected": str(pipe["transform"]), "observed": "existing_log1p" if already_log1p else "log1p", "status": "PASS"},
    ]).to_csv(target / "effective_pipeline.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "study_id": config["study_id"], "scope_id": args.scope,
        "source_membership_order_sha256": source_membership_order_sha256,
        "membership_cell_set_sha256": membership_set_sha256,
        "execution_cell_order_sha256": execution_cell_order_sha256,
        "frozen_execution_object": scope_cfg.get("frozen_execution_object", "NOT_APPLICABLE"),
        "frozen_execution_object_sha256": frozen_object_sha256,
        "execution_contract": config.get("gse171_execution_contract", "source_membership_order"),
        "n_cells": len(membership), "cell_set_status": "PASS", "execution_order_status": "PASS",
    }]).to_csv(target / "cell_order_contract.tsv", sep="\t", index=False)
    metadata = {
        "study_id": config["study_id"], "scope_id": args.scope, "n_cells": len(membership), "n_features": len(features),
        "source_state_count": int(membership["source_state_id"].nunique()),
        "source_states_used_for_clustering": False, "forbidden_biological_fields_read": False,
        "pas_apa_accessed": False, "source_membership_order_sha256": source_membership_order_sha256,
        "membership_cell_set_sha256": membership_set_sha256,
        "execution_cell_order_sha256": execution_cell_order_sha256,
        "frozen_execution_object_sha256": frozen_object_sha256,
        "execution_contract": config.get("gse171_execution_contract", "source_membership_order"),
        "feature_universe_sha256": ordered_sha256(features["gene_id"]), "bootstrap_strata_n": membership["library_id"].nunique(),
        "bootstrap_strata_contains_na": False,
    }
    (target / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary = target / ".ok.tmp"
    temporary.write_text(f"study_id={config['study_id']}\nscope_id={args.scope}\nstatus=PASS\n", encoding="utf-8")
    os.replace(temporary, target / ".ok")


if __name__ == "__main__":
    main()
