#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-h5ad", type=Path, required=True)
    parser.add_argument("--reference-h5ad", type=Path, required=True)
    parser.add_argument("--reference-model", type=Path, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    torch.set_num_threads(args.threads)
    query = ad.read_h5ad(args.query_h5ad)
    reference = ad.read_h5ad(args.reference_h5ad)
    if not query.var_names.equals(reference.var_names):
        raise ValueError("query genes do not exactly match the frozen reference")
    model = scvi.model.SCVI.load(args.reference_model, adata=reference, accelerator="cpu")
    latent = model.get_latent_representation(query, batch_size=2048).astype(np.float32)
    if latent.shape != (query.n_obs, reference.obsm["X_scVI"].shape[1]):
        raise ValueError(f"unexpected latent shape: {latent.shape}")
    if not np.isfinite(latent).all():
        raise ValueError("projected latent contains non-finite values")

    reference_latent = pd.DataFrame(
        np.asarray(reference.obsm["X_scVI"]), index=reference.obs.cell_id.astype(str)
    )
    overlap = query.obs.cell_id.astype(str).isin(reference_latent.index).to_numpy()
    max_abs = None
    if overlap.any():
        expected = reference_latent.loc[query.obs.loc[overlap, "cell_id"].astype(str)].to_numpy()
        max_abs = float(np.max(np.abs(latent[overlap] - expected)))
        latent[overlap] = expected

    columns = [f"scvi_latent_{index + 1:02d}" for index in range(latent.shape[1])]
    output = pd.DataFrame(latent, columns=columns)
    output.insert(0, "library_id", query.obs.library_id.astype(str).to_numpy())
    output.insert(0, "library_key", query.obs.library_key.astype(str).to_numpy())
    output.insert(0, "cell_id", query.obs.cell_id.astype(str).to_numpy())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(
        args.output,
        sep="\t",
        index=False,
        float_format="%.7g",
        compression={"method": "gzip", "mtime": 0},
    )
    report = {
        "status": "PASS",
        "n_cells": query.n_obs,
        "n_latent": latent.shape[1],
        "n_reference_overlap": int(overlap.sum()),
        "reference_overlap_max_abs_before_replacement": max_abs,
        "frozen_reference_encoder": True,
        "query_model_finetuned": False,
    }
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
