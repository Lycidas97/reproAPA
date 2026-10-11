#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from cluster_utils import write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", required=True)
    parser.add_argument("--shards", nargs="+", required=True)
    parser.add_argument("--manifests", nargs="+", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    project = yaml.safe_load(Path(args.project_config).read_text(encoding="utf-8"))
    n_resamples = int(project["selection"]["n_resamples"])
    with np.load(args.full, allow_pickle=False) as archive:
        resolutions = archive["resolutions"].astype(np.float32)
        full_labels = archive["full_labels"].astype(np.int16)
    sampled = np.full((len(resolutions), n_resamples, full_labels.shape[1]), -1, dtype=np.int16)
    seen = np.zeros(n_resamples, dtype=np.int8)
    for path in args.shards:
        with np.load(path, allow_pickle=False) as archive:
            if not np.array_equal(archive["resolutions"].astype(np.float32), resolutions):
                raise ValueError(f"resolution drift in shard {path}")
            indices = archive["replicate_indices"].astype(int)
            values = archive["sampled_labels"].astype(np.int16)
        if values.shape != (len(resolutions), len(indices), full_labels.shape[1]):
            raise ValueError(f"shape drift in shard {path}")
        sampled[:, indices] = values
        seen[indices] += 1
    if not np.all(seen == 1):
        raise ValueError("merged shards do not cover each replicate exactly once")
    manifests = pd.concat([pd.read_csv(path, sep="\t") for path in args.manifests], ignore_index=True)
    manifests = manifests.sort_values("replicate_index", kind="stable")
    if manifests["replicate_index"].tolist() != list(range(n_resamples)):
        raise ValueError("manifest replicate order/coverage drift")
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(target / "partitions.npz", resolutions=resolutions, full_labels=full_labels,
                        sampled_labels=sampled)
    manifests.to_csv(target / "resample_manifest.tsv", sep="\t", index=False)
    write_json(target / "metadata.json", {"n_cells": full_labels.shape[1], "n_resamples": n_resamples,
                                           "n_shards": len(args.shards),
                                           "resolutions": resolutions.astype(float).tolist(),
                                           "replicate_coverage_status": "PASS"})
    write_marker(target / ".ok", status="PASS", n_resamples=n_resamples, n_shards=len(args.shards))


if __name__ == "__main__":
    main()
