#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from threadpoolctl import threadpool_limits

from cluster_utils import graph_from_coordinates, leiden, stratified_indices, write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    study = yaml.safe_load(Path(args.study_config).read_text(encoding="utf-8"))
    project = yaml.safe_load(Path(args.project_config).read_text(encoding="utf-8"))
    pipe = study["scopes"][args.scope]["pipeline"]
    selection = project["selection"]
    resolutions = np.asarray(sorted(set(float(x) for x in selection["resolutions"])), dtype=np.float32)
    root = Path(args.output_root) / "scopes" / args.scope
    target = root / "sweep"
    target.mkdir(parents=True, exist_ok=True)
    coordinates = np.load(root / "embedding" / "pca.npy", mmap_mode="r")[:, :int(pipe["pca_n_components_used"])]
    cohort = pd.read_csv(root / "cohort" / "evaluation_cells.tsv.gz", sep="\t", dtype={"full_row_index": int, "cell_id": str, "library_id": str})
    eval_coords = np.asarray(coordinates[cohort["full_row_index"].to_numpy(int)], dtype=np.float32)
    n_resamples = int(selection["n_resamples"])
    fraction = float(selection["subsample_fraction"])
    seed = int(selection["seed"])
    n_cells = len(cohort)
    full_labels = np.full((len(resolutions), n_cells), -1, dtype=np.int16)
    sampled_labels = np.full((len(resolutions), n_resamples, n_cells), -1, dtype=np.int16)
    manifest_rows = []
    implementation = str(pipe["leiden_implementation"])
    with threadpool_limits(limits=args.threads):
        full_graph = graph_from_coordinates(eval_coords, int(pipe["neighbor_n_neighbors"]), seed, args.threads, str(pipe["neighbor_backend"]))
        for r_index, resolution in enumerate(resolutions):
            full_labels[r_index] = leiden(full_graph, float(resolution), seed, int(pipe["leiden_iterations"]), implementation)
        for replicate in range(n_resamples):
            replicate_seed = seed + replicate + 1
            sample_size = int(np.ceil(fraction * n_cells))
            chosen = stratified_indices(cohort["library_id"], sample_size, replicate_seed)
            graph = graph_from_coordinates(eval_coords[chosen], int(pipe["neighbor_n_neighbors"]), replicate_seed,
                                           args.threads, str(pipe["neighbor_backend"]))
            for r_index, resolution in enumerate(resolutions):
                sampled_labels[r_index, replicate, chosen] = leiden(
                    graph, float(resolution), replicate_seed, int(pipe["leiden_iterations"]), implementation,
                )
            manifest_rows.append({"replicate": replicate + 1, "seed": replicate_seed, "sampled_n": len(chosen),
                                  "sampled_fraction": len(chosen) / n_cells})
    np.savez_compressed(target / "partitions.npz", resolutions=resolutions, full_labels=full_labels,
                        sampled_labels=sampled_labels)
    pd.DataFrame(manifest_rows).to_csv(target / "resample_manifest.tsv", sep="\t", index=False)
    pd.DataFrame({"cell_id": cohort["cell_id"]}).join(
        pd.DataFrame({f"res_{float(resolution):g}": full_labels[index] for index, resolution in enumerate(resolutions)})
    ).to_csv(target / "full_partitions.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    write_json(target / "metadata.json", {"study_id": study["study_id"], "scope_id": args.scope,
                                           "n_cells": n_cells, "n_resamples": n_resamples,
                                           "subsample_fraction": fraction, "resolutions": resolutions.tolist()})
    write_marker(target / ".ok", status="PASS", scope_id=args.scope, n_resamples=n_resamples, n_resolutions=len(resolutions))


if __name__ == "__main__":
    main()
