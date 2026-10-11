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
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--n-shards", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    study = yaml.safe_load(Path(args.study_config).read_text(encoding="utf-8"))
    project = yaml.safe_load(Path(args.project_config).read_text(encoding="utf-8"))
    pipe = study["scopes"][args.scope]["pipeline"]
    selection = project["selection"]
    n_resamples = int(selection["n_resamples"])
    if not 0 <= args.shard_index < args.n_shards:
        raise ValueError("invalid shard index")
    replicate_indices = np.arange(args.shard_index, n_resamples, args.n_shards, dtype=np.int16)
    if not len(replicate_indices):
        raise ValueError("sweep shard has no assigned replicates")
    resolutions = np.asarray(sorted(set(float(x) for x in selection["resolutions"])), dtype=np.float32)
    root = Path(args.output_root) / "scopes" / args.scope
    target = root / "sweep" / "shards" / f"{args.shard_index:03d}"
    target.mkdir(parents=True, exist_ok=True)
    coordinates = np.load(root / "embedding" / "pca.npy", mmap_mode="r")[:, :int(pipe["pca_n_components_used"])]
    cohort = pd.read_csv(root / "cohort" / "evaluation_cells.tsv.gz", sep="\t",
                         dtype={"full_row_index": int, "cell_id": str, "library_id": str})
    eval_coords = np.asarray(coordinates[cohort["full_row_index"].to_numpy(int)], dtype=np.float32)
    labels = np.full((len(resolutions), len(replicate_indices), len(cohort)), -1, dtype=np.int16)
    manifest = []
    seed = int(selection["seed"])
    sample_size = int(np.ceil(float(selection["subsample_fraction"]) * len(cohort)))
    with threadpool_limits(limits=args.threads):
        for local_index, replicate_index in enumerate(replicate_indices.astype(int)):
            replicate_seed = seed + replicate_index + 1
            chosen = stratified_indices(cohort["library_id"], sample_size, replicate_seed)
            graph = graph_from_coordinates(eval_coords[chosen], int(pipe["neighbor_n_neighbors"]), replicate_seed,
                                           args.threads, str(pipe["neighbor_backend"]))
            for resolution_index, resolution in enumerate(resolutions):
                labels[resolution_index, local_index, chosen] = leiden(
                    graph, float(resolution), replicate_seed, int(pipe["leiden_iterations"]),
                    str(pipe["leiden_implementation"]),
                )
            manifest.append({"replicate": replicate_index + 1, "replicate_index": replicate_index,
                             "seed": replicate_seed, "sampled_n": len(chosen),
                             "sampled_fraction": len(chosen) / len(cohort),
                             "shard_index": args.shard_index})
    np.savez_compressed(target / "partitions.npz", resolutions=resolutions,
                        replicate_indices=replicate_indices, sampled_labels=labels)
    pd.DataFrame(manifest).to_csv(target / "manifest.tsv", sep="\t", index=False)
    write_json(target / "metadata.json", {"shard_index": args.shard_index, "n_shards": args.n_shards,
                                           "replicate_indices": replicate_indices.astype(int).tolist()})
    write_marker(target / ".ok", status="PASS", shard_index=args.shard_index,
                 n_replicates=len(replicate_indices))


if __name__ == "__main__":
    main()
