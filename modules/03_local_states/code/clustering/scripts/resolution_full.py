#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from threadpoolctl import threadpool_limits

from cluster_utils import graph_from_coordinates, leiden, write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-config", required=True)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    study = yaml.safe_load(Path(args.study_config).read_text(encoding="utf-8"))
    project = yaml.safe_load(Path(args.project_config).read_text(encoding="utf-8"))
    pipe = study["scopes"][args.scope]["pipeline"]
    resolutions = np.asarray(sorted(set(float(x) for x in project["selection"]["resolutions"])), dtype=np.float32)
    root = Path(args.output_root) / "scopes" / args.scope
    target = root / "sweep" / "full"
    target.mkdir(parents=True, exist_ok=True)
    coordinates = np.load(root / "embedding" / "pca.npy", mmap_mode="r")[:, :int(pipe["pca_n_components_used"])]
    cohort = pd.read_csv(root / "cohort" / "evaluation_cells.tsv.gz", sep="\t",
                         dtype={"full_row_index": int, "cell_id": str})
    selected = cohort["full_row_index"].to_numpy(int)
    eval_coords = np.asarray(coordinates[selected], dtype=np.float32)
    seed = int(project["selection"]["seed"])
    labels = np.full((len(resolutions), len(cohort)), -1, dtype=np.int16)
    with threadpool_limits(limits=args.threads):
        graph = graph_from_coordinates(eval_coords, int(pipe["neighbor_n_neighbors"]), seed, args.threads,
                                       str(pipe["neighbor_backend"]))
        for index, resolution in enumerate(resolutions):
            labels[index] = leiden(graph, float(resolution), seed, int(pipe["leiden_iterations"]),
                                   str(pipe["leiden_implementation"]))
    np.savez_compressed(target / "full_labels.npz", resolutions=resolutions, full_labels=labels)
    pd.DataFrame({"cell_id": cohort["cell_id"]}).join(
        pd.DataFrame({f"res_{float(resolution):g}": labels[index]
                      for index, resolution in enumerate(resolutions)})
    ).to_csv(root / "sweep" / "full_partitions.tsv.gz", sep="\t", index=False,
             compression={"method": "gzip", "mtime": 0})
    write_json(target / "metadata.json", {"n_cells": len(cohort), "resolutions": resolutions.tolist(),
                                           "seed": seed})
    write_marker(target / ".ok", status="PASS", n_cells=len(cohort), n_resolutions=len(resolutions))


if __name__ == "__main__":
    main()
