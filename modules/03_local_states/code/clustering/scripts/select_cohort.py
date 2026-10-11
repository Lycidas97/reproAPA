#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from cluster_utils import ordered_sha256, stratified_indices, write_json, write_marker


def choose_evaluation_size(
    full_n: int,
    minimum_cells: int,
    target_fraction: float,
    maximum_cells: int,
) -> int:
    if full_n <= 0:
        raise ValueError("full_n must be positive")
    if minimum_cells <= 0 or maximum_cells <= 0 or minimum_cells > maximum_cells:
        raise ValueError("evaluation cell bounds are invalid")
    if not 0 < target_fraction <= 1:
        raise ValueError("target_fraction must be in (0, 1]")
    target = max(minimum_cells, math.ceil(full_n * target_fraction))
    return min(full_n, maximum_cells, target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", required=True)
    parser.add_argument("--max-cells", type=int)
    parser.add_argument("--minimum-cells", type=int)
    parser.add_argument("--target-fraction", type=float)
    parser.add_argument("--maximum-cells", type=int)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    cells = pd.read_csv(args.cells, sep="\t", dtype=str)
    adaptive = (args.minimum_cells, args.target_fraction, args.maximum_cells)
    if args.max_cells is not None and any(value is not None for value in adaptive):
        raise ValueError("legacy --max-cells cannot be combined with adaptive evaluation bounds")
    if args.max_cells is not None:
        evaluation_n = min(args.max_cells, len(cells))
        policy = "fixed_maximum"
    elif all(value is not None for value in adaptive):
        evaluation_n = choose_evaluation_size(
            len(cells), args.minimum_cells, args.target_fraction, args.maximum_cells
        )
        policy = "bounded_fraction"
    else:
        raise ValueError(
            "provide --max-cells or all of --minimum-cells, --target-fraction and --maximum-cells"
        )
    selected = stratified_indices(cells["library_id"], evaluation_n, args.seed)
    cohort = cells.iloc[selected].copy()
    cohort.insert(0, "full_row_index", selected)
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    cohort.to_csv(target / "evaluation_cells.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    strata = cohort.groupby("library_id", sort=True).size().rename("evaluation_n").reset_index()
    strata = strata.merge(cells.groupby("library_id", sort=True).size().rename("full_n").reset_index(), on="library_id")
    strata["evaluation_fraction"] = strata["evaluation_n"] / strata["full_n"]
    strata.to_csv(target / "strata.tsv", sep="\t", index=False)
    metadata = {
        "full_n": len(cells),
        "evaluation_n": len(cohort),
        "evaluation_fraction": len(cohort) / len(cells),
        "evaluation_size_policy": policy,
        "seed": args.seed,
        "cell_order_sha256": ordered_sha256(cohort["cell_id"]),
    }
    if policy == "bounded_fraction":
        metadata.update({
            "minimum_cells": args.minimum_cells,
            "target_fraction": args.target_fraction,
            "maximum_cells": args.maximum_cells,
        })
    else:
        metadata["maximum_cells"] = args.max_cells
    write_json(target / "metadata.json", metadata)
    write_marker(target / ".ok", status="PASS", full_n=len(cells), evaluation_n=len(cohort))


if __name__ == "__main__":
    main()
