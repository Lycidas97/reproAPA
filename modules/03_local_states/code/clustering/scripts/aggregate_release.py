#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from cluster_utils import write_json, write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--dataset", default=None)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    root = Path(args.output_root)
    paths = sorted((root / args.dataset / "scopes").glob("*/release/scope_release.json")) if args.dataset else sorted(root.glob("*/release/dataset_release.tsv"))
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    if args.dataset:
        rows = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
        if not rows:
            raise ValueError(f"no scope releases found for {args.dataset}")
        frame = pd.DataFrame(rows)
        frame.insert(0, "dataset", args.dataset)
        frame.to_csv(target / "dataset_release.tsv", sep="\t", index=False)
        status = "REVIEW" if frame["status"].eq("REVIEW").any() else "PASS"
        write_json(target / "dataset_release.json", {"dataset": args.dataset, "n_scopes": len(frame),
                                                      "n_cells": int(frame["n_cells"].sum()), "status": status})
    else:
        frames = [pd.read_csv(path, sep="\t") for path in paths]
        if not frames:
            raise ValueError("no dataset releases found")
        frame = pd.concat(frames, ignore_index=True)
        frame.to_csv(target / "project_release.tsv", sep="\t", index=False)
        status = "REVIEW" if frame["status"].eq("REVIEW").any() else "PASS"
        write_json(target / "project_release.json", {"n_datasets": frame["dataset"].nunique(),
                                                      "n_scopes": len(frame), "n_cells": int(frame["n_cells"].sum()),
                                                      "status": status})
    write_marker(target / ".ok", status=status)


if __name__ == "__main__":
    main()
