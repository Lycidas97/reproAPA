"""Summarize route-level WPP effects without changing the frozen estimator."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def metrics(values: pd.Series) -> dict:
    x = np.sort(values.to_numpy(float))
    cut = int(np.floor(len(x) * 0.1))
    trimmed = x[cut : len(x) - cut]
    result = {
        "genes": len(x),
        "mean_delta_WPP": float(x.mean()),
        "median_delta_WPP": float(np.median(x)),
        "trimmed_mean_10_each_tail": float(trimmed.mean()),
    }
    for threshold in [0.02, 0.05]:
        label = str(threshold).replace(".", "")
        result[f"proximal_{label}"] = int((x <= -threshold).sum())
        result[f"distal_{label}"] = int((x >= threshold).sum())
        result[f"subthreshold_{label}"] = int((np.abs(x) < threshold).sum())
    return result


def main(input_path: Path, output: Path, provisional: bool) -> None:
    route = pd.read_csv(input_path, sep="\t")
    rows = []
    keys = ["species", "study_id", "route", "view", "mode"]
    for key, frame in route.groupby(keys):
        rows.append(dict(zip(keys, key)) | metrics(frame.delta))
    summary = pd.DataFrame(rows)
    summary.to_csv(output, sep="\t", index=False)
    audit = {
        "status": "PROVISIONAL_EXCLUDED_LIBRARY" if provisional else "PASS",
        "rows": len(summary),
        "routes": int(summary.route.nunique()),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
    }
    output.with_suffix(".audit.json").write_text(json.dumps(audit, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--provisional", action="store_true")
    args = parser.parse_args()
    main(args.input, args.output, args.provisional)
