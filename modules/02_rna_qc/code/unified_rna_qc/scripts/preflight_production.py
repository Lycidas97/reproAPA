#!/usr/bin/env python3
"""Validate every production matrix and shared QC contract before fan-out."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py

from common import load_config, read_manifest, resolve, sha256_text


def count_mtx_barcodes(path: Path) -> int:
    for name in ("barcodes.tsv.gz", "barcodes.tsv"):
        candidate = path / name
        if candidate.is_file():
            import gzip
            opener = gzip.open if candidate.suffix == ".gz" else open
            with opener(candidate, "rt", encoding="utf-8") as handle:
                return sum(1 for line in handle if line.strip())
    raise FileNotFoundError(f"Missing barcodes file under {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--marker", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    with args.registry.open(newline="", encoding="utf-8") as handle:
        studies = list(csv.DictReader(handle, delimiter="\t"))
    if len(studies) != 6:
        raise ValueError(f"Expected six study scopes, observed {len(studies)}")

    rows = []
    method_contract = None
    seen_libraries = set()
    total_cells = 0
    for study in studies:
        config_path = resolve(root, study["config_path"])
        config = load_config(config_path)
        manifest = read_manifest(root, config)
        current_contract = {
            "adapters": config["adapters"],
            "basic_qc": config["basic_qc"],
            "sample_gates": config["sample_gates"],
            "scrublet": config["scrublet"],
            "global_cell_id_template": config["global_cell_id_template"],
        }
        if method_contract is None:
            method_contract = current_contract
        elif current_contract != method_contract:
            raise ValueError(f"Shared QC contract drift in {study['scope']}")
        for record in manifest.to_dict("records"):
            library = str(record["library_id"])
            if library in seen_libraries:
                raise ValueError(f"Duplicate library across studies: {library}")
            seen_libraries.add(library)
            path = resolve(root, str(record["input_path"]))
            adapter = str(record["input_adapter"])
            if adapter == "cellranger_filtered_h5":
                if not path.is_file():
                    raise FileNotFoundError(path)
                with h5py.File(path, "r") as handle:
                    if "matrix/shape" not in handle:
                        raise ValueError(f"Missing matrix/shape: {path}")
                    shape = tuple(int(value) for value in handle["matrix/shape"][:])
                    if len(shape) != 2 or shape[0] <= 0:
                        raise ValueError(f"Invalid H5 matrix shape: {path}: {shape}")
                    observed = shape[1]
            elif adapter == "cellranger_filtered_mtx":
                observed = count_mtx_barcodes(path)
            else:
                raise ValueError(f"Unsupported production adapter: {adapter}")
            expected = int(record["expected_called_cells"])
            if observed != expected:
                raise ValueError(f"Called-cell mismatch for {library}: {observed} != {expected}")
            total_cells += observed
            rows.append({
                "scope": study["scope"], "study_id": config["study_id"],
                "library_id": library, "sample_accession": record["sample_accession"],
                "experimental_time": record["experimental_time"], "condition": record["condition"],
                "input_adapter": adapter, "observed_called_cells": observed,
                "expected_called_cells": expected, "status": "PASS",
            })

    if len(rows) != 193 or total_cells != 778207:
        raise ValueError(f"Production scope drift: {len(rows)} libraries, {total_cells} cells")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0])
    with args.report.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "status": "PASS", "study_scopes": len(studies), "libraries": len(rows),
        "called_cells": total_cells,
        "method_contract_sha256": sha256_text(json.dumps(method_contract, sort_keys=True)),
    }
    args.marker.write_text(json.dumps(summary, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
