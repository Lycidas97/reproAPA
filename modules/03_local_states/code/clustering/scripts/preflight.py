#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from cluster_utils import ordered_sha256, write_json, write_marker


FORBIDDEN_FIELDS = {
    "time", "stage", "marker", "program_score", "author_label", "route_outcome",
    "pas", "apa", "cross_study_state", "cross_study_mapping",
}


def set_sha256(values) -> str:
    return ordered_sha256(sorted((str(value) for value in values), key=lambda value: value.encode("utf-8")))


def read_tier(path: str, mapping: dict, expected_n: int, expected_hash: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    cell_column = mapping["cell_id"]
    if cell_column not in frame or frame[cell_column].duplicated().any() or frame[cell_column].isna().any():
        raise ValueError(f"invalid tier membership: {path}")
    if len(frame) != expected_n or set_sha256(frame[cell_column]) != expected_hash:
        raise ValueError(f"tier membership contract drift: {path}")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--marker", required=True)
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    if config["study_id"] == "GSE225142" or args.scope not in config["scopes"]:
        raise ValueError("scope is outside the six-study contract")
    if FORBIDDEN_FIELDS & set(config.get("partition_selection_fields", [])):
        raise ValueError("biological labels/programs cannot be used to select the partition")
    scope = config["scopes"][args.scope]
    required = [
        Path(config["membership_path"]), Path(config["sensitivity_membership_path"]),
        Path(config["technical_membership_path"]), Path(config["library_metadata_path"]),
        Path(config["qc_contract"]["independent_validation_path"]),
    ]
    if scope.get("expression_sources"):
        required.extend(Path(source["path"]) for source in scope["expression_sources"])
    else:
        required.append(Path(scope["expression_path"]))
    for optional in ["expression_cells_path", "expression_features_path"]:
        if config.get(optional):
            required.append(Path(config[optional]))
    if scope.get("frozen_execution_object"):
        required.append(Path(scope["frozen_execution_object"]))
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing inputs: {missing}")
    mapping = config["membership_columns"]
    tiers = config["qc_contract"]["tiers"]
    membership = read_tier(
        config["membership_path"], mapping, int(scope["n_cells"]), scope["cell_set_sha256"]
    )
    sensitivity = read_tier(
        config["sensitivity_membership_path"], mapping, int(scope["sensitivity_n_cells"]),
        scope["sensitivity_cell_set_sha256"],
    )
    technical = read_tier(
        config["technical_membership_path"], mapping, int(tiers["technical_exclude"]["n_cells"]),
        tiers["technical_exclude"]["cell_set_sha256"],
    )
    all_ids = pd.concat([membership[mapping["cell_id"]], sensitivity[mapping["cell_id"]], technical[mapping["cell_id"]]])
    if all_ids.duplicated().any():
        raise ValueError("QC tiers overlap")
    validation = pd.read_csv(config["qc_contract"]["independent_validation_path"], sep="\t", dtype=str)
    if "status" not in validation or not validation["status"].eq("PASS").all():
        raise ValueError("independent QC validation is not all PASS")
    library_metadata = pd.read_csv(config["library_metadata_path"], sep="\t", dtype=str)
    if library_metadata["library_id"].duplicated().any():
        raise ValueError("library metadata is not one-to-one")
    membership_libraries = set(pd.concat([membership, sensitivity, technical])[mapping["library_id"]])
    if membership_libraries != set(library_metadata["library_id"]):
        raise ValueError("QC memberships do not close against library metadata")
    if scope.get("expression_sources"):
        source_libraries = {str(source["value"]) for source in scope["expression_sources"]}
        if source_libraries != set(library_metadata["library_id"]):
            raise ValueError("expression sources do not close against library metadata")
    observed = len(membership)
    report = {"study_id": config["study_id"], "scope_id": args.scope, "n_cells": observed,
              "n_sensitivity_cells": len(sensitivity), "n_technical_cells": len(technical),
              "n_libraries": len(library_metadata),
              "qc_method_contract_sha256": config["qc_contract"]["method_contract_sha256"],
              "biological_fields_used_for_clustering": False,
              "expression_mode": config["expression_mode"], "inputs": [str(path) for path in required], "status": "PASS"}
    write_json(Path(args.report), report)
    write_marker(Path(args.marker), status="PASS", study_id=config["study_id"], scope_id=args.scope)


if __name__ == "__main__":
    main()
