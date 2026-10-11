#!/usr/bin/env python3
"""Materialize deterministic tier memberships from an independently accepted QC ledger."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


FILENAMES = {
    "primary_pass": "primary_membership.tsv.gz",
    "sensitivity_only": "sensitivity_membership.tsv.gz",
    "technical_exclude": "technical_membership.tsv.gz",
}
MEMBERSHIP_COLUMNS = ["study_id", "library_id", "barcode", "global_cell_id"]


def set_hash(values) -> str:
    payload = "".join(f"{value}\n" for value in sorted((str(value) for value in values)))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--decisions", required=True)
    parser.add_argument("--tier-counts", required=True)
    parser.add_argument("--independent-validation", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    validation = pd.read_csv(args.independent_validation, sep="\t", dtype=str)
    if "status" not in validation or not validation["status"].eq("PASS").all():
        raise ValueError("independent QC validation is not all PASS")
    decisions = pd.read_parquet(args.decisions)
    missing = set([*MEMBERSHIP_COLUMNS, "final_qc_tier"]) - set(decisions)
    if missing:
        raise ValueError(f"QC decision ledger is missing columns: {sorted(missing)}")
    if decisions["global_cell_id"].duplicated().any() or decisions["global_cell_id"].isna().any():
        raise ValueError("invalid global_cell_id in QC decision ledger")
    expected = pd.read_csv(args.tier_counts, sep="\t", dtype=str).set_index("tier")
    if set(expected.index) != set(FILENAMES):
        raise ValueError("tier count contract drift")
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    rows = []
    for tier, filename in FILENAMES.items():
        membership = decisions.loc[
            decisions["final_qc_tier"].eq(tier), MEMBERSHIP_COLUMNS
        ].sort_values("global_cell_id", kind="mergesort")
        observed_hash = set_hash(membership["global_cell_id"])
        if len(membership) != int(expected.loc[tier, "n_cells"]):
            raise ValueError(f"{tier}: membership count mismatch")
        if observed_hash != str(expected.loc[tier, "membership_sha256"]):
            raise ValueError(f"{tier}: membership hash mismatch")
        membership.to_csv(
            target / filename, sep="\t", index=False,
            compression={"method": "gzip", "mtime": 0},
        )
        rows.append({"tier": tier, "n_cells": len(membership), "membership_sha256": observed_hash})
    pd.DataFrame(rows).to_csv(target / "materialization.tsv", sep="\t", index=False)
    (target / "materialization.ok").write_text("status=PASS\n", encoding="utf-8")


if __name__ == "__main__":
    main()
