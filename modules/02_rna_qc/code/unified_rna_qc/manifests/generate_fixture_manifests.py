#!/usr/bin/env python3
"""Materialize the two explicitly named task fixtures from retained local evidence."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT.parent / "reproapa-unified-rna-qc-v3-builder-pilot-gse178325"


def main() -> None:
    inventory = pd.read_csv(PRIOR / "input_inventory.tsv", sep="\t")
    samples = pd.read_csv(PRIOR / "sample_qc_decisions.tsv", sep="\t")
    joined = samples.merge(inventory[["library_id", "filtered_h5", "n_called_cells"]], on="library_id", suffixes=("", "_inventory"), validate="one_to_one")
    rows = []
    for row in joined.itertuples(index=False):
        prior = f"../reproapa-unified-rna-qc-v3-builder-pilot-gse178325/per_library/{row.library_id}"
        rows.append({
            "study_id": "GSE178325", "species": "Homo sapiens", "dataset_id": row.dataset_id,
            "library_id": row.library_id, "sample_id": row.sample_id,
            "sample_accession": row.sample_accession, "experimental_time": row.experimental_time,
            "condition": row.condition, "author_state": row.author_state,
            "input_adapter": "cellranger_filtered_h5", "input_path": row.filtered_h5,
            "expected_called_cells": int(row.n_called_cells),
            "basic_qc_evidence_adapter": "validated_prior_run", "prior_basic_qc_path": prior,
            "diagnostic_evidence_adapter": "validated_prior_run", "prior_diagnostic_path": prior,
        })
    pd.DataFrame(rows).to_csv(ROOT / "manifests/GSE178325_regression.tsv", sep="\t", index=False)

    matrix_root = Path("/root/nfsdata/apadetect/data/reprogramming_differentiation_apa/analyses/raw_realignment_apa_2026-06-07/matrices/gse103221_pluripotency_reprogramming_core")
    policy = "contract_inputs/historical_qc/report.md"
    diagnostic = "contract_inputs/historical_qc/production_policy_recommendation.tsv"
    small = [
        ("gse103221_pluripotency_reprogramming_core_GSM3629847", "GSM3629847", 4678, "D0", "starting_state"),
        ("gse103221_pluripotency_reprogramming_core_GSM3629849", "GSM3629849", 9072, "D8", "reprogramming"),
    ]
    rows = []
    for library, accession, cells, timepoint, condition in small:
        rows.append({
            "study_id": "GSE103221", "species": "Mus musculus", "dataset_id": "gse103221_two_library_real_pilot",
            "library_id": library, "sample_id": accession, "sample_accession": accession,
            "experimental_time": timepoint, "condition": condition, "author_state": "NOT_USED",
            "input_adapter": "cellranger_filtered_mtx", "input_path": str(matrix_root / library / "filtered_feature_bc_matrix"),
            "expected_called_cells": cells,
            "basic_qc_evidence_adapter": "policy_document_recompute", "prior_basic_qc_path": policy,
            "diagnostic_evidence_adapter": "scrublet_policy_recompute", "prior_diagnostic_path": diagnostic,
        })
    pd.DataFrame(rows).to_csv(ROOT / "manifests/GSE103221_small.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
