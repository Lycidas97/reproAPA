#!/usr/bin/env python3
"""Build frozen study-wide clustering configs from an accepted QC production run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml


DATASET_KEYS = {
    "WOT": "wot",
    "GSE103221": "gse103221",
    "GSE118258": "gse118258",
    "GSE150311": "gse150311",
    "GSE171127": "gse171127",
    "GSE178325": "gse178325",
}

PIPELINE = {
    "hvg_n_top_genes": 3000,
    "pca_n_components_fitted": 50,
    "pca_n_components_used": 30,
    "pca_seed": 20260825,
    "neighbor_n_neighbors": 30,
    "leiden_resolution": 0.6,
    "leiden_seed": 20260825,
    "filter_genes_min_cells": 10,
    "normalization_method": "normalize_total",
    "normalization_target_sum": 10000.0,
    "transform": "log1p",
    "hvg_flavor": "scanpy_seurat",
    "scale_max_value": 10,
    "pca_solver": "sklearn_randomized_iterated_power_4",
    "neighbor_backend": "pynndescent.NNDescent",
    "neighbor_metric": "euclidean",
    "leiden_implementation": "leidenalg.find_partition",
    "leiden_iterations": -1,
}


def ordered_sha256(values) -> str:
    payload = "".join(f"{value}\n" for value in values)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def set_sha256(values) -> str:
    return ordered_sha256(sorted((str(value) for value in values), key=lambda value: value.encode("utf-8")))


def read_membership(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t", dtype=str)
    required = {"study_id", "library_id", "barcode", "global_cell_id"}
    if missing := required - set(frame):
        raise ValueError(f"{path}: missing membership columns {sorted(missing)}")
    if frame["global_cell_id"].duplicated().any() or frame["global_cell_id"].isna().any():
        raise ValueError(f"{path}: invalid global_cell_id")
    return frame


def tier_contract(run_dir: Path) -> tuple[dict, dict[str, pd.DataFrame]]:
    validation_path = run_dir / "independent_validation.tsv"
    if not validation_path.is_file():
        raise FileNotFoundError(validation_path)
    validation = pd.read_csv(validation_path, sep="\t", dtype=str)
    if "status" not in validation or not validation["status"].eq("PASS").all():
        raise ValueError(f"{run_dir}: independent QC validation is not all PASS")
    counts = pd.read_csv(run_dir / "tier_counts_and_hashes.tsv", sep="\t", dtype=str).set_index("tier")
    filenames = {
        "primary_pass": "primary_membership.tsv.gz",
        "sensitivity_only": "sensitivity_membership.tsv.gz",
        "technical_exclude": "technical_membership.tsv.gz",
    }
    memberships = {tier: read_membership(run_dir / filename) for tier, filename in filenames.items()}
    contract = {}
    for tier, frame in memberships.items():
        expected_n = int(counts.loc[tier, "n_cells"])
        expected_hash = str(counts.loc[tier, "membership_sha256"])
        observed_hash = set_sha256(frame["global_cell_id"])
        if len(frame) != expected_n or observed_hash != expected_hash:
            raise ValueError(f"{run_dir}: {tier} membership does not match the accepted QC ledger")
        contract[tier] = {"n_cells": len(frame), "cell_set_sha256": observed_hash}
    return contract, memberships


def load_baseline_hash(path: Path) -> str:
    if not path.is_file():
        return "NOT_AVAILABLE"
    baseline = yaml.safe_load(path.read_text(encoding="utf-8"))
    return str(baseline["scopes"]["studywide"]["cell_set_sha256"])


def validate_qc_acceptance(qc_root: Path) -> tuple[Path, dict]:
    acceptance_path = qc_root / "runs/_production_acceptance/acceptance.json"
    if not acceptance_path.is_file():
        raise FileNotFoundError(acceptance_path)
    report = json.loads(acceptance_path.read_text(encoding="utf-8"))
    expected = {"n_studies": 6, "expected_libraries": 193, "expected_called_cells": 778207}
    if report.get("status") != "PASS" or any(report.get(key) != value for key, value in expected.items()):
        raise ValueError("QC production acceptance is not the frozen six-study PASS contract")
    statuses = {row.get("scope"): row.get("status") for row in report.get("studies", [])}
    if statuses != {scope: "PASS" for scope in DATASET_KEYS}:
        raise ValueError("QC production acceptance does not contain six study-level PASS records")
    review = pd.read_csv(qc_root / "production/manual_visual_review.tsv", sep="\t", dtype=str).fillna("")
    if set(review["scope"]) != set(DATASET_KEYS) or not review["status"].eq("PASS").all():
        raise ValueError("QC manual visual review ledger is not all PASS")
    registry = pd.read_csv(qc_root / "production/studies.tsv", sep="\t", dtype=str)
    sources = [qc_root / "production/manual_visual_review.tsv", qc_root / "production/studies.tsv"]
    for item in registry.itertuples(index=False):
        run = qc_root / "runs" / item.output_dir
        sources.extend([
            run / "independent_validation.tsv",
            run / "tier_counts_and_hashes.tsv",
            run / "figures/qc_contact_sheets/contact_sheet_index.tsv",
        ])
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"QC acceptance inputs are missing: {missing}")
    if acceptance_path.stat().st_mtime_ns < max(path.stat().st_mtime_ns for path in sources):
        raise ValueError("QC acceptance report is stale relative to a validation or visual-review input")
    return acceptance_path, report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-root", required=True)
    parser.add_argument("--clustering-root", required=True)
    parser.add_argument("--baseline-config-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    qc_root = Path(args.qc_root).resolve()
    clustering_root = Path(args.clustering_root).resolve()
    baseline_dir = Path(args.baseline_config_dir).resolve()
    target = Path(args.output_dir).resolve()
    target.mkdir(parents=True, exist_ok=True)
    configs_dir = target / "configs"
    configs_dir.mkdir(exist_ok=True)

    acceptance_path, acceptance = validate_qc_acceptance(qc_root)
    preflight = json.loads((qc_root / "runs/_production_preflight/preflight.ok").read_text(encoding="utf-8"))
    if preflight.get("status") != "PASS" or int(preflight.get("study_scopes", 0)) != 6:
        raise ValueError("QC production preflight is not the accepted six-study contract")
    registry = pd.read_csv(qc_root / "production/studies.tsv", sep="\t", dtype=str)
    if set(registry["scope"]) != set(DATASET_KEYS):
        raise ValueError("QC registry study scope drift")

    project_configs = {}
    rows = []
    for item in registry.itertuples(index=False):
        dataset = DATASET_KEYS[item.scope]
        run_dir = qc_root / "runs" / item.output_dir
        contract, memberships = tier_contract(run_dir)
        manifest_path = qc_root / item.config_path.replace("configs/", "manifests/").replace(".yaml", ".tsv")
        manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
        if len(manifest) != int(item.expected_libraries):
            raise ValueError(f"{item.scope}: QC manifest library count drift")
        if set(manifest["library_id"]) != set(pd.concat(memberships.values())["library_id"]):
            raise ValueError(f"{item.scope}: tier memberships do not close against the QC manifest")
        if not manifest["input_adapter"].eq("cellranger_filtered_h5").all():
            raise ValueError(f"{item.scope}: unsupported expression input adapter")
        missing_inputs = [path for path in manifest["input_path"] if not Path(path).is_file()]
        if missing_inputs:
            raise FileNotFoundError(f"{item.scope}: missing expression inputs: {missing_inputs[:3]}")

        primary = memberships["primary_pass"]
        sensitivity = memberships["sensitivity_only"]
        config = {
            "version": "3.0",
            "study_id": "WOT_GSE115943" if item.scope == "WOT" else item.scope,
            "species": sorted(manifest["species"].unique().tolist()),
            "membership_path": str(run_dir / "primary_membership.tsv.gz"),
            "sensitivity_membership_path": str(run_dir / "sensitivity_membership.tsv.gz"),
            "technical_membership_path": str(run_dir / "technical_membership.tsv.gz"),
            "library_metadata_path": str(manifest_path),
            "membership_columns": {
                "cell_id": "global_cell_id",
                "state_id": None,
                "scope_id": None,
                "library_id": "library_id",
                "barcode": "barcode",
            },
            "expression_mode": "multi_source",
            "expression_already_log1p": False,
            "partition_selection_fields": [],
            "library_use": "stratified_resampling_and_postfreeze_audit_only",
            "qc_contract": {
                "method_contract_sha256": preflight["method_contract_sha256"],
                "production_acceptance_path": str(acceptance_path),
                "independent_validation_path": str(run_dir / "independent_validation.tsv"),
                "tiers": contract,
            },
            "scopes": {
                "studywide": {
                    "role": "primary",
                    "n_cells": len(primary),
                    "source_membership_order_sha256": ordered_sha256(primary["global_cell_id"]),
                    "cell_set_sha256": set_sha256(primary["global_cell_id"]),
                    "sensitivity_n_cells": len(sensitivity),
                    "sensitivity_source_membership_order_sha256": ordered_sha256(sensitivity["global_cell_id"]),
                    "sensitivity_cell_set_sha256": set_sha256(sensitivity["global_cell_id"]),
                    "expression_sources": [
                        {
                            "selector": "library_id",
                            "value": row.library_id,
                            "mode": "10x_h5",
                            "path": row.input_path,
                        }
                        for row in manifest.sort_values("library_id", kind="mergesort").itertuples(index=False)
                    ],
                    "pipeline": dict(PIPELINE),
                }
            },
        }
        config_path = configs_dir / f"{dataset}.yaml"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        project_configs[dataset] = str(config_path.relative_to(clustering_root))
        baseline_hash = load_baseline_hash(baseline_dir / f"{dataset}.yaml")
        rows.append({
            "dataset": dataset,
            "study_scope": item.scope,
            "libraries": len(manifest),
            "primary_cells": len(primary),
            "sensitivity_cells": len(sensitivity),
            "technical_cells": len(memberships["technical_exclude"]),
            "primary_cell_set_sha256": contract["primary_pass"]["cell_set_sha256"],
            "baseline_primary_cell_set_sha256": baseline_hash,
            "baseline_exact_match": baseline_hash == contract["primary_pass"]["cell_set_sha256"],
            "clustering_action": "reuse_only_after_new_release_gates" if baseline_hash == contract["primary_pass"]["cell_set_sha256"] else "rerun",
            "config_path": str(config_path),
        })

    project = {
        "python": "/root/miniconda3/envs/single-cell-qc/bin/python",
        "project_config_path": str(target / "config.production.yaml"),
        "output_root": "results_qc_v3_2_studywide_v1",
        "dataset_configs": project_configs,
        "selection": {
            "min_evaluation_cells": 10000,
            "target_evaluation_fraction": 0.20,
            "max_evaluation_cells": 50000,
            "subsample_fraction": 0.8,
            "n_resamples": 100,
            "n_shards": 20,
            "seed": 20260825,
            "ci_resamples": 25000,
            "resolutions": [0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0],
        },
        "diagnostics": {
            "minimum_cluster_cells": 50,
            "minimum_evaluation_cluster_cells": 20,
            "minimum_transfer_ari": 0.70,
            "minimum_cluster_transfer_jaccard": 0.50,
            "minimum_cluster_transfer_purity": 0.70,
            "minimum_cluster_transfer_recall": 0.70,
            "maximum_cluster_count_drift_fraction": 0.25,
            "classifier_max_cells": 50000,
            "marker_top_n": 30,
            "umap_max_points_per_panel": 120000,
            "maximum_cell_cycle_polarization_without_review": 0.80,
            "maximum_positive_program_mean_without_review": 0.50,
            "projection_neighbors": 30,
            "projection_support_quantile": 0.99,
        },
        "resources": {
            "preflight_mem_mb": 4000,
            "prepare_mem_mb": 64000,
            "embedding_mem_mb": 128000,
            "sweep_mem_mb": 64000,
            "sweep_shard_mem_mb": 8000,
            "sweep_merge_mem_mb": 64000,
            "summarize_mem_mb": 32000,
            "diagnostic_mem_mb": 128000,
        },
    }
    project_path = target / "config.production.yaml"
    project_path.write_text(yaml.safe_dump(project, sort_keys=False), encoding="utf-8")
    pd.DataFrame(rows).sort_values("dataset", kind="mergesort").to_csv(
        target / "qc_to_clustering_registry.tsv", sep="\t", index=False
    )
    (target / "build.ok").write_text(
        f"status=PASS\nstudies={len(rows)}\nmethod_contract_sha256={preflight['method_contract_sha256']}\n"
        f"qc_acceptance_status={acceptance['status']}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
