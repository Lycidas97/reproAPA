from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


ROOT = Path(__file__).parents[1]
PREPARE = ROOT / "scripts" / "prepare_scope.py"
EMBED = ROOT / "scripts" / "fit_embedding.py"
PREPARE_SENSITIVITY = ROOT / "scripts" / "prepare_sensitivity.py"


def ordered_sha256(values) -> str:
    return hashlib.sha256("".join(f"{value}\n" for value in values).encode("utf-8")).hexdigest()


def set_sha256(values) -> str:
    return ordered_sha256(sorted((str(value) for value in values), key=lambda value: value.encode("utf-8")))


def write_10x_h5(path: Path, counts: np.ndarray, barcodes: list[str]) -> None:
    matrix = sparse.csc_matrix(counts.T)
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=matrix.data)
        group.create_dataset("indices", data=matrix.indices)
        group.create_dataset("indptr", data=matrix.indptr)
        group.create_dataset("shape", data=np.asarray(matrix.shape, dtype=np.int64))
        group.create_dataset("barcodes", data=np.asarray(barcodes, dtype="S"))
        features = group.create_group("features")
        features.create_dataset("id", data=np.asarray([f"gene_{index}" for index in range(counts.shape[1])], dtype="S"))
        features.create_dataset("name", data=np.asarray([f"GENE{index}" for index in range(counts.shape[1])], dtype="S"))
        features.create_dataset("feature_type", data=np.asarray(["Gene Expression"] * counts.shape[1], dtype="S"))


class PrimarySensitivityPreparationTest(unittest.TestCase):
    def test_primary_fit_and_sensitivity_preparation_share_feature_and_pca_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = np.random.default_rng(11)
            counts = rng.poisson(lam=np.linspace(1.0, 6.0, 12), size=(120, 12)).astype(np.int32)
            barcodes = [f"BC{index:03d}-1" for index in range(120)]
            h5_path = root / "matrix.h5"
            write_10x_h5(h5_path, counts, barcodes)
            primary_ids = [f"lib1|{barcode}" for barcode in barcodes[:100]]
            sensitivity_ids = [f"lib1|{barcode}" for barcode in barcodes[100:]]
            primary_path = root / "primary.tsv.gz"
            sensitivity_path = root / "sensitivity.tsv.gz"
            pd.DataFrame({
                "study_id": ["TEST"] * 100,
                "library_id": ["lib1"] * 100,
                "barcode": barcodes[:100],
                "global_cell_id": primary_ids,
            }).to_csv(primary_path, sep="\t", index=False, compression="gzip")
            pd.DataFrame({
                "study_id": ["TEST"] * 20,
                "library_id": ["lib1"] * 20,
                "barcode": barcodes[100:],
                "global_cell_id": sensitivity_ids,
            }).to_csv(sensitivity_path, sep="\t", index=False, compression="gzip")
            metadata_path = root / "metadata.tsv"
            pd.DataFrame([{
                "library_id": "lib1", "species": "Homo sapiens", "sample_accession": "GSMTEST",
                "sample_title": "test", "experimental_time": "D0", "condition": "source",
            }]).to_csv(metadata_path, sep="\t", index=False)
            config = {
                "version": "3.0",
                "study_id": "GSE171127",
                "membership_path": str(primary_path),
                "sensitivity_membership_path": str(sensitivity_path),
                "library_metadata_path": str(metadata_path),
                "membership_columns": {
                    "cell_id": "global_cell_id", "state_id": None, "scope_id": None,
                    "library_id": "library_id", "barcode": "barcode",
                },
                "expression_mode": "multi_source",
                "expression_already_log1p": False,
                "partition_selection_fields": [],
                "scopes": {"studywide": {
                    "n_cells": 100,
                    "source_membership_order_sha256": ordered_sha256(primary_ids),
                    "cell_set_sha256": set_sha256(primary_ids),
                    "sensitivity_n_cells": 20,
                    "sensitivity_source_membership_order_sha256": ordered_sha256(sensitivity_ids),
                    "sensitivity_cell_set_sha256": set_sha256(sensitivity_ids),
                    "expression_sources": [{
                        "selector": "library_id", "value": "lib1", "mode": "10x_h5", "path": str(h5_path),
                    }],
                    "pipeline": {
                        "hvg_n_top_genes": 5, "pca_n_components_fitted": 3,
                        "pca_n_components_used": 3, "pca_seed": 11,
                        "neighbor_n_neighbors": 5, "leiden_resolution": 0.6,
                        "leiden_seed": 11, "filter_genes_min_cells": 1,
                        "normalization_method": "normalize_total", "normalization_target_sum": 10000.0,
                        "transform": "log1p", "hvg_flavor": "scanpy_seurat", "scale_max_value": 10,
                        "pca_solver": "sklearn_randomized_iterated_power_4",
                        "neighbor_backend": "pynndescent.NNDescent", "neighbor_metric": "euclidean",
                        "leiden_implementation": "leidenalg.find_partition", "leiden_iterations": -1,
                    },
                }},
            }
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            output_root = root / "output"
            subprocess.run([
                sys.executable, str(PREPARE), "--config", str(config_path), "--scope", "studywide",
                "--output-root", str(output_root),
            ], check=True)
            subprocess.run([
                sys.executable, str(EMBED), "--config", str(config_path), "--scope", "studywide",
                "--output-root", str(output_root), "--threads", "1",
            ], check=True)
            scope_root = output_root / "scopes" / "studywide"
            subprocess.run([
                sys.executable, str(PREPARE_SENSITIVITY), "--config", str(config_path),
                "--scope", "studywide", "--scope-root", str(scope_root),
            ], check=True)
            primary_features = pd.read_csv(scope_root / "prepared/features.tsv.gz", sep="\t")
            sensitivity_matrix = sparse.load_npz(scope_root / "sensitivity/prepared/log1p_gex.npz")
            sensitivity_cells = pd.read_csv(scope_root / "sensitivity/prepared/cells.tsv.gz", sep="\t")
            preparation_metadata = json.loads((scope_root / "prepared/metadata.json").read_text(encoding="utf-8"))
            pca_model_exists = (scope_root / "embedding/pca_model.npz").is_file()
        self.assertTrue(pca_model_exists)
        self.assertEqual(sensitivity_matrix.shape, (20, len(primary_features)))
        self.assertEqual(sensitivity_cells["cell_id"].tolist(), sensitivity_ids)
        self.assertIn("experimental_time", sensitivity_cells)
        self.assertEqual(preparation_metadata["execution_contract"], "source_membership_order")
        self.assertEqual(preparation_metadata["frozen_execution_object_sha256"], "NOT_APPLICABLE")


if __name__ == "__main__":
    unittest.main()
