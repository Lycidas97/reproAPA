from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


SCRIPT = Path(__file__).parents[1] / "scripts" / "project_sensitivity.py"
SCRIPTS = SCRIPT.parent
sys.path.insert(0, str(SCRIPTS))

from project_sensitivity import cluster_color_map  # noqa: E402


class SensitivityProjectionTest(unittest.TestCase):
    def test_cluster_palette_does_not_repeat_above_tab20_limit(self) -> None:
        levels = [f"C{index:03d}" for index in range(1, 28)]
        colors = cluster_color_map(levels)
        self.assertEqual(len(colors), 27)
        self.assertEqual(len(set(colors.values())), 27)

    def test_projection_retains_cells_and_recovers_primary_neighborhoods(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = np.random.default_rng(7)
            primary_pca = np.vstack([
                rng.normal(-2.0, 0.12, size=(60, 5)),
                rng.normal(2.0, 0.12, size=(60, 5)),
            ]).astype(np.float32)
            np.save(root / "primary_pca.npy", primary_pca)
            primary = pd.DataFrame({
                "cell_id": [f"primary_{index:03d}" for index in range(120)],
                "cluster_id": ["C001"] * 60 + ["C002"] * 60,
                "library_id": ["lib_a"] * 60 + ["lib_b"] * 60,
            })
            primary.to_csv(root / "primary_clusters.tsv.gz", sep="\t", index=False, compression="gzip")
            primary.assign(umap_1=primary_pca[:, 0], umap_2=primary_pca[:, 1]).to_csv(
                root / "primary_umap.tsv.gz", sep="\t", index=False, compression="gzip"
            )
            np.savez_compressed(
                root / "pca_model.npz",
                selected_feature_indices=np.arange(5, dtype=np.int64),
                feature_mean=np.full(5, 3.0, dtype=np.float32),
                feature_scale=np.ones(5, dtype=np.float32),
                components=np.eye(5, dtype=np.float32),
                scale_max_value=np.asarray(10.0, dtype=np.float32),
            )
            sensitivity_matrix = sparse.csr_matrix(np.vstack([
                np.full((10, 5), 1.0, dtype=np.float32),
                np.full((10, 5), 5.0, dtype=np.float32),
            ]))
            sparse.save_npz(root / "sensitivity.npz", sensitivity_matrix)
            pd.DataFrame({
                "study_id": ["TEST"] * 20,
                "scope_id": ["studywide"] * 20,
                "cell_id": [f"sensitivity_{index:03d}" for index in range(20)],
                "library_id": ["lib_a"] * 10 + ["lib_b"] * 10,
            }).to_csv(root / "sensitivity_cells.tsv.gz", sep="\t", index=False, compression="gzip")
            output = root / "output"
            subprocess.run([
                sys.executable, str(SCRIPT),
                "--sensitivity-matrix", str(root / "sensitivity.npz"),
                "--sensitivity-cells", str(root / "sensitivity_cells.tsv.gz"),
                "--pca-model", str(root / "pca_model.npz"),
                "--primary-pca", str(root / "primary_pca.npy"),
                "--primary-clusters", str(root / "primary_clusters.tsv.gz"),
                "--primary-umap", str(root / "primary_umap.tsv.gz"),
                "--neighbors", "10",
                "--support-quantile", "0.99",
                "--pca-components-used", "5",
                "--seed", "7",
                "--threads", "1",
                "--max-plot-cells", "1000",
                "--output-dir", str(output),
            ], check=True)
            projection = pd.read_csv(output / "sensitivity_projection.tsv.gz", sep="\t")
        self.assertEqual(len(projection), 20)
        self.assertEqual(projection.iloc[:10]["assigned_cluster"].unique().tolist(), ["C001"])
        self.assertEqual(projection.iloc[10:]["assigned_cluster"].unique().tolist(), ["C002"])
        self.assertTrue((projection["maximum_cluster_probability"] > 0.99).all())


if __name__ == "__main__":
    unittest.main()
