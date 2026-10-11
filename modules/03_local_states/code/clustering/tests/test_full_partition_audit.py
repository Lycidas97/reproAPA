from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "scripts" / "full_partition_audit.py"


def run_audit(tmp_path: Path, frozen_labels: list[str]) -> tuple[dict, pd.DataFrame]:
    cell_ids = [f"cell_{index:03d}" for index in range(100)]
    evaluation_labels = ["0"] * 50 + ["1"] * 50
    evaluation = tmp_path / "full_partitions.tsv.gz"
    selected = tmp_path / "selected_resolution.json"
    frozen = tmp_path / "cell_clusters.tsv.gz"
    output = tmp_path / "audit"
    pd.DataFrame({"cell_id": cell_ids, "res_0.6": evaluation_labels}).to_csv(
        evaluation, sep="\t", index=False, compression="gzip"
    )
    selected.write_text(json.dumps({"selected_resolution": 0.6}), encoding="utf-8")
    pd.DataFrame({"cell_id": cell_ids, "cluster_id": frozen_labels}).to_csv(
        frozen, sep="\t", index=False, compression="gzip"
    )
    subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--evaluation-partitions", str(evaluation),
            "--selected-json", str(selected),
            "--frozen-clusters", str(frozen),
            "--minimum-cluster-cells", "10",
            "--minimum-evaluation-cluster-cells", "5",
            "--minimum-transfer-ari", "0.70",
            "--minimum-cluster-jaccard", "0.50",
            "--minimum-cluster-purity", "0.70",
            "--minimum-cluster-recall", "0.70",
            "--maximum-cluster-count-drift-fraction", "0.25",
            "--output-dir", str(output),
        ],
        check=True,
    )
    return (
        json.loads((output / "summary.json").read_text(encoding="utf-8")),
        pd.read_csv(output / "cluster_transfer_audit.tsv", sep="\t"),
    )


class FullPartitionAuditTest(unittest.TestCase):
    def test_exact_evaluation_to_full_transfer_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            summary, clusters = run_audit(Path(directory), ["C001"] * 50 + ["C002"] * 50)
        self.assertEqual(summary["status"], "PASS")
        self.assertEqual(summary["adjusted_rand_index"], 1.0)
        self.assertEqual(summary["unresolved_cluster_count"], 0)
        self.assertEqual(clusters["full_n_cells"].sum(), 100)

    def test_full_data_split_is_retained_but_requires_review(self) -> None:
        labels = ["C001"] * 25 + ["C003"] * 25 + ["C002"] * 50
        with tempfile.TemporaryDirectory() as directory:
            summary, clusters = run_audit(Path(directory), labels)
        self.assertEqual(summary["status"], "REVIEW")
        self.assertIs(summary["all_full_cells_retained"], True)
        self.assertGreaterEqual(summary["unresolved_cluster_count"], 1)
        self.assertEqual(clusters["full_n_cells"].sum(), 100)
        self.assertEqual(
            set(clusters["cluster_release_class"]), {"stable_candidate", "rare_or_unresolved"}
        )


if __name__ == "__main__":
    unittest.main()
