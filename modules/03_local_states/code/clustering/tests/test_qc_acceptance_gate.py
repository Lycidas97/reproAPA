from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_qc_clustering_inputs import DATASET_KEYS, validate_qc_acceptance  # noqa: E402


class QCAcceptanceGateTest(unittest.TestCase):
    def test_qc_scope_contract_uses_wot_registry_name(self) -> None:
        self.assertEqual(
            set(DATASET_KEYS),
            {"WOT", "GSE103221", "GSE118258", "GSE150311", "GSE171127", "GSE178325"},
        )
        self.assertEqual(DATASET_KEYS["WOT"], "wot")

    def test_non_pass_acceptance_is_rejected_before_clustering_config_generation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / "runs/_production_acceptance"
            target.mkdir(parents=True)
            report = {
                "status": "VISUAL_PENDING",
                "n_studies": 6,
                "expected_libraries": 193,
                "expected_called_cells": 778207,
                "studies": [{"scope": scope, "status": "VISUAL_PENDING"} for scope in DATASET_KEYS],
            }
            (target / "acceptance.json").write_text(json.dumps(report), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not the frozen six-study PASS contract"):
                validate_qc_acceptance(root)


if __name__ == "__main__":
    unittest.main()
