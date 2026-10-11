from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd


SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from program_audit import add_program_review_flags  # noqa: E402


class ProgramAuditTest(unittest.TestCase):
    def test_program_polarization_requires_review_without_dropping_clusters(self) -> None:
        clusters = pd.DataFrame({
            "cluster_id": ["C001", "C002", "C003"],
            "phase_fraction_G1_like": [0.40, 0.05, 0.40],
            "phase_fraction_S_like": [0.30, 0.45, 0.30],
            "phase_fraction_G2M_like": [0.30, 0.50, 0.30],
            "immediate_early_stress_mean": [0.10, 0.10, 0.60],
            "interferon_response_mean": [0.10, 0.10, 0.20],
        })

        audited = add_program_review_flags(clusters, 0.80, 0.50)

        self.assertEqual(audited["cluster_id"].tolist(), clusters["cluster_id"].tolist())
        self.assertEqual(audited["requires_program_review"].tolist(), [False, True, True])
        self.assertEqual(audited.loc[1, "program_review_reasons"], "cell_cycle_polarized")
        self.assertEqual(audited.loc[2, "program_review_reasons"], "stress_or_interferon_high")


if __name__ == "__main__":
    unittest.main()
