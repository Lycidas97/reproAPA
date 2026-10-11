from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from select_cohort import choose_evaluation_size  # noqa: E402


class EvaluationCohortPolicyTest(unittest.TestCase):
    def test_large_studies_receive_twenty_percent_up_to_fifty_thousand(self) -> None:
        self.assertEqual(choose_evaluation_size(240829, 10000, 0.20, 50000), 48166)
        self.assertEqual(choose_evaluation_size(145072, 10000, 0.20, 50000), 29015)

    def test_smaller_studies_preserve_the_ten_thousand_floor_without_using_all_cells(self) -> None:
        self.assertEqual(choose_evaluation_size(29000, 10000, 0.20, 50000), 10000)
        self.assertEqual(choose_evaluation_size(12484, 10000, 0.20, 50000), 10000)

    def test_study_smaller_than_floor_uses_every_available_cell(self) -> None:
        self.assertEqual(choose_evaluation_size(8000, 10000, 0.20, 50000), 8000)

    def test_invalid_policy_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            choose_evaluation_size(100, 10000, 0.0, 50000)


if __name__ == "__main__":
    unittest.main()
