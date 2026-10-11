from __future__ import annotations

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


PRODUCTION = Path(__file__).parents[1] / "production"
sys.path.insert(0, str(PRODUCTION))

from check_acceptance import membership_ids, set_hash  # noqa: E402


class ProductionAcceptanceTest(unittest.TestCase):
    def write_membership(self, path: Path, values: list[str]) -> None:
        with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["global_cell_id"], delimiter="\t")
            writer.writeheader()
            writer.writerows({"global_cell_id": value} for value in values)

    def test_membership_hash_is_order_independent(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "membership.tsv.gz"
            self.write_membership(path, ["cell_b", "cell_a"])
            values = membership_ids(path)
        self.assertEqual(values, {"cell_a", "cell_b"})
        self.assertEqual(set_hash(values), set_hash({"cell_b", "cell_a"}))

    def test_duplicate_membership_ids_fail_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "membership.tsv.gz"
            self.write_membership(path, ["cell_a", "cell_a"])
            with self.assertRaisesRegex(ValueError, "duplicate global_cell_id"):
                membership_ids(path)


if __name__ == "__main__":
    unittest.main()
