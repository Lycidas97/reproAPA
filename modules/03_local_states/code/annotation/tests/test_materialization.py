from __future__ import annotations

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "release"


def rows(name: str) -> list[dict[str, str]]:
    with (RELEASE / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class MaterializationTest(unittest.TestCase):
    def test_validation_passes(self) -> None:
        validation = rows("validation_summary.tsv")
        self.assertTrue(validation, "Validation contract must contain explicit checks")
        self.assertEqual({row["status"] for row in validation}, {"PASS"})

    def test_cluster_and_cell_closure(self) -> None:
        master = rows("cluster_annotation_master.tsv")
        self.assertEqual(len(master), 94, "Frozen cluster universe must not change")
        self.assertEqual(sum(int(row["n_primary_cells"]) for row in master), 606_186)
        self.assertEqual(sum(int(row["n_sensitivity_cells"]) for row in master), 76_435)

    def test_wot_c001_uses_r10b(self) -> None:
        master = rows("cluster_annotation_master.tsv")
        row = next(
            row
            for row in master
            if row["study_id"] == "WOT_GSE115943" and row["local_cluster_id"] == "C001"
        )
        expected = {
            "authority_release": "R10B",
            "local_annotation_label": "early_source_context",
            "local_working_state": "source_like",
            "transition_mode": "stable_source",
            "route_relation": "not_applicable",
            "hard_identity": "NA",
            "manual_hold": "true",
            "public_state_eligible": "false",
            "cross_study_vote_allowed": "false",
        }
        self.assertEqual({key: row[key] for key in expected}, expected)

    def test_only_strict_anchors_are_initially_eligible(self) -> None:
        adjudication = rows("anchor_adjudication.tsv")
        observed = {
            (row["study_id"], row["local_cluster_id"])
            for row in adjudication
            if row["alignment_anchor_eligible"] == "true"
        }
        expected = {
            ("GSE171127", "C006"),
            ("WOT_GSE115943", "C007"),
            ("WOT_GSE115943", "C009"),
            ("WOT_GSE115943", "C010"),
            ("WOT_GSE115943", "C011"),
        }
        self.assertEqual(observed, expected)
        for row in adjudication:
            if row["alignment_anchor_eligible"] == "true":
                self.assertEqual(row["manual_hold"], "false")
                self.assertEqual(row["library_context_triggered"], "false")

    def test_reviewed_anchor_gate(self) -> None:
        adjudication = rows("anchor_adjudication_final.tsv")
        promoted = [
            row for row in adjudication if row["alignment_training_anchor_eligible"] == "true"
        ]
        self.assertEqual(
            len(promoted),
            15,
            "Neural-like remains a study-specific identity; four existing and eleven promoted anchors train",
        )
        for row in promoted:
            self.assertEqual(row["manual_hold"], "false")
            self.assertEqual(row["library_context_triggered"], "false")

    def test_species_pilot_gate(self) -> None:
        readiness = {row["species"]: row for row in rows("species_pilot_readiness.tsv")}
        self.assertEqual(readiness["Mus_musculus"]["pilot_status"], "READY")
        self.assertEqual(readiness["Homo_sapiens"]["pilot_status"], "BLOCKED")


if __name__ == "__main__":
    unittest.main()
