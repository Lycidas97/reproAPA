#!/usr/bin/env python3

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from common import read_manifest, sha256_file
from path_contract import managed_run_dir
from process_library import evidence_guard


class PathContractTests(unittest.TestCase):
    def test_managed_fresh_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = managed_run_dir(root, "runs", "study", "fresh")
            self.assertEqual(run, root / "runs" / "study")
            run.mkdir(parents=True)
            (run / "partial.txt").write_text("partial", encoding="utf-8")
            with self.assertRaises(ValueError):
                managed_run_dir(root, "runs", "study", "fresh")
            self.assertEqual(managed_run_dir(root, "runs", "study", "resume"), run)

    def test_escape_and_symlink_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as outside:
            root = Path(tmp)
            for output_root, output_dir in (("/tmp/runs", "study"), ("runs", "../study")):
                with self.assertRaises(ValueError):
                    managed_run_dir(root, output_root, output_dir, "fresh")
            (root / "runs").symlink_to(Path(outside), target_is_directory=True)
            with self.assertRaises(ValueError):
                managed_run_dir(root, "runs", "study", "fresh")


class MetadataContractTests(unittest.TestCase):
    def fixture(self, root: Path) -> dict:
        registry = root / "registry.csv"
        registry.write_text(
            "series_accession,sample_accession,sample_title,organism,source_url\n"
            "GSE1,GSM1,10x_osk_d1,Mus musculus,https://example/GSM1\n",
            encoding="utf-8",
        )
        contract = root / "contract.tsv"
        contract.write_text(
            "sample_accession\texperimental_time\tcondition\tmetadata_derivation\n"
            "GSM1\tD1\tOSK_reprogramming\tGEO title identifies day 1\n",
            encoding="utf-8",
        )
        fields = [
            "study_id", "species", "dataset_id", "library_id", "sample_id",
            "sample_accession", "sample_title", "experimental_time", "condition",
            "author_state", "metadata_derivation", "source_url", "input_adapter",
            "input_path", "expected_called_cells", "basic_qc_evidence_adapter",
            "prior_basic_qc_path", "diagnostic_evidence_adapter", "prior_diagnostic_path",
        ]
        row = dict.fromkeys(fields, "unused")
        row.update({
            "study_id": "GSE1", "species": "Mus musculus", "dataset_id": "d",
            "library_id": "lib1", "sample_id": "GSM1", "sample_accession": "GSM1",
            "sample_title": "10x_osk_d1", "experimental_time": "D1",
            "condition": "OSK_reprogramming", "metadata_derivation": "GEO title identifies day 1",
            "source_url": "https://example/GSM1", "input_adapter": "cellranger_filtered_h5",
            "expected_called_cells": "10", "basic_qc_evidence_adapter": "validated_prior_run",
            "diagnostic_evidence_adapter": "validated_prior_run",
        })
        manifest = root / "manifest.tsv"
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        return {
            "library_manifest": "manifest.tsv", "metadata_registry": "registry.csv",
            "metadata_contract": "contract.tsv", "expected_library_count": 1,
            "expected_called_cell_count": 10, "study_id": "GSE1", "species": "Mus musculus",
            "library_id_pattern": "^lib[0-9]+$",
            "adapters": {"raw_input": ["cellranger_filtered_h5"],
                         "basic_qc_evidence": ["validated_prior_run"],
                         "diagnostic_evidence": ["validated_prior_run"]},
        }

    def test_wrong_stage_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self.fixture(root)
            self.assertEqual(read_manifest(root, config).iloc[0]["experimental_time"], "D1")
            manifest = (root / "manifest.tsv").read_text(encoding="utf-8").replace("\tD1\t", "\tD8\t")
            (root / "manifest.tsv").write_text(manifest, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "experimental_time"):
                read_manifest(root, config)


class PolicyEvidenceTests(unittest.TestCase):
    def test_generic_policy_can_govern_a_new_study(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            policy = root / "policy.txt"
            policy.write_text("ReproAPA RNA QC v1\nScrublet\n", encoding="utf-8")
            config = {
                "study_id": "GSE_NEW",
                "evidence_sha256": {"policy.txt": sha256_file(policy)},
                "evidence_required_text": ["ReproAPA RNA QC v1", "Scrublet"],
            }
            observed, digest = evidence_guard(root, "policy.txt", config)
            self.assertEqual(observed, policy)
            self.assertEqual(digest, sha256_file(policy))

    def test_policy_file_drift_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            policy = root / "policy.txt"
            policy.write_text("ReproAPA RNA QC v1\nScrublet\n", encoding="utf-8")
            config = {
                "evidence_sha256": {"policy.txt": "0" * 64},
                "evidence_required_text": ["ReproAPA RNA QC v1", "Scrublet"],
            }
            with self.assertRaisesRegex(ValueError, "SHA256 drift"):
                evidence_guard(root, "policy.txt", config)


if __name__ == "__main__":
    unittest.main()
