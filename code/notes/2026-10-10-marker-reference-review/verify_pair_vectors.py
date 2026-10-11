"""Check every delivered pair vector and its original/star cohort identity."""
import csv
import gzip
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "notes/2026-10-10-marker-reference-sensitivity"
HERE = Path(__file__).resolve().parent
CONTRASTS = {"D00": ("original_P", "original_M"), "D01": ("original_P", "Mstar"),
             "D10": ("Pstar", "original_M"), "D11": ("Pstar", "Mstar")}


def main():
    total, maximum = 0, 0.
    with gzip.open(SOURCE / "pair_effects.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            position = np.asarray(json.loads(row["position"]))
            a = np.asarray(json.loads(row["target_usage"]))
            b = np.asarray(json.loads(row["reference_usage"]))
            assert len(position) == len(a) == len(b) == len(json.loads(row["PAS_ids"]))
            assert abs(a.sum() - 1) < 1e-12 and abs(b.sum() - 1) < 1e-12
            error = abs(float((a - b) @ position) - float(row["delta_DWUI"]))
            assert error < 1e-12
            maximum = max(maximum, error); total += 1
            for side, cohort in zip(["target", "reference"], CONTRASTS[row["contrast"]]):
                expected = row["version"] + "_" + cohort if cohort in ["Pstar", "Mstar"] else cohort
                assert row[side + "_cohort"] == expected
                assert row[side + "_group"] == "|".join([row["study_id"], row["route"], row["library_key"], row["setting"], row["gene_id"], expected])
    report = dict(status="PASS", pair_rows=total, max_pair_DWUI_vector_abs_difference=maximum,
                  all_usage_mass_and_contrast_group_IDs_checked=True)
    (HERE / "frozen_pair_vectors_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
