"""Verify every gradient keeps its original P and whole-M source profiles."""
import json

import numpy as np

from verify_APA_delivery import HERE, read


def main():
    keys = ["study_id", "route", "gene_id", "library_key_a", "library_key_b"]
    baseline = read("baseline_pair_effects.tsv.gz", usecols=keys + ["target_usage", "reference_usage", "PAS_ids", "delta_DWUI"]).set_index(keys)
    gradient = read("gradient_pair_effects.tsv.gz")
    maximum = 0.
    for row in gradient.to_dict("records"):
        key = (row["study_id"], row["route"], row["gene_id"], row["library_key"], row["library_key"])
        original = baseline.loc[key]
        assert json.loads(row["PAS_ids"]) == json.loads(original.PAS_ids)
        for field, original_field in [("P_usage", "target_usage"), ("M_whole_usage", "reference_usage")]:
            error = float(np.max(np.abs(np.asarray(json.loads(row[field])) - np.asarray(json.loads(original[original_field])))))
            assert error < 1e-12
            maximum = max(maximum, error)
        assert abs(row["delta_whole"] - original.delta_DWUI) < 1e-12
        assert row["P_group"] == "|".join([row["study_id"], row["route"], row["library_key"], "supported50", row["gene_id"], "original_P"])
    report = dict(status="PASS", gradient_pairs=len(gradient), fixed_original_P_and_whole_M_max_vector_error=maximum)
    (HERE / "gradient_reference_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
