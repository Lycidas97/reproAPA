"""Independent observation, null-tail and correction replay; no production imports."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE / "analysis"
LABELS = ["proximal", "distal", "nonrecurrent"]


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="stable")
    ranked = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty(len(values))
    result[order] = np.minimum(ranked, 1)
    return result


def main():
    read = lambda p: pd.read_csv(p, sep="\t", float_precision="round_trip")
    independent = read(HERE / "root_observed_labels.tsv")
    production = read(ANALYSIS / "measured_ortholog_members.tsv.gz")
    key = ["definition", "threshold", "human_gene_id", "mouse_gene_id"]
    cols = key + ["human_direction", "mouse_direction"]
    pd.testing.assert_frame_equal(independent[cols].sort_values(key).reset_index(drop=True),
                                  production[cols].sort_values(key).reset_index(drop=True), check_exact=True)
    results = read(ANALYSIS / "ortholog_3x3_threshold_scan.tsv")
    assert len(results) == 972 and int(results.tested.sum()) == 432
    raw_checks = 0
    with np.load(ANALYSIS / "null_draws.npz") as archive:
        assert len(archive.files) == 72
        for context, rows in results.groupby(["definition", "contrast", "threshold", "model"], sort=True):
            definition, contrast, threshold, model = context
            parts = ["ortholog", contrast, str(threshold), model]
            if contrast == "P_M":
                parts.append(definition)
            draws = archive["|".join(parts)]
            assert draws.shape == (100000, 3, 3)
            assert (draws.sum((1, 2)) == (728 if contrast == "P_M" else 1656)).all()
            for row in rows.itertuples():
                i, j = LABELS.index(row.human_direction), LABELS.index(row.mouse_direction)
                values = draws[:, i, j]
                tail = int((values >= row.observed).sum())
                p = (1 + tail) / 100001
                sd = float(values.std(ddof=1))
                assert tail == row.tail_count
                assert abs(float(values.mean()) - row.null_mean) < 1e-12
                assert abs(sd - row.null_sd) < 1e-12
                assert abs(float(values.mean()) - row.analytic_expected) <= 6 * sd / np.sqrt(100000) + 1e-9
                if row.tested:
                    assert abs(row.p_for_correction - (p if sd > 0 else 1)) < 1e-15
                if contrast == "P_M":
                    source = independent[(independent.definition == definition) & (independent.threshold == threshold)]
                    observed = int(((source.human_direction == row.human_direction) &
                                    (source.mouse_direction == row.mouse_direction)).sum())
                    assert observed == row.observed
                    assert abs(row.raw_p_enrichment - p) < 1e-15
                    if sd == 0:
                        assert row.inferential_status == "not_estimable_zero_null_variance" and pd.isna(row.p_enrichment)
                raw_checks += 1
    formal = results[results.tested].copy()
    family_counts = {}
    for name, fields, n in [("q8", ["definition", "threshold", "model"], 8),
                            ("q72", ["definition", "model"], 72),
                            ("q144", ["definition"], 144)]:
        count = 0
        for _, rows in formal.groupby(fields, sort=True):
            assert len(rows) == n
            np.testing.assert_allclose(bh(rows.p_for_correction), rows[name], rtol=0, atol=1e-15)
            count += 1
        family_counts[name] = count
    ps = formal[formal.contrast == "P_S"].drop_duplicates(["contrast", "threshold", "model", "human_direction", "mouse_direction"])
    unique = pd.concat([formal[formal.contrast == "P_M"], ps], ignore_index=True)
    assert len(unique) == 288
    np.testing.assert_allclose(bh(unique.p_for_correction), unique.q288, rtol=0, atol=1e-15)
    for _, rows in formal[formal.contrast == "P_S"].groupby(["threshold", "model", "human_direction", "mouse_direction"]):
        assert len(rows) == 3 and rows.p_for_correction.nunique() == rows.q288.nunique() == 1
    inputs = [HERE / "root_observed_labels.tsv", ANALYSIS / "measured_ortholog_members.tsv.gz",
              ANALYSIS / "ortholog_3x3_threshold_scan.tsv", ANALYSIS / "null_draws.npz"]
    result = {"status": "PASS", "independent_member_rows": len(independent), "raw_cell_checks": raw_checks,
              "formal_rows": len(formal), "correction_families_checked": family_counts,
              "unique_global_family": 288, "all_null_draw_totals_verified": True,
              "zero_variance_handling_verified": True,
              "input_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (HERE / "root_result_QA.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
