"""Independent coverage, parent-contrast, conditional-null and hash checks."""
from pathlib import Path
import hashlib
import json
import math

import numpy as np
import pandas as pd

from verify_summary import bh


BASE = Path(__file__).resolve().parent
SITE = BASE / "sites"


def main():
    seal = json.loads((SITE / "event_sequence_sealed_manifest.json").read_text())
    assert seal["status"] == "SEALED_EVENT_SEQUENCE_FOR_INDEPENDENT_REVIEW"
    for name, expected in seal["files"].items():
        assert hashlib.sha256((SITE / name).read_bytes()).hexdigest() == expected, name
    assert hashlib.sha256((BASE / "SITE-STUDY-VOTE-CLARIFICATION.zh.md").read_bytes()).hexdigest() == seal["root_study_vote_clarification_sha256"]
    events = pd.read_csv(SITE / "site_event_contexts.tsv.gz", sep="\t")
    membership = pd.read_csv(BASE / "function/comparison_membership_features.tsv.gz", sep="\t")
    expected_keys = set()
    keys = ["source_population", "gene_threshold", "human_gene_id", "contrast", "tolerance"]
    for row in membership.to_dict("records"):
        if row["set_name"] in ["main48", "list45"]:
            contrasts = [c for c in ["P_S", "P_M"] if row[c + "_measurable"]]
        else:
            contrasts = ["hN_P__mouse_P_M"]
        for contrast in contrasts:
            for tolerance in [25, 50]:
                expected_keys.add((row["set_name"], row["threshold"], row["human_gene_id"], contrast, tolerance))
    assert set(events[keys].itertuples(index=False, name=None)) == expected_keys
    assert not events.duplicated(keys).any()
    for sp in ["human", "mouse"]:
        driver_ok = (events[sp + "_gain_id"] != events[sp + "_loss_id"]) & events[sp + "_gain_delta"].ge(.02) & events[sp + "_loss_delta"].le(-.02)
        assert np.array_equal(events[sp + "_driver_threshold"], driver_ok)
    assert np.array_equal(events.strict_UTR, events.human_strict_UTR & events.mouse_strict_UTR)
    assert np.array_equal(events.mass_eligible, events.order_eligible & events.human_min_retained_mass.ge(.8) & events.mouse_min_retained_mass.ge(.8))
    assert np.array_equal(events.strict_double_driver, events.strict_UTR & events.mass_eligible & events.double_driver_correspondence)
    assert np.array_equal(events.two_study_exact_event_both, events.strict_double_driver & events.human_gain_support_studies.ge(2) & events.mouse_gain_support_studies.ge(2))
    assert not (events.sequence_candidate & events.parent_case & events.context_role.ne("shared_recurrent")).any()
    assert not (events.sequence_candidate & ~events.parent_case & events.context_role.ne("single_species_recurrent")).any()

    group_columns = ["source_population", "gene_threshold", "set_name", "tolerance", "parent_case", "context_role"]
    stages = ["all_contexts", "strict_UTR", "order_eligible", "mass_eligible", "driver_judgeable", "double_driver_correspondence", "two_study_exact_event_both"]
    rebuilt = []
    for key, group in events.groupby(group_columns, sort=True):
        for stage in stages:
            if stage != "all_contexts":
                group = group[group[stage]]
            rebuilt.append(dict(zip(group_columns, key), stage=stage, contexts=len(group), unique_orthologs=group.human_gene_id.nunique()))
    attrition = pd.read_csv(SITE / "event_attrition.tsv", sep="\t")
    sort_columns = group_columns + ["stage"]
    pd.testing.assert_frame_equal(pd.DataFrame(rebuilt).sort_values(sort_columns).reset_index(drop=True), attrition.sort_values(sort_columns).reset_index(drop=True), check_dtype=False)

    gene_features = pd.read_csv(SITE / "sequence_gene_features_and_strata.tsv", sep="\t")
    results = pd.read_csv(SITE / "sequence_tests.tsv", sep="\t")
    draws = np.load(SITE / "sequence_null_draws.npz")
    assert len(results) == 8
    checked = []
    for row in results.itertuples(index=False):
        frame = gene_features[(gene_features.set_name == row.set_name) & (gene_features.tolerance == row.tolerance)]
        assert not frame.human_gene_id.duplicated().any()
        assert len(frame) == row.genes and int(frame.case.sum()) == row.cases
        score = frame.loc[frame.case, "score"].mean()
        expected = variance = 0.0
        mobile_cases = 0
        for _, layer in frame.groupby("sequence_" + row.model):
            n, k = len(layer), int(layer.case.sum())
            expected += k * layer.score.mean()
            if 0 < k < n:
                mobile_cases += k
            if n > 1:
                variance += k * (n - k) / (n * (n - 1)) * ((layer.score - layer.score.mean()) ** 2).sum()
        expected /= row.cases
        variance /= row.cases ** 2
        assert math.isclose(score, row.observed_mean, abs_tol=1e-12)
        assert math.isclose(expected, row.conditional_expected, abs_tol=1e-12)
        assert math.isclose(variance, row.null_sd ** 2, abs_tol=1e-12)
        assert mobile_cases == row.mobile_cases
        null = draws[f"{row.set_name}|{row.tolerance}|{row.model}"]
        assert len(null) == 20000
        tail = int(np.sum(null >= score - 1e-12))
        assert tail == row.upper_tail_count
        estimable = row.cases >= 3 and variance > 1e-24
        assert estimable == row.estimable
        expected_p = (tail + 1) / 20001 if estimable else 1.0
        assert math.isclose(expected_p, row.p, abs_tol=1e-12)
        checked.append({"set_name": row.set_name, "tolerance": row.tolerance, "model": row.model, "status": row.status})
    np.testing.assert_allclose(bh(results.p), results.q8, atol=1e-12)
    shared = events[events.source_population.isin(["main48", "hN79"]) & events.gene_threshold.eq(.02) & events.parent_case & events.context_role.eq("shared_recurrent")]
    shared.to_csv(BASE / "shared_event_summary.tsv", sep="\t", index=False)
    outcome = {
        "status": "PASS_SITE_SUMMARY_LOGIC_AND_HASHES",
        "sealed_files_verified": len(seal["files"]),
        "all_measurable_event_contexts_verified": len(events),
        "full_attrition_rows_rebuilt": len(attrition),
        "sequence_tests": checked,
        "scope": "All key/flag/attrition/null summaries. Source-vector and PAS-coordinate audit is separate in review/.",
    }
    (BASE / "root_sites_validation.json").write_text(json.dumps(outcome, indent=2) + "\n")
    print(json.dumps({k: v for k, v in outcome.items() if k not in ["sequence_tests", "scope"]}))


if __name__ == "__main__":
    main()
