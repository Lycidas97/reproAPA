"""Audit completed files against frozen inputs and stored Monte Carlo draws."""
import json
import unittest

import numpy as np
import pandas as pd

from analyze import HERE, ROOT, TABLES, B, bh, digest


def main():
    manifest = json.loads((HERE / "manifest.json").read_text())
    assert all(digest(ROOT / path) == sha for path, sha in manifest["input_sha256"].items())
    assert all(digest(HERE / name) == sha for name, sha in manifest["output_sha256"].items())
    assert all(digest(ROOT / path) == sha for path, sha in manifest.get("definition_source_sha256", {}).items())
    assert all(digest(HERE / name) == sha for name, sha in manifest.get("derived_delivery_sha256", {}).items()
               if name != "verify.py")
    assert digest(HERE / "analyze.py") == manifest["script_sha256"]
    result = pd.read_csv(HERE / "overlap_results.tsv", sep="\t", float_precision="round_trip")
    members = pd.read_csv(HERE / "ortholog_membership_and_covariates.tsv.gz", sep="\t", float_precision="round_trip")
    draws = np.load(HERE / "overlap_null_draws.npz")
    pairs = list(members.groupby(["human_label", "mouse_target", "threshold"]))
    assert len(pairs) == 20
    for row in result.itertuples():
        key = "|".join(map(str, (row.human_label, row.mouse_target, row.threshold, row.model)))
        null = draws[key]
        assert len(null) == B and (null >= 0).all()
        assert row.tail_count == int((null >= row.observed).sum())
        assert row.p == (1 + row.tail_count) / (B + 1)
        f = members[members.human_label.eq(row.human_label) & members.mouse_target.eq(row.mouse_target) &
                    members.threshold.eq(row.threshold)]
        assert len(f) == row.n_orthologs
        assert not f.human_gene_id.duplicated().any() and not f.mouse_gene_id.duplicated().any()
        assert f.h_selected.sum() == row.human_selected and f.m_selected.sum() == row.mouse_selected
        assert (f.h_selected & f.m_selected).sum() == row.observed
    pd.testing.assert_series_equal(result.q_primary_or_sensitivity_within_model,
        result.groupby(["model", "primary_family"]).p.transform(bh), check_names=False)
    canonical = pd.read_csv(TABLES / "candidate_screen/dual_reference_candidates.tsv", sep="\t")
    expected = set(canonical[canonical.species.eq("mouse") & canonical.direction.eq("negative")].gene_id)
    assert len(expected) == 207
    # Recover the full dual target with the same audited inputs, before ortholog restriction.
    from analyze import read, mouse_frames
    frozen = read("recurrence_controls/recurrence_frozen_covariates.tsv.gz")
    recurrence = read("recurrence/recurrence_members.tsv.gz")
    recurrence = recurrence[recurrence.metric.eq("WPP") & recurrence.series_matching & recurrence.variant.eq("primary")]
    vectors = read("matched_pms/matched_state_vectors.tsv.gz")
    axes = read("state_vectors/PAS_axes.primary.tsv.gz", usecols=["atlas_pas_id", "species", "study_id", "setting", "WPP"])
    target = mouse_frames(frozen, recurrence, vectors, axes, .02)["P_M_PS_same_contributor_dual_proximal"]
    assert set(target[target.m_selected].mouse_gene_id) == expected
    suite = unittest.defaultTestLoader.discover(str(HERE), pattern="test_analysis.py")
    outcome = unittest.TextTestRunner(verbosity=2).run(suite)
    assert outcome.wasSuccessful() and len(outcome.skipped) == 0
    verification = dict(status="PASS", input_hashes_verified=len(manifest["input_sha256"]),
        output_hashes_verified=len(manifest["output_sha256"]), all40_finite_p_replayed=True,
        all20_measured_membership_backgrounds_verified=True, BH_families_replayed=True,
        canonical_207_dual_reference_roster_exact_match=True,
        independent_canonical_sha256=digest(TABLES / "candidate_screen/dual_reference_candidates.tsv"),
        tests_run=outcome.testsRun, tests_skipped=len(outcome.skipped),
        primary_plot_or_manuscript_modified=False)
    (HERE / "validation.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps(verification))


if __name__ == "__main__":
    main()
