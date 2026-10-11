"""Tests for inferential intent: measured cohorts, state specificity, conditional null."""
import itertools
import json
import unittest

import numpy as np
import pandas as pd

from analyze import binary_draw, bh, human_frames, membership_probabilities, mouse_frames, simulate


class AnalysisIntentTests(unittest.TestCase):
    def test_missing_naive_condition_is_not_a_nonresponse(self):
        axes = pd.DataFrame(dict(atlas_pas_id=["a", "b"], study_id=["GSE150311"] * 2,
                                 setting=["supported50_p08"] * 2, WPP=[.2, 1.]))
        records = []
        for gene, conditions in [("complete", ["6", "8"]), ("only_one_measured", ["6"])]:
            for condition in conditions:
                records.append(dict(version="v2_primary", setting="supported50_p08", gene_id=gene,
                    gene_symbol=gene, PAS_number=2, tandem_utr_compatible=True,
                    library_a="lib" + condition, node_a="N", library_b="lib7", node_b="P",
                    n_mc_a=10, n_mc_b=10, pas_counts_a=300, pas_counts_b=400,
                    PAS_ids=json.dumps(["a", "b"]), usage_a=json.dumps([.8, .2]),
                    usage_b=json.dumps([.5, .5]), delta_DWUI=-.24))
        frames, _ = human_frames(pd.DataFrame(records), axes)
        for name, frame in frames.items():
            if name != "RSeT_same_library_N_vs_P":
                self.assertEqual(set(frame.human_gene_id), {"complete"})
                # The supported axis starts at .2, so capacity is .6-.2=.4.
                self.assertAlmostEqual(frame.h_capacity.iloc[0], .4)

    def test_PM_recurrence_alone_does_not_establish_source_relative_specificity(self):
        frozen, recurrence, vectors, axes = [], [], [], []
        for number, gene in enumerate(["fails_source_reference", "passes_both_references"]):
            for unit, route in [("GSE103221", "CIP"), ("GSE115943", "WOT")]:
                frozen.append(dict(species="mouse", contrast="P_M", threshold=.05, unit=unit,
                    gene_id=gene, gene_symbol=gene, event=-1, mean_GEX=30., APA_coverage=30.,
                    PAS_number=2, baseline=.6, capacity_negative=.6))
                pas = [unit + str(number) + "a", unit + str(number) + "b"]
                axes.extend(dict(species="mouse", setting="supported50_p08", atlas_pas_id=p, WPP=w)
                            for p, w in zip(pas, [0., 1.]))
                vectors.append(dict(species="mouse", setting="supported50_p08", route=route,
                    gene_id=gene, PAS_ids=json.dumps(pas), M_DWUI=.6, S_DWUI=.6,
                    DWUI_PM=-.1, DWUI_PS=.1 if number == 0 and route == "CIP" else -.1))
            recurrence.append(dict(species="mouse", contrast="P_M", threshold=.05,
                                   direction="negative", gene_id=gene))
        frames = mouse_frames(pd.DataFrame(frozen), pd.DataFrame(recurrence),
                              pd.DataFrame(vectors), pd.DataFrame(axes), .05)
        self.assertEqual(frames["P_M_proximal_recurrent"].m_selected.sum(), 2)
        dual = frames["P_M_PS_same_contributor_dual_proximal"].set_index("mouse_gene_id")
        self.assertFalse(dual.loc["fails_source_reference", "m_selected"])
        self.assertTrue(dual.loc["passes_both_references", "m_selected"])

    def test_independent_species_null_matches_exact_conditional_enumeration(self):
        h = np.array([True, False, True, False])
        m = np.array([True, True, False, False])
        hs = [np.array([0, 1]), np.array([2, 3])]
        ms = [np.array([0, 2]), np.array([1, 3])]
        def possibilities(strata):
            result = []
            for first, second in itertools.product(*[list(ii) for ii in strata]):
                label = np.zeros(4, bool)
                label[[first, second]] = True
                result.append(label)
            return result
        exact = np.array([(a & b).sum() for a in possibilities(hs) for b in possibilities(ms)])
        analytic = membership_probabilities(h, hs) @ membership_probabilities(m, ms)
        self.assertAlmostEqual(analytic, exact.mean())
        null = simulate(h, m, hs, ms, 20000, np.random.default_rng(47))
        self.assertLess(abs(null.mean() - exact.mean()), .03)
        self.assertLess(abs((null >= 2).mean() - (exact >= 2).mean()), .02)
        for labels, strata in [(h, hs), (m, ms)]:
            draws = binary_draw(labels, strata, 100, np.random.default_rng(19))
            for ii in strata:
                self.assertTrue(np.all(draws[:, ii].sum(1) == labels[ii].sum()))

    def test_BH_keeps_primary_family_correction_instead_of_uncorrected_significance(self):
        self.assertTrue(np.allclose(bh([.03, .9]), [.06, .9]))
        self.assertTrue(np.allclose(bh([.001, .02, .5]), [.003, .03, .5]))


if __name__ == "__main__":
    unittest.main()
