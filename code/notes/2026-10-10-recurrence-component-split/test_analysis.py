"""Tests encode the sign interpretation and the conditional null's purpose."""
import unittest
from itertools import permutations

import numpy as np
import pandas as pd

from analyze import add_bh, bh, bins, dual_reference, inference, label_null


class ComponentIntent(unittest.TestCase):
    def test_alignment_distinguishes_distal_MS_from_proximal_PM(self):
        # Both patterns require positive raw A, but have opposite aligned allocations.
        distal_MS, distal_PM = .08, .01
        proximal_MS, proximal_PM = .01, -.08
        raw = np.array([distal_MS - distal_PM, proximal_MS - proximal_PM])
        aligned = np.array([1, -1]) * raw
        self.assertTrue((raw > 0).all())
        self.assertGreater(aligned[0], 0)
        self.assertLess(aligned[1], 0)
        self.assertAlmostEqual(aligned[0] - aligned[1], raw.sum())

    def test_labels_preserve_both_directions_and_stratum_expectations(self):
        values = np.array([-.1, .1, .3, -.4, .2, .6])
        labels = np.array([-1, 0, 1, -1, 0, 1])
        cells = [np.array([0, 1, 2]), np.array([3, 4, 5])]
        draws, expected = label_null(values, labels, cells, 10000, np.random.default_rng(17))
        exact = []
        for a in permutations(labels[:3]):
            for b in permutations(labels[3:]):
                p = np.array(a + b)
                self.assertEqual(np.sum(p == -1), 2)
                self.assertEqual(np.sum(p == 1), 2)
                exact.append([values[p == -1].mean(), values[p == 1].mean()])
        exact = np.array(exact)
        np.testing.assert_allclose([expected[-1], expected[1]], exact.mean(0), atol=1e-12)
        np.testing.assert_allclose(draws.mean(0), exact.mean(0), atol=.006)
        self.assertTrue(all(np.any(np.all(np.isclose(exact, row), axis=1)) for row in draws))

    def test_missing_route_cannot_create_dual_reference_gene(self):
        members = pd.DataFrame(dict(species=["mouse"] * 3, contrast=["P_M"] * 3,
            threshold=[.02] * 3, gene_id=["complete", "missing", "wrong_PS"],
            gene_symbol=["a", "b", "c"], direction=["proximal"] * 3))
        matched = pd.DataFrame([
            dict(species="mouse", route=route, gene_id=gene, DWUI_PS=ps, DWUI_PM=-.05)
            for gene, route, ps in [("complete", "CIP", -.03), ("complete", "WOT", -.04),
                                   ("missing", "CIP", -.08), ("wrong_PS", "CIP", -.05),
                                   ("wrong_PS", "WOT", .03)]])
        out = dual_reference(matched, members, .02)
        self.assertEqual(set(out.gene_id), {"complete"})

    def test_ties_stay_together_and_constant_null_has_no_p_value(self):
        values = np.array([0, 0, 0, 1, 1, 1])
        assigned, _ = bins(values, 5)
        self.assertEqual(len(set(assigned[:3])), 1)
        self.assertEqual(len(set(assigned[3:])), 1)
        values = np.array([-.1, .2, .3])
        labels = np.array([-1, 0, 1])
        draws, _ = label_null(values, labels, [np.array([0]), np.array([1]), np.array([2])],
                              20, np.random.default_rng(1))
        result = inference(-.1, draws[:, 0])
        self.assertFalse(result["informative"])
        self.assertTrue(np.isnan(result["p_upper"]))

    def test_BH_uses_entire_defined_family(self):
        np.testing.assert_allclose(bh([.01, .04, .03, .9]), [.04, .053333333333, .053333333333, .9])
        rows = pd.DataFrame(dict(null_model=["a"]*4, threshold=[.02]*4, p_upper=[.01,.04,.03,np.nan]))
        result = add_bh(rows)
        self.assertEqual(result.BH_family_size_per_threshold.iloc[0], 4)
        np.testing.assert_allclose(result.q_BH_per_threshold.iloc[:3], [.04,.053333333333,.053333333333])
        self.assertTrue(np.isnan(result.q_BH_per_threshold.iloc[3]))
        result = inference(.2, np.array([.2 - 1e-15, .1, .3]))
        self.assertEqual(result["tail_count"], 2)


if __name__ == "__main__":
    unittest.main()
