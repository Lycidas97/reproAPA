"""Intent checks: preserve biological contrast dependence and honest inference."""
import importlib.util
import itertools
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("analysis", Path(__file__).with_name("analyze.py"))
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class SpecificityIntent(unittest.TestCase):
    def test_joint_patterns_are_not_independently_reassigned(self):
        # Overlapping state responses must keep their dependence under the null.
        labels = np.array([[1, 0, 1, 0], [1, 0, 0, 1], [0, 1, 1, 0], [0, 0, 0, 0]], bool)
        strata = [np.arange(4)]
        out = a.joint_draw(labels, strata, 500, np.random.default_rng(8))
        reference = sorted(map(tuple, labels))
        for draw in out:
            self.assertEqual(sorted(map(tuple, draw)), reference)
            self.assertEqual(int((draw[:, 0] & draw[:, 2]).sum()), 1)

    def test_matched_expectation_agrees_with_complete_enumeration(self):
        # Independent species permutations do not mean independent human labels.
        h = np.array([[1, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 1], [0, 0, 1, 0]], bool)
        m = np.array([[1, 0], [0, 1], [1, 0], [0, 0]], bool)
        hs, ms = [np.array([0, 1]), np.array([2, 3])], [np.array([0, 2]), np.array([1, 3])]
        expectation = a.probabilities(h, hs).T @ a.probabilities(m, ms)
        results = []
        for hp in itertools.product(*[list(itertools.permutations(ii)) for ii in hs]):
            hd = h.copy()
            for ii, perm in zip(hs, hp):
                hd[ii] = h[list(perm)]
            for mp in itertools.product(*[list(itertools.permutations(ii)) for ii in ms]):
                md = m.copy()
                for ii, perm in zip(ms, mp):
                    md[ii] = m[list(perm)]
                results.append(hd.astype(int).T @ md.astype(int))
        np.testing.assert_allclose(np.mean(results, axis=0), expectation)
        null = a.simulate(h, m, hs, ms, 20000, np.random.default_rng(31))
        np.testing.assert_allclose(null.mean(0), expectation, atol=.025)

    def test_identical_control_cannot_establish_specificity(self):
        # An identical non-naive response must give zero specificity, even when overlap enriches.
        one = np.array([[1, 0], [0, 0], [0, 1], [1, 0]], bool)
        h = np.column_stack([one, one])
        m = np.array([[1, 0], [0, 1], [0, 0], [1, 0]], bool)
        ii = [np.arange(4)]
        exp = a.probabilities(h, ii).T @ a.probabilities(m, ii)
        null = a.simulate(h, m, ii, ii, 500, np.random.default_rng(1))
        diff = null[:, 0, 0]/exp[0, 0]-null[:, 2, 0]/exp[2, 0]
        np.testing.assert_array_equal(diff, 0)
        self.assertEqual(diff.std(), 0)

    def test_missing_ortholog_measurement_is_not_a_negative_label(self):
        # The comparison denominator is measured in every contrast, not the full ortholog map.
        orth = pd.DataFrame({"human_gene_id": ["h1", "h2", "h3"], "mouse_gene_id": ["m1", "m2", "m3"]})
        naive = pd.DataFrame({"human_gene_id": ["h1", "h2"], "naive": [True, False]})
        control = pd.DataFrame({"human_gene_id": ["h1", "h3"], "control": [False, True]})
        measured = orth.merge(naive, on="human_gene_id").merge(control, on="human_gene_id")
        self.assertEqual(measured.human_gene_id.tolist(), ["h1"])

    def test_fixed_thresholds_and_disjoint_signs(self):
        values = np.array([-.05, -.02, -.019, 0, .019, .02, .05])
        labels = a.direction_labels(values, .02)
        self.assertFalse((labels[:, 0] & labels[:, 1]).any())
        self.assertEqual(labels[:, 0].sum(), 2)
        self.assertEqual(labels[:, 1].sum(), 2)
        self.assertEqual(a.direction_labels(values, .05)[:, 0].sum(), 1)


if __name__ == "__main__":
    unittest.main()
