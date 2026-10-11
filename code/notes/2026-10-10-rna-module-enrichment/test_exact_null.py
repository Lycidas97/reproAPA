"""Tests encode why matching preserves size and can render inference impossible."""
import itertools
import math
import unittest
import numpy as np
from exact_null import component, distribution, upper_tail, rejection_probability, detectable_effect, bh


class ExactNullTests(unittest.TestCase):
    def test_fixed_target_stratum_size_not_independent_bernoulli(self):
        N, K, n = 7, 3, 2
        counts = [sum(i < K for i in subset) for subset in itertools.combinations(range(N), n)]
        lo, p = component(N, K, n)
        observed = np.bincount(counts, minlength=len(p)) / math.comb(N, n)
        np.testing.assert_allclose(p, observed, atol=1e-15)

    def test_joint_stratum_convolution_preserves_every_margin(self):
        params = [(5, 2, 2), (4, 1, 1)]
        samples = []
        for left in itertools.combinations(range(5), 2):
            for right in itertools.combinations(range(4), 1):
                samples.append(sum(i < 2 for i in left) + sum(i < 1 for i in right))
        lo, p = distribution(params)
        np.testing.assert_allclose(p, np.bincount(samples, minlength=len(p)) / len(samples), atol=1e-15)
        self.assertEqual(lo, 0)

    def test_singletons_do_not_supply_inferential_power(self):
        params = [(1, 1, 1), (1, 0, 1), (1, 1, 0)]
        lo, p = distribution(params)
        self.assertEqual(lo, 1)
        np.testing.assert_array_equal(p, [1.0])
        self.assertEqual(upper_tail((lo, p), 1), 1.0)
        self.assertFalse(detectable_effect(params, .05)['attainable'])

    def test_zero_annotation_category_remains_in_prespecified_family(self):
        self.assertEqual(upper_tail(distribution([(100, 0, 14)]), 0), 1.0)
        np.testing.assert_allclose(bh([.01, 1.0]), [.02, 1.0])

    def test_power_is_exact_under_null_and_increases_under_annotation_tilt(self):
        params = [(80, 20, 15)]
        null = distribution(params)
        powers = [rejection_probability(null, distribution(params, odds), .05)
                  for odds in [1, 2, 4, 8]]
        self.assertLessEqual(powers[0], .05)
        self.assertTrue(np.all(np.diff(powers) > 0))
        mde = detectable_effect(params, .05)
        self.assertTrue(mde['attainable'])
        self.assertGreaterEqual(rejection_probability(null, distribution(params, mde['minimum_odds']), .05), .8 - 1e-12)

    def test_disjoint_rare_terms_can_be_unattainable_after_family_correction(self):
        self.assertFalse(detectable_effect([(100, 1, 14)], .05 / 14)['attainable'])


if __name__ == '__main__':
    unittest.main()
