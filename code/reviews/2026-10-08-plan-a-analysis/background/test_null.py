"""Signed randomization must not count an opposite-direction event as support."""
import itertools
import unittest
import numpy as np
from analyze import signed_pmf


class SignedNullTests(unittest.TestCase):
    def test_matches_exhaustive_signed_assignments(self):
        for npos, nneg, nzero, kpos, kneg in [(2, 2, 1, 1, 2), (1, 2, 2, 2, 1), (0, 2, 1, 0, 2)]:
            population = np.array([1] * npos + [-1] * nneg + [0] * nzero)
            counts = np.zeros(kpos + kneg + 1)
            for positive in itertools.combinations(range(len(population)), kpos):
                remaining = set(range(len(population))) - set(positive)
                for negative in itertools.combinations(remaining, kneg):
                    support = int((population[list(positive)] == 1).sum() + (population[list(negative)] == -1).sum())
                    counts[support] += 1
            np.testing.assert_allclose(signed_pmf(npos, nneg, nzero, kpos, kneg), counts / counts.sum(), atol=1e-12)

    def test_no_candidates_has_no_support(self):
        np.testing.assert_array_equal(signed_pmf(3, 2, 1, 0, 0), [1.])

    def test_all_positive_candidates_cannot_match_negative_events(self):
        np.testing.assert_array_equal(signed_pmf(1, 2, 0, 3, 0), [0., 1., 0., 0.])


if __name__ == '__main__':
    unittest.main()
