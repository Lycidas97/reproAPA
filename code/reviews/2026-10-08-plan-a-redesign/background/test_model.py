"""Small exhaustive checks for the scientific invariants of A1 randomization."""
import importlib.util
import itertools
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

spec = importlib.util.spec_from_file_location('conditional_reference', Path(__file__).with_name('analyze.py'))
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class ConditionalReferenceTests(unittest.TestCase):
    def test_bounded_transform_is_identity_without_clipping_including_boundaries(self):
        m = np.array([0., 0., 1., 1., .2, .2, .8, .8])
        s = np.array([0., .7, 1., .3, 0., 1., .2, .9])
        z = analysis.bounded_residual(m, s)
        np.testing.assert_allclose(analysis.inverse_residual(m, z), s, atol=1e-15)
        self.assertTrue(np.all((z >= -1) & (z <= 1)))
        # Only identical boundary status can share donors; no latent direction collapses.
        for indices in [[0, 1], [2, 3], [4, 5, 6, 7]]:
            for order in itertools.permutations(indices):
                recipient = analysis.inverse_residual(m[indices], z[list(order)])
                np.testing.assert_allclose(analysis.bounded_residual(m[indices], recipient), z[list(order)], atol=1e-15)

    def test_assignment_moments_match_every_permutation_not_independent_bernoulli(self):
        matrices = [np.array([[1]]), np.array([[1, 0], [0, 1]]),
                    np.array([[1, 1, 0], [0, 0, 0], [0, 1, 1]]),
                    np.array([[1, 0, 0, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])]
        for hits in matrices:
            values = np.array([hits[np.arange(len(hits)), order].sum() for order in itertools.permutations(range(len(hits)))])
            mean, variance = analysis.assignment_moments(hits)
            self.assertAlmostEqual(mean, values.mean(), places=14)
            self.assertAlmostEqual(variance, values.var(), places=14)
            i, j = linear_sum_assignment(hits)
            self.assertEqual(int(hits[i, j].sum()), int(values.min()))
            i, j = linear_sum_assignment(hits, maximize=True)
            self.assertEqual(int(hits[i, j].sum()), int(values.max()))

    def test_fixed_shared_P_PM_selection_and_joint_route_donor(self):
        p = np.array([[.25, .30], [.55, .50], [.35, .75], [.70, .65]])
        m = np.array([[.35, .40], [.45, .40], [.45, .65], [.60, .55]])
        s = np.array([[.40, .45], [.35, .30], [.20, .80], [.80, .70]])
        direction = np.array([-1, 1, -1, 1])
        eligible = np.array([True, True, False, False])
        z = analysis.bounded_residual(m, s)
        f = pd.DataFrame({'direction': direction, 'pm_eligible': eligible,
                          'dual': eligible & ((p - s) * direction[:, None] >= analysis.TAU).all(axis=1)})
        for index, route in enumerate(analysis.ROUTES['mouse']):
            for name, vector in [('P_DWUI', p), ('M_DWUI', m), ('S_DWUI', s), ('bounded_S_minus_M', z)]:
                f[route + '__' + name] = vector[:, index]
        for model in ['bounded_MS', 'S_vector']:
            hits = analysis.hit_matrix(f, 'mouse', model)
            self.assertFalse(hits[~eligible].any())
            for order in itertools.permutations(range(len(f))):
                reference = analysis.inverse_residual(m, z[list(order)]) if model == 'bounded_MS' else s[list(order)]
                rebuilt_PS = p - reference
                np.testing.assert_allclose(rebuilt_PS, (p - m) + (m - reference), atol=1e-15)
                direct = eligible & (rebuilt_PS * direction[:, None] >= analysis.TAU).all(axis=1)
                np.testing.assert_array_equal(hits[np.arange(len(f)), order], direct)
        # A PM-ineligible gene cannot gain membership, even when its PS can align.
        self.assertFalse(analysis.hit_matrix(f, 'mouse', 'S_vector')[2:].any())


if __name__ == '__main__':
    unittest.main()
