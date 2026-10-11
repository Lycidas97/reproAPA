"""Display smoothing must not change temporal topology or invent displacement."""
import unittest
import numpy as np
from matplotlib.path import Path

from smooth_trajectory_display import smooth_chain


class TrajectorySmoothingTests(unittest.TestCase):
    def test_small_zigzags_are_removed_without_moving_the_source_or_target(self):
        xy = np.array([[0., 0.], [1., .15], [2., -.2], [3., .1], [4., 0.]])
        original = xy.copy()
        path, audit = smooth_chain(xy)
        np.testing.assert_array_equal(path.vertices[[0, -1]], xy[[0, -1]])
        np.testing.assert_array_equal(xy, original)
        self.assertEqual(audit['retained_waypoints'], 2)

    def test_major_turns_are_curved_and_remain_near_the_original_path(self):
        path, audit = smooth_chain([[0., 0.], [10., 0.], [10., 10.]])
        self.assertIn(Path.CURVE3, path.codes)
        self.assertLessEqual(audit['deviation_upper_bound_mm'], 1.5)
        self.assertTrue(((path.vertices >= 0) & (path.vertices <= 10)).all())

    def test_shared_fork_is_exactly_the_same_point_in_all_incident_chains(self):
        trunk, _ = smooth_chain([[0., 0.], [2., .2], [4., 1.]])
        a, _ = smooth_chain([[4., 1.], [5., 2.5], [7., 3.]])
        b, _ = smooth_chain([[4., 1.], [5., -.2], [7., -2.]])
        np.testing.assert_array_equal(trunk.vertices[-1], a.vertices[0])
        np.testing.assert_array_equal(trunk.vertices[-1], b.vertices[0])

    def test_single_edge_stays_straight_and_keeps_forward_endpoint_order(self):
        path, _ = smooth_chain([[7., 3.], [-1., 5.]])
        np.testing.assert_array_equal(path.vertices, [[7., 3.], [-1., 5.]])
        self.assertEqual(list(path.codes), [Path.MOVETO, Path.LINETO])

    def test_large_corner_rounding_is_reduced_to_respect_the_bound(self):
        _, audit = smooth_chain([[0., 0.], [30., 0.], [30., 30.]], cut_mm=20.)
        self.assertLess(audit['rounding_cut_mm'], 20.)
        self.assertLessEqual(audit['deviation_upper_bound_mm'], 1.5)

    def test_invalid_geometry_fails_instead_of_silently_dropping_a_route(self):
        for points in ([[0., 0.], [np.nan, 1.]], [[0., 0.], [0., 0.]]):
            with self.assertRaises(ValueError):
                smooth_chain(points)


if __name__ == '__main__':
    unittest.main()
