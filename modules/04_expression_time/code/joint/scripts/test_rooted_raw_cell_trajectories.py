"""Tests guard why a continuous display path is a supported temporal path."""
import unittest
import pandas as pd

from build_rooted_raw_cell_trajectories import best_rooted_paths, timed_node


def table(rows):
    return pd.DataFrame(rows, columns=['source', 'target', 'day_from', 'day_to', 'mass', 'target_fraction'])


class RootedTrajectoryTests(unittest.TestCase):
    def test_same_state_continuation_restores_the_missing_time_interval(self):
        a, b, c = [timed_node('state_R', day) for day in (0, 1, 2)]
        edges = table([(a, b, 0, 1, .8, .9), (b, c, 1, 2, .8, .9)])
        chosen, paths = best_rooted_paths(edges, a, [c], .005, .1)
        self.assertEqual(len(chosen), 2)
        self.assertEqual(paths[0]['edge_indices'], [0, 1])

    def test_a_gap_must_not_be_reconnected_by_state_name(self):
        edges = table([('S0', 'R1', 0, 1, .9, .9), ('R2', 'E3', 2, 3, .9, .9)])
        with self.assertRaisesRegex(ValueError, 'cannot be reached'):
            best_rooted_paths(edges, 'S0', ['E3'], .005, .1)

    def test_a_strong_incoming_edge_from_an_unrooted_source_is_not_a_new_root(self):
        edges = table([('S0', 'R1', 0, 1, .7, .7), ('R1', 'E2', 1, 2, .2, .3),
                       ('other1', 'E2', 1, 2, .8, .7)])
        chosen, _ = best_rooted_paths(edges, 'S0', ['E2'], .005, .1)
        self.assertEqual(list(chosen.index), [0, 1])

    def test_supported_branches_share_their_trunk_without_duplicate_roots(self):
        edges = table([('S0', 'R1', 0, 1, .9, .9), ('R1', 'E2', 1, 2, .6, .8),
                       ('R1', 'M2', 1, 2, .3, .7)])
        chosen, paths = best_rooted_paths(edges, 'S0', ['E2', 'M2'], .005, .1)
        self.assertEqual(len(chosen), 3)
        self.assertEqual(set(chosen.source)-set(chosen.target), {'S0'})
        self.assertEqual(len(paths), 2)

    def test_path_choice_uses_full_path_support_not_just_last_edge(self):
        edges = table([('S0', 'A1', 0, 1, .8, .8), ('S0', 'B1', 0, 1, .2, .2),
                       ('A1', 'E2', 1, 2, .4, .4), ('B1', 'E2', 1, 2, .6, .6)])
        chosen, _ = best_rooted_paths(edges, 'S0', ['E2'], .005, .1)
        self.assertEqual(list(chosen.index), [0, 2])

    def test_weak_or_nonforward_connections_cannot_be_used_as_bridges(self):
        with self.assertRaisesRegex(ValueError, 'cannot be reached'):
            best_rooted_paths(table([('S0', 'E1', 0, 1, .001, .9)]), 'S0', ['E1'], .005, .1)
        with self.assertRaisesRegex(ValueError, 'strictly forward'):
            best_rooted_paths(table([('S1', 'E0', 1, 0, .8, .9)]), 'S1', ['E0'], .005, .1)


if __name__ == '__main__':
    unittest.main()
