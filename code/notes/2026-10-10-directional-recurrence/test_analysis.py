"""Checks for the scientific meaning of directional recurrence and its null."""
import itertools
import unittest

import numpy as np
import pandas as pd

from analyze import COLS, SIM_SOURCE, bh, numerical_function, observed_members, prepare_units, simulate


def frame_from_votes(votes):
    rows = []
    for gene, labels in votes.items():
        for study, label in enumerate(labels):
            if label is None:
                continue
            row = dict(gene_id=gene, gene_symbol=gene, unit=f"s{study}", event=label)
            row.update({c: 0 for c in COLS})
            row.update(can_negative=True, can_positive=True, eligible_units=sum(x is not None for x in labels))
            rows.append(row)
    return pd.DataFrame(rows)


class DirectionalRecurrenceTest(unittest.TestCase):
    def test_two_studies_determine_direction_despite_opposing_third(self):
        frame = frame_from_votes({"prox": [-1, -1, 1], "dist": [-1, 1, 1],
                                  "missing": [-1, None, 0], "two": [1, 1, None],
                                  "zero": [0, 0, 0]})
        result = observed_members(frame).set_index("gene_id")
        self.assertEqual(result.direction.to_dict(), {"dist": "distal", "prox": "proximal", "two": "distal"})
        self.assertEqual(result.loc["prox", "distal_support"], 1)
        self.assertEqual(result.loc["two", "eligible_studies"], 2)

    def test_four_studies_require_a_different_double_recurrence_rule(self):
        with self.assertRaisesRegex(AssertionError, "four or more"):
            observed_members(frame_from_votes({"both": [-1, -1, 1, 1]}))

    def test_directional_null_retains_original_draws_and_event_margins(self):
        frame = frame_from_votes({"a": [-1, -1], "b": [-1, 0], "c": [1, 1], "d": [0, 1]})
        genes, units, analytic, margins = prepare_units(frame)
        self.assertEqual(margins[["proximal", "non_event", "distal"]].values.tolist(), [[2, 1, 1], [1, 1, 2]])
        original = numerical_function(SIM_SOURCE, "simulate")
        directional = simulate(units, len(genes), 10000, np.random.default_rng(83))
        total = original(units, len(genes), 10000, np.random.default_rng(83))
        np.testing.assert_array_equal(directional.sum(axis=1), total)
        # Enumerating label assignments checks the null rather than mirroring the simulator.
        left = set(itertools.permutations([-1, -1, 1, 0]))
        right = set(itertools.permutations([-1, 0, 1, 1]))
        exact = []
        for a, b in itertools.product(left, right):
            exact.append([sum(x == y == s for x, y in zip(a, b)) for s in [-1, 1]])
        exact_mean = np.mean(exact, axis=0)
        np.testing.assert_allclose([analytic["proximal"], analytic["distal"]], exact_mean)
        np.testing.assert_allclose(directional.mean(axis=0), exact_mean, atol=0.025)

    def test_bh_correction_uses_the_whole_family(self):
        np.testing.assert_allclose(bh([0.001, 0.01, 0.04, 0.2]), [0.004, 0.02, 0.16 / 3, 0.2])


if __name__ == "__main__":
    unittest.main()
