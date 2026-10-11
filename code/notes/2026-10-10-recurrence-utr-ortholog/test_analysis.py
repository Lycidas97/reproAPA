"""Intent checks for annotation gates and conditional ortholog inference."""
import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

SPEC = importlib.util.spec_from_file_location('utr_ortholog_analysis', Path(__file__).with_name('analyze.py'))
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)


class Intent(unittest.TestCase):
    def test_entire_axis_must_be_tandem_in_every_measured_study(self):
        frame = pd.DataFrame(dict(species=['human']*4, unit=['a', 'b', 'a', 'b'],
                                  gene_id=['g', 'g', 'h', 'h'], event=[1, 1, -1, -1]))
        eligibility = pd.DataFrame(dict(species=['human']*4, study_id=['a', 'b', 'a', 'b'],
                                         gene_id=['g', 'g', 'h', 'h'],
                                         tandem_utr_compatible=[True, False, True, True]))
        self.assertEqual(A.eligible_tandem(frame, eligibility), {'h'})
        frame.event *= -1
        self.assertEqual(A.eligible_tandem(frame, eligibility), {'h'})

    def test_missing_study_is_not_a_measured_nonresponse(self):
        frame = pd.DataFrame(dict(species=['human', 'human', 'human', 'mouse', 'mouse', 'mouse', 'mouse'],
                                  gene_id=['h1', 'h1', 'h2', 'm1', 'm1', 'm2', 'm2'],
                                  unit=['a', 'b', 'a', 'c', 'd', 'c', 'd']))
        orth = pd.DataFrame(dict(human_gene_id=['h1', 'h2'], mouse_gene_id=['m1', 'm2']))
        result = A.common_pairs(frame, orth)
        self.assertEqual(result.human_gene_id.tolist(), ['h1'])

    def test_class_null_keeps_recurrence_margins(self):
        labels = np.array([0, 0, 1, 2, 2, 1, 2], np.int8)
        cells = [np.arange(5), np.arange(5, 7)]
        draws = A.draw_classes(labels, cells, 300, np.random.default_rng(1))
        for ii in cells:
            for label in range(3):
                self.assertTrue(((draws[:, ii] == label).sum(1) == (labels[ii] == label).sum()).all())
        self.assertGreater(np.unique(draws, axis=0).shape[0], 1)

    def test_ortholog_expectation_conditions_on_fixed_species_margins(self):
        h = np.array([0, 0, 1, 2], np.int8)
        m = np.array([0, 1, 2, 2], np.int8)
        cells = [np.arange(4)]
        expected = A.class_probabilities(h, cells).T @ A.class_probabilities(m, cells)
        np.testing.assert_allclose(expected, np.outer(np.bincount(h, minlength=3),
                                                    np.bincount(m, minlength=3))/4)
        draws = A.cross_draws(h, m, cells, cells, np.random.default_rng(2), n=2000)
        self.assertTrue((draws.sum((1, 2)) == 4).all())
        np.testing.assert_allclose(draws.mean(0), expected, atol=.07, rtol=0)

    def test_two_supports_can_allow_opposing_third_study(self):
        frame = pd.DataFrame(dict(gene_id=['g']*3, gene_symbol=['G']*3,
                                  unit=['a', 'b', 'c'], event=[-1, -1, 1]))
        self.assertEqual(A.class_labels(frame, ['g']).tolist(), [0])


if __name__ == '__main__':
    unittest.main()
