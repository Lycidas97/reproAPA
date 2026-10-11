"""Intent: time curves must preserve contributors and keep dates as test units."""
import importlib.util
import itertools
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('analysis', Path(__file__).with_name('analyze.py'))
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class AnalysisIntent(unittest.TestCase):
    def test_missing_planned_library_excludes_gene(self):
        frame = pd.DataFrame({'gene_id': ['complete', 'complete', 'missing'], 'library_key': ['a', 'b', 'a']})
        self.assertEqual(a.fixed_library_genes(frame, ['a', 'b']), {'complete'})
        self.assertEqual(a.fixed_library_genes(frame, ['a', 'b', 'unobserved']), set())

    def test_sample_units_then_dates_not_library_or_gene_pseudoreplication(self):
        frame = pd.DataFrame({'gene_id': ['g'] * 4, 'time': ['D1'] * 3 + ['D2'], 'sample_unit': ['unit1', 'unit1', 'unit2', 'unit3'], 'value': [0., 0., 1., 1.]})
        daily = a.daily_values(frame).set_index('time').value
        self.assertEqual(daily['D1'], .5)
        self.assertEqual(daily['D2'], 1.)
        self.assertEqual(daily.mean(), .75)

    def test_exact_date_null_keeps_values_and_can_detect_both_signs(self):
        x, y = np.arange(4), np.array([.1, .2, .3, .4])
        estimate, draws, p = a.date_permutation(x, y, 1, exact=True)
        expected = np.array([a.slope(x, v) for v in itertools.permutations(y)])
        np.testing.assert_allclose(draws, expected)
        self.assertAlmostEqual(p, 2 / 24)
        negative, negative_draws, negative_p = a.date_permutation(x, -y, 1, exact=True)
        self.assertAlmostEqual(negative, -estimate)
        self.assertEqual(negative_p, p)
        np.testing.assert_allclose(negative_draws, -draws)

    def test_only_documented_numeric_days_can_be_temporal_clock(self):
        self.assertEqual(a.numeric_day('D8.25'), 8.25)
        self.assertTrue(np.isnan(a.numeric_day('endpoint')))
        self.assertTrue(np.isnan(a.numeric_day('Stage_IV_day_10')))
        self.assertTrue(np.isnan(a.numeric_day('D0_D7_pool')))

    def test_bh_controls_two_predefined_primary_curves(self):
        np.testing.assert_allclose(a.bh([.03, .2]), [.06, .2])


if __name__ == '__main__':
    unittest.main(verbosity=2)
