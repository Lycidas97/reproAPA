"""Intent checks for sign, full axes, qualified scores and hierarchy."""
import importlib.util
from pathlib import Path
import unittest
import json
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('sm_analysis', Path(__file__).with_name('analyze.py'))
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class IntentTests(unittest.TestCase):
    def test_equal_target_dates_despite_capture_counts(self):
        # A date with three captures must not outweigh another date with one.
        q = pd.DataFrame(dict(time_a=['D1', 'D1', 'D1', 'D2'], sample_unit_a=['A', 'B', 'C', 'D'],
                              time_b=['S'] * 4, sample_unit_b=['source'] * 4))
        w = a.measurement.weights(q)
        np.testing.assert_allclose(w, [1/6, 1/6, 1/6, 1/2], atol=1e-15)
        np.testing.assert_allclose(a.hierarchical_weights(q), w, atol=1e-15)

    def test_multiple_reference_dates_units_and_captures(self):
        q = pd.DataFrame(dict(time_a=['D1'] * 4 + ['D2'] * 2,
                              sample_unit_a=['A'] * 3 + ['B', 'C', 'C'],
                              time_b=['D0', 'D0', 'D0', 'D0', 'D0', 'D0.5'],
                              sample_unit_b=['source1', 'source1', 'source2', 'source1', 'source1', 'source2']))
        expected = a.measurement.weights(q)
        np.testing.assert_allclose(a.hierarchical_weights(q), expected, rtol=0, atol=1e-15)

    def test_qualified_score_excludes_low_depth_MC(self):
        counts = np.array([[20, 0], [0, 20], [10, 10], [0, 0]])
        meta = pd.DataFrame(dict(n_cells=[10] * 4, Cycle=[1., 2., 3., 999.],
                                 Adjusted_Cycle=[.5, 1., 1.5, 999.]))
        result = a.profile(counts, np.arange(4), meta)
        self.assertEqual(result['n_mc'], 3)
        self.assertEqual(result['n_cells'], 30)
        self.assertEqual(result['Cycle'], 2.)
        self.assertEqual(result['Adjusted_Cycle'], 1.)
        np.testing.assert_allclose(result['usage'], [.5, .5])

    def test_S_M_reverses_APA_and_cycle_together(self):
        # Higher-cycle S is proximal relative to M; M target hierarchy is reversed.
        # A third unsupported PAS completes the full axis. It must not make .4 become 1.
        features = pd.DataFrame(dict(primary_gene_id=['g'] * 3, primary_gene_symbol=['G'] * 3,
                                     coordinate_eligible=[True] * 3, matched_source_count_50nt=[1, 1, 0],
                                     atlas_pas_id=['p1', 'p2', 'p3'], position=[0., .4, 1.]))
        z = {'GSE103221_features_json': np.array(features.to_json(orient='records')),
             'GSE103221_pas_data': np.array([20] * 7 + [30]),
             'GSE103221_pas_indices': np.array([4, 5, 6, 7, 0, 1, 2, 3]),
             'GSE103221_pas_indptr': np.array([0, 4, 7, 8]), 'GSE103221_pas_shape': np.array([8, 3])}
        meta = pd.DataFrame(dict(metacell_id=['m' + str(i) for i in range(8)], n_cells=[10] * 8,
                                 Cycle=[1., 1., 1., 999.] + [2.] * 4,
                                 Adjusted_Cycle=[.5, .5, .5, 999.] + [1.] * 4))
        m = dict(library_key='m', time='D1', sample_unit='M', replicate_status='unknown', rows=np.arange(4),
                 Cycle=250.5, Adjusted_Cycle=250.125)
        s = dict(library_key='s', time='D0', sample_unit='S', replicate_status='unknown', rows=np.arange(4, 8),
                 Cycle=2., Adjusted_Cycle=1.)
        group = {('CIP', 'm', 'M'): m, ('CIP', 's', 'S'): s}
        d = a.pair_design('GSE103221', 'CIP', m, s, 'main', 'M', 'S', 'documented_source')
        result, eligibility = a.effects(z, 'GSE103221', meta, group, [d], {'g'})
        self.assertEqual(len(result), 1)
        self.assertTrue(eligibility.eligible.all())
        self.assertEqual(result.delta_DWUI.iloc[0], -.4)
        self.assertEqual(result.delta_Cycle_exact.iloc[0], 1.)
        self.assertEqual(result.delta_Adjusted_Cycle_exact.iloc[0], .5)
        self.assertEqual(result.n_qualified_MC_a.iloc[0], 3)
        self.assertNotIn('m3', json.loads(result.qualified_MC_ids_a.iloc[0]))

    def test_failure_is_not_zero_measurement(self):
        q = a.summary(dict(route='missing'), pd.DataFrame())
        self.assertEqual(q['status'], 'not_estimable')
        self.assertEqual(q['N'], 0)
        self.assertTrue(np.isnan(q['mean_DWUI']))
        self.assertTrue(np.isnan(q['proximal']))

    def test_gene_threshold_boundaries(self):
        q = pd.DataFrame(dict(delta_DWUI=[-.02, .02, -.019, .019, 0.]))
        for field in a.VALUE_COLS[1:]:
            q[field] = 0.
        row = a.summary({}, q)
        self.assertEqual(row['proximal'], 1)
        self.assertEqual(row['distal'], 1)
        self.assertEqual(row['below_threshold'], 3)

    def test_two_PAS_require_distinct_state_support(self):
        first = dict(n_mc=3, n_cells=30, molecules=np.array([30, 0]), nonzero_mc=np.array([3, 0]))
        second = dict(n_mc=3, n_cells=30, molecules=np.array([0, 30]), nonzero_mc=np.array([0, 3]))
        self.assertEqual(a.measurement.reason(first, second), '')
        unsupported = dict(n_mc=3, n_cells=30, molecules=np.array([0, 30]), nonzero_mc=np.array([0, 1]))
        self.assertEqual(a.measurement.reason(first, unsupported), 'fewer_than_two_joint_supported_PAS')


if __name__ == '__main__':
    unittest.main()
