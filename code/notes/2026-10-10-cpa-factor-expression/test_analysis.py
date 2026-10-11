"""Tests encode physical-library and reference restrictions, not RNA conclusions."""
import unittest
import numpy as np
import pandas as pd
from analyze import legacy_library_rows, comparison, documented_source_flags, expression_contrasts


class ExpressionContractTests(unittest.TestCase):
    def test_all_gene_denominator_gate_and_zero_expression_retained(self):
        meta = pd.DataFrame(dict(route=['r']*5, library_key=['a']*3+['b']*2, node=['P']*5, experimental_time=['D10']*5, n_cells=[10]*5))
        counts = np.array([[10, 0], [10, 0], [10, 0], [5, 0], [5, 0]])
        totals = np.array([1000]*5)
        result = legacy_library_rows()(meta, counts, totals, ['NUDT21', 'PCF11'], 'study', 'human')
        a = result[result.library_key.eq('a') & result.factor.eq('NUDT21')].iloc[0]
        self.assertTrue(a.qualified)
        self.assertAlmostEqual(a.mean_log2CPM1, np.log2(10001))
        self.assertEqual(a.total_gex, 3000)
        self.assertFalse(result[result.library_key.eq('b')].qualified.any())
        zero = result[result.factor.eq('PCF11') & result.library_key.eq('a')].iloc[0]
        self.assertTrue(zero.mapped)
        self.assertEqual(zero.mean_log2CPM1, 0)

    def test_physical_library_pairing_excludes_unmatched_extremes(self):
        target = pd.DataFrame(dict(library_key=['shared', 'P_only'], time=['D1', 'D1'], mean_log2CPM1=[3., 100.], pooled_CPM=[7., 1000.]))
        reference = pd.DataFrame(dict(library_key=['shared', 'M_only'], time=['D1', 'D1'], mean_log2CPM1=[1., -100.], pooled_CPM=[3., 0.]))
        effect, pairs = comparison(('human', 's', 'r', 'NUDT21'), 'P_M', target, reference, paired=True)
        self.assertEqual(effect['delta_log2CPM1'], 2)
        self.assertEqual(effect['n_pairs'], 1)
        self.assertEqual(len(pairs), 1)

    def test_library_equal_weight_and_time_equal_sensitivity_are_distinct(self):
        target = pd.DataFrame(dict(library_key=['a', 'b', 'c'], time=['D1', 'D1', 'D2'], mean_log2CPM1=[1., 1., 4.], pooled_CPM=[1., 1., 4.], n_metacells=[100, 3, 3]))
        ref = pd.DataFrame(dict(library_key=['s'], time=['D0'], mean_log2CPM1=[0.], pooled_CPM=[0.]))
        effect, _ = comparison(('human', 's', 'r', 'NUDT21'), 'P_S', target, ref)
        self.assertEqual(effect['delta_log2CPM1'], 2)
        self.assertEqual(effect['time_equal_delta'], 2.5)

    def test_observed_source_is_not_documented_starting_source(self):
        libraries = pd.DataFrame(dict(study_id=['Xing', 'Xing', 'Liu'], route=['r']*3, library_key=['D0', 'D2', 'pool'], state=['S']*3))
        design = libraries.assign(time=['D0', 'D2', 'D0_D7_pool'])
        z = documented_source_flags(libraries, design)
        self.assertEqual(z.documented_source.tolist(), [True, False, False])

    def test_naive_primed_reference_cannot_be_replaced_by_RSeT_P(self):
        z = pd.DataFrame(dict(species=['human']*4, study_id=['GSE150311']*4, route=['GSE150311']*4, factor=['NUDT21']*4, library_key=['GSM4546546', 'GSM4546548', 'GSM4546547', 'GSM4546548'], state=['N', 'N', 'P', 'P'], time=['endpoint']*4, mean_log2CPM1=[2., 4., 5., 100.], pooled_CPM=[2., 4., 5., 100.], qualified=[True]*4, mapped=[True]*4, documented_source=[False]*4))
        effects, pairs = expression_contrasts(z)
        effect = effects[effects.contrast.eq('hN_P_primed_condition')].iloc[0]
        self.assertEqual(effect.delta_log2CPM1, -2)
        self.assertEqual(effect.n_target_libraries, 2)
        self.assertEqual(effect.n_reference_libraries, 1)
        same = effects[effects.contrast.eq('hN_P_same_library')].iloc[0]
        self.assertEqual(same.n_pairs, 1)
        self.assertEqual(same.delta_log2CPM1, -96)

    def test_missing_one_naive_condition_invalidates_pooled_complete_contrast(self):
        z = pd.DataFrame(dict(species=['human']*2, study_id=['GSE150311']*2, route=['GSE150311']*2, factor=['NUDT21']*2, library_key=['GSM4546548', 'GSM4546547'], state=['N', 'P'], time=['endpoint']*2, mean_log2CPM1=[4., 5.], pooled_CPM=[4., 5.], qualified=[True]*2, mapped=[True]*2, documented_source=[False]*2))
        effects, _ = expression_contrasts(z)
        effect = effects[effects.contrast.eq('hN_P_primed_condition')].iloc[0]
        self.assertEqual(effect.status, 'unavailable')
        self.assertTrue(np.isnan(effect.delta_log2CPM1))


if __name__ == '__main__':
    unittest.main()
