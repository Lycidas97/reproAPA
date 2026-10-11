"""Intent tests for reference definitions, support, and measurement identity."""
import unittest

import numpy as np
import pandas as pd

from analyze import aggregate_bridge, aggregate_new, bridge_identities, classify, dense_gene, gradient_aggregate, original_groups, programme, reason, remove_failed_support_pairs, state, summaries, weights


class IntentTests(unittest.TestCase):
    def test_unique_support_loss_rematches_each_background_without_reusing_common_flags(self):
        rows = []
        profiles = {'D00': ([.2, .8], [.8, .2]), 'D01': ([.2, .8], [.5, .5]),
                    'D10': ([.3, .7], [.8, .2]), 'D11': ([.3, .7], [.5, .5])}
        for library in ['A', 'B']:
            for contrast, (a,b) in profiles.items():
                rows.append(dict(study_id='S', route='R', gene_id='G', library_key=library, setting='supported50', version='core_third_strict',
                    contrast=contrast, time_a=library, time_b=library, sample_unit_a=library, sample_unit_b=library,
                    target_usage=np.array(a), reference_usage=np.array(b), position=np.array([0.,1.]), PAS_ids='["a","b"]',
                    delta_DWUI=a[1]-b[1], full_weight=.5, common_weight=.5, common_four=True))
        pairs = pd.DataFrame(rows)
        lost = pairs[pairs.library_key.eq('B') & pairs.contrast.eq('D01')][['study_id','route','library_key','setting','version','contrast','gene_id']]
        kept = remove_failed_support_pairs(pairs, lost)
        self.assertEqual(len(kept),7)
        self.assertFalse(kept.common_four.any())
        self.assertTrue(kept.common_weight.isna().all())
        full = aggregate_new(kept,False); common = aggregate_new(kept,True)
        self.assertEqual(full.set_index('contrast').loc['D00','n_libraries'],2)
        self.assertTrue(common.n_libraries.eq(1).all())
        self.assertFalse(kept.loc[kept.library_key.eq('B'),'common_four'].any())
        m,_,_=aggregate_bridge(kept,'common_M_bridge',['D00','D01'])
        p,_,_=aggregate_bridge(kept,'common_P_bridge',['D00','D10'])
        self.assertTrue(m.n_libraries.eq(1).all())
        self.assertTrue(p.n_libraries.eq(2).all())

    def test_route_summary_uses_route_contributors_without_library_key_column(self):
        q = pd.DataFrame([dict(study_id='GSE103221', route='CIP', contrast='D00', delta_DWUI=-.03,
                               n_libraries=2, contributor_libraries='["A","B"]'),
                          dict(study_id='GSE103221', route='CIP', contrast='D00', delta_DWUI=.03,
                               n_libraries=2, contributor_libraries='["B","C"]')])
        rows = summaries(q, 'common_four', 'core_third_strict', 'supported50')
        main = next(r for r in rows if r['route']=='CIP' and r['contrast']=='D00' and r['threshold']==.02)
        self.assertEqual(main['N'], 2)
        self.assertEqual(main['n_contributor_libraries'], 3)
        self.assertEqual(main['max_libraries_per_gene'], 2)
        self.assertEqual(main['n_proximal'], 1)
        absent = next(r for r in rows if r['route']=='WOT' and r['contrast']=='D00' and r['threshold']==.02)
        self.assertEqual(absent['status'], 'not_estimable')
        self.assertEqual(absent['N'], 0)

    def test_marker_constants_keep_the_frozen_panel_denominator(self):
        cp = np.array([[0., 0., 2., 5.], [1., 2., 2., 5.], [3., 4., 2., 5.]])
        score, nvariable, expressed = programme(cp, [0, 1, 2, 3])
        values = np.log1p(cp[:, :2])
        expected = ((values - values.mean(0)) / values.std(0, ddof=0)).sum(1) / 4
        np.testing.assert_allclose(score, expected)
        self.assertEqual(nvariable, 2)
        np.testing.assert_array_equal(expressed, [2, 4, 4])

    def test_tied_programmes_cannot_be_split_to_manufacture_states(self):
        result, reasons, cuts, available = classify(np.ones(12), np.arange(12),
                                                  np.full(12, 3), np.zeros(12), 1/3, True)
        self.assertFalse(available)
        self.assertTrue((result == 'unassigned/mixed').all())
        self.assertTrue((reasons == 'programme_quantiles_not_separated').all())

    def test_quantile_separation_does_not_replace_absolute_expression(self):
        scores = np.arange(12, dtype=float)
        result, reasons, cuts, available = classify(scores, -scores, np.ones(12),
                                                  np.ones(12), 1/3, True)
        self.assertTrue(available)
        self.assertFalse(np.isin(result, ['Pstar', 'Mstar']).any())
        self.assertTrue((reasons == 'absolute_marker_support_not_passed').any())

    def test_literal_two_thirds_boundary_is_inclusive_without_float_complement(self):
        scores = np.arange(13, dtype=float)
        result, reasons, cuts, available = classify(scores, -scores, np.full(13, 2),
                                                  np.zeros(13), 1/3, True)
        self.assertEqual(cuts[1], 8.)
        self.assertEqual(result[8], 'Pstar')

    def test_usage_is_mc_equal_and_low_depth_cannot_change_the_mean(self):
        counts = np.array([[100, 0], [10, 10], [0, 19]])
        meta = pd.DataFrame({'n_cells': [20, 20, 20]})
        value = state(counts, np.arange(3), meta)
        np.testing.assert_allclose(value['usage'], [.75, .25])
        self.assertEqual(value['n_mc'], 2)
        self.assertEqual(value['n_cells'], 40)
        self.assertNotAlmostEqual(value['usage'][0], 110/120)

    def test_joint_support_is_a_gate_and_does_not_trim_usage_or_rescale_axis(self):
        counts = np.array([[10, 10, 1]] * 6)
        meta = pd.DataFrame({'n_cells': [20] * 6})
        a = state(counts, np.arange(3), meta)
        b = state(counts, np.arange(3, 6), meta)
        self.assertEqual(reason(a, b), '')
        self.assertEqual(len(a['usage']), 3)
        self.assertAlmostEqual(a['usage'] @ np.array([0., .5, 1.]), 6/21)
        # A retained subset still uses the original complete primary positions.
        self.assertAlmostEqual(np.array([.5, .5]) @ np.array([.5, 1.]), .75)

    def test_date_unit_weighting_does_not_reward_extra_captures(self):
        q = pd.DataFrame({'time_a': ['D1', 'D1', 'D1', 'D1', 'D2'],
                          'time_b': ['D1', 'D1', 'D1', 'D1', 'D2'],
                          'sample_unit_a': ['A', 'A', 'A', 'B', 'C'],
                          'sample_unit_b': ['A', 'A', 'A', 'B', 'C']})
        np.testing.assert_allclose(weights(q), [1/12, 1/12, 1/12, 1/4, 1/2])

    def test_full_pool_labels_do_not_add_original_baseline_contributors(self):
        meta = pd.DataFrame({'in_planned_route': [True] * 6, 'original_selected': [True] * 3 + [False] * 3,
                             'original_purity': [1.] * 6, 'original_node': ['P'] * 6,
                             'route': ['X'] * 6, 'library_key': ['L'] * 6, 'n_cells': [20] * 6})
        np.testing.assert_array_equal(original_groups(meta)[('X', 'L', 'original_P')], [0, 1, 2])

    def test_sparse_pas_replay_preserves_integer_counts_and_mc_row_order(self):
        result = dense_gene(np.array([7, 3, 11]), np.array([0, 2, 1]), np.array([0, 2, 3]), 3, [1, 0])
        np.testing.assert_array_equal(result, [[0, 7], [11, 0], [0, 3]])
        self.assertTrue(np.issubdtype(result.dtype, np.integer))

    def test_four_comparisons_share_exact_libraries_and_weights_in_common_background(self):
        profiles = {'D00': ([.2, .8], [.8, .2]), 'D01': ([.2, .8], [.5, .5]),
                    'D10': ([.3, .7], [.8, .2]), 'D11': ([.3, .7], [.5, .5])}
        rows = []
        for contrast, (a, b) in profiles.items():
            rows.append(dict(study_id='S', route='R', gene_id='G', library_key='A', contrast=contrast,
                time_a='D1', time_b='D1', sample_unit_a='U1', sample_unit_b='U1',
                target_usage=np.array(a), reference_usage=np.array(b), position=np.array([0., 1.]),
                PAS_ids='["a","b"]', delta_DWUI=a[1]-b[1], full_weight=np.nan, common_weight=np.nan, common_four=False))
        rows.append(dict(rows[0], library_key='B', time_a='D2', time_b='D2', sample_unit_a='U2', sample_unit_b='U2'))
        q = pd.DataFrame(rows)
        aggregate_new(q, common=False)
        common = aggregate_new(q, common=True)
        self.assertEqual(len(common), 4)
        self.assertTrue(common.n_libraries.eq(1).all())
        self.assertTrue(q.loc[q.library_key.eq('A'), 'common_four'].all())
        self.assertFalse(q.loc[q.library_key.eq('B'), 'common_four'].any())
        np.testing.assert_allclose(q.loc[q.library_key.eq('A'), 'common_weight'], 1.)
        self.assertAlmostEqual(q.loc[0, 'full_weight'], .5)
        identities = bridge_identities(q, common, 'test', 'supported50')
        self.assertTrue(all(row['abs_error'] < 1e-12 for row in identities))

    def test_M_gradient_changes_only_reference_usage_with_shared_P(self):
        rows = []
        for library, p, time in [('A', [.2, .8], 'D1'), ('B', [.1, .9], 'D2')]:
            rows.append(dict(study_id='S', route='R', panel='core', gene_id='G', library_key=library,
                time_a=time, time_b=time, sample_unit_a=library, sample_unit_b=library,
                P_usage=np.array(p), M_low_usage=np.array([.8, .2]), M_middle_usage=np.array([.65, .35]),
                M_high_usage=np.array([.5, .5]), M_whole_usage=np.array([.65, .35]), position=np.array([0., 1.]), weight=np.nan))
        result = gradient_aggregate(pd.DataFrame(rows))
        self.assertAlmostEqual(result.difference_high_minus_low.iloc[0], -.3)
        self.assertAlmostEqual(result.direct_M_low_minus_high.iloc[0], -.3)
        self.assertLess(result.identity_abs_error.iloc[0], 1e-12)

    def test_depth_and_group_size_cannot_replace_two_supported_PAS(self):
        counts = np.array([[100, 0]] * 6)
        meta = pd.DataFrame({'n_cells': [50] * 6})
        self.assertEqual(reason(state(counts, np.arange(3), meta), state(counts, np.arange(3, 6), meta)),
                         'fewer_than_two_joint_supported_PAS')

    def test_M_bridge_retains_matching_libraries_without_Pstar(self):
        rows = []
        for library, contrasts in [('A', ['D00', 'D01', 'D10', 'D11']), ('B', ['D00', 'D01'])]:
            profiles = {'D00': ([.2, .8], [.8, .2]), 'D01': ([.2, .8], [.5, .5]),
                        'D10': ([.3, .7], [.8, .2]), 'D11': ([.3, .7], [.5, .5])}
            for contrast in contrasts:
                a, b = profiles[contrast]
                rows.append(dict(study_id='S', route='R', gene_id='G', library_key=library, contrast=contrast,
                    time_a=library, time_b=library, sample_unit_a=library, sample_unit_b=library,
                    target_usage=np.array(a), reference_usage=np.array(b), position=np.array([0., 1.]),
                    PAS_ids='["a","b"]', delta_DWUI=a[1]-b[1]))
        pairs = pd.DataFrame(rows)
        m, mw, mc = aggregate_bridge(pairs, 'common_M_bridge', ['D00', 'D01'])
        p, pw, pc = aggregate_bridge(pairs, 'common_P_bridge', ['D00', 'D10'])
        self.assertTrue(m.n_libraries.eq(2).all())
        self.assertTrue(p.n_libraries.eq(1).all())
        self.assertEqual(len(mw), 4)
        self.assertEqual(len(pw), 2)
        self.assertTrue(all(abs(row['shared_weight'] - .5) < 1e-12 for row in mw))
        self.assertTrue(all(row['identity_abs_error'] < 1e-12 for row in mc + pc))


if __name__ == '__main__':
    unittest.main()
