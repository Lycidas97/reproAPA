"""Intent tests: threshold agreement, fixed background and whole-scan inference."""
import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

SPEC = importlib.util.spec_from_file_location('threshold_scan', Path(__file__).with_name('analyze.py'))
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)


class Intent(unittest.TestCase):
    def base(self, low, high, negative_capacity=None, positive_capacity=None):
        return pd.DataFrame(dict(effect_min=low, effect_max=high,
            capacity_negative=negative_capacity if negative_capacity is not None else [.5]*len(low),
            capacity_positive=positive_capacity if positive_capacity is not None else [.5]*len(low)))

    def test_all_series_must_cross_the_threshold_even_when_the_mean_does(self):
        f = self.base([-.08, .005, -.07], [.001, .08, -.025])
        self.assertEqual(A.threshold_events(f, .02).event.tolist(), [0, 0, -1])

    def test_event_boundary_is_inclusive_without_borrowing_capacity_tolerance(self):
        f = self.base([-.02, .02, -.0199999995, .0199999995],
                      [-.02, .02, -.0199999995, .0199999995],
                      [.02, .02, .0199999995, .0199999995], [.02]*4)
        out = A.threshold_events(f, .02)
        self.assertEqual(out.event.tolist(), [-1, 1, 0, 0])
        self.assertTrue(out.can_negative.all())

    def test_a_capacity_guard_cannot_allow_an_infeasible_label(self):
        with self.assertRaises(AssertionError):
            A.threshold_events(self.base([-.05], [-.05], [.01], [.99]), .02)

    def test_high_threshold_does_not_remove_a_measured_nonrecurrent_pair(self):
        f = pd.DataFrame(dict(species=['human']*2+['mouse']*2,
            gene_id=['h']*2+['m']*2, gene_symbol=['H']*2+['M']*2, unit=['a', 'b', 'c', 'd'],
            effect_min=[.02]*4, effect_max=[.02]*4,
            capacity_negative=[.5]*4, capacity_positive=[.5]*4))
        orth = pd.DataFrame(dict(human_gene_id=['h'], mouse_gene_id=['m']))
        low, high = A.threshold_events(f, .01), A.threshold_events(f, .05)
        self.assertTrue(A.A.common_pairs(low, orth).equals(A.A.common_pairs(high, orth)))
        self.assertEqual(A.A.class_labels(high[high.species.eq('human')], ['h']).tolist(), [2])

    def test_missing_study_is_excluded_from_the_common_background(self):
        f = pd.DataFrame(dict(species=['human', 'mouse', 'mouse'],
            gene_id=['h', 'm', 'm'], unit=['a', 'c', 'd']))
        orth = pd.DataFrame(dict(human_gene_id=['h'], mouse_gene_id=['m']))
        self.assertEqual(len(A.A.common_pairs(f, orth)), 0)

    def test_two_same_direction_votes_allow_an_opposing_third_study(self):
        f = pd.DataFrame(dict(gene_id=['g']*3, gene_symbol=['G']*3,
                              unit=['a', 'b', 'c'], event=[-1, -1, 1]))
        self.assertEqual(A.A.class_labels(f, ['g']).tolist(), [0])

    def test_class_shuffle_keeps_every_species_stratum_margin(self):
        labels = np.array([0, 0, 1, 2, 2, 1, 2], np.int8)
        cells = [np.arange(5), np.arange(5, 6), np.arange(6, 7)]
        draws = A.A.draw_classes(labels, cells, 300, np.random.default_rng(1))
        for ii in cells:
            for label in range(3):
                self.assertTrue(((draws[:, ii] == label).sum(1) == (labels[ii] == label).sum()).all())
        self.assertTrue((draws[:, 5:] == labels[5:]).all())

    def test_additional_draws_keep_the_original_chunk_prefix_for_replay(self):
        h = np.array([0, 0, 1, 2], np.int8)
        m = np.array([0, 1, 2, 2], np.int8)
        cells = [np.arange(4)]
        short = A.A.cross_draws(h, m, cells, cells, np.random.default_rng(2), n=500)
        long = A.A.cross_draws(h, m, cells, cells, np.random.default_rng(2), n=1000)
        np.testing.assert_array_equal(short, long[:500])

    def test_scan_correction_cannot_use_only_the_selected_threshold(self):
        rows = []
        for model in ['measurement_PAS_capacity', 'full_assay']:
            for threshold in [.01, .015, .02, .025, .03, .035, .04, .045, .05]:
                for contrast in ['P_M', 'P_S']:
                    for h in ['proximal', 'distal']:
                        for m in ['proximal', 'distal']:
                            p = .001 if (model, threshold, contrast, h, m) == (
                                'full_assay', .02, 'P_M', 'proximal', 'proximal') else 1.
                            rows.append(dict(model=model, threshold=threshold, contrast=contrast,
                                human_direction=h, mouse_direction=m, tested=True, p_enrichment=p))
        out = A.apply_correction(pd.DataFrame(rows))
        hit = out[out.p_enrichment.eq(.001)].iloc[0]
        self.assertAlmostEqual(hit.q_eight_per_threshold_model, .008)
        self.assertAlmostEqual(hit.q_scan_72_per_model, .072)
        self.assertAlmostEqual(hit.q_scan_144_across_models, .144)

    def test_zero_null_tails_do_not_imply_zero_sampling_uncertainty(self):
        out = A.tail_precision(0, 100000)
        self.assertGreater(out['p_tail_Wilson95_hi'], 0)
        self.assertAlmostEqual(out['permutation_p_floor'], 1/100001)

    def test_constant_expression_bins_do_not_stand_in_for_threshold_capacity(self):
        f = pd.DataFrame(dict(gene_id=['a', 'b'], measured_pattern=['x;y']*2,
            proximal_capacity=[1, 2], distal_capacity=[2, 2], PAS_number_bin=[0, 0],
            mean_GEX_bin=[0, 0], APA_coverage_bin=[0, 0], baseline_bin=[0, 0]))
        self.assertEqual(len(A.A.partitions(f, 'full_assay')), 2)
        f['proximal_capacity'] = [2, 2]
        self.assertEqual(len(A.A.partitions(f, 'full_assay')), 1)


if __name__ == '__main__':
    unittest.main()
