"""Intent: fixed margins, frozen ties, conservative cores and reciprocal gap-free coordinates."""
import unittest

import numpy as np
import pandas as pd

from analyze_intervals import core_interval, inclusion
from analyze_sequence import bh, bins, conditional_draws
from chain_mapper import map_chain


class SequenceIntent(unittest.TestCase):
    def test_identical_covariates_cannot_be_split_to_create_exchangeability(self):
        result = bins(pd.Series([0., 0., 0., 1., 1., 1.]))
        self.assertEqual(result.iloc[:3].nunique(), 1)
        self.assertEqual(result.iloc[3:].nunique(), 1)

    def test_null_preserves_case_counts_and_excludes_frozen_stratum_variance(self):
        frame = pd.DataFrame(dict(score=[0., 1., 10.], case=[True, False, True], layer=['mobile', 'mobile', 'frozen']))
        draws, expected, variance, margins = conditional_draws(frame, 'layer', np.random.default_rng(7), 1000)
        self.assertEqual(set(draws), {5., 5.5})
        self.assertAlmostEqual(expected, 5.25)
        self.assertAlmostEqual(variance, .0625)
        self.assertEqual(margins.cases.sum(), 2)

    def test_bh_keeps_ne_hypotheses_as_p_one_in_full_family(self):
        q = bh([.01, 1., 1., 1., 1., 1., 1., 1.])
        self.assertAlmostEqual(q[0], .08)
        self.assertTrue(np.all(q[1:] == 1))


class IntervalIntent(unittest.TestCase):
    def test_core_is_shared_across_native_boundary_uncertainty(self):
        self.assertEqual(core_interval('+', [100, 110], [200, 210]), (110, 200))
        self.assertEqual(core_interval('-', [200, 210], [100, 110]), (110, 200))

    def test_full_core_inclusion_uses_native_cuts_without_matched_renormalization(self):
        self.assertAlmostEqual(inclusion('+', 110, 200, [100, 150, 200, 210], [.3, .2, .4, .1]), .5)
        self.assertAlmostEqual(inclusion('-', 110, 200, [100, 110, 150, 210], [.2, .3, .1, .4]), .5)

    def test_minus_chain_flips_both_interval_and_feature_strand(self):
        lines = ['chain 1 chrH 1000 + 100 120 chrM 2000 - 200 220 1', '20']
        q = dict(query_id='seed', chrom='chrH', start0=105, end0=113, strand='+')
        result = map_chain(lines, [q])['seed']
        self.assertEqual(result['status'], 'unique')
        self.assertEqual(result['hits'][0], dict(chrom='chrM', strand='-', start0=1787, end0=1795, chain_id='1'))

    def test_seed_straddling_any_alignment_gap_is_unjudgeable(self):
        lines = ['chain 1 chrH 1000 + 100 122 chrM 2000 + 200 220 1', '10 2 0', '10']
        q = dict(query_id='seed', chrom='chrH', start0=105, end0=113, strand='+')
        self.assertEqual(map_chain(lines, [q])['seed']['status'], 'gap_or_partial')

    def test_multiple_chain_hits_cannot_be_best_score_selected(self):
        lines = ['chain 100 chrH 1000 + 100 120 chrM 2000 + 200 220 1', '20', '',
                 'chain 1 chrH 1000 + 100 120 chrM 2000 + 300 320 2', '20']
        q = dict(query_id='seed', chrom='chrH', start0=105, end0=113, strand='+')
        self.assertEqual(map_chain(lines, [q])['seed']['status'], 'ambiguous')

    def test_adjacent_blocks_with_no_gap_are_one_valid_alignment(self):
        lines = ['chain 1 chrH 1000 + 100 120 chrM 2000 + 200 220 1', '10 0 0', '10']
        q = dict(query_id='seed', chrom='chrH', start0=105, end0=113, strand='+')
        self.assertEqual(map_chain(lines, [q])['seed']['status'], 'unique')


if __name__ == '__main__':
    unittest.main()
