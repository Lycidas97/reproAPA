"""Intent checks for the interpretation limits of site correspondence."""
import unittest
from types import SimpleNamespace
import pandas as pd

from analyze_sites import average_profiles, characterize, matched_mass, mean_supports_drivers, merge_native, strict_utr, supports_drivers, top_drivers


class SiteIntent(unittest.TestCase):
    def test_unmapped_maximum_cannot_be_replaced_by_mappable_runner_up(self):
        gain, loss = top_drivers({'UNMAPPED:native': .4, 'mapped_gain': .2, 'mapped_loss': -.6})
        self.assertEqual(gain, 'UNMAPPED:native')
        self.assertEqual(loss, 'mapped_loss')

    def test_merged_native_sites_do_not_gain_extra_event_votes(self):
        self.assertEqual(merge_native(['a', 'b', 'c'], [.2, .3, .5], {'a': 'shared', 'b': 'shared'}),
                         {'shared': .5, 'UNMAPPED:c': .5})

    def test_studies_are_equal_despite_number_of_measured_series(self):
        profiles = [dict(study_id=study, target={'gain': effect, 'loss': 1 - effect}, reference={'gain': 0., 'loss': 1.})
                    for study, effect in [('Xing', .9), ('Guan', .1), ('Guan', .2)]]
        delta, studies = average_profiles(profiles)
        self.assertAlmostEqual(delta['gain'], (.9 + (.1 + .2) / 2) / 2)
        self.assertEqual(len(studies), 2)

    def test_two_conditions_or_series_are_one_independent_study(self):
        profiles = [dict(study_id='same_study', target={'gain': .1, 'loss': .9}, reference={'gain': 0., 'loss': 1.}) for _ in range(2)]
        self.assertEqual(supports_drivers(profiles, 'gain', 'loss'), 1)

    def test_opposing_series_prevents_exact_study_event_support(self):
        profiles = [dict(study_id='Guan', target={'gain': effect, 'loss': 1 - effect}, reference={'gain': .1, 'loss': .9}) for effect in [.3, .0]]
        self.assertEqual(supports_drivers(profiles, 'gain', 'loss'), 0)

    def test_retained_mass_checks_each_state_not_only_pooled_mean(self):
        profiles = [dict(target={'matched': .81, 'other': .19}, reference={'matched': .79, 'other': .21})]
        self.assertAlmostEqual(matched_mass(profiles, ['matched']), .79)

    def test_all_contributing_full_axes_must_be_utr_compatible(self):
        e = {('human', study, 'gene'): SimpleNamespace(tandem_utr_compatible=ok, compatible_transcripts='tx', native_axis_span_nt=1000)
             for study, ok in [('a', True), ('b', False)]}
        result = strict_utr('human', 'gene', [dict(study_id='a'), dict(study_id='b')], e)
        self.assertFalse(result[0])
        self.assertEqual(result[1], 'non_tandem_full_primary_axis')

    def test_common_transcript_cannot_be_inferred_from_separate_transcripts(self):
        e = {('human', study, 'gene'): SimpleNamespace(tandem_utr_compatible=True, compatible_transcripts=study, native_axis_span_nt=1000) for study in ['a', 'b']}
        result = strict_utr('human', 'gene', [dict(study_id='a'), dict(study_id='b')], e)
        self.assertFalse(result[0])
        self.assertEqual(result[1], 'no_common_terminal_UTR_transcript')

    def test_missing_naive_qualification_is_unavailable_not_inherited(self):
        result = strict_utr('human', 'gene', [dict(study_id='Liu', pair_tandem=True)], {})
        self.assertFalse(result[0])
        self.assertEqual(result[1], 'missing_frozen_full_axis_UTR_qualification')

    def test_study_mean_support_does_not_replace_original_series_vote(self):
        profiles = [dict(study_id='Guan', target={'gain': effect, 'loss': 1 - effect}, reference={'gain': .1, 'loss': .9}) for effect in [.15, .10]]
        self.assertEqual(supports_drivers(profiles, 'gain', 'loss'), 0)
        self.assertEqual(mean_supports_drivers(profiles, 'gain', 'loss'), 1)

    def test_naive_other_state_key_cannot_replace_unverified_actual_axis(self):
        e = {('human', 'Liu', 'gene'): SimpleNamespace(tandem_utr_compatible=True, compatible_transcripts='tx', native_axis_span_nt=1000)}
        result = strict_utr('human', 'gene', [dict(study_id='Liu', axis_verified=False)], e)
        self.assertFalse(result[0])
        self.assertEqual(result[1], 'actual_naive_axis_not_verified')

    def test_parent_case_cannot_use_nonshared_contrast_as_sequence_control(self):
        catalog = pd.DataFrame([dict(id=f'{sp}{n}', cut_site=1000 + 1000*n, strand='+') for sp in ['h', 'm'] for n in [0, 1]]).set_index('id')
        pairs = pd.DataFrame([dict(human_gene_id='hg', mouse_gene_id='mg', human_id=f'h{n}', mouse_id=f'm{n}', reciprocal25=True) for n in [0, 1]])
        profiles = {sp: [dict(study_id=study, series=study, target={prefix+'0': .6, prefix+'1': .4}, reference={prefix+'0': .4, prefix+'1': .6}, original_delta=-.2) for study in ['a', 'b']] for sp, prefix in [('human', 'h'), ('mouse', 'm')]}
        eligibility = {(sp, study, gene): SimpleNamespace(tandem_utr_compatible=True, compatible_transcripts='tx', native_axis_span_nt=1000) for sp, gene in [('human', 'hg'), ('mouse', 'mg')] for study in ['a', 'b']}
        frame = dict(human_gene_id='hg', mouse_gene_id='mg', case=True)
        row, _ = characterize(frame, 'main', 'P_M', 'single_species_recurrent', 'human', profiles, eligibility, catalog, pairs, 25)
        self.assertTrue(row['strict_double_driver'])
        self.assertFalse(row['sequence_candidate'])


if __name__ == '__main__':
    unittest.main()
