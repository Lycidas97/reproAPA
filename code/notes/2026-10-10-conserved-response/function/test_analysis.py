"""Intent checks: annotation unknowns, opportunity, ties and conditional inference."""
import itertools
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import analyze as a


class IntentTests(unittest.TestCase):
    def test_gene_gets_one_vote_across_shared_contrasts(self):
        f = pd.DataFrame(dict(P_S_measurable=[True], P_M_measurable=[True],
            P_S_human_direction=['proximal'], P_S_mouse_direction=['proximal'],
            P_M_human_direction=['proximal'], P_M_mouse_direction=['proximal']))
        g = a.classify_same_contrast(f)
        self.assertEqual(int(g.case_all_four.sum()), 1)
        self.assertFalse(g.control.any())

    def test_recurrence_in_different_contrasts_is_not_single_species(self):
        f = pd.DataFrame(dict(P_S_measurable=[True], P_M_measurable=[True],
            P_S_human_direction=['proximal'], P_S_mouse_direction=['nonrecurrent'],
            P_M_human_direction=['nonrecurrent'], P_M_mouse_direction=['proximal']))
        g = a.classify_same_contrast(f)
        self.assertFalse(g.case_all_four.any())
        self.assertFalse(g.control.any())

    def test_unmeasured_contrast_cannot_become_negative_evidence(self):
        f = pd.DataFrame(dict(P_S_measurable=[True], P_M_measurable=[False],
            P_S_human_direction=['distal'], P_S_mouse_direction=['distal'],
            P_M_human_direction=['proximal'], P_M_mouse_direction=['nonrecurrent']))
        g = a.classify_same_contrast(f)
        self.assertEqual(g.direction_class.iloc[0], 'distal')
        self.assertTrue(pd.isna(g.P_M_human_direction.iloc[0]))
        self.assertFalse(g.control.any())

    def test_opposing_bilateral_response_is_not_control(self):
        f = pd.DataFrame(dict(P_S_measurable=[True], P_M_measurable=[False],
            P_S_human_direction=['proximal'], P_S_mouse_direction=['distal'],
            P_M_human_direction=[None], P_M_mouse_direction=[None]))
        g = a.classify_same_contrast(f)
        self.assertFalse(g.case_all_four.any()); self.assertFalse(g.control.any())

    def test_first_three_sensitivity_excludes_only_PM_distal(self):
        f = pd.DataFrame(dict(P_S_measurable=[False], P_M_measurable=[True],
            P_S_human_direction=[None], P_S_mouse_direction=[None],
            P_M_human_direction=['distal'], P_M_mouse_direction=['distal']))
        g = a.classify_same_contrast(f)
        self.assertTrue(g.case_all_four.iloc[0]); self.assertFalse(g.case_first_three.iloc[0])

    def test_ties_at_cut_stay_together_and_missing_is_explicit(self):
        vals = pd.Series([0., .3, .3, .8, 1., np.nan])
        b = a.bins(vals)
        self.assertEqual(b.iloc[1], b.iloc[2]); self.assertEqual(b.iloc[1], 0)
        self.assertEqual(b.iloc[-1], -1)
        self.assertEqual(a.bins(pd.Series([.5]*9)).tolist(), [0]*9)

    def test_ontology_includes_part_of_not_regulates_and_resolves_alias(self):
        text = 'format-version: 1.2\n\n[Term]\nid: GO:0000001\nname: root\nnamespace: biological_process\n\n[Term]\nid: GO:0000002\nname: child\nnamespace: biological_process\nalt_id: GO:9999999\nrelationship: part_of GO:0000001 ! root\n\n[Term]\nid: GO:0000003\nname: regulator\nnamespace: biological_process\nrelationship: regulates GO:0000001 ! root\n\n[Term]\nid: GO:0000004\nname: obsolete\nnamespace: biological_process\nis_obsolete: true\n'
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'ontology.obo'; p.write_text(text)
            terms, aliases, ancestors = a.ontology(p)
        self.assertIn('GO:0000001', ancestors('GO:0000002'))
        self.assertNotIn('GO:0000001', ancestors('GO:0000003'))
        self.assertEqual(aliases['GO:9999999'], 'GO:0000002')
        self.assertNotIn('GO:0000004', terms)

    def test_unknown_namespace_does_not_count_as_nonmember(self):
        f = pd.DataFrame(dict(case=[True, True, False, False], coarse_stratum=['a']*4))
        row, qualified, _, _, _ = a.test_category(f, np.array([True, False, True, True]),
            np.array([True, False, True, False]), 'coarse')
        self.assertEqual(row['annotation_unknown_cases'], 1)
        self.assertEqual((row['N'], row['cases'], row['K'], row['hits']), (3, 1, 2, 1))
        self.assertEqual(qualified.coarse_stratum.tolist(), ['a']*3)
        self.assertAlmostEqual(row['p'], 2/3)

    def test_exact_stratified_distribution_matches_full_enumeration(self):
        populations = [[1, 1, 0, 0], [1, 0, 0]]; sizes = [2, 1]
        counts = []
        for draw1 in itertools.combinations(range(4), 2):
            for draw2 in itertools.combinations(range(3), 1):
                counts.append(sum(populations[0][i] for i in draw1)+sum(populations[1][i] for i in draw2))
        offset, pmf = a.exact_distribution([(4, 2, 2), (3, 1, 1)])
        empirical = np.bincount(counts, minlength=4)/len(counts)
        np.testing.assert_allclose(pmf, empirical[offset:offset+len(pmf)], rtol=0, atol=1e-14)

    def test_frozen_labels_report_NE_and_occupy_family_p_one(self):
        f = pd.DataFrame(dict(case=[True, True, False], coarse_stratum=['case_only', 'case_only', 'control_only']))
        row, _, _, _, _ = a.test_category(f, np.ones(3, bool), np.array([True, False, False]), 'coarse')
        self.assertTrue(row['status'].startswith('NE')); self.assertEqual(row['p'], 1.)
        self.assertEqual(row['mobile_cases'], 0)

    def test_BH_keeps_uninformative_planned_hypothesis_in_family(self):
        np.testing.assert_allclose(a.bh([.001, .03, .04, 1.]), [.004, .05333333333333334, .05333333333333334, 1.])

    def test_full_tail_is_a_valid_probability_despite_roundoff(self):
        f = pd.DataFrame(dict(case=[True]+[False]*8, coarse_stratum=['a']*9))
        row, _, _, _, _ = a.test_category(f, np.ones(9, bool), np.array([False]+[True]*4+[False]*4), 'coarse')
        self.assertEqual(row['hits'], 0)
        self.assertEqual(row['p'], 1.)
        self.assertEqual(a.bh([row['p']]).tolist(), [1.])


if __name__ == '__main__':
    unittest.main(verbosity=2)
