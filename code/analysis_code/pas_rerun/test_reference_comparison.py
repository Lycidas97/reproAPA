import unittest

import numpy as np
import pandas as pd

from compare_reference_runs import effect_comparison, member_comparison


class ReferenceComparisonIntent(unittest.TestCase):
    def row(self,gene,setting='supported50_p08',effect=.03):
        return dict(species='human',study_id='study',route='route',gene_id=gene,gene_symbol=gene,
            contrast='P_S',design='true_source',setting=setting,delta_DWUI=effect,delta_RWUI=effect,
            delta_PDUI=effect,delta_PPUI=-effect,MPRO=abs(effect))

    def test_missing_in_control_stays_missing_and_unsigned_metric_has_no_direction(self):
        a=pd.DataFrame([self.row('a'),self.row('b')]);b=pd.DataFrame([self.row('a',effect=-.03)])
        joined,summary=effect_comparison(a,b)
        self.assertTrue(np.isnan(joined.set_index('gene_id').loc['b','delta_DWUI_v3']))
        self.assertTrue(summary.unmeasured_in_v3.eq(1).all())
        self.assertTrue(summary[summary.metric.eq('MPRO')].direction_agreement.isna().all())
        self.assertTrue(summary[summary.metric.eq('delta_DWUI')].direction_agreement.eq(0).all())

    def test_expanded_primary_change_is_a_pipeline_error_not_reference_sensitivity(self):
        a=pd.DataFrame([self.row('a','expanded_primary_p08')])
        b=pd.DataFrame([self.row('a','expanded_primary_p08',.04)])
        with self.assertRaises(AssertionError):effect_comparison(a,b)
        effect_comparison(a,a.copy())

    def test_membership_keeps_metric_and_threshold_denominators_separate(self):
        a=pd.DataFrame([dict(species='mouse',metric='DWUI',threshold=.02,gene_id='a'),
                        dict(species='mouse',metric='MPRO',threshold=.1,gene_id='a')])
        b=pd.DataFrame([dict(species='mouse',metric='DWUI',threshold=.05,gene_id='a')])
        members,summary=member_comparison(a,b,['species','metric','threshold','gene_id'],['species','metric','threshold'])
        self.assertEqual(len(members),3)
        self.assertEqual(int(summary.shared_members.sum()),0)

    def test_empty_planned_sets_remain_zero_without_inventing_perfect_agreement(self):
        empty=pd.DataFrame(columns=['species','gene_id'])
        planned=pd.DataFrame({'species':['human','mouse']})
        members,summary=member_comparison(empty,empty,['species','gene_id'],['species'],planned)
        self.assertEqual(len(members),0)
        self.assertEqual(len(summary),2)
        self.assertTrue(summary['union'].eq(0).all())
        self.assertTrue(summary.jaccard.isna().all())


if __name__=='__main__':unittest.main()
