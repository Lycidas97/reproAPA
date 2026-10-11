"""Independent statistical intent fixtures, without production imports."""
import unittest

import numpy as np
import pandas as pd


def study_event(effects, threshold):
    """A study votes only when every observed series exceeds the same bound."""
    values=np.asarray(effects,dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError('Only actual finite, evaluable series may be supplied')
    return -1 if values.max() <= -threshold else 1 if values.min() >= threshold else 0


def recurrent_class(study_events):
    """Missing studies do not vote; independent studies, not series, recur."""
    available=[v for v in study_events if pd.notna(v)]
    if any(v not in [-1,0,1] for v in available):
        raise ValueError('Invalid study vote')
    negative=available.count(-1)>=2
    positive=available.count(1)>=2
    if negative and positive:raise ValueError('Both directions recur: design exceeds reviewed three-study rule')
    return 'proximal' if negative else 'distal' if positive else 'nonrecurrent'


def independent_bh(values):
    p=np.asarray(values,dtype=float)
    if not np.isfinite(p).all() or ((p<0)|(p>1)).any():
        raise ValueError('Formal family requires explicit finite p, including NE=1')
    order=np.argsort(p,kind='stable')
    ranked=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty(len(p));out[order]=np.minimum(1,ranked)
    return out


class StatisticalIntent(unittest.TestCase):
    def test_inclusive_threshold_boundary(self):
        self.assertEqual(study_event([-.02,-.03],.02),-1)
        self.assertEqual(study_event([.02,.03],.02),1)

    def test_same_sign_but_incomplete_threshold_agreement_is_not_vote(self):
        self.assertEqual(study_event([-.01,-.05],.02),0)
        self.assertEqual(study_event([.01,.05],.02),0)

    def test_average_does_not_replace_series_agreement(self):
        self.assertLess(np.mean([-.06,.01]),-.02)
        self.assertEqual(study_event([-.06,.01],.02),0)

    def test_third_study_opposition_does_not_cancel_two_study_recurrence(self):
        self.assertEqual(recurrent_class([-1,-1,1]),'proximal')
        self.assertEqual(recurrent_class([1,1,-1]),'distal')

    def test_missing_study_does_not_create_vote(self):
        self.assertEqual(recurrent_class([-1,np.nan]),'nonrecurrent')
        with self.assertRaises(ValueError):study_event([-.03,np.nan],.02)

    def test_repeated_series_are_one_study(self):
        event=study_event([-.03,-.04,-.05],.02)
        self.assertEqual(recurrent_class([event]),'nonrecurrent')

    def test_declared_ne_item_stays_in_family(self):
        without=independent_bh([.01,.1])
        with_ne=independent_bh([.01,.1,1])
        self.assertAlmostEqual(without[0],.02)
        self.assertAlmostEqual(with_ne[0],.03)
        self.assertEqual(with_ne[2],1)

    def test_nonrec_is_not_missing_qualification(self):
        records=pd.DataFrame({'id':['a','b'],'eligible':[True,False],
                              'classification':['nonrecurrent',None]})
        self.assertEqual(records.loc[records.eligible,'classification'].tolist(),['nonrecurrent'])
        self.assertTrue(records.loc[~records.eligible,'classification'].isna().all())


if __name__=='__main__':unittest.main(verbosity=2)
