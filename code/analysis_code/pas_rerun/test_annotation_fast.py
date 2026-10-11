import unittest
import numpy as np

from polyapipe_annotation_fast import CandidateIndex


class IntervalLookupIntentTests(unittest.TestCase):
    def item(self, starts, ends):
        starts=np.array(starts);ends=np.array(ends)
        return starts,ends,np.maximum.accumulate(ends),[{'row':i} for i in range(len(starts))]

    def test_long_enclosing_interval_is_not_lost(self):
        item=self.item([0,100,200,300],[10000,110,210,310])
        self.assertEqual(list(CandidateIndex.hits(item,500,551)),[{'row':0}])

    def test_half_open_boundaries_and_original_order(self):
        item=self.item([10,20,20,30],[20,21,40,41])
        self.assertEqual(list(CandidateIndex.hits(item,20,30)),[{'row':1},{'row':2}])
        self.assertEqual(list(CandidateIndex.hits(item,41,42)),[])

    def test_nested_overlaps_match_independent_exhaustive_search(self):
        rng=np.random.default_rng(20261003)
        starts=np.sort(rng.integers(0,10000,2000));ends=starts+rng.integers(1,3000,2000)
        item=self.item(starts,ends)
        for start in rng.integers(-100,13000,500):
            end=start+51
            expected=[{'row':i} for i,(a,b) in enumerate(zip(starts,ends)) if a<end and b>start]
            self.assertEqual(list(CandidateIndex.hits(item,start,end)),expected)


if __name__=='__main__':
    unittest.main()
