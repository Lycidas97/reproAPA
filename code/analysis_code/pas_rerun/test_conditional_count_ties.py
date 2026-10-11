from math import comb
import unittest
import numpy as np
from scipy.stats import random_table
from conditional_count_test import probability_tail_mask,fisher

class ConditionalTieTests(unittest.TestCase):
    def test_real_two_by_three_enumeration_and_batching_preserve_equal_masses(self):
        observed=np.array([[5,114,76],[2,37,24]])
        margins=observed.sum(0);n=observed[0].sum();tables=[];weights=[]
        for x in range(max(0,n-margins[1]-margins[2]),min(margins[0],n)+1):
            for y in range(max(0,n-x-margins[2]),min(margins[1],n-x)+1):
                v=np.array([x,y,n-x-y]);tables.append([v,margins-v]);weights.append(comb(int(margins[0]),x)*comb(int(margins[1]),y)*comb(int(margins[2]),int(v[2])))
        tables=np.array(tables);obs=comb(7,5)*comb(151,114)*comb(100,76)
        expected=np.array([w<=obs for w in weights]);actual=probability_tail_mask(tables,observed)
        np.testing.assert_array_equal(actual,expected)
        exact=sum(w for w,keep in zip(weights,actual) if keep)/comb(258,195)
        self.assertAlmostEqual(exact,.9622042499127195,places=14)
        # A tie's classification must not depend on the shape or order of a vectorized batch.
        for q in [observed,np.array([[6,113,76],[1,38,24]])]:
            self.assertTrue(probability_tail_mask(q[None],observed)[0])
            self.assertTrue(probability_tail_mask(np.tile(q,(10000,1,1)),observed).all())
    def test_current_MC_estimate_contains_the_exact_enumerated_probability(self):
        r=fisher((0,np.array([5,114,76]),np.array([2,37,24]),1234))
        self.assertLess(abs(r[1]-.9622042499127195),.004)
    def test_small_other_margins_match_integer_probability_order(self):
        rng=np.random.default_rng(31)
        for observed in [np.array([[1,9,3],[5,2,7]]),np.array([[0,2,5,1],[4,1,0,8]])]:
            samples=random_table.rvs(observed.sum(1),observed.sum(0),size=1000,random_state=rng)
            from math import factorial,prod
            denominator=lambda x:prod(factorial(int(v)) for v in x.ravel())
            expected=np.array([denominator(q)>=denominator(observed) for q in samples])
            np.testing.assert_array_equal(probability_tail_mask(samples,observed),expected)

if __name__=='__main__':unittest.main()
