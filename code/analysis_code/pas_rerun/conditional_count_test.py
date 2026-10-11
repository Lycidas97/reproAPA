"""Conditional Fisher tests with exact probability ordering near floating ties."""
from functools import lru_cache
from math import factorial,prod
import numpy as np
from scipy.special import gammaln
from scipy.stats import fisher_exact,random_table

B=131071


def probability_tail_mask(tables,observed):
    """For fixed margins probability is proportional to 1 / product(cell!)."""
    observed=np.asarray(observed,dtype=np.int64);tables=np.asarray(tables,dtype=np.int64)
    assert tables.ndim==3 and observed.shape==tables.shape[1:]
    observed_score=gammaln(observed+1).sum();scores=gammaln(tables+1).sum(axis=(1,2))
    # Bound roundoff on the factorial terms, not on the much smaller log mass
    # obtained by subtracting large constants. Resolve the boundary exactly.
    tolerance=100*np.finfo(float).eps*max(1.,float(gammaln(observed.sum()+1)))
    result=scores>observed_score+tolerance
    near=np.flatnonzero(abs(scores-observed_score)<=tolerance)
    if len(near):
        @lru_cache(None)
        def denominator(values):return prod(factorial(int(v)) for v in values)
        base=denominator(tuple(sorted(observed.ravel())))
        for i in near:result[i]=denominator(tuple(sorted(tables[i].ravel())))>=base
    return result


def fisher(job):
    index,a,b,random_seed=job
    table=np.stack((a,b)).astype(np.int64)
    assert np.all(table>=0) and np.all(table.sum(1)>0)
    table=table[:,table.sum(0)>0];assert table.shape[1]>=2
    if table.shape[1]==2:
        return index,float(fisher_exact(table,alternative='two-sided').pvalue),'Fisher_exact_2x2',0,-1,random_seed,2
    distribution=random_table(table.sum(1),table.sum(0),seed=np.random.default_rng(random_seed));extreme=0
    for start in range(0,B,10000):
        simulated=distribution.rvs(size=min(10000,B-start))
        extreme+=int(probability_tail_mask(simulated,table).sum())
    return index,(extreme+1)/(B+1),'Fisher_MC_2xn',B,extreme,random_seed,table.shape[1]
