"""Small shared statistical and provenance helpers for the fixed protocol."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = {}

def read(path, **kwargs):
    path = Path(path)
    INPUTS[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, sep='\t', **kwargs)

def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.15g',
        compression={'method':'gzip', 'mtime':0} if name.endswith('.gz') else None)

def manifest(name, **kwargs):
    (HERE / name).write_text(json.dumps(dict(input_sha256=INPUTS, **kwargs), indent=2)+'\n')

def rng(*keys):
    return np.random.default_rng(20261002 + int(hashlib.sha256('|'.join(map(str,keys)).encode()).hexdigest()[:12],16))

def bh(p):
    p = np.asarray(p, float)
    assert np.isfinite(p).all() and ((p>=0)&(p<=1)).all()
    order = np.argsort(p); out = np.empty(len(p))
    out[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return out

def probability(null, observed):
    return (1 + len(null)-np.searchsorted(np.sort(null), observed, side='left'))/(len(null)+1)

def rate(rejected):
    k, n = int(np.sum(rejected)), len(rejected)
    return dict(rejected=k, simulations=n, rejection_rate=k/n,
        lo=0. if k==0 else beta.ppf(.025,k,n-k+1),
        hi=1. if k==n else beta.ppf(.975,k+1,n-k))

def groups(frame, columns):
    if not columns:
        return [np.arange(len(frame))]
    return list(frame.groupby(columns, dropna=False, sort=True).indices.values())

def binary_draw(labels, cells, n, random):
    out = np.broadcast_to(labels,(n,len(labels))).copy()
    for ii in cells:
        if 0 < labels[ii].sum() < len(ii):
            out[:,ii] = random.permuted(np.broadcast_to(labels[ii],(n,len(ii))), axis=1)
    return out

def mobility(labels, cells):
    probabilities = np.zeros(len(labels)); mixed = 0; changed = 0.; singletons = 0
    for ii in cells:
        p = labels[ii].mean(); probabilities[ii] = p
        mixed += len(ii)*(0<p<1); changed += len(ii)*2*p*(1-p)
        singletons += len(ii)==1
    return probabilities, dict(strata=len(cells), singleton_strata=singletons,
        movable_fraction=mixed/len(labels), expected_changed_fraction=changed/len(labels))
