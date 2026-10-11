"""Verify metadata-index caching against the frozen pre-cache implementation."""
import argparse
import importlib.util
import json
from pathlib import Path
import time

import numpy as np

from downstream_runtime import Context,Study


def compare(left,right):
    assert set(left)==set(right)
    for key in left:
        if isinstance(left[key],dict):compare(left[key],right[key])
        else:np.testing.assert_array_equal(left[key],right[key],err_msg=str(key))


def main(args):
    spec=importlib.util.spec_from_file_location('frozen_runtime',args.original)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    ctx=Context(args.config,'membership_cache_validation',args.out)
    before=old.Study(ctx,args.study);after=Study.__new__(Study);after.__dict__.update(before.__dict__)
    times={'original':0.,'cached':0.};cases=0
    for gene in list(before.genes)[:args.genes]:
        for setting in ['supported50_p08','supported50_p09','expanded_primary_p08']:
            axis=before.axis(gene,setting)
            if len(axis)<2:continue
            data=before.gene_counts(axis)
            for depth in [20,40,80]:
                start=time.monotonic();a=before.profile_states(data,setting,depth);times['original']+=time.monotonic()-start
                start=time.monotonic();b=after.profile_states(data,setting,depth);times['cached']+=time.monotonic()-start
                compare(a,b);cases+=1
    report=dict(study=args.study,genes=args.genes,profile_cases=cases,exact_array_equality=True,
        seconds=times,speedup=times['original']/times['cached'])
    (ctx.out/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--original',required=True)
    p.add_argument('--out',required=True);p.add_argument('--study',default='GSE103221');p.add_argument('--genes',type=int,default=20)
    main(p.parse_args())
