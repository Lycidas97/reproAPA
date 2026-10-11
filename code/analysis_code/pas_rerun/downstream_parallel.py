"""Ordered, process-isolated tasks with complete input provenance."""
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
import os
import time

from downstream_runtime import Context


def _task(function, config, node, out, item):
    ctx=Context(config,node,out)
    start=time.monotonic()
    value=function(ctx,item)
    return value,ctx.inputs,dict(task=item,pid=os.getpid(),seconds=time.monotonic()-start)


def ordered_map(ctx,function,items):
    items=list(items)
    workers=ctx.config.get('node_workers',{}).get(ctx.node,1)
    assert isinstance(workers,int) and workers>=1
    if not items:return []
    if workers==1:
        return [function(ctx,item) for item in items]
    # Each child reads immutable inputs and writes only its task-owned paths.
    # Aggregation, joint fits and multiple-testing correction stay in the caller.
    values=[];records=[]
    with ProcessPoolExecutor(max_workers=min(workers,len(items)),
                             mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(_task,function,str(ctx.config_path),ctx.node,str(ctx.out),item) for item in items]
        for future in futures:
            value,inputs,record=future.result()
            for path,sha in inputs.items():
                assert path not in ctx.inputs or ctx.inputs[path]==sha, f'Conflicting worker input: {path}'
                ctx.inputs[path]=sha
            values.append(value);records.append(record)
    (ctx.out/'parallel-execution.json').write_text(json.dumps(dict(workers=workers,tasks=records),indent=2)+'\n')
    return values
