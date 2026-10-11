"""One-command source switch: bind prepared input(s), compute statistics, export tables.

Rendering is intentionally a separate command. No polling or future scheduling.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from pathlib import Path
import re
from apadetect_prepare import prepare,derive_route
from prepare_pipeline import prepare as bind
from pipeline import execute,read,require,write,digest
from result_records import collect_sources,source_identity


def execute_sources(runs,out,workers,source_workers=1,prepare_source=None,prepare_only=False):
    """Independent runs may overlap, but never share output directories or identities."""
    require(source_workers>0 and workers>0,'Worker counts must be positive')
    require(len({r['source_id'] for r in runs})==len(runs),'Duplicate source IDs')
    require(len({str(Path(r['config']).resolve()) for r in runs})==len(runs),'Sources share a pipeline config')
    states={r['source_id']:dict(r,status='pending') for r in runs};active={};pending=list(runs)
    def one(row):
        bound=prepare_source(dict(row)) if prepare_source else row
        if not prepare_only:execute(bound['config'],workers)
        return bound
    def status(phase,**extra):
        write(Path(out)/'source-runs.json',dict(status=phase,routes=list(states.values()),
            source_workers=source_workers,node_workers_per_source=workers,production_replaced=False,figures_rendered=False,**extra))
    status('running')
    with ThreadPoolExecutor(max_workers=source_workers) as pool:
        while pending or active:
            while pending and len(active)<source_workers:
                row=pending.pop(0);states[row['source_id']]['status']='running'
                active[pool.submit(one,row)]=row['source_id']
            status('running');finished,_=wait(active,return_when=FIRST_COMPLETED)
            for future in finished:
                key=active.pop(future)
                try:states[key].update(future.result(),status='prepared' if prepare_only else 'statistics_complete')
                except Exception as error:states[key].update(status='failed',error=repr(error))
            status('running')
    failed=[k for k,r in states.items() if r['status']=='failed']
    status('partial_failure' if failed else 'prepared' if prepare_only else 'statistics_complete')
    # Each successful route retains a complete bundle even if another fails.
    complete=[dict(source_id=r['source_id'],bundle=r['bundle'],identity=r['identity'])
              for r in states.values() if r['status']=='statistics_complete']
    if complete:
        try:collect_sources(complete,Path(out)/('completed-result-series.json' if failed else 'result-series.json'))
        except Exception as error:
            status('result_index_failed',error=repr(error));raise
    require(not failed,'Failed sources (other sources were allowed to finish): '+', '.join(failed))


def run(a):
    out=Path(a.out).resolve();require(not out.exists(),'Use a new run-series directory');out.mkdir(parents=True)
    if getattr(a,'inputs',None):
        require(not a.prepared_native and not a.prepared_base and not a.intake_config and a.routes=='base',
                '--inputs cannot be combined with a single native input or derived routes')
        items=read(a.inputs)['inputs'];require(bool(items),'Empty input-source list');sources={}
        for item in items:
            key=item['source_id'];require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',key) is not None,'Unsafe source_id')
            require(key not in sources,'Duplicate source_id');sources[key]=Path(item['prepared_native']).resolve()
    elif a.prepared_native:
        require(not a.intake_config and not a.prepared_base and a.routes=='base','Prepared inputs already define one route')
        sources={'prepared':Path(a.prepared_native).resolve()}
    else:
        require(a.intake_config,'Provide --prepared-native or --intake-config')
        intake=read(a.intake_config);known={r['route'] for r in intake['analysis_routes']}
        routes=sorted(known,key=lambda r:(r!='base',r)) if a.routes=='all' else list(dict.fromkeys(a.routes.split(',')))
        require(set(routes)<=known,'Unknown analysis route')
        if a.prepared_base:
            base=Path(a.prepared_base).resolve()
            receipt=read(base/'intake-complete.json');manifest=read(base/'intake.native.json')
            require(receipt['status']=='READY_FOR_CORE' and manifest['analysis_route']=='base','A complete base intake is required')
            require(manifest['run_id']==intake['run_id']+'__base','Prepared base belongs to a different intake')
        else:
            base=out/'native/base';prepare(a.intake_config,base)
        sources={}
        for route in routes:
            p=base if route=='base' else out/'native'/route
            sources[route]=p
    runs=[]
    for route,source in sources.items():
        path=out/'runs'/route
        runs.append(dict(route=route,source_id=route,config=str(path/'pipeline.config.json'),
                         native_source=str(source),bundle=str(path/'statistics/bundle')))
    def prepare_source(row):
        source=Path(row['native_source']);stage=out/'source-stages'/(row['source_id']+'.json')
        def state(phase,**extra):write(stage,dict(source_id=row['source_id'],phase=phase,**extra))
        try:
            if not getattr(a,'inputs',None) and not a.prepared_native and row['route']!='base':
                state('deriving_native');derive_route(base,row['route'],source)
            state('binding_pipeline')
            cfg=bind(a.environment,source,Path(row['config']).parent,python=a.python)
            row.update(run_id=cfg['run_id'],identity=source_identity(cfg))
            state('prepared' if a.prepare_only else 'statistics_dispatched',config=row['config'])
            return row
        except Exception as error:state('preparation_failed',error=repr(error));raise
    execute_sources(runs,out,a.workers,getattr(a,'source_workers',1),prepare_source,a.prepare_only)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--environment',required=True);p.add_argument('--out',required=True)
    p.add_argument('--prepared-native');p.add_argument('--prepared-base');p.add_argument('--intake-config');p.add_argument('--routes',default='base')
    p.add_argument('--inputs',help='JSON inputs: [{source_id, prepared_native}], including different PAS callers')
    p.add_argument('--source-workers',type=int,default=1,help='Concurrent input-source DAGs; default 1 bounds memory')
    p.add_argument('--prepare-only',action='store_true');p.add_argument('--workers',type=int,default=2);p.add_argument('--python')
    a=p.parse_args();require(a.workers>0 and a.source_workers>0,'workers must be positive');run(a)
