"""Reserve concurrent lanes for base, ASNR and ASNR plus terminal protection.

An adopted current-run route is observed to completion, never recalculated.
New routes use the same scientific modules with independent study/route outputs.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import time
import re
from pipeline import read,write,digest,require,execute
from prepare_pipeline import prepare as bind
from result_records import ResultBundle,collect_sources,source_identity
from parallel_processing import derive_route_parallel

FAMILIES=('uniform_count_base','ASNR','ASNR_plus_terminal')


def process_identity(pid):
    p=Path('/proc')/str(pid)
    if not p.exists():return None
    try:
        st=(p/'stat').read_text().rsplit(')',1)[1].split()
        if st[0]=='Z':return None
        return dict(pid=pid,start_ticks=st[19],argv=(p/'cmdline').read_bytes().decode().rstrip('\0').split('\0'))
    except FileNotFoundError:return None


def family_plan(intake,adoptions):
    rows=intake['analysis_routes'];require(len({r['route'] for r in rows})==len(rows),'Duplicate route')
    require({r['family'] for r in rows}==set(FAMILIES),'Unexpected processing family')
    require(set(adoptions)<={r['route'] for r in rows},'Unknown adopted route')
    require(all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',r['route']) for r in rows),'Unsafe route name')
    groups={f:[] for f in FAMILIES}
    for row in rows:groups[row['family']].append(dict(row))
    for f in groups:groups[f].sort(key=lambda r:r['route'] not in adoptions)
    return groups


def wait_for_adopted(route,adoption,controller,stage,poll_seconds=15):
    config=Path(adoption['config']);bundle=Path(adoption['bundle'])
    def verify_done():
        cfg=read(config);identity=source_identity(cfg)
        require(cfg['analysis_route']==route,'Adoption route mismatch')
        require(identity==source_identity(adoption['identity']) if 'identity' in adoption else True,'Adoption identity changed')
        result=ResultBundle(bundle,identity)
        state=read(Path(cfg['analysis_root'])/'pipeline-status.json')
        require(state['phase']=='statistics_complete','Bundle exists before whole DAG acceptance')
        return dict(source_id=route,bundle=str(bundle),identity=result.identity,config=str(config),adopted=True)
    while True:
        if config.is_file():
            cfg=read(config);require(cfg['analysis_route']==route,'Adoption config belongs to another route')
            p=Path(cfg['analysis_root'])/'pipeline-status.json'
            if p.exists():
                s=read(p);require(s['phase']!='failed','Adopted route failed: '+str(s.get('error')))
                if s['phase']=='statistics_complete' and (bundle/'bundle-complete.json').is_file():return verify_done()
                write(stage,dict(source_id=route,phase='adopted_statistics_running',
                      completed_nodes=len(s.get('completed',[])),bundle=str(bundle)))
        prep=Path(adoption['source_stage'])
        if prep.exists():require(read(prep)['phase']!='preparation_failed','Adopted input preparation failed')
        actual=process_identity(controller['pid'])
        if actual!=controller:
            if config.is_file() and (bundle/'bundle-complete.json').is_file():return verify_done()
            raise ValueError('Adopted controller exited or changed before route completion: '+route)
        time.sleep(poll_seconds)


def run_group(rows,capacity,work,status_path):
    states={r['route']:dict(route=r['route'],status='pending') for r in rows}
    results=[]
    def status(phase):write(status_path,dict(phase=phase,capacity=capacity,routes=list(states.values())))
    pending=list(rows);active={};status('running')
    from concurrent.futures import wait,FIRST_COMPLETED
    with ThreadPoolExecutor(max_workers=capacity) as pool:
        while pending or active:
            while pending and len(active)<capacity:
                row=pending.pop(0);states[row['route']]['status']='running';active[pool.submit(work,row)]=row['route']
            status('running');done,_=wait(active,return_when=FIRST_COMPLETED)
            for future in done:
                route=active.pop(future)
                try:
                    result=future.result();results.append(result);states[route].update(status='statistics_complete',bundle=result['bundle'])
                except Exception as error:states[route].update(status='failed',error=repr(error))
            status('running')
    failures=[r for r in states.values() if r['status']=='failed'];status('partial_failure' if failures else 'statistics_complete')
    return dict(results=results,failures=failures)


def run(a):
    out=Path(a.out).resolve();base=Path(a.prepared_base).resolve()
    require(not out.exists() and not out.is_relative_to(base),'Use a fresh family-series directory')
    require(not any(p in out.parts for p in ['Fig','main','supplement','production-v15']),'Never write into production')
    require(read(base/'intake-complete.json')['status']=='READY_FOR_CORE','Full base intake required')
    adoption=read(a.adoption) if a.adoption else dict(routes={},controller=None)
    groups=family_plan(read(a.intake_config),adoption['routes']);out.mkdir(parents=True)
    write(out/'family-plan.json',dict(groups=groups,family_workers=a.family_workers,core_workers=a.core_workers,
          preparation_workers=a.preparation_workers,node_workers=a.workers,adoption=adoption,
          frozen_modules=True,production_replaced=False,new_BAM_recount=False))
    def work(row):
        route=row['route'];stage=out/'source-stages'/(route+'.json')
        def state(phase,**extra):write(stage,dict(source_id=route,family=row['family'],phase=phase,**extra))
        try:
            if route in adoption['routes']:
                state('adopted_running',bundle=adoption['routes'][route]['bundle'])
                return wait_for_adopted(route,adoption['routes'][route],adoption['controller'],stage)
            native=base if route=='base' else out/'native'/route
            if route!='base':
                state('deriving_native',study_workers=a.preparation_workers)
                derive_route_parallel(base,route,native,a.preparation_workers)
            state('binding_pipeline');folder=out/'runs'/route
            cfg=bind(a.environment,native,folder,core_workers=a.core_workers)
            state('statistics_running',config=str(folder/'pipeline.config.json'),core_workers=a.core_workers,node_workers=a.workers)
            execute(folder/'pipeline.config.json',a.workers)
            result=dict(source_id=route,bundle=str(folder/'statistics/bundle'),identity=source_identity(cfg),adopted=False)
            state('statistics_complete',bundle=result['bundle']);return result
        except Exception as error:state('failed',error=repr(error));raise
    write(out/'family-status.json',dict(phase='running',families=list(groups),routes=sum(map(len,groups.values())),production_replaced=False))
    results=[];failures=[]
    with ThreadPoolExecutor(max_workers=len(groups)) as pool:
        futures={pool.submit(run_group,rows,min(a.family_workers,len(rows)),work,
                            out/'family-status'/(family+'.json')):family for family,rows in groups.items()}
        for future in as_completed(futures):
            r=future.result();results+=r['results'];failures+=r['failures']
    if results:
        try:collect_sources(results,out/('completed-result-series.json' if failures else 'result-series.json'))
        except Exception as error:
            write(out/'family-status.json',dict(phase='result_index_failed',error=repr(error)));raise
    write(out/'family-status.json',dict(phase='partial_failure' if failures else 'statistics_complete',
          complete=len(results),failures=failures,production_replaced=False))
    require(not failures,'Some processing routes failed; independent results were preserved')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['environment','intake-config','prepared-base','out']:p.add_argument('--'+key,required=True)
    p.add_argument('--adoption');p.add_argument('--family-workers',type=int,default=2)
    p.add_argument('--core-workers',type=int,default=7);p.add_argument('--preparation-workers',type=int,default=4)
    p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    require(min(a.family_workers,a.core_workers,a.preparation_workers,a.workers)>0,'All worker counts must be positive')
    run(a)
