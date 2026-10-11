"""Attach supplementary completion to each finished source in an active series."""
import argparse
import fcntl
from pathlib import Path
import subprocess
import sys
import time
import re

from pipeline import digest,read,write,require
from result_records import ResultBundle,collect_sources


def source_rows(series):
    """Resolve family adoptions from their plan, never from stale pending rows."""
    plan=series/'family-plan.json'
    if not plan.exists():
        rows=read(series/'source-runs.json')['routes']
        # The common-catalogue controller owns preparation and the main DAG;
        # it publishes these terminal names without starting supplementary work.
        terminal={'MAIN_STATISTICS_COMPLETE':'statistics_complete','FAILED':'failed'}
        return [dict(r,status=terminal.get(r['status'],r['status'])) for r in rows]
    plan=read(plan);groups=plan['groups'];adopted=plan['adoption']['routes'];rows=[]
    require(set(adopted)<={r['route'] for group in groups.values() for r in group},'Unknown adopted source')
    for family,group in groups.items():
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',family),'Invalid family name')
        require(all(r['family']==family for r in group),'Family plan mismatch')
        current=read(series/'family-status'/f'{family}.json')['routes']
        states={r['route']:r for r in current}
        require(len(states)==len(current)==len(group) and set(states)=={r['route'] for r in group},'Family route membership changed')
        for item in group:
            key=item['route'];status=states[key];folder=series/'runs'/key
            binding=adopted.get(key,dict(config=str(folder/'pipeline.config.json'),bundle=str(folder/'statistics/bundle')))
            require('bundle' not in status or status['bundle']==binding['bundle'],'Family bundle binding changed')
            rows.append(dict(source_id=key,config=binding['config'],bundle=binding['bundle'],
                family=family,adopted=key in adopted,status=status['status'],error=status.get('error')))
    return rows


def source_bindings(rows):
    ids=[r['source_id'] for r in rows]
    require(bool(ids) and len(set(ids))==len(ids) and all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',k) for k in ids),'Invalid or duplicate source IDs')
    for field in ['config','bundle']:
        require(all(Path(r[field]).is_absolute() for r in rows),'Source paths must be absolute')
        require(len({str(Path(r[field]).resolve()) for r in rows})==len(rows),'Sources share '+field)
    return {r['source_id']:{k:r[k] for k in ['config','bundle','family','adopted'] if k in r} for r in rows}


def run(series,contracts,out,poll_seconds=60,once=False):
    series=Path(series).resolve();contracts=Path(contracts).resolve();out=Path(out).resolve()
    require(not out.is_relative_to(series),'Use a separate continuation directory')
    out.mkdir(parents=True,exist_ok=True)
    with (out/'continuation.lock').open('a') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        bindings=source_bindings(source_rows(series));expected=set(bindings)
        state=read(out/'status.json') if (out/'status.json').exists() else dict(status='waiting_for_source_bundles',routes={})
        hashes={p.name:digest(p) for p in contracts.glob('supplement-*.tsv')}
        if 'contracts_sha256' in state:
            require(state['contracts_sha256']==hashes and state['source_series']==str(series),'Continuation source or contracts changed')
        if 'source_bindings' in state:require(state['source_bindings']==bindings,'Continuation source bindings changed')
        state.update(source_series=str(series),contracts_sha256=hashes,
            source_bindings=bindings,source_count=len(expected),production_replaced=False,source_jobs_interrupted=False)
        while True:
            rows=source_rows(series);require(source_bindings(rows)==bindings,'Source series bindings changed')
            for row in rows:
                key=row['source_id'];previous=state['routes'].get(key,{})
                if previous.get('status') in ['complete','failed','upstream_failed']:continue
                state['routes'][key]=dict(source_bundle=row['bundle'],source_config=row['config'],upstream_status=row['status'],
                    **{k:row[k] for k in ['family','adopted'] if k in row})
                if row['status']=='failed':state['routes'][key].update(status='upstream_failed',error=row.get('error','Source pipeline failed'));continue
                bundle=Path(row['bundle']);config=Path(row['config'])
                # A bundle can be written before final input/DAG acceptance.
                if row['status']!='statistics_complete' or not (bundle/'bundle-complete.json').exists():
                    state['routes'][key].update(status='waiting_for_source_bundle');continue
                cfg=read(config)
                require(read(Path(cfg['analysis_root'])/'pipeline-status.json')['phase']=='statistics_complete','Source bundle precedes final DAG acceptance')
                count_audit=bundle/'tables/count_inference/audit.json'
                if count_audit.exists() and read(count_audit)['checks'].get('robust_probability_tie_comparison') is not True:
                    state['routes'][key].update(status='waiting_for_corrected_count_inference')
                    continue
                if 'family' in row:require(cfg['analysis_route']==key,'Family config belongs to another route')
                destination=out/'runs'/key
                if destination.exists():
                    state['routes'][key].update(status='failed',error='Unfinished previous continuation directory; inspect rather than overwrite');continue
                state['routes'][key].update(status='running',output=str(destination))
                state['status']='running';write(out/'status.json',state)
                log=out/'logs'/f'{key}.log';log.parent.mkdir(exist_ok=True)
                command=[sys.executable,str(Path(__file__).with_name('complete_supplement_bundle.py')),
                    '--bundle',str(bundle),'--config',str(config),'--contracts',str(contracts),'--out',str(destination)]
                with log.open('w') as handle:result=subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT)
                if result.returncode:
                    state['routes'][key].update(status='failed',exit_code=result.returncode,log=str(log))
                else:
                    ready=destination/'bundle';reader=ResultBundle(ready)
                    state['routes'][key].update(status='complete',bundle=str(ready),source_key=reader.identity['source_key'],log=str(log))
                write(out/'status.json',state)
            terminal={'complete','failed','upstream_failed'}
            done=all(state['routes'].get(k,{}).get('status') in terminal for k in expected)
            state['status']=('complete' if all(state['routes'][k]['status']=='complete' for k in expected) else 'complete_with_failures') if done else 'waiting_for_source_bundles'
            state['checked_at_unix']=time.time();write(out/'status.json',state)
            if done:
                completed=[dict(source_id=k,bundle=r['bundle']) for k,r in state['routes'].items() if r['status']=='complete']
                target=out/('result-series.json' if state['status']=='complete' else 'completed-result-series.json')
                if completed and not target.exists():collect_sources(completed,target)
                return 0 if state['status']=='complete' else 1
            if once:return 0
            time.sleep(poll_seconds)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['series','contracts','out']:p.add_argument('--'+name,required=True)
    p.add_argument('--poll-seconds',type=int,default=60);p.add_argument('--once',action='store_true')
    a=p.parse_args();require(a.poll_seconds>=10,'Poll interval must be at least ten seconds')
    sys.exit(run(a.series,a.contracts,a.out,a.poll_seconds,a.once))
