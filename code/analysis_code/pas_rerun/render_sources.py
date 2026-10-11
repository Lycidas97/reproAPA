"""Dispatch the same registered figure builders over one or many result sources.

This only routes immutable result records. A legacy/unadapted builder is a hard
preflight failure, never an excuse to draw an old result or publish a partial set.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import subprocess
import sys
from pipeline import digest,read,write,require
from result_records import ResultBundle,contained


def plan(series,project,out,source_ids=None,figures=None,python=None):
    series=read(series);require(series['schema']=='pas-result-series/1','Unknown result-series format')
    project=Path(project).resolve();out=Path(out).resolve()
    require(not out.exists(),'Use a new candidate render directory')
    require(not any(p in out.parts for p in ['Fig','main','supplement','production-v15']),'Never render into production')
    sources=series['sources'];known={s['source_id'] for s in sources}
    require(len(known)==len(sources),'Duplicate source IDs in series')
    if source_ids is not None:require(set(source_ids)<=known,'Unknown requested source');sources=[s for s in sources if s['source_id'] in source_ids]
    require(bool(sources),'No sources selected')
    jobs=[];schemas=set();keys=set();registry_ids=set()
    for source in sources:
        name=source['source_id'];folder=contained(out,name)
        require(folder.parent==out,'Source ID must be a single directory name')
        result=ResultBundle(source['bundle'],source['identity']);key=result.identity['source_key']
        require(key not in keys,'Same source listed twice');keys.add(key)
        require(digest(result.root/'result-records.json')==source['result_records_sha256'],'Series points to a changed result index')
        require(result.record['statistics_schema_id']==source['statistics_schema_id'],'Series statistical schema mismatch')
        schemas.add(source['statistics_schema_id'])
        require(not out.is_relative_to(result.root),'Do not write figures inside an immutable result bundle')
        registry=result.record['figures'];available={r['figure']:r for r in registry}
        require(len(available)==len(registry),'Duplicate figure registration')
        from result_records import content_id
        registry_ids.add(content_id(registry))
        requested=list(available) if figures is None else figures
        require(len(set(requested))==len(requested),'Duplicate requested figure')
        require(set(requested)<=set(available),'Unknown requested figure')
        for figure in requested:
            builder=available[figure];script=contained(project,builder['script'])
            require(script.is_file() and digest(script)==builder['sha256'],'Registered drawing script changed: '+figure)
            require(set(builder['panels'])=={r['panel'] for r in result.record['components'] if r['figure']==figure},
                    'Figure registration does not cover its panels: '+figure)
            # Verify all records before any renderer starts, including cross-source/table hashes.
            required_tables={}
            for c in result.record['components']:
                if c['figure']==figure:
                    result.table(c['panel'],c['component_id'])
                    required_tables[c['table']]=c['table_sha256']
            dest=contained(folder,figure);require(dest.parent==folder,'Figure ID must be a single directory name')
            ready=builder['interface']=='pas-panel-result-records/1' and builder['status']=='ready'
            jobs.append(dict(source_id=name,source=result.identity,figure=figure,panels=builder['panels'],
                statistics_schema_id=source['statistics_schema_id'],bundle=str(result.root),
                result_records_sha256=source['result_records_sha256'],out=str(dest),
                script_sha256=builder['sha256'],required_tables=required_tables,ready=ready,
                blocked_reason='' if ready else builder.get('blocked_reason','Builder has not been adapted to result records'),
                command=[python or sys.executable,str(script),'--records',str(result.root/'result-records.json'),'--out',str(dest)]))
    require(len(schemas)==1,'Selected sources have different statistical schemas')
    require(len(registry_ids)==1,'Selected sources have different figure registrations')
    require(bool(jobs),'No figures selected')
    return dict(schema='pas-render-plan/1',status='ready' if all(j['ready'] for j in jobs) else 'blocked_by_unadapted_builders',
                jobs=jobs,statistics_recomputed=False,production_replaced=False)


def execute(plan,out,workers):
    out=Path(out);require(workers>0,'workers must be positive')
    require(not out.exists(),'Use a fresh render directory');out.mkdir(parents=True)
    write(out/'render-plan.json',plan)
    require(plan['status']=='ready','Figure adapters are incomplete; see render-plan.json. No renderer was executed.')
    statuses=[]
    def job(row):
        dest=Path(row['out']);dest.mkdir(parents=True)
        try:
            require(digest(Path(row['command'][1]))==row['script_sha256'],'Drawing script changed after planning')
            current=ResultBundle(row['bundle'],row['source'])
            require(digest(current.root/'result-records.json')==row['result_records_sha256'],'Results changed after planning')
            with (dest/'command.log').open('w') as log:
                subprocess.run(row['command'],stdout=log,stderr=subprocess.STDOUT,check=True)
            receipt=read(dest/'render-receipt.json')
            require(receipt['source']==row['source'],'Cross-source render receipt')
            require(receipt['result_records_sha256']==row['result_records_sha256'],'Renderer used a different result index')
            require(receipt['figure']==row['figure'] and set(receipt['panels'])==set(row['panels']),'Incomplete figure render')
            require(receipt.get('tables')==row['required_tables'],'Renderer did not account for every registered result table')
            require(receipt.get('statistics_recomputed') is False and receipt.get('production_replaced') is False,'Invalid renderer scope')
            artifacts=receipt.get('artifacts',[]);require(bool(artifacts),'Renderer produced no artifacts')
            require(any(a['path'].endswith('.pdf') for a in artifacts),'Missing candidate PDF')
            require(any(a['path'].endswith('.png') for a in artifacts),'Missing candidate preview')
            for a in artifacts:require(digest(contained(dest,a['path']))==a['sha256'],'Changed render artifact')
            return dict(source_id=row['source_id'],figure=row['figure'],status='rendered',
                        receipt_sha256=digest(dest/'render-receipt.json'),publication_review='pending')
        except Exception as error:
            return dict(source_id=row['source_id'],figure=row['figure'],status='failed',error=repr(error))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(job,row) for row in plan['jobs']]
        for f in as_completed(futures):
            statuses.append(f.result());write(out/'render-status.json',dict(status='running',jobs=statuses,production_replaced=False))
    failed=[s for s in statuses if s['status']=='failed']
    write(out/'render-status.json',dict(status='failed' if failed else 'candidate_render_complete',jobs=statuses,
                                      production_replaced=False,publication_review='pending'))
    require(not failed,'Some figure builds failed; successful candidates were preserved')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',required=True);p.add_argument('--out',required=True)
    p.add_argument('--project',default=str(Path(__file__).resolve().parents[2]));p.add_argument('--source-ids',default='all')
    p.add_argument('--figures',default='all');p.add_argument('--workers',type=int,default=1);p.add_argument('--python')
    p.add_argument('--plan-only',action='store_true');a=p.parse_args()
    plan_=plan(a.sources,a.project,a.out,None if a.source_ids=='all' else a.source_ids.split(','),
               None if a.figures=='all' else a.figures.split(','),a.python)
    if a.plan_only:
        write(Path(a.out)/'render-plan.json',plan_)
        print(plan_['status'])
    else:execute(plan_,a.out,a.workers)
