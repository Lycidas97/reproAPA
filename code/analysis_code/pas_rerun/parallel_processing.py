"""Parallelize independent study exports without changing APADETECT eligibility."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing
from pathlib import Path
import pandas as pd
from apadetect_prepare import apply_view,checked,truth
from polyapipe_prepare import save
from pipeline import read,write,digest,require


def derive_study(source,view,out,row):
    source=Path(source);out=Path(out);sid=row['study_id'];folder=source/'studies'/sid
    old=read(source/'intake.native.json');inputs={};dest=out/'studies'/sid;dest.mkdir(parents=True)
    for name,h in old['outputs'][sid].items():checked(folder/name,inputs,h)
    f=pd.read_csv(folder/'features.tsv.gz',sep='\t',low_memory=False)
    result=apply_view(f,view)
    require(result.atlas_pas_id.equals(f.atlas_pas_id),'Parallel route changed stored PAS order')
    save(result,dest/'features.tsv.gz')
    for name in old['outputs'][sid]:
        if name!='features.tsv.gz':(dest/name).symlink_to(folder/name)
    summary=dict(row,primary_pas=int(result.keep_for_apa_primary.sum()),
                 supported50_primary_pas=int((result.keep_for_apa_primary & result.matched_source_count_50nt.ge(1)).sum()))
    return dict(study_id=sid,summary=summary,inputs=inputs,outputs={p.name:digest(p) for p in dest.iterdir()})


def derive_route_parallel(source,view,out,workers=4):
    source=Path(source).resolve();out=Path(out).resolve()
    require(workers>0 and view!='base','Use a positive study concurrency and a filtered route')
    require(not out.exists() and not out.is_relative_to(source),'Use a fresh sibling route directory')
    old=read(source/'intake.native.json')
    require(old.get('analysis_route',old.get('primary_view'))=='base','Derive from base, never from another filter')
    require(read(source/'intake-complete.json')['status']=='READY_FOR_CORE','Complete base required')
    require(read(source/'annotation.complete.json')['support_reference_variant']=='v2_primary','Only v2 support')
    inputs={str(source/'intake.native.json'):digest(source/'intake.native.json'),str(Path(__file__).resolve()):digest(__file__)}
    # Record the exact unchanged eligibility implementation as well as this executor.
    import apadetect_prepare
    inputs[str(Path(apadetect_prepare.__file__).resolve())]=digest(apadetect_prepare.__file__)
    out.mkdir(parents=True);records={}
    def status(phase,**extra):
        write(out/'study-preparation-status.json',dict(phase=phase,analysis_route=view,
              study_workers=workers,completed=sorted(records),total=len(old['studies']),**extra))
    status('running')
    try:
        # Spawn is safe when several processing-family threads create study pools.
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
            futures={pool.submit(derive_study,str(source),view,str(out),r):r['study_id'] for r in old['studies']}
            for future in as_completed(futures):
                result=future.result();records[result['study_id']]=result;status('running')
        summaries=[];products={}
        for row in old['studies']:
            r=records[row['study_id']];summaries.append(r['summary']);products[r['study_id']]=r['outputs']
            for p,h in r['inputs'].items():
                require(p not in inputs or inputs[p]==h,'Study input fingerprints disagree');inputs[p]=h
        for name in ['frozen_native_membership_audit.json','frozen_native_membership_audit.tsv']:
            checked(source/name,inputs);(out/name).symlink_to(source/name)
        for p,h in inputs.items():require(digest(p)==h,'Base changed while deriving route: '+p)
        record=dict(old,run_id=old['run_id'].removesuffix('__base')+'__'+view,analysis_route=view,
                    analysis_route_policy='parallel_primary_sources',studies=summaries,inputs=inputs,outputs=products,
                    base_intake_sha256=digest(source/'intake.native.json'))
        record.pop('primary_view',None)
        write(out/'intake.native.json',record)
        write(out/'annotation.complete.json',dict(status='complete',rows=sum(s['atlas_pas'] for s in summaries),
              studies=products,support_reference_variant='v2_primary',analysis_route=view))
        write(out/'intake-complete.json',dict(status='READY_FOR_CORE',inputs=inputs,studies=summaries,
              outputs=products,analysis_route=view,shared_immutable_counts=True,new_BAM_recount=False,
              core_or_statistics_run=False,figures_replaced=False,full_manuscript_rerun_complete=False))
        status('complete')
    except Exception as error:status('failed',error=repr(error));raise
