"""Native/annotation site sensitivity and fresh R/E terminal-family backtracking."""
import argparse
import copy
import json
import shutil

import numpy as np
import pandas as pd

from downstream_parallel import ordered_map
from downstream_runtime import Context,Study,local_state_labels
from downstream_regulatory import numerical,effect,F5,ROUTES,strict_interval_pas

KEY=['species','route','gene_id','category','chrom','strand','start0','end0','contrast']
COORDINATE_TABLES=['pas_coordinate_mapping.tsv.gz','pas_representative_distance_summary.tsv',
    'family_migration.tsv','site_classification_sensitivity.tsv.gz','site_classification_summary.tsv']


def changed(frame,prefix):
    lo20=frame['eq20_lo_'+prefix];hi20=frame['eq20_hi_'+prefix]
    lo50=frame['eq50_lo_'+prefix];hi50=frame['eq50_hi_'+prefix]
    return (lo20.ge(.1)&lo50.ge(.1))|(hi20.le(-.1)&hi50.le(-.1))


def prepare_family(study,raw,annotation,mapping,strict):
    f=study.f;lookup=f.set_index('atlas_pas_id').feature_index_1based
    indexes=lookup.loc[[p['atlas_pas_id'] for p in raw['pas']]].to_numpy()-1
    family=study.counts[:,indexes].toarray()
    gene=f[strict_interval_pas(f,strict)&f.primary_gene_id.eq(raw['gene_id'])&f.matched_source_count_50nt.ge(1)]
    total=np.asarray(study.counts[:,gene.feature_index_1based.to_numpy()-1].sum(1)).ravel()
    amap=mapping.set_index('native_pas_id').annotation_id
    projection=np.array([[int(amap[p['atlas_pas_id']]==q['atlas_pas_id']) for q in annotation['pas']] for p in raw['pas']])
    assert (projection.sum(1)==1).all()
    return family,total,projection


def rebuild(study,raw,annotation,mapping,meta,strict,prepared=None):
    native=copy.deepcopy(raw);annotated=copy.deepcopy(annotation)
    family,total,projection=prepare_family(study,raw,annotation,mapping,strict) if prepared is None else prepared
    native['profiles']={};annotated['profiles']={}
    for state,m in meta.groupby('node'):
        libraries=[];denominators=[]
        for lib,g in m.groupby('library_key'):
            if len(g)<3 or g.n_cells.sum()<30:continue
            counts=family[g.index].sum(0);denom=int(total[g.index].sum());value=int(counts.sum())
            libraries.append(dict(library=lib,total=value,coverage=value/denom if denom else 0.,counts=counts,
                usage=counts/value if value else np.zeros(len(counts))))
            denominators.append(denom)
        if not libraries:continue
        counts=sum(l['counts'] for l in libraries);value=int(counts.sum());denom=sum(denominators)
        native['profiles'][state]=dict(total=value,coverage=value/denom if denom else 0.,counts=counts,
            usage=counts/value if value else np.zeros(len(counts)),libraries=libraries)
        profile=copy.deepcopy(native['profiles'][state])
        for obj in [profile]+profile['libraries']:
            before=np.asarray(obj['counts']);after=before@projection;assert before.sum()==after.sum()
            obj['counts']=after;obj['usage']=after/after.sum() if after.sum() else np.zeros(len(after))
        annotated['profiles'][state]=profile
    return native,annotated


def coordinate_inventory(ctx):
    cls=ctx.functions('notes/2026-09-29-core-sensitivity/figure5_coordinate_audit.py',['structural_status'])
    mapping=ctx.read(ctx.analysis/'intervals/pas_mapping.tsv.gz')
    mapping['representative_distance']=(mapping.native_cut_interbase-mapping.annotation_cut).abs()
    ctx.table(mapping,'pas_coordinate_mapping.tsv.gz')
    summary=mapping.groupby(['species','route']).agg(PAS=('native_pas_id','size'),genes=('gene_id','nunique'),
        representative_distance_median=('representative_distance','median'),representative_distance_max=('representative_distance','max'),
        source_distance_max=('distance_nt','max'),representative_gt50=('representative_distance',lambda v:int(v.gt(50).sum()))).reset_index()
    ctx.table(summary,'pas_representative_distance_summary.tsv')
    effects=ctx.read(ctx.analysis/'intervals/model_effects.tsv');effects['single_pas_after_merge']=effects.annotated_pas.lt(2)
    ctx.table(effects,'family_migration.tsv')
    joins=[]
    for sp in ROUTES:
        a=ctx.read(ctx.analysis/f'site_scan/native/full-cohort-v2/data/{sp}/site_effects.tsv.gz')
        b=ctx.read(ctx.analysis/f'site_scan/full-cohort-v2/data/{sp}/site_effects.tsv.gz')
        for frame in [a,b]:
            assert not frame.duplicated(KEY).any()
            frame['structure']=frame.apply(cls.structural_status,axis=1) if len(frame) else pd.Series(dtype=str)
        z=a.merge(b,on=KEY,how='outer',suffixes=('_native','_annotation'),indicator=True,validate='one_to_one')
        z['strict_changed_native']=changed(z,'native');z['strict_changed_annotation']=changed(z,'annotation')
        z['structure_changed']=z._merge.eq('both')&z.structure_native.ne(z.structure_annotation)
        z['strict_change_changed']=z.strict_changed_native.ne(z.strict_changed_annotation);joins.append(z)
    joined=pd.concat(joins,ignore_index=True);ctx.table(joined,'site_classification_sensitivity.tsv.gz')
    ctx.table(joined.groupby(['species','route','contrast','category']).agg(native_rows=('_merge',lambda v:int(v.ne('right_only').sum())),
        annotation_rows=('_merge',lambda v:int(v.ne('left_only').sum())),shared_rows=('_merge',lambda v:int(v.eq('both').sum())),
        structure_changed=('structure_changed','sum'),strict_change_changed=('strict_change_changed','sum')).reset_index(),'site_classification_summary.tsv')


def backtrack_route(ctx,item):
    sp,sid,route=item
    _,migration,_=numerical(ctx)
    onset=ctx.functions(F5+'site-recurrence-v1/supplement-compressed-v1/prepare.py',['first_observed'])
    prepare=ctx.functions('notes/2026-10-01-mpro-biology/legacy_prepare_pilot.py',['strict_pas'])
    mapping=ctx.read(ctx.analysis/'intervals/pas_mapping.tsv.gz')
    regions=ctx.read(ctx.analysis/'intervals/all_coordinate_regions.tsv.gz');claims=ctx.read(ctx.analysis/'dwui_pool/claim_index.tsv')
    claimkeys=set(claims[['species','gene_id','contrast']].itertuples(index=False,name=None))
    measurements=[];first=[];coverage=[];endpoint_checks=0
    study=Study(ctx,sid)
    if sp=='mouse':local=local_state_labels(ctx.read(study.root/'metacells.tsv.gz'),study.meta.metacell_id)
    native=json.loads(ctx.register(ctx.analysis/f'intervals/native/{sp}/{route}/models.json').read_text())
    annotated=json.loads(ctx.register(ctx.analysis/f'intervals/annotation/{sp}/{route}/models.json').read_text())
    for gid,old in native.items():
        if len(annotated[gid]['pas'])<2:continue
        rr=regions[regions.route.eq(route)&regions.gene_id.eq(gid)].drop_duplicates(['chrom','strand','start0','end0'])
        if rr.empty:continue
        base=study.meta[study.meta.route.eq(route)&study.meta.purity.ge(.8)]
        gene_mapping=mapping[mapping.route.eq(route)&mapping.gene_id.eq(gid)]
        prepared=prepare_family(study,old,annotated[gid],gene_mapping,prepare.strict_pas)
        variants=['baseline','exclude_C012','exclude_L00','exclude_both'] if sp=='mouse' else ['baseline']
        for variant in variants:
            meta=base.copy()
            if sid=='GSE103221' and variant in ['exclude_C012','exclude_both']:
                meta=meta[~(meta.node.eq('E')&(local[meta.index]=='C012'))]
            if sid=='GSE115943' and variant in ['exclude_L00','exclude_both']:
                meta=meta[~(meta.node.eq('R')&(local[meta.index]=='L00'))]
            _,model=rebuild(study,old,annotated[gid],gene_mapping,meta,prepare.strict_pas,prepared)
            if variant=='baseline':
                for state,profile in model['profiles'].items():
                    prior=annotated[gid]['profiles'][state]
                    np.testing.assert_array_equal(profile['counts'],prior['counts'])
                    assert abs(profile['coverage']-prior['coverage'])<1e-12
            else:
                for state in ['S','P','M']:
                    if state in model['profiles']:
                        np.testing.assert_array_equal(model['profiles'][state]['counts'],annotated[gid]['profiles'][state]['counts']);endpoint_checks+=1
            stats=migration.state_stats(model);cuts=np.array([p['cut_site'] for p in model['pas']])
            for r in rr.itertuples():
                mask=migration.mask_for(r.start0,r.end0,r.strand,cuts);values={}
                for cutoff in [20,50]:
                    values[cutoff]=np.array([float(stats[cutoff,s][0]@mask) for s in 'SREP'])
                    for state in 'SREPM':
                        usage,n,count=stats[cutoff,state];value=float(usage@mask)
                        measurements.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=old['gene_symbol'],region_id=r.region_id,
                            chrom=r.chrom,strand=r.strand,start0=r.start0,end0=r.end0,variant=variant,cutoff=cutoff,state=state,inclusion=value,
                            qualified_libraries=n,include_count=float(count@mask),exclude_count=float(count@(1-mask)),
                            status='measured' if np.isfinite(value) else 'not_measurable'))
                for threshold in [.05,.10,.15,.20]:
                    first.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=old['gene_symbol'],region_id=r.region_id,variant=variant,
                        threshold=threshold,first_observed=onset.first_observed(values,threshold),
                        DWUI_SP_pool=(sp,gid,'SP') in claimkeys,union_pool=True))
            for state in 'SREPM':
                p=model['profiles'].get(state,{});coverage.append(dict(species=sp,route=route,gene_id=gid,variant=variant,state=state,
                    pooled_UMI=p.get('total',0),family_coverage=p.get('coverage',np.nan),qualified20=stats[20,state][1],qualified50=stats[50,state][1]))
    return measurements,first,coverage,endpoint_checks


def run(ctx):
    if 'coordinate_resume' in ctx.config.get('resources',{}):
        manifest=json.loads(ctx.resource('coordinate_resume').read_text())
        assert manifest['run_id']==ctx.run_id and manifest['input_fingerprint']==ctx.config['input_fingerprint']
        assert {r['name'] for r in manifest['files']}==set(COORDINATE_TABLES)
        for row in manifest['files']:
            source=ctx.register(row['source']);assert ctx.inputs[str(source)]==row['sha256']
            destination=ctx.out/row['name'];assert not destination.exists();shutil.copy2(source,destination)
        for path,value in manifest['upstream_inputs'].items():
            source=ctx.register(path);assert ctx.inputs[str(source)]==value
    else:coordinate_inventory(ctx)
    measurements=[];first=[];coverage=[];endpoint_checks=0
    for m,f,c,n in ordered_map(ctx,backtrack_route,[(sp,sid,route) for sp,routes in ROUTES.items() for sid,route in routes]):
        measurements.extend(m);first.extend(f);coverage.extend(c);endpoint_checks+=n
    ctx.table(measurements,'region_state_inclusion.tsv.gz',columns=['species','route','gene_id','region_id','variant','cutoff','state','inclusion'])
    ctx.table(first,'first_crossing.tsv.gz',columns=['species','route','gene_id','region_id','variant','threshold','first_observed'])
    ctx.table(coverage,'state_family_coverage.tsv.gz',columns=['species','route','gene_id','variant','state','pooled_UMI','family_coverage'])
    ctx.finish(dict(native_annotation_site_inventory_outer_join=True,source_and_representative_distances_separate=True,
        many_to_one_and_single_PAS_audited=True,baseline_rebuilt_from_native_counts=True,C012_L00_exclusions_requalified=True,
        first_crossing_all_four_cutoffs=True,missing_intermediates_remain_incomplete=True),unchanged_endpoint_count_checks=endpoint_checks)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'coordinate_states',a.out))
