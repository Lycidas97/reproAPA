"""Adapt reviewed measurement controls to the current run's vectors and cells."""
import argparse
import json
import shutil
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

from downstream_parallel import ordered_map
from downstream_runtime import Context,Study,STUDIES,KEY,seed
from downstream_recurrence import numerical,events_from_routes,recurrence_null
from polyapipe_effects import weights

SOURCE='notes/2026-10-02-reviewer-followup/apa_followup.py'


def collapse_common_axis(axis,counts,mapping):
    keys=['study_id','atlas_pas_id','primary_gene_id']
    chosen=axis[keys].reset_index().merge(mapping,on=keys,validate='one_to_one')
    regions=chosen[['primary_gene_id','cluster_id','common_position','common_cut_site']].drop_duplicates().sort_values('cluster_id')
    assert regions.cluster_id.is_unique
    ci=pd.Index(regions.cluster_id).get_indexer(chosen.cluster_id)
    mat=sparse.csr_matrix((np.ones(len(chosen),np.int64),(chosen['index'],ci)),shape=(len(axis),len(regions)))
    collapsed=(counts@mat).tocsr()
    assert np.array_equal(np.asarray(collapsed.sum(1)),np.asarray(counts[:,chosen['index']].sum(1)))
    regions=regions.rename(columns={'cluster_id':'atlas_pas_id','common_position':'position','common_cut_site':'cut_site'}).reset_index(drop=True)
    regions['rank']=regions.groupby('primary_gene_id').position.rank(method='first')
    regions['n']=regions.groupby('primary_gene_id').position.transform('size')
    return regions,collapsed


def restore_completed_blocks(ctx):
    if 'measurement_resume' not in ctx.config.get('resources',{}):return set()
    manifest=json.loads(ctx.resource('measurement_resume').read_text())
    assert manifest['run_id']==ctx.run_id and manifest['input_fingerprint']==ctx.config['input_fingerprint']
    completed=set()
    for block in manifest['blocks']:
        key=(block['kind'],block['study_id']);assert key not in completed
        for item in block['files']:
            source=ctx.register(item['source']);assert ctx.inputs[str(source)]==item['sha256']
            destination=ctx.out/item['relative_path']
            assert destination.is_relative_to(ctx.out) and not destination.exists()
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
        if block['kind']!='same-state':
            audit=json.loads((ctx.out/block['kind']/block['study_id']/'audit.json').read_text())
            assert audit['status']=='PASS' and audit['max_baseline_error']<1e-8
            assert (audit['bootstrap'],audit['split_half'],audit['thinning'])==(500,100,20)
            for path,value in audit['input_sha256'].items():
                source=ctx.register(path);assert ctx.inputs[str(source)]==value
        completed.add(key)
    return completed


def common_axes(axes):
    rows=[]
    for species,studies in [('mouse',['GSE103221','GSE115943']),('human',['GSE118258','GSE178325'])]:
        for gene,group in axes[axes.study_id.isin(studies)].groupby('primary_gene_id'):
            if group.study_id.nunique()!=2 or group.strand.nunique()!=1 or group.chrom.nunique()!=1:continue
            clusters=[]
            for row in group.sort_values(['cut_site','study_id','atlas_pas_id']).to_dict('records'):
                if not clusters or row['cut_site']-clusters[-1][0]['cut_site']>50:clusters.append([])
                clusters[-1].append(row)
            clusters=[c for c in clusters if {r['study_id'] for r in c}==set(studies)]
            if len(clusters)<2:continue
            cut=np.array([np.mean([r['cut_site'] for r in c]) for c in clusters])
            signed=cut*(1 if group.strand.iloc[0]=='+' else -1)
            if np.ptp(signed)==0:continue
            position=(signed-signed.min())/np.ptp(signed)
            for j,cluster in enumerate(clusters):
                for row in cluster:rows.append(dict(species=species,study_id=row['study_id'],primary_gene_id=gene,
                    atlas_pas_id=row['atlas_pas_id'],cluster_id=f'{species}|{gene}|{j}',common_position=position[j],
                    common_cut_site=cut[j],cluster_span=max(r['cut_site'] for r in cluster)-min(r['cut_site'] for r in cluster)))
    return pd.DataFrame(rows,columns=['species','study_id','primary_gene_id','atlas_pas_id','cluster_id',
                                     'common_position','common_cut_site','cluster_span'])


def measurement_block(ctx,item):
    sid,kind=item
    base=ctx.analysis/'state_vectors'
    pairs=ctx.read(base/'pair_effects.primary.tsv.gz');routes=ctx.read(base/'route_effects.primary.tsv.gz')
    states=ctx.read(base/'library_states.primary.tsv.gz');states['WPP']=states.APA
    axes=ctx.read(base/'PAS_axes.primary.tsv.gz')
    mapping=common_axes(axes)
    cycle=ctx.resource('expression_cycle_inputs')
    virtual=Path('/__current_pas_run__')
    tables={virtual/'pair_effects.tsv.gz':pairs,virtual/'route_effects.tsv.gz':routes,virtual/'library_state_metrics.tsv.gz':states}
    def read(path,**kw):
        if Path(path) in tables:return tables[Path(path)].copy()
        return ctx.read(path,**kw)
    def save(frame,path):
        path=Path(path);name=str(path.relative_to(ctx.out))
        if frame.empty and not len(frame.columns):
            schemas={'route_variants.tsv.gz':KEY+['variant','replicate','delta_WPP','MPRO'],
                'conditional_intervals.tsv.gz':KEY+['interval_excludes_zero','mpro_exceeds_noise95'],
                'pseudo_route_effects.tsv.gz':KEY+['replicate','delta_WPP','MPRO','null_kind'],
                'depth_associations.tsv.gz':KEY+['delta_WPP','delta_log_pas_depth'],
                'coverage.tsv.gz':['study_id','gene_id','status'],
                'same_state_split_summary.tsv.gz':['study_id','gene_id','route','library_key','state','split_estimable']}
            columns=schemas.get(path.name,KEY+['pseudo_kind','replicate','delta_WPP','MPRO','n_pairs'])
            frame=pd.DataFrame(columns=columns)
        return ctx.table(frame,name)
    def load(sid,kind='native'):
        study=Study(ctx,sid);a=axes[axes.study_id.eq(sid)].reset_index(drop=True)
        counts=study.counts[:,a.feature_index_1based.to_numpy()-1].tocsr()
        if kind=='common':
            a,counts=collapse_common_axis(a,counts,mapping[mapping.study_id.eq(sid)])
        return study.meta,a,counts
    def cycle_covariates(sid,meta):
        with np.load(cycle,allow_pickle=False) as archive:
            if sid+'_meta_json' not in archive:return None
            original=pd.DataFrame(json.loads(str(archive[sid+'_meta_json'])))
            species='mouse' if sid in ['GSE103221','GSE115943','GSE171127'] else 'human'
            symbols=archive[species+'_marker_symbols'].tolist();programs=json.loads(str(archive[species+'_program_json']))
            marker=archive[sid+'_marker_logcp10k'];cov=original[['metacell_id']].copy()
            for name in ['Cycle_S','Cycle_G2M']:
                cov[name]=marker[:,[symbols.index(g.upper()) for g in programs[name].split()]].mean(1)
        out=meta.merge(cov,on='metacell_id',validate='one_to_one')
        assert len(out)==len(meta)
        for route,part in out.groupby('route'):
            for name in ['Cycle_S','Cycle_G2M']:
                edges=np.unique(np.quantile(part[name],[1/3,2/3]));out.loc[part.index,name+'_bin']=np.searchsorted(edges,part[name],side='right')
        out['cycle_bin']=(out.Cycle_S_bin.astype(int)*3+out.Cycle_G2M_bin.astype(int)).astype(int)
        return out
    def pair_weights(frame):
        f=frame.copy();f['weight']=weights(f);return f
    ns=dict(OLD=virtual,HERE=ctx.out,BASE=KEY,BOOT=500,SPLIT=100,THIN=20,INPUTS=ctx.inputs,
        time=time,read=read,save=save,load_study=load,cycle_covariates=cycle_covariates,
        legacy=SimpleNamespace(pair_weights=pair_weights),rng_for=lambda *k:np.random.default_rng(seed(*k)))
    math=ctx.functions(SOURCE,['hypergeom_thin','state_record','joint_good','cycle_pair','make_route','measure'],ns)
    halves=ctx.functions('notes/2026-10-02-reviewer-followup/same_state.py',['extract'],ns)
    if kind=='same-state':halves.extract(sid)
    else:math.measure(sid,kind)
    return dict(study_id=sid,kind=kind,status='measurement_complete')


def measurement(ctx):
    base=ctx.analysis/'state_vectors'
    pairs=ctx.read(base/'pair_effects.primary.tsv.gz')
    routes=ctx.read(base/'route_effects.primary.tsv.gz')
    axes=ctx.read(base/'PAS_axes.primary.tsv.gz')
    mapping=common_axes(axes);ctx.table(mapping,'common_axes.tsv.gz')
    tasks=[(sid,kind) for sid in STUDIES if pairs.study_id.eq(sid).any()
           for kind in (['native','same-state','common'] if mapping.study_id.eq(sid).any() else ['native','same-state'])]
    ordered_map(ctx,measurement_block,tasks)
    status=[dict(study_id=sid,status='measurement_complete' if pairs.study_id.eq(sid).any() else 'no_qualified_endpoint_pairs',
                 common_axis_genes=mapping[mapping.study_id.eq(sid)].primary_gene_id.nunique()) for sid in STUDIES]
    same_state_summary(ctx,routes)
    measurement_nulls(ctx)
    ctx.table(status,'study_status.tsv')
    ctx.finish(dict(current_counts_and_candidates=True,source_bootstrap_draw_reused=True,
        literal_same_state_halves=True,without_replacement_thinning=True,new_common_axes=True,
        exact_series_nulls=True,planned_112_positions=True),bootstrap=500,split_halves=100,thinning=20)


def measurement_nulls(ctx):
    math=numerical(ctx)
    routes=pd.concat([ctx.read(p) for p in sorted(ctx.out.glob('*/*/route_variants.tsv.gz'))],ignore_index=True)
    intervals=pd.concat([ctx.read(p) for p in sorted(ctx.out.glob('native/*/conditional_intervals.tsv.gz'))],ignore_index=True)
    confident=routes[routes.variant.eq('native')&routes.replicate.eq(0)].merge(
        intervals[KEY+['interval_excludes_zero','mpro_exceeds_noise95']],on=KEY,validate='one_to_one')
    summaries=[];members=[];archives={};thresholds={'WPP':[.02,.05],'MPRO':[.1,.2]}
    for variant in ['native','interval_supported','depth20','depth40','cycle_matched','common']:
        f=confident.copy() if variant=='interval_supported' else routes[routes.variant.eq(variant)&routes.replicate.eq(0)]
        if f.empty:
            template=events_from_routes(confident,math).iloc[:0]
            events=template
        else:events=events_from_routes(f,math)
        if variant=='interval_supported':
            index=pd.MultiIndex.from_frame(events[['species','contrast','unit','gene_id']])
            for metric,column in [('WPP','interval_excludes_zero'),('MPRO','mpro_exceeds_noise95')]:
                good=f.groupby(['species','contrast','study_id','gene_id'])[column].all().reindex(index).to_numpy()
                events.loc[events.metric.eq(metric)&~good,'event']=0
        ctx.table(events,f'nulls/{variant}_study_events.tsv.gz')
        for exact in [False,True] if variant=='native' else [True]:
            s,m,st,di,a,lo=recurrence_null(events,math,exact_series=exact,tag=variant,thresholds=thresholds)
            summaries.append(s);members.append(m);archives.update(a)
    out=pd.concat(summaries,ignore_index=True)
    assert len(out)==112
    from downstream_runtime import bh
    out['q_all_followup_tests']=bh(out.empirical_p)
    out['q_within_variant_metric']=out.groupby(['variant','metric','series_matching']).empirical_p.transform(bh)
    ctx.table(out,'nulls/recurrence_summary.tsv');ctx.table(pd.concat(members),'nulls/recurrence_members.tsv.gz')
    measurement_summaries(ctx,routes,out,pd.concat(members))
    np.savez_compressed(ctx.out/'nulls/null_draws.npz',**archives)


def measurement_summaries(ctx,routes,summary,members):
    keys=['variant','metric','species','contrast','threshold'];base=members[members.variant.eq('native')&members.series_matching]
    rows=[]
    for r in summary[summary.series_matching].itertuples():
        f=members[members.series_matching];ref=base
        for k in keys:
            f=f[f[k].eq(getattr(r,k))]
            if k!='variant':ref=ref[ref[k].eq(getattr(r,k))]
        a,b=set(zip(f.gene_id,f.direction)),set(zip(ref.gene_id,ref.direction))
        rows.append({k:getattr(r,k) for k in keys}|dict(original=len(b),current=len(a),same_direction_intersection=len(a&b),status=r.status))
    ctx.table(rows,'nulls/member_retention.tsv')
    associations=pd.concat([ctx.read(p) for p in sorted(ctx.out.glob('native/*/depth_associations.tsv.gz'))])
    rows=[]
    for key,f in associations.groupby(['species','study_id','route','contrast']):
        ok=np.isfinite(f.delta_WPP)&np.isfinite(f.delta_log_pas_depth);d=f[ok]
        estimable=len(d)>=3 and d.delta_WPP.nunique()>1 and d.delta_log_pas_depth.nunique()>1
        rows.append(dict(zip(['species','study_id','route','contrast'],key))|dict(n=len(d),status='estimable' if estimable else 'insufficient_variation',
            spearman_dwui_depth=float(spearmanr(d.delta_WPP,d.delta_log_pas_depth).statistic) if estimable else np.nan))
    ctx.table(rows,'nulls/depth_associations.tsv',columns=['species','study_id','route','contrast','n','status','spearman_dwui_depth'])
    original=routes[routes.variant.eq('native')&routes.replicate.eq(0)][KEY+['delta_WPP','MPRO']]
    d=routes[routes.variant.isin(['depth20','depth40'])].merge(original,on=KEY,suffixes=('','_original'),validate='many_to_one')
    rows=[]
    for key,f in d.groupby(['variant','species','study_id','route','contrast']):
        for metric,cuts in [('WPP',[.02,.05]),('MPRO',[.1,.2])]:
            col='delta_WPP' if metric=='WPP' else 'MPRO'
            for cut in cuts:
                calls=np.where(f[col]>=cut,1,np.where((metric=='WPP')&(f[col]<=-cut),-1,0))
                ref=np.where(f[col+'_original']>=cut,1,np.where((metric=='WPP')&(f[col+'_original']<=-cut),-1,0))
                rows.append(dict(zip(['variant','species','study_id','route','contrast'],key))|dict(metric=metric,threshold=cut,
                    n_gene_draws=len(f),genes=f.gene_id.nunique(),draws=f.replicate.nunique(),classification_agreement=float((calls==ref).mean()),
                    absolute_effect_difference_median=float((f[col]-f[col+'_original']).abs().median())))
    ctx.table(rows,'nulls/thinning_stability.tsv',columns=['variant','species','study_id','route','contrast','metric','threshold','n_gene_draws','genes','draws','classification_agreement','absolute_effect_difference_median'])
    recurrence=[];math=numerical(ctx)
    for (variant,replicate),frame in routes[routes.variant.isin(['depth20','depth40'])].groupby(['variant','replicate']):
        events=events_from_routes(frame,math)
        for key,g in events[events.metric.isin(['WPP','MPRO'])].groupby(['metric','species','contrast','threshold']):
            pivot=g.pivot(index='gene_id',columns='unit',values='event')
            votes=np.maximum(pivot.eq(1).sum(1),pivot.eq(-1).sum(1))
            recurrence.append(dict(zip(['metric','species','contrast','threshold'],key))|dict(variant=variant,replicate=replicate,
                eligible_genes=int(pivot.notna().sum(1).ge(2).sum()),recurrent_genes=int(votes.ge(2).sum())))
    ctx.table(recurrence,'nulls/thinning_recurrence_all_draws.tsv',columns=['metric','species','contrast','threshold','variant','replicate','eligible_genes','recurrent_genes'])



def same_state_summary(ctx,observed):
    HERE=ctx.out;read=ctx.read;BASE=KEY
    def save(f,p):return ctx.table(f,str(p.relative_to(HERE)))
    pseudo=pd.concat([read(p) for p in sorted((HERE / 'same-state').glob('GSE*.tsv.gz'))],ignore_index=True)
    rows,noise=[] ,[]
    for kind in ['reference','target']:
        for (species,contrast),p in pseudo[pseudo.pseudo_kind.eq(kind)].groupby(['species','contrast']):
            for metric,thresholds in [('WPP',[.02,.05]),('MPRO',[.1,.2])]:
                field='delta_WPP' if metric=='WPP' else 'MPRO'
                for rep,frame in p.groupby('replicate'):
                    frame=frame[np.isfinite(frame[field])]
                    # All original series must have eligible halves; no partial-series recollapse.
                    expected=observed[observed.species.eq(species)&observed.contrast.eq(contrast)].groupby(['study_id','gene_id']).size()
                    available=frame.groupby(['study_id','gene_id']).size();good=available[available.eq(expected.reindex(available.index))].index
                    idx=pd.MultiIndex.from_frame(frame[['study_id','gene_id']]);frame=frame[idx.isin(good)]
                    sub=observed.merge(frame[BASE],on=BASE,validate='one_to_one')
                    for threshold in thresholds:
                        def recurrent(d):
                            f=d.groupby(['study_id','gene_id'])[field].agg(['min','max'])
                            plus=f['min'].ge(threshold);minus=f['max'].le(-threshold) if metric=='WPP' else pd.Series(False,index=f.index)
                            votes=pd.DataFrame({'plus':plus,'minus':minus}).groupby('gene_id').sum();nunits=f.reset_index().groupby('gene_id').study_id.nunique()
                            return int(((votes.plus>=2)|(votes.minus>=2)).sum()),int((nunits>=2).sum())
                        n,eligible=recurrent(frame);obs,oe=recurrent(sub);assert eligible==oe
                        rows.append(dict(pseudo_kind=kind,species=species,contrast=contrast,metric=metric,threshold=threshold,replicate=rep,
                                         eligible_cross_study_genes=eligible,pseudo_recurrent=n,observed_on_same_background=obs))
                values=p[field].dropna().to_numpy();noise.append(dict(pseudo_kind=kind,species=species,contrast=contrast,metric=metric,
                    n_pseudo_routes=len(values),absolute_median=float(np.median(abs(values))),absolute_p90=float(np.quantile(abs(values),.9)) if len(values) else np.nan,absolute_p95=float(np.quantile(abs(values),.95)) if len(values) else np.nan,
                    lower_threshold_fraction=float((abs(values)>=(.02 if metric=='WPP' else .1)).mean())))
    result=pd.DataFrame(rows,columns=["pseudo_kind","species","contrast","metric","threshold","replicate","eligible_cross_study_genes","pseudo_recurrent","observed_on_same_background"]);save(result,HERE / 'same-state/recurrence_draws.tsv.gz');save(pd.DataFrame(noise,columns=["pseudo_kind","species","contrast","metric","n_pseudo_routes","absolute_median","absolute_p90","absolute_p95","lower_threshold_fraction"]),HERE / 'same-state/noise_thresholds.tsv')
    summary=[]
    keys=['pseudo_kind','species','contrast','metric','threshold']
    for key,f in result.groupby(keys):
        summary.append(dict(zip(keys,key)) | dict(eligible_median=float(f.eligible_cross_study_genes.median()),pseudo_median=float(f.pseudo_recurrent.median()),pseudo_p95=float(f.pseudo_recurrent.quantile(.95)),
            observed_same_background_median=float(f.observed_on_same_background.median()),observed_exceeds_pseudo_fraction=float((f.observed_on_same_background>f.pseudo_recurrent).mean())))
    save(pd.DataFrame(summary,columns=keys+["eligible_median","pseudo_median","pseudo_p95","observed_same_background_median","observed_exceeds_pseudo_fraction"]),HERE / 'same-state/recurrence_summary.tsv');print(pd.DataFrame(summary).to_string(index=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();measurement(Context(a.config,'measurement',a.out))
