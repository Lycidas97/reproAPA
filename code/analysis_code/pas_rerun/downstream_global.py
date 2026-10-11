"""Global-direction backgrounds and a full-primary-axis tandem UTR audit."""
import argparse
import gzip
import re

import numpy as np
import pandas as pd

from downstream_runtime import Context,Study,STUDIES,SETTINGS,KEY

GROUP=['species','study_id','route','contrast','design']
SOURCE='notes/2026-09-26-global-direction-audit/'


def describe(values):
    x=np.asarray(values,float);assert np.isfinite(x).all()
    result=dict(n=len(x),status='measured' if len(x) else 'empty_background')
    for name,fn in [('mean',np.mean),('median',np.median),('q10',lambda a:np.quantile(a,.1)),
        ('q25',lambda a:np.quantile(a,.25)),('q75',lambda a:np.quantile(a,.75)),('q90',lambda a:np.quantile(a,.9)),
        ('positive_contribution',lambda a:np.maximum(a,0).mean()),('negative_contribution',lambda a:np.minimum(a,0).mean())]:
        result[name]=fn(x) if len(x) else np.nan
    for t in [.02,.05]:
        suffix=str(t).replace('.','');p=int((x>=t).sum());n=int((x<=-t).sum())
        result.update({f'proximal_{suffix}':n,f'distal_{suffix}':p,f'distal_excess_{suffix}':(p-n)/len(x) if len(x) else np.nan})
    return result


def collapse(ctx,pairs):
    func=ctx.functions(SOURCE+'analyze.py',['collapse_pairs'],{'KEYS':KEY+['setting']})
    return func.collapse_pairs(pairs)


def backgrounds(route,orth,planned):
    frames={k:q for k,q in route.groupby(GROUP,dropna=False)};fixed={};species={};rows=[];members=[]
    for key in planned:
        q=frames.get(key,route.iloc[:0]);fixed[key]=set.intersection(*(set(q.loc[q.setting.eq(s),'gene_id']) for s in SETTINGS))
    for sp,co in sorted({(k[0],k[3]) for k in planned if k[-1] in ['true_source','same_library']}):
        keys=[k for k in planned if k[0]==sp and k[3]==co and k[-1] in ['true_source','same_library']]
        species[sp,co]=set.intersection(*(fixed[k] for k in keys))
    for key in planned:
        sp,_,_,co,design=key;q=frames.get(key,route.iloc[:0])
        paired=orth[orth.human_gene_id.isin(species.get(('human',co),set()))&orth.mouse_gene_id.isin(species.get(('mouse',co),set()))]
        choices={'fixed_settings_within_route':fixed[key]}
        if design in ['true_source','same_library']:
            choices.update(common_within_species_all_settings=species[sp,co],common_one_to_one_ortholog_all_settings=set(paired[sp+'_gene_id']))
        for setting in SETTINGS:
            z=q[q.setting.eq(setting)]
            for label,ids in dict(available=set(z.gene_id),**choices).items():
                b=z[z.gene_id.isin(ids)];r=dict(zip(GROUP,key),setting=setting,background=label,**describe(b.delta))
                if 'delta_nt' in b:r.update(mean_delta_nt=b.delta_nt.mean(),median_delta_nt=b.delta_nt.median())
                rows.append(r)
                members.extend(dict(zip(GROUP,key),setting=setting,background=label,gene_id=g) for g in sorted(ids))
    return pd.DataFrame(rows),pd.DataFrame(members,columns=GROUP+['setting','background','gene_id'])


def observed_audit(ctx):
    effects=[];inventories=[];references=[]
    for sid in STUDIES:
        study=Study(ctx,sid);meta=study.meta
        for (route,node,time),q in meta.groupby(['route','node','experimental_time']):
            inventories.append(dict(species=q.species.iloc[0],study_id=sid,route=route,state=node,time=time,
                n_mc=len(q),n_cells=int(q.n_cells.sum()),libraries=q.library_key.nunique(),
                true_source=node=='S' and str(time) in ['D0','source','Stage_I_day_0']))
        for setting in ['supported50_p08','expanded_primary_p08']:
            libs=study.core('library_states',setting);nodes=study.core('node_positions',setting)
            for route,q in nodes.groupby('route'):
                w=q.pivot(index='gene_id',columns='node',values='position_index').reindex(columns=['S','P','M'])
                for contrast,reference in [('P_S','S'),('P_M','M')]:
                    z=(w.P-w[reference]).dropna()
                    effects.extend(dict(species=q.species.iloc[0],study_id=sid,route=route,setting=setting,gene_id=g,
                        contrast=contrast,reference='all_observed_'+reference,delta=d) for g,d in z.items())
                design=study.design[study.design.route.eq(route)][['library_key','time']].drop_duplicates()
                z=libs[libs.route.eq(route)].merge(design,on='library_key',validate='many_to_one')
                source=z[z.node.eq('S')]
                references.extend(dict(species=q.species.iloc[0],study_id=sid,route=route,setting=setting,
                    library_key=lib,time=t,genes=len(b),documented_source=str(t) in ['D0','source','Stage_I_day_0'])
                    for (lib,t),b in source.groupby(['library_key','time']))
                true=source[source.time.astype(str).isin(['D0','source','Stage_I_day_0'])].groupby('gene_id').position_index.mean()
                d=(w.P-true).dropna()
                effects.extend(dict(species=q.species.iloc[0],study_id=sid,route=route,setting=setting,gene_id=g,
                    contrast='P_S',reference='documented_S_profile_library_equal',delta=v) for g,v in d.items())
                both=pd.concat([w.P.rename('P'),w.S.rename('S'),true.rename('source')],axis=1).dropna()
                for ref in ['S','source']:
                    effects.extend(dict(species=q.species.iloc[0],study_id=sid,route=route,setting=setting,gene_id=g,
                        contrast='P_S',reference='matched_genes_'+ref,delta=v) for g,v in (both.P-both[ref]).items())
    ctx.table(effects,'observed_state_effects.tsv.gz',columns=['species','study_id','route','setting','gene_id','contrast','reference','delta'])
    ctx.table(inventories,'state_time_inventory.tsv');ctx.table(references,'source_library_inventory.tsv')


def global_direction(ctx):
    pairs=ctx.read(ctx.analysis/'state_vectors/pair_effects.tsv.gz');route,dates=collapse(ctx,pairs)
    assert np.isfinite(route.delta).all()
    original=ctx.read(ctx.analysis/'state_vectors/route_effects.tsv.gz')
    compare=route.merge(original[KEY+['setting','delta_DWUI']],on=KEY+['setting'],validate='one_to_one')
    assert len(compare)==len(route)==len(original) and np.allclose(compare.delta,compare.delta_DWUI,rtol=0,atol=1e-12)
    orth=ctx.read(ctx.resource('ortholog'))
    planned=[]
    for sid in STUDIES:
        study=Study(ctx,sid)
        planned.extend((study.meta.species.iloc[0],sid,*key) for key in study.pairs[['route','contrast','design']].drop_duplicates().itertuples(index=False,name=None))
    planned=sorted(set(planned));summary,members=backgrounds(route,orth,planned)
    ctx.table(route,'qualified_route_effects.tsv.gz');ctx.table(dates,'qualified_date_effects.tsv.gz')
    ctx.table(summary,'global_direction_summary.tsv');ctx.table(members,'background_gene_membership.tsv.gz')
    ctx.table([dict(zip(GROUP,k)) for k in planned],'planned_route_contrasts.tsv')
    leave=[];sources=[];decomposition=[]
    for key,q in dates[dates.setting.eq('supported50_p08')].groupby(GROUP):
        w=q.pivot(index='gene_id',columns='time_a',values='APA_delta').dropna()
        for omitted in ['none',*w.columns] if len(w.columns)>1 else ['none']:
            z=w if omitted=='none' else w.drop(columns=omitted)
            leave.append(dict(zip(GROUP,key),omitted_date=omitted,available_dates=len(w.columns),**describe(z.mean(1))))
    for key,q in pairs[pairs.setting.eq('supported50_p08')&pairs.design.eq('true_source')].groupby(GROUP):
        design=Study(ctx,key[1]).pairs
        libraries=sorted(design[design.route.eq(key[2])&design.contrast.eq(key[3])&design.design.eq(key[4])].library_key_b.unique())
        joint=q.groupby(['gene_id','library_key_a']).library_key_b.transform('nunique').eq(len(libraries))
        fixed=set.intersection(*(set(q.loc[q.library_key_b.eq(lib),'gene_id']) for lib in libraries))
        for lib in libraries:
            for background,z in [('available',q[q.library_key_b.eq(lib)]),('fixed_genes_all_sources',q[q.library_key_b.eq(lib)&q.gene_id.isin(fixed)]),
                ('matched_target_libraries_all_sources',q[q.library_key_b.eq(lib)&joint])]:
                r,_=collapse(ctx,z);sources.append(dict(zip(GROUP,key),source_library=lib,background=background,**describe(r.delta)))
    for key,q in route.groupby(GROUP):
        a=q[q.setting.eq('supported50_p08')].set_index('gene_id').delta;b=q[q.setting.eq('expanded_primary_p08')].set_index('gene_id').delta
        ids=a.index.intersection(b.index)
        for label,z in [('supported_common',a.loc[ids]),('expanded_common',b.loc[ids]),('supported_only',a.loc[a.index.difference(ids)]),('expanded_only',b.loc[b.index.difference(ids)])]:
            decomposition.append(dict(zip(GROUP,key),subset=label,**describe(z)))
    source=pd.DataFrame(sources,columns=GROUP+['source_library','background',*describe([])])
    stable=[]
    for key,q in summary.groupby(GROUP):
        r=dict(zip(GROUP,key),scenarios=len(q),empty_scenarios=int(q.n.eq(0).sum()),mean_min=q['mean'].min(),mean_max=q['mean'].max())
        for metric in ['mean','median','distal_excess_002','distal_excess_005']:
            values=q[metric];sign=int(np.sign(values.iloc[0])) if values.notna().all() and np.sign(values).nunique()==1 else 0
            r[metric+'_stable_sign']=sign
            z=source[source.route.eq(key[2])&source.contrast.eq(key[3])&source.background.ne('available')]
            r[metric+'_source_checked_sign']=sign if key[4]!='true_source' or (len(z) and z[metric].notna().all() and np.sign(z[metric]).eq(sign).all()) else 0
        stable.append(r)
    ctx.table(leave,'leave_target_date_summary.tsv');ctx.table(source,'single_source_library_summary.tsv')
    ctx.table(decomposition,'pas_layer_coverage_decomposition.tsv');ctx.table(stable,'direction_stability.tsv')
    observed_audit(ctx)
    ctx.finish(dict(new_qualified_pairs=True,route_weights_independently_recovered=True,planned_four_settings_in_intersections=True,
        empty_backgrounds_explicit=True,source_and_observed_S_separate=True,native_counts_reused_without_BAM_recount=True),
        counting_rule_sensitivity='separate legacy-method dependency; not run or substituted into native results')


def fits_full_primary(axis,utrs):
    if len(axis)<2 or axis.chrom.nunique()!=1 or axis.strand.nunique()!=1:return utrs.iloc[:0]
    cut=axis.cut_site.to_numpy()+int(axis.strand.iloc[0]=='-')
    return utrs[utrs.chrom.astype(str).str.removeprefix('chr').eq(str(axis.chrom.iloc[0]).removeprefix('chr'))&
        utrs.strand.eq(axis.strand.iloc[0])&utrs.utr_start.le(cut.min())&utrs.utr_end.ge(cut.max())]


def tandem(ctx):
    effects=ctx.read(ctx.analysis/'global/qualified_route_effects.tsv.gz');planned=ctx.read(ctx.analysis/'global/planned_route_contrasts.tsv')
    function=ctx.functions(SOURCE+'tandem_utr.py',['terminal_utrs'],dict(gzip=gzip,re=re,INPUTS={},sha=lambda p:ctx.inputs[str(ctx.register(p))]))
    utrs=pd.concat([function.terminal_utrs(ctx.resource(sp+'_ensembl_gtf'),set(effects.loc[effects.species.eq(sp),'gene_id']),sp) for sp in ['human','mouse']],ignore_index=True)
    ctx.table(utrs,'annotation_terminal_utrs.tsv.gz');bygene=dict(tuple(utrs.groupby('gene_id')));rows=[]
    for sid in STUDIES:
        study=Study(ctx,sid)
        for gene,axis in study.genes.items():
            u=bygene.get(gene,utrs.iloc[:0]);fit=fits_full_primary(axis,u)
            rows.append(dict(species=study.meta.species.iloc[0],study_id=sid,gene_id=gene,primary_pas=len(axis),
                native_axis_span_nt=axis.cut_site.max()-axis.cut_site.min(),tandem_utr_compatible=not fit.empty,
                compatible_transcripts=';'.join(sorted(fit.transcript_id)),coordinate_rule='interbase cut to 1-based retained nucleotide: plus +0, minus +1'))
    membership=pd.DataFrame(rows);ctx.table(membership,'tandem_utr_eligibility.tsv.gz')
    joined=effects.merge(membership,on=['species','study_id','gene_id'],validate='many_to_one');assert len(joined)==len(effects)
    x=joined[joined.tandem_utr_compatible].copy();x['delta_nt']=x.delta*x.native_axis_span_nt
    assert np.sign(x.delta).eq(np.sign(x.delta_nt)).all()
    ctx.table(x,'tandem_utr_effects.tsv.gz')
    summary,genes=backgrounds(x,ctx.read(ctx.resource('ortholog')),list(planned[GROUP].itertuples(index=False,name=None)))
    ctx.table(summary,'tandem_utr_summary.tsv');ctx.table(genes,'tandem_background_gene_membership.tsv.gz')
    ctx.finish(dict(full_primary_axis_defines_tandem=True,same_tandem_membership_across_layers=True,stop_codon_in_terminal_exon=True,
        minus_retained_base_converted=True,delta_nt_uses_native_span=True,empty_planned_backgrounds_preserved=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--node',choices=['global','tandem'],required=True)
    a=p.parse_args();ctx=Context(a.config,a.node,a.out);(global_direction if a.node=='global' else tandem)(ctx)
