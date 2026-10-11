"""Rebuild named-locus inclusion and influence diagnostics from current families."""
import json
import numpy as np
import pandas as pd

from downstream_runtime import Study
from downstream_regulatory import ROUTES,numerical

AUDITED={'mouse':['Usp3','Tgfbr2','Rac1','Dnajc5','Cds2','Sbno1'],
         'human':['FBXW7','SMAD2','SPARC','RAB21']}
FEATURED={'HMGN2','EIF4EBP2','RAD23B'}


def named_sites(sites,gene):
    if gene=='EIF4EBP2':return sites.factor.str.upper().eq('PCBP1')
    motif=sites.factor.eq('AAGTGCT') if gene=='RAD23B' else sites.names.str.contains('miR-23a|miR-23b',case=False,regex=True,na=False)
    return sites.category.eq('miRNA')&motif


def mean_summary(frame,value):
    x=frame[value].to_numpy(float);n=len(x)
    return dict(n_libraries=n,equal_library=float(x.mean()) if n else np.nan,
        time_equal=float(frame.groupby('time')[value].mean().mean()) if n else np.nan,
        minimum=float(x.min()) if n else np.nan,maximum=float(x.max()) if n else np.nan,
        leave_one_min=float(((x.sum()-x)/(n-1)).min()) if n>1 else np.nan,
        leave_one_max=float(((x.sum()-x)/(n-1)).max()) if n>1 else np.nan,
        interval_type='observed_range_and_leave_one_library_influence_not_CI')


def qualified(model,state,cutoff):
    return {a['library']:a for a in model['profiles'].get(state,{}).get('libraries',[]) if a['total']>=cutoff and a['coverage']>=.5}


def weighted_state(frame):
    if frame.empty:return np.nan
    return float(frame.groupby(['time','sample_unit']).value.mean().groupby('time').mean().mean())


def fixed_source_influence(ctx):
    library=ctx.read(ctx.analysis/'intermediate/library_state_values.tsv.gz')
    library=library[library.gene_symbol.str.upper().eq('EIF4EBP2')&library.setting.eq('supported50_p08')&library.route.isin(['CIP','WOT'])]
    rows=[];contributors=[]
    for route,group in library.groupby('route'):
        sid=group.study_id.iloc[0];study=Study(ctx,sid)
        assert group.gene_id.nunique()==1
        axis=study.axis(group.gene_id.iloc[0]);data=study.gene_counts(axis)
        source=set(study.pairs.loc[study.pairs.route.eq(route)&study.pairs.design.eq('true_source'),'library_key_b'])
        group=group[group.state.isin(list('REP'))|(group.state.eq('S')&group.library_key.isin(source))]
        support={(a.library_key,a.state):study.state(data,route,a.library_key,a.state) for a in group.itertuples()}
        assert all(q['endpoint_eligible'] for q in support.values())
        contributors.append(group)
        for kind,column in [('library','library_key'),('observed_sampling_unit','sample_unit')]:
            for omitted in [None,*sorted(group[column].unique())]:
                q=group if omitted is None else group[group[column].ne(omitted)]
                a=q[q.state.eq('S')];reference=weighted_state(a)
                for state in 'REP':
                    b=q[q.state.eq(state)];target=weighted_state(b)
                    participating=[support[a.library_key,a.state] for a in pd.concat([a,b]).itertuples()]
                    counts=sum((v['site_counts'] for v in participating),np.zeros(len(axis),np.int64))
                    metacells=sum((v['site_mc'] for v in participating),np.zeros(len(axis),np.int64))
                    joint=int(((counts>=10)&(metacells>=2)).sum())
                    status='source_or_target_empty' if not np.isfinite(reference+target) else 'joint_PAS_support_failed' if joint<2 else 'measured'
                    rows.append(dict(study_id=sid,route=route,gene='EIF4EBP2',state=state,omission_kind=kind,
                        omitted='none' if omitted is None else omitted,reference_DWUI=reference,target_DWUI=target,
                        delta_DWUI=target-reference if status=='measured' else np.nan,n_joint_pas=joint,
                        n_source_libraries=len(a),n_target_libraries=len(b),status=status))
    ctx.table(rows,'fixed_source_DWUI_influence.tsv',columns=['route','gene','state','omission_kind','omitted','delta_DWUI','status'])
    ctx.table(pd.concat(contributors,ignore_index=True) if contributors else library.iloc[:0],'fixed_source_DWUI_contributors.tsv.gz')


def feature_libraries(model,mask,times,cutoff):
    rows=[];paired=[]
    for state in 'SREPM':
        for lib,a in qualified(model,state,cutoff).items():
            counts=np.asarray(a['counts']);assert counts.sum()==a['total']
            rows.append(dict(library_key=lib,state=state,time=times[lib],inclusion=float((counts/counts.sum())@mask),family_UMI=a['total'],coverage=a['coverage']))
    m,p=qualified(model,'M',cutoff),qualified(model,'P',cutoff)
    for lib in sorted(set(m)&set(p)):
        mi=float(np.asarray(m[lib]['counts'])@mask/m[lib]['total']);pi=float(np.asarray(p[lib]['counts'])@mask/p[lib]['total'])
        paired.append(dict(library_key=lib,time=times[lib],reference_inclusion=mi,target_inclusion=pi,delta=pi-mi))
    return pd.DataFrame(rows,columns=['library_key','state','time','inclusion','family_UMI','coverage']),pd.DataFrame(paired,columns=['library_key','time','reference_inclusion','target_inclusion','delta'])


def run(ctx):
    old=ctx.read(ctx.resource('legacy_ten_loci'))
    assert {(a.species,a.gene) for a in old.itertuples()}=={(sp,g) for sp,genes in AUDITED.items() for g in genes}
    fixed_source_influence(ctx)
    _,migration,_=numerical(ctx)
    regions=ctx.read(ctx.analysis/'intervals/all_coordinate_regions.tsv.gz')
    sites=ctx.read(ctx.analysis/'site_statistics/site_factor_inventory.tsv.gz')
    changed=ctx.read(ctx.analysis/'site_statistics/all_changed_region_pairs.tsv.gz')
    partial=ctx.read(ctx.analysis/'site_statistics/cross_study_partial_retention.tsv')
    selected=ctx.read(ctx.analysis/'candidate_screen/example_selection.tsv')
    featured=FEATURED|set(selected.gene_symbol.str.upper())
    librows=[];pairedrows=[];states=[];effects=[];features=[];axes=[];eligibility=[];audit=[]
    for sp,routes in ROUTES.items():
        for gene in AUDITED[sp]:
            q=changed[changed.species.eq(sp)&changed.gene_symbol_a.str.upper().eq(gene.upper())]
            z=partial[partial.species.eq(sp)&partial.gene_symbol_a.str.upper().eq(gene.upper())] if 'gene_symbol_a' in partial else partial.iloc[:0]
            for co in ['SP','MP']:
                f=q[q.contrast.eq(co)];p=z[z.contrast.eq(co)]
                audit.append(dict(species=sp,gene=gene,contrast=co,matched_region_pairs=len(f),genes=f.gene_id.nunique(),
                    effect_a_min=f.effect20_a.min(),effect_a_max=f.effect20_a.max(),effect_b_min=f.effect20_b.min(),effect_b_max=f.effect20_b.max(),
                    exact_partial_recurrence=bool(len(p) and p.exact_partial_recurrence.any()),status='measured' if len(f) else 'no_new_shared_changed_region'))
        for sid,route in routes:
            study=Study(ctx,sid);meta=study.meta[study.meta.route.eq(route)]
            assert meta.groupby('library_key').experimental_time.nunique(dropna=False).eq(1).all()
            times=meta.groupby('library_key').experimental_time.first().astype(str).to_dict()
            models=json.loads(ctx.register(ctx.analysis/f'intervals/annotation/{sp}/{route}/models.json').read_text())
            matched={m['gene_symbol'].upper():m for m in models.values() if m['gene_symbol'].upper() in featured}
            for name in sorted(featured):
                model=matched.get(name);eligibility.append(dict(species=sp,study_id=sid,route=route,gene=name,
                    status='new_terminal_family' if model else 'not_in_new_eligible_terminal_family'))
                if model is None:continue
                gid=model['gene_id'];cuts=np.array([a['cut_site'] for a in model['pas']]);strand=model['family']['strand']
                base=dict(species=sp,study_id=sid,route=route,gene_id=gid,gene=name)
                axes.extend(dict(base,atlas_pas_id=p['atlas_pas_id'],cut_interbase=p['cut_interbase'],strand=strand,chrom=p['chrom']) for p in model['pas'])
                rr=regions[regions.route.eq(route)&regions.gene_id.eq(gid)].drop_duplicates(['start0','end0'])
                selected=[dict(start0=int(a.start0),end0=int(a.end0),kind='interval',label='interval') for a in rr.itertuples()]
                ss=sites[sites.route.eq(route)&sites.gene_id.eq(gid)]
                for a in ss.drop_duplicates(['start0','end0','factor']).itertuples():
                    selected.append(dict(start0=int(a.start0),end0=int(a.end0),kind='site',label=a.factor))
                for f in selected:
                    feature_id=f'{route}|{gid}|{f["kind"]}|{f["start0"]}-{f["end0"]}|{f["label"]}'
                    info=dict(base,**f,feature_id=feature_id,strand=strand,chrom=model['family']['chrom'])
                    mask=migration.mask_for(f['start0'],f['end0'],strand,cuts);features.append(dict(info,mask=json.dumps(mask.tolist())))
                    for cutoff in [20,50]:
                        libs,pairs=feature_libraries(model,mask,times,cutoff);key=dict(info,cutoff=cutoff)
                        librows.extend(libs.assign(**key).to_dict('records'));pairedrows.extend(pairs.assign(**key).to_dict('records'))
                        for state in 'SREPM':
                            q=libs[libs.state.eq(state)];states.append(dict(key,state=state,**mean_summary(q,'inclusion')))
                        for co in ['SP','MP']:
                            if co=='MP':row=mean_summary(pairs,'delta')
                            else:
                                a=libs[libs.state.eq('S')];b=libs[libs.state.eq('P')]
                                base_effect=b.inclusion.mean()-a.inclusion.mean();time_effect=b.groupby('time').inclusion.mean().mean()-a.groupby('time').inclusion.mean().mean()
                                values=[]
                                for omitted in sorted(set(a.library_key)|set(b.library_key)):
                                    aa=a[a.library_key.ne(omitted)];bb=b[b.library_key.ne(omitted)]
                                    if len(aa) and len(bb):values.append(bb.inclusion.mean()-aa.inclusion.mean())
                                row=dict(n_libraries=len(a)+len(b),n_reference=len(a),n_target=len(b),equal_library=base_effect,time_equal=time_effect,
                                    leave_one_min=min(values) if values else np.nan,leave_one_max=max(values) if values else np.nan,
                                    interval_type='leave_one_physical_library_influence_not_CI')
                            effects.append(dict(key,contrast=co,**row))
    for rows,name,cols in [(features,'inclusion_features.tsv',['feature_id','gene','route','kind','start0','end0']),
        (axes,'annotated_example_PAS.tsv',['gene','route','cut_interbase']),
        (librows,'inclusion_library_states.tsv.gz',['feature_id','gene','route','cutoff','state','library_key','inclusion']),
        (pairedrows,'inclusion_paired_libraries.tsv.gz',['feature_id','gene','route','cutoff','library_key','reference_inclusion','target_inclusion','delta']),
        (states,'inclusion_state_summary.tsv',['feature_id','gene','route','cutoff','state','equal_library']),
        (effects,'inclusion_effect_sensitivity.tsv',['feature_id','gene','route','cutoff','contrast','equal_library']),
        (eligibility,'inclusion_eligibility.tsv',None),(audit,'ten_locus_migration.tsv',None)]:ctx.table(rows,name,columns=cols)
