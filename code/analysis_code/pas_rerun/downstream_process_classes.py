"""Rebuild the common-reference R/E process cohort and its fixed-cohort sensitivities."""
import numpy as np
import pandas as pd

SOURCE='analysis_code/recurrence_correction_2026-09-25/sources/MOUSE/analysis/figure3-recurrence-common-reference-2026-09-25/analyze.py'
AXES={'mouse':'GSE103221__GSE115943','human':'GSE118258__GSE178325'}
STUDIES={'mouse':('GSE103221','GSE115943'),'human':('GSE118258','GSE178325')}
ROUTES={'mouse':{'GSE103221':'CIP','GSE115943':'WOT'},'human':{'GSE118258':'GSE118258','GSE178325':'178325_0618'}}


def process_classes(ctx,states):
    math=ctx.functions(SOURCE,['weighted_window','aggregate_fixed_reference','endpoint_table','signed_class','process_tables'],
        dict(AXES=AXES,STUDIES=STUDIES,PROCESS_ROUTES=ROUTES,SOURCE_TIMES={'D0','source','Stage_I_day_0'}))
    members=ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz')
    members=members[members.metric.eq('WPP')&members.threshold.eq(.02)&members.contrast.eq('P_S')&members.series_matching].copy()
    members['direction_sign']=members.direction.map({'negative':-1,'positive':1})
    frozen={(sp,'P_S'):members[members.species.eq(sp)] for sp in AXES}
    ctx.table(members,'process/new_primary_candidate_roster.tsv')
    all_effects=[];endpoint_rows=[];classes=[];profiles=[];primary={};coverage=[]
    native=ctx.read(ctx.out/'library_state_values.tsv.gz');native=native[native.setting.eq('supported50_p08')].rename(columns={'value':'APA'})
    variants={setting:f[f.axis_id.isin(AXES.values())] for setting,f in states.groupby('setting')}
    variants['native_PAS_coordinates']=native
    effects_columns=['species','study_id','route','setting','gene_id','gene_symbol','variant','source_kind','state','source_WPP','target_WPP','delta_WPP']
    for variant,frame in variants.items():
        frame=frame.merge(members[['species','gene_id']],on=['species','gene_id'],validate='many_to_one')
        valid=frame[frame.state.eq('S')&frame.time.astype(str).isin(['D0','source','Stage_I_day_0'])]
        effects=math.aggregate_fixed_reference(frame,frozen,variant,'native_axis' if variant=='native_PAS_coordinates' else 'shared_axis') if len(valid) else pd.DataFrame(columns=effects_columns)
        all_effects.append(effects)
    # Admission always comes from this run's baseline shared-axis endpoints.
    effect=pd.concat(all_effects,ignore_index=True);ctx.table(effect,'process/state_effects.tsv.gz')
    for sp in AXES:
        genes=frozen[sp,'P_S']
        if genes.empty:continue
        base=effect[effect.variant.eq('shared_p08')]
        primary[sp]=math.endpoint_table(base,genes,sp,.02,'new_shared_axis_0.02')
        for variant in variants:
            if sp=='human' and variant.startswith('exclude_'):continue
            current=effect[effect.variant.eq(variant)]
            for endpoint_cut,intermediate_cut,scope in [(.02,.02,'requalified'),(.02,.05,'fixed_primary_admission'),(.05,.05,'stronger_endpoint')]:
                ep=primary[sp].copy() if scope=='fixed_primary_admission' else math.endpoint_table(current,genes,sp,endpoint_cut,'current_variant')
                label='|'.join([variant,scope,str(endpoint_cut),str(intermediate_cut)])
                ep['analysis']=label;ep['variant']=variant;ep['scope']=scope
                pair,profile=math.process_tables(current,genes,sp,ep,intermediate_cut,label)
                for f in [pair,profile]:f['variant']=variant;f['scope']=scope
                endpoint_rows.append(ep);classes.append(pair);profiles.append(profile)
                eligible=pair[pair.endpoint_retained&pair.pair_measurable]
                required=['R','E'] if sp=='mouse' else ['R'];complete=eligible.groupby('gene_id').state.nunique().eq(len(required)).sum()
                coverage.append(dict(species=sp,variant=variant,scope=scope,endpoint_threshold=endpoint_cut,intermediate_threshold=intermediate_cut,
                    original_candidates=len(genes),endpoint_measurable=int(ep.endpoint_measurable.sum()),endpoint_retained=int(ep.endpoint_retained.sum()),
                    R=int(eligible.state.eq('R').sum()),E=int(eligible.state.eq('E').sum()) if sp=='mouse' else np.nan,
                    complete_required_states=int(complete)))
    endpoints=pd.concat(endpoint_rows,ignore_index=True) if endpoint_rows else pd.DataFrame(columns=['analysis','species','gene_id','variant','scope','endpoint_retained'])
    pair=pd.concat(classes,ignore_index=True) if classes else pd.DataFrame(columns=['analysis','species','gene_id','state','variant','scope','endpoint_retained','pair_measurable','pair_class'])
    profile=pd.concat(profiles,ignore_index=True) if profiles else pd.DataFrame(columns=['analysis','species','gene_id','variant','scope','profile_class'])
    ctx.table(endpoints,'process/endpoint_eligibility.tsv');ctx.table(pair,'process/same_series_process_classes.tsv.gz');ctx.table(profile,'process/profiles.tsv')
    ctx.table(coverage,'process/coverage_flow.tsv',columns=['species','variant','scope','endpoint_threshold','intermediate_threshold','original_candidates','endpoint_measurable','endpoint_retained','R','E','complete_required_states'])
    primary_pairs=pair[pair.variant.eq('shared_p08')&pair.scope.eq('requalified')]
    rows=[]
    for (sp,state),p in primary_pairs.groupby(['species','state']):
        cohort=p[p.endpoint_retained&p.pair_measurable][['gene_id','pair_class']].rename(columns={'pair_class':'primary_class'})
        for (variant,scope),q in pair[pair.species.eq(sp)&pair.state.eq(state)].groupby(['variant','scope']):
            z=cohort.merge(q,on='gene_id',how='left',validate='one_to_one')
            z['outcome']=np.select([~z.pair_measurable.fillna(False),~z.endpoint_retained.fillna(False),z.pair_class.eq(z.primary_class)],
                ['unavailable','endpoint_no_longer_passes','stable'],default='category_changed')
            rows.append(dict(species=sp,state=state,variant=variant,scope=scope,fixed_cohort=len(z),
                **{key:int(z.outcome.eq(key).sum()) for key in ['stable','category_changed','endpoint_no_longer_passes','unavailable']}))
    ctx.table(rows,'process/fixed_cohort_outcomes.tsv',columns=['species','state','variant','scope','fixed_cohort','stable','category_changed','endpoint_no_longer_passes','unavailable'])
    support=ctx.read(ctx.resource('state_expression_support'))
    ctx.table(support.groupby(['species','study_id','node']).agg(metacells=('metacell_id','size'),support_q95=('support_q95','sum'),support_q90=('support_q90','sum')).reset_index(),'process/frozen_expression_support_summary.tsv')
