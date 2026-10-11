"""Requalify intermediate states and all seven paired state contrasts."""
import argparse
import json
import shutil

import numpy as np
import pandas as pd
from scipy import stats

from downstream_parallel import ordered_map
from downstream_runtime import Context,Study,STUDIES,bh

PAIRS=[('R','S'),('E','S'),('P','S'),('E','R'),('P','E'),('P','M'),('M','S')]
SETTINGS=['supported50_p08','supported50_p09','expanded_primary_p08']
BASE=['species','study_id','route','setting','gene_id','gene_symbol']
STATE_TABLES=['library_state_values.tsv.gz','paired_library_effects.tsv.gz','global_contrast_eligibility.tsv.gz',
    'state_library_inventory.tsv','state_date_values.tsv.gz','state_means.tsv.gz','state_mean_contrasts.tsv.gz',
    'paired_gene_effects.tsv.gz','WOT_date_stratified_tests.tsv','all_gene_state_patterns.tsv.gz',
    'cross_study_patterns.tsv.gz','cross_study_summary.tsv','paired_cross_study_effects.tsv.gz','paired_cross_study_summary.tsv']


def qualified_state_groups(meta,purity):
    selected=meta[meta.purity.ge(purity)&meta.node.isin(list('SREPM'))]
    return {key:g for key,g in selected.groupby(['route','library_key','node']) if len(g)>=3 and g.n_cells.sum()>=30}


def extract_study(ctx,item):
    sid,setting=item
    states=[];pairs=[];eligibility=[];inventory=[]
    study=Study(ctx,sid);design=study.design[['route','library_key','sample_unit','replicate_status','time']].drop_duplicates()
    assert not design.duplicated(['route','library_key']).any();design=design.set_index(['route','library_key']).to_dict('index')
    purity=.9 if setting.endswith('p09') else .8
    groups=qualified_state_groups(study.meta,purity)
    assert all((route,lib) in design for route,lib,node in groups),'Qualified named states need a frozen sampling-design row'
    for key,g in groups.items():inventory.append(dict(study_id=sid,setting=setting,route=key[0],library_key=key[1],state=key[2],n_mc=len(g),n_cells=int(g.n_cells.sum())))
    for gene in study.genes:
        axis=study.axis(gene,setting)
        if len(axis)<2:continue
        data=study.gene_counts(axis);library={};route_totals={}
        for (route,lib,node),group in groups.items():
            q=study.state(data,route,lib,node,setting)
            if not q['endpoint_eligible']:continue
            d=design[route,lib];key=dict(species=axis.species.iloc[0],study_id=sid,route=route,setting=setting,gene_id=gene,gene_symbol=axis.primary_gene_symbol.iloc[0])
            row=key|dict(library_key=lib,sample_unit=d['sample_unit'],replicate_status=d['replicate_status'],time=d['time'],state=node,
                value=float(q['usage']@axis.position),n_mc=q['n_mc'],n_cells=q['n_cells'],pas_counts=int(q['site_counts'].sum()),available_pas=len(axis))
            states.append(row);library[route,lib,node]=(row,q['site_counts'],q['site_mc'])
            total=route_totals.setdefault((route,node),[np.zeros(len(axis),np.int64),np.zeros(len(axis),np.int64)])
            total[0]+=q['site_counts'];total[1]+=q['site_mc']
        for route in {k[0] for k in route_totals}:
            for a,b in PAIRS:
                if (route,a) not in route_totals or (route,b) not in route_totals:continue
                ca,ma=route_totals[route,a];cb,mb=route_totals[route,b]
                eligibility.append(dict(species=axis.species.iloc[0],study_id=sid,route=route,setting=setting,gene_id=gene,
                    contrast=a+'_'+b,n_joint_pas=int(((ca+cb>=10)&(ma+mb>=2)).sum())))
        for route,lib in {(k[0],k[1]) for k in library}:
            for a,b in PAIRS:
                if (route,lib,a) not in library or (route,lib,b) not in library:continue
                da,ca,ma=library[route,lib,a];db,cb,mb=library[route,lib,b];n=int(((ca+cb>=10)&(ma+mb>=2)).sum())
                if n<2:continue
                assert da['time']==db['time']
                pairs.append({k:da[k] for k in BASE+['library_key','sample_unit','replicate_status','time']}|dict(contrast=a+'_'+b,
                    delta=da['value']-db['value'],n_joint_pas=n,n_mc_a=da['n_mc'],n_mc_b=db['n_mc']))
    print(sid,setting,len(states),len(pairs),flush=True)
    return (states,pairs,eligibility,inventory)


def extract(ctx):
    states=[];pairs=[];eligibility=[];inventory=[]
    for result in ordered_map(ctx,extract_study,[(sid,setting) for sid in STUDIES for setting in SETTINGS]):
        for target,values in zip([states,pairs,eligibility,inventory],result):target.extend(values)
    s=pd.DataFrame(states,columns=BASE+['library_key','sample_unit','replicate_status','time','state','value','n_mc','n_cells','pas_counts','available_pas'])
    p=pd.DataFrame(pairs,columns=BASE+['library_key','sample_unit','replicate_status','time','contrast','delta','n_joint_pas','n_mc_a','n_mc_b'])
    e=pd.DataFrame(eligibility,columns=BASE[:-1]+['contrast','n_joint_pas'])
    ctx.table(s,'library_state_values.tsv.gz');ctx.table(p,'paired_library_effects.tsv.gz');ctx.table(e,'global_contrast_eligibility.tsv.gz');ctx.table(inventory,'state_library_inventory.tsv')
    return s,p,e


def analyze(ctx,s,p,e):
    k=BASE+['state'];u=s.groupby(k+['time','sample_unit']).value.mean().reset_index();daily=u.groupby(k+['time']).value.mean().reset_index()
    state=daily.groupby(k).value.agg(value='mean',n_dates='size',date_min='min',date_max='max').reset_index()
    state=state.merge(u.groupby(k).size().rename('n_units').reset_index(),on=k)
    ctx.table(daily,'state_date_values.tsv.gz');ctx.table(state,'state_means.tsv.gz')
    wide=state.pivot(index=BASE,columns='state',values='value').reindex(columns=list('SREPM'));effects=[]
    for a,b in PAIRS:
        z=wide[[a,b]].dropna().reset_index();z['contrast']=a+'_'+b;z['delta']=z[a]-z[b];z=z.drop(columns=[a,b])
        z=z.merge(e,on=BASE[:-1]+['contrast'],validate='one_to_one');effects.append(z[z.n_joint_pas>=2])
    effect=pd.concat(effects,ignore_index=True);ctx.table(effect,'state_mean_contrasts.tsv.gz')
    kp=BASE+['contrast'];pu=p.groupby(kp+['sample_unit','time']).delta.mean().reset_index();pdaily=pu.groupby(kp+['time']).delta.mean().reset_index()
    pg=pdaily.groupby(kp).delta.agg(delta='mean',n_dates='size',date_min='min',date_max='max').reset_index();pg=pg.merge(pu.groupby(kp).size().rename('n_units').reset_index(),on=kp)
    pg['leave_date_min']=np.where(pg.n_dates>1,(pg.delta*pg.n_dates-pg.date_max)/(pg.n_dates-1),np.nan)
    pg['leave_date_max']=np.where(pg.n_dates>1,(pg.delta*pg.n_dates-pg.date_min)/(pg.n_dates-1),np.nan);ctx.table(pg,'paired_gene_effects.tsv.gz')
    tests=[]
    for key,z in pu[pu.route.eq('WOT')&pu.setting.eq('supported50_p08')].groupby(kp):
        v=z.groupby('time').delta.agg(['mean','size']);df=len(z)-len(v);estimate=v['mean'].mean();se=np.nan;pval=1.;status='insufficient_units_or_df'
        if len(z)>=6 and df>=3:
            residual=z.delta-z.time.map(v['mean']);variance=(residual**2).sum()/df;se=np.sqrt(variance*(1/v['size']).sum()/len(v)**2)
            if se>0:pval=2*stats.t.sf(abs(estimate/se),df);status='tested'
            elif estimate==0:pval=1.;status='zero_effect_zero_variance'
            else:status='nonzero_effect_zero_variance_descriptive'
        tests.append(dict(zip(kp,key))|dict(effect=estimate,se=se,p_value=pval,n_units=len(z),n_dates=len(v),df=df,status=status))
    t=pd.DataFrame(tests,columns=kp+['effect','se','p_value','n_units','n_dates','df','status'])
    t['q_within_contrast']=t.groupby('contrast').p_value.transform(bh);t['q_all_requested_contrasts']=bh(t.p_value);ctx.table(t,'WOT_date_stratified_tests.tsv')
    patterns=[]
    for key,q in state.groupby(['species','study_id','route','setting']):
        vals=q.pivot(index=['gene_id','gene_symbol'],columns='state',values='value').reindex(columns=list('SREPM'))
        ev=effect
        for col,value in zip(['species','study_id','route','setting'],key):ev=ev[ev[col].eq(value)]
        valid={c:set(z.gene_id) for c,z in ev.groupby('contrast')}
        for inter in ['R','E']:
            z=vals[['S',inter,'P']].dropna();z=z[z.index.get_level_values('gene_id').isin(valid.get('P_S',set())&valid.get(inter+'_S',set()))]
            for threshold,small in [(.02,.01),(.05,.02)]:
                for (gid,symbol),r in z.iterrows():
                    ps=r.P-r.S;i_s=r[inter]-r.S;i_p=r[inter]-r.P;endpoint=abs(ps)>=threshold
                    patterns.append(dict(zip(['species','study_id','route','setting'],key))|dict(gene_id=gid,gene_symbol=symbol,intermediate=inter,
                        threshold=threshold,small_endpoint_cutoff=small,P_S=ps,I_S=i_s,I_P=i_p,endpoint=endpoint,
                        intermediate_excursion=abs(ps)<small and abs(i_s)>=threshold and abs(i_p)>=threshold and i_s*i_p>0,
                        intermediate_already_shifted=endpoint and ps*i_s>0 and abs(i_s)>=threshold,
                        endpoint_concentrated=endpoint and abs(i_s)<small and abs(i_p)>=threshold,intermediate_fraction=i_s/ps if endpoint else np.nan))
    b=pd.DataFrame(patterns,columns=BASE+['intermediate','threshold','small_endpoint_cutoff','P_S','I_S','I_P','endpoint','intermediate_excursion','intermediate_already_shifted','endpoint_concentrated','intermediate_fraction'])
    ctx.table(b,'all_gene_state_patterns.tsv.gz');cross=[];summary=[]
    for name,sp,r1,r2,inter in [('mouse_R','mouse','CIP','WOT','R'),('mouse_E','mouse','CIP','WOT','E'),('human_R','human','GSE118258','178325_0618','R'),('human_E_mixed','human','GSE118258','GSE150311','E')]:
        for setting in SETTINGS:
            for threshold in [.02,.05]:
                q=b[b.species.eq(sp)&b.intermediate.eq(inter)&b.setting.eq(setting)&b.threshold.eq(threshold)]
                w=q[q.route.eq(r1)].merge(q[q.route.eq(r2)],on='gene_id',suffixes=('_1','_2'),validate='one_to_one')
                same=w.endpoint_1.astype(bool)&w.endpoint_2.astype(bool)&(w.P_S_1*w.P_S_2>0)
                exc=w.intermediate_excursion_1.astype(bool)&w.intermediate_excursion_2.astype(bool)&(w.I_S_1*w.I_S_2>0)
                w['cohort']=name;w['setting']=setting;w['threshold']=threshold;w['repeated_endpoint']=same;w['repeated_excursion']=exc
                w['endpoint_both_early']=same&w.intermediate_already_shifted_1.astype(bool)&w.intermediate_already_shifted_2.astype(bool)
                w['endpoint_both_late']=same&w.endpoint_concentrated_1.astype(bool)&w.endpoint_concentrated_2.astype(bool)
                w['early_late_mismatch']=same&((w.intermediate_already_shifted_1.astype(bool)&w.endpoint_concentrated_2.astype(bool))|(w.intermediate_already_shifted_2.astype(bool)&w.endpoint_concentrated_1.astype(bool)))
                cross.append(w);summary.append(dict(cohort=name,setting=setting,threshold=threshold,background=len(w),repeated_endpoint=int(same.sum()),repeated_intermediate_excursion=int(exc.sum()),
                    endpoint_both_intermediate_shifted=int(w.endpoint_both_early.sum()),endpoint_both_endpoint_concentrated=int(w.endpoint_both_late.sum()),endpoint_intermediate_pattern_mismatch=int(w.early_late_mismatch.sum())))
    ctx.table(pd.concat(cross,ignore_index=True),'cross_study_patterns.tsv.gz');ctx.table(summary,'cross_study_summary.tsv')
    paired=[];psummary=[]
    for sp,r1,r2 in [('mouse','CIP','WOT'),('human','GSE118258','GSE150311'),('human','GSE118258','GSE221739')]:
        for setting in SETTINGS:
            for a,b in PAIRS:
                co=a+'_'+b;q=pg[pg.species.eq(sp)&pg.setting.eq(setting)&pg.contrast.eq(co)]
                w=q[q.route.eq(r1)].merge(q[q.route.eq(r2)],on='gene_id',suffixes=('_1','_2'),validate='one_to_one')
                for threshold in [.02,.05]:
                    same=w.delta_1.abs().ge(threshold)&w.delta_2.abs().ge(threshold)&(w.delta_1*w.delta_2>0)
                    paired.append(w.assign(repeated=same,threshold=threshold,setting=setting,contrast=co,species=sp,comparison=r1+'__'+r2))
                    psummary.append(dict(species=sp,comparison=r1+'__'+r2,setting=setting,contrast=co,threshold=threshold,background=len(w),repeated=int(same.sum())))
    ctx.table(pd.concat(paired,ignore_index=True),'paired_cross_study_effects.tsv.gz');ctx.table(psummary,'paired_cross_study_summary.tsv')


def main(ctx):
    if 'intermediate_resume' in ctx.config.get('resources',{}):
        manifest=json.loads(ctx.resource('intermediate_resume').read_text())
        assert manifest['run_id']==ctx.run_id and manifest['input_fingerprint']==ctx.config['input_fingerprint']
        assert {r['name'] for r in manifest['files']}==set(STATE_TABLES)
        for row in manifest['files']:
            source=ctx.register(row['source']);assert ctx.inputs[str(source)]==row['sha256']
            destination=ctx.out/row['name'];assert not destination.exists();shutil.copy2(source,destination)
        for path,value in manifest['upstream_inputs'].items():
            source=ctx.register(path);assert ctx.inputs[str(source)]==value
    else:
        states,pairs,eligibility=extract(ctx);analyze(ctx,states,pairs,eligibility)
    from downstream_shared_states import extract as shared_extract,compare as shared_compare
    shared=shared_extract(ctx);shared_compare(ctx,shared)
    from downstream_process_classes import process_classes
    process_classes(ctx,shared)
    ctx.finish(dict(new_gene_qualified_states=True,new_joint_support=True,fixed_physical_libraries=True,dates_and_units_balanced=True,
        all_gene_patterns_reselected=True,all_planned_comparisons_retained=True,new_shared_axes=True,
        shared_contributor_GEX=True,conditional_cycle_expression_depth_models=True,local_state_exclusions_recomputed=True,
        common_window_triangle_identity=True,new_endpoint_admission_and_fixed_cohort_classes=True,
        frozen_expression_support_applied_to_new_PAS_counts=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();main(Context(a.config,'intermediate',a.out))
