"""Rebuild study events, conditional nulls and nested-control families."""
import argparse
import hashlib
import json
from itertools import product

import numpy as np
import pandas as pd
from scipy.stats import beta

from downstream_runtime import Context, bh, seed

SOURCE='notes/2026-09-30-multimetric-statistics/nulls/analyze.py'
THRESHOLDS={'WPP':[.02,.05],'RWUI':[.02,.05],'MPRO':[.1,.2]}


def numerical(ctx):
    return ctx.functions(SOURCE,['qbin','collapse','stratify','shuffle','summarize','cross_species'],
        dict(THRESHOLDS=THRESHOLDS,HUMAN=['GSE118258','GSE178325'],MOUSE=['GSE103221','GSE115943']))


def events_from_routes(routes,math):
    assert not routes.duplicated(['species','study_id','route','gene_id','contrast']).any()
    events=math.collapse(routes)
    patterns=routes.groupby(['species','contrast','study_id','gene_id']).route.agg(lambda x:';'.join(sorted(x)))
    keys=pd.MultiIndex.from_frame(events[['species','contrast','unit','gene_id']])
    events['series_pattern']=patterns.reindex(keys).to_numpy()
    assert events.series_pattern.notna().all()
    return events


def recurrence_null(events,math,n=10000,exact_series=True,tag='primary',thresholds=None):
    summaries=[];members=[];strata=[];diagnostics=[];archive={};loso=[]
    for metric,cuts in (THRESHOLDS if thresholds is None else thresholds).items():
        for species,contrast,threshold in product(['human','mouse'],['P_S','P_M'],cuts):
            key=dict(metric=metric,species=species,contrast=contrast,threshold=threshold,series_matching=exact_series,variant=tag)
            f=events[events.metric.eq(metric)&events.species.eq(species)&events.contrast.eq(contrast)&events.threshold.eq(threshold)]
            if f.unit.nunique()<2:
                summaries.append(key|dict(observed=0,empirical_p=1.,status='fewer_than_two_studies',n_genes=f.gene_id.nunique()))
                continue
            genes=sorted(f.gene_id.unique());index=pd.Index(genes)
            positive=np.zeros((n,len(genes)),np.int8);negative=np.zeros_like(positive)
            op=np.zeros(len(genes),int);om=np.zeros(len(genes),int)
            rng=np.random.default_rng(seed(tag,metric,species,contrast,threshold,exact_series))
            labels=f.pivot(index='gene_id',columns='unit',values='event').reindex(genes)
            for unit,part in f.groupby('unit'):
                sf=math.stratify(part,metric,threshold)
                if exact_series:sf['stratum']=sf.stratum+'|'+sf.series_pattern
                draws,prob,diag=math.shuffle(sf,rng,n)
                ix=index.get_indexer(sf.gene_id)
                positive[:,ix]+=draws==1;negative[:,ix]+=draws==-1
                op[ix]+=sf.event.eq(1).to_numpy();om[ix]+=sf.event.eq(-1).to_numpy()
                strata.append(sf.assign(**{k:v for k,v in key.items() if k not in sf}))
                diagnostics.append(diag.assign(unit=unit,**key))
            support=np.maximum(op,om);null=(np.maximum(positive,negative)>=2).sum(1)
            record=key|math.summarize((support>=2).sum(),null)|dict(status='tested',n_genes=len(genes),n_studies=f.unit.nunique())
            for n_support in range(1,f.unit.nunique()+1):
                exact=(np.maximum(positive,negative)==n_support).sum(1)
                record.update({f'exact_{n_support}_observed':int((support==n_support).sum()),f'exact_{n_support}_null_median':float(np.median(exact)),
                    f'exact_{n_support}_null_p025':float(np.quantile(exact,.025)),f'exact_{n_support}_null_p975':float(np.quantile(exact,.975))})
            summaries.append(record)
            archive['|'.join(map(str,key.values()))]=null
            symbols=f.drop_duplicates('gene_id').set_index('gene_id').gene_symbol
            for i in np.flatnonzero(support>=2):
                members.append(key|dict(gene_id=genes[i],gene_symbol=symbols.loc[genes[i]],
                    eligible_units=int(labels.iloc[i].notna().sum()),positive_support=op[i],negative_support=om[i],
                    direction='unsigned' if metric=='MPRO' else 'positive' if op[i]>=om[i] else 'negative'))
            for unit in labels.columns:
                sub=labels.drop(columns=unit)
                retained=np.maximum(sub.eq(1).sum(axis=1),sub.eq(-1).sum(axis=1)).ge(2)
                loso.append(key|dict(omitted_study=unit,remaining_studies=sub.shape[1],recurrent_genes=int(retained.sum()),
                    retained_original=int((retained.to_numpy()&(support>=2)).sum())))
    summary=pd.DataFrame(summaries)
    summary['q_metric_8']=summary.groupby('metric').empirical_p.transform(bh)
    summary['q_all_24']=bh(summary.empirical_p)
    member_columns=['metric','species','contrast','threshold','series_matching','variant','gene_id',
                    'gene_symbol','eligible_units','positive_support','negative_support','direction']
    return summary,pd.DataFrame(members,columns=member_columns),strata,diagnostics,archive,pd.DataFrame(loso)


def omit_study(events,species,unit):
    """Rebuild measurable-study counts; missing studies never become zero events."""
    out=events[events.species.eq(species)&events.unit.ne(unit)].copy()
    out['eligible_units']=out.groupby(['species','contrast','metric','threshold','gene_id']).unit.transform('nunique')
    return out


def leave_one_study_nulls(events,math,n=10000):
    summaries=[];strata=[];diagnostics=[];archive={}
    for species,source in events.groupby('species',sort=True):
        for unit in sorted(source.unit.unique()):
            contrasts=set(source[source.unit.eq(unit)].contrast)
            reduced=omit_study(events,species,unit)
            reduced=reduced[reduced.contrast.isin(contrasts)]
            for exact in [False,True]:
                s,m,st,d,a,l=recurrence_null(reduced,math,n=n,exact_series=exact,tag='omit_'+unit)
                s=s[s.species.eq(species)].drop(columns=['q_metric_8','q_all_24']).copy()
                s.loc[~s.contrast.isin(contrasts),'status']='omitted_study_not_evaluable'
                s['omitted_study']=unit
                # Four planned contrast/threshold tests per metric within each omission.
                assert s.groupby('metric').size().eq(4).all() and len(s)==12
                s['q_metric_4']=s.groupby('metric').empirical_p.transform(bh)
                s['q_all_12']=bh(s.empirical_p)
                summaries.append(s)
                strata.extend(f.assign(omitted_study=unit) for f in st)
                diagnostics.extend(f.assign(omitted_study=unit) for f in d)
                assert not set(archive)&set(a);archive.update(a)
    return pd.concat(summaries,ignore_index=True),strata,diagnostics,archive


def recurrence(ctx):
    math=numerical(ctx)
    routes=ctx.read(ctx.analysis/'state_vectors/route_effects.primary.tsv.gz')
    events=events_from_routes(routes,math);ctx.table(events,'study_event_covariates.tsv.gz')
    outputs=[];members=[];strata=[];diags=[];loso=[];archive={}
    for exact in [False,True]:
        s,m,st,d,a,l=recurrence_null(events,math,exact_series=exact)
        outputs.append(s);members.append(m);strata.extend(st);diags.extend(d);archive.update(a);loso.append(l)
    ctx.table(pd.concat(outputs),'recurrence_summary.tsv')
    ctx.table(pd.concat(members),'recurrence_members.tsv.gz',columns=['metric','species','contrast','threshold','series_matching','variant','gene_id','gene_symbol','eligible_units','positive_support','negative_support','direction'])
    ctx.table(pd.concat(strata) if strata else pd.DataFrame(columns=['unit','gene_id','stratum']),'recurrence_strata.tsv.gz')
    ctx.table(pd.concat(diags) if diags else pd.DataFrame(columns=['unit','stratum','genes']),'null_mobility.tsv.gz')
    ctx.table(pd.concat(loso),'leave_one_study_out.tsv',columns=['metric','species','contrast','threshold',
        'series_matching','omitted_study','remaining_studies','recurrent_genes','retained_original'])
    np.savez_compressed(ctx.out/'null_draws.npz',**archive)
    ls,lst,ld,la=leave_one_study_nulls(events,math)
    ctx.table(ls,'leave_one_study_out_nulls.tsv')
    ctx.table(pd.concat(lst) if lst else pd.DataFrame(columns=['omitted_study','unit','gene_id','stratum']),
        'leave_one_study_out_strata.tsv.gz')
    ctx.table(pd.concat(ld) if ld else pd.DataFrame(columns=['omitted_study','unit','stratum','genes']),
        'leave_one_study_out_mobility.tsv.gz')
    np.savez_compressed(ctx.out/'leave_one_study_out_null_draws.npz',**la)
    ctx.finish(dict(new_study_event_universe=True,all_evaluable_series_agree=True,
        one_vote_per_GEO=True,exact_series_control=True,planned_BH_families_rebuilt=True,
        omission_eligibility_strata_and_nulls_rebuilt=True),
        primary='series_matching=True; earlier number-only matching retained as sensitivity',permutations=10000)


def cross_sharing(ctx):
    math=numerical(ctx);events=ctx.read(ctx.analysis/'recurrence/study_event_covariates.tsv.gz')
    ortholog=ctx.read(ctx.resource('ortholog'))
    assert not ortholog.human_gene_id.duplicated().any() and not ortholog.mouse_gene_id.duplicated().any()
    diagnostics,archive=cross_species_current(ctx,math,events,ortholog,10000,20261003)
    ctx.table(pd.concat(diagnostics) if diagnostics else pd.DataFrame(columns=['species','metric','human_study','mouse_study','threshold','stratum','genes']),'null_mobility.tsv.gz')
    np.savez_compressed(ctx.out/'null_draws.npz',**archive)
    ctx.finish(dict(new_ortholog_backgrounds=True,native_species_metrics=True,
        distinct_8_24_56_families=True,documented_source_scope=True),permutations=10000)


def cross_strata(math,frame,metric,threshold,model):
    sf=math.stratify(frame,metric,threshold,pair=True)
    if model=='Range10':
        assert metric=='WPP','The prespecified range sensitivity applies to DWUI'
        for side in ['min','max']:
            value=sf['baseline_'+side].to_numpy(float)
            assert np.isfinite(value).all() and ((value>=-1e-12)&(value<=1+1e-12)).all()
            sf['baseline_'+side+'_bin']=np.floor(value*10).clip(0,9).astype(int)
        fields=['eligible_units','n_routes','mean_GEX_bin','APA_coverage_bin','PAS_number_bin',
            'baseline_min_bin','baseline_max_bin','can_negative','can_positive']
        sf['stratum']=sf[fields].astype(str).agg('|'.join,axis=1)
    else:assert model=='Q5'
    return sf


def cross_species_current(ctx, math, events, ortholog, n, seed,model='Q5'):
    HUMAN=['GSE118258','GSE178325'];MOUSE=['GSE103221','GSE115943']
    stratify=math.stratify;shuffle=math.shuffle;summarize=math.summarize;save=ctx.table
    summaries, backgrounds, members, strata, diagnostics, archive = [], [], [], [], [], {}
    for metric, thresholds in (THRESHOLDS if model=='Q5' else {'WPP':THRESHOLDS['WPP']}).items():
        for h in HUMAN:
            for m in MOUSE:
                for threshold in thresholds:
                    context = dict(metric=metric, human_study=h, mouse_study=m, threshold=threshold,model=model)
                    selected = events[events.metric.eq(metric) & events.threshold.eq(threshold) & events.contrast.eq("P_S")]
                    hs = selected[selected.unit.eq(h)]
                    ms = selected[selected.unit.eq(m)]
                    f = ortholog.merge(hs.add_prefix("h_"), left_on="human_gene_id", right_on="h_gene_id", validate="one_to_one")
                    f = f.merge(ms.add_prefix("m_"), left_on="mouse_gene_id", right_on="m_gene_id", validate="one_to_one")
                    f = f.sort_values("human_gene_id").reset_index(drop=True)
                    key = "sharing__" + "__".join(map(str, context.values()))
                    run_seed = seed + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
                    rng = np.random.default_rng(run_seed)
                    draws, probabilities = [], []
                    for species, prefix in ([("human", "h_"), ("mouse", "m_")] if len(f) else []):
                        sub = f[[column for column in f if column.startswith(prefix)]].rename(columns=lambda x: x[len(prefix):])
                        sf = cross_strata(math,sub,metric,threshold,model)
                        shuffled, prob, diag = shuffle(sf, rng, n)
                        draws.append(shuffled)
                        probabilities.append(prob)
                        strata.append(sf.assign(human_gene_id=f.human_gene_id, mouse_gene_id=f.mouse_gene_id, analysis="sharing", **context))
                        diagnostics.append(diag.assign(analysis="sharing", species=species, **context))
                    if not len(f):
                        draws=[np.empty((n,0),int)]*2;probabilities=[np.empty((0,3))]*2
                    tests = ["shared_unsigned"] if metric == "MPRO" else ["same_direction", "shared_negative", "shared_positive"]
                    backgrounds.append(f.assign(**context))
                    shared = f.h_event.eq(f.m_event) & f.h_event.ne(0)
                    members.append(f[shared].assign(**context))
                    for test in tests:
                        values = [1] if test == "shared_unsigned" or test == "shared_positive" else [-1] if test == "shared_negative" else [-1, 1]
                        observed = sum(((f.h_event == value) & (f.m_event == value)).sum() for value in values)
                        null = sum(((draws[0] == value) & (draws[1] == value)).sum(axis=1) for value in values)
                        analytic = sum((probabilities[0][:, value + 1] * probabilities[1][:, value + 1]).sum() for value in values)
                        assert abs(null.mean() - analytic) <= 6 * null.std(ddof=1) / np.sqrt(n) + 1e-10
                        summaries.append(context | {"test": test, "common_orthologs": len(f), "analytic_null_mean": analytic,
                                                    "seed": run_seed, "permutations": n, "status": "tested" if len(f) else "no_common_measurable_orthologs"} | summarize(observed, null))
                        archive[key + "__" + test] = null
                    print(key, len(f), summaries[-len(tests)], flush=True)
    summary = pd.DataFrame(summaries)
    summary["q_metric_all"] = summary.groupby("metric").empirical_p.transform(bh)
    for metric, f in summary.groupby("metric"):
        primary = f[f.test.isin(["same_direction", "shared_unsigned"])]
        assert len(primary) == 8
        summary.loc[primary.index, "q_primary_8"] = bh(primary.empirical_p)
        directional = f[~f.test.isin(["same_direction", "shared_unsigned"])]
        if len(directional):
            assert len(directional) == 16
            summary.loc[directional.index, "q_directional_16"] = bh(directional.empirical_p)
    primary = summary[summary.test.isin(["same_direction", "shared_unsigned"])]
    if model=='Q5':
        assert len(primary)==24 and len(summary)==56
        summary.loc[primary.index, "q_joint_primary_24"] = bh(primary.empirical_p)
        summary["q_joint_all_56"] = bh(summary.empirical_p)
    else:
        assert len(primary)==8 and len(summary)==24
        summary['q_DWUI_all_24']=bh(summary.empirical_p)
    save(summary, "sharing_summary.tsv")
    save(pd.concat(backgrounds, ignore_index=True), "sharing_backgrounds.tsv.gz")
    save(pd.concat(members, ignore_index=True), "sharing_members.tsv.gz")
    save(pd.concat(strata, ignore_index=True) if strata else pd.DataFrame(columns=["species","gene_id","stratum"]), "sharing_strata.tsv.gz")
    return diagnostics, archive


def controls(ctx):
    math=numerical(ctx)
    common=ctx.functions('notes/2026-10-02-null-model-controls/common.py',
        ['probability','rate','groups'],dict(beta=beta))
    sim=ctx.functions('notes/2026-10-02-null-model-controls/recurrence.py',['simulate'])
    base=['eligible_units','series_pattern','can_negative','can_positive'];meas=['mean_GEX_bin','APA_coverage_bin']
    models={'L0_feasible_series':base,'L1_measurement':base+meas,'L2_baseline':base+meas+['baseline_bin'],
        'L2_PAS_first':base+meas+['PAS_number_bin'],'L3_full':base+meas+['baseline_bin','PAS_number_bin']}
    events=ctx.read(ctx.analysis/'recurrence/study_event_covariates.tsv.gz')
    rows=[];strata=[];power=[];frozen=[];archive={}
    for key,f in events[events.metric.eq('WPP')].groupby(['species','contrast','threshold']):
        context=dict(zip(['species','contrast','threshold'],key));genes=pd.Index(sorted(f.gene_id.unique()))
        lab=f.pivot(index='gene_id',columns='unit',values='event')
        observed=int(np.maximum(lab.eq(1).sum(1),lab.eq(-1).sum(1)).ge(2).sum())
        prepared=[]
        for unit,part in f.groupby('unit'):
            sf=math.stratify(part,'WPP',key[2]);prepared.append((sf,genes.get_indexer(sf.gene_id)));frozen.append(sf)
        nulls={};units_by_model={}
        for model,cols in models.items():
            units=[];fixed=0;changed=0.;nstrata=0
            for sf,ix in prepared:
                cells=common.groups(sf,cols);units.append((sf,ix,cells));nstrata+=len(cells)
                for j,ii in enumerate(cells):
                    freq=sf.event.iloc[ii].value_counts().reindex([-1,0,1],fill_value=0).to_numpy()
                    moving=1-(freq/len(ii))@(freq/len(ii));changed+=len(ii)*moving
                    fixed+=len(ii)*(np.count_nonzero(freq)==1)
                    strata.append(context|dict(model=model,unit=sf.unit.iloc[0],stratum=j,n=len(ii),negative=freq[0],zero=freq[1],positive=freq[2],expected_changed=len(ii)*moving))
            units_by_model[model]=units
            null=sim.simulate(units,len(genes),10000,np.random.default_rng(seed('recurrence_control',*key,model)))
            nulls[model]=null;archive['|'.join(map(str,(*key,model)))]=null
            rows.append(context|dict(model=model,n_genes=len(genes),n_study_genes=len(f),n_studies=f.unit.nunique(),
                observed=observed,expected=null.mean(),null_sd=null.std(ddof=1),p=float(common.probability(null,observed)),
                strata=nstrata,frozen_label_fraction=fixed/len(f),expected_changed_fraction=changed/len(f)))
        for strength in [0,.5,1,2]:
            draw=sim.simulate(units_by_model['L3_full'],len(genes),1000,np.random.default_rng(seed('recurrence_inject',*key,strength)),strength)
            archive['|'.join(map(str,(*key,'inject',strength)))]=draw
            for model,reference in nulls.items():
                for alpha in [.05,.05/8]:
                    power.append(context|dict(generator='L3_full',tested_model=model,shared_logweight_sd=strength,alpha=alpha,
                        simulated_mean=draw.mean(),**common.rate(common.probability(reference,draw)<=alpha)))
    summary=pd.DataFrame(rows);summary['q_eight']=summary.groupby('model').p.transform(bh);summary['q_all_40']=bh(summary.p)
    ctx.table(summary,'recurrence_ladder.tsv');ctx.table(strata,'recurrence_strata.tsv.gz')
    ctx.table(pd.concat(frozen),'recurrence_frozen_covariates.tsv.gz');ctx.table(power,'recurrence_power.tsv')
    np.savez_compressed(ctx.out/'recurrence_draws.npz',**archive)
    series_unit_sensitivity(ctx,math)
    ctx.finish(dict(new_conditional_covariates=True,exact_series_identity=True,event_margins_preserved=True,
        signal_injection_calibrated=True,fresh_BH=True),reference_draws=10000,injection_draws=1000)


def series_unit_sensitivity(ctx,math):
    routes=ctx.read(ctx.analysis/'state_vectors/route_effects.primary.tsv.gz')
    routes=routes[routes.species.eq('human')&routes.contrast.eq('P_S')].copy()
    mask=routes.study_id.eq('GSE178325')
    routes.loc[mask,'study_id']=routes.loc[mask,'study_id']+'_'+routes.loc[mask,'route']
    events=events_from_routes(routes,math)
    summary,members,strata,diags,archive,loso=recurrence_null(events,math,exact_series=True,
        tag='Guan_series_as_separate_pseudounits',thresholds={'WPP':[.02,.05]})
    summary['scope']='sensitivity only: Guan series treated as separate pseudounits, not independent studies'
    ctx.table(summary,'guan_series_unit_sensitivity.tsv')
    ctx.table(members,'guan_series_unit_members.tsv.gz')
    np.savez_compressed(ctx.out/'guan_series_unit_nulls.npz',**archive)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--node',choices=['recurrence','cross_sharing','recurrence_controls'],required=True)
    a=p.parse_args();{'recurrence':recurrence,'cross_sharing':cross_sharing,'recurrence_controls':controls}[a.node](Context(a.config,a.node,a.out))
