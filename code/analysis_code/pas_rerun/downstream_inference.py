"""Fresh-cohort biological-unit and conditional-count inference."""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

import numpy as np
import pandas as pd
from scipy.stats import t, beta, fisher_exact, random_table

from downstream_runtime import Context, Study, bh, seed
from downstream_replicate_diagnostics import confidence_intervals

BASE='notes/2026-09-30-multimetric-statistics/'
METRICS=['DWUI','RWUI','PDUI','PPUI','MPRO']
FAMILY=['species','study_id','route','setting','contrast','design','library_key_a','library_key_b']


def load(ctx):
    p=ctx.read(ctx.analysis/'state_vectors/pair_effects.primary.tsv.gz')
    a=ctx.read(ctx.analysis/'state_vectors/PAS_axes.primary.tsv.gz')
    axes={(s,g):f for (s,g),f in a.groupby(['study_id','primary_gene_id'],sort=False)}
    return p,axes


def vectors(frame):
    return np.stack(frame.target_usage.map(json.loads))-np.stack(frame.reference_usage.map(json.loads))


def full_study_axes(study):
    # The exclusion denominator must include genes with no qualified endpoint pair.
    return {(study.id,gene):study.axis(gene) for gene in study.genes}


def source_states(study,axes):
    rows=[]
    groups=list(study.meta[study.meta.node.isin(['S','P'])&study.meta.purity.ge(.8)].groupby(['route','library_key','node']).groups)
    for (sid,gene),axis in axes.items():
        if sid!=study.id or len(axis)<2:continue
        data=study.gene_counts(axis)
        for route,library,node in groups:
            s=study.state(data,route,library,node)
            if not s['endpoint_eligible']:continue
            rows.append(dict(gene_id=gene,library_key=library,state=node,APA=float(s['usage']@axis.position),
                site_counts=json.dumps(s['site_counts'].tolist()),site_mc=json.dumps(s['site_mc'].tolist())))
    return pd.DataFrame(rows,columns=['gene_id','library_key','state','APA','site_counts','site_mc'])


def finish_tests(ctx,rows,name='gene_metric_tests.tsv.gz'):
    f=pd.DataFrame(rows,columns=['gene_id','gene_symbol','metric','effect','p','n_units','n_dates',
        'PAS_number','weighted_delta_usage','sign_patterns','exact_sign_enumeration',
        'minimum_two_sided_p','statistic','SE','df','CI95_low','CI95_high','date_min','date_max',
        'leave_date_min','leave_date_max'])
    f=confidence_intervals(f)
    f['q_BH_joint']=bh(f.p)
    f['q_BH_metric']=f.groupby('metric').p.transform(bh)
    f['joint_q05']=f.q_BH_joint.lt(.05)
    ctx.table(f,name)
    ctx.table(f.groupby('metric').agg(genes=('gene_id','size'),significant_joint=('joint_q05','sum'),
        minimum_p=('p','min'),minimum_joint_q=('q_BH_joint','min')).reset_index(),'detection_summary.tsv')
    return f


def mc_stability(tests):
    """Bound finite sign-sampling error while leaving the primary tests intact."""
    out=tests[['gene_id','gene_symbol','metric','p','q_BH_joint',
        'exact_sign_enumeration','sign_patterns']].copy()
    mc=~out.exact_sign_enumeration.astype(bool)
    n=out.loc[mc,'sign_patterns'].to_numpy(dtype=int)
    k=np.rint(out.loc[mc,'p'].to_numpy()*(n+1)-1).astype(int)
    assert ((k>=0)&(k<=n)).all()
    assert np.allclose(out.loc[mc,'p'],(k+1)/(n+1),atol=1e-14,rtol=0)
    alpha=.05/int(mc.sum()) if mc.any() else .05
    lower=np.zeros(len(k));upper=np.ones(len(k))
    low=k>0;high=k<n
    lower[low]=beta.ppf(alpha/2,k[low],n[low]-k[low]+1)
    upper[high]=beta.ppf(1-alpha/2,k[high]+1,n[high]-k[high])
    out['MC_extreme_patterns']=pd.Series(pd.NA,index=out.index,dtype='Int64')
    out.loc[mc,'MC_extreme_patterns']=k
    out['p_lower_simultaneous95']=out.p;out['p_upper_simultaneous95']=out.p
    out.loc[mc,'p_lower_simultaneous95']=lower;out.loc[mc,'p_upper_simultaneous95']=upper
    assert np.allclose(out.q_BH_joint,bh(out.p),atol=1e-14,rtol=0)
    out['q_BH_joint_lower']=bh(out.p_lower_simultaneous95)
    out['q_BH_joint_upper']=bh(out.p_upper_simultaneous95)
    assert out.q_BH_joint_lower.le(out.q_BH_joint_upper+1e-14).all()
    out['stable_joint_q05_reject']=out.q_BH_joint_upper.lt(.05)
    out['stable_joint_q05_nonreject']=out.q_BH_joint_lower.ge(.05)
    out['joint_q05_boundary_uncertain']=~(out.stable_joint_q05_reject|out.stable_joint_q05_nonreject)
    return out


def well(ctx):
    pairs,axes=load(ctx)
    study=Study(ctx,'GSE115943');axes=full_study_axes(study)
    p=pairs[pairs.study_id.eq('GSE115943') & pairs.contrast.eq('P_M')].copy()
    assert p.library_key_a.eq(p.library_key_b).all()
    assert p.replicate_status_a.eq('biological_independent').all()
    core=ctx.functions(BASE+'inference/analyze.py',['signs','sign_tests','nested_summary'],
        dict(SIGN_CACHE={},SEED=20260930,B=131071,t=t))
    multi=ctx.functions(BASE+'metric_inference/analyze.py',['statistics'],dict(CORE=core))
    rows=[];coverage=[];used=[];omnibus=[]
    for (sid,gene),axis in axes.items():
        if sid!='GSE115943':continue
        w=p[p.gene_id.eq(gene)].copy()
        w['date_replicates']=w.groupby('time_a').sample_unit_a.transform('nunique')
        q=w[w.date_replicates.ge(2)].sort_values(['time_a','library_key_a'])
        eligible=len(q)>=8 and q.time_a.nunique()>=4
        coverage.append(dict(gene_id=gene,qualified_wells=len(w),qualified_dates=w.time_a.nunique(),
            duplicate_date_wells=len(q),duplicate_dates=q.time_a.nunique(),tested=eligible,
            exclusion_reason='none' if eligible else 'fewer_than_4_duplicate_dates'))
        if not eligible:continue
        assert not q.sample_unit_a.duplicated().any()
        weights=1/q.time_a.nunique()/q.date_replicates.to_numpy()
        delta=vectors(q);linear=axis[['WPP','RWUI','PDUI','PPUI']].to_numpy()
        effect,pvalue,stat,metadata,mean=multi.statistics(delta,weights,linear)
        for j,metric in enumerate(METRICS):
            rec=dict(gene_id=gene,gene_symbol=q.gene_symbol.iloc[0],metric=metric,effect=effect[j],p=pvalue[j],
                statistic=stat[j],n_units=len(q),n_dates=q.time_a.nunique(),PAS_number=len(axis),
                weighted_delta_usage=json.dumps(mean.tolist()),**metadata)
            if j<4:
                summary=core.nested_summary(delta@linear[:,j],q.time_a)
                assert abs(summary.pop('estimate')-effect[j])<1e-12
                rec.update(summary)
            rows.append(rec)
        ca=np.stack(q.molecules_a.map(json.loads));cb=np.stack(q.molecules_b.map(json.loads))
        count_delta=ca/ca.sum(1)[:,None]-cb/cb.sum(1)[:,None]
        combined=np.column_stack([count_delta,count_delta@linear[:,0],delta@linear[:,0]])
        pv,meta=core.sign_tests(combined,weights)
        omnibus.append(dict(gene_id=gene,omnibus_p=pv[0],count_WPP_p=pv[1],equal_MC_WPP_p=pv[2],**meta))
        q['inference_weight']=weights;used.append(q)
    tests=finish_tests(ctx,rows)
    stability=mc_stability(tests);ctx.table(stability,'mc_stability.tsv')
    ctx.table(stability.groupby('metric').agg(tests=('gene_id','size'),
        MC_tests=('exact_sign_enumeration',lambda x:int((~x.astype(bool)).sum())),
        stable_reject=('stable_joint_q05_reject','sum'),
        stable_nonreject=('stable_joint_q05_nonreject','sum'),
        boundary_uncertain=('joint_q05_boundary_uncertain','sum')).reset_index(),'mc_stability_summary.tsv')
    ctx.table(coverage,'coverage_and_exclusions.tsv')
    ctx.table(pd.concat(used) if used else p.iloc[:0],'qualified_well_pairs.tsv.gz')
    om=pd.DataFrame(omnibus,columns=['gene_id','omnibus_p','count_WPP_p','equal_MC_WPP_p',
        'exact_sign_enumeration','sign_patterns','omnibus_statistic','count_WPP_statistic','equal_MC_WPP_statistic'])
    pc=['omnibus_p','count_WPP_p','equal_MC_WPP_p']
    combined=bh(om[pc].to_numpy().ravel()).reshape(len(om),3)
    for i,c in enumerate(pc):om[c+'_q_family']=bh(om[c]);om[c+'_q_combined']=combined[:,i]
    ctx.table(om,'gene_tests.tsv.gz')
    source_endpoint(ctx,axes,study)
    ctx.finish(dict(fresh_well_eligibility=True,shared_signs_on_full_vectors=True,
        fresh_five_metric_BH=True,source_variance_once=True,simultaneous_MC_bounds=True,
        primary_p_q_unchanged_by_MC_audit=True),tests=len(tests),joint_family_size=len(tests),
        MC_tests=int((~tests.exact_sign_enumeration.astype(bool)).sum()),MC_confidence=.95)


def source_endpoint(ctx,axes,study):
    states=source_states(study,axes)
    ctx.table(states,'all_qualified_source_target_states.tsv.gz')
    design=ctx.read(ctx.native/'studies/WOT_GSE115943/sampling_design.tsv')
    states=states.merge(design[['library_key','state','time','sample_unit','replicate_status']],
        on=['library_key','state'],validate='many_to_one')
    rows=[];coverage=[]
    by_gene=dict(tuple(states.groupby('gene_id')))
    for sid,gene in axes:
        if sid!='GSE115943':continue
        w=by_gene.get(gene,states.iloc[:0])
        source=w[w.state.eq('S') & w.time.eq('D0')].sort_values('library_key')
        dates=[]
        if len(source)==2:
            sc=np.stack(source.site_counts.map(json.loads));sm=np.stack(source.site_mc.map(json.loads))
            for date,part in w[w.state.eq('P')].groupby('time'):
                if len(part)!=2:continue
                cc=np.stack(part.site_counts.map(json.loads));cm=np.stack(part.site_mc.map(json.loads))
                support=((cc[:,None,:]+sc[None,:,:]>=10)&(cm[:,None,:]+sm[None,:,:]>=2)).sum(2)
                if (support>=2).all():dates.append(part)
        ok=len(dates)>=4
        coverage.append(dict(gene_id=gene,source_units=len(source),duplicate_target_dates=len(dates),tested=ok))
        if not ok:continue
        axis=axes['GSE115943',gene]
        for metric in ['equal_MC_WPP','count_WPP']:
            def values(frame):
                if metric=='equal_MC_WPP':return frame.APA.to_numpy()
                c=np.stack(frame.site_counts.map(json.loads));return c@axis.position/c.sum(1)
            s=values(source);ts=[values(d) for d in dates]
            parts=np.r_[s.var(ddof=1)/2,[v.var(ddof=1)/len(v)/len(dates)**2 for v in ts]]
            dfs=np.r_[1,[len(v)-1 for v in ts]];variance=parts.sum();denom=(parts**2/dfs).sum()
            df=variance**2/denom if denom else np.nan;se=np.sqrt(variance)
            effect=np.mean([v.mean() for v in ts])-s.mean()
            pvalue=2*t.sf(abs(effect/se),df) if se else 1.
            margin=t.ppf(.975,df)*se if se else np.nan
            rows.append(dict(gene_id=gene,metric=metric,effect=effect,p=pvalue,SE=se,df=df,
                CI95_low=effect-margin,CI95_high=effect+margin,source_variance=s.var(ddof=1),
                source_variance_fraction=parts[0]/variance if variance else np.nan,
                inference_status='estimable' if se else 'degenerate_observed_between_well_variance',
                source_units=2,target_units=sum(map(len,dates)),n_dates=len(dates)))
    f=pd.DataFrame(rows,columns=['gene_id','metric','effect','p','SE','df','CI95_low','CI95_high',
        'source_variance','source_variance_fraction','inference_status','source_units','target_units','n_dates'])
    f['q_BH_joint']=bh(f.p);f['q_BH_metric']=f.groupby('metric').p.transform(bh)
    f['interpretation']='D0_to_late_2i_endpoint_combined_state_date_medium_contrast'
    ctx.table(f,'source_endpoint_WPP_tests.tsv.gz')
    ctx.table(coverage,'source_endpoint_coverage_and_exclusions.tsv',columns=['gene_id','source_units','duplicate_target_dates','tested'])


def pool(ctx):
    pairs,axes=load(ctx)
    axes=full_study_axes(Study(ctx,'GSE221739'))
    pairs=pairs[pairs.study_id.eq('GSE221739') & pairs.contrast.eq('P_M')].copy()
    assert pairs.library_key_a.eq(pairs.library_key_b).all()
    assert pairs.sample_unit_a.eq(pairs.sample_unit_b).all()
    assert pairs.replicate_status_a.eq('pooled_chip_sample').all()
    math=ctx.functions(BASE+'pool_inference/analyze.py',['sign_statistics'])
    rows=[];coverage=[];unit_rows=[]
    for (sid,gene),axis in axes.items():
        if sid!='GSE221739':continue
        f=pairs[pairs.gene_id.eq(gene)]
        units=f[['time_a','sample_unit_a']].drop_duplicates()
        duplicate_dates=int(units.groupby('time_a').size().ge(2).sum())
        ok=len(units)>=4 and duplicate_dates>=2
        coverage.append(dict(gene_id=gene,captures=len(f),independent_pools=len(units),
            duplicate_pool_dates=duplicate_dates,tested=ok))
        if not ok:continue
        delta=[];meta=[]
        for (date,unit),part in f.groupby(['time_a','sample_unit_a'],sort=True):
            d=vectors(part).mean(0);delta.append(d)
            meta.append(dict(gene_id=gene,time=date,sample_unit=unit,n_captures=len(part),
                captures=';'.join(sorted(part.library_key_a)),delta_usage=json.dumps(d.tolist())))
        u=pd.DataFrame(meta);weights=1/u.time.nunique()/u.groupby('time').sample_unit.transform('size').to_numpy()
        assert not u.sample_unit.duplicated().any()
        effect,pvalue,mean,null,extreme=math.sign_statistics(np.stack(delta),weights,axis[['WPP','RWUI','PDUI','PPUI']].to_numpy())
        u['weight']=weights;unit_rows.extend(u.to_dict('records'))
        for j,metric in enumerate(METRICS):
            rows.append(dict(gene_id=gene,gene_symbol=f.gene_symbol.iloc[0],metric=metric,effect=effect[j],p=pvalue[j],
                n_units=len(u),n_dates=u.time.nunique(),PAS_number=len(axis),weighted_delta_usage=json.dumps(mean.tolist()),
                sign_patterns=2**len(u),exact_sign_enumeration=True,minimum_two_sided_p=2/2**len(u)))
    tests=finish_tests(ctx,rows)
    ctx.table(coverage,'coverage_and_exclusions.tsv')
    ctx.table(unit_rows,'independent_pool_vectors.tsv.gz',columns=['gene_id','time','sample_unit','n_captures','captures','delta_usage','weight'])
    ctx.finish(dict(capture_aliases_aggregated_first=True,fresh_pool_eligibility=True,
        all_sign_patterns=True,fresh_five_metric_BH=True),tests=len(tests))


def initialize_count_worker(config,out):
    global COUNT_FISHER
    from conditional_count_test import fisher
    COUNT_FISHER=fisher


def count_worker(job):
    return COUNT_FISHER(job)


def count_results(ctx,jobs,workers):
    """Isolate Python-heavy conditional tests while retaining per-record RNGs."""
    assert isinstance(workers,int) and workers>=1
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
                             initializer=initialize_count_worker,initargs=(str(ctx.config_path),str(ctx.out))) as executor:
        yield from executor.map(count_worker,jobs,chunksize=16)


def count(ctx):
    pairs,axes=load(ctx)
    math=ctx.functions(BASE+'pair_count_tests/analyze.py',['count_effects'],
        dict(fisher_exact=fisher_exact,random_table=random_table,B=131071))
    jobs=[];rows=[]
    for i,r in enumerate(pairs.to_dict('records')):
        a=np.array(json.loads(r['molecules_a']),dtype=np.int64);b=np.array(json.loads(r['molecules_b']),dtype=np.int64)
        assert ((a+b)>0).sum()>=2
        identity=[r[k] for k in FAMILY+['gene_id']]
        jobs.append((i,a,b,seed(*identity)))
        row={k:r[k] for k in FAMILY+['gene_id','gene_symbol']}
        row.update(math.count_effects(a,b,axes[r['study_id'],r['gene_id']]))
        row.update(equal_MC_DWUI=r['delta_DWUI'],PAS_number=len(a));rows.append(row)
    out=pd.DataFrame(rows)
    for i,pvalue,method,draws,extreme,random_seed,n_pas in count_results(ctx,jobs,int(ctx.config.get('count_workers',4))):
        for key,value in dict(p=pvalue,method=method,draws=draws,extreme=extreme,seed=random_seed,test_PAS=n_pas).items():out.loc[i,key]=value
        if (i+1)%100==0:print('conditional count tests',i+1,'/',len(out),flush=True)
    mc=out.method.eq('Fisher_MC_2xn');alpha=.05/int(mc.sum()) if mc.any() else .05
    out['p_MC_lower']=out.p;out['p_MC_upper']=out.p
    x=out.loc[mc,'extreme'].to_numpy();n=out.loc[mc,'draws'].to_numpy()
    out.loc[mc,'p_MC_lower']=np.minimum(out.loc[mc,'p'],np.where(x==0,0,beta.ppf(alpha/2,x,n-x+1)))
    out.loc[mc,'p_MC_upper']=np.maximum(out.loc[mc,'p'],np.where(x==n,1,beta.ppf(1-alpha/2,x+1,n-x)))
    for column,label in [('p','q_BH'),('p_MC_lower','q_BH_MC_lower'),('p_MC_upper','q_BH_MC_upper')]:
        out[label]=out.groupby(FAMILY,dropna=False)[column].transform(bh)
    out['MC_classification_uncertain']=out.q_BH_MC_lower.lt(.05)&out.q_BH_MC_upper.ge(.05)
    out['q_point_supported']=out.q_BH.lt(.05);out['q_MC_upper_supported']=out.q_BH_MC_upper.lt(.05)
    for metric in METRICS:
        for threshold in ([.1,.2] if metric=='MPRO' else [.02,.05]):
            tag=metric+'_'+str(threshold).replace('.','p')
            effect=out['count_'+metric].abs().ge(threshold)
            out[tag+'_point_call']=out.q_point_supported&effect
            out[tag+'_MC_stable_call']=out.q_MC_upper_supported&effect
    ctx.table(out,'gene_pair_Fisher_tests.tsv.gz')
    ctx.table(out.groupby(FAMILY,dropna=False).agg(genes=('gene_id','size'),
        q_point_supported=('q_point_supported','sum'),q_MC_stable_supported=('q_MC_upper_supported','sum'),
        MC_classification_uncertain=('MC_classification_uncertain','sum')).reset_index(),'FDR_families.tsv')
    ctx.finish(dict(integer_count_input=True,full_axis_effects=True,double_zero_removed_only_for_test=True,
        fresh_comparison_BH=True,simultaneous_MC_intervals=True,robust_probability_tie_comparison=True),
        tests=len(out),MC_draws=131071,probability_ordering='fixed-margin factorial products; exact integer comparison near floating roundoff')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--node',choices=['well_inference','pool_inference','count_inference'],required=True)
    a=p.parse_args();{'well_inference':well,'pool_inference':pool,'count_inference':count}[a.node](Context(a.config,a.node,a.out))
