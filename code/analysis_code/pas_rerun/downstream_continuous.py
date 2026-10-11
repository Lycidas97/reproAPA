"""Rebuild continuous ortholog comparisons and reciprocal discovery validation."""
import argparse
import itertools

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from downstream_runtime import Context,bh,seed
from downstream_cross import MOUSE,HUMAN,qcodes

METRICS=['spearman','pearson','direction_concordance','direction_score_all']
PAIRS={s:list(itertools.product(MOUSE,h)) for s,h in
       [('R',['H1','H3']),('E',['H1','H2']),('P',HUMAN),('M',HUMAN)]}
SOURCE='github/reproapa-analysis/modules/09_cross_species/code/figure4/'


def average(x,axis=None):
    a=np.asarray(x,float);n=np.isfinite(a).sum(axis=axis)
    return np.divide(np.nansum(a,axis=axis),n,out=np.full(np.shape(n),np.nan),where=n>0)


def standardized(x):
    x=np.asarray(x,float)
    return (x-x.mean())/x.std() if len(x)>1 and x.std()>0 else np.full(len(x),np.nan)


def prepare(x,y):
    return [np.asarray(x),np.asarray(y),standardized(rankdata(x)),standardized(rankdata(y)),standardized(x),standardized(y)]


def metrics(a,ix=None):
    x,y,rx,ry,zx,zy=a
    if ix is not None:y=y[ix];ry=ry[ix];zy=zy[ix]
    strong=(abs(x)>=.02)&(abs(y)>=.02);same=strong&(np.sign(x)==np.sign(y));opposite=strong&~same
    return np.array([average(rx*ry),average(zx*zy),same.sum()/strong.sum() if strong.any() else np.nan,
                     (same.sum()-opposite.sum())/len(x) if len(x) else np.nan])


def effects_frame(nodes,layer,mouse,human,states):
    frames=[]
    for code,prefix in [(mouse,'m'),(human,'h')]:
        n=nodes[nodes.layer.eq(layer)&nodes.code.eq(code)&nodes.measurement_pass&nodes.human_gene_id.notna()]
        w=n.pivot(index='human_gene_id',columns='node',values='position_index').reindex(columns=['S',*states]).dropna()
        frames.append(pd.DataFrame({prefix+'_'+s:w[s]-w.S for s in states}))
    return frames[0].join(frames[1],how='inner')


def null_record(value,null):
    valid=np.isfinite(null)
    estimable=bool(np.isfinite(value) and valid.all() and len(null)>0)
    return dict(observed=value,p_greater=(1+int((null>=value).sum()))/(len(null)+1) if estimable else 1.,
        null_mean=average(null),null_lo=np.nanquantile(null,.025) if valid.any() else np.nan,
        null_hi=np.nanquantile(null,.975) if valid.any() else np.nan,draws=len(null),
        status='estimated' if estimable else 'not_estimable',finite_null_draws=int(valid.sum()))


def state_comparisons(ctx,nodes):
    rows=[];pair_rows=[];obs=[];boot_rows=[];boot_pairs=[]
    for layer in ['supported50','expanded_primary']:
        for state,pairs in PAIRS.items():
            for matched in [False,True] if state!='P' else [False]:
                states=[state,'P'] if matched else [state];label='P_vs_'+state if matched else state
                frames=[effects_frame(nodes,layer,m,h,states) for m,h in pairs]
                ready=[[prepare(q['m_'+s],q['h_'+s]) for s in states] for q in frames]
                observed=average([[metrics(a) for a in block] for block in ready],axis=0)
                value=observed[1]-observed[0] if matched else observed[0]
                draws=10000 if state=='M' else 1000;rng=np.random.default_rng(seed('continuous',layer,label))
                null=np.empty((draws,4))
                for b in range(draws):
                    estimates=[]
                    for block in ready:
                        ix=rng.permutation(len(block[0][0]));estimates.append([metrics(a,ix) for a in block])
                    z=average(estimates,axis=0);null[b]=z[1]-z[0] if matched else z[0]
                for i,metric in enumerate(METRICS):
                    rows.append(dict(layer=layer,comparison=label,matched=matched,metric=metric,
                        **null_record(value[i],null[:,i]),state_value=observed[0,i],P_value=observed[-1,i] if matched else np.nan))
                for (mouse,human),q,block in zip(pairs,frames,ready):
                    for s,a in zip(states,block):
                        pair_rows.append(dict(layer=layer,comparison=label,state=s,mouse=mouse,human=human,n=len(q),
                            **dict(zip(METRICS,metrics(a)))))
                        if matched:
                            for gene,r in q.iterrows():
                                x,y=r['m_'+s],r['h_'+s]
                                obs.append(dict(layer=layer,comparison=label,state=s,pair=mouse+'-'+human,mouse_code=mouse,human_code=human,
                                    human_gene_id=gene,mouse_effect=x,human_effect=y,relation='weak_or_mixed' if min(abs(x),abs(y))<.02 else
                                    'concordant' if np.sign(x)==np.sign(y) else 'discordant'))
                if matched:continue
                genes=sorted(set().union(*(set(q.index) for q in frames)));arrays=[q.reindex(genes)[['m_'+state,'h_'+state]].to_numpy() for q in frames]
                boot=np.full((1000,len(frames)),np.nan)
                if genes:
                    for b in range(len(boot)):
                        ix=rng.integers(len(genes),size=len(genes))
                        for i,a in enumerate(arrays):
                            v=a[ix];v=v[np.isfinite(v).all(1)]
                            boot[b,i]=metrics(prepare(v[:,0],v[:,1]))[0]
                def bounds(v):
                    return dict(lo=np.nanquantile(v,.025),hi=np.nanquantile(v,.975)) if np.isfinite(v).any() else dict(lo=np.nan,hi=np.nan)
                boot_rows.append(dict(layer=layer,state=state,union_genes=len(genes),pairs=len(pairs),rho=value[0],**bounds(average(boot,axis=1))))
                for i,((m,h),q) in enumerate(zip(pairs,frames)):
                    boot_pairs.append(dict(layer=layer,state=state,mouse=m,human=h,n=len(q),rho=metrics(ready[i][0])[0],**bounds(boot[:,i])))
    result=pd.DataFrame(rows);result['q_family']=np.nan
    # The historical R/E/P family and later M family retain distinct identities.
    result['family']=np.where(result.matched,'matched_P_vs_R_E_M_12',np.where(result.comparison.eq('M'),'M_only_4','R_E_P_12'))
    result['q_within_layer']=result.groupby(['layer','family']).p_greater.transform(bh)
    result['q_both_layers']=result.groupby('family').p_greater.transform(bh)
    assert result[result.matched].groupby('layer').size().eq(12).all()
    ctx.table(result,'state_and_matched_permutation.tsv');ctx.table(pair_rows,'state_pair_metrics.tsv')
    ctx.table(boot_rows,'continuous_gene_block_summary.tsv');ctx.table(boot_pairs,'continuous_pair_summary.tsv')
    observations=pd.DataFrame(obs,columns=['layer','comparison','state','pair','mouse_code','human_code','human_gene_id','mouse_effect','human_effect','relation'])
    return observations


def technical(ctx,nodes,libs,obs):
    valid=nodes[nodes.measurement_pass&nodes.human_gene_id.notna()]
    columns=['n_libraries','qualified_metacells','supported_pas_in_route','supported_span']
    n=valid.set_index(['layer','human_gene_id','code','node'])[columns]
    sd=libs[libs.human_gene_id.notna()].groupby(['layer','human_gene_id','code','node']).position_index.std()
    rows=[]
    for r in obs.itertuples(index=False):
        keys=[(r.layer,r.human_gene_id,c,s) for c in [r.mouse_code,r.human_code] for s in ['S',r.state]]
        a=n.loc[keys];disp=sd.reindex(keys)
        rows.append(dict(r._asdict(),min_libraries=a.n_libraries.min(),min_metacells=a.qualified_metacells.min(),
            min_pas=a.supported_pas_in_route.min(),min_span=a.supported_span.min(),max_library_sd=disp.max(),all_nodes_multi_library=a.n_libraries.ge(2).all()))
    tech=pd.DataFrame(rows,columns=[*obs.columns,'min_libraries','min_metacells','min_pas','min_span','max_library_sd','all_nodes_multi_library'])
    keys=['comparison','state','pair','human_gene_id']
    rel=tech.pivot(index=keys,columns='layer',values='relation').reindex(columns=['supported50','expanded_primary'])
    rel['relation_stable']=rel.supported50.eq(rel.expanded_primary).where(rel.notna().all(1))
    tech=tech.merge(rel[['relation_stable']].reset_index(),on=keys,how='left',validate='many_to_one')
    ctx.table(tech,'technical_observation_audit.tsv.gz');tests=[]
    for cov in ['min_libraries','min_metacells','min_pas','min_span','max_library_sd','all_nodes_multi_library','relation_stable']:
        binary=cov in ['all_nodes_multi_library','relation_stable']
        q=tech[tech.layer.eq('supported50')&tech.relation.isin(['concordant','discordant'])].dropna(subset=[cov]).reset_index(drop=True)
        labels=q.relation.eq('discordant').to_numpy();values=q[cov].to_numpy(float)
        if not binary:values=np.log1p(values)
        def difference(l):return average(values[l])-average(values[~l])
        value=difference(labels);groups=list(q.groupby(['comparison','state','pair']).indices.values())
        rng=np.random.default_rng(seed('technical',cov));null=np.empty(1000)
        for b in range(len(null)):
            perm=labels.copy()
            for ix in groups:perm[ix]=labels[rng.permutation(ix)]
            null[b]=difference(perm)
        rec=null_record(abs(value),abs(null))
        tests.append(dict(covariate=cov,family='binary_2' if binary else 'continuous_5',n=len(q),difference=value,**rec))
    tests=pd.DataFrame(tests);tests['q_family']=tests.groupby('family').p_greater.transform(bh);ctx.table(tests,'technical_covariate_permutation.tsv')


def add_bins(x):
    x=x.copy();x['effect_bin']=qcodes(x.disc_meta.abs(),min(10,max(1,x.disc_meta.abs().nunique())))
    x['stratum']=x[['disc_direction','disc_eligible','val_eligible','effect_bin']].astype(str).agg('|'.join,axis=1) if len(x) else pd.Series(dtype=str)
    return x


def reciprocal(ctx,nodes,pm):
    ns=dict(nodes=nodes,pm=pm,codes_mouse=MOUSE,codes_human=HUMAN)
    parts=[]
    for layer in ['supported50','expanded_primary']:
        q=nodes[nodes.layer.eq(layer)&nodes.measurement_pass&nodes.code.isin(MOUSE+HUMAN)].copy()
        w=q.pivot(index=['human_gene_id','code'],columns='node',values='position_index').reindex(columns=['P','S'])
        ps=(w.P-w.S).rename('effect').dropna().reset_index()
        cov=q[q.node.isin(['P','S'])].groupby(['human_gene_id','code']).agg(pas=('supported_pas_in_route','min'),
            support=('qualified_metacells','min'),libraries=('n_libraries','min'),span=('supported_span','min')).reset_index()
        ps=ps.merge(cov,on=['human_gene_id','code'],validate='one_to_one').assign(contrast='P-S',layer=layer)
        z=pm[pm.layer.eq(layer)&pm.quality_pass&pm.code.isin(MOUSE+HUMAN)].copy()
        z=z[(z.code.isin(MOUSE)&z.n_sample_units.ge(2))|(z.code.isin(HUMAN)&z.n_sample_units.ge(1))]
        z=z[['human_gene_id','code','mean_delta','supported_sites','n_sample_units','n_libraries','supported_span']].rename(
            columns=dict(mean_delta='effect',supported_sites='pas',n_sample_units='support',n_libraries='libraries',supported_span='span'))
        parts.extend([ps,z.assign(contrast='P-M',layer=layer)])
    effects=pd.concat(parts,ignore_index=True)
    funcs=ctx.functions(SOURCE+'recurrence_predicts_cross_species_v26.py',['species_summary','validation_table','match_sets'],dict(ns,effects=effects,add_bins=add_bins))
    audit=ctx.functions(SOURCE+'audit_figure4c_reused_controls_v53.py',['odds','analysis'])
    ctx.table(effects,'study_effects_and_covariates.tsv.gz');universes=[];matches=[];tests=[];balance=[];doses=[];trends=[]
    for layer,contrast,discovery,tau in itertools.product(['supported50','expanded_primary'],['P-S','P-M'],['human','mouse'],[.01,.02,.03,.05]):
        key=dict(layer=layer,contrast=contrast,discovery=discovery,tau=tau);primary=layer=='supported50' and tau==.02
        x=funcs.validation_table(layer,contrast,discovery,tau);universes.append(x)
        m=funcs.match_sets(x) if len(x) else pd.DataFrame();n=m.match_id.nunique() if len(m) else 0
        record=dict(key,primary=primary,n_sets=n,status='not_estimable_fewer_than_5_matched_cases',
            exact_permutation_p_two_sided=1.,effect_bin_permutation_p_two_sided=1.)
        if n:
            m['analysis']='|'.join(map(str,key.values()));matches.append(m)
            for cov in ['disc_meta','disc_pas','disc_support','disc_libraries','disc_span']:
                z=m.copy();z['value']=z[cov].abs() if cov=='disc_meta' else z[cov]
                case=z[z.role.eq('recurrent')].set_index('match_id').value
                z=z[z.role.eq('control')];ctrl=(z.value*z.match_weight).groupby(z.match_id).sum().reindex(case.index)
                pooled=np.sqrt((case.var()+ctrl.var())/2)
                balance.append(dict(key,covariate=cov,case_mean=case.mean(),control_mean=ctrl.mean(),standardized_difference=(case.mean()-ctrl.mean())/pooled if pooled>0 else np.nan))
        if n>=5:record.update(audit.analysis(m,seed('reciprocal',*key.values()),20000 if primary else 5000),status='estimated')
        tests.append(record)
        for label,mask in [('low',x.disc_score.lt(.49)),('medium',x.disc_score.ge(.49)&x.disc_score.lt(.99)),('high',x.disc_score.ge(.99))]:
            q=x[mask];doses.append(dict(key,score_group=label,n=len(q),validation_rate=q.validation_same_direction.mean(),mean_signed_validation=q.signed_validation.mean()))
        if primary:
            z=add_bins(x).reset_index(drop=True);groups=list(z.groupby('stratum').indices.values())
            rng=np.random.default_rng(seed('score_trend',contrast,discovery))
            for outcome in ['signed_validation','validation_same_direction']:
                a=prepare(z.disc_score.to_numpy(),z[outcome].to_numpy(float));null=np.empty(500)
                for b in range(len(null)):
                    ix=np.arange(len(z))
                    for g in groups:ix[g]=rng.permutation(g)
                    null[b]=metrics(a,ix)[0]
                trends.append(dict(key,outcome=outcome,n=len(z),**null_record(metrics(a)[0],null)))
    tests=pd.DataFrame(tests)
    for mode in ['exact','effect_bin']:
        p=mode+'_permutation_p_two_sided';tests[mode+'_q_planned_32']=bh(tests[p]);tests[mode+'_q_primary_four']=np.nan
        tests.loc[tests.primary,mode+'_q_primary_four']=bh(tests.loc[tests.primary,p])
    ctx.table(pd.concat(universes,ignore_index=True),'reciprocal_gene_universe.tsv.gz')
    ctx.table(pd.concat(matches,ignore_index=True) if matches else [],'matched_comparators.tsv.gz',columns=['analysis','human_gene_id','match_id','role','match_weight'])
    ctx.table(tests,'reciprocal_unique_gene_tests.tsv');ctx.table(balance,'matching_balance.tsv',columns=['layer','contrast','discovery','tau','covariate'])
    ctx.table(doses,'recurrence_score_dose_response.tsv');trends=pd.DataFrame(trends);trends['q_eight']=bh(trends.p_greater);ctx.table(trends,'recurrence_score_trend_tests.tsv')
    leave=[]
    for contrast in ['P-S','P-M']:
        f=effects[effects.layer.eq('supported50')&effects.contrast.eq(contrast)]
        for discovery in ['human','mouse']:
            dc,vc=(HUMAN,MOUSE) if discovery=='human' else (MOUSE,HUMAN)
            for side,codes in [('discovery',dc),('validation',vc)]:
                for omitted in codes:
                    dd=[c for c in dc if not(side=='discovery' and c==omitted)];vv=[c for c in vc if not(side=='validation' and c==omitted)]
                    ds=funcs.species_summary(f,dd,'disc',.02);vs=funcs.species_summary(f,vv,'val',.02)
                    z=ds.join(vs,how='inner');z=z[z.disc_recurrent&z.val_eligible.ge(min(2,len(vv)))&z.disc_direction.notna()]
                    signed=z.val_meta*z.disc_direction
                    leave.append(dict(contrast=contrast,discovery=discovery,side=side,omitted=omitted,n=len(z),
                        validation_rate=signed.ge(.02).mean(),mean_signed_validation=signed.mean(),
                        status='not_estimable_single_discovery_study' if len(dd)<2 else 'descriptive'))
    ctx.table(leave,'leave_one_study_sensitivity.tsv')


def run(ctx):
    nodes=ctx.read(ctx.analysis/'cross_overlap/all_mapped_nodes.tsv.gz');libs=ctx.read(ctx.analysis/'cross_overlap/all_mapped_library_nodes.tsv.gz')
    pm=ctx.read(ctx.analysis/'cross_overlap/paired_pm_study_effects.tsv.gz')
    # Unmapped genes must never form a spurious shared NaN ortholog.
    nodes=nodes[nodes.human_gene_id.notna()];libs=libs[libs.human_gene_id.notna()];pm=pm[pm.human_gene_id.notna()]
    observations=state_comparisons(ctx,nodes);technical(ctx,nodes,libs,observations);reciprocal(ctx,nodes,pm)
    ctx.finish(dict(new_state_qualified_orthologs=True,matched_state_joint_permutations=True,matched_twelve_test_families=True,
        shared_ortholog_block_bootstrap=True,new_reciprocal_matches=True,reused_controls_resampled_once=True,
        all_32_planned_reciprocal_positions=True,unmeasurable_positions_explicit=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'cross_continuous',a.out))
