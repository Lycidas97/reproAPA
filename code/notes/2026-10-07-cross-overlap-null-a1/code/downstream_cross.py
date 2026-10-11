"""Fresh observed-state cross-species cohorts and conditional overlap nulls."""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import beta,hypergeom,nchypergeom_fisher

from downstream_runtime import Context,Study,bh,seed

CODES={'W':('GSE115943','WOT'),'C':('GSE103221','CIP'),'H1':('GSE118258','GSE118258'),
       'H2':('GSE150311','GSE150311'),'H3':('GSE178325','178325_0618'),'H4':('GSE221739','GSE221739')}
MOUSE=['W','C'];HUMAN=['H1','H2','H3','H4']
MATCH='notes/2026-09-29-core-sensitivity/figure4_matched_overlap.py'
CONTROL='notes/2026-10-02-null-model-controls/'
MODELS={'L0_measurement':['effect_n'],'L1_capacity':['effect_n','prox_capacity'],
        'L2_baseline':['effect_n','prox_capacity','ref_q5'],'L2_PAS_first':['effect_n','prox_capacity','pas_q3'],
        'L3_full':['effect_n','prox_capacity','ref_q5','pas_q3']}


def qcodes(values,bins=5):
    out=pd.Series(-1,index=values.index,dtype=int);mask=values.notna()
    if mask.any():out.loc[mask]=pd.qcut(values[mask],bins,labels=False,duplicates='drop').fillna(0).astype(int)
    return out


def proximal(frame,codes,tau=.02):
    x=frame.reindex(columns=codes)
    if len(codes)==2:return x.le(-tau).all(axis=1)
    return x.le(-tau).sum(axis=1).ge(2)&x.ge(tau).sum(axis=1).eq(0)


def collect(ctx):
    orth=ctx.read(ctx.resource('ortholog'));assert not orth.human_gene_id.duplicated().any() and not orth.mouse_gene_id.duplicated().any()
    parts=[];libs=[];pm=[];sid_cache={}
    for code,(sid,route) in CODES.items():
        if sid not in sid_cache:sid_cache[sid]=Study(ctx,sid)
        study=sid_cache[sid]
        for setting,layer in [('supported50_p08','supported50'),('expanded_primary_p08','expanded_primary')]:
            nodes=study.core('node_positions',setting);nodes=nodes[nodes.route.eq(route)].copy()
            geometry={};support={}
            groups=[g.index.to_numpy() for _,g in study.meta[study.meta.route.eq(route)&study.meta.purity.ge(.8)].groupby(['library_key','node'])
                    if len(g)>=3 and g.n_cells.sum()>=30]
            for gene in nodes.gene_id.unique():
                axis=study.axis(gene,setting);supported=study.axis(gene,'supported50_p08')
                geometry[gene]=(len(supported),supported.cut_site.max()-supported.cut_site.min())
                data=study.gene_counts(axis);rows=np.concatenate(groups) if groups else np.array([],int)
                counts=data[rows];counts=counts[counts.sum(1)>=20]
                support[gene]=int(((counts.sum(0)>=10)&((counts>0).sum(0)>=2)).sum())
            nodes['supported_pas_in_route']=nodes.gene_id.map(support)
            nodes['supported_sites']=nodes.gene_id.map({g:v[0] for g,v in geometry.items()})
            nodes['supported_span']=nodes.gene_id.map({g:v[1] for g,v in geometry.items()})
            nodes['code']=code;nodes['layer']=layer
            nodes['measurement_pass']=nodes.position_index.notna()&nodes.qualified_metacells.ge(10)&nodes.supported_pas_in_route.ge(2)&nodes.supported_span.ge(200)
            parts.append(nodes)
            lib=study.core('library_states',setting);lib=lib[lib.route.eq(route)].copy();libs.append(lib.assign(code=code,layer=layer))
            design=study.design[study.design.route.eq(route)][['library_key','time','sample_unit']].drop_duplicates()
            assert not design.library_key.duplicated().any()
            for gid,q in lib.groupby('gene_id'):
                geo=geometry[gid]
                if geo[0]<2 or geo[1]<200 or support[gid]<2:continue
                wide=q.pivot(index='library_key',columns='node',values='position_index').reindex(columns=['P','M']).dropna()
                if not len(wide):continue
                wide['delta']=wide.P-wide.M
                matched=wide.reset_index().merge(design,on='library_key',validate='one_to_one')
                units=matched.groupby(['time','sample_unit']).delta.mean()
                if len(units)>=(2 if code in MOUSE else 1):
                    pm.append(dict(layer=layer,code=code,gene_id=gid,mean_delta=units.groupby(level='time').mean().mean(),
                        n_sample_units=len(units),n_libraries=len(wide),supported_span=geo[1],supported_sites=geo[0],quality_pass=True))
    nodes=pd.concat(parts,ignore_index=True);libraries=pd.concat(libs,ignore_index=True)
    mapper=orth.set_index('mouse_gene_id').human_gene_id
    nodes['human_gene_id']=np.where(nodes.species.eq('human'),nodes.gene_id,nodes.gene_id.map(mapper))
    libraries['human_gene_id']=np.where(libraries.species.eq('human'),libraries.gene_id,libraries.gene_id.map(mapper))
    paired=pd.DataFrame(pm,columns=['layer','code','gene_id','mean_delta','n_sample_units','n_libraries','supported_span','supported_sites','quality_pass'])
    paired['human_gene_id']=paired.gene_id.map(dict(zip(orth.mouse_gene_id,orth.human_gene_id))|dict(zip(orth.human_gene_id,orth.human_gene_id)))
    return orth,nodes,libraries,paired


def matrices(nodes,pm,layer):
    f=nodes[nodes.layer.eq(layer)&nodes.measurement_pass];wide=f.pivot(index='gene_id',columns=['code','node'],values='position_index')
    refs={};effects={}
    for contrast,state in [('PS','S'),('PM','M')]:
        refs[contrast]=pd.DataFrame({c:wide.get((c,state),pd.Series(np.nan,index=wide.index)) for c in CODES})
        effects[contrast]=pd.DataFrame({c:wide.get((c,'P'),pd.Series(np.nan,index=wide.index))-refs[contrast][c] for c in CODES})
    effects['PM']=pm[pm.layer.eq(layer)].pivot(index='gene_id',columns='code',values='mean_delta').reindex(columns=CODES)
    return effects,refs


def feature_frame(orth,nodes,pm,layer,tau=.02):
    f=orth.copy().reset_index(drop=True);effects,refs=matrices(nodes,pm,layer);n=nodes[nodes.layer.eq(layer)&nodes.measurement_pass]
    for prefix,codes,idcol in [('h',HUMAN,'human_gene_id'),('m',MOUSE,'mouse_gene_id')]:
        z=n[n.code.isin(codes)]
        for label,values in [('n',z.groupby('gene_id').code.nunique()),('pas',z.groupby('gene_id').supported_pas_in_route.median())]:
            f[prefix+'_'+label]=f[idcol].map(values).fillna(-1).round().astype(int)
        for contrast in ['PS','PM']:
            e=effects[contrast].reindex(index=f[idcol],columns=codes).reset_index(drop=True)
            r=refs[contrast].reindex(index=f[idcol],columns=codes).reset_index(drop=True).where(e.notna())
            stem=f'{prefix}_{contrast}'
            f[stem]=proximal(e,codes,tau).to_numpy();f[stem+'_effect_n']=r.notna().sum(1)
            # Paired effects can remain measurable when the node-reference proxy is missing.
            f[stem+'_actual_effect_n']=e.notna().sum(1)
            assert not (f[stem]&f[stem+'_actual_effect_n'].lt(2)).any()
            for op in ['median','min','max']:f[stem+'_ref_'+op]=getattr(r,op)(axis=1)
            f[stem+'_prox_capacity']=r.ge(tau).sum(1);f[stem+'_ref_q5']=qcodes(f[stem+'_ref_median'])
            f[stem+'_pas_q3']=qcodes(f[prefix+'_pas'].replace(-1,np.nan),3)
            for op in ['min','max']:
                value=f[stem+'_ref_'+op];f[stem+'_ref_'+op+'10']=np.floor(value*10).clip(upper=9).fillna(-1).astype(int)
            f[stem+'_mean_effect']=e.mean(1)
    return f


def overlap(ctx,cached=False):
    if cached:
        frozen=ctx.read(ctx.analysis/'cross_overlap/ortholog_matched_features.tsv.gz')
        orth=ctx.read(ctx.resource('ortholog'))
        nodes=ctx.read(ctx.analysis/'cross_overlap/all_mapped_nodes.tsv.gz')
        pm=ctx.read(ctx.analysis/'cross_overlap/paired_pm_study_effects.tsv.gz')
    else:
        orth,nodes,libs,pm=collect(ctx)
        ctx.table(nodes,'all_mapped_nodes.tsv.gz');ctx.table(libs,'all_mapped_library_nodes.tsv.gz');ctx.table(pm,'paired_pm_study_effects.tsv.gz')
    joint=ctx.functions(MATCH,['analytic_stratum','matched_null'],{'hypergeom':hypergeom,'SEED':seed('cross_joint')%(2**32),'B':5000})
    numerical=ctx.functions(MATCH,['analytic_stratum','matched_null','independent_null'],{'hypergeom':hypergeom,'SEED':seed('cross_overlap')%(2**32),'B':20000})
    primary=feature_frame(orth,nodes,pm,'supported50')
    if cached:pd.testing.assert_frame_equal(primary[frozen.columns],frozen,check_dtype=False,rtol=0,atol=1e-12)
    ctx.table(primary,'ortholog_matched_features.tsv.gz')
    archive={}
    rows=[];diagnostics=[];members=[];sens=[]
    for layer in ['supported50','expanded_primary']:
        for tau in [.01,.02,.03,.05]:
            f=primary if layer=='supported50' and tau==.02 else feature_frame(orth,nodes,pm,layer,tau)
            for co in ['PS','PM']:
                hit=f[f['h_'+co]&f['m_'+co]].copy();hit['contrast']=co;hit['layer']=layer;hit['threshold']=tau;members.append(hit)
                row,diagnostic=joint.matched_null(f,co,'joint_measurement_PAS',['h_n','m_n','h_pas','m_pas',f'h_{co}_actual_effect_n',f'm_{co}_actual_effect_n'],draws=archive if layer=='supported50' and tau==.02 else None)
                sens.append(row|dict(layer=layer,threshold=tau));diagnostics.append(diagnostic.assign(layer=layer,threshold=tau,side='joint'))
                if layer!='supported50' or tau!=.02:continue
                rows.append(row)
                for model,cols in [('independent_q5_reference_capacity',['effect_n','actual_effect_n','pas_q3','prox_capacity','ref_q5']),
                                   ('independent_range10_reference_capacity',['effect_n','actual_effect_n','pas_q3','prox_capacity','ref_min10','ref_max10'])]:
                    row,diagnostic=numerical.independent_null(f,co,model,[f'h_{co}_{x}' for x in cols],[f'm_{co}_{x}' for x in cols],draws=archive)
                    rows.append(row);diagnostics.append(diagnostic.assign(layer=layer,threshold=tau))
    results=pd.DataFrame(rows);results['q_within_model_two_contrasts']=results.groupby('model').empirical_p.transform(bh)
    sensitivity=pd.DataFrame(sens);sensitivity['q_16']=bh(sensitivity.empirical_p)
    ctx.table(results,'overlap_null_results.tsv');ctx.table(sensitivity,'overlap_threshold_layer_sensitivity.tsv')
    ctx.table(pd.concat(diagnostics,ignore_index=True),'stratum_diagnostics.tsv.gz')
    ctx.table(pd.concat(members,ignore_index=True),'new_shared_members.tsv.gz')
    np.savez_compressed(ctx.out/'overlap_null_draws.npz',**archive)
    ctx.finish(dict(complete_ortholog_reference_universe=True,new_proximal_sets=True,observed_S_and_paired_PM_separate=True,
        PM_node_reference_proxy_explicit=True,actual_effect_study_count_in_all_null_strata=True,node10MC_and_supported_span200=True,PM_group_qualified_route_support=True,new_covariates_and_strata=True),threshold=.02,null_draws=20000)


def controls(ctx):
    f=ctx.read(ctx.analysis/'cross_overlap/ortholog_matched_features.tsv.gz')
    common=ctx.functions(CONTROL+'common.py',['groups','binary_draw','mobility','probability','rate'],{'beta':beta})
    num=ctx.functions(CONTROL+'cross_species.py',['independent','injected'],dict(binary_draw=common.binary_draw,hypergeom=hypergeom,nchypergeom_fisher=nchypergeom_fisher))
    rows=[];strata=[];balances=[];powers=[];archives={}
    for co in ['PS','PM']:
        x=f['h_'+co].to_numpy(bool);y=f['m_'+co].to_numpy(bool);observed=int((x&y).sum());cells={};nulls={}
        for model,cols in MODELS.items():
            hg=common.groups(f,[f'h_{co}_{c}' for c in cols]);mg=common.groups(f,[f'm_{co}_{c}' for c in cols]);cells[model]=(hg,mg)
            hp,hd=common.mobility(x,hg);mp,md=common.mobility(y,mg)
            null=num.independent(x,y,hg,mg,20000,np.random.default_rng(seed('cross',co,model)));nulls[model]=null;archives[co+'|'+model]=null
            row=dict(contrast=co,model=model,n_orthologs=len(f),human_events=int(x.sum()),mouse_events=int(y.sum()),
                observed=observed,expected=float(hp@mp),mc_mean=null.mean(),null_sd=null.std(ddof=1),excess=observed-float(hp@mp),
                p=common.probability(null,observed),null_lo=np.quantile(null,.025),null_hi=np.quantile(null,.975),
                **{'human_'+k:v for k,v in hd.items()},**{'mouse_'+k:v for k,v in md.items()})
            for side,lab,prob,groups in [('h',x,hp,hg),('m',y,mp,mg)]:
                for j,ii in enumerate(groups):strata.append(dict(contrast=co,model=model,side=side,stratum=j,n=len(ii),events=int(lab[ii].sum())))
                for col in ['effect_n','prox_capacity','ref_median','pas_q3']:
                    value=f[f'{side}_{co}_{col}'].to_numpy(float);valid=np.isfinite(value);selected=lab&valid;weight=prob*valid
                    balances.append(dict(contrast=co,model=model,side=side,covariate=col,selected_mean=value[selected].mean() if selected.any() else np.nan,
                        null_selected_mean=np.nansum(prob*value)/weight.sum() if weight.sum() else np.nan,selected_missing=int((lab&~valid).sum()),expected_missing=float(prob[~valid].sum())))
                row[side+'_expected_events_below_capacity']=float(prob[f[f'{side}_{co}_prox_capacity'].lt(2)].sum())
            rows.append(row)
        for odds in [1,2,4,8]:
            sim=num.injected(x,y,*cells['L3_full'],3000,odds,np.random.default_rng(seed('cross_inject',co,odds)));archives[co+'|OR'+str(odds)]=sim
            for model,null in nulls.items():
                p=common.probability(null,sim)
                for alpha in [.05,.025]:powers.append(dict(contrast=co,generator='L3_full',tested_model=model,odds_ratio=odds,alpha=alpha,
                    experiment='conditional_injection',simulated_mean=sim.mean(),simulated_sd=sim.std(ddof=1),**common.rate(p<=alpha)))
        for model,(hg,mg) in cells.items():
            sim=num.independent(x,y,hg,mg,3000,np.random.default_rng(seed('cross_calibration',co,model)))
            powers.append(dict(contrast=co,generator=model,tested_model=model,odds_ratio=1,alpha=.05,experiment='self_calibration',
                simulated_mean=sim.mean(),simulated_sd=sim.std(ddof=1),**common.rate(common.probability(nulls[model],sim)<=.05)))
    out=pd.DataFrame(rows);out['q_two_contrasts']=out.groupby('model').p.transform(bh);out['q_all_10']=bh(out.p)
    ctx.table(out,'cross_species_ladder.tsv');ctx.table(strata,'cross_species_strata.tsv.gz');ctx.table(balances,'cross_species_balance.tsv');ctx.table(powers,'cross_species_power.tsv')
    np.savez_compressed(ctx.out/'cross_species_draws.npz',**archives)
    ctx.finish(dict(current_cohort_fixed_across_models=True,new_nulls=True,planned_ten_test_family=True,event_margins_preserved=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--node',choices=['cross_overlap','cross_controls'],required=True)
    p.add_argument('--cached-overlap',action='store_true',help='Rerun overlap nulls using frozen measured effects')
    a=p.parse_args();ctx=Context(a.config,a.node,a.out)
    if a.node=='cross_overlap':overlap(ctx,cached=a.cached_overlap)
    else:
        assert not a.cached_overlap
        controls(ctx)
