"""Descriptive numerical products formerly calculated in publication renderers."""
import json
import numpy as np
import pandas as pd


def rank_bins(frame,bins=8):
    rows=[];members=[];summaries=[]
    f=frame[frame.design.eq('PMS')&frame.tolerance.eq(25)&frame.mass_min.ge(.8)]
    for pair in ['H1_C','H1_W']:
        for contrast in ['PM','PS','MS']:
            q=f[f.pair.eq(pair)].dropna(subset=[f'rank_{contrast}_h',f'rank_{contrast}_m']).copy()
            x=q[f'rank_{contrast}_h'];y=q[f'rank_{contrast}_m'];n=len(q)
            key=dict(pair=pair,contrast=contrast)
            rho=x.rank().corr(y.rank()) if n>1 else np.nan
            status='estimated' if n>=bins else 'fewer_genes_than_planned_bins'
            summaries.append(dict(key,n=n,rho=rho,bins=bins,status=status))
            # Preserve the source-table order to resolve equal human effects.
            labels=pd.qcut(x.rank(method='first'),bins,labels=False) if n>=bins else pd.Series(-1,index=q.index)
            for row,bin_id in zip(q.itertuples(),labels):
                members.append(dict(key,human_gene_id=row.human_gene_id,bin=int(bin_id)))
            for b in range(bins):
                xx=x[labels.eq(b)];yy=y[labels.eq(b)];nn=len(yy)
                rows.append(dict(key,bin=b,n=nn,x_mean=xx.mean(),y_mean=yy.mean(),
                    y_se=yy.std(ddof=1)/np.sqrt(nn) if nn>1 else np.nan,status=status,
                    uncertainty='within_bin_standard_error_not_gene_bootstrap_CI'))
    return pd.DataFrame(rows),pd.DataFrame(members,columns=['pair','contrast','human_gene_id','bin']),pd.DataFrame(summaries)


def distributions(frame,groups,measure,weights=None,edges=None):
    """Export raw observations, ECDF, quantiles and histogram; retain empty groups."""
    points=[];hist=[];summary=[]
    for key,g in frame.groupby(groups,dropna=False):
        key=key if isinstance(key,tuple) else (key,);base=dict(zip(groups,key))
        z=g[np.isfinite(g[measure])].sort_values(measure,kind='stable')
        x=z[measure].to_numpy(float);w=z[weights].to_numpy(float) if weights else np.ones(len(x))
        assert np.isfinite(w).all() and (w>=0).all()
        total=w.sum();valid=len(x)>0 and total>0
        cdf=np.cumsum(w)/total if valid else np.full(len(x),np.nan)
        quant=np.interp([.25,.5,.75],cdf,x) if weights and valid else np.quantile(x,[.25,.5,.75]) if valid else [np.nan]*3
        summary.append(dict(base,n=len(x),weight_sum=total,mean=np.average(x,weights=w) if valid else np.nan,
            q25=quant[0],q50=quant[1],q75=quant[2],minimum=x.min() if valid else np.nan,maximum=x.max() if valid else np.nan,
            status='estimated' if valid else 'no_finite_values',weighted=bool(weights)))
        points.extend(dict(base,order=i,value=v,weight=ww,cumulative=c) for i,(v,ww,c) in enumerate(zip(x,w,cdf)))
        if valid:
            narrow=x.max()-x.min()<64*np.finfo(float).eps*max(1,abs(x).max())
            bounds=edges if edges is not None else np.array([np.nextafter(x.min(),-np.inf),np.nextafter(x.max(),np.inf)]) if narrow else np.histogram_bin_edges(x,bins=min(32,max(8,int(np.sqrt(len(x))))))
            counts,_=np.histogram(x,bounds,weights=w)
            assert np.isclose(counts.sum(),total),'Histogram bounds exclude observations'
            hist.extend(dict(base,bin=i,left=a,right=b,weight=count,peak_fraction=count/counts.max())
                for i,(a,b,count) in enumerate(zip(bounds[:-1],bounds[1:],counts)))
    return (pd.DataFrame(points,columns=[*groups,'order','value','weight','cumulative']),
        pd.DataFrame(hist,columns=[*groups,'bin','left','right','weight','peak_fraction']),pd.DataFrame(summary,
        columns=[*groups,'n','weight_sum','mean','q25','q50','q75','minimum','maximum','status','weighted']))


def heatmap_long(f,members):
    """Centred route x node profiles of the selected k3 genes, with both heatmap orders.

    within_class_order: per route x class (mean centred position of that route).
    class_order: one gene order per species x class across every route x R/E/P/M cell, as drawn in the
    Figure 2A heatmap (mean of the available cells; ties keep the membership-table order)."""
    f=f[f.setting.eq('supported50_p08')&f.depth.eq(20)]
    f=f.merge(members[['species','gene_id','cluster']],on=['species','gene_id'],validate='many_to_one')
    w=f.pivot(index=['species','route','cluster','gene_id'],columns='node',values='position_index')
    w=w.reindex(columns=list('SREPM'));d=w.sub(w.S,axis=0)
    # Heatmap order is a statistic too: within each class, mean centred position.
    order=d.mean(axis=1).rename('ordering_mean').reset_index().sort_values(
        ['species','route','cluster','ordering_mean','gene_id'],kind='stable')
    order['within_class_order']=order.groupby(['species','route','cluster']).cumcount()
    long=d.reset_index().melt(id_vars=['species','route','cluster','gene_id'],var_name='node',value_name='delta_from_S')
    long=long.merge(order,on=['species','route','cluster','gene_id'],validate='many_to_one')
    cells=long[long.node.isin(list('REPM'))]
    cls=cells.groupby(['species','cluster','gene_id']).delta_from_S.mean().rename('class_ordering_mean').reset_index()
    rank=members[['species','gene_id']].reset_index(drop=True).reset_index().rename(columns={'index':'member_index'})
    cls=cls.merge(rank,on=['species','gene_id'],validate='one_to_one').sort_values(
        ['species','cluster','class_ordering_mean','member_index'],kind='stable',na_position='last')
    cls['class_order']=cls.groupby(['species','cluster']).cumcount()
    return long.merge(cls.drop(columns='member_index'),on=['species','cluster','gene_id'],validate='many_to_one')


def profiles(ctx,read):
    members=read('profiles/gene_membership_k3.tsv.gz')
    long=heatmap_long(read('profiles/all_route_node_positions.tsv.gz'),members)
    ctx.table(long,'profile-centred-heatmap.tsv.gz')
    ctx.table(long.groupby(['species','route','cluster','node'],as_index=False).agg(
        mean_delta=('delta_from_S','mean'),n=('delta_from_S','count')),'profile-centroid-means.tsv')
    a=read('profiles/architecture_gene_metrics.tsv.gz');rows=[]
    for key,g in a.groupby(['species','route','cluster']):
        base=dict(zip(['species','route','cluster'],key));x=np.sort(g.transient_fraction_T.dropna())
        for t in np.linspace(0,1,401):rows.append(dict(base,threshold=t,n=len(x),survival=(len(x)-np.searchsorted(x,t,side='left'))/len(x) if len(x) else np.nan))
    ctx.table(rows,'transient-survival.tsv',columns=['species','route','cluster','threshold','n','survival'])
    for label,left,right in [('SR','S','R'),('RE','R','E'),('EP','E','P')]:a['signed_'+label]=a['delta_'+right]-a['delta_'+left]
    ctx.table(a.groupby(['species','route','cluster'],as_index=False).agg(n=('gene_id','size'),
        signed_SR=('signed_SR','mean'),signed_RE=('signed_RE','mean'),signed_EP=('signed_EP','mean')),'signed-stage-means.tsv')
    from downstream_panel_statistics import distribution
    ctx.table(distribution(a,['cluster'],['endpoint_retention_D_over_L','transient_fraction_T','path_length_L']),
        'pooled-class-geometry.tsv')
    z=read('profiles/trajectory_amplitude_sensitivity.tsv')
    ctx.table(z.groupby(['cluster','measure','cutoff'],as_index=False).agg(median_of_route_medians=('median_T','median'),
        n_routes=('median_T','count')),'amplitude-sensitivity-display.tsv')
    rows=[]
    for cl,g in a[a.route.isin(['CIP','WOT'])].groupby('cluster'):
        w=g.pivot(index='gene_id',columns='route',values='dominant_stage').reindex(columns=['CIP','WOT']).dropna()
        left=w.CIP.value_counts(normalize=True);right=w.WOT.value_counts(normalize=True)
        rows.append(dict(cluster=cl,n_genes=len(w),observed=w.CIP.eq(w.WOT).mean(),
            marginal_expected=(left*right).fillna(0).sum() if len(w) else np.nan,scope='descriptive_independent_marginals_no_test'))
    ctx.table(rows,'dominant-stage-marginal-expectation.tsv',columns=['cluster','n_genes','observed','marginal_expected','scope'])


def cis(ctx,read):
    cv=read('cis/recomputed_cv.tsv');keys=['layer','contrast','discovery','tau'];group=keys+['model']
    ctx.table(cv.groupby(group,as_index=False).agg(auc=('auc','mean'),logloss=('logloss','mean'),
        n_repeats=('auc','count')),'cis-repeat-means.tsv')
    base=cv[cv.model.eq('M1')].drop(columns='model').rename(columns={'auc':'baseline_auc','logloss':'baseline_logloss','status':'baseline_status'})
    paired=cv.merge(base,on=keys+['repeat'],validate='many_to_one');paired['auc_gain']=paired.auc-paired.baseline_auc
    paired['logloss_reduction']=paired.baseline_logloss-paired.logloss
    ctx.table(paired,'cis-repeat-paired-gains.tsv')
    ctx.table(paired.groupby(group,as_index=False)[['auc_gain','logloss_reduction']].mean(),'cis-gain-means.tsv')
    f=read('cis/paired_cv_uncertainty.tsv').copy();reverse=f.model.eq('M0') ^ f.metric.eq('logloss')
    for col,source in [('delta','delta'),('lo','hi'),('hi','lo')]:
        f['display_'+col]=np.where(reverse,-f[source],f[col])
    f['display_comparison']=np.where(f.model.eq('M0'),'M1 minus M0',f.model+' minus M1')
    f['display_metric']=np.where(f.metric.eq('logloss'),'logloss_reduction','auc_gain')
    ctx.table(f,'cis-oriented-gains-and-intervals.tsv')
    x=read('cis/descriptive_matched_features.tsv.gz');x=x[x.layer.eq('supported50')&x.tau.eq(.02)]
    for key in ['A_g','phastcons_driver_mean']:
        for name,f in zip(['ecdf','histogram','quantiles'],distributions(x,['discovery','contrast','class'],key,'match_weight',np.linspace(0,1,13))):
            ctx.table(f,f'cis-{key}-{name}.tsv')


def regulatory(ctx,read):
    claims=read('union_pool/figure5_union_claims.tsv');models=read('intervals/model_effects.tsv');rows=[]
    stats=read('site_statistics/region_consistency_statistics.tsv');stats=stats[stats.category.eq('all')&stats.unit.eq('region_pair')]
    for sp in ['human','mouse']:
        for co in ['SP','MP']:
            genes=set(claims.loc[claims.species.eq(sp)&claims.contrast.eq(co),'gene_id'])
            q=models[models.species.eq(sp)&models.contrast.eq(co)&models.eligible&models.annotated_pas.ge(2)&models.n_ref50.ge(1)&models.n_target50.ge(1)]
            n=q.groupby('gene_id').route.nunique();valid=set(n[n.eq(2)].index)&genes
            q=stats[stats.species.eq(sp)&stats.contrast.eq(co)]
            assert len(q)==1
            rows.append(dict(species=sp,contrast=co,union_genes=len(genes),measurable_family_genes=len(valid),
                selected_interval_pairs=int(q.n.iloc[0]),selected_interval_genes=int(q.genes.iloc[0])))
    ctx.table(rows,'interval-admission-counts.tsv')
    region=read('site_statistics/all_changed_region_pairs.tsv.gz')
    ctx.table(region.groupby(['species','contrast','effect_mean_a','effect_mean_b']).size().rename('multiplicity').reset_index(),'interval-effect-ties.tsv')
    x=read('weighted_retention/gene_retention.tsv');x['cutoff_difference']=x.weighted_sites50-x.weighted_sites20
    for name,f in zip(['ecdf','histogram','quantiles'],distributions(x,['category'],'cutoff_difference')):ctx.table(f,'retention-cutoff-'+name+'.tsv')
    a=read('motif_background/motif_attrition.tsv')
    ctx.table(a.groupby(['species','contrast','status']).size().rename('windows').reset_index(),'motif-admission-counts.tsv')
    a=read('motif_controls/motif_power.tsv')
    ctx.table(a.groupby(['species','category','injected_sites_per_kb'],as_index=False).agg(
        minimum=('rejection_rate','min'),median=('rejection_rate','median'),maximum=('rejection_rate','max'),n_scenarios=('rejection_rate','size')),'motif-power-display.tsv')
    a=read('examples/inclusion_paired_libraries.tsv.gz')
    # Every exported current feature; no hard-coded 445-nt interval or privileged seed.
    keys=['feature_id','gene','route','cutoff']
    ctx.table(a.groupby(keys,as_index=False).agg(n=('library_key','size'),reference_mean=('reference_inclusion','mean'),
        target_mean=('target_inclusion','mean'),delta_mean=('delta','mean')),'inclusion-paired-means.tsv')
    a=read('examples/inclusion_library_states.tsv.gz')
    for name,f in zip(['ecdf','histogram','quantiles'],distributions(a,keys+['state'],'inclusion')):ctx.table(f,'inclusion-'+name+'.tsv')
    a=read('rbp_expression/cross_study_summary.tsv')
    for sp,g in a.groupby('species'):
        edges=np.histogram_bin_edges(g.effect_b.dropna(),bins=18) if g.effect_b.notna().any() else np.linspace(-1,1,19)
        for name,f in zip(['ecdf','histogram','quantiles'],distributions(g,['species','contrast'],'effect_b',edges=edges)):
            ctx.table(f,f'rbp-{sp}-marginal-{name}.tsv')


def cohort_summaries(ctx,read):
    x=read('recurrence/recurrence_members.tsv.gz');x=x[x.series_matching&x.variant.eq('primary')]
    ctx.table(x.groupby(['species','contrast','metric','threshold','direction']).size().rename('genes').reset_index(),'recurrence-direction-counts.tsv')
    rows=[]
    for (sp,co),g in x.groupby(['species','contrast']):
        a=set(g.loc[g.metric.eq('WPP')&g.threshold.eq(.02),'gene_id']);b=set(g.loc[g.metric.eq('MPRO')&g.threshold.eq(.1),'gene_id'])
        rows.append(dict(species=sp,contrast=co,DWUI_only=len(a-b),both=len(a&b),MPRO_only=len(b-a)))
    ctx.table(rows,'recurrence-metric-overlap.tsv',columns=['species','contrast','DWUI_only','both','MPRO_only'])
    x=read('functional/nested_universe_membership.tsv.gz')
    ctx.table(x.groupby(['species','contrast'],as_index=False)[['in_G0','in_G1_nested','in_G2']].sum(),'functional-universe-counts.tsv')
    x=pd.concat([read('functional/go_ladder.tsv.gz'),read('functional/kegg_ladder.tsv.gz')]);x['significant']=x.q_database_model.lt(.05)
    ctx.table(x.groupby(['model','comparison','database'],as_index=False).agg(tests=('significant','size'),q05=('significant','sum')),'functional-hit-counts.tsv')
    x=read('functional/functional_power.tsv.gz');x=x[x.model.eq('L3_full')]
    from downstream_panel_statistics import distribution
    ctx.table(distribution(x,['database','threshold','odds_ratio'],['rejection_probability']),'functional-power-display.tsv')
    x=read('profile_gex/rank_enrichment_all_tests.tsv.gz');x['significant']=x.q_both_databases.lt(.05)
    ctx.table(x.groupby(['species','route','layer','measure'],as_index=False).agg(tests=('significant','size'),q05=('significant','sum')),'profile-enrichment-counts.tsv')
    x=read('global/single_source_library_summary.tsv');x=x[x.contrast.eq('P_S')]
    ctx.table(x.groupby('route',as_index=False).agg(median=('mean','median'),minimum=('mean','min'),maximum=('mean','max'),n_settings=('mean','count')),'source-choice-ranges.tsv')
    x=read('apa_gex/four_class_nulls.tsv.gz');x=x[x.cell.eq('B')&x.baseline_bins.eq(5)&x.null_model.eq('direction_specific')]
    ctx.table(x.groupby(['species','threshold'],as_index=False).agg(comparable=('comparable_genes','max'),movable=('shuffleable_genes','min')),'apa-gex-mobility-counts.tsv')
    x=read('intermediate/process/frozen_expression_support_summary.tsv')
    ctx.table(x.groupby('node',as_index=False)[['support_q90','support_q95']].sum(),'intermediate-expression-support-counts.tsv')
    x=read('intermediate/process/same_series_process_classes.tsv.gz')
    x=x[x.scope.eq('requalified')&x.intermediate_threshold.eq(.02)&x.endpoint_retained&x.pair_measurable].copy()
    x['aligned']=x.pair_class.eq('endpoint_aligned')
    ctx.table(x.groupby(['species','state','variant','pair_class']).size().rename('genes').reset_index(),'intermediate-class-counts.tsv')
    ctx.table(x.groupby(['species','state','variant'],as_index=False).agg(n=('aligned','size'),aligned_fraction=('aligned','mean')),'intermediate-alignment-fractions.tsv')


def empirical_displays(ctx,read):
    specifications=[
        ('naive-same-library','naive_state/RSeT_same_library_triangle_gene_effects.tsv.gz',['version','setting','metric'],['hN_P','hN_S','P_S']),
        ('motif-density','motif_background/gene_density.tsv',['species','contrast','category'],['density_difference']),
        ('external-mirna-density','external_mirna/gene_family_densities.tsv',['species','contrast','group'],['difference']),
        ('apa-gex-profile','profile_gex/APA_GEX_gene_correlations.tsv.gz',['route','cluster'],['APA_GEX_r']),
        ('geometry','profiles/architecture_gene_metrics.tsv.gz',['species','route','cluster'],['endpoint_retention_D_over_L','transient_fraction_T']),
        ('profile-ora','profile_gex/enrichment_all_tests.tsv.gz',['species'],['q_value']),
        ('profile-rank','profile_gex/rank_enrichment_all_tests.tsv.gz',['layer','measure'],['q_both_databases']),
        ('functional-mobility','functional/functional_mobility.tsv',['model'],['expected_changed_fraction']),
        ('well-source-variance','well_inference/source_endpoint_WPP_tests.tsv.gz',['metric'],['source_variance_fraction']),
        ('same-state-recurrence','measurement/same-state/recurrence_draws.tsv.gz',['species','contrast','metric','threshold'],['pseudo_recurrent']),
    ]
    for prefix,path,groups,values in specifications:
        x=read(path)
        for value in values:
            for name,f in zip(['ecdf','histogram','quantiles'],distributions(x,groups,value)):
                ctx.table(f,f'{prefix}-{value}-{name}.tsv.gz')


def pas_coordinate_lookup(frame):
    keys=['study_id','primary_gene_id','atlas_pas_id']
    axis=frame[keys+['chrom','strand','cut_site']].drop_duplicates()
    assert not axis.duplicated(keys).any(),'Ambiguous PAS coordinate mapping'
    assert axis.strand.isin(['+','-']).all() and np.isfinite(axis.cut_site).all()
    assert np.equal(axis.cut_site,np.floor(axis.cut_site)).all()
    axis['oriented_cut']=axis.cut_site*np.where(axis.strand.eq('+'),1,-1)
    return axis.set_index(keys).oriented_cut


def coordinates_for_row(lookup,row):
    ids=json.loads(row.PAS_ids)
    keys=pd.MultiIndex.from_tuples([(row.study_id,row.gene_id,p) for p in ids])
    values=lookup.reindex(keys)
    assert values.notna().all(),'PAS absent from the frozen study/gene coordinate axis'
    return values.to_numpy(dtype=np.int64)


def pas_curves(ctx,read):
    selected=read('examples/example_eligibility.tsv');genes=set(selected.gene_id.dropna())
    curves=[];lookup=pas_coordinate_lookup(read('state_vectors/PAS_axes.primary.tsv.gz'))
    for source,path in [('source_endpoint','state_vectors/route_effects.primary.tsv.gz'),('same_contributors_PMS','matched_pms/matched_state_vectors.tsv.gz')]:
        x=read(path);x=x[x.gene_id.isin(genes)&x.setting.eq('supported50_p08')].copy()
        if source=='source_endpoint':x=x[x.contrast.eq('P_S')&x.design.eq('true_source')]
        axes={i:coordinates_for_row(lookup,row) for i,row in x.iterrows()}
        for gid,g in x.groupby('gene_id'):
            origin=min(axes[i].min() for i in g.index)
            for i,row in g.iterrows():
                coords=axes[i];order=np.argsort(coords);ids=json.loads(row.PAS_ids)
                states=[('S','reference_usage'),('P','target_usage')] if source=='source_endpoint' else [(s,s+'_usage') for s in 'SPM']
                for state,column in states:
                    usage=np.array(json.loads(row[column]));assert len(usage)==len(coords) and np.isclose(usage.sum(),1)
                    cumulative=np.cumsum(usage[order])
                    curves.extend(dict(source=source,species=row.species,gene_id=gid,gene_symbol=row.gene_symbol,route=row.route,
                        state=state,PAS_id=ids[j],oriented_position_nt=int(coords[j]-origin),usage=usage[j],cumulative=cumulative[k]) for k,j in enumerate(order))
    ctx.table(curves,'example-cumulative-PAS-curves.tsv.gz',columns=['source','species','gene_id','gene_symbol','route','state','PAS_id','oriented_position_nt','usage','cumulative'])


def run(ctx):
    read=lambda p:ctx.read(ctx.analysis/p)
    profiles(ctx,read);cis(ctx,read);regulatory(ctx,read);cohort_summaries(ctx,read);empirical_displays(ctx,read);pas_curves(ctx,read)
    # Default parser reproduces the original figure's rank tie handling exactly.
    f=pd.read_csv(ctx.register(ctx.analysis/'same_pas/v2_primary_same_axis_gene_effects.tsv.gz'),sep='\t')
    for name,frame in zip(['same-pas-rank-bins.tsv','same-pas-rank-bin-members.tsv','same-pas-rank-bin-summary.tsv'],rank_bins(f)):ctx.table(frame,name)
