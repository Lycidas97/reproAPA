"""Numerical display summaries computed before rendering; no graphical imports."""
import argparse
import numpy as np
import pandas as pd
from downstream_runtime import Context


def distribution(frame,groups,columns):
    rows=[]
    for keys,q in frame.groupby(groups,dropna=False):
        if len(groups)==1:keys=(keys,) if not isinstance(keys,tuple) else keys
        for col in columns:
            x=q[col].dropna().to_numpy(float)
            rec=dict(zip(groups,keys));rec.update(measure=col,n=len(x),mean=float(x.mean()) if len(x) else np.nan)
            rec.update({f'q{int(p*100):02d}':float(np.quantile(x,p)) if len(x) else np.nan for p in [0,.1,.25,.5,.75,.9,1]})
            rows.append(rec)
    return pd.DataFrame(rows)


def run(ctx):
    read=lambda x:ctx.read(ctx.analysis/x)
    f=read('profiles/trajectory_geometry.tsv.gz')
    keys=[x for x in ['species','route','cluster'] if x in f]
    columns=[x for x in ['path_length_L','endpoint_distance_D','D_over_L','transient_fraction_T','endpoint_retention','L','D','T'] if x in f]
    # Source schemas are pinned; preserve every numerical geometry field when aliases differ.
    if not columns:columns=[c for c in f.select_dtypes('number') if c not in ['depth','k']]
    ctx.table(distribution(f,keys,columns),'profile-geometry-quantiles.tsv')
    a=read('profiles/architecture_gene_metrics.tsv.gz')
    columns=[c for c in ['contribution_SR','contribution_RE','contribution_EP','endpoint_retention','endpoint_retention_D_over_L','transient_fraction_T'] if c in a]
    ctx.table(distribution(a,['species','route','cluster'],columns),'stage-contribution-means-and-quantiles.tsv')
    members=read('profiles/gene_membership_k3.tsv.gz')
    counts=members.groupby(['species','cluster']).size().rename('genes').reset_index()
    counts['species_total']=counts.groupby('species').genes.transform('sum');counts['fraction']=counts.genes/counts.species_total
    ctx.table(counts,'profile-species-composition.tsv')
    robust=read('profile_controls/clustering_robustness.tsv')
    ctx.table(distribution(robust,[c for c in ['kind','method','k','reference_method','status'] if c in robust],
        [c for c in ['ARI','centroid_RMSE','assignment_agreement','mean_Jaccard','minimum_Jaccard'] if c in robust]),'clustering-robustness-quantiles.tsv')
    well=read('well_inference/gene_metric_tests.tsv.gz');r=read('recurrence/recurrence_members.tsv.gz')
    r=r[r.species.eq('mouse')&r.contrast.eq('P_M')&r.metric.eq('WPP')&r.threshold.eq(.02)&r.series_matching]
    rows=[]
    for metric,q in well.groupby('metric'):
        for bg,z in [('all_tested',q),('recurrent_P_M_intersection',q[q.gene_id.isin(r.gene_id)])]:
            rows.append(dict(metric=metric,background=bg,tested=len(z),joint_q05=int(z.q_BH_joint.lt(.05).sum()),correction='original all-gene all-metric joint family; not recomputed in subset'))
    ctx.table(rows,'well-denominators.tsv')
    qc=read('pas_qc/catalog_summary.tsv');qc=qc[qc.layer.eq('primary')].copy()
    qc['supported50_PAS_fraction']=qc.supported50_sites/qc.sites
    qc['supported50_UMI_fraction']=qc.supported50_UMI/qc.UMI
    ctx.table(qc,'primary-PAS-and-UMI-support.tsv')
    frames=[]
    for subset,path in [('All genes','global/global_direction_summary.tsv'),('Tandem-compatible','tandem/tandem_utr_summary.tsv')]:
        x=read(path);x=x[x.setting.eq('supported50_p08')&x.background.eq('available')&x.contrast.eq('P_S')&x.design.eq('true_source')].copy();x['subset']=subset;frames.append(x)
    ctx.table(pd.concat(frames),'tandem-comparison.tsv')
    interval=read('site_statistics/region_consistency_statistics.tsv')
    cols=[c for c in ['species','contrast','category','unit','n','genes','cobs','cobs_ci_low','cobs_ci_high','rho_mean','rho_mean_ci_low','rho_mean_ci_high'] if c in interval]
    ctx.table(interval[cols],'observed-interval-agreement.tsv')
    cv=read('cis/paired_cv_uncertainty.tsv');ctx.table(cv,'cis-model-increments.tsv')
    inventory=read('global/state_time_inventory.tsv');planned=read('global/design_qualified_planned_pairs.tsv') if (ctx.analysis/'global/design_qualified_planned_pairs.tsv').exists() else None
    ctx.table(inventory,'state-time-inventory.tsv')
    # Separate numerical arrays retain full test/background scope. Renderers may
    # filter labels and marks, never recalculate these tests or choose examples.
    from downstream_display import run as display
    display(ctx)
    ctx.finish(dict(display_statistics_outside_rendering=True,no_old_numbers=True,well_denominators_separate=True,
                    observed_agreement_only=True,quantiles_explicit=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'panel_statistics',a.out))
