"""Effect, precision, conditioning and design-resolution diagnostics."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'notes/2026-10-02-reviewer-followup'))
from apa_followup import read,save,rng_for,nullmod,INPUTS
PROJECTS=Path('/Users/lycidas/work/科研/projects')
M=PROJECTS/'2026-09-02_mouse-time-aware-integration-pilot'
C=PROJECTS/'2026-09-22_reproapa-cis-supplement'
F=PROJECTS/'2026-09-13_reproapa-figure4-real'


def sharing():
    base=ROOT/'notes/2026-09-30-multimetric-statistics/nulls'
    d=read(base/'sharing_summary.tsv')
    d['excess']=d.observed-d.null_mean
    d['standardized_null_excess']=d.excess/d.null_sd
    d['raw_supported']=d.empirical_p<.05
    d['primary_joint_supported']=d.q_joint_primary_24<.05
    d['excess_per_1000_orthologs']=1000*d.excess/d.common_orthologs
    save(d,HERE/'sharing_diagnostics.tsv')
    mobility=read(base/'null_mobility.tsv')
    save(mobility[mobility.analysis.eq('sharing')],HERE/'sharing_mobility.tsv')
    fixed=read(ROOT/'notes/2026-09-29-core-sensitivity/figure4-matched-overlap/overlap_null_results.tsv')
    save(fixed,HERE/'fixed_proximal_overlap.tsv')
    # The newer published input fixes responder status and Guan series patterns.
    # The older strengthening output has different human nulls and is superseded.
    p=ROOT/'table/Supplementary-Source-Data-20260928/figure-data/Figure3/panel_B'
    d=read(p/'D_four_class_nulls.tsv.gz')
    d['excess']=d.observed-d.null_mean
    save(d[d.cell.eq('B')],HERE/'apa_gex_opposite_nulls.tsv')
    save(read(p/'D_prespecified_eight_test_family.tsv'),HERE/'apa_gex_primary_and_q3.tsv')
    mobility=d[['species','threshold','null_model','baseline_bins','shuffled_study','comparable_genes','shuffleable_genes']].drop_duplicates()
    mobility['shuffleable_fraction']=mobility.shuffleable_genes/mobility.comparable_genes
    save(mobility,HERE/'apa_gex_mobility.tsv')


def weighted_auc(y,p,weights):
    order=np.argsort(p,kind='stable'); ps=p[order];yy=y[order];w=weights[:,order]
    starts=np.r_[0,np.flatnonzero(np.diff(ps))+1]
    positives=np.add.reduceat(w*yy,starts,axis=1)
    negatives=np.add.reduceat(w*(1-yy),starts,axis=1)
    wins=(positives*(np.cumsum(negatives,axis=1)-.5*negatives)).sum(1)
    return wins/(positives.sum(1)*negatives.sum(1))


def cis():
    d=read(C/'data/primary_oof_predictions.tsv.gz')
    features=read(C/'data/primary_observed_features.tsv.gz')
    summaries=[]; missing=[]
    for keys,part in d.groupby(['contrast','discovery']):
        con,disc=keys
        genes=sorted(part.human_gene_id.unique())
        rng=rng_for('cis_gene_bootstrap',*keys)
        y=part.drop_duplicates('human_gene_id').set_index('human_gene_id').loc[genes,'Y'].to_numpy()
        weights=np.vstack([np.ones(len(genes)),rng.multinomial(len(genes),np.ones(len(genes))/len(genes),size=2000)])
        estimates={}
        for model,block in part.groupby('model'):
            pred=block.pivot(index='human_gene_id',columns='repeat',values='oof_probability').loc[genes].to_numpy()
            assert pred.shape==(len(genes),10)
            clipped=np.clip(pred,1e-15,1-1e-15)
            loss=-(y[:,None]*np.log(clipped)+(1-y[:,None])*np.log1p(-clipped))
            auc=np.mean([weighted_auc(y,pred[:,r],weights) for r in range(10)],axis=0)
            estimates[model]=(auc,weights@loss.mean(1)/len(genes),weights@((pred-y[:,None])**2).mean(1)/len(genes))
        for model,baseline in [('M1','M0'),('PAS','M1'),('Architecture','M1'),('Sequence','M1'),('All cis','M1')]:
            auc,ll,brier=estimates[model];ba,bl,bb=estimates[baseline]
            row=dict(contrast=con,discovery=disc,model=model,baseline=baseline,n=len(genes),positives=int(y.sum()),
                positive_fraction=y.mean(),baseline_auc=ba[0],model_auc=auc[0],baseline_logloss=bl[0],model_logloss=ll[0],
                constant_prevalence_logloss=-(y.mean()*np.log(y.mean())+(1-y.mean())*np.log1p(-y.mean())))
            for label,values in [('delta_auc',auc-ba),('logloss_gain',bl-ll),('brier_gain',bb-brier)]:
                row.update({label:values[0],label+'_lo':np.quantile(values[1:],.025),label+'_hi':np.quantile(values[1:],.975)})
            summaries.append(row)
        f=features[features.contrast.eq(con)&features.discovery.eq(disc)]
        for col in [c for c in f if c.endswith('_missing')]:
            missing.append(dict(contrast=con,discovery=disc,feature=col,missing=int(f[col].sum()),n=len(f)))
    save(pd.DataFrame(summaries),HERE/'cis_oof_gene_bootstrap.tsv')
    save(pd.DataFrame(missing),HERE/'cis_missingness.tsv')
    calibration=read(C/'analysis/D_calibration.tsv')
    calibration['raw_supported']=calibration.p_boot<.05
    calibration['joint_supported']=calibration.q_boot_32<.05
    save(calibration,HERE/'cis_calibration_diagnostics.tsv')


def continuous():
    nodes=read(F/'data/figure4-program-v9/all_mapped_nodes.tsv.gz')
    valid=nodes[nodes.measurement_pass&nodes.code.isin(['W','C','H1','H2','H3','H4'])]
    wide=valid.pivot(index=['layer','human_gene_id','code'],columns='node',values='position_index')
    originals=read(F/'data/figure4-state-recurrence-v25/state_specific_permutation.tsv')
    summaries=[];pairsout=[]
    for layer in ['supported50','expanded_primary']:
        for state,hcodes in [('R',['H1','H3']),('E',['H1','H2']),('P',['H1','H2','H3','H4'])]:
            frames=[]
            for m in ['W','C']:
                for h in hcodes:
                    mm=wide.xs((layer,m),level=('layer','code')).dropna(subset=['S',state])
                    hh=wide.xs((layer,h),level=('layer','code')).dropna(subset=['S',state])
                    ids=mm.index.intersection(hh.index)
                    frames.append(pd.DataFrame(dict(gene_id=ids,x=(mm.loc[ids,state]-mm.loc[ids,'S']).values,y=(hh.loc[ids,state]-hh.loc[ids,'S']).values,m=m,h=h)))
            union=sorted(set().union(*(set(f.gene_id) for f in frames)))
            arrays=[]
            for f in frames:
                a=f.set_index('gene_id')[['x','y']].reindex(union).to_numpy();arrays.append(a)
            obs=np.array([spearmanr(f.x,f.y).statistic for f in frames])
            old=originals[originals.layer.eq(layer)&originals.state.eq(state)&originals.metric.eq('spearman')].iloc[0]
            assert np.isclose(obs.mean(),old.observed,atol=1e-12)
            rng=rng_for('cross_species_gene_bootstrap',layer,state)
            boot=[]
            for _ in range(1000):
                ix=rng.integers(len(union),size=len(union));values=[]
                for a in arrays:
                    v=a[ix];v=v[np.isfinite(v).all(1)]
                    values.append(spearmanr(v[:,0],v[:,1]).statistic)
                boot.append(values)
            boot=np.array(boot)
            for i,f in enumerate(frames):
                pairsout.append(dict(layer=layer,state=state,mouse=f.m.iloc[0],human=f.h.iloc[0],n=len(f),rho=obs[i],lo=np.quantile(boot[:,i],.025),hi=np.quantile(boot[:,i],.975)))
            summaries.append(dict(layer=layer,state=state,union_genes=len(union),pairs=len(frames),rho=obs.mean(),
                lo=np.quantile(boot.mean(1),.025),hi=np.quantile(boot.mean(1),.975),min_pair_rho=obs.min(),max_pair_rho=obs.max(),
                positive_pairs=int((obs>0).sum()),original_p=old.p_greater,original_q=old.q_all_24))
            print('Continuous',layer,state,flush=True)
    save(pd.DataFrame(summaries),HERE/'continuous_gene_block_summary.tsv')
    save(pd.DataFrame(pairsout),HERE/'continuous_pair_summary.tsv')


def intervals():
    # All eligible pairs include opposite directions; a recurrent-only table
    # would condition on the direction agreement being evaluated.
    d=read(ROOT/'notes/2026-10-01-figure5-union-update/Figure.5B.display-region-pairs.tsv')
    original=read(ROOT/'notes/2026-10-01-figure5-union-update/Figure.5B.display-statistics.tsv')
    d['a']=d[['effect20_a','effect50_a']].mean(1);d['b']=d[['effect20_b','effect50_b']].mean(1)
    out=[]
    for (sp,con),f in d.groupby(['species','contrast']):
        groups=[b for _,b in f.groupby('gene_id')];n=len(groups)
        ref=original[original.species.eq(sp)&original.contrast.eq(con)].iloc[0]
        assert len(f)==ref.n and n==ref.genes
        stats=[]
        for b in groups:
            stats.append([(b.a*b.b>0).mean(),(b.a>0).mean(),(b.b>0).mean(),b.a.median(),b.b.median()])
        stats=np.array(stats)
        def direction(v):
            concordance,pa,pb=v[:,:3].mean(0)
            return concordance-(pa*pb+(1-pa)*(1-pb))
        rng=rng_for('interval_gene_bootstrap',sp,con)
        obs=spearmanr(f.a,f.b).statistic
        assert np.isclose(obs,ref.rho_mean)
        boot=[];loo=[]
        if n>=10:
            for _ in range(2000):
                ix=rng.integers(n,size=n)
                v=pd.concat([groups[i] for i in ix])
                boot.append([spearmanr(v.a,v.b).statistic,direction(stats[ix])])
            for i in range(n):
                v=pd.concat([b for j,b in enumerate(groups) if j!=i])
                loo.append(spearmanr(v.a,v.b).statistic)
        boot=np.array(boot)
        pa,pb=stats[:,1:3].mean(0)
        out.append(dict(species=sp,contrast=con,genes=n,regions=len(f),rho=obs,
            rho_lo=np.nanquantile(boot[:,0],.025) if n>=10 else np.nan,rho_hi=np.nanquantile(boot[:,0],.975) if n>=10 else np.nan,
            loo_rho_min=np.nanmin(loo) if loo else np.nan,loo_rho_max=np.nanmax(loo) if loo else np.nan,direction_agreement=stats[:,0].mean(),
            direction_expectation=pa*pb+(1-pa)*(1-pb),direction_excess=direction(stats),
            direction_lo=np.quantile(boot[:,1],.025) if n>=10 else np.nan,direction_hi=np.quantile(boot[:,1],.975) if n>=10 else np.nan,
            largest_gene_region_fraction=max(map(len,groups))/len(f),constant_direction_margins=bool(pa in (0,1) and pb in (0,1))))
    save(pd.DataFrame(out),HERE/'interval_gene_block_diagnostics.tsv')


def coverage():
    d=read(ROOT/'notes/2026-09-26-state-comparability-audit/metacell_expression_support.tsv.gz')
    d['label_pass']=d.same_label_vote>=.6;d['distance_pass']=d.cross_same_mean5_distance<=d.caliper_q95
    assert ((d.label_pass&d.distance_pass)==d.support_q95).all()
    d['reason']=np.select([d.support_q95,~d.label_pass&~d.distance_pass,~d.label_pass,d.label_pass&~d.distance_pass],['supported','both_fail','label_only_fail','distance_only_fail'],default='unknown')
    out=d.groupby(['species','study_id','node','calibration','reason']).size().rename('metacells').reset_index()
    save(out,HERE/'expression_support_failure_decomposition.tsv')
    rows=[]
    for p in sorted((ROOT/'notes/2026-10-02-reviewer-followup/native').glob('*/coverage.tsv.gz')):
        f=read(p);rows.append(f.groupby(['study_id','contrast','status']).size().rename('genes').reset_index())
    save(pd.concat(rows),HERE/'measurement_control_coverage.tsv')
    a=read(ROOT/'notes/2026-09-30-multimetric-statistics/inference/original_candidate_statistical_support.tsv')
    a['diagnosis']=np.select([~a.new_tested,a.equal_MC_WPP_estimate*a.direction<0,
        a.equal_MC_WPP_q_BH_combined>=.05,a.equal_MC_WPP_estimate.abs()<.02],
        ['uncovered','opposite_direction','q_not_passed','below_effect_threshold'],default='supported_direction_q_effect')
    save(a.groupby(['contrast','diagnosis']).size().rename('genes').reset_index(),HERE/'independent_well_candidate_failures.tsv')
    p=read(ROOT/'notes/2026-09-30-multimetric-statistics/pool_inference/gene_metric_tests.tsv.gz')
    save(p,HERE/'pool_design_resolution.tsv.gz')
    external=read(ROOT/'notes/2026-10-02-gse99922-pilot/vector_results/vector_summary.tsv')
    save(external,HERE/'external_full_design_grid.tsv')
    e=read(ROOT/'notes/2026-10-02-gse99922-pilot/vector_results/vector_candidate_status.tsv')
    save(e.groupby(['configuration','set_name','match_status']).size().rename('genes').reset_index(),HERE/'external_candidate_attrition.tsv')


if __name__=='__main__':
    for function in [sharing,cis,coverage,intervals,continuous]:
        function();print('DONE',function.__name__,flush=True)
        (HERE/'other_input_sha256.json').write_text(json.dumps(INPUTS,indent=2)+'\n')
