"""Exact site recurrence and region inference with newly computed denominators."""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr

from downstream_runtime import Context,bh,seed
from downstream_regulatory import F5,ROUTES,observed_regions

KEY=['species','gene_id','category','factor','chrom','strand','start0','end0','contrast']
SITE_COLUMNS=['species','route','gene_id','gene_symbol','category','factor','chrom','strand','start0','end0','contrast','names','provenance','mask','status',
              'region_id','region_start0','region_end0','effect20','effect50','quantified','observed','changed','strict_changed','retained','evaluable']


def site_inventory(data,regions):
    lookup={k:g.sort_values('start0') for k,g in regions.groupby(['species','route','gene_id','contrast'])};rows=[]
    for key,g in data.groupby(['species','route','gene_id','contrast'],sort=False):
        rr=lookup.get(key);starts=rr.start0.to_numpy() if rr is not None else np.array([])
        for x in g.itertuples():
            mask=np.asarray(json.loads(x.mask_json));exact=bool(x.exact);finite=np.isfinite(x.eq20_lo+x.eq50_lo)
            constant=exact and np.all(mask==1);excluded=exact and np.all(mask==0);region=None
            ix=np.searchsorted(starts,x.start0,side='right')-1
            if ix>=0:
                candidate=rr.iloc[ix]
                if x.end0<=candidate.end0 and exact and np.array_equal(mask[:,0],json.loads(candidate['mask'])):
                    region=candidate
                    for t in [20,50]:np.testing.assert_allclose(getattr(x,f'eq{t}_lo'),candidate[f'effect{t}'],atol=1e-12,equal_nan=True)
            observed=region is not None and bool(region.observed_both5);quantified=finite and exact and region is not None
            change=quantified and min(abs(x.eq20_lo),abs(x.eq50_lo))>=.1 and x.eq20_lo*x.eq50_lo>0
            status=('ambiguous' if not exact else 'constitutive_retained' if constant and finite else 'constitutive_excluded' if excluded else
                    'coverage_insufficient' if not finite else 'boundary_unresolved' if region is None else 'variable_observed' if observed else 'variable_low_observation')
            identity=x.feature_groups if x.category=='miRNA' else x.names
            for factor in sorted(set(str(identity).split(';'))):
                rows.append(dict(species=x.species,route=x.route,gene_id=x.gene_id,gene_symbol=x.gene_symbol,category='miRNA' if x.category=='miRNA' else 'RBP',factor=factor,
                    chrom=x.chrom,strand=x.strand,start0=x.start0,end0=x.end0,contrast=x.contrast,names=x.names,provenance=x.category,mask=x.mask_json,status=status,
                    region_id=region.region_id if region is not None else '',region_start0=region.start0 if region is not None else np.nan,
                    region_end0=region.end0 if region is not None else np.nan,effect20=x.eq20_lo,effect50=x.eq50_lo,quantified=quantified,observed=observed,
                    changed=change,strict_changed=change and observed,retained=constant and finite,evaluable=(constant and finite) or (quantified and observed)))
    f=pd.DataFrame(rows,columns=SITE_COLUMNS)
    for col in ['quantified','observed','changed','strict_changed','retained','evaluable']:f[col]=f[col].astype(bool)
    keys=['route']+KEY
    provenance=f.groupby(keys).provenance.agg(lambda v:';'.join(sorted(set(v))))
    f=f.drop_duplicates(keys).set_index(keys);f['provenance']=provenance
    return f.reset_index()


def partial(ctx,sites,recurrent,math):
    keys=['species','route','contrast','category','gene_id','factor']
    f=sites.assign(positive_change=sites.strict_changed&sites.effect20.gt(0),negative_change=sites.strict_changed&sites.effect20.lt(0))
    f['coord']=f.chrom+':'+f.start0.astype(str)+'-'+f.end0.astype(str)
    pair=f.groupby(keys).agg(gene_symbol=('gene_symbol','first'),n_sites=('start0','size'),n_evaluable=('evaluable','sum'),
        n_changed=('strict_changed','sum'),n_retained=('retained','sum'),n_positive=('positive_change','sum'),n_negative=('negative_change','sum'))
    pair['changed_sites']=f[f.strict_changed].groupby(keys).coord.agg(';'.join)
    pair['retained_sites']=f[f.retained].groupby(keys).coord.agg(';'.join);pair[['changed_sites','retained_sites']]=pair[['changed_sites','retained_sites']].fillna('')
    pair['n_unresolved']=pair.n_sites-pair.n_evaluable;pair['has_partial']=pair.n_changed.gt(0)&pair.n_retained.gt(0)
    pair['classification']=[math.classify([1]*r.n_positive+[-1]*r.n_negative,r.n_retained,r.n_evaluable,r.n_sites)[0] for r in pair.itertuples()]
    pair=pair.reset_index();summary=[]
    for sp in ROUTES:
        for route in [r for s,r in ROUTES[sp]]:
            for co in ['SP','MP']:
                for cat in ['miRNA','RBP']:
                    q=pair[pair.species.eq(sp)&pair.route.eq(route)&pair.contrast.eq(co)&pair.category.eq(cat)]
                    eligible=q[q.n_evaluable>=2];active=eligible[eligible.n_changed>0]
                    blocks=active.groupby('gene_id').agg(success=('has_partial','sum'),total=('has_partial','size'))
                    ci=math.bootstrap_ratio(blocks,np.random.default_rng(seed('partial',sp,route,co,cat)),4000)
                    summary.append(dict(species=sp,route=route,contrast=co,category=cat,multisite_pairs=len(eligible),changed_multisite_pairs=len(active),
                        changed_genes=active.gene_id.nunique(),partial_including_mixed=int(active.has_partial.sum()),partial_genes=active[active.has_partial].gene_id.nunique(),
                        partial_fraction=active.has_partial.mean(),ci_low=ci[0],ci_high=ci[1],
                        **{c:int(active.classification.eq(c).sum()) for c in ['partial','mixed','all_annotated_sites_changed','changed_with_other_or_unresolved']}))
    cols=['species','contrast','category','gene_id','factor'];r=recurrent[recurrent.recurrent_strict].copy()
    r['coord']=r.chrom+':'+r.start0.astype(str)+'-'+r.end0.astype(str);lookup=r.groupby(cols).coord.agg(set).to_dict();repeated=[]
    for sp,routes in ROUTES.items():
        a,b=[r for s,r in routes]
        q=pair[pair.route.eq(a)&pair.has_partial].merge(pair[pair.route.eq(b)&pair.has_partial],on=cols,suffixes=('_a','_b'))
        for x in q.to_dict('records'):
            change=set(x['changed_sites_a'].split(';'))&set(x['changed_sites_b'].split(';'))&lookup.get(tuple(x[k] for k in cols),set())
            retained=(set(x['retained_sites_a'].split(';'))&set(x['retained_sites_b'].split(';')))-{''}
            repeated.append(x|dict(shared_recurrent_changed_sites=';'.join(sorted(change)),shared_retained_sites=';'.join(sorted(retained)),exact_partial_recurrence=bool(change and retained)))
    ctx.table(pair,'gene_factor_classification.tsv.gz');ctx.table(summary,'partial_retention_summary.tsv')
    frame=pd.DataFrame(repeated,columns=list(q.columns)+['shared_recurrent_changed_sites','shared_retained_sites','exact_partial_recurrence'])
    ctx.table(frame,'cross_study_partial_retention.tsv')
    exact=frame[frame.exact_partial_recurrence.astype(bool)]
    ctx.table(exact.groupby(['species','contrast','category']).agg(exact_partial_pairs=('gene_id','size'),genes=('gene_id','nunique')).reset_index(),'exact_partial_recurrence_summary.tsv')


def consistency(ctx,rec):
    num=ctx.functions(F5+'site-recurrence-v1/consistency.py',['safe_spearman','stats_frame','block_bootstrap','gene_label_test'],
        dict(rankdata=rankdata,spearmanr=spearmanr,RNG_SEED=seed('region')%(2**32)))
    x=rec[rec.both_observed&rec.both_changed].copy()
    x['effect_mean_a']=(x.effect20_a+x.effect50_a)/2;x['effect_mean_b']=(x.effect20_b+x.effect50_b)/2
    keys=['species','contrast','gene_id','gene_symbol_a','chrom','strand','region_id_a','region_id_b']
    checks=x.groupby(keys)[['effect20_a','effect50_a','effect20_b','effect50_b']].nunique();assert checks.le(1).all().all()
    cat=x.groupby(keys).category.agg(lambda z:';'.join(sorted(set(z)))).rename('categories')
    d=x.drop_duplicates(keys).set_index(keys).join(cat).reset_index();ctx.table(d,'all_changed_region_pairs.tsv.gz')
    rows=[];equal=[];loo=[]
    for sp in ROUTES:
        for co in ['SP','MP']:
            for category in ['all','miRNA','RBP']:
                source=x[x.species.eq(sp)&x.contrast.eq(co)]
                base=d[d.species.eq(sp)&d.contrast.eq(co)]
                if category!='all':source=source[source.category.eq(category)];base=base[base.categories.str.split(';').apply(lambda v:category in v)]
                for unit in ['region_pair','mask_pair']:
                    z=base if unit=='region_pair' else source.drop_duplicates(['gene_id','mask_a','mask_b'])
                    key=dict(species=sp,contrast=co,category=category,unit=unit);rngseed=seed('region',*key.values())
                    stats=num.stats_frame(z);n=z.gene_id.nunique();row=key|stats
                    if n:
                        row.update(num.block_bootstrap(z,5000,rngseed))
                        gm=z.groupby('gene_id')[['effect_mean_a','effect_mean_b']].median()
                        if len(gm)>=3 and gm.nunique().min()>1:row.update(num.gene_label_test(z,19999,rngseed))
                        else:row.update(gene_median_rho=np.nan,gene_label_p=np.nan,permutations=19999)
                    else:row.update(gene_median_rho=np.nan,gene_label_p=np.nan,permutations=19999)
                    row['status']='inference' if n>=3 and np.isfinite(row['gene_label_p']) else 'descriptive_only'
                    row['p_for_family']=row['gene_label_p'] if row['status']=='inference' else 1.;rows.append(row)
                    if not n:continue
                    ge=z.assign(same=np.sign(z.effect20_a)==np.sign(z.effect20_b),pa=z.effect20_a.gt(0),pb=z.effect20_b.gt(0)).groupby('gene_id')[['same','pa','pb']].mean()
                    values=ge.to_numpy();r=np.random.default_rng(rngseed);boot=values[r.integers(n,size=(2000,n))].mean(1)
                    c0=boot[:,1]*boot[:,2]+(1-boot[:,1])*(1-boot[:,2]);mean=values.mean(0)
                    equal.append(key|dict(genes=n,Cobs=mean[0],C0=mean[1]*mean[2]+(1-mean[1])*(1-mean[2]),
                        Cobs_lo=np.quantile(boot[:,0],.025),Cobs_hi=np.quantile(boot[:,0],.975),delta_C_lo=np.quantile(boot[:,0]-c0,.025),delta_C_hi=np.quantile(boot[:,0]-c0,.975)))
                    for gene in ge.index:
                        sub=z[z.gene_id.ne(gene)];loo.append(key|dict(omitted_gene=gene,remaining_genes=sub.gene_id.nunique(),**num.stats_frame(sub)))
    result=pd.DataFrame(rows);result['primary_test']=result.unit.eq('region_pair')&result.category.eq('all')
    result['gene_label_q_primary']=np.nan
    # Separate-state DWUI design has four primary tests. Union paired-PM results
    # retain nominal P; broad correction remains explicitly a sensitivity family.
    if ctx.node=='original_pool_sites':
        ix=result.primary_test;result.loc[ix,'gene_label_q_primary']=bh(result.loc[ix,'p_for_family'])
    result['gene_label_q_all_rows']=bh(result.p_for_family)
    ctx.table(result,'region_consistency_statistics.tsv');ctx.table(equal,'gene_equal_consistency.tsv',columns=['species','contrast','category','unit','genes','Cobs','C0'])
    ctx.table(loo,'region_leave_one_gene_out.tsv',columns=['species','contrast','category','unit','omitted_gene','remaining_genes'])


def run(ctx):
    regions=ctx.read(ctx.analysis/'intervals/all_coordinate_regions.tsv.gz');frames=[]
    original=ctx.node=='original_pool_sites'
    claims=ctx.read(ctx.analysis/'dwui_pool/claim_index.tsv') if original else None
    for sp in ROUTES:
        f=ctx.read(ctx.analysis/f'site_scan/full-cohort-v2/data/{sp}/site_effects.tsv.gz')
        if original:
            f=f.merge(claims[['species','gene_id','contrast']].drop_duplicates(),on=['species','gene_id','contrast'],validate='many_to_one')
            for t in [20,50]:
                for col in ['lo','hi','n_reference','n_target']:f[f'eq{t}_{col}']=f[f'separate_eq{t}_{col}']
        frames.append(f)
    if original:
        regions=regions.merge(claims[['species','gene_id','contrast']].drop_duplicates(),on=['species','gene_id','contrast'],validate='many_to_one')
    regions=observed_regions(regions,common=not original)
    sites=site_inventory(pd.concat(frames,ignore_index=True),regions);ctx.table(sites,'site_factor_inventory.tsv.gz')
    math=ctx.functions(F5+'site-recurrence-v1/analyze.py',['classify','bootstrap_ratio'])
    recs=[];summaries=[];tests=[]
    for sp,routes in ROUTES.items():
        a,b=[r for s,r in routes];r=sites[sites.route.eq(a)].merge(sites[sites.route.eq(b)],on=KEY,suffixes=('_a','_b'),validate='one_to_one')
        r['both_quantified']=r.quantified_a&r.quantified_b;r['both_observed']=r.both_quantified&r.observed_a&r.observed_b
        r['both_changed']=r.changed_a&r.changed_b;r['same_direction']=r.effect20_a*r.effect20_b>0
        r['recurrent']=r.both_changed&r.same_direction;r['recurrent_strict']=r.recurrent&r.both_observed;recs.append(r)
        for co in ['SP','MP']:
            for cat in ['miRNA','RBP']:
                allg=r[r.contrast.eq(co)&r.category.eq(cat)]
                for scope,col in [('quantified','both_quantified'),('observed','both_observed')]:
                    g=allg[allg[col]];phy=g.drop_duplicates(['gene_id','chrom','strand','start0','end0']);yes=phy[phy.recurrent];both=phy[phy.both_changed]
                    blocks=phy.groupby('gene_id').agg(success=('recurrent','sum'),total=('recurrent','size'))
                    ci=math.bootstrap_ratio(blocks,np.random.default_rng(seed('site_ratio',sp,co,cat,scope)),4000)
                    summaries.append(dict(species=sp,contrast=co,category=cat,scope=scope,shared_site_factor_links=len(g),shared_physical_sites=len(phy),shared_genes=phy.gene_id.nunique(),
                        recurrent_links=int(g.recurrent.sum()),recurrent_sites=len(yes),recurrent_genes=yes.gene_id.nunique(),recurrent_fraction=len(yes)/len(phy) if len(phy) else np.nan,
                        ci_low=ci[0],ci_high=ci[1],both_changed_sites=len(both),opposite_sites=int((both.effect20_a*both.effect20_b<0).sum()),
                        concordance=len(yes)/len(both) if len(both) else np.nan,gains=int(yes.effect20_a.gt(0).sum()),losses=int(yes.effect20_a.lt(0).sum())))
                    if scope!='observed':continue
                    block=both.assign(sa=np.sign(both.effect20_a),sb=np.sign(both.effect20_b)).groupby('gene_id').agg(sa=('sa','mean'),sb=('sb','mean'),same=('same_direction','mean'))
                    va,vb=block.sa.to_numpy(),block.sb.to_numpy();observed=float(np.mean(va*vb)) if len(block) else np.nan
                    random=np.random.default_rng(seed('site_sign',sp,co,cat));null=np.array([np.mean(va*random.permutation(vb)) for i in range(9999)]) if len(block)>=3 else np.array([])
                    tests.append(dict(species=sp,contrast=co,category=cat,genes=len(block),gene_equal_site_concordance=block.same.mean(),paired_gene_direction_score=observed,
                        shuffled_expectation=float(va.mean()*vb.mean()) if len(block) else np.nan,p=(1+(null>=observed-1e-12).sum())/10000 if len(null) else 1.,
                        status='inference' if len(null) else 'fewer_than_3_genes_descriptive_only'))
    rec=pd.concat(recs,ignore_index=True);ctx.table(rec,'exact_site_cross_study.tsv.gz');ctx.table(summaries,'recurrence_summary.tsv')
    tests=pd.DataFrame(tests);tests['q_8']=bh(tests.p);ctx.table(tests,'gene_block_direction_tests.tsv')
    partial(ctx,sites,rec,math);consistency(ctx,rec)
    recurrent=rec[rec.recurrent_strict];ctx.table(recurrent,'recurrent_sites.tsv.gz')
    union=recurrent.drop_duplicates(['species','category','gene_id','chrom','strand','start0','end0']).groupby(['species','category']).agg(unique_physical_sites=('gene_id','size'),genes=('gene_id','nunique')).reset_index()
    ctx.table(union,'recurrence_union_across_contrasts.tsv')
    candidates=recurrent[recurrent.category.eq('RBP')].groupby(['species','factor']).agg(recurrent_genes=('gene_id','nunique'),site_contrast_links=('start0','size')).reset_index()
    ctx.table(candidates,'recurrent_rbp_candidates.tsv')
    ctx.finish(dict(new_site_inventory=True,physical_sites_deduplicated=True,new_evaluability_and_partial_denominators=True,
        gene_block_inference=True,small_samples_explicit=True,planned_families_preserved=True),pool='new_DWUI' if original else 'new_DWUI_MPRO_union',
        estimand='separate-state means' if original else 'SP separate-state; MP common physical libraries',site_ratio_bootstrap=4000,region_bootstrap=5000,gene_equal_bootstrap=2000)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--node',choices=['site_statistics','original_pool_sites'],required=True)
    a=p.parse_args();run(Context(a.config,a.node,a.out))
