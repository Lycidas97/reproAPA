"""Validate legacy hits and run the frozen site/region summaries on the union."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
from scan_union import HERE, ROOT, UNION, F5, ROUTES, module, migration, paired, digest

OLD=ROOT/'notes/2026-09-29-core-sensitivity/figure5-paired-pm'

def replay():
    rows=[]
    keys=['route','gene_id','category','chrom','strand','start0','end0','contrast']
    exact=['names','feature_groups','site_types','mask_json','exact','variable']
    numeric=['pooled_lo','pooled_hi','reference_UMI','target_UMI']+[f'eq{t}_{x}' for t in [20,50] for x in ['lo','hi','n_reference','n_target']]
    for sp in ROUTES:
        old=pd.read_csv(OLD/f'full-cohort-v2/data/{sp}/site_effects.tsv.gz',sep='\t')
        old=old[old.contrast.isin(['SP','MP'])].set_index(keys).sort_index()
        new=pd.read_csv(HERE/f'full-cohort-v2/data/{sp}/site_effects.tsv.gz',sep='\t').set_index(keys).sort_index()
        assert old.index.is_unique and new.index.is_unique
        assert old.index.isin(new.index).all(),(sp,'old recognition intervals lost')
        same=new.loc[old.index]
        for c in exact:
            assert old[c].fillna('').equals(same[c].fillna('')),(sp,c)
        for c in numeric: np.testing.assert_allclose(old[c],same[c],rtol=1e-10,atol=1e-12,equal_nan=True,err_msg=sp+' '+c)
        rows.append(dict(species=sp,original_site_contrast_rows=len(old),union_site_contrast_rows=len(new),replay='PASS'))
    pd.DataFrame(rows).to_csv(HERE/'legacy_site_replay.tsv',sep='\t',index=False)
    print('Original recognition-site/effect replay PASS',rows,flush=True)

def regions():
    d=pd.read_csv(UNION/'annotation/interval_effects.tsv.gz',sep='\t')
    d['region_id']=d.route+'|'+d.gene_id+'|annotation:'+d.start0.astype(str)+'-'+d.end0.astype(str)
    d.to_csv(HERE/'all_coordinate_regions.paired_support.tsv.gz',sep='\t',index=False)
    paired_support=d.observed_both5.copy()
    for (sp,route),indices in d[d.contrast.eq('MP')].groupby(['species','route']).groups.items():
        models=json.loads((UNION/f'annotation/{sp}/{route}/models.json').read_text())
        cache={}
        for i in indices:
            row=d.loc[i]; gid=row.gene_id
            if gid not in cache: cache[gid]=migration.state_stats(models[gid])
            stats=cache[gid]; mask=np.array(json.loads(row['mask']))
            d.at[i,'observed_both5']=all(min(float(stats[50,s][2]@mask),float(stats[50,s][2]@(1-mask)))>=5 for s in 'MP')
    d.to_csv(HERE/'all_coordinate_regions.tsv.gz',sep='\t',index=False)
    d[d.observed_both5.ne(paired_support)].to_csv(HERE/'support_gate_differences.tsv',sep='\t',index=False)
    # Frozen original regions with finite effects must retain the same geometry,
    # effects, and observation gate within the original candidate pool.
    old=pd.read_csv(OLD/'all_coordinate_regions.paired_pm_original_support.tsv.gz',sep='\t')
    old=old[np.isfinite(old.effect20+old.effect50)]
    keys=['species','route','gene_id','contrast','start0','end0']
    check=old.set_index(keys); new=d.set_index(keys).loc[check.index]
    for c in ['effect20','effect50']: np.testing.assert_allclose(check[c],new[c],rtol=1e-10,atol=1e-12)
    assert np.array_equal(check.observed_both5,new.observed_both5)
    print('Original region replay PASS',len(check),'support differences',int(d.observed_both5.ne(paired_support).sum()),flush=True)

def summaries(out):
    # The frozen summarizer has only deterministic top-level table transforms.
    path=F5/'site-recurrence-v1/summarize.py'
    source=path.read_text().replace("OUT=ROOT/'analysis'",'OUT=Path('+repr(str(out))+')')
    exec(compile(source,str(path),'exec'),{'__file__':str(path),'__name__':'union_summary'})

def run_analysis(out,reg,consistency=True):
    out.mkdir(exist_ok=True)
    rec=module('union_recurrence_'+out.name,F5/'site-recurrence-v1/analyze.py')
    rec.P=HERE;rec.OUT=out;rec.REG=reg;rec.main()
    if consistency:
        con=module('union_consistency_'+out.name,F5/'site-recurrence-v1/consistency.py')
        con.ROOT=HERE;con.OUT=out;con.SOURCE=out/'exact_site_cross_study.tsv.gz';con.main()
    summaries(out)

def main():
    replay();regions()
    run_analysis(HERE/'analysis',HERE/'all_coordinate_regions.tsv.gz')
    run_analysis(HERE/'analysis-paired-support',HERE/'all_coordinate_regions.paired_support.tsv.gz',False)
    gate=[]
    for name in ['recurrence_union_across_contrasts.tsv','exact_partial_recurrence_summary.tsv']:
        a=pd.read_csv(HERE/'analysis'/name,sep='\t');b=pd.read_csv(HERE/'analysis-paired-support'/name,sep='\t')
        gate.append(dict(table=name,identical=a.equals(b)))
    cols=['species','contrast','gene_id','category','factor','chrom','strand','start0','end0','region_id_a','region_id_b','effect20_a','effect20_b','effect50_a','effect50_b']
    frames=[pd.read_csv(HERE/name/'exact_site_cross_study.tsv.gz',sep='\t') for name in ['analysis','analysis-paired-support']]
    picked=[d[d.both_observed&d.both_changed][cols].sort_values(cols[:9]).reset_index(drop=True) for d in frames]
    fractions=[pd.read_csv(HERE/name/'partial_retention_summary.tsv',sep='\t').partial_fraction for name in ['analysis','analysis-paired-support']]
    (HERE/'verification.json').write_text(json.dumps(dict(status='PASS',recognition_sites_recomputed=True,legacy_site_replay=True,legacy_region_replay=True,estimator='P-S separate state means; P-M common physical library paired differences',primary_support='original marginal qualified counts',paired_support_comparison=gate,paired_support_region_statistic_inputs_identical=picked[0].equals(picked[1]),partial_fraction_max_absolute_difference_pp=float((fractions[0]-fractions[1]).abs().max()*100),scope='DWUI union MPRO candidates; sequence predictions, no occupancy or regulatory activity inference'),indent=2)+'\n')

if __name__=='__main__':main()
