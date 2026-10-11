"""Region-level directional and continuous cross-study correspondence."""
import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis'
SOURCE=OUT/'exact_site_cross_study.tsv.gz'
RNG_SEED=20260920


def safe_spearman(a,b):
    a=np.asarray(a); b=np.asarray(b)
    if len(a)<2 or np.ptp(a)==0 or np.ptp(b)==0:
        return np.nan
    return float(spearmanr(a,b).statistic)


def stats_frame(d):
    same=np.sign(d.effect20_a)==np.sign(d.effect20_b)
    pa=float((d.effect20_a>0).mean()); pb=float((d.effect20_b>0).mean())
    cobs=float(same.mean()); c0=pa*pb+(1-pa)*(1-pb)
    out=dict(n=len(d),genes=d.gene_id.nunique(),cobs=cobs,
        p_a_positive=pa,p_a_negative=1-pa,p_b_positive=pb,p_b_negative=1-pb,
        c0=c0,delta_c=cobs-c0)
    for label,a,b in [('eq20','effect20_a','effect20_b'),('eq50','effect50_a','effect50_b'),
                      ('mean','effect_mean_a','effect_mean_b')]:
        out['rho_'+label]=safe_spearman(d[a],d[b])
    return out


def block_bootstrap(d,n=5000,seed=RNG_SEED):
    rng=np.random.default_rng(seed)
    cols=['effect20_a','effect20_b','effect50_a','effect50_b','effect_mean_a','effect_mean_b']
    blocks=[g[cols].to_numpy() for _,g in d.groupby('gene_id',sort=True)]
    vals=np.empty((n,6),dtype=float)
    for _ in range(n):
        z=np.concatenate([blocks[i] for i in rng.integers(0,len(blocks),len(blocks))])
        pa=np.mean(z[:,0]>0); pb=np.mean(z[:,1]>0)
        cobs=np.mean(np.sign(z[:,0])==np.sign(z[:,1])); c0=pa*pb+(1-pa)*(1-pb)
        vals[_,:]=[cobs,c0,cobs-c0,safe_spearman(z[:,0],z[:,1]),
                    safe_spearman(z[:,2],z[:,3]),safe_spearman(z[:,4],z[:,5])]
    names=['cobs','c0','delta_c','rho_eq20','rho_eq50','rho_mean']
    return {f'{name}_{tail}':float(np.nanquantile(vals[:,i],q)) for i,name in enumerate(names)
            for tail,q in [('ci_low',.025),('ci_high',.975)]}


def gene_label_test(d,n=19999,seed=RNG_SEED):
    """Test gene identity using one robust effect summary per gene."""
    g=d.groupby('gene_id').agg(a=('effect_mean_a','median'),b=('effect_mean_b','median')).dropna()
    if len(g)>2:
        # Spearman correlation is Pearson correlation of ranks.  Standardizing
        # the ranks once makes the permutation test exact and avoids repeatedly
        # ranking the same marginal distributions.
        a=rankdata(g.a.to_numpy()); b=rankdata(g.b.to_numpy())
        a=(a-a.mean())/np.linalg.norm(a-a.mean())
        b=(b-b.mean())/np.linalg.norm(b-b.mean())
        observed=float(a@b)
        rng=np.random.default_rng(seed)
        perm=np.fromiter((a@rng.permutation(b) for _ in range(n)),dtype=float,count=n)
    else:
        observed=np.nan
        perm=np.array([])
    p=float((1+np.count_nonzero(np.abs(perm)>=abs(observed)-1e-12))/(n+1)) if len(perm) else np.nan
    return dict(gene_median_rho=observed,gene_label_p=p,permutations=n,
        null='shuffle study-B gene median effects; preserves both marginal effect distributions')


def bh(p):
    p=np.asarray(p,dtype=float)
    order=np.argsort(p); adjusted=np.empty(len(p))
    adjusted[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
    return adjusted


class Tests(unittest.TestCase):
    def test_global_bias_reference(self):
        d=pd.DataFrame(dict(effect20_a=[-1,-1,1,1],effect20_b=[-2,2,-2,2],
            effect50_a=[-1,-1,1,1],effect50_b=[-2,2,-2,2],effect_mean_a=[-1,-1,1,1],
            effect_mean_b=[-2,2,-2,2],gene_id=list('abcd')))
        s=stats_frame(d)
        self.assertEqual(s['cobs'],.5)
        self.assertEqual(s['c0'],.5)
        self.assertEqual(s['delta_c'],0)
    def test_perfect_correspondence(self):
        d=pd.DataFrame(dict(effect20_a=[-4,-2,1,3],effect20_b=[-8,-4,2,6],
            effect50_a=[-4,-2,1,3],effect50_b=[-8,-4,2,6],effect_mean_a=[-4,-2,1,3],
            effect_mean_b=[-8,-4,2,6],gene_id=list('abcd')))
        self.assertAlmostEqual(stats_frame(d)['rho_mean'],1.)
    def test_gene_block_resampling_keeps_whole_blocks(self):
        d=pd.DataFrame(dict(gene_id=['a','a','b'],effect20_a=[-1,-2,1],effect20_b=[-1,-2,1],
            effect50_a=[-1,-2,1],effect50_b=[-1,-2,1],effect_mean_a=[-1,-2,1],effect_mean_b=[-1,-2,1]))
        ci=block_bootstrap(d,n=20)
        self.assertEqual(ci['cobs_ci_low'],1.)
        self.assertEqual(ci['cobs_ci_high'],1.)


def main():
    t=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    assert t.wasSuccessful() and not t.skipped
    x=pd.read_csv(SOURCE,sep='\t',low_memory=False)
    x=x[x.both_observed & x.both_changed].copy()
    x['effect_mean_a']=(x.effect20_a+x.effect50_a)/2
    x['effect_mean_b']=(x.effect20_b+x.effect50_b)/2
    pairkeys=['species','contrast','gene_id','gene_symbol_a','chrom','strand','region_id_a','region_id_b']
    checks=x.groupby(pairkeys).agg(n_a20=('effect20_a','nunique'),n_a50=('effect50_a','nunique'),
        n_b20=('effect20_b','nunique'),n_b50=('effect50_b','nunique'))
    assert checks.max().max()==1
    categories=x.groupby(pairkeys).category.agg(lambda z:';'.join(sorted(set(z)))).rename('categories')
    factors=x.groupby(pairkeys).factor.agg(lambda z:';'.join(sorted(set(z)))).rename('factors')
    d=x.drop_duplicates(pairkeys).set_index(pairkeys).join([categories,factors]).reset_index()
    d['same_direction']=np.sign(d.effect20_a)==np.sign(d.effect20_b)
    d['effect_mean_a']=(d.effect20_a+d.effect50_a)/2
    d['effect_mean_b']=(d.effect20_b+d.effect50_b)/2
    d.to_csv(OUT/'all_changed_region_pairs.tsv.gz',sep='\t',index=False)
    results=[]
    detail=[]
    for (sp,contrast),base in d.groupby(['species','contrast']):
        for category in ['all','miRNA','RBP']:
            q=base if category=='all' else base[base.categories.str.split(';').apply(lambda z:category in z)]
            for unit in ['region_pair','mask_pair']:
                if unit=='region_pair':
                    z=q.copy()
                else:
                    # Effects are determined by masks; collapse repeated physical
                    # regions carrying the same paired PAS-retention structure.
                    source=x[(x.species==sp)&(x.contrast==contrast)]
                    source=source if category=='all' else source[source.category.eq(category)]
                    mk=['gene_id','mask_a','mask_b']
                    z=source.drop_duplicates(mk).copy()
                    z['effect_mean_a']=(z.effect20_a+z.effect50_a)/2
                    z['effect_mean_b']=(z.effect20_b+z.effect50_b)/2
                seed=RNG_SEED+sum(map(ord,sp+contrast+category+unit))
                result=dict(species=sp,contrast=contrast,category=category,unit=unit,
                    **stats_frame(z),**block_bootstrap(z,seed=seed),**gene_label_test(z,seed=seed))
                results.append(result)
                for r in z.itertuples():
                    detail.append(dict(species=sp,contrast=contrast,category=category,unit=unit,
                        gene_id=r.gene_id,effect20_a=r.effect20_a,effect50_a=r.effect50_a,
                        effect20_b=r.effect20_b,effect50_b=r.effect50_b,
                        effect_mean_a=r.effect_mean_a,effect_mean_b=r.effect_mean_b))
    result=pd.DataFrame(results)
    # The four category-agnostic physical-region tests are primary. Category
    # strata overlap, and mask-pair rows are sensitivity analyses, so including
    # all 24 rows in the primary FDR family would count the same APA effects
    # repeatedly.  A second all-row adjustment is retained for transparency.
    result['primary_test']=result.unit.eq('region_pair') & result.category.eq('all')
    result['gene_label_q_primary']=np.nan
    primary=result.primary_test
    result.loc[primary,'gene_label_q_primary']=bh(result.loc[primary,'gene_label_p'])
    result['gene_label_q_all_rows']=bh(result.gene_label_p)
    result.to_csv(OUT/'region_consistency_statistics.tsv',sep='\t',index=False)
    pd.DataFrame(detail).to_csv(OUT/'region_consistency_plot_data.tsv.gz',sep='\t',index=False)
    (OUT/'consistency_verification.json').write_text(json.dumps(dict(status='PASS',unit_tests=t.testsRun,
        skipped=0,source_rows=len(x),unique_region_pairs=len(d),results=len(result),seed=RNG_SEED),indent=2))
    print(result[result.unit.eq('region_pair')].to_string(index=False))


if __name__=='__main__':main()
