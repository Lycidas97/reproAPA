"""Exact-coordinate recurrence and within-pair partial retention, without imputation."""
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
P=ROOT.parent
OUT=ROOT/'analysis'
REG=P/'cis-regulatory-assessment-v1/analysis/all_coordinate_regions.tsv.gz'
KEY=['species','gene_id','category','factor','chrom','strand','start0','end0','contrast']


def changed(a,b):
    return np.isfinite(a+b) and min(abs(a),abs(b))>=.1 and a*b>0


def classify(changes,retained,n_eval,n_all):
    signs=set(changes)
    partial=bool(signs and retained)
    if n_eval<2:return 'insufficient_evaluable_sites',partial
    if len(signs)>1:return 'mixed',partial
    if partial:return 'partial',partial
    if signs and len(changes)==n_eval==n_all:return 'all_annotated_sites_changed',partial
    if signs:return 'changed_with_other_or_unresolved',partial
    return 'no_measured_change',partial


class Tests(unittest.TestCase):
    def test_constant_plus_change_is_partial(self):
        self.assertEqual(classify([-1],1,2,2),('partial',True))
    def test_missing_prevents_complete(self):
        self.assertEqual(classify([-1,-1],0,2,3)[0],'changed_with_other_or_unresolved')
    def test_mixed_retains_partial_flag(self):
        self.assertEqual(classify([-1,1],1,3,3),('mixed',True))
    def test_no_observation_is_not_change(self):
        self.assertFalse(changed(np.nan,.5))
        self.assertFalse(changed(.2,-.3))


def bootstrap_ratio(g,rng,n=4000):
    if not len(g):return [np.nan,np.nan]
    vals=g[['success','total']].to_numpy(float)
    sums=vals[rng.integers(0,len(vals),size=(n,len(vals)))].sum(axis=1)
    ratios=np.divide(sums[:,0],sums[:,1],out=np.full(n,np.nan),where=sums[:,1]>0)
    return np.nanquantile(ratios,[.025,.975]).tolist()


def main():
    t=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    assert t.wasSuccessful() and not t.skipped
    OUT.mkdir(exist_ok=True)
    regs=pd.read_csv(REG,sep='\t')
    rlookup={k:g.sort_values('start0') for k,g in regs.groupby(['species','route','gene_id','contrast'])}
    sources=[Path(__file__),REG]; parts=[]
    for sp in ['mouse','human']:
        path=P/f'full-cohort-v2/data/{sp}/site_effects.tsv.gz';sources.append(path)
        d=pd.read_csv(path,sep='\t',low_memory=False)
        d=d[d.contrast.isin(['SP','MP'])].copy()
        rows=[]
        for k,g in d.groupby(['species','route','gene_id','contrast'],sort=False):
            rr=rlookup.get(k)
            starts=rr.start0.to_numpy() if rr is not None else np.array([])
            for x in g.itertuples():
                mask=np.array(json.loads(x.mask_json)); exact=bool(x.exact)
                finite=np.isfinite(x.eq20_lo+x.eq50_lo)
                constant=exact and np.all(mask==1)
                excluded=exact and np.all(mask==0)
                region=None
                ix=np.searchsorted(starts,x.start0,side='right')-1
                if ix>=0:
                    candidate=rr.iloc[ix]
                    if x.end0<=candidate.end0 and exact and np.array_equal(mask[:,0],json.loads(candidate['mask'])):
                        region=candidate
                        for threshold in [20,50]:
                            np.testing.assert_allclose(getattr(x,f'eq{threshold}_lo'),candidate[f'effect{threshold}'],atol=1e-12,equal_nan=True)
                observed=region is not None and bool(region.observed_both5)
                quantified=finite and exact and region is not None
                change=quantified and changed(x.eq20_lo,x.eq50_lo)
                status=('ambiguous' if not exact else 'constitutive_retained' if constant and finite else
                        'constitutive_excluded' if excluded else 'coverage_insufficient' if not finite else
                        'boundary_unresolved' if region is None else 'variable_observed' if observed else 'variable_low_observation')
                identity=x.feature_groups if x.category=='miRNA' else x.names
                for factor in sorted(set(str(identity).split(';'))):
                    rows.append(dict(species=sp,route=x.route,gene_id=x.gene_id,gene_symbol=x.gene_symbol,
                        category='miRNA' if x.category=='miRNA' else 'RBP',factor=factor,
                        chrom=x.chrom,strand=x.strand,start0=x.start0,end0=x.end0,contrast=x.contrast,
                        names=x.names,provenance=x.category,mask=x.mask_json,status=status,
                        region_id=region.region_id if region is not None else '',
                        region_start0=region.start0 if region is not None else np.nan,
                        region_end0=region.end0 if region is not None else np.nan,
                        effect20=x.eq20_lo,effect50=x.eq50_lo,quantified=quantified,
                        observed=observed,changed=change,strict_changed=change and observed,
                        retained=constant and finite,evaluable=(constant and finite) or (quantified and observed)))
        s=pd.DataFrame(rows)
        groupkey=['route']+KEY
        provenance=s.groupby(groupkey).provenance.agg(lambda x:';'.join(sorted(set(x))))
        # Direct and orthology-transferred motifs for the same factor and exact
        # span are one relationship, not independent regulatory observations.
        s=s.drop_duplicates(groupkey).set_index(groupkey)
        s['provenance']=provenance
        parts.append(s.reset_index())
        print(sp,'site-factor rows',len(s),flush=True)
    sites=pd.concat(parts,ignore_index=True)
    sites.to_csv(OUT/'site_factor_inventory.tsv.gz',sep='\t',index=False)
    recs=[]; summaries=[]; tests=[]; rng=np.random.default_rng(20260919)
    for sp,a,b in [('mouse','CIP','WOT'),('human','GSE118258','GSE150311')]:
        sa=sites[sites.route.eq(a)]; sb=sites[sites.route.eq(b)]
        r=sa.merge(sb,on=KEY,suffixes=('_a','_b'),validate='one_to_one')
        r['both_quantified']=r.quantified_a&r.quantified_b
        r['both_observed']=r.both_quantified&r.observed_a&r.observed_b
        r['both_changed']=r.changed_a&r.changed_b
        r['same_direction']=r.effect20_a*r.effect20_b>0
        r['recurrent']=r.both_changed&r.same_direction
        r['recurrent_strict']=r.recurrent&r.both_observed
        recs.append(r)
        for (contrast,cat),allg in r.groupby(['contrast','category']):
            for scope,col in [('quantified','both_quantified'),('observed','both_observed')]:
                g=allg[allg[col]].copy()
                # Physical sites are counted once per category, even if they
                # match several factors. Factor identities remain in the table.
                phy=g.drop_duplicates(['gene_id','chrom','strand','start0','end0'])
                yes=phy[phy.recurrent]; both=phy[phy.both_changed]
                block=phy.groupby('gene_id').agg(success=('recurrent','sum'),total=('recurrent','size'))
                ci=bootstrap_ratio(block,rng)
                summaries.append(dict(species=sp,contrast=contrast,category=cat,scope=scope,
                    shared_site_factor_links=len(g),shared_physical_sites=len(phy),shared_genes=phy.gene_id.nunique(),
                    recurrent_links=int(g.recurrent.sum()),recurrent_sites=len(yes),recurrent_genes=yes.gene_id.nunique(),
                    recurrent_fraction=len(yes)/len(phy) if len(phy) else np.nan,ci_low=ci[0],ci_high=ci[1],
                    both_changed_sites=len(both),opposite_sites=int((both.effect20_a*both.effect20_b<0).sum()),
                    concordance=len(yes)/len(both) if len(both) else np.nan,
                    gains=int((yes.effect20_a>0).sum()),losses=int((yes.effect20_a<0).sum())))
                if scope=='observed' and len(both):
                    block=both.assign(sa=np.sign(both.effect20_a),sb=np.sign(both.effect20_b)).groupby('gene_id').agg(sa=('sa','mean'),sb=('sb','mean'),same=('same_direction','mean'))
                    va=block.sa.to_numpy(); vb=block.sb.to_numpy()
                    obs=float(np.mean(va*vb)); expected=float(va.mean()*vb.mean())
                    sim=np.array([np.mean(va*rng.permutation(vb)) for _ in range(9999)])
                    p=(1+np.count_nonzero(sim>=obs-1e-12))/10000
                    tests.append(dict(species=sp,contrast=contrast,category=cat,genes=len(block),
                        gene_equal_site_concordance=block.same.mean(),paired_gene_direction_score=obs,
                        shuffled_expectation=expected,excess_score=obs-expected,p=p,
                        null='shuffle study B gene-mean directions; preserves marginal direction bias; conditional on both changed'))
    pd.concat(recs,ignore_index=True).to_csv(OUT/'exact_site_cross_study.tsv.gz',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'recurrence_summary.tsv',sep='\t',index=False)
    tests=pd.DataFrame(tests); ix=np.argsort(tests.p.to_numpy()); q=np.empty(len(ix))
    q[ix]=np.minimum.accumulate((tests.p.to_numpy()[ix]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1].clip(0,1)
    tests['q']=q; tests.to_csv(OUT/'gene_block_direction_tests.tsv',sep='\t',index=False)
    finish_partial(sites,pd.concat(recs,ignore_index=True))
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)
    print(tests.to_string(index=False),flush=True)


def finish_partial(sites,recurrence):
    rng=np.random.default_rng(20260920)
    keys=['species','route','contrast','category','gene_id','factor']
    sites=sites.copy()
    sites['positive_change']=sites.strict_changed & sites.effect20.gt(0)
    sites['negative_change']=sites.strict_changed & sites.effect20.lt(0)
    sites['coord']=sites.chrom+':'+sites.start0.astype(str)+'-'+sites.end0.astype(str)
    pair=sites.groupby(keys,sort=False).agg(gene_symbol=('gene_symbol','first'),n_sites=('start0','size'),
        n_evaluable=('evaluable','sum'),n_changed=('strict_changed','sum'),n_retained=('retained','sum'),
        n_positive=('positive_change','sum'),n_negative=('negative_change','sum'))
    pair['changed_sites']=sites[sites.strict_changed].groupby(keys).coord.agg(';'.join)
    pair['retained_sites']=sites[sites.retained].groupby(keys).coord.agg(';'.join)
    pair[['changed_sites','retained_sites']]=pair[['changed_sites','retained_sites']].fillna('')
    pair['n_unresolved']=pair.n_sites-pair.n_evaluable
    pair['has_partial']=pair.n_changed.gt(0)&pair.n_retained.gt(0)
    pair['classification']=np.select([
        pair.n_evaluable.lt(2),pair.n_positive.gt(0)&pair.n_negative.gt(0),pair.has_partial,
        pair.n_changed.gt(0)&pair.n_changed.eq(pair.n_evaluable)&pair.n_evaluable.eq(pair.n_sites),
        pair.n_changed.gt(0)],['insufficient_evaluable_sites','mixed','partial',
        'all_annotated_sites_changed','changed_with_other_or_unresolved'],default='no_measured_change')
    # Check the batch implementation against the independently specified rules.
    for x in pair.sample(min(2000,len(pair)),random_state=1).itertuples():
        expected=classify([1]*x.n_positive+[-1]*x.n_negative,x.n_retained,x.n_evaluable,x.n_sites)
        assert expected==(x.classification,x.has_partial)
    pair=pair.reset_index()
    pair.to_csv(OUT/'gene_factor_classification.tsv.gz',sep='\t',index=False)
    ps=[]
    for key,g in pair.groupby(['species','route','contrast','category']):
        eligible=g[g.n_evaluable>=2]; active=eligible[eligible.n_changed>0]
        block=active.groupby('gene_id').agg(success=('has_partial','sum'),total=('has_partial','size'))
        ci=bootstrap_ratio(block,rng)
        ps.append(dict(zip(['species','route','contrast','category'],key))|dict(
            multisite_pairs=len(eligible),changed_multisite_pairs=len(active),changed_genes=active.gene_id.nunique(),
            partial_including_mixed=int(active.has_partial.sum()),partial_genes=active[active.has_partial].gene_id.nunique(),
            partial_fraction=active.has_partial.mean(),ci_low=ci[0],ci_high=ci[1],
            **{c:int(active.classification.eq(c).sum()) for c in ['partial','mixed','all_annotated_sites_changed','changed_with_other_or_unresolved']}))
    pd.DataFrame(ps).to_csv(OUT/'partial_retention_summary.tsv',sep='\t',index=False)
    repeated=[]
    cols=['species','contrast','category','gene_id','factor']
    rr=recurrence[recurrence.recurrent_strict].copy()
    rr['coord']=rr.chrom+':'+rr.start0.astype(str)+'-'+rr.end0.astype(str)
    rec_lookup=rr.groupby(cols).coord.agg(set).to_dict()
    for sp,a,b in [('mouse','CIP','WOT'),('human','GSE118258','GSE150311')]:
        pp=pair[pair.route.eq(a)&pair.has_partial].merge(pair[pair.route.eq(b)&pair.has_partial],on=cols,suffixes=('_a','_b'))
        for x in pp.to_dict('records'):
            changed_intersection=set(x['changed_sites_a'].split(';'))&set(x['changed_sites_b'].split(';'))
            retained_intersection=set(x['retained_sites_a'].split(';'))&set(x['retained_sites_b'].split(';'))
            recurrent_coords=rec_lookup.get(tuple(x[k] for k in cols),set())
            repeated.append(x|dict(shared_recurrent_changed_sites=';'.join(sorted(changed_intersection&recurrent_coords)),
                shared_retained_sites=';'.join(sorted(retained_intersection)),
                exact_partial_recurrence=bool(changed_intersection&recurrent_coords) and bool(retained_intersection)))
    pd.DataFrame(repeated).to_csv(OUT/'cross_study_partial_retention.tsv',sep='\t',index=False)
    sources=[Path(__file__),REG,P/'full-cohort-v2/data/mouse/site_effects.tsv.gz',P/'full-cohort-v2/data/human/site_effects.tsv.gz']
    (OUT/'run_verification.json').write_text(json.dumps(dict(unit_tests=4,skipped=0,
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}),indent=2))
    print(pd.DataFrame(ps).to_string(index=False),flush=True)


if __name__=='__main__':
    import sys
    if '--partial-only' in sys.argv:
        tests=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
        assert tests.wasSuccessful() and not tests.skipped
        finish_partial(pd.read_csv(OUT/'site_factor_inventory.tsv.gz',sep='\t',low_memory=False),
                       pd.read_csv(OUT/'exact_site_cross_study.tsv.gz',sep='\t',low_memory=False))
    else:main()
