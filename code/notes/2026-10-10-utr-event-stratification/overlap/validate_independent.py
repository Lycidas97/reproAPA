"""Independent arithmetic and source-roster audit, without production imports."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'notes/2026-10-10-ortholog-threshold-scan'
CLASSES = HERE.parent / 'classification/event_class_pairs.tsv'
check_count = 0


def check(condition, meaning):
    global check_count
    check_count += 1
    if not condition:
        raise AssertionError(meaning)


def bh(p):
    x = np.asarray(p, float)
    order = np.argsort(x, kind='stable')
    ordered = x[order]
    corrected = np.minimum.accumulate((ordered * len(x) / np.arange(1, len(x)+1))[::-1])[::-1]
    out = np.empty(len(x))
    out[order] = np.minimum(1, corrected)
    return out


def near(a, b, message, atol=1e-12):
    check(np.allclose(a, b, atol=atol, rtol=0, equal_nan=True), message)


def main():
    cls = pd.read_csv(CLASSES, sep='\t')
    members = pd.read_csv(SOURCE / 'measured_ortholog_members.tsv.gz', sep='\t', float_precision='round_trip')
    members = members[members.contrast.eq('P_M')]
    features = pd.read_csv(SOURCE / 'ortholog_assay_features.tsv.gz', sep='\t', float_precision='round_trip')
    features = features[features.contrast.eq('P_M')]
    r = pd.read_csv(HERE / 'class_direction_3x3.tsv', sep='\t', float_precision='round_trip')
    diag = pd.read_csv(HERE / 'stratum_diagnostics.tsv.gz', sep='\t', float_precision='round_trip')
    mobility = pd.read_csv(HERE / 'mobility_summary.tsv', sep='\t', float_precision='round_trip')
    support = pd.read_csv(HERE / 'null_support.tsv', sep='\t', float_precision='round_trip')
    thresholds = [.01,.015,.02,.025,.03,.035,.04,.045,.05]
    labels = ['proximal','distal','nonrecurrent']
    contexts = ['pair_class','threshold','model']
    check(len(r)==324, 'Every planned 3x3 cell is retained in every context.')
    check(len(r.groupby(contexts))==36, 'No post-result omission of a class, threshold, or model.')
    check(not cls.duplicated(['human_gene_id','mouse_gene_id']).any(), 'One ortholog unit per pair.')
    joined = members.merge(cls, on=['human_gene_id','mouse_gene_id'], validate='many_to_one')
    independent_random = []
    max_expected_difference = 0.
    with np.load(HERE / 'null_draws.npz') as archive:
        for (event_class,threshold,model), part in r.groupby(contexts, sort=True):
            pair = joined[joined.pair_class.eq(event_class) & joined.threshold.eq(threshold)].reset_index(drop=True)
            check(pair.human_gene_id.is_unique and pair.mouse_gene_id.is_unique,
                  'Class restriction preserves independent ortholog identities.')
            check(part.measured_orthologs.eq(len(pair)).all(), 'Foreground and background use the same class.')
            cols = ['measured_pattern','proximal_capacity','distal_capacity','PAS_number_bin']
            if model=='full_assay':
                cols += ['mean_GEX_bin','APA_coverage_bin','baseline_bin']
            probabilities, species_labels, groups = {}, {}, {}
            for species in ['human','mouse']:
                ff = features[features.species.eq(species)&features.threshold.eq(threshold)]
                ff = ff.set_index('gene_id').loc[pair[species+'_gene_id'].tolist()]
                check(set(ff.index)==set(pair[species+'_gene_id']), 'No other-class member enters a stratum.')
                ll = pair[species+'_direction'].map(dict(zip(labels, range(3)))).to_numpy()
                groups[species] = list(ff.groupby(cols, sort=True).indices.values())
                probabilities[species] = np.zeros((len(pair),3))
                species_labels[species] = ll
                singleton_n = frozen_n = 0
                expected_change = 0.
                for idx in groups[species]:
                    counts = np.array([int((ll[idx]==a).sum()) for a in range(3)])
                    probabilities[species][idx] = counts/len(idx)
                    singleton_n += len(idx) if len(idx)==1 else 0
                    frozen_n += len(idx) if (counts>0).sum()==1 else 0
                    expected_change += len(idx)*(1-np.square(counts/len(idx)).sum())
                dd = diag[diag.pair_class.eq(event_class)&diag.threshold.eq(threshold)&
                          diag.model.eq(model)&diag.species.eq(species)]
                mm = mobility[mobility.pair_class.eq(event_class)&mobility.threshold.eq(threshold)&
                              mobility.model.eq(model)&mobility.species.eq(species)].iloc[0]
                check(len(dd)==len(groups[species]), 'Reported strata match old bins after class restriction.')
                check(int(dd.n.sum())==len(pair), 'Diagnostics cover the entire restricted universe.')
                near(mm.singleton_gene_fraction,singleton_n/len(pair),'Singleton fraction is gene-weighted.')
                near(mm.frozen_label_fraction,frozen_n/len(pair),'Homogeneous fraction is gene-weighted.')
                near(mm.expected_changed_fraction,expected_change/len(pair),'Mobility is not counted by strata.')
            expected = probabilities['human'].T @ probabilities['mouse']
            key = '|'.join(map(str,(event_class,threshold,model)))
            null = archive[key]
            check(null.shape==(100000,3,3), 'Exactly 100000 null draws in every planned context.')
            check(np.all(null.sum((1,2))==len(pair)), 'Null never escapes the class-qualified background.')
            for a in range(3):
                check(np.all(null.sum(2)[:,a]==(species_labels['human']==a).sum()),'Human margins remain fixed.')
                check(np.all(null.sum(1)[:,a]==(species_labels['mouse']==a).sum()),'Mouse margins remain fixed.')
            for _, row in part.iterrows():
                a,b = labels.index(row.human_direction),labels.index(row.mouse_direction)
                observed = int(((species_labels['human']==a)&(species_labels['mouse']==b)).sum())
                values = null[:,a,b]
                check(row.observed==observed, 'Overlap count independently recovered from fixed source labels.')
                difference = abs(row.analytic_expected-expected[a,b])
                max_expected_difference=max(max_expected_difference,difference)
                near(row.analytic_expected,expected[a,b],'Expected count uses restricted strata, not full 751.')
                near(row.null_mean,values.mean(),'Archived draws reproduce null mean.')
                near(row.null_sd,values.std(ddof=1),'Archived draws reproduce null variation.')
                near(row.excess,observed-expected[a,b],'Excess uses the correct class-conditional expectation.')
                check(row.tail_count==(values>=observed).sum(),'Upper tail includes equality.')
                check(row.null_support_points==len(np.unique(values)),'Null support is fully reported.')
                if a<2 and b<2:
                    near(row.p_enrichment,(1+(values>=observed).sum())/100001,'Permutation add-one correction.')
                    ss = support[support.pair_class.eq(event_class)&support.threshold.eq(threshold)&
                        support.model.eq(model)&support.human_direction.eq(labels[a])&
                        support.mouse_direction.eq(labels[b])]
                    v,n=np.unique(values,return_counts=True)
                    check(np.array_equal(ss.overlap_count.to_numpy(),v),'Support points are exact.')
                    check(np.array_equal(ss.draws.to_numpy(),n),'Support counts are exact.')
                    near(ss.probability.to_numpy(),n/100000,'Support mass uses all draws.')
                    if values.std()==0:
                        check(not row.inferentially_mobile,'Frozen null is explicitly flagged.')
                        near(row.p_enrichment,1.,'Frozen comparison retained with p=1, not deleted from family.')
                else:
                    check(pd.isna(row.p_enrichment),'Non-directional cells are not inferential tests.')
            # Independent small replay checks the null, not only output arithmetic.
            if threshold in [.02,.05]:
                draws=5000
                rng=np.random.default_rng(90217+len(independent_random))
                simulated={}
                for species in ['human','mouse']:
                    ll=species_labels[species]
                    arr=np.tile(ll,(draws,1))
                    for idx in groups[species]:
                        if np.unique(ll[idx]).size>1:
                            # Independent ranking of IID continuous draws generates uniform permutations.
                            order=np.argsort(rng.random((draws,len(idx))),axis=1)
                            arr[:,idx]=ll[idx][order]
                    simulated[species]=arr
                values=((simulated['human']==0)&(simulated['mouse']==0)).sum(1)
                old=null[:,0,0]
                margin=6*np.sqrt(old.var()/len(old)+values.var()/draws)+1e-12
                check(abs(old.mean()-values.mean())<=margin,'Independent permutation construction agrees in mean.')
                observed=int(((species_labels['human']==0)&(species_labels['mouse']==0)).sum())
                pnew=float((values>=observed).mean());pold=float((old>=observed).mean())
                margin=6*np.sqrt(pnew*(1-pnew)/draws+pold*(1-pold)/100000)+1/draws+1e-12
                check(abs(pnew-pold)<=margin,'Independent permutation construction agrees in tail.')
                independent_random.append(dict(pair_class=event_class,threshold=threshold,model=model,
                    draws=draws,seed=90217+len(independent_random),mean=float(values.mean()),tail=pnew))
    primary=r.human_direction.eq('proximal')&r.mouse_direction.eq('proximal')
    for threshold,column in [(.02,'q_primary4'),(.05,'q_sensitivity4')]:
        mask=primary&r.threshold.eq(threshold)
        check(mask.sum()==4,'Primary family retains both classes and both models.')
        near(r.loc[mask,column],bh(r.loc[mask,'p_enrichment']),'Primary correction independently recovered.')
    check(primary.sum()==36,'Threshold scan family includes every planned threshold.')
    near(r.loc[primary,'q_scan36'],bh(r.loc[primary,'p_enrichment']),'No favorable-threshold correction.')
    directional=r.human_direction.isin(labels[:2])&r.mouse_direction.isin(labels[:2])
    check(directional.sum()==144,'Secondary scan includes all four directional cells.')
    for threshold in thresholds:
        mask=directional&r.threshold.eq(threshold)
        check(mask.sum()==16,'Secondary family includes near-near too.')
        near(r.loc[mask,'q_secondary16_per_threshold'],bh(r.loc[mask,'p_enrichment']),'Secondary threshold BH.')
    near(r.loc[directional,'q_secondary_scan144'],bh(r.loc[directional,'p_enrichment']),'Secondary entire scan BH.')
    strict=pd.read_csv(HERE/'full_axis_tandem_members.tsv.gz',sep='\t',float_precision='round_trip')
    expected_ids=cls.loc[cls.pair_class.eq('both_tandem')&cls.both_full_axis_tandem,'human_gene_id']
    for threshold in thresholds:
        ss=strict[strict.threshold.eq(threshold)]
        check(set(ss.human_gene_id)==set(expected_ids),'Strict sensitivity is tandem AND both full axes.')
    manifest=json.loads((HERE/'analysis_manifest.json').read_text())
    for name,value in manifest['input_sha256'].items():
        check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==value,'Frozen inputs remained unchanged.')
    for name,value in manifest['output_sha256'].items():
        check(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==value,'Production output hash is confirmed.')
    result=dict(status='PASS',checks=check_count,skipped=0,
        maximum_expected_difference=max_expected_difference,
        independent_permutation_replays=independent_random,
        no_production_functions_imported=True,
        source_sha256={str(CLASSES.relative_to(ROOT)):hashlib.sha256(CLASSES.read_bytes()).hexdigest()})
    (HERE/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
