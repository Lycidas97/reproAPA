import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score
from common import *


def unit_ranks(a):
    r=rankdata(a,axis=1,method='average').astype(np.float64)
    r-=r.mean(axis=1,keepdims=True)
    norm=np.linalg.norm(r,axis=1,keepdims=True)
    return np.divide(r,norm,out=np.zeros_like(r),where=norm>0)


def run():
    out=RUN/'A04';x=obs()
    wait_file('A04',RUN/'shared/features_ready.json','A03')
    counts=sparse.load_npz(RUN/'shared/ortholog_counts.npz')
    features=np.load(RUN/'shared/feature_sets.npz')['holdout']
    f=pd.read_csv(RUN/'shared/feature_selection.tsv',sep='\t').iloc[features].copy()
    totals=np.load(RUN/'shared/full_gene_umi.npy')
    a=counts[:,features].tocsr();del counts
    detection={study:np.asarray(a[idx].getnnz(axis=0)).ravel()>=10
               for study,idx in x.groupby('true_study',sort=True).indices.items()}
    a=(sparse.diags(10000/totals)@a).tocsr();a.data=np.log1p(a.data)
    # Library means, then equally weighted libraries within each original state.
    group=x.true_study+'|'+x.state_key+'|'+x.sample_key
    code,levels=pd.factorize(group,sort=True)
    sizes=np.bincount(code)
    mapping=sparse.csr_matrix((1/sizes[code],(code,np.arange(len(x)))),shape=(len(levels),len(x)))
    profile=np.asarray((mapping@a).toarray())
    meta=x[['true_study','species_key','state_key','display_group','sample_key']].copy()
    meta['group']=group
    meta=meta.drop_duplicates('group').set_index('group').reindex(levels).reset_index()
    meta['n_cells']=sizes
    meta.to_csv(out/'state_library_profile_registry.tsv',sep='\t',index=False)
    np.save(out/'heldout_state_library_logCPM_means.npy',profile)
    state_meta=meta.drop_duplicates('state_key').set_index('state_key')[['true_study','species_key','display_group']]
    state_mean=pd.DataFrame(profile,index=meta.state_key).groupby(level=0).mean().reindex(state_meta.index).to_numpy()
    rng=np.random.default_rng(SEED)
    query=np.sort(np.concatenate([rng.choice(ix,min(100,len(ix)),replace=False)
                   for ix in x.groupby(['state_key','sample_key'],sort=True).indices.values()]))
    x.iloc[query][['cell_id','state_key','sample_key']].assign(atlas_index=query).to_csv(out/'query_registry.tsv',sep='\t',index=False)
    results=[];aucs=[];folds=[]
    for test_study,idx in x.groupby('true_study',sort=True).indices.items():
        species=x.iloc[idx[0]].species_key
        for mode in ['within_species_leave_study_out','cross_species']:
            train=state_meta.species_key.eq(species)&state_meta.true_study.ne(test_study) if mode.startswith('within') else state_meta.species_key.ne(species)
            train_meta=state_meta[train]
            training_studies=sorted(train_meta.true_study.unique())
            eligible=np.stack([detection[s] for s in training_studies]).sum(axis=0)>=2
            assert eligible.sum()>=100, (test_study,mode,int(eligible.sum()))
            folds.append(dict(test_study=test_study,mode=mode,training_studies=training_studies,genes=int(eligible.sum())))
            pd.DataFrame(dict(human_gene_id=f.human_gene_id.to_numpy(),training_detection_eligible=eligible)).to_csv(
                out/f'{test_study}_{mode}_features.tsv',sep='\t',index=False)
            # Score each training state separately, then average state scores within a study-family,
            # and equally across training studies. This avoids giving a study extra votes for more clusters.
            reference=unit_ranks(state_mean[train.to_numpy()][:,eligible])
            state_names=train_meta.index.to_numpy()
            group_keys=train_meta.true_study+'|'+train_meta.display_group
            sg,sg_levels=pd.factorize(group_keys,sort=True)
            sg_count=np.bincount(sg)
            s_to_group=np.zeros((len(train_meta),len(sg_levels)))
            s_to_group[np.arange(len(train_meta)),sg]=1/sg_count[sg]
            sg_meta=train_meta.assign(sg=group_keys).drop_duplicates('sg').set_index('sg').reindex(sg_levels)
            family_names=np.sort(sg_meta.display_group.unique())
            group_to_family=np.stack([(sg_meta.display_group.eq(v)/sg_meta.display_group.eq(v).sum()).to_numpy() for v in family_names],axis=1)
            tq=query[x.iloc[query].true_study.eq(test_study).to_numpy()]
            fold_scores=[];fold_rows=[]
            progress('A04','RUNNING',stage='heldout_expression_prediction',test_study=test_study,mode=mode,genes=int(eligible.sum()),queries=len(tq))
            for start in range(0,len(tq),512):
                qi=tq[start:start+512]
                values=a[qi][:,eligible].toarray()
                ranks=unit_ranks(values)
                state_scores=ranks@reference.T
                scores=state_scores@s_to_group@group_to_family
                order=np.argsort(-scores,axis=1)
                r=x.iloc[qi][['cell_id','species_key','true_study','state_key','sample_key','display_group']].copy()
                r['mode']=mode;r['heldout_features']=int(eligible.sum())
                r['best_reference_state']=state_names[state_scores.argmax(axis=1)]
                r['best_reference_family']=family_names[order[:,0]]
                r['best_family_spearman']=scores[np.arange(len(qi)),order[:,0]]
                r['family_score_margin']=scores[np.arange(len(qi)),order[:,0]]-scores[np.arange(len(qi)),order[:,1]]
                r['reference_contains_expected_family']=r.display_group.isin(family_names)
                r['query_evaluable']=np.ptp(values,axis=1)>0
                r['descriptive_label_agreement']=r.best_reference_family.eq(r.display_group)
                r['interpretation']=np.where(r.reference_contains_expected_family,'candidate_expression_support','no_expected_family_in_reference')
                r.loc[~r.query_evaluable,'best_reference_family']='NOT_EVALUABLE'
                r.loc[~r.query_evaluable,'descriptive_label_agreement']=False
                r.loc[~r.query_evaluable,'interpretation']='constant_expression_in_heldout_features'
                # A direct SciPy calculation validates the actual rank-correlation operation.
                if start==0:
                    ri=np.flatnonzero(np.ptp(state_mean[train.to_numpy()][:,eligible],axis=1)>0)
                    vi=np.flatnonzero(r.query_evaluable.to_numpy())
                    if len(ri) and len(vi):
                        direct=spearmanr(values[vi[0]],state_mean[train.to_numpy()][ri[0],eligible]).statistic
                        assert np.isclose(direct,state_scores[vi[0],ri[0]],atol=1e-10)
                fold_scores.append(scores);fold_rows.append(r)
            joined=pd.concat(fold_rows,ignore_index=True);scores=np.vstack(fold_scores)
            for j,family in enumerate(family_names):
                valid=joined.query_evaluable.to_numpy()
                y=joined.loc[valid,'display_group'].eq(family).to_numpy()
                auc=float(roc_auc_score(y,scores[valid,j])) if y.any() and (~y).any() else None
                aucs.append(dict(test_study=test_study,mode=mode,family=family,queries=len(y),positive_queries=int(y.sum()),AUROC=auc))
            results.append(joined)
    result=pd.concat(results,ignore_index=True)
    result.to_csv(out/'heldout_cell_predictions.tsv.gz',sep='\t',index=False)
    pd.DataFrame(aucs).to_csv(out/'heldout_family_AUROC.tsv',sep='\t',index=False)
    result.groupby(['species_key','true_study','state_key','display_group','sample_key','mode']).agg(
        queries=('cell_id','size'),label_agreement=('descriptive_label_agreement','mean'),
        mean_best_family_correlation=('best_family_spearman','mean'),mean_margin=('family_score_margin','mean')).reset_index().to_csv(
            out/'state_library_expression_support.tsv',sep='\t',index=False)
    write_json(out/'folds.json',folds)
    finish('A04','feature_holdout_LOSO',queries=len(query),prediction_rows=len(result),heldout_candidates=len(features),
        folds=len(folds),method='Spearman expression similarity to library-balanced state profiles, study-balanced family scores',
        next_steps=['Original-marker conflicts and reference-label evidence','Doublet exclusion profile refits','Manual evidence-tier review'],
        boundary='Feature holdout in the same cells; label agreement and AUROC do not alone establish homology or independent experiments.')


if __name__=='__main__': execute('A04',run)
