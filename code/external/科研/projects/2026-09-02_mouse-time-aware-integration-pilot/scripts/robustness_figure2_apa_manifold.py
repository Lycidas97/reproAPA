"""Audit fixed-gene panels and graph sensitivity for Figure 2 APA coherence."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

import project_k3_group_wpp_umap as projection
import test_figure2_apa_manifold_coherence as pilot
from robustness_figure2_apa_knn import graph_for_study


OUT = pilot.PROJECT / 'outputs/figure2-apa-manifold-coherence-robustness-v1'


def gene_wpp(source, support, members):
    meta = pd.read_csv(source / 'metacells.tsv.gz', sep='\t', usecols=['metacell_id'])
    features = pd.read_csv(source / 'features.tsv.gz', sep='\t',
                           usecols=projection.FEATURE_COLUMNS, low_memory=False)
    assert np.array_equal(features.feature_index_1based.to_numpy(), np.arange(1, len(features)+1))
    features = features.loc[features.keep_for_apa_primary &
                            features.primary_gene_id.isin(members.gene_id)].copy()
    features = features.sort_values(['primary_gene_id','cut_site','atlas_pas_id']).reset_index(drop=True)
    reference = pd.read_csv(support / source.name / 'primary_pas_reference_matches.tsv.gz',
                            sep='\t',usecols=['atlas_pas_id','matched_source_count_50nt'])
    features = features.merge(reference,on='atlas_pas_id',validate='one_to_one')
    genes = pd.Index(features.primary_gene_id.unique())
    gene_index = genes.get_indexer(features.primary_gene_id)
    starts = np.r_[0,np.flatnonzero(np.diff(gene_index))+1]
    signed=features.cut_site.to_numpy(np.float64)*np.where(features.strand.eq('+'),1.,-1.)
    low=np.minimum.reduceat(signed,starts)
    high=np.maximum.reduceat(signed,starts)
    matched=features.matched_source_count_50nt.ge(1).to_numpy()
    valid=(features.groupby('primary_gene_id',sort=False).chrom.nunique().to_numpy()==1)
    valid &= features.groupby('primary_gene_id',sort=False).strand.nunique().to_numpy()==1
    valid &= np.isin(features.strand.iloc[starts],['+','-']) & (high>low)
    valid &= np.add.reduceat(matched.astype(np.int32),starts)>=2
    relative=np.divide(signed-low[gene_index],(high-low)[gene_index],
                       out=np.zeros(len(features)),where=(high-low)[gene_index]>0)
    counts=sparse.load_npz(source/'metacell_by_pas_counts.npz').tocsr()
    assert len(meta)==counts.shape[0]
    counts=counts[:,features.feature_index_1based.to_numpy(np.int64)-1].multiply(matched).tocsr()
    mapping=sparse.csr_matrix((np.ones(len(features),np.float32),
                               (np.arange(len(features)),gene_index)),shape=(len(features),len(genes)))
    weighting=sparse.csr_matrix((relative.astype(np.float32),
                                 (np.arange(len(features)),gene_index)),shape=mapping.shape)
    totals=(counts@mapping).toarray()
    numerators=(counts@weighting).toarray()
    observable=(totals>=20)&valid
    wpp=np.divide(numerators,totals,out=np.full(totals.shape,np.nan,np.float32),where=observable)
    class_ids=members.set_index('gene_id').cluster.reindex(genes).to_numpy()
    original=pd.read_csv(projection.OUT/'k3_group_wpp_metacells.tsv.gz',sep='\t')
    original=original.set_index('metacell_id').reindex(meta.metacell_id)
    observed=np.isfinite(wpp)
    gene_mean=np.divide(np.nansum(wpp,axis=0),observed.sum(axis=0),
                        out=np.full(wpp.shape[1],np.nan),where=observed.sum(axis=0)>0)
    for klass in projection.CLASSES:
        take=class_ids==klass
        count=observed[:,take].sum(axis=1)
        actual=np.divide(np.nansum(wpp[:,take]-gene_mean[take],axis=1),count,
                         out=np.full(len(meta),np.nan),where=count>0)
        expected=original[f'{klass}_centered'].to_numpy()
        assert np.allclose(actual,expected,atol=2e-5,equal_nan=True), (source.name,klass)
    return meta.metacell_id.to_numpy(),genes.to_numpy(),class_ids,wpp


def main(permutations):
    OUT.mkdir(parents=True,exist_ok=True)
    members=pd.read_csv(projection.MEMBERSHIP,sep='\t',usecols=['species','gene_id','cluster'])
    nodes=pd.read_csv(pilot.INTEGRATION/'common_metacell_nodes.tsv.gz',sep='\t',
                      usecols=['metacell_id','study_id','species_key','library_key','study_group'],low_memory=False)
    metadata=nodes.merge(pilot.state_assignment(),on='metacell_id',validate='one_to_one').merge(
        pilot.pas_support(),on='metacell_id',validate='one_to_one').set_index('metacell_id')
    assert len(metadata)==11104
    original=pd.read_csv(projection.OUT/'k3_group_wpp_metacells.tsv.gz',sep='\t').set_index('metacell_id')
    latent=np.load(pilot.INTEGRATION/'species_balanced_scvi_latent.npy')
    assert latent.shape==(len(nodes),30)
    original_index=pd.Index(nodes.metacell_id)
    audits=[];tests=[]
    for species,root,support in [('human',projection.HUMAN,projection.HUMAN_SUPPORT),
                                 ('mouse',projection.MOUSE,projection.MOUSE_SUPPORT)]:
        specific=members[members.species.eq(species)]
        for source in sorted((root/'studies').iterdir()):
            ids,genes,classes,wpp=gene_wpp(source,support,specific)
            state=metadata.loc[ids]
            assert state.species_key.eq(species).all()
            latent_index=original_index.get_indexer(ids)
            assert (latent_index>=0).all()
            gene_means=np.divide(np.nansum(wpp,axis=0),np.isfinite(wpp).sum(axis=0),
                                 out=np.full(wpp.shape[1],np.nan),where=np.isfinite(wpp).sum(axis=0)>0)
            print('STUDY',source.name,'metacells',len(ids),flush=True)
            for klass in projection.CLASSES:
                group_values=wpp[:,classes==klass]
                coverage=np.isfinite(group_values).mean(axis=0)
                for threshold in (.7,.8,.9):
                    selected=coverage>=threshold
                    complete=np.isfinite(group_values[:,selected]).all(axis=1) if selected.any() else np.zeros(len(ids),bool)
                    audits.append(dict(species=species,study=source.name,cluster=klass,
                                       threshold=threshold,candidate_genes=int((classes==klass).sum()),
                                       selected_genes=int(selected.sum()),
                                       complete_metacells=int(complete.sum()),
                                       complete_fraction=float(complete.mean()),
                                       selected_gene_ids=';'.join(genes[classes==klass][selected])))
                    eligible=complete & state.state_purity.ge(.8).to_numpy()
                    subset_ids=ids[eligible]
                    if len(subset_ids)<30:
                        for metric in ('original_same_nodes','fixed_panel'):
                            tests.append(dict(species=species,study=source.name,cluster=klass,
                                              threshold=threshold,metric=metric,status='INSUFFICIENT_N',
                                              selected_genes=int(selected.sum()),n_metacells=len(subset_ids)))
                        continue
                    subset_meta=state.loc[subset_ids]
                    graph=graph_for_study(latent[latent_index[eligible]],20,
                                          subset_meta.study_group)
                    fixed_genes=np.flatnonzero(classes==klass)[selected]
                    fixed_score=np.mean(wpp[eligible][:,fixed_genes]-gene_means[fixed_genes],axis=1)
                    original_score=original.loc[subset_ids,f'{klass}_centered'].to_numpy(float)
                    assert np.isfinite(fixed_score).all() and np.isfinite(original_score).all()
                    original_n=original.loc[subset_ids,f'{klass}_n'].to_numpy(float)
                    design=pd.get_dummies(subset_meta.state_group,dtype=float).to_numpy()
                    design=np.column_stack([design,np.log1p(original_n),
                                            np.log1p(subset_meta.primary_pas_umi.to_numpy(float))])
                    blocks=(subset_meta.library_key_x+'|'+subset_meta.state_group).to_numpy()
                    for mi,(metric,score) in enumerate((('original_same_nodes',original_score),
                                                        ('fixed_panel',fixed_score))):
                        residual=score-design@np.linalg.lstsq(design,score,rcond=None)[0]
                        stat=pilot.permuted_moran(graph,residual,blocks,permutations,
                              seed=280924+int(threshold*100)*1000+len(tests)*3+mi)
                        tests.append(dict(species=species,study=source.name,cluster=klass,
                                          threshold=threshold,metric=metric,status='PASS',
                                          selected_genes=len(fixed_genes),n_metacells=len(subset_ids),
                                          n_undirected_edges=graph.nnz//2,
                                          n_blocks=len(pd.unique(blocks)),permutations=permutations,**stat))
    audit=pd.DataFrame(audits)
    audit.to_csv(OUT/'fixed_panel_feasibility.tsv',sep='\t',index=False)
    results=pd.DataFrame(tests)
    for (threshold,metric),index in results.loc[results.status.eq('PASS')].groupby(
            ['threshold','metric']).groups.items():
        results.loc[index,'q_bh']=pilot.bh(results.loc[index,'empirical_p'].to_numpy())
    results.to_csv(OUT/'fixed_panel_coherence.tsv',sep='\t',index=False,float_format='%.9g')
    print(audit.drop(columns='selected_gene_ids').to_string(index=False),flush=True)
    print(results.groupby(['threshold','metric']).agg(n=('status','size'),
          estimable=('status',lambda x:(x=='PASS').sum()),
          significant=('q_bh',lambda x:(x<.05).sum()),
          median_excess=('excess_i','median')).to_string(),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--permutations',type=int,default=9999)
    args=parser.parse_args()
    main(args.permutations)
