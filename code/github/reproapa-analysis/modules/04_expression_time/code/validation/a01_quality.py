import faiss
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from common import *


def evaluate(z,x,queries,name):
    faiss.omp_set_num_threads(2)
    index=faiss.IndexFlatL2(z.shape[1]);index.add(np.ascontiguousarray(z,dtype=np.float32))
    pieces=[]
    for start in range(0,len(queries),256):
        q=queries[start:start+256]
        distance,ids=index.search(np.ascontiguousarray(z[q],dtype=np.float32),51)
        ids=np.stack([row[row!=cell][:50] for row,cell in zip(ids,q)])
        for k in [15,30,50]:
            nn=ids[:,:k]
            r=x.iloc[q][['cell_id','species_key','true_study','sample_key','state_key','display_group']].copy()
            r['representation']=name;r['k']=k
            for label in ['species_key','true_study','sample_key','state_key','display_group']:
                v=x[label].to_numpy()
                r[label+'_same_fraction']=(v[nn]==v[q,None]).mean(axis=1)
            shared=r.display_group.str.startswith('shared_').to_numpy()
            same_family=x.display_group.to_numpy()[nn]==r.display_group.to_numpy()[:,None]
            opposite_species=x.species_key.to_numpy()[nn]!=r.species_key.to_numpy()[:,None]
            r['shared_family_query']=shared
            r['cross_species_same_family_fraction']=(same_family&opposite_species).mean(axis=1)
            pieces.append(r)
        if start%2560==0:progress('A01','RUNNING',stage='exact_full_reference_neighbors',representation=name,queries_done=start,queries_total=len(queries))
    result=pd.concat(pieces,ignore_index=True)
    result.to_csv(RUN/'A01'/f'{name}_query_metrics.tsv.gz',sep='\t',index=False)
    values=[c for c in result if c.endswith('_fraction')]
    result.groupby(['representation','k','species_key','true_study','state_key','display_group'])[values].mean().reset_index().to_csv(
        RUN/'A01'/f'{name}_state_metrics.tsv',sep='\t',index=False)
    np.save(RUN/'A01'/f'{name}_latent.npy',z)


def run():
    x=obs();rng=np.random.default_rng(SEED)
    queries=np.sort(np.concatenate([rng.choice(ix,min(200,len(ix)),replace=False)
                  for ix in x.groupby(['species_key','true_study','state_key'],sort=True).indices.values()]))
    q=x.iloc[queries][['cell_id','state_key']].copy(); q['atlas_index']=queries
    q['inclusion_probability']=q.state_key.map(q.state_key.value_counts())/q.state_key.map(x.state_key.value_counts())
    q.to_csv(RUN/'A01/query_registry.tsv',sep='\t',index=False)
    if not (RUN/'A01/existing_representations.complete.json').exists():
        z=np.load(BASE/'species_balanced_raw_cell_scvi_latent.npy')
        evaluate(z,x,queries,'joint_scvi_30D')
        evaluate(x[['UMAP1','UMAP2']].to_numpy(np.float32),x,queries,'joint_UMAP_2D_display_only')
        write_json(RUN/'A01/existing_representations.complete.json',dict(queries=len(queries),reference_cells=len(x)))
    wait_file('A01',RUN/'shared/counts_ready.json','A03')
    counts=sparse.load_npz(RUN/'shared/ortholog_counts.npz')
    f=pd.read_csv(RUN/'shared/ortholog_features.tsv',sep='\t')
    original=pd.read_csv(ROOT/'inputs/S6-integration-orthologs.tsv',sep='\t')
    ix=pd.Index(f.human_gene_id).get_indexer(original.human_gene_id)
    assert (ix>=0).all() and len(ix)==2128
    total=np.load(RUN/'shared/full_gene_umi.npy')
    data=counts[:,ix].astype(np.float32);del counts
    data=sparse.diags(10000/np.maximum(total,1))@data;data=data.tocsr();data.data=np.log1p(data.data)
    train=np.load(RUN/'shared/pilot_indices.npy')
    a=data[train].toarray();mean=a.mean(0,dtype=np.float64);std=a.std(0,dtype=np.float64);std=np.maximum(std,1e-6)
    a=np.clip((a-mean)/std,-10,10).astype(np.float32)
    pca=PCA(n_components=30,svd_solver='randomized',random_state=SEED)
    pca.fit(a);del a
    z=np.empty((len(x),30),np.float32)
    for start in range(0,len(x),2048):
        a=np.clip((data[start:start+2048].toarray()-mean)/std,-10,10).astype(np.float32)
        z[start:start+len(a)]=pca.transform(a)
    np.savez(RUN/'A01/uncorrected_PCA_parameters.npz',features=ix,training_indices=train,
             mean=mean,scale=std,components=pca.components_,explained_variance_ratio=pca.explained_variance_ratio_)
    evaluate(z,x,queries,'uncorrected_ortholog_PCA30')
    finish('A01','neighbor_baselines',queries=len(queries),reference_cells=len(x),
        next_steps=['A04 independently supported labels','biological preservation and formal batch metrics'],
        interpretation='Exact full-reference neighborhoods; existing-label preservation is descriptive, not annotation truth.')


if __name__=='__main__': execute('A01',run)
