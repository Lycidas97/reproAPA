"""Re-evaluate strict Figure 2 APA coherence across expression kNN sizes."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

import test_figure2_apa_manifold_coherence as pilot


OUT = pilot.PROJECT / 'outputs/figure2-apa-manifold-coherence-robustness-v1'


def graph_for_study(latent, k, groups=None):
    n=len(latent)
    groups=np.zeros(n,dtype=int) if groups is None else np.asarray(groups)
    rows=[];cols=[]
    for group in pd.unique(groups):
        chosen=np.flatnonzero(groups==group)
        neighbors=min(k,len(chosen)-1)
        if neighbors<1:
            continue
        found=cKDTree(latent[chosen]).query(latent[chosen],k=neighbors+1,workers=-1)[1][:,1:]
        rows.append(np.repeat(chosen,neighbors))
        cols.append(chosen[found.ravel()])
    graph=sparse.csr_matrix((np.ones(sum(len(x) for x in rows),np.float32),
                             (np.concatenate(rows),np.concatenate(cols))),shape=(n,n))
    graph=graph.maximum(graph.T)
    graph.setdiag(0)
    graph.eliminate_zeros()
    assert (graph-graph.T).nnz==0
    return graph


def main(permutations):
    OUT.mkdir(parents=True,exist_ok=True)
    nodes=pd.read_csv(pilot.INTEGRATION/'common_metacell_nodes.tsv.gz',sep='\t',
                      usecols=['metacell_id','study_id','species_key','library_key','study_group'],low_memory=False)
    apa=pd.read_csv(pilot.APA,sep='\t')
    labels=pilot.state_assignment()
    frame=nodes.merge(apa,on='metacell_id',validate='one_to_one').merge(
        labels,on='metacell_id',validate='one_to_one').merge(
        pilot.pas_support(),on='metacell_id',validate='one_to_one')
    assert len(frame)==len(nodes)==11104 and frame.metacell_id.eq(nodes.metacell_id).all()
    latent=np.load(pilot.INTEGRATION/'species_balanced_scvi_latent.npy')
    assert latent.shape==(len(frame),30) and np.isfinite(latent).all()
    frozen=sparse.load_npz(pilot.INTEGRATION/'local_expression_kernel.npz').tocsr()
    frozen=frozen.maximum(frozen.T)
    frozen.data[:]=1
    frozen.setdiag(0)
    frozen.eliminate_zeros()
    records=[]
    for si,study in enumerate(pd.unique(frame.study_id)):
        indices=np.flatnonzero(frame.study_id.eq(study))
        within=frame.iloc[indices].reset_index(drop=True)
        print('STUDY',study,flush=True)
        for k in (10,15,20,30):
            full_graph=graph_for_study(latent[indices],k,within.study_group)
            if k==20:
                assert (full_graph-frozen[indices][:,indices]).nnz==0,study
            for ci,klass in enumerate(pilot.CLASSES):
                eligible=within.state_purity.ge(.8)&within[f'{klass}_centered'].notna()
                chosen=np.flatnonzero(eligible)
                local=full_graph[chosen][:,chosen].tocsr()
                connected=np.diff(local.indptr)>0
                chosen=chosen[connected]
                local=local[connected][:,connected].tocsr()
                subset=within.iloc[chosen]
                assert np.diff(local.indptr).min()>0
                value=subset[f'{klass}_centered'].to_numpy(float)
                design=pd.get_dummies(subset.state_group,dtype=float).to_numpy()
                design=np.column_stack([design,np.log1p(subset[f'{klass}_n'].to_numpy(float)),
                                        np.log1p(subset.primary_pas_umi.to_numpy(float))])
                residual=value-design@np.linalg.lstsq(design,value,rcond=None)[0]
                blocks=(subset.library_key_x+'|'+subset.state_group).to_numpy()
                result=pilot.permuted_moran(local,residual,blocks,permutations,
                                             seed=270924+si*1000+ci*100+k)
                records.append(dict(study_id=study,species=subset.species_key.iloc[0],
                                    cluster=klass,k=k,n_metacells=len(subset),
                                    n_edges=local.nnz//2,n_isolated=int((~connected).sum()),
                                    permutations=permutations,**result))
    output=pd.DataFrame(records)
    for k in (10,15,20,30):
        take=output.k.eq(k)
        output.loc[take,'q_bh']=pilot.bh(output.loc[take,'empirical_p'].to_numpy())
    output.to_csv(OUT/'knn_sensitivity.tsv',sep='\t',index=False,float_format='%.9g')
    print(output.groupby('k').agg(significant=('q_bh',lambda x:(x<.05).sum()),
                                  median_excess=('excess_i','median'),
                                  min_excess=('excess_i','min')).to_string(),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--permutations',type=int,default=9999)
    args=parser.parse_args()
    main(args.permutations)
