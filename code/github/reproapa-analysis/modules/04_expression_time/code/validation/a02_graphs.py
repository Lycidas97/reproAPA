import ast
import gc
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.linalg import orthogonal_procrustes
from scipy.sparse.csgraph import connected_components
from common import *


def normalize(g):
    g=g.tocsr().astype(np.float32)
    total=np.asarray(g.sum(1)).ravel()
    return (sparse.diags(np.divide(1,total,out=np.zeros_like(total),where=total>0))@g).tocsr()


def sym_scale(g):
    g=(g+g.T)*.5;g=g.tocsr();g/=g.data.max();g.eliminate_zeros();return g


def run():
    from umap.umap_ import simplicial_set_embedding, find_ab_params
    x=obs();out=RUN/'A02';n_h=int(x.species_key.eq('human').sum())
    assert x.iloc[:n_h].species_key.eq('human').all()
    hroot=HUMAN/'results_temporal_v1'
    mroot=MOUSE/'results_full_single_cell_time_umap_v1'
    H=sparse.load_npz(hroot/'layout/time_graph.npz').tocsr()
    M=sparse.load_npz(mroot/'layout_graph.npz').tocsr()
    hc=pd.read_csv(hroot/'expression/node_order.tsv.gz',sep='\t')
    mc=pd.read_csv(mroot/'coordinates.tsv.gz',sep='\t')
    assert np.array_equal(hc.cell_id,x.iloc[:n_h].cell_id)
    assert np.array_equal(mc.cell_id,x.iloc[n_h:].cell_id)
    W=sparse.block_diag([H,M],format='csr')
    B=sparse.load_npz(BASE/'cross_species_expression_mnn.npz').tocsr()
    original=sparse.load_npz(BASE/'integrated_time_graph.npz').tocsr()
    reconstructed=sym_scale(normalize(.65*W+.35*B))
    diff=original-reconstructed
    maximum=float(np.max(np.abs(diff.data))) if diff.nnz else 0.0
    assert maximum<2e-6, ('Current graph reconstruction failed',maximum)
    del diff,reconstructed,original
    write_json(out/'reconstruction.json',dict(status='PASS',maximum_difference=maximum,
        native_graph_cell_orders_verified=True))
    sources=[hroot/'layout/time_graph.npz',mroot/'layout_graph.npz',BASE/'cross_species_expression_mnn.npz']
    # Remove the temporal kernel before normalization; do not just erase trajectory arrows.
    no_time=[]
    for name,local,cross,temporal,native in [
        ('human',hroot/'expression/local_kernel.npz',hroot/'expression/cross_kernel.npz',hroot/'layout/temporal_kernel.npz',H),
        ('mouse',MOUSE/'results_full_single_cell_joint_expression_v1/single_cells/expression_kernel.npz',
         MOUSE/'results_full_single_cell_joint_expression_v1/single_cells/cross_kernel.npz',
         mroot/'sampled_bidirectional_time_kernel.npz',M)]:
        L=sparse.load_npz(local).tocsr();C=sparse.load_npz(cross).tocsr();T=sparse.load_npz(temporal).tocsr()
        recovered=sym_scale(normalize(.4*L+.4*T+.2*C))
        error=native-recovered
        delta=float(np.max(np.abs(error.data))) if error.nnz else 0.0
        write_json(out/(name+'_native_reconstruction.json'),dict(maximum_difference=delta))
        assert delta<2e-6, (name,'Native layers fail to reconstruct frozen graph',delta)
        no_time.append(sym_scale(normalize(.4*L+.2*C)))
        sources.extend([local,cross,temporal]);del L,C,T,recovered,error
    W0=sparse.block_diag(no_time,format='csr');del no_time,H,M
    receipt(sources,out/'input_receipt.json')
    latent=np.load(BASE/'species_balanced_raw_cell_scvi_latent.npy')
    human_xy=np.load(hroot/'layout/time_coordinates.npy')
    mouse_xy=mc[['UMAP1','UMAP2']].to_numpy(np.float32)
    def standardize(z):
        c=np.median(z,axis=0);return (z-c)/max(float(np.sqrt(np.mean((z-c)**2))),1e-6)
    h=standardize(human_xy);m=standardize(mouse_xy)
    anchors=B[:n_h,n_h:].tocoo();rng=np.random.default_rng(20260910)
    take=rng.choice(anchors.nnz,min(100000,anchors.nnz),replace=False)
    rotation,_=orthogonal_procrustes(m[anchors.col[take]],h[anchors.row[take]])
    initial=np.vstack([h,m@rotation]).astype(np.float32)
    np.save(out/'frozen_native_aligned_initialization.npy',initial)
    variants=[('baseline_065_035',.65,False),('within_080_cross_020',.8,False),
              ('within_050_cross_050',.5,False),('no_cross_species',1.,False),('no_temporal',.65,True)]
    summaries=[];a,b=find_ab_params(1.10,.35)
    for name,weight,no_temporal in variants:
        progress('A02','RUNNING',stage='construct_ablation',variant=name)
        graph=sym_scale(normalize(weight*(W0 if no_temporal else W)+(1-weight)*B))
        dest=out/name;dest.mkdir(exist_ok=True)
        sparse.save_npz(dest/'graph.npz',graph)
        cross_mass=float(graph[:n_h,n_h:].sum(dtype=np.float64))*2
        total=float(graph.sum(dtype=np.float64))
        summary=dict(variant=name,within_coefficient=weight,cross_coefficient=1-weight,
            directed_edges=graph.nnz,total_weight=total,cross_weight_fraction=cross_mass/total,
            components=int(connected_components(graph,directed=False,return_labels=False)))
        summaries.append(summary);write_json(out/'graph_comparisons.json',summaries)
        seeds=[SEED,SEED+1,SEED+2]
        jobs=[('native_aligned',s) for s in seeds]+([('random',s) for s in seeds] if name=='baseline_065_035' else [])
        for init_name,seed in jobs:
            path=dest/f'layout_{init_name}_{seed}.npy'
            if path.exists():continue
            progress('A02','RUNNING',stage='full_cell_UMAP_ablation',variant=name,initialization=init_name,seed=seed,cells=len(x),epochs=250)
            xy,_=simplicial_set_embedding(data=latent,graph=graph.copy(),n_components=2,initial_alpha=1.,
                a=a,b=b,gamma=1.,negative_sample_rate=5,n_epochs=250,
                init=initial.copy() if init_name=='native_aligned' else 'random',
                random_state=np.random.RandomState(seed),metric='euclidean',metric_kwds={},
                densmap=False,densmap_kwds={},output_dens=False,parallel=True,verbose=False)
            assert xy.shape==(len(x),2) and np.isfinite(xy).all()
            np.save(path,xy)
        del graph;gc.collect()
    finish('A02','weight_time_initialization_ablations',variants=5,layouts=18,
        next_steps=['Expression-only graph baseline','Time-label permutation with OT refitting','A01/A04 independent evaluation'],
        boundary='No temporal-permutation result is implied by no-temporal ablation.')


if __name__=='__main__': execute('A02',run)
