"""New APA class scores on the unchanged expression graph and embeddings."""
import argparse

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from downstream_runtime import Context,Study,STUDIES,bh,seed

SOURCE='../科研/projects/2026-09-02_mouse-time-aware-integration-pilot/scripts/'


def scores(study,members):
    ctx=study.ctx;meta=ctx.read(study.root/'metacells.tsv.gz')
    counts=sparse.load_npz(ctx.register(study.root/'metacell_by_pas_counts.npz')).tocsc()
    axes=[study.axis(g) for g in members.gene_id if g in study.genes]
    axes=[a for a in axes if len(a)>=2]
    if axes:
        f=pd.concat(axes,ignore_index=True);genes=pd.Index(f.primary_gene_id.unique());gi=genes.get_indexer(f.primary_gene_id)
        c=counts[:,f.feature_index_1based.to_numpy()-1].tocsr()
        mapping=sparse.csr_matrix((np.ones(len(f)),(np.arange(len(f)),gi)),shape=(len(f),len(genes)))
        weighted=sparse.csr_matrix((f.position,(np.arange(len(f)),gi)),shape=mapping.shape)
        total=(c@mapping).toarray();numerator=(c@weighted).toarray()
        values=np.divide(numerator,total,out=np.full(total.shape,np.nan),where=total>=20)
    else:genes=pd.Index([],dtype=str);values=np.empty((len(meta),0))
    primary=study.f.keep_for_apa_primary.astype(str).str.lower().eq('true')
    primary_umi=np.asarray(counts[:,np.flatnonzero(primary)].sum(axis=1)).ravel()
    return meta,genes,values,primary_umi


def moran(ctx):
    math=ctx.functions(SOURCE+'test_figure2_apa_manifold_coherence.py',['permuted_moran'])
    graphmath=ctx.functions(SOURCE+'robustness_figure2_apa_knn.py',['graph_for_study'],dict(cKDTree=cKDTree))
    nodes=ctx.read(ctx.resource('expression_graph_nodes'));latent=np.load(ctx.resource('expression_latent'))
    frozen=sparse.load_npz(ctx.resource('expression_graph')).tocsr();frozen=frozen.maximum(frozen.T)
    frozen.data[:]=1;frozen.setdiag(0);frozen.eliminate_zeros()
    states=ctx.read(ctx.resource('frozen_metacell_states'))
    metadata=nodes.merge(states,on=['metacell_id','library_key'],validate='one_to_one').set_index('metacell_id')
    assert len(metadata)==len(nodes)==len(latent) and frozen.shape==(len(nodes),len(nodes))
    members=ctx.read(ctx.analysis/'profiles/gene_membership_k3.tsv.gz');index=pd.Index(nodes.metacell_id)
    records=[];panels=[];scoreparts=[]
    def test(frame,graph,value,block,testname,context):
        if len(frame)<=20 or not graph.nnz or np.var(value)<=1e-20:
            return context|dict(test=testname,status='not_estimable',n_metacells=len(frame),empirical_p=1.)
        return context|dict(test=testname,status='tested',n_metacells=len(frame),
            **math.permuted_moran(graph,value,block,9999,seed('moran',str(context),testname)))
    for sid in STUDIES:
        study=Study(ctx,sid);species=str(study.f.species.iloc[0]);m=members[members.species.eq(species)]
        meta,genes,values,primary=scores(study,m);ids=meta.metacell_id.to_numpy();state=metadata.loc[ids].copy()
        state['primary_pas_umi']=primary;ix=index.get_indexer(ids);assert (ix>=0).all()
        means=np.divide(np.nansum(values,axis=0),np.isfinite(values).sum(axis=0),out=np.full(len(genes),np.nan),where=np.isfinite(values).sum(axis=0)>0)
        classes=m.set_index('gene_id').cluster.reindex(genes).to_numpy();score=pd.DataFrame(dict(metacell_id=ids,species=species,study_id=sid))
        graphs={k:graphmath.graph_for_study(latent[ix],k,state.study_group) for k in [10,15,20,30]}
        assert (graphs[20]-frozen[ix][:,ix]).nnz==0,'Frozen expression graph changed'
        for klass in ['C1','C2','C3']:
            take=classes==klass;observed=np.isfinite(values[:,take]);n=observed.sum(axis=1)
            centered=np.divide(np.nansum(values[:,take]-means[take],axis=1),n,out=np.full(len(ids),np.nan),where=n>0)
            raw=np.divide(np.nansum(values[:,take],axis=1),n,out=np.full(len(ids),np.nan),where=n>0)
            score[klass+'_centered']=centered;score[klass+'_wpp']=raw;score[klass+'_n']=n
            context=dict(study_id=sid,species=species,cluster=klass)
            for k,graph in graphs.items():
                chosen=np.flatnonzero(state.state_purity.ge(.8).to_numpy()&np.isfinite(centered))
                local=graph[chosen][:,chosen].tocsr();connected=np.diff(local.indptr)>0;chosen=chosen[connected]
                local=local[connected][:,connected];frame=state.iloc[chosen];value=centered[chosen]
                labels=frame.state_group.to_numpy();blocks=(frame.library_key+'|'+frame.state_group).to_numpy()
                design=pd.get_dummies(frame.state_group,dtype=float).to_numpy()
                design=np.column_stack([design,np.log1p(n[chosen]),np.log1p(primary[chosen])])
                residual=value-design@np.linalg.lstsq(design,value,rcond=None)[0] if len(frame) else value
                records.append(test(frame,local,residual,blocks,'state_coverage_library_state',context|dict(k=k)))
                if k==20:
                    state_residual=value-pd.Series(value).groupby(labels).transform('mean').to_numpy()
                    for label,x,block in [('total_unrestricted',value,None),('state_residual_within_state',state_residual,labels),
                                          ('state_residual_within_library_state',state_residual,blocks)]:
                        records.append(test(frame,local,x,block,label,context|dict(k=k)))
            group=values[:,take];coverage=np.isfinite(group).mean(axis=0)
            for threshold in [.7,.8,.9]:
                selected=coverage>=threshold
                complete=np.isfinite(group[:,selected]).all(axis=1) if selected.any() else np.zeros(len(ids),bool)
                chosen=np.flatnonzero(complete&state.state_purity.ge(.8).to_numpy())
                panels.append(context|dict(threshold=threshold,selected_genes=int(selected.sum()),complete_metacells=len(chosen),
                    selected_gene_ids=';'.join(genes[take][selected])))
                if len(chosen)<30:
                    for label in ['original_same_nodes','fixed_panel']:
                        records.append(context|dict(k=20,threshold=threshold,test=label,status='insufficient_nodes',n_metacells=len(chosen),empirical_p=1.))
                    continue
                frame=state.iloc[chosen];graph=graphmath.graph_for_study(latent[ix[chosen]],20,frame.study_group)
                columns=np.flatnonzero(take)[selected];fixed=np.mean(values[chosen][:,columns]-means[columns],axis=1)
                design=np.column_stack([pd.get_dummies(frame.state_group,dtype=float),np.log1p(n[chosen]),np.log1p(primary[chosen])])
                blocks=(frame.library_key+'|'+frame.state_group).to_numpy()
                for label,value in [('original_same_nodes',centered[chosen]),('fixed_panel',fixed)]:
                    residual=value-design@np.linalg.lstsq(design,value,rcond=None)[0]
                    records.append(test(frame,graph,residual,blocks,label,context|dict(k=20,threshold=threshold)))
        scoreparts.append(score);print('moran',sid,flush=True)
    result=pd.DataFrame(records);result['family']=result.test+'|k'+result.k.astype(str)+'|panel'+result.threshold.fillna(0).astype(str)
    core=result.test.isin(['total_unrestricted','state_residual_within_state']);result.loc[core,'family']='primary_total_and_state'
    result['q_BH']=result.groupby('family').empirical_p.transform(bh)
    ctx.table(result,'moran_results.tsv');ctx.table(panels,'fixed_panel_feasibility.tsv')
    ctx.table(pd.concat(scoreparts),'k3_group_wpp_metacells.tsv.gz')
    ctx.finish(dict(frozen_expression_graph_recovered=True,new_class_scores=True,new_primary_PAS_coverage=True,
        library_state_exchangeability=True,separate_BH_families=True,fixed_gene_panel_controls=True),permutations=9999)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();moran(Context(a.config,'moran',a.out))
