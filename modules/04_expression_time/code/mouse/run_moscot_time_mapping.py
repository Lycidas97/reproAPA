#!/usr/bin/env python3
import argparse, json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from moscot.problems.generic import SinkhornProblem
from moscot.problems.time import TemporalProblem

QUERIES=["gse103221_osk","gse103221_cip","gse171127_d3_getm","gse171127_d3_oskm","gse171127_d6_pooled","gse171127_d12_pooled"]
REFERENCE="wot_common_to_2i"

def progress(obs):
    out=np.zeros(len(obs)); ds=obs.dataset.astype(str).to_numpy(); day=obs.day.to_numpy(float)
    for route in np.unique(ds):
        ix=np.flatnonzero(ds==route)
        denom=12.0 if route.startswith("gse171127") else max(day[ix].max(),1.0)
        out[ix]=day[ix]/denom
    return np.clip(out,0,1)

def choose_lambda(z,obs,p,candidates=(0,.25,1,4,9)):
    target=np.flatnonzero(obs.dataset.to_numpy()==REFERENCE); rows=[]
    base={}
    for route in QUERIES:
        src=np.flatnonzero(obs.dataset.to_numpy()==route)
        for lam in candidates:
            aug_t=np.c_[z[target],np.sqrt(lam)*p[target,None]]; aug_s=np.c_[z[src],np.sqrt(lam)*p[src,None]]
            ind=NearestNeighbors(n_neighbors=1).fit(aug_t).kneighbors(aug_s,return_distance=False)[:,0]; match=target[ind]
            expr=np.sqrt(((z[src]-z[match])**2).sum(1)); pm=np.abs(p[src]-p[match])
            if lam==0: base[route]=np.median(expr)
            rows.append({"route":route,"lambda":lam,"progress_mae":np.mean(pm),"expression_distance_median":np.median(expr)})
    tab=pd.DataFrame(rows); agg=tab.groupby("lambda").agg(progress_mae=("progress_mae","mean"),expression_distance_median=("expression_distance_median","mean")).reset_index()
    baseline=agg.loc[agg["lambda"].eq(0),"expression_distance_median"].iloc[0]
    ok=agg[(agg.progress_mae<=.15)&(agg.expression_distance_median<=1.25*baseline)]
    if len(ok): chosen=float(ok.sort_values("lambda").iloc[0]["lambda"])
    else:
        best=(agg.progress_mae/agg.progress_mae.max()+agg.expression_distance_median/baseline).idxmin()
        chosen=float(agg.loc[best,"lambda"])
    return chosen,tab,agg

def ood_scores(z,obs,p):
    target=np.flatnonzero(obs.dataset.to_numpy()==REFERENCE); ref=[]
    for i in target:
        cand=target[(np.abs(p[target]-p[i])<=.08)&(target!=i)]
        if len(cand): ref.append(np.sqrt(((z[cand]-z[i])**2).sum(1)).min())
    ref=np.sort(ref); result=np.full(len(obs),np.nan)
    for route in QUERIES:
        for i in np.flatnonzero(obs.dataset.to_numpy()==route):
            cand=target[np.abs(p[target]-p[i])<=.15]
            if not len(cand): cand=target
            d=np.sqrt(((z[cand]-z[i])**2).sum(1)).min()
            result[i]=np.searchsorted(ref,d,side="right")/len(ref)
    return result,float(np.quantile(ref,.99))

def add_edges(mat,src,tgt,rows,cols,data,weight,topk=5):
    x=np.array(mat,copy=True); x/=np.maximum(x.sum(1,keepdims=True),1e-300)
    for i in range(x.shape[0]):
        take=np.argpartition(x[i],-min(topk,x.shape[1]))[-min(topk,x.shape[1]):]
        scale=max(x[i,take].max(),1e-12)
        for j in take: rows.extend([src[i],tgt[j]]); cols.extend([tgt[j],src[i]]); data.extend([weight*x[i,j]/scale]*2)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--h5ad",required=True); ap.add_argument("--out",required=True); args=ap.parse_args(); out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    a=sc.read_h5ad(args.h5ad); a.obs["progress"]=progress(a.obs); z=a.obsm["X_scVI"].astype(float); z=(z-z.mean(0))/np.maximum(z.std(0),1e-6)
    lam,grid,agg=choose_lambda(z,a.obs,a.obs.progress.to_numpy()); grid.to_csv(out/"time_lambda_grid_by_route.tsv",sep="\t",index=False); agg.to_csv(out/"time_lambda_grid_summary.tsv",sep="\t",index=False)
    a.obsm["X_time_aug"]=np.c_[z,np.sqrt(lam)*a.obs.progress.to_numpy()[:,None]]
    pairs=[(r,REFERENCE) for r in QUERIES]
    cross=SinkhornProblem(a).prepare(key="dataset",joint_attr="X_time_aug",policy="explicit",subset=pairs).solve(epsilon=.05,tau_a=.9,tau_b=.9,max_iterations=2000,threshold=1e-3)
    rows=[]; cols=[]; vals=[]; summaries=[]; mapped_reference_progress=np.full(a.n_obs,np.nan)
    for route,target_name in pairs:
        src=np.flatnonzero(a.obs.dataset.to_numpy()==route); tgt=np.flatnonzero(a.obs.dataset.to_numpy()==target_name); t=np.asarray(cross[route,target_name].solution.transport_matrix)
        add_edges(t,src,tgt,rows,cols,vals,.25)
        rn=t.astype(np.float64,copy=False)/np.maximum(t.sum(1,keepdims=True,dtype=np.float64),1e-300); mapped=rn@a.obs.iloc[tgt].progress.to_numpy()
        mapped_reference_progress[src]=mapped
        entropy=-(rn*np.log(np.clip(rn,1e-300,None))).sum(1)/np.log(rn.shape[1])
        summaries.append({"route":route,"n_cells":len(src),"mapped_progress_mae":float(np.mean(np.abs(mapped-a.obs.iloc[src].progress))),"mean_mapping_entropy":float(np.mean(entropy))})
    temporal_rows=[]
    for route in [REFERENCE,"gse103221_osk","gse103221_cip"]:
        ix=np.flatnonzero(a.obs.dataset.to_numpy()==route); sub=a[ix].copy()
        tp=TemporalProblem(sub).prepare(time_key="day",joint_attr="X_scVI").solve(epsilon=.05,tau_a=.95,tau_b=.95,max_iterations=2000,threshold=1e-3)
        for (d0,d1),prob in tp.problems.items():
            aa=ix[np.flatnonzero(sub.obs.day.to_numpy()==d0)]; bb=ix[np.flatnonzero(sub.obs.day.to_numpy()==d1)]; t=np.asarray(prob.solution.transport_matrix)
            add_edges(t,aa,bb,rows,cols,vals,.15,topk=3); temporal_rows.append({"route":route,"day_from":d0,"day_to":d1,"mass":float(t.sum())})
    pd.DataFrame(summaries).to_csv(out/"cross_study_mapping_summary.tsv",sep="\t",index=False); pd.DataFrame(temporal_rows).to_csv(out/"moscot_temporal_edges.tsv",sep="\t",index=False)
    a.obs["mapped_reference_progress"]=mapped_reference_progress
    a.obs["ood_probability"],threshold=ood_scores(z,a.obs,a.obs.progress.to_numpy()); a.obs["mapping_status"]=np.where(a.obs.ood_probability>=.99,"OOD",np.where(a.obs.ood_probability.notna(),"mapped","reference"))
    query=a.obs[a.obs.dataset.astype(str).isin(QUERIES)].copy(); query["is_ood"]=query.mapping_status.eq("OOD")
    route_day=(query.groupby(["dataset","day"],observed=True).agg(n_cells=("mapped_reference_progress","size"),observed_progress=("progress","mean"),mapped_progress_mean=("mapped_reference_progress","mean"),mapped_progress_sd=("mapped_reference_progress","std"),mapping_mae=("mapped_reference_progress",lambda x: float(np.mean(np.abs(x-query.loc[x.index,"progress"])))),ood_fraction=("is_ood","mean")).reset_index())
    route_day.to_csv(out/"route_day_mapping_summary.tsv",sep="\t",index=False)
    graph=a.obsp["connectivities"].maximum(sparse.csr_matrix((vals,(rows,cols)),shape=(a.n_obs,a.n_obs)))
    a.obsp["time_connectivities"]=graph; a.uns["time_neighbors"]={"connectivities_key":"time_connectivities","distances_key":"distances","params":{"n_neighbors":30,"method":"umap","metric":"euclidean","use_rep":"X_scVI"}}
    a.write_h5ad(out/"time_graph_intermediate.h5ad",compression="gzip")
    sc.tl.umap(a,neighbors_key="time_neighbors",key_added="time_regularized",random_state=20260902)
    xy=a.obsm["time_regularized"]; coords=a.obs.copy(); coords[["UMAP1","UMAP2"]]=xy; coords.to_csv(out/"time_regularized_umap.tsv.gz",sep="\t",index=False)
    fig,axs=plt.subplots(2,2,figsize=(14,11))
    for ax,col in zip(axs.flat,["study_id","dataset","day","mapping_status"]):
        if col=="day": q=ax.scatter(xy[:,0],xy[:,1],c=a.obs.day,s=3,cmap="viridis",rasterized=True); fig.colorbar(q,ax=ax,label="day")
        else: sns.scatterplot(x=xy[:,0],y=xy[:,1],hue=a.obs[col],s=5,linewidth=0,ax=ax,legend="brief")
        ax.set_title(col)
    fig.tight_layout(); fig.savefig(out/"time_regularized_overview.png",dpi=180); plt.close(fig)
    selected=agg.loc[agg["lambda"].eq(lam)].iloc[0]; baseline=agg.loc[agg["lambda"].eq(0),"expression_distance_median"].iloc[0]
    summary={"selected_time_lambda":lam,"lambda_acceptance_met":bool(selected.progress_mae<=.15 and selected.expression_distance_median<=1.25*baseline),"ood_reference_99pct_distance":threshold,"ood_cells":int((a.obs.mapping_status=="OOD").sum()),"moscot_version":"0.5.0","time_used_in_scvi":False,"time_used_in_ot":True}
    (out/"time_mapping_summary.json").write_text(json.dumps(summary,indent=2)+"\n")

if __name__=="__main__": main()
