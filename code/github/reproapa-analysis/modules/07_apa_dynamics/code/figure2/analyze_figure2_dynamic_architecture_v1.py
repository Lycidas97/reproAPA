"""Quantify trajectory architecture for Figure 2 without changing its visualization."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


P=Path(__file__).resolve().parents[1]
OUT=P/'outputs/figure2-dynamic-architecture-v1'
APA=P/'data/apa-human-mouse-joint-v1/all_route_node_positions.tsv.gz'
MATCHED=P/'data/figure2-cde-v1/matched_node_values.tsv.gz'
FIT=P/'data/apa-human-mouse-joint-v1/joint_fit.npz'
MEMBERS=P/'outputs/figure2-k3-srelative-v2/gene_membership_k3.tsv.gz'
STATES=list('SREP')
TRANSITIONS=[('S','R'),('R','E'),('E','P')]
ROUTES={'mouse':['WOT','CIP'],'human':['GSE118258']}
CLUSTERS=['C1','C2','C3']
SEED=91327


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def unit_centered(values):
    centered=values-values.mean(axis=1,keepdims=True)
    norm=np.linalg.norm(centered,axis=1,keepdims=True)
    assert np.all(norm>1e-10)
    return centered/norm


def summarize(values,keys,metrics):
    rows=[]
    for key,block in values.groupby(keys,dropna=False):
        if not isinstance(key,tuple):key=(key,)
        row=dict(zip(keys,key));row['n']=len(block)
        for metric in metrics:
            x=block[metric].dropna()
            row.update({f'{metric}_mean':float(x.mean()),f'{metric}_median':float(x.median()),f'{metric}_q25':float(x.quantile(.25)),f'{metric}_q75':float(x.quantile(.75))})
        rows.append(row)
    return pd.DataFrame(rows)


def trajectory_architecture(members,source):
    d=source[source.layer.eq('supported50')&source.node.isin(STATES)].merge(members[['species','gene_id','cluster']],on=['species','gene_id'],validate='many_to_one')
    wide=d.pivot(index=['species','route','cluster','gene_id'],columns='node',values='position_index').reindex(columns=STATES).dropna().reset_index()
    delta=wide[STATES].sub(wide.S,axis=0)
    steps=pd.DataFrame({name:(delta[b]-delta[a]).abs() for name,(a,b) in zip(['SR','RE','EP'],TRANSITIONS)})
    L=steps.sum(axis=1);D=delta.P.abs();valid=L.gt(1e-12)
    assert valid.all() and np.all(D<=L+1e-12)
    result=wide[['species','route','cluster','gene_id']].copy()
    for node in STATES:result[f'delta_{node}']=delta[node]
    for stage in ['SR','RE','EP']:result[f'abs_{stage}']=steps[stage]
    result['path_length_L']=L;result['endpoint_displacement_D']=D;result['hidden_remodeling_H']=L-D
    result['endpoint_retention_D_over_L']=D/L;result['transient_fraction_T']=1-D/L
    for stage in ['SR','RE','EP']:result[f'contribution_{stage}']=steps[stage]/L
    result['dominant_stage']=result[[f'contribution_{s}' for s in ['SR','RE','EP']]].idxmax(axis=1).str.removeprefix('contribution_')
    assert np.allclose(result[[f'contribution_{s}' for s in ['SR','RE','EP']]].sum(axis=1),1)
    metrics=['path_length_L','endpoint_displacement_D','hidden_remodeling_H','endpoint_retention_D_over_L','transient_fraction_T','contribution_SR','contribution_RE','contribution_EP']
    route=summarize(result,['species','route','cluster'],metrics)
    species=summarize(result,['species','cluster'],metrics)
    cluster=summarize(result,['cluster'],metrics)
    for frame,keys in [(route,['species','route','cluster']),(species,['species','cluster']),(cluster,['cluster'])]:
        dom=result.groupby(keys+['dominant_stage']).size().rename('dominant_n').reset_index()
        totals=result.groupby(keys).size().rename('total').reset_index();dom=dom.merge(totals,on=keys);dom['dominant_fraction']=dom.dominant_n/dom.total
        pivot=dom.pivot(index=keys,columns='dominant_stage',values='dominant_fraction').fillna(0).add_prefix('dominant_fraction_').reset_index()
        frame[pivot.columns.difference(keys)]=frame[keys].merge(pivot,on=keys)[pivot.columns.difference(keys)].to_numpy()
    return result,route,species,cluster


def loo_route_validation(members,source):
    z=np.load(FIT,allow_pickle=True);ids=z['gene_ids'].astype(str);species=z['species'].astype(str);features=z['features'];weights=z['weights']
    lookup={(s,g):i for i,(s,g) in enumerate(zip(species,ids))}
    ordered=members.set_index(['species','gene_id']).loc[list(zip(species,ids))].cluster.to_numpy()
    full_centers=np.array([np.average(features[ordered==c],axis=0,weights=weights[ordered==c]) for c in CLUSTERS])
    route_features={}
    for sp,routes in ROUTES.items():
        spids=ids[species==sp]
        for route in routes:
            raw=source[source.species.eq(sp)&source.route.eq(route)&source.layer.eq('supported50')].pivot(index='gene_id',columns='node',values='position_index').reindex(index=spids,columns=STATES).to_numpy()
            assert np.isfinite(raw).all()
            route_features[(sp,route)]=unit_centered(raw)
    reconstructed=np.concatenate([np.mean([route_features[(sp,r)] for r in ROUTES[sp]],axis=0) for sp in ['mouse','human']])
    assert np.allclose(reconstructed,features,atol=1e-14)
    rows=[];cluster_rows=[]
    for held_sp,held_route in [(sp,r) for sp in ['mouse','human'] for r in ROUTES[sp]]:
        train_parts=[];train_species=[];train_global=[]
        for sp in ['mouse','human']:
            remaining=[r for r in ROUTES[sp] if not (sp==held_sp and r==held_route)]
            if not remaining:continue
            block=np.mean([route_features[(sp,r)] for r in remaining],axis=0)
            train_parts.append(block);train_species.extend([sp]*len(block));train_global.extend(np.flatnonzero(species==sp))
        train=np.concatenate(train_parts);train_species=np.array(train_species);train_global=np.array(train_global)
        if len(set(train_species))==2:
            w=np.array([.5/np.sum(train_species==sp) for sp in train_species])
        else:w=np.repeat(1/len(train),len(train))
        model=KMeans(n_clusters=3,random_state=SEED,n_init=50,max_iter=500,tol=1e-7).fit(train,sample_weight=w)
        aa,bb=linear_sum_assignment(((full_centers[:,None]-model.cluster_centers_[None])**2).sum(2));old_to_main=dict(zip(bb,aa))
        centers=np.array([model.cluster_centers_[bb[np.where(aa==j)[0][0]]] for j in range(3)])
        train_assigned=np.array([old_to_main[x] for x in model.labels_])
        held=route_features[(held_sp,held_route)];projected=((held[:,None]-centers[None])**2).sum(2).argmin(1)
        full_labels=np.array([CLUSTERS.index(c) for c in ordered[species==held_sp]])
        same_species_train=np.flatnonzero(train_species==held_sp)
        agreement_training=np.nan;ari_training=np.nan
        if len(same_species_train):
            training_for_held=train_assigned[same_species_train]
            agreement_training=float(np.mean(projected==training_for_held));ari_training=float(adjusted_rand_score(projected,training_for_held))
        rows.append(dict(heldout_species=held_sp,heldout_route=held_route,training_species='+'.join(sorted(set(train_species))),n_heldout=len(held),
            assignment_agreement_vs_full=float(np.mean(projected==full_labels)),ARI_vs_full=float(adjusted_rand_score(projected,full_labels)),
            assignment_agreement_vs_remaining_route=agreement_training,ARI_vs_remaining_route=ari_training,
            training_inertia=float(model.inertia_)))
        for j,c in enumerate(CLUSTERS):
            use=full_labels==j;mean=held[use].mean(0);r=float(np.corrcoef(mean,centers[j])[0,1]);rmse=float(np.sqrt(np.mean((mean-centers[j])**2)))
            nearest=int(((mean[None]-centers)**2).sum(1).argmin())
            cluster_rows.append(dict(heldout_species=held_sp,heldout_route=held_route,cluster=c,n=int(use.sum()),centroid_correlation=r,centroid_RMSE=rmse,
                nearest_LOO_cluster=CLUSTERS[nearest],shape_recovered=nearest==j))
    # Direct support using each route against the full taxonomy centers.
    support=[]
    for (sp,route),held in route_features.items():
        labels=np.array([CLUSTERS.index(c) for c in ordered[species==sp]])
        projected=((held[:,None]-full_centers[None])**2).sum(2).argmin(1)
        for j,c in enumerate(CLUSTERS):
            mean=held[labels==j].mean(0);nearest=int(((mean[None]-full_centers)**2).sum(1).argmin())
            support.append(dict(species=sp,route=route,cluster=c,n=int((labels==j).sum()),centroid_correlation=float(np.corrcoef(mean,full_centers[j])[0,1]),
                centroid_RMSE=float(np.sqrt(np.mean((mean-full_centers[j])**2))),nearest_full_cluster=CLUSTERS[nearest],shape_supported=nearest==j,
                gene_assignment_agreement=float(np.mean(projected[labels==j]==j))))
    return pd.DataFrame(rows),pd.DataFrame(cluster_rows),pd.DataFrame(support)


def apa_gex_correspondence(members,matched):
    d=matched[matched.layer.eq('supported50')&matched.node.isin(STATES)].drop(columns=['cluster','row_order'],errors='ignore').merge(members[['species','gene_id','cluster']],on=['species','gene_id'],validate='many_to_one')
    gene_rows=[];transition_rows=[]
    for key,block in d.groupby(['species','route','cluster','gene_id']):
        sp,route,cluster,gene=key
        a=block.pivot(index='gene_id',columns='node',values='position_index').reindex(columns=STATES).iloc[0]
        e=block.pivot(index='gene_id',columns='node',values='gex_log2_cpm1').reindex(columns=STATES).iloc[0]
        ok=a.notna()&e.notna()
        nodes=[s for s in STATES if ok[s]]
        if 'S' in nodes and len(nodes)>=3:
            x=a[nodes].to_numpy(float)-a.S;y=e[nodes].to_numpy(float)-e.S
            valid=np.ptp(x)>1e-12 and np.ptp(y)>1e-12
            gene_rows.append(dict(species=sp,route=route,cluster=cluster,gene_id=gene,n_states=len(nodes),states=''.join(nodes),r=float(np.corrcoef(x,y)[0,1]) if valid else np.nan,status='estimable' if valid else 'constant'))
        for stage,(left,right) in zip(['SR','RE','EP'],TRANSITIONS):
            if ok[left] and ok[right]:
                da=float(a[right]-a[left]);de=float(e[right]-e[left]);product=da*de
                transition_rows.append(dict(species=sp,route=route,cluster=cluster,gene_id=gene,transition=stage,delta_APA=da,delta_GEX=de,
                    concordance='same' if product>0 else 'opposite' if product<0 else 'zero'))
    genes=pd.DataFrame(gene_rows);transitions=pd.DataFrame(transition_rows)
    summaries=[]
    for keys in [['species'],['species','cluster'],['species','route'],['species','route','cluster']]:
        for key,block in genes.groupby(keys):
            if not isinstance(key,tuple):key=(key,)
            valid=block.r.dropna();row=dict(zip(keys,key));row.update(n_total=len(block),n_estimable=len(valid),median_r=float(valid.median()),q25_r=float(valid.quantile(.25)),q75_r=float(valid.quantile(.75)),
                fraction_r_ge_0_5=float(valid.ge(.5).mean()),fraction_r_le_minus_0_5=float(valid.le(-.5).mean()),fraction_abs_r_lt_0_5=float(valid.abs().lt(.5).mean()))
            summaries.append(row)
    transition_summary=transitions.groupby(['species','cluster','transition','concordance']).size().rename('n').reset_index()
    totals=transition_summary.groupby(['species','cluster','transition']).n.transform('sum');transition_summary['fraction']=transition_summary.n/totals
    return genes,pd.DataFrame(summaries),transitions,transition_summary


def main():
    OUT.mkdir(exist_ok=True)
    members=pd.read_csv(MEMBERS,sep='\t')
    source=pd.read_csv(APA,sep='\t')
    matched=pd.read_csv(MATCHED,sep='\t')
    architecture,route_summary,species_summary,cluster_summary=trajectory_architecture(members,source)
    endpoint_thresholds=architecture.groupby(['species','route','cluster']).agg(n=('gene_id','size'),
        fraction_T_ge_0_5=('transient_fraction_T',lambda x:float(x.ge(.5).mean())),fraction_T_ge_0_75=('transient_fraction_T',lambda x:float(x.ge(.75).mean())),
        median_L=('path_length_L','median'),median_D=('endpoint_displacement_D','median'),median_H=('hidden_remodeling_H','median')).reset_index()
    loo,loo_clusters,support=loo_route_validation(members,source)
    gene_corr,corr_summary,transitions,transition_summary=apa_gex_correspondence(members,matched)
    outputs=[(architecture,'trajectory_gene_route_metrics.tsv.gz'),(route_summary,'trajectory_route_cluster_summary.tsv'),(species_summary,'trajectory_species_cluster_summary.tsv'),(cluster_summary,'trajectory_cluster_summary.tsv'),(endpoint_thresholds,'endpoint_information_loss_threshold_summary.tsv'),
        (loo,'leave_one_route_out_summary.tsv'),(loo_clusters,'leave_one_route_out_cluster_metrics.tsv'),(support,'route_archetype_support.tsv'),
        (gene_corr,'APA_GEX_gene_route_correlations.tsv.gz'),(corr_summary,'APA_GEX_correlation_summary.tsv'),(transitions,'APA_GEX_transition_concordance.tsv.gz'),(transition_summary,'APA_GEX_transition_summary.tsv')]
    for frame,name in outputs:frame.to_csv(OUT/name,sep='\t',index=False)
    qa=dict(status='PASS',taxonomy='fixed original centered-feature species-balanced joint k=3',complete_gene_route_trajectories=len(architecture),complete_routes=sorted(architecture.route.unique()),
        triangle_inequality=True,stage_contributions_sum_to_one=True,LOO_routes=len(loo),APA_GEX_gene_routes=len(gene_corr),APA_GEX_estimable=int(gene_corr.r.notna().sum()),APA_GEX_transitions=len(transitions),
        sources={str(f.relative_to(P)):sha(f) for f in [APA,MATCHED,FIT,MEMBERS]},script_sha256=sha(Path(__file__)))
    (OUT/'analysis_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
    print('TRAJECTORY CLUSTER SUMMARY')
    print(cluster_summary.to_string(index=False))
    print('\nLEAVE ONE ROUTE OUT')
    print(loo.to_string(index=False))
    print('\nROUTE SHAPE SUPPORT')
    print(support.groupby('cluster').agg(routes=('shape_supported','size'),supported=('shape_supported','sum'),minimum_r=('centroid_correlation','min'),maximum_RMSE=('centroid_RMSE','max')).to_string())
    print('\nAPA-GEX')
    print(corr_summary[corr_summary.columns.intersection(['species','cluster','n_total','n_estimable','median_r','q25_r','q75_r','fraction_r_ge_0_5','fraction_r_le_minus_0_5'])].dropna(subset=['cluster']).to_string(index=False))


if __name__=='__main__':main()
