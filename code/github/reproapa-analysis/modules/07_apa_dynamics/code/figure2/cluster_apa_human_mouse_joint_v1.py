"""Fit one species-balanced shape dictionary on complete native SREP profiles."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_samples
from scipy.optimize import linear_sum_assignment
from compute_apa_node_position import sha

SPECIES=['mouse','human']
FIT_ROUTES={'mouse':['WOT','CIP'],'human':['GSE118258']}
STATES=list('SREP')
SEED=91327


def standardize(v):
    assert np.isfinite(v).all(), 'Complete observed SREP is required; no imputation'
    z=v-v.mean(axis=1,keepdims=True);norm=np.linalg.norm(z,axis=1)
    assert (norm>1e-10).all(), 'A fitted profile is flat'
    return z/norm[:,None]


def features(views):
    return np.concatenate([np.mean([standardize(v) for v in views[sp]],axis=0) for sp in SPECIES])


def weights(species):
    return np.array([.5/np.sum(species==s) for s in species])


def fit(x,species,k,seed=SEED):
    model=KMeans(n_clusters=k,random_state=seed,n_init=50,max_iter=500,tol=1e-7)
    model.fit(x,sample_weight=weights(species))
    return model


def match(reference,other,k):
    table=np.zeros((k,k),int);np.add.at(table,(reference,other),1)
    a,b=linear_sum_assignment(-table);mapping=dict(zip(b,a));mapped=np.array([mapping[v] for v in other])
    jac=[np.sum((reference==j)&(mapped==j))/np.sum((reference==j)|(mapped==j)) for j in range(k)]
    return mapped,dict(ari=float(adjusted_rand_score(reference,other)),minimum_jaccard=float(min(jac)),jaccard=jac)


def views_from_frame(frame,members,layer):
    return {sp:[frame[frame.species.eq(sp)&frame.layer.eq(layer)&frame.route.eq(r)].pivot(index='gene_id',columns='node',values='position_index').reindex(index=members[sp].gene_id,columns=STATES).to_numpy() for r in FIT_ROUTES[sp]] for sp in SPECIES}


def main(root,old,out,ortholog):
    assert not out.exists();out.mkdir(parents=True)
    source=pd.read_csv(root/'all_route_node_positions.tsv.gz',sep='\t')
    members={sp:pd.read_csv(root/f'{sp}_gene_membership.tsv',sep='\t').query("cluster != 'U'").sort_values('gene_id').reset_index(drop=True) for sp in SPECIES}
    species=np.concatenate([np.repeat(sp,len(members[sp])) for sp in SPECIES]);ids=np.concatenate([members[sp].gene_id for sp in SPECIES])
    views=views_from_frame(source,members,'supported50');x=features(views);n=len(x)
    expanded=features(views_from_frame(source,members,'expanded_primary'))
    assert np.allclose(x.sum(1),0) and np.isclose(weights(species).sum(),1)
    assert all(np.isclose(weights(species)[species==sp].sum(),.5) for sp in SPECIES)
    models={k:fit(x,species,k) for k in range(2,7)};parts={k:[] for k in models};libparts={k:[] for k in models}
    metrics={k:[] for k in models};libmetrics={k:[] for k in models};inputs={};library_weights={};hashes={}
    for sp in SPECIES:
        path=old/sp/'library_inputs.npz';data=np.load(path);assert np.array_equal(data['gene_ids'],members[sp].gene_id)
        arrays=[data[f'values_{j}'] for j in range(len(FIT_ROUTES[sp]))];indices=[data[f'indices_{j}'] for j in range(len(arrays))]
        for a,v in zip(arrays,views[sp]):assert np.allclose(np.nanmean(a,axis=0),v,equal_nan=True,atol=1e-12)
        inputs[sp]=(arrays,indices);library_weights[sp]=np.load(old/sp/'validation.npz')['weights'];hashes[str(path)]=sha(path);hashes[str(old/sp/'validation.npz')]=sha(old/sp/'validation.npz')
    rng=np.random.default_rng(SEED+1)
    for rep in range(24):
        ix=np.sort(np.concatenate([rng.choice(np.flatnonzero(species==sp),int(np.ceil(np.sum(species==sp)*.8)),replace=False) for sp in SPECIES]))
        for k,model in models.items():
            lab=fit(x[ix],species[ix],k).labels_;mapped,metric=match(model.labels_[ix],lab,k)
            full=np.full(n,-1,int);full[ix]=mapped;parts[k].append(full);metrics[k].append(metric)
    for rep in range(32):
        vv={}
        for sp in SPECIES:
            arrays,indices=inputs[sp];vv[sp]=[]
            for a,ix in zip(arrays,indices):
                w=library_weights[sp][rep,ix,None,None];finite=np.isfinite(a);den=(finite*w).sum(0)
                vv[sp].append(np.divide((np.where(finite,a,0)*w).sum(0),den,out=np.full(den.shape,np.nan),where=den>0))
        xx=features(vv)
        for k,model in models.items():
            lab=fit(xx,species,k).labels_;mapped,metric=match(model.labels_,lab,k);libparts[k].append(mapped);libmetrics[k].append(metric)
    rows=[];extra={}
    for k,model in models.items():
        lab=model.labels_;exmap,exmetric=match(lab,fit(expanded,species,k).labels_,k);extra[k]=exmap
        _,startmetric=match(lab,fit(x,species,k,SEED+100).labels_,k)
        sil=silhouette_samples(x,lab)
        rows.append(dict(k=k,selected_display=k==4,species_balanced_silhouette=float(np.average(sil,weights=weights(species))),inertia=float(model.inertia_),minimum_size=int(np.bincount(lab).min()),
            minimum_mouse_size=int(np.bincount(lab[species=='mouse'],minlength=k).min()),minimum_human_size=int(np.bincount(lab[species=='human'],minlength=k).min()),
            subsample_median_ari=float(np.median([v['ari'] for v in metrics[k]])),subsample_median_minimum_jaccard=float(np.median([v['minimum_jaccard'] for v in metrics[k]])),
            library_median_ari=float(np.median([v['ari'] for v in libmetrics[k]])),library_median_minimum_jaccard=float(np.median([v['minimum_jaccard'] for v in libmetrics[k]])),expanded_ari=exmetric['ari'],different_seed_ari=startmetric['ari']))
    pd.DataFrame(rows).to_csv(out/'resolution_comparison.tsv',sep='\t',index=False)
    model=models[4];order=np.argsort(model.cluster_centers_[:,3]-model.cluster_centers_[:,0]);newlab=np.empty(4,int);newlab[order]=np.arange(4);labels=newlab[model.labels_];centers=model.cluster_centers_[order]
    center_rows=[];agreement=[];offset=0
    for sp in SPECIES:
        m=members[sp].copy();nn=len(m);sl=slice(offset,offset+nn);labs=labels[sl];m['previous_cluster']=m.cluster;m['cluster']=[f'C{j+1}' for j in labs]
        m['joint_center_distance']=np.linalg.norm(x[sl]-centers[labs],axis=1)
        sub=np.stack(parts[4])[:,sl];m['joint_subsample_retention']=np.sum(sub==model.labels_[sl],axis=0)/np.sum(sub>=0,axis=0)
        m['joint_library_retention']=np.mean(np.stack(libparts[4])[:,sl]==model.labels_[sl],axis=0);m['joint_expanded_agreement']=extra[4][sl]==model.labels_[sl]
        # Old sensitivity fields describe a different fit and must not be mistaken for new evidence.
        m=m.drop(columns=['algorithm_retention','library_cluster_retention','expanded_group_agreement','old_cluster','k'],errors='ignore')
        m=m.sort_values(['cluster','joint_center_distance','gene_id']);m['row_order']=np.arange(nn);m.to_csv(out/f'{sp}_gene_membership.tsv',sep='\t',index=False)
        for j in range(4):
            mu=x[sl][labs==j].mean(0)
            center_rows.append(dict(species=sp,cluster=f'C{j+1}',genes=int(np.sum(labs==j)),**dict(zip(STATES,mu))))
        offset+=nn
    for j,c in enumerate(centers):center_rows.append(dict(species='joint_center',cluster=f'C{j+1}',genes=int(np.sum(labels==j)),**dict(zip(STATES,c))))
    center_table=pd.DataFrame(center_rows);center_table.to_csv(out/'shape_centers.tsv',sep='\t',index=False)
    for j in range(4):
        a=center_table.query("species == 'mouse'").iloc[j][STATES].to_numpy(float);b=center_table.query("species == 'human'").iloc[j][STATES].to_numpy(float)
        agreement.append(dict(cluster=f'C{j+1}',human_mouse_mean_shape_r=float(np.corrcoef(a,b)[0,1])))
    pd.DataFrame(agreement).to_csv(out/'human_mouse_shape_agreement.tsv',sep='\t',index=False)
    # Score native >=3-node profiles against the same learned centers; do not fill a missing node.
    evidence=[];offset=0
    for sp in SPECIES:
        mm=members[sp];ll=labels[offset:offset+len(mm)];offset+=len(mm)
        for route,part in source[source.species.eq(sp)&source.layer.eq('supported50')&source.node.isin(STATES)].groupby('route'):
            nodes=[s for s in STATES if part[part.node.eq(s)].position_index.notna().any()]
            if len(nodes)<3:continue
            v=part.pivot(index='gene_id',columns='node',values='position_index').reindex(index=mm.gene_id,columns=nodes).to_numpy();ok=np.isfinite(v).all(1)&(np.ptp(v,axis=1)>=.01);ix=np.flatnonzero(ok)
            zz=standardize(v[ok]);cc=standardize(centers[:,[STATES.index(s) for s in nodes]]);rho=zz@cc.T;best=rho.argmax(1)
            for t,i in enumerate(ix):evidence.append(dict(species=sp,route=route,gene_id=mm.gene_id.iloc[i],cluster=f'C{ll[i]+1}',native_states=''.join(nodes),used_for_shape_fit=route in FIT_ROUTES[sp],assigned_center_r=float(rho[t,ll[i]]),nearest_center=f'C{best[t]+1}',nearest_agrees=bool(best[t]==ll[i]),positive_shape_support=bool(rho[t,ll[i]]>=.6)))
    pd.DataFrame(evidence).to_csv(out/'native_shape_evidence.tsv.gz',sep='\t',index=False)
    e=pd.DataFrame(evidence);summary=e.groupby(['species','route','cluster','native_states','used_for_shape_fit']).agg(eligible_genes=('gene_id','size'),median_assigned_center_r=('assigned_center_r','median'),nearest_agreement=('nearest_agrees','mean'),fraction_r_ge_06=('positive_shape_support','mean')).reset_index();summary.to_csv(out/'native_shape_summary.tsv',sep='\t',index=False)
    ortho=pd.read_csv(ortholog,sep='\t');ortho=ortho[ortho.documented_one_to_one.eq(True)]
    for sp in SPECIES:
        m=pd.read_csv(out/f'{sp}_gene_membership.tsv',sep='\t');ortho=ortho.merge(m[['gene_id','cluster']].rename(columns={'gene_id':f'{sp}_gene_id','cluster':f'{sp}_cluster'}),on=f'{sp}_gene_id')
    assert ortho.mouse_gene_id.is_unique and ortho.human_gene_id.is_unique
    ortho['same_joint_shape']=ortho.mouse_cluster.eq(ortho.human_cluster);ortho.to_csv(out/'ortholog_shape_pairs.tsv',sep='\t',index=False)
    # Independent assignment / center / objective checks encode the balanced fitting contract.
    dist=((x[:,None,:]-model.cluster_centers_[None,:,:])**2).sum(2);assert np.array_equal(dist.argmin(1),model.labels_)
    for j in range(4):assert np.allclose(model.cluster_centers_[j],np.average(x[model.labels_==j],axis=0,weights=weights(species)[model.labels_==j]),atol=1e-8)
    direct=0.;offset=0
    for sp in SPECIES:
        nn=len(members[sp]);cc=model.cluster_centers_[model.labels_[offset:offset+nn]]
        for v in views[sp]:direct+=.5/nn/len(views[sp])*np.sum((standardize(v)-cc)**2)
        offset+=nn
    constant=sum(.5/len(members[sp])*sum(np.sum((standardize(v)-np.mean([standardize(t) for t in views[sp]],axis=0))**2) for v in views[sp])/len(views[sp]) for sp in SPECIES)
    assert np.isclose(direct,model.inertia_+constant,atol=1e-10)
    assert np.array_equal(fit(x,species,4).labels_,model.labels_)
    np.savez_compressed(out/'joint_fit.npz',gene_ids=ids,species=species,features=x,centers=centers,labels=labels,weights=weights(species),
        **{f'labels_k{k}':m.labels_ for k,m in models.items()},**{f'subsamples_k{k}':np.stack(v) for k,v in parts.items()},**{f'library_labels_k{k}':np.stack(v) for k,v in libparts.items()},**{f'expanded_labels_k{k}':v for k,v in extra.items()})
    for name in ['config.json','all_route_node_positions.tsv.gz']:shutil.copy2(root/name,out/name)
    for name in ['human_unassigned_genes.tsv','mouse_unassigned_genes.tsv']:
        if (root/name).exists():shutil.copy2(root/name,out/name)
    source_hashes={str(root/name):sha(root/name) for name in ['all_route_node_positions.tsv.gz','mouse_gene_membership.tsv','human_gene_membership.tsv','config.json']};source_hashes.update(hashes);source_hashes[str(ortholog)]=sha(ortholog)
    result=dict(status='PASS',method='species-balanced KMeans of mean study-normalized complete SREP profiles',selected_k=4,k_selection='Display resolution retained from previous figure; not claimed optimal',seed=SEED,n_init=50,fit_routes=FIT_ROUTES,counts={sp:len(members[sp]) for sp in SPECIES},shape_agreement=agreement,ortholog_pairs=len(ortho),ortholog_same_shape=int(ortho.same_joint_shape.sum()),source_sha256=source_hashes,script_sha256=sha(Path(__file__)),assignment_and_center_check=True,study_level_objective_check=True,deterministic_replay=True,node_imputation=False,library_reweighting_is_biological_bootstrap=False)
    (out/'analysis_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['root','old','out','ortholog']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.root,a.old,a.out,a.ortholog)
