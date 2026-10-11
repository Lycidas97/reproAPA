"""Evaluate fixed Figure 2 members after prespecified state-member exclusions."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import sys
PYTHON='/root/miniconda3/envs/apasim_analysis/bin/python'
if os.path.realpath(sys.executable)!=os.path.realpath(PYTHON):
    os.execv(PYTHON,[PYTHON,__file__]+sys.argv[1:])
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn import __version__ as sklearn_version

ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data'
BASE=Path('/root/nfsdata/apadetect/analysis_projects');PRIOR=BASE/'reproapa_pas_reliability_2026_09_13'
MOUSE=BASE/'mouse_time_aware_integration_pilot_v1';JOINT=MOUSE/'results_apa_human_mouse_joint_v1'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(d,name):d.to_csv(D/name,sep='\t',index=False)
def sign(a):return np.where(a>=.01,1,np.where(a<=-.01,-1,0))
SPEC=json.loads((D/'apa_frozen_variants.json').read_text());VARIANTS=list(SPEC['variants'])
assert json.loads((D/'apa_reaggregation_verification.json').read_text())['status']=='PASS'
receipts={str(p):sha(p) for p in [JOINT/'joint_fit.npz',JOINT/'mouse_gene_membership.tsv',JOINT/'human_gene_membership.tsv',D/'affected_library_positions.tsv.gz',D/'apa_frozen_variants.json']}
source=np.load(JOINT/'joint_fit.npz',allow_pickle=True)
ids=source['gene_ids'];species=source['species'];oldx=source['features'];centers=source['centers'];labels=source['labels']
members={sp:read(JOINT/f'{sp}_gene_membership.tsv') for sp in ['mouse','human']}
assert {sp:len(m) for sp,m in members.items()}=={'mouse':540,'human':897}
for sp,m in members.items():
    assert set(ids[species==sp])==set(m.gene_id)
    assert (m.set_index('gene_id').reindex(ids[species==sp]).cluster.to_numpy()==np.array(['C'+str(x+1) for x in labels[species==sp]])).all()
# Original estimates from all seven studies preserve every fixed study denominator.
old=[]
for folder in sorted((PRIOR/'revision/data').iterdir()):
    p=folder/'library_positions.tsv.gz'
    if not folder.is_dir() or not p.exists():continue
    receipts[str(p)]=sha(p);z=read(p);old.append(z[z['mode'].eq('original')].drop(columns='mode'))
old=pd.concat(old,ignore_index=True);assert old.study.nunique()==7
affected=read(D/'affected_library_positions.tsv.gz')
keys=['species','study','route','node','gene_id']
original_nodes=old.groupby(keys).position_index.mean().reset_index()
nodeparts=[]
for variant in VARIANTS:
    for policy in ['reapply_purity','frozen_purity']:
        z=affected[affected.variant.eq(variant)&affected.policy.eq(policy)]
        current=pd.concat([old[~old.study.isin(SPEC['studies'])],z[old.columns]],ignore_index=True)
        n=current.groupby(keys).position_index.mean().reset_index();n['variant']=variant;n['policy']=policy;nodeparts.append(n)
nodes=pd.concat(nodeparts,ignore_index=True);save(nodes,'all_variant_node_positions.tsv.gz')
# Repeat the no-exclusion readback at node level, including the untouched studies.
b=nodes[nodes.variant.eq('baseline')&nodes.policy.eq('reapply_purity')].set_index(keys).position_index.sort_index()
e=original_nodes.set_index(keys).position_index.sort_index()
assert b.index.equals(e.index) and np.allclose(b,e,atol=1e-10,rtol=0,equal_nan=True)
routes={'mouse':['WOT','CIP'],'human':['GSE118258']}
def features(variant,policy):
    parts=[]
    for sp in ['mouse','human']:
        gs=ids[species==sp];vectors=[]
        for route in routes[sp]:
            t=nodes[nodes.species.eq(sp)&nodes.variant.eq(variant)&nodes.policy.eq(policy)&nodes.route.eq(route)]
            a=t.pivot(index='gene_id',columns='node',values='position_index').reindex(index=gs,columns=list('SREP')).to_numpy()
            centered=a-a.mean(1,keepdims=True);norm=np.linalg.norm(centered,axis=1)
            vectors.append(np.divide(centered,norm[:,None],out=np.full(a.shape,np.nan),where=norm[:,None]>1e-10))
        parts.append(np.mean(vectors,axis=0))
    return np.concatenate(parts)
baseline_features=features('baseline','reapply_purity')
assert np.allclose(baseline_features,oldx,rtol=0,atol=1e-10)
def fit(x,sp):
    w=np.array([.5/np.sum(sp==s) for s in sp])
    return KMeans(n_clusters=4,random_state=91327,n_init=50,max_iter=500,tol=1e-7).fit(x,sample_weight=w).labels_
def matched(ref,new):
    cross=np.zeros((4,4),int);np.add.at(cross,(ref,new),1)
    a,b=linear_sum_assignment(-cross);mapping=dict(zip(b,a));return np.array([mapping[x] for x in new])
shape_rows=[];shape_summary=[]
for variant in VARIANTS:
    for policy in ['reapply_purity','frozen_purity']:
        x=features(variant,policy);ok=np.isfinite(x).all(1)
        nearest=np.full(len(ids),-1,int);refit=np.full(len(ids),-1,int)
        eligible_baseline_refit=np.full(len(ids),-1,int);eligible_new_refit=np.full(len(ids),-1,int)
        nearest[ok]=((x[ok,None,:]-centers[None,:,:])**2).sum(2).argmin(1)
        if ok.sum()>=4:
            newfit=fit(x[ok],species[ok]);refit[ok]=matched(labels[ok],newfit)
            eligible_baseline_refit[ok]=fit(oldx[ok],species[ok])
            eligible_new_refit[ok]=matched(eligible_baseline_refit[ok],newfit)
        if variant=='baseline':assert ok.all() and np.array_equal(refit,labels) and np.array_equal(nearest,labels)
        for i,gene in enumerate(ids):
            shape_rows.append(dict(species=species[i],gene_id=gene,variant=variant,policy=policy,
                original_cluster='C'+str(labels[i]+1),evaluable=bool(ok[i]),
                fixed_center_cluster='C'+str(nearest[i]+1) if ok[i] else None,
                refitted_cluster='C'+str(refit[i]+1) if ok[i] else None,
                baseline_refit_on_same_eligible_set=int(eligible_baseline_refit[i]) if ok[i] else None,
                new_refit_matched_to_same_eligible_baseline=int(eligible_new_refit[i]) if ok[i] else None,
                feature_distance=float(np.linalg.norm(x[i]-oldx[i])) if ok[i] else np.nan))
        for sp in members:
            mask=species==sp;use=mask&ok
            shape_summary.append(dict(species=sp,variant=variant,policy=policy,fixed_genes=int(mask.sum()),
                evaluable=int(use.sum()),unavailable=int((mask&~ok).sum()),same_fixed_center=int((use&(nearest==labels)).sum()),
                same_refitted_cluster=int((use&(refit==labels)).sum()),
                same_refit_on_identical_eligible_set=int((use&(eligible_baseline_refit==eligible_new_refit)).sum()),
                original_features_vs_new_features_refit_ARI=float(adjusted_rand_score(eligible_baseline_refit[use],eligible_new_refit[use])) if use.sum()>1 else np.nan,
                refit_species_present=';'.join(sorted(set(species[ok]))),
                original_vs_refit_ARI=float(adjusted_rand_score(labels[use],refit[use])) if use.sum()>1 else np.nan,
                p90_feature_distance=float(np.quantile(np.linalg.norm(x[use]-oldx[use],axis=1),.9)) if use.any() else np.nan))
        print('SHAPE_COMPLETE',variant,policy,'available',int(ok.sum()),flush=True)
save(pd.DataFrame(shape_rows),'gene_shape_assignments.tsv.gz');save(pd.DataFrame(shape_summary),'shape_retention.tsv')
contrasts=['SR','RE','EP','SP','SM','MP'];effect_rows=[];effect_summary=[];profile_rows=[]
for sp,m in members.items():
    gs=m.gene_id;cluster=m.set_index('gene_id').cluster
    local=nodes[nodes.species.eq(sp)]
    for route in sorted(local.route.unique()):
        z=local[local.route.eq(route)]
        tables={(v,p):q.pivot(index='gene_id',columns='node',values='position_index').reindex(index=gs,columns=list('SREPM'))
            for (v,p),q in z.groupby(['variant','policy'])}
        base=tables['baseline','reapply_purity']
        for (variant,policy),new in tables.items():
            for contrast in contrasts:
                a=base[contrast[1]]-base[contrast[0]];b=new[contrast[1]]-new[contrast[0]]
                pa=a.notna()&b.notna();sa=sign(a);sb=sign(b)
                outcome=np.select([a.isna(),b.isna(),(sa!=0)&(sa==sb),(sa*sb==-1),(sa!=0)&(sb==0),(sa==0)&(sb!=0)],
                    ['baseline_unavailable','new_unavailable','same_direction','opposite_direction','weakened','new_effect'],default='both_weak')
                q=pd.DataFrame(dict(species=sp,route=route,variant=variant,policy=policy,contrast=contrast,
                    gene_id=gs.to_numpy(),original_cluster=gs.map(cluster).to_numpy(),baseline_delta=a.to_numpy(),new_delta=b.to_numpy(),
                    baseline_direction=sa,new_direction=sb,outcome=outcome,abs_difference=(b-a).abs().to_numpy()))
                effect_rows.append(q)
                selected=a.notna()&(sa!=0);available=pa&selected
                rho=float(spearmanr(a[pa],b[pa]).statistic) if pa.sum()>2 and a[pa].nunique()>1 and b[pa].nunique()>1 else np.nan
                effect_summary.append(dict(species=sp,route=route,variant=variant,policy=policy,contrast=contrast,
                    fixed_genes=len(gs),baseline_available=int(a.notna().sum()),paired=int(pa.sum()),
                    baseline_effects_ge001=int(selected.sum()),same_direction=int((selected&(outcome=='same_direction')).sum()),
                    opposite_direction=int((selected&(outcome=='opposite_direction')).sum()),weakened=int((selected&(outcome=='weakened')).sum()),
                    new_unavailable=int((selected&b.isna()).sum()),spearman=rho,
                    median_abs_difference=float((b[pa]-a[pa]).abs().median()),p90_abs_difference=float((b[pa]-a[pa]).abs().quantile(.9))))
            # Paired original-gene profiles: each displayed curve uses exactly the same genes at all four nodes.
            if route in routes[sp]:
                eligible=base[list('SREP')].notna().all(1)&new[list('SREP')].notna().all(1)
                for group in ['C1','C2','C3','C4']:
                    take=eligible & (gs.map(cluster).to_numpy()==group)
                    for estimate,table in [('baseline',base),('variant',new)]:
                        for node in 'SREP':profile_rows.append(dict(species=sp,route=route,variant=variant,policy=policy,
                            cluster=group,estimate=estimate,node=node,paired_genes=int(take.sum()),
                            fixed_cluster_genes=int((gs.map(cluster)==group).sum()),mean_position=table.loc[take,node].mean()))
save(pd.concat(effect_rows,ignore_index=True),'gene_local_effects.tsv.gz');save(pd.DataFrame(effect_summary),'local_effect_summary.tsv')
save(pd.DataFrame(profile_rows),'paired_fixed_group_profiles.tsv')
# Preserve the original local-comparison eligibility and independent-study denominator.
H={'GSE118258':['GSE118258'],'GSE150311':['GSE150311'],'GSE178325':['178325_0618'],'GSE221739':['GSE221739']}
B={'GSE115943':['WOT'],'GSE103221':['CIP','OSK']};C={'GSE115943':['WOT'],'GSE103221':['CIP']}
COMPARISONS={'human':{'SP':H,'SM':H,'MP':H},'mouse':{'SR':B,'RE':B,'EP':C,'SP':C,'MP':C,'SM':{**B,'GSE171127':['171127_GETM','171127_OSKM']}}}
effect=pd.concat(effect_rows,ignore_index=True).set_index(['species','gene_id','route','contrast','variant','policy'])
majority=[]
for sp,m in members.items():
    for gene in m.gene_id:
        for contrast,studymap in COMPARISONS[sp].items():
            need=len(studymap)//2+1;calls={}
            for variant in VARIANTS:
                for policy in ['reapply_purity','frozen_purity']:
                    votes=[];available=0
                    for study,rs in studymap.items():
                        vals=np.array([effect.loc[(sp,gene,r,contrast,variant,policy),'new_delta'] for r in rs])
                        good=np.isfinite(vals).all();available+=good
                        votes.append(1 if good and (vals>=.01).all() else -1 if good and (vals<=-.01).all() else 0)
                    direction=1 if votes.count(1)>=need else -1 if votes.count(-1)>=need else 0
                    calls[variant,policy]=(direction,available)
            od,oa=calls['baseline','reapply_purity']
            for (variant,policy),(nd,na) in calls.items():
                status='same_majority' if od and od==nd else 'opposite_majority' if od and od==-nd else 'insufficient_coverage' if od and na<need else 'no_majority' if od else 'new_majority' if nd else 'neither_selected'
                majority.append(dict(species=sp,gene_id=gene,contrast=contrast,variant=variant,policy=policy,
                    fixed_studies=len(studymap),required_studies=need,baseline_direction=od,new_direction=nd,
                    baseline_available_studies=oa,new_available_studies=na,outcome=status))
majority=pd.DataFrame(majority);save(majority,'gene_majority_outcomes.tsv.gz')
save(majority.groupby(['species','variant','policy','contrast','outcome']).size().rename('genes').reset_index(),'majority_summary.tsv')
assert all(sha(p)==h for p,h in receipts.items())
(D/'apa_analysis_verification.json').write_text(json.dumps(dict(status='PASS',fixed_gene_counts={sp:len(m) for sp,m in members.items()},
    baseline_joint_features_max_abs_difference=float(np.max(np.abs(baseline_features-oldx))),baseline_refit_exact=True,
    source_sha256=receipts,sklearn_version=sklearn_version,python_executable=sys.executable,
    comparison_registry=COMPARISONS,shape_rows=len(shape_rows),local_effect_rows=sum(len(q) for q in effect_rows),
    majority_rows=len(majority),all_sources_unchanged=True),indent=2)+'\n')
print('APA_ANALYSIS_COMPLETE',flush=True)
