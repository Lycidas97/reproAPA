"""Re-select and refit trajectories; preserve the frozen expression scaffold."""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_samples

from downstream_parallel import ordered_map
from downstream_runtime import Context, Study, STUDIES, seed

NOISE='notes/2026-10-02-reviewer-followup/trajectory_noise.py'
ROUTES={'CIP':'GSE103221','WOT':'GSE115943','GSE118258':'GSE118258','178325_0618':'GSE178325','GSE150311':'GSE150311'}
STATES={'CIP':'SREP','WOT':'SREP','GSE118258':'SREP','178325_0618':'SRP','GSE150311':'SEP'}
FIT={'mouse':['WOT','CIP'],'human':['GSE118258']}


def transform(raw,method='centered'):
    z=raw-raw.mean(1,keepdims=True) if method=='centered' else raw-raw[:,[0]]
    norm=np.linalg.norm(z,axis=1,keepdims=True)
    assert np.isfinite(z).all() and (norm>1e-10).all()
    return z/norm


def feature(data,members,method='centered',omit=None):
    x=[];species=[];genes=[]
    for sp,routes in FIT.items():
        views=[];ids=members[sp]
        for route in routes:
            if route==omit:continue
            g,v=data[route];ix=pd.Index(g).get_indexer(ids)
            assert (ix>=0).all()
            views.append(transform(v[ix],method))
        assert views
        x.append(np.mean(views,axis=0));species.extend([sp]*len(ids));genes.extend(ids)
    return np.concatenate(x),np.array(species),np.array(genes)


def fit(x,species,k=3,random_seed=91327):
    weights=np.array([.5/np.sum(species==s) for s in species])
    assert set(species)=={'mouse','human'}
    model=KMeans(n_clusters=k,n_init=50,max_iter=500,tol=1e-7,random_state=random_seed).fit(x,sample_weight=weights)
    order=np.argsort(model.cluster_centers_[:,-1]-model.cluster_centers_[:,0])
    mapping=np.empty(k,int);mapping[order]=np.arange(k)
    return model,mapping[model.labels_],model.cluster_centers_[order],weights


def select_data(nodes,setting='supported50_p08',depth=20):
    data={}
    for route,states in STATES.items():
        q=nodes[nodes.route.eq(route)&nodes.setting.eq(setting)&nodes.depth.eq(depth)]
        f=q.pivot(index='gene_id',columns='node',values='position_index').reindex(columns=list(states)).dropna()
        data[route]=(f.index.to_numpy(str),f.to_numpy())
    return data


def sensitivity_members(data,members):
    """A fixed candidate can lose its classifiable shape at another depth/layer."""
    chosen={};audit=[]
    for sp,routes in FIT.items():
        ids=members[sp];keep=np.ones(len(ids),bool)
        for route in routes:
            genes,values=data[route];indices=pd.Index(genes).get_indexer(ids)
            for j,(gene,index) in enumerate(zip(ids,indices)):
                status='missing_complete_profile'
                if index>=0:
                    v=values[index];norm=np.linalg.norm(v-v.mean())
                    status='classifiable' if np.isfinite(v).all() and norm>1e-10 else 'constant_or_nonfinite_profile'
                keep[j]&=status=='classifiable'
                audit.append(dict(species=sp,route=route,gene_id=gene,status=status))
        chosen[sp]=ids[keep]
    return chosen,pd.DataFrame(audit)


def reconstruct_profile_setting(ctx,item):
    sid,setting=item
    libs=[];mc_audit=[]
    study=Study(ctx,sid);study.load_gex()
    current=study.core('library_states',setting)
    if current.empty:return (libs,mc_audit)
    for gene,group in current.groupby('gene_id',sort=False):
        axis=study.axis(gene,setting);data=study.gene_counts(axis)
        for depth in [20,40,80]:
            eligible=study.profile_states(data,setting,depth)
            for row in group.to_dict('records'):
                s=eligible.get((row['route'],row['library_key'],row['node']))
                if s is None:continue
                mean=float(s['usage']@axis.position)
                if depth==20:assert abs(mean-row['position_index'])<1e-12
                expr=study.expression_summary(gene,s['rows'])
                libs.append(dict(species=row['species'],study_id=sid,route=row['route'],setting=setting,
                    gene_id=gene,gene_symbol=row['gene_symbol'],library_key=row['library_key'],node=row['node'],depth=depth,
                    position_index=mean,GEX_CPM=expr['GEX_CPM'],gex_log2_cpm1=expr['GEX'],qualified_metacells=s['n_mc'],
                    qualified_cells=s['n_cells'],PAS_number=len(axis),mean_pas_umi=float(s['counts'].sum(1).mean()),
                    within_library_apa_range=float(np.ptp(s['counts']@axis.position/s['counts'].sum(1))),
                    metacell_ids=json.dumps(study.meta.metacell_id.iloc[s['rows']].tolist())))
    mc_audit.append(dict(study_id=sid,setting=setting,core_library_gene_states=len(current),
        core_genes=current.gene_id.nunique(),status='recovered_at_depth20'))
    print('profiles',sid,flush=True)
    return (libs,mc_audit)


def profiles(ctx):
    libs=[];mc_audit=[]
    for values,audit in ordered_map(ctx,reconstruct_profile_setting,[(sid,setting) for sid in STUDIES for setting in ['supported50_p08','expanded_primary_p08']]):
        libs.extend(values);mc_audit.extend(audit)
    libraries=pd.DataFrame(libs)
    nodes=libraries.groupby(['species','study_id','route','setting','depth','gene_id','gene_symbol','node'],dropna=False).agg(
        position_index=('position_index','mean'),GEX_CPM=('GEX_CPM','mean'),gex_log2_cpm1=('gex_log2_cpm1','mean'),n_libraries=('library_key','nunique'),qualified_metacells=('qualified_metacells','sum'),
        mean_pas_umi=('mean_pas_umi','mean'),within_library_apa_range=('within_library_apa_range','max')).reset_index()
    ctx.table(libraries,'all_library_positions.tsv.gz');ctx.table(nodes,'all_route_node_positions.tsv.gz')
    math=ctx.functions(NOISE,['selection','choose','transient'])
    data=select_data(nodes);members=math.choose(data)
    assert min(map(len,members.values()))>=3,'Too few newly reproducible profiles to fit the planned joint taxonomy'
    x,species,ids=feature(data,members)
    all_labels=[];centers=[];metrics=[];archive=dict(features=x,species=species,gene_ids=ids,
        weights=np.array([.5/np.sum(species==s) for s in species]))
    for k in range(2,7):
        assert len(x)>k
        model,labels,centroid,weights=fit(x,species,k)
        all_labels.append(pd.DataFrame(dict(k=k,species=species,gene_id=ids,cluster=[f'C{j+1}' for j in labels])))
        centers.extend(dict(k=k,cluster=f'C{i+1}',**dict(zip('SREP',c))) for i,c in enumerate(centroid))
        sil=silhouette_samples(x,labels)
        metrics.append(dict(k=k,species_balanced_silhouette=np.average(sil,weights=weights),inertia=model.inertia_,
            minimum_mouse_size=np.bincount(labels[species=='mouse'],minlength=k).min(),
            minimum_human_size=np.bincount(labels[species=='human'],minlength=k).min()))
        archive[f'labels_k{k}']=labels
    membership=pd.concat(all_labels);primary=membership[membership.k.eq(3)].drop(columns='k')
    ctx.table(primary,'gene_membership_k3.tsv.gz');ctx.table(membership,'all_resolution_memberships.tsv.gz')
    ctx.table(centers,'shape_centers.tsv');ctx.table(metrics,'resolution_comparison.tsv')
    np.savez_compressed(ctx.out/'joint_fit.npz',**archive)
    geometry=[]
    for route,(genes,raw) in data.items():
        sp='mouse' if route in ['CIP','WOT'] else 'human';m=primary[primary.species.eq(sp)]
        for rec in m[m.gene_id.isin(genes)].itertuples():
            values=raw[pd.Index(genes).get_loc(rec.gene_id)];length=abs(np.diff(values)).sum();displacement=abs(values[-1]-values[0])
            geometry.append(dict(species=sp,route=route,gene_id=rec.gene_id,cluster=rec.cluster,L=length,D=displacement,
                T=(length-displacement)/length if length>1e-12 else 0.,endpoint_delta=values[-1]-values[0]))
    ctx.table(geometry,'trajectory_geometry.tsv.gz');ctx.table(mc_audit,'state_reconstruction.tsv')
    architecture(ctx,nodes,primary)
    ctx.finish(dict(full_background_candidate_selection=True,new_k3_fit=True,
        species_balanced=True,expanded_and_20_40_80_depth_profiles=True,matched_GEX_contributors=True),
        selected_genes={s:len(g) for s,g in members.items()},seed=91327)


def architecture(ctx,nodes,members):
    source='github/reproapa-analysis/modules/07_apa_dynamics/code/figure2/analyze_figure2_dynamic_architecture_v1.py'
    math=ctx.functions(source,['unit_centered','summarize','trajectory_architecture','loo_route_validation'],
        dict(STATES=list('SREP'),TRANSITIONS=[('S','R'),('R','E'),('E','P')],ROUTES=FIT,
            CLUSTERS=['C1','C2','C3'],SEED=91327,FIT=ctx.out/'joint_fit.npz',KMeans=KMeans,
            linear_sum_assignment=linear_sum_assignment,adjusted_rand_score=adjusted_rand_score))
    n=nodes[nodes.depth.eq(20)].copy();n['layer']=n.setting.map({'supported50_p08':'supported50','expanded_primary_p08':'expanded_primary'})
    gene,route,species,cluster=math.trajectory_architecture(members,n)
    for f,name in [(gene,'architecture_gene_metrics.tsv.gz'),(route,'architecture_route_summary.tsv'),
        (species,'architecture_species_summary.tsv'),(cluster,'architecture_cluster_summary.tsv')]:ctx.table(f,name)
    loo,centroid,support=math.loo_route_validation(members,n)
    ctx.table(loo,'leave_route_projection.tsv');ctx.table(centroid,'leave_route_centroids.tsv');ctx.table(support,'route_shape_support.tsv')
    sensitivity=[]
    for measure in ['path_length_L','SREP_range']:
        q=gene.copy();q['SREP_range']=q[['delta_'+s for s in 'SREP']].max(axis=1)-q[['delta_'+s for s in 'SREP']].min(axis=1)
        for cutoff in np.round(np.arange(0,.201,.02),2):
            for key,g in q.groupby(['species','route','cluster']):
                retained=g[g[measure].ge(cutoff)]
                sensitivity.append(dict(zip(['species','route','cluster'],key))|dict(measure=measure,cutoff=cutoff,
                    original_n=len(g),retained_n=len(retained),median_T=retained.transient_fraction_T.median()))
    ctx.table(sensitivity,'trajectory_amplitude_sensitivity.tsv')


def profile_gex(ctx):
    nodes=ctx.read(ctx.analysis/'profiles/all_route_node_positions.tsv.gz')
    members=ctx.read(ctx.analysis/'profiles/gene_membership_k3.tsv.gz')
    q=nodes[nodes.setting.eq('supported50_p08')&nodes.depth.eq(20)].merge(members,on=['species','gene_id'],validate='many_to_one')
    rows=[];transitions=[]
    for key,f in q.groupby(['species','study_id','route','gene_id','cluster']):
        states=[s for s in 'SREP' if s in set(f.node)]
        if len(states)<3:continue
        d=f.set_index('node').loc[states];a=d.position_index.to_numpy();g=d.gex_log2_cpm1.to_numpy()
        finite=np.isfinite(g).all();corr=np.corrcoef(a,g)[0,1] if finite and np.ptp(a)>0 and np.ptp(g)>0 else np.nan
        base=dict(zip(['species','study_id','route','gene_id','cluster'],key))
        rank_corr=pd.Series(a).corr(pd.Series(g),method='spearman') if finite and np.ptp(a)>0 and np.ptp(g)>0 else np.nan
        rows.append(base|dict(n_states=len(states),APA_GEX_r=corr,APA_GEX_r_spearman=rank_corr,missing_GEX=not finite))
        for i in range(len(states)-1):
            da=a[i+1]-a[i];dg=g[i+1]-g[i]
            transitions.append(base|dict(transition=states[i]+'_'+states[i+1],APA_delta=da,GEX_log_delta=dg,
                same_direction=bool(da*dg>0),zero_change=bool(da==0 or dg==0),GEX_observed=bool(np.isfinite(dg))))
    ctx.table(rows,'APA_GEX_gene_correlations.tsv.gz');correlation_summary(ctx,pd.DataFrame(rows));tr=pd.DataFrame(transitions);ctx.table(tr,'adjacent_transitions.tsv.gz')
    ctx.table(tr[tr.GEX_observed].groupby(['species','route','cluster','transition']).agg(
        genes=('gene_id','size'),same_direction=('same_direction','sum'),zero_changes=('zero_change','sum')).reset_index(),
        'adjacent_direction_summary.tsv')
    from downstream_profile_functional import functional_ranges
    functional_ranges(ctx,nodes,members)
    ctx.finish(dict(full_SREP_ranges_and_functional_families=True,same_PAS_qualified_MC_GEX=True,observed_states_only=True,zero_changes_in_denominator=True))


def correlation_summary(ctx,frame):
    rows=[]
    for keys in [['species','cluster'],['species','cluster','n_states'],['species','route','cluster']]:
        for key,g in frame.groupby(keys):
            for method,column in [('Pearson','APA_GEX_r'),('Spearman','APA_GEX_r_spearman')]:
                x=g[column].dropna()
                rows.append(dict(zip(keys,key))|dict(grouping=';'.join(keys),method=method,n_total=len(g),n_estimable=len(x),
                    median_r=x.median(),q25_r=x.quantile(.25),q75_r=x.quantile(.75)))
    ctx.table(rows,'APA_GEX_correlation_summary.tsv')


def controls(ctx):
    nodes=ctx.read(ctx.analysis/'profiles/all_route_node_positions.tsv.gz')
    libraries=ctx.read(ctx.analysis/'profiles/all_library_positions.tsv.gz')
    member=ctx.read(ctx.analysis/'profiles/gene_membership_k3.tsv.gz')
    members={s:f.gene_id.to_numpy(str) for s,f in member.groupby('species')}
    data=select_data(nodes);x,species,ids=feature(data,members)
    baselines={m:{k:fit(feature(data,members,m)[0],species,k) for k in range(2,7)} for m in ['centered','S_zero']}
    models=baselines['centered'];rows=[];rng=np.random.default_rng(91328)
    def compare(kind,rep,alternate,selected=None,method='centered',reference=None):
        if selected is None:selected=np.arange(len(x))
        for k,(model,labels,centroid,w) in baselines[reference or method].items():
            if len(selected)<=k:
                rows.append(dict(kind=kind,replicate=rep,method=method,k=k,status='insufficient_genes',genes=len(selected)));continue
            fitted,ll,cc,ww=fit(alternate,species[selected],k)
            aa,bb=linear_sum_assignment(((centroid[:,None]-cc[None])**2).sum(2));mapping=dict(zip(bb,aa));aligned=np.array([mapping[v] for v in ll])
            overlaps=[((labels[selected]==j)&(aligned==j)).sum()/((labels[selected]==j)|(aligned==j)).sum()
                for j in range(k) if ((labels[selected]==j)|(aligned==j)).any()]
            rows.append(dict(kind=kind,replicate=rep,method=method,k=k,ARI=adjusted_rand_score(labels[selected],ll),
                assignment_agreement=float((labels[selected]==aligned).mean()),mean_Jaccard=np.mean(overlaps),minimum_Jaccard=min(overlaps),
                centroid_RMSE=float(np.sqrt(np.mean((centroid[aa]-cc[bb])**2))),reference_method=reference or method,genes=len(selected),status='refitted'))
    for method in ['centered','S_zero']:
        xx,_,_=feature(data,members,method)
        if method=='S_zero':compare('S_origin',0,xx,method=method,reference='centered')
        for rep in range(24):
            ix=np.sort(np.concatenate([rng.choice(np.flatnonzero(species==sp),int(np.ceil((species==sp).sum()*.8)),replace=False) for sp in FIT]))
            compare('80pct_gene_subsample',rep,xx[ix],ix,method)
        for rep in range(32):
            parts=[]
            for route in set(sum(FIT.values(),[])):
                f=libraries[libraries.route.eq(route)&libraries.setting.eq('supported50_p08')&libraries.depth.eq(20)].copy()
                lw={lib:rng.exponential() for lib in f.library_key.unique()};f['weight']=f.library_key.map(lw)
                f['weighted']=f.position_index*f.weight
                v=f.groupby(['gene_id','node']).agg(numerator=('weighted','sum'),denominator=('weight','sum')).reset_index()
                v['position_index']=v.numerator/v.denominator;v['route']=route;v['setting']='supported50_p08';v['depth']=20;parts.append(v)
            alt=select_data(pd.concat(parts))
            compare('library_reweight',rep,feature(alt,members,method)[0],method=method)
        for new_seed in [91427,42,20261003]:
            for k,(_,labels,_,_) in baselines[method].items():
                _,ll,_,_=fit(xx,species,k,new_seed)
                rows.append(dict(kind='seed',replicate=new_seed,method=method,k=k,ARI=adjusted_rand_score(labels,ll),genes=len(x),status='refitted'))
    availability=[]
    for setting,depth in [('expanded_primary_p08',20),('supported50_p08',40),('supported50_p08',80)]:
        alt=select_data(nodes,setting,depth)
        chosen,eligibility=sensitivity_members(alt,members)
        eligibility['setting']=setting;eligibility['depth']=depth;availability.append(eligibility)
        if min(map(len,chosen.values()))<1:
            rows.append(dict(kind=f'{setting}_depth{depth}',status='no_classifiable_two_species_cohort',genes=0));continue
        xx,ss,gg=feature(alt,chosen);lookup={(s,g):i for i,(s,g) in enumerate(zip(species,ids))}
        ix=np.array([lookup[s,g] for s,g in zip(ss,gg)])
        compare(f'{setting}_depth{depth}',0,xx,ix)
    for route in ['CIP','WOT']:
        compare('leave_route_'+route,0,feature(data,members,omit=route)[0])
    ctx.table(rows,'clustering_robustness.tsv')
    ctx.table(pd.concat(availability,ignore_index=True),'profile_sensitivity_eligibility.tsv')
    trajectory_noise(ctx,libraries,member)
    ctx.finish(dict(new_candidate_perturbations=True,centered_and_S_origin=True,
        full_background_monotone_selection=True,reclustering_each_draw=True),bootstrap_draws=500)


def trajectory_noise(ctx,libraries,member):
    math=ctx.functions(NOISE,['selection','choose','transient'])
    banks={};audits=[]
    for route,sid in ROUTES.items():
        study=Study(ctx,sid);states=list(STATES[route])
        libs=libraries[libraries.route.eq(route)&libraries.setting.eq('supported50_p08')&libraries.depth.eq(20)]
        genes=pd.Index(sorted(libs.gene_id.unique()));raw=np.full((len(genes),len(states)),np.nan)
        boot=np.full((500,len(genes),len(states)),np.nan)
        for gene,f in libs.groupby('gene_id',sort=False):
            i=genes.get_loc(gene);axis=study.axis(gene);data=study.gene_counts(axis)
            for j,node in enumerate(states):
                means=[];draws=[]
                for r in f[f.node.eq(node)].itertuples():
                    s=study.state(data,route,r.library_key,node)
                    vals=s['counts']@axis.position/s['counts'].sum(1)
                    assert len(vals)==r.qualified_metacells and abs(vals.mean()-r.position_index)<1e-12
                    random=np.random.default_rng(seed('trajectory',sid,route,r.library_key,node,gene))
                    means.append(vals.mean());draws.append(vals[random.integers(len(vals),size=(500,len(vals)))].mean(1))
                if means:raw[i,j]=np.mean(means);boot[:,i,j]=np.mean(draws,axis=0)
        ok=np.isfinite(raw).all(1);raw=raw[ok];boot=boot[:,ok];gid=genes[ok].to_numpy(str)
        monotone=raw[:,[0]]+(raw[:,[-1]]-raw[:,[0]])*np.linspace(0,1,len(states))
        null=np.clip(monotone[None]+boot-raw[None],0,1)
        banks[route]=dict(gene_ids=gid,observed=raw,bootstrap=boot,null=null,states=np.array(states))
        np.savez_compressed(ctx.out/f'trajectory_{route}.npz',**banks[route])
        audits.append(dict(route=route,complete_background=len(gid)))
    original={r:(d['gene_ids'],d['observed']) for r,d in banks.items()};selected=math.choose(original)
    for sp,g in selected.items():assert set(g)==set(member[member.species.eq(sp)].gene_id)
    rows=[];fixed=[]
    for route in ['WOT','CIP','GSE118258']:
        sp='mouse' if route in ['WOT','CIP'] else 'human';d=banks[route]
        m=member[member.species.eq(sp)];ix=pd.Index(d['gene_ids']).get_indexer(m.gene_id)
        for klass in ['C1','C2','C3']:
            use=m.cluster.eq(klass).to_numpy()
            if not use.any():continue
            observed=np.median(math.transient(d['observed'][ix[use]]));null=np.median(math.transient(d['null'][:,ix[use]]),axis=1)
            fixed.append(dict(species=sp,route=route,cluster=klass,genes=int(use.sum()),observed_median_T=observed,
                null_median_T=np.median(null),null_lo=np.quantile(null,.025),null_hi=np.quantile(null,.975)))
    for rep in range(500):
        data={r:(d['gene_ids'],d['null'][rep]) for r,d in banks.items()};cohort=math.choose(data)
        if min(map(len,cohort.values()))<3:
            rows.append(dict(replicate=rep,status='insufficient_selected_genes'));continue
        x,sp,ids=feature(data,cohort);_,labels,_,_=fit(x,sp)
        for species,routes in FIT.items():
            for route in routes:
                g,v=data[route];raw=v[pd.Index(g).get_indexer(cohort[species])]
                for j in range(3):
                    use=labels[sp==species]==j
                    rows.append(dict(replicate=rep,species=species,route=route,cluster=f'C{j+1}',n=int(use.sum()),
                        selected_species_genes=len(cohort[species]),median_T=np.median(math.transient(raw[use])) if use.any() else np.nan,
                        status='full_background_reselection_reclustering'))
        if rep%50==0:print('monotone refit',rep,flush=True)
    ctx.table(rows,'trajectory_reselected_draws.tsv.gz');ctx.table(fixed,'trajectory_fixed_class_T.tsv')
    ctx.table(audits,'trajectory_noise_background.tsv')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--node',choices=['profiles','profile_controls','profile_gex'],required=True)
    a=p.parse_args();{'profiles':profiles,'profile_controls':controls,'profile_gex':profile_gex}[a.node](Context(a.config,a.node,a.out))
