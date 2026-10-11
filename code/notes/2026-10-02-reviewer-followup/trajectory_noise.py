"""Empirical state measurement error and full-cohort selection-aware monotone null."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from apa_followup import HERE, M, H, SPECS, read, save, rng_for, INPUTS

ROUTES = {'GSE103221': ['CIP'], 'GSE115943': ['WOT'], 'GSE118258': ['GSE118258'],
          'GSE178325': ['178325_0618'], 'GSE150311': ['GSE150311']}
STATES = {'CIP':'SREP','WOT':'SREP','GSE118258':'SREP','178325_0618':'SRP','GSE150311':'SEP'}
OUT = HERE / 'trajectory-noise'
B = 500


def extract(study):
    species, folder = SPECS[study]
    base = (M / 'data/metacell-pas-v7-v1/studies' if species == 'mouse' else H / 'data/results_metacell_pas_v2/studies') / folder
    support = (M / 'data/pas-annotation-support-v1/studies' if species == 'mouse' else M / 'data/apa-majority-shared-v1/human_annotation') / folder
    f = read(base / 'features.tsv.gz', low_memory=False)
    f = f[f.keep_for_apa_primary & f.primary_gene_id.notna()].copy()
    shape = f.groupby('primary_gene_id').agg(n=('cut_site','nunique'),c=('chrom','nunique'),s=('strand','nunique'))
    f = f[f.primary_gene_id.isin(shape.query('n>=2 and c==1 and s==1').index)]
    coord = f.cut_site * np.where(f.strand.eq('+'),1,-1)
    f['position'] = (coord - coord.groupby(f.primary_gene_id).transform('min')) / (coord.groupby(f.primary_gene_id).transform('max')-coord.groupby(f.primary_gene_id).transform('min'))
    f = f.merge(read(support / 'primary_pas_reference_matches.tsv.gz',usecols=['atlas_pas_id','matched_source_count_50nt']),on='atlas_pas_id',validate='one_to_one')
    f = f[f.matched_source_count_50nt.ge(1)].reset_index(drop=True)
    meta = read(M / 'data/apa-majority-shared-v1/profiles' / folder / 'selected_metacells.tsv.gz')
    libs = read(M / 'data/apa-majority-shared-v1/profiles' / folder / 'supported50_library_positions.tsv.gz')
    libs = libs[libs.route.isin(ROUTES[study]) & libs.node.isin(list('SREP'))]
    wanted = set(libs.gene_id)
    f = f[f.primary_gene_id.isin(wanted)].reset_index(drop=True)
    genes = pd.Index(sorted(wanted))
    gi = genes.get_indexer(f.primary_gene_id)
    c = sparse.load_npz(base / 'metacell_by_pas_counts.npz').tocsr()[meta.matrix_row.to_numpy()][:,f.feature_index_1based.to_numpy()-1]
    mapping = sparse.csr_matrix((np.ones(len(f)),(np.arange(len(f)),gi)),shape=(len(f),len(genes)))
    weighted = sparse.csr_matrix((f.position,(np.arange(len(f)),gi)),shape=mapping.shape)
    totals = (c @ mapping).toarray()
    numerators = (c @ weighted).toarray()
    x = np.divide(numerators,totals,out=np.full(totals.shape,np.nan),where=totals>=20)
    pieces, audit = [], []
    for route in ROUTES[study]:
        states = list(STATES[route])
        raw = np.full((len(genes),len(states)),np.nan)
        boot = np.full((B,len(genes),len(states)),np.nan)
        for j,node in enumerate(states):
            sums = np.zeros(len(genes)); nlib = np.zeros(len(genes),int); draws = np.zeros((B,len(genes)))
            for lib, block in libs[libs.route.eq(route) & libs.node.eq(node)].groupby('library_key'):
                rows = np.flatnonzero(meta.route.eq(route) & meta.library_key.eq(lib) & meta.node.eq(node) & meta.purity.ge(.8))
                for rec in block.itertuples():
                    ix = genes.get_loc(rec.gene_id)
                    values = x[rows,ix];values = values[np.isfinite(values)]
                    assert len(values)==rec.qualified_metacells
                    mean = values.mean()
                    assert abs(mean-rec.position_index)<1e-8,(study,lib,node,rec.gene_id)
                    rng = rng_for('trajectory',study,route,lib,node,rec.gene_id)
                    v = values[rng.integers(len(values),size=(B,len(values)))].mean(1)
                    sums[ix] += mean;draws[:,ix] += v;nlib[ix] += 1
            raw[:,j] = np.divide(sums,nlib,out=np.full(len(genes),np.nan),where=nlib>0)
            boot[:,:,j] = np.divide(draws,nlib[None],out=np.full(draws.shape,np.nan),where=nlib[None]>0)
        complete = np.isfinite(raw).all(1)
        gid = genes[complete].to_numpy(str); raw = raw[complete]; boot = boot[:,complete]
        fraction = np.linspace(0,1,len(states))
        monotone = raw[:,[0]] + (raw[:,[-1]]-raw[:,[0]])*fraction
        simulated = np.clip(monotone[None] + (boot-raw[None]),0,1)
        OUT.mkdir(exist_ok=True,parents=True)
        np.savez_compressed(OUT / f'{route}.npz',gene_ids=gid,observed=raw,bootstrap=boot,null=simulated,states=np.array(states))
        audit.append(dict(study_id=study,route=route,complete_background=len(gid),library_state_recovery='PASS'))
        print('TRAJECTORY',route,len(gid),flush=True)
    (OUT / f'{study}-audit.json').write_text(json.dumps(dict(records=audit,input_sha256=INPUTS),indent=2)+'\n')


def transient(x):
    length = abs(np.diff(x,axis=-1)).sum(-1)
    displacement = abs(x[...,-1]-x[...,0])
    return np.divide(length-displacement,length,out=np.zeros_like(length),where=length>1e-12)


def selection(a,b,ia,ib):
    u,v=a[:,ia],b[:,ib]
    u=u-u.mean(1,keepdims=True);v=v-v.mean(1,keepdims=True)
    norm=np.linalg.norm(u,axis=1)*np.linalg.norm(v,axis=1)
    corr=np.divide((u*v).sum(1),norm,out=np.zeros(len(u)),where=norm>1e-12)
    return (corr>=.6)&(np.ptp(a[:,ia],axis=1)>=.01)&(np.ptp(b[:,ib],axis=1)>=.01)


def choose(data):
    a,b=data['WOT'],data['CIP'];g=sorted(set(a[0])&set(b[0]))
    aa=a[1][pd.Index(a[0]).get_indexer(g)];bb=b[1][pd.Index(b[0]).get_indexer(g)]
    mouse=np.array(g)[selection(aa,bb,np.arange(4),np.arange(4))]
    human=set();x=data['GSE118258']
    for route,indices in [('178325_0618',[0,1,3]),('GSE150311',[0,2,3])]:
        y=data[route];g=sorted(set(x[0])&set(y[0]));xx=x[1][pd.Index(x[0]).get_indexer(g)];yy=y[1][pd.Index(y[0]).get_indexer(g)]
        use=selection(xx,yy,np.array(indices),np.arange(3))&(np.ptp(xx,axis=1)>=.01)
        human.update(np.array(g)[use])
    return {'mouse':mouse,'human':np.array(sorted(human))}


def summarize():
    paths={r:OUT / f'{r}.npz' for r in STATES}
    banks={r:dict(np.load(p,allow_pickle=False)) for r,p in paths.items()}
    members=read(M / 'outputs/figure2-k3-srelative-v2/gene_membership_k3.tsv.gz')
    original={r:(d['gene_ids'],d['observed']) for r,d in banks.items()}
    selected=choose(original)
    for species,ids in selected.items():
        assert set(ids)==set(members[members.species.eq(species)].gene_id),(species,len(ids),'candidate replay mismatch')
    fixed=[]
    for route,species in [('WOT','mouse'),('CIP','mouse'),('GSE118258','human')]:
        d=banks[route];m=members[members.species.eq(species)].copy();ix=pd.Index(d['gene_ids']).get_indexer(m.gene_id)
        assert (ix>=0).all()
        obs=transient(d['observed'][ix]);null=transient(d['null'][:,ix])
        for cluster in ['C1','C2','C3']:
            mask=m.cluster.eq(cluster).to_numpy();observed=float(np.median(obs[mask]));draws=np.median(null[:,mask],axis=1)
            fixed.append(dict(species=species,route=route,cluster=cluster,n=int(mask.sum()),observed_median_T=observed,
                null_median_T=float(np.median(draws)),null_lo=float(np.quantile(draws,.025)),null_hi=float(np.quantile(draws,.975)),
                observed_minus_null_median=observed-float(np.median(draws)),scope='fixed_members_and_labels_conditional'))
    save(pd.DataFrame(fixed),OUT / 'fixed_class_T.tsv')
    rows=[]
    for rep in range(B):
        data={r:(d['gene_ids'],d['null'][rep]) for r,d in banks.items()}
        cohort=choose(data);features=[];species=[];geneids=[];rawroutes={}
        for sp,routes in [('mouse',['WOT','CIP']),('human',['GSE118258'])]:
            ids=cohort[sp];views=[]
            for route in routes:
                gids,x=data[route];raw=x[pd.Index(gids).get_indexer(ids)];rawroutes[sp,route]=raw
                z=raw-raw.mean(1,keepdims=True);views.append(z/np.linalg.norm(z,axis=1,keepdims=True))
            features.append(np.mean(views,axis=0));species.extend([sp]*len(ids));geneids.extend(ids)
        if min(len(cohort['mouse']),len(cohort['human']))<3:
            rows.append(dict(replicate=rep,status='insufficient_selected_genes'));continue
        x=np.concatenate(features);sp=np.array(species);weights=np.array([.5/np.sum(sp==s) for s in sp])
        fit=KMeans(n_clusters=3,n_init=50,tol=1e-7,max_iter=500,random_state=91327).fit(x,sample_weight=weights)
        centers=fit.cluster_centers_;order=np.argsort(centers[:,3]-centers[:,0]);label={k:f'C{j+1}' for j,k in enumerate(order)}
        labs=np.array([label[k] for k in fit.labels_])
        for (s,r),raw in rawroutes.items():
            for cluster in ['C1','C2','C3']:
                mask=labs[sp==s]==cluster
                rows.append(dict(replicate=rep,species=s,route=r,cluster=cluster,n=int(mask.sum()),selected_species_genes=len(cohort[s]),
                                 median_T=float(np.median(transient(raw[mask]))) if mask.any() else np.nan,status='full_background_reselection_reclustering'))
        if rep%50==0:print('RESELECT',rep,{k:len(v) for k,v in cohort.items()},flush=True)
    frame=pd.DataFrame(rows);save(frame,OUT / 'reselected_draws.tsv.gz')
    summary=[]
    for key,f in frame[frame.status.eq('full_background_reselection_reclustering')].groupby(['species','route','cluster']):
        d=f.median_T.dropna();summary.append(dict(zip(['species','route','cluster'],key)) | dict(valid_draws=len(d),median_T=float(d.median()),T_lo=float(d.quantile(.025)),T_hi=float(d.quantile(.975)),selected_n_median=float(f.n.median()),selected_species_n_median=float(f.selected_species_genes.median())))
    save(pd.DataFrame(summary),OUT / 'reselected_summary.tsv')
    print(pd.DataFrame(fixed).to_string(index=False));print(pd.DataFrame(summary).to_string(index=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['extract','summarize']);p.add_argument('--study',choices=list(ROUTES));args=p.parse_args()
    extract(args.study) if args.command=='extract' else summarize()
