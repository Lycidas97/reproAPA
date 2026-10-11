import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import subprocess
import sys
from common import *


def decode(a): return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])


def run():
    import anndata as ad
    import scanpy as sc
    x=obs();out=RUN/'A03';shared=RUN/'shared'
    libraries=pd.read_csv(shared/'raw_count_manifest.tsv',sep='\t')
    orthologs=JOINT/'sources/local/ensembl-116-human-mouse-one-to-one-orthologs.tsv'
    f=pd.read_csv(orthologs,sep='\t',dtype=str).fillna('')
    f=f[f.documented_one_to_one.str.lower().eq('true')].reset_index(drop=True)
    assert f.human_gene_id.is_unique and f.mouse_gene_id.is_unique
    present=np.ones(len(f),bool);headers={};input_receipts=[]
    for row in libraries.itertuples():
        with h5py.File(row.input_path,'r') as h:
            g=decode(h['matrix/features/id'][:]);t=decode(h['matrix/features/feature_type'][:])
        assert pd.Index(g).is_unique
        genes=g[t=='Gene Expression']
        species=x.loc[x.raw_library_id.eq(row.library_id),'species_key'].unique();assert len(species)==1
        present &= f[species[0]+'_gene_id'].isin(genes).to_numpy()
        headers[row.library_id]=(g,t,species[0])
    f=f[present].reset_index(drop=True)
    assert len(f)>2128
    f.to_csv(shared/'ortholog_features.tsv',sep='\t',index=False)
    nparts=[];orders=[];total=np.zeros(len(x),np.float64)
    for number,row in enumerate(libraries.sort_values('sample_accession').itertuples(),1):
        progress('A03','RUNNING',stage='extract_full_ortholog_counts',library=number,total_libraries=161,sample=row.sample_accession,features=len(f))
        source=Path(row.input_path);before=source.stat()
        digest=sha(source)
        with h5py.File(source,'r') as h:
            m=h['matrix'];b=decode(m['barcodes'][:])
            raw=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr()
        g,t,species=headers[row.library_id]
        selected=x.index[x.raw_library_id.eq(row.library_id)].to_numpy()
        ci=pd.Index(b).get_indexer(x.loc[selected,'barcode']);assert (ci>=0).all()
        fi=pd.Index(g).get_indexer(f[species+'_gene_id']);assert (fi>=0).all()
        whole=raw[ci]; assert np.all(whole.data>=0) and np.all(whole.data==np.rint(whole.data))
        total[selected]=np.asarray(whole[:,t=='Gene Expression'].sum(1)).ravel()
        part=whole[:,fi].astype(np.float32).tocsr()
        nparts.append(part);orders.append(selected)
        assert source.stat().st_mtime_ns==before.st_mtime_ns and source.stat().st_size==before.st_size
        input_receipts.append(dict(library=row.library_id,path=str(source),sha256=digest,bytes=before.st_size,
                                   mtime_ns=before.st_mtime_ns,cells=len(selected),features=len(f)))
        write_json(out/'raw_count_receipt.json',input_receipts)
        del raw,whole,part
    counts=sparse.vstack(nparts,format='csr');order=np.concatenate(orders);del nparts
    assert len(np.unique(order))==len(x) and np.array_equal(np.sort(order),np.arange(len(x)))
    counts=counts[np.argsort(order)]
    assert counts.shape==(len(x),len(f)) and (total>0).all()
    sparse.save_npz(shared/'ortholog_counts.npz',counts)
    np.save(shared/'full_gene_umi.npy',total)
    # Equal allocation across species, then studies; no state label or outcome selection.
    rng=np.random.default_rng(SEED);pilot=[]
    for species,d in x.groupby('species_key',sort=True):
        groups=list(d.groupby('true_study',sort=True).indices.keys())
        cap=30000//len(groups)
        for study in groups:
            idx=x.index[x.species_key.eq(species)&x.true_study.eq(study)].to_numpy()
            pilot.extend(rng.choice(idx,min(cap,len(idx)),replace=False))
    pilot=np.sort(np.asarray(pilot,int));np.save(shared/'pilot_indices.npy',pilot)
    write_json(shared/'counts_ready.json',dict(cells=len(x),features=len(f),pilot_cells=len(pilot),integer_counts=True,
        shared_counts_sha256=sha(shared/'ortholog_counts.npz'),full_gene_UMI_normalization=True))
    progress('A03','RUNNING',stage='freeze_holdout_and_symmetric_features')
    original=pd.read_csv(ROOT/'inputs/S6-integration-orthologs.tsv',sep='\t')
    excluded=set(original.human_gene_id)|set(original.mouse_gene_id)
    human_hvg=HUMAN/'results_v1/model/hvg_features.tsv'
    hd=pd.read_csv(human_hvg,sep='\t');excluded.update(hd.astype(str).to_numpy().ravel())
    prior=pd.read_csv(ROOT/'prior_evidence/mouse_state_matching/inputs_snapshot/features.tsv',sep='\t')
    excluded.update(prior.loc[prior.panel.eq('original_hvg'),'gene_id'])
    scopes=json.loads((RUN/'runtime/scopes.json').read_text())
    exclusion_sources=[human_hvg]
    for unit,path in scopes.items():
        if unit=='wot':continue
        source=Path(path)/'embedding/selected_features.tsv'
        part=pd.read_csv(source,sep='\t');excluded.update(part.gene_id.astype(str));exclusion_sources.append(source)
    wot=Path('/root/nfsdata/multica/worktrees/apadetect/DOC-54-wot-readiness/analyses/wot-unsupervised-first-topology-audit-v1/unsupervised_feature_manifest.tsv')
    w=pd.read_csv(wot,sep='\t');excluded.update(w.loc[w.primary_HVG,'gene_id']);exclusion_sources.append(wot)
    marker=set()
    labels=pd.read_csv(ROOT/'inputs/S5-local-state-annotations.tsv',sep='\t').fillna('')
    for value in labels.marker_evidence: marker.update(v.strip().upper() for v in value.split(';') if v.strip())
    programs=pd.read_csv(ROOT/'inputs/S7-mouse-annotation-programs.tsv',sep='\t')
    for col in programs:
        if 'symbol' in col:marker.update(programs[col].dropna().str.upper())
    marker.update(['POU5F1','SOX2','KLF4','MYC','NANOG','LIN28A','LIN28B'])
    invalid=(f.human_gene_id.isin(excluded)|f.mouse_gene_id.isin(excluded)|
             f.human_gene_symbol.str.upper().isin(marker)|f.mouse_gene_symbol.str.upper().isin(marker)|
             f.human_gene_symbol.str.match(r'^(MT-|RPL|RPS)')|f.mouse_gene_symbol.str.match(r'^(mt-|Rp[ls])'))
    f['holdout_eligible']=~invalid
    f['holdout_hash']=[hashlib.sha256(('reproAPA-holdout-20260913|'+g).encode()).hexdigest() for g in f.human_gene_id]
    hold=f.loc[~invalid].sort_values('holdout_hash').head(2000).index.to_numpy()
    assert len(hold)>=200, 'Too few unused genes for meaningful feature-holdout validation.'
    f['heldout']=f.index.isin(hold)
    # These genes will never be used to train the new symmetric-feature model.
    normalized=(sparse.diags(10000/total)@counts).tocsr();normalized.data=np.log1p(normalized.data)
    rank_tables=[]
    for (species,study),idx in x.groupby(['species_key','true_study'],sort=True).indices.items():
        progress('A03','RUNNING',stage='study_balanced_HVG_ranks',species=species,study=study)
        a=ad.AnnData(normalized[idx]);a.var_names=f.human_gene_id
        hv=sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=None,inplace=False)
        scores=hv.dispersions_norm.replace([np.inf,-np.inf],np.nan).fillna(-np.inf)
        rank=scores.rank(method='average',ascending=False,pct=True).to_numpy()
        rank_tables.append(pd.DataFrame(dict(feature_index=np.arange(len(f)),species=species,study=study,percentile_rank=rank)))
    ranks=pd.concat(rank_tables,ignore_index=True)
    means=ranks.groupby(['feature_index','species']).percentile_rank.mean().unstack('species')
    f['symmetric_rank']=means[['human','mouse']].mean(axis=1).to_numpy()
    selected=f.loc[~f.heldout].sort_values(['symmetric_rank','human_gene_id']).head(2128).index.to_numpy()
    assert len(selected)==2128 and not np.intersect1d(selected,hold).size
    original_idx=pd.Index(f.human_gene_id).get_indexer(original.human_gene_id)
    assert (original_idx>=0).all() and not np.intersect1d(original_idx,hold).size
    np.savez(shared/'feature_sets.npz',original=original_idx,symmetric=selected,holdout=hold)
    f.to_csv(shared/'feature_selection.tsv',sep='\t',index=False)
    ranks.to_csv(out/'study_HVG_ranks.tsv.gz',sep='\t',index=False)
    receipt([orthologs,*exclusion_sources,shared/'feature_sets.npz'],out/'feature_source_receipt.json')
    write_json(shared/'features_ready.json',dict(features=len(f),heldout=len(hold),symmetric=2128,original=2128,
        heldout_selection='hash ordered after exclusions; expression detection filtered separately within training folds'))
    del counts,normalized
    # A separate process imports torch only after a device is available.
    subprocess.run([sys.executable,str(RUN/'runtime/a03_train_queue.py')],check=True)


if __name__=='__main__': execute('A03',run)
