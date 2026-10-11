"""Read frozen partitions and QC; quantify structure without new identity labels."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2')
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data';OUT.mkdir(exist_ok=True)
BASE=Path('/root/nfsdata/apadetect/analysis_projects/reproapa_integration_validation_20260913')
RUN=BASE/'execution_v1'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(d,n):d.to_csv(OUT/n,sep='\t',index=False)
S=read(ROOT/'inputs/S5-local-state-annotations.tsv').set_index('original_state_id')
M=read(RUN/'A07/partition_metrics.tsv');receipts={str(RUN/'A07/partition_metrics.tsv'):sha(RUN/'A07/partition_metrics.tsv')}
records=[];states=[];overlaps=[]
for r in M.itertuples():
    name=f'{r.role}_r{r.resolution:.8f}_seed{r.seed}'
    p=RUN/'A07'/r.unit/(name+'_overlap.tsv');receipts[str(p)]=sha(p)
    t=pd.read_csv(p,sep='\t',index_col=0);c=t.to_numpy(int)
    assert t.index.isin(S.index).all() and c.sum()>0
    n=c.sum();nr=c.sum(1);nc=c.sum(0)
    # Independent ARI computation from pair counts, without rerunning Leiden.
    choose=lambda x:np.asarray(x)*(np.asarray(x)-1)/2
    total=choose(n);a=choose(nr).sum();b=choose(nc).sum();expected=a*b/total
    denom=(a+b)/2-expected
    ari=float((choose(c).sum()-expected)/denom) if denom else 1.
    assert np.isclose(ari,r.ARI,atol=1e-12)
    family=S.loc[t.index,'display_group']
    fc=t.groupby(family.to_numpy()).sum().to_numpy()
    precision=float((fc**2/nc).sum()/n)
    recall=float((fc**2/fc.sum(1)[:,None]).sum()/n)
    records.append(dict(unit=r.unit,resolution=r.resolution,seed=r.seed,role=r.role,
        n_cells=int(n),original_states=len(t),alternative_clusters=len(t.columns),ARI=ari,
        family_Bcubed_precision=precision,family_Bcubed_recall=recall))
    if r.role=='baseline_replay':assert ari>1-1e-12
    for i,state in enumerate(t.index):
        famcounts=c[family.to_numpy()==family.iloc[i]].sum(0)
        states.append(dict(state_id=state,unit=r.unit,resolution=r.resolution,seed=r.seed,role=r.role,
            n_cells=int(nr[i]),best_cluster_jaccard=float(np.max(c[i]/(nr[i]+nc-c[i]))),
            state_Bcubed_precision=float((c[i]**2/nc).sum()/nr[i]),
            state_Bcubed_recall=float((c[i]**2).sum()/nr[i]**2),
            original_family_neighbor_fraction=float((c[i]*famcounts/nc).sum()/nr[i]),
            alternative_clusters_with_5percent=int((c[i]/nr[i]>=.05).sum())))
        for j,count in enumerate(c[i]):
            if count:overlaps.append(dict(unit=r.unit,resolution=r.resolution,seed=r.seed,role=r.role,
                original_state_id=state,alternative_cluster=t.columns[j],cells=int(count)))
save(pd.DataFrame(records),'partition_runs.tsv');st=pd.DataFrame(states);save(st,'state_partition_runs.tsv')
save(pd.DataFrame(overlaps),'partition_overlaps.tsv.gz')
summary=st[st.role.eq('sensitivity')].groupby(['unit','state_id']).agg(
    runs=('seed','size'),worst_jaccard=('best_cluster_jaccard','min'),
    median_jaccard=('best_cluster_jaccard','median'),
    worst_family_neighbor_fraction=('original_family_neighbor_fraction','min'),
    median_family_neighbor_fraction=('original_family_neighbor_fraction','median')).reset_index()
assert len(summary)==111 and summary.runs.eq(9).all();save(summary,'state_partition_summary.tsv')
qpath=BASE/'A06/matched_cell_diagnostics.parquet';receipts[str(qpath)]=sha(qpath)
q=pd.read_parquet(qpath)
assert len(q)==551201 and q.cell_id.is_unique
assert q.groupby('state_key').size().sort_index().equals(S.n_cells.sort_index().astype('int64'))
metrics=['rna_total_umi','rna_detected_genes','rna_pct_mito']
for metric in metrics:
    n=q.groupby('sample_key')[metric].transform('count')
    q[metric+'_within_library_percentile']=(q.groupby('sample_key')[metric].rank(method='average')-.5)/n
keys=['species_key','true_study','state_key','sample_key']
agg={'n_cells':('cell_id','size'),'automatic_call_fraction':('automatic_call','mean')}
for metric in metrics:
    agg[metric+'_median']=(metric,'median')
    agg[metric+'_percentile_mean']=(metric+'_within_library_percentile','mean')
lib=q.groupby(keys).agg(**agg).reset_index()
adjpath=RUN/'A06/UMI_adjusted_state_library_scores.tsv';receipts[str(adjpath)]=sha(adjpath)
adj=read(adjpath);lib=lib.merge(adj[keys+['valid_cells','UMI_stratified_percentile_difference']],on=keys,how='left',validate='one_to_one')
assert len(lib)==875;save(lib,'state_library_QC.tsv')
qc=[]
for (sp,study,state),z in lib.groupby(keys[:3]):
    row=dict(species=sp,study_id=study,state_id=state,n_cells=int(z.n_cells.sum()),n_libraries=len(z),
        n_libraries_ge20=int(z.n_cells.ge(20).sum()),largest_library_cell_fraction=float(z.n_cells.max()/z.n_cells.sum()))
    for subset,v in [('all_libraries',z),('ge20',z[z.n_cells.ge(20)])]:
        for col in [m+'_percentile_mean' for m in metrics]+['automatic_call_fraction','UMI_stratified_percentile_difference']:
            row[subset+'_'+col]=v[col].mean()
    qc.append(row)
qc=pd.DataFrame(qc);assert len(qc)==111;save(qc,'state_QC_summary.tsv')
assert all(sha(Path(p))==h for p,h in receipts.items())
(OUT/'structure_qc_verification.json').write_text(json.dumps(dict(status='PASS',runs=len(records),
    sensitivity_runs=int(M.role.eq('sensitivity').sum()),baseline_replays=int(M.role.eq('baseline_replay').sum()),
    original_states=111,cells=551201,state_library_strata=875,ARI_independently_recomputed=True,
    source_sha256=receipts,interpretation='Structural compatibility with original families; not new biological labels.'),indent=2)+'\n')
print('STRUCTURE_QC_COMPLETE',len(records),len(st),len(lib),flush=True)
