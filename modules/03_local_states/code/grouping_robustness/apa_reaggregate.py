"""Remove predeclared cells, retaining PAS axes and the original metacell framework."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
from pathlib import Path
import json,hashlib,sys,time
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data'
BASE=Path('/root/nfsdata/apadetect/analysis_projects')
PRIOR=BASE/'reproapa_pas_reliability_2026_09_13'
MOUSE=BASE/'mouse_time_aware_integration_pilot_v1'
RUN=BASE/'reproapa_integration_validation_20260913/execution_v1'
sys.path.insert(0,str(ROOT/'inputs/estimator_source'))
from metrics import Estimator

TARGETS={'WOT_L00':'wot_common_to_2i|WOT::common_to_2i::L00',
    'CIP_C012':'gse103221_cip|C012','human_C005':'GSE118258::studywide::C005',
    'human_C011':'GSE118258::studywide::C011','human_C004':'GSE118258::studywide::C004',
    'human_C006':'GSE118258::studywide::C006'}
VARIANTS={'baseline':[], 'exclude_WOT_L00':['WOT_L00'], 'exclude_CIP_C012':['CIP_C012'],
    'exclude_human_C005':['human_C005'], 'exclude_three':['WOT_L00','CIP_C012','human_C005'],
    'exclude_human_candidates':['human_C011','human_C004','human_C005','human_C006']}
STUDIES=['gse115943_pluripotency_reprogramming_core','gse103221_pluripotency_reprogramming_core',
    'gse118258_human_pluripotency_reprogramming']
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False)
def receipt(p,records):
    records[str(p)]=dict(sha256=sha(p),bytes=p.stat().st_size)

def main():
    OUT.mkdir(exist_ok=True)
    frozen=dict(targets=TARGETS,variants=VARIANTS,studies=STUDIES,
        contract_sha256=sha(ROOT/'analysis-contract.zh.md'),script_sha256=sha(Path(__file__)))
    (OUT/'apa_frozen_variants.json').write_text(json.dumps(frozen,indent=2)+'\n')
    atlaspath=RUN/'shared/atlas.parquet';atlas=pd.read_parquet(atlaspath)
    annpath=ROOT/'inputs/S5-local-state-annotations.tsv';ann=read(annpath).set_index('original_state_id')
    atlas['family']=atlas.state_key.map(ann.original_family_code.replace({'T':'P'}))
    assert atlas.family.notna().all() and set(TARGETS.values())<=set(atlas.state_key)
    receipts={};receipt(atlaspath,receipts);receipt(annpath,receipts)
    for p in (ROOT/'inputs/estimator_source').glob('*.py'):receipt(p,receipts)
    audit=[];cell_audit=[];mc_audit=[];position_parts=[]
    for study in STUDIES:
        dst=OUT/study;dst.mkdir(exist_ok=True)
        source=PRIOR/'revision/data'/study;original=PRIOR/'data'/study
        fp=source/'features.tsv.gz';mp=original/'selected_metacells.tsv.gz';rp=original/'matrix_rows.tsv.gz';lp=source/'libraries.tsv'
        for p in [fp,mp,rp,lp]:receipt(p,receipts)
        f=read(fp);meta=read(mp);rows=read(rp);libs=read(lp)
        species=str(f.species.iloc[0]);nmc=len(rows);nf=len(f)
        state_universe=set(atlas.loc[atlas.sample_key.isin(libs.library_key)|atlas.library_key.isin(libs.library_key),'state_key'])
        local_sets={name:tuple(sorted(set(TARGETS[k] for k in targets)&state_universe)) for name,targets in VARIANTS.items()}
        sets=sorted(set(local_sets.values()))
        counts={s:sparse.csr_matrix((nmc,nf),dtype=np.int64) for s in sets}
        ncells={s:np.zeros(nmc,dtype=int) for s in sets};same={s:np.zeros(nmc,dtype=int) for s in sets}
        node=np.full(nmc,'_',dtype=object);node[meta.matrix_row]=meta.node
        for j,lib in enumerate(libs.itertuples(),1):
            folder=source/'recount'/lib.library_key
            paths=[folder/'cells.tsv.gz',folder/'original_cells.npz',folder/'original_mc.npz']
            for p in paths:receipt(p,receipts)
            bind=read(paths[0]);raw=sparse.load_npz(paths[1]).tocsr();expected=sparse.load_npz(paths[2]).tocsr()
            assert raw.shape==(len(bind),nf) and bind.barcode.is_unique
            local=atlas[(atlas.library_key.eq(lib.library_key))|atlas.sample_key.eq(lib.library_key)].set_index('barcode')
            assert local.index.is_unique and len(bind)==len(local) and bind.barcode.isin(local.index).all()
            local=local.loc[bind.barcode].reset_index(drop=True)
            mci=pd.Index(rows.metacell_id).get_indexer(bind.metacell_id);assert (mci>=0).all()
            ag=sparse.csr_matrix((np.ones(len(bind),int),(mci,np.arange(len(bind)))),shape=(nmc,len(bind)))
            baseline=(ag@raw).tocsr();diff=baseline-expected;diff.eliminate_zeros();assert diff.nnz==0
            for states in sets:
                keep=~local.state_key.isin(states).to_numpy()
                counts[states]=counts[states]+(baseline if not states else ag[:,keep]@raw[keep])
                ncells[states]+=np.bincount(mci[keep],minlength=nmc)
                correct=keep & (local.family.to_numpy()==node[mci])
                same[states]+=np.bincount(mci[correct],minlength=nmc)
                if states:
                    cell_audit.append(dict(study=study,library_key=lib.library_key,excluded_states=';'.join(states),
                        all_cells=len(bind),removed_cells=int((~keep).sum()),remaining_cells=int(keep.sum()),
                        original_UMI=int(raw.sum()),removed_UMI=int(raw[~keep].sum())))
            print('CELL_REAGGREGATION',study,j,len(libs),lib.library_key,flush=True)
        ix=meta.matrix_row.to_numpy(int)
        assert np.array_equal(ncells[()][ix],meta.n_cells),'Baseline metacell cell-count mismatch'
        purity=np.divide(same[()][ix],ncells[()][ix])
        purity_error=float(np.max(np.abs(purity-meta.purity.to_numpy())))
        assert np.allclose(purity,meta.purity,atol=1e-12,rtol=0),('Baseline family-purity mismatch',study,purity_error)
        estimates={}
        for states in sets:
            for policy in ['reapply_purity','frozen_purity']:
                current=meta.copy();current['n_cells']=ncells[states][ix]
                current['retained_purity']=np.divide(same[states][ix],ncells[states][ix],out=np.zeros(len(ix)),where=ncells[states][ix]>0)
                eligible=(current.retained_purity.ge(.8) if policy=='reapply_purity' else np.ones(len(current),bool))&current.n_cells.gt(0)
                selected=current[eligible].copy()
                value=Estimator(f,selected).calculate(counts[states]);estimates[states,policy]=value
                entry=current[['metacell_id','library_key','route','node','matrix_row','purity','n_cells','retained_purity']].copy()
                entry['study']=study;entry['excluded_states']=';'.join(states);entry['policy']=policy;entry['passes_cell_and_purity']=eligible.to_numpy()
                entry['original_n_cells']=meta.n_cells.to_numpy();mc_audit.append(entry)
        base=estimates[(),'reapply_purity']
        ep=MOUSE/'results_apa_majority_shared_v1/profiles'/study/'supported50_library_positions.tsv.gz';receipt(ep,receipts)
        expected=read(ep);keys=['route','node','library_key','gene_id']
        a=base.dropna(subset=['position_index']).set_index(keys).position_index.sort_index()
        b=expected[expected.gene_id.isin(f.primary_gene_id)].dropna(subset=['position_index']).set_index(keys).position_index.sort_index()
        assert a.index.equals(b.index) and np.allclose(a,b,atol=1e-10,rtol=0),('APA baseline mismatch',study)
        for variant,states in local_sets.items():
            for policy in ['reapply_purity','frozen_purity']:
                value=estimates[states,policy].copy();value['variant']=variant;value['policy']=policy
                value['study']=study;value['species']=species;position_parts.append(value)
        audit.append(dict(study=study,libraries=len(libs),cells=int(ncells[()].sum()),features=nf,
            baseline_purity_max_abs_difference=purity_error,baseline_PA_max_abs_difference=float(np.max(np.abs(a-b))),
            baseline_raw_MC_counts_exact=True,baseline_library_PA_rows=len(a)))
        save(pd.concat([v.assign(excluded_states=';'.join(k[0]),policy=k[1]) for k,v in estimates.items()],ignore_index=True),dst/'library_positions.tsv.gz')
        (OUT/'apa_input_receipt.partial.json').write_text(json.dumps(receipts,indent=2)+'\n')
        del counts,estimates
    save(pd.concat(position_parts,ignore_index=True),OUT/'affected_library_positions.tsv.gz')
    save(pd.DataFrame(cell_audit),OUT/'removed_cell_library_counts.tsv');save(pd.concat(mc_audit,ignore_index=True),OUT/'metacell_retention.tsv.gz')
    for p,h in receipts.items():assert sha(Path(p))==h['sha256'],('Source changed',p)
    (OUT/'apa_reaggregation_verification.json').write_text(json.dumps(dict(status='PASS',studies=audit,
        exact_cell_level_exclusion=True,original_PAS_and_metacell_framework_preserved=True,
        source_receipts=receipts,variants=frozen,all_source_hashes_unchanged=True),indent=2)+'\n')
    print('APA_REAGGREGATION_COMPLETE',flush=True)

if __name__=='__main__':main()
