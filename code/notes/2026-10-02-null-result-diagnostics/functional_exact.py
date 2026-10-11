"""Exact convolution for every term in the frozen stratified sampling tests."""
from pathlib import Path
import os, sys, json, importlib.util, hashlib
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'notes/2026-10-02-reviewer-followup'))
from apa_followup import save, read, INPUTS, nullmod

G=Path('/Users/lycidas/work/科研/projects/2026-09-26_reproapa-go-analysis')
os.environ['REPROAPA_GEX_ROOT']=str(Path('/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot/analysis/figure3-apa-gex-state-analysis-2026-09-19/raw'))
spec=importlib.util.spec_from_file_location('frozen_kegg',G/'scripts/run_kegg_nested.py')
K=importlib.util.module_from_spec(spec);spec.loader.exec_module(K)
N=K.NESTED;B=N.BASE


@lru_cache(maxsize=200000)
def pmf(npop, success, draws):
    lo=max(0,draws+success-npop); hi=min(draws,success)
    values=hypergeom.pmf(np.arange(lo,hi+1),npop,success,draws)
    assert abs(values.sum()-1)<1e-9
    return lo,values/values.sum()


def exact(strata,hits):
    distribution=np.array([1.]); offset=0; mean=0.; variance=0.
    for npop,success,draws in strata:
        if not success or not draws:
            continue
        mean += draws*success/npop
        if npop>1:
            variance += draws*(success/npop)*(1-success/npop)*(npop-draws)/(npop-1)
        lo,values=pmf(npop,success,draws)
        offset+=lo
        distribution=np.convolve(distribution,values)
    assert abs(distribution.sum()-1)<1e-8
    p=float(distribution[max(0,hits-offset):].sum())
    return min(1.,p),mean,np.sqrt(variance)


def main():
    # qcut boundary behavior changed between the installed pandas versions.
    # Require the environment that exactly replays the archived null strata.
    assert pd.__version__ == '2.2.3', 'Use the frozen pandas 2.2.3 runtime documented in RESULTS.zh.md'
    (HERE/'functional_environment.json').write_text(json.dumps(dict(pandas=pd.__version__,numpy=np.__version__,python=sys.version),indent=2)+'\n')
    for p in [G/'scripts/run_go.py',G/'scripts/run_go_nested.py',G/'scripts/run_kegg_nested.py']:
        INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    for p in Path(os.environ['REPROAPA_GEX_ROOT']).glob('*/full_gex_states.tsv.gz'):
        INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    for p in list((G/'sources').glob('*'))+list((G/'sources_kegg').glob('*')):
        if p.is_file(): INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    g0=N.add_expression_strata(N.build_g0())
    events=read(G/'data/event_baseline_matrix.tsv.gz')
    membership=read(ROOT/'table/Supplementary-Source-Data-20260928/functional/GO/nested_universe_membership.tsv.gz')
    terms,aliases,ancestors=B.ontology()
    annotations={'GO':{},'KEGG':{}}
    for species in ['mouse','human']:
        annotations['GO'][species]=B.load_go_genes(species,terms,aliases,ancestors)
        annotations['KEGG'][species]=K.load_kegg(species)[0]
    summaries=[]; frozen=[]
    for database, filename in [('GO','nested_go_enrichment.tsv.gz'),('KEGG','kegg_enrichment.tsv.gz')]:
        old=read(ROOT/f'table/Supplementary-Source-Data-20260928/functional/{database}/{filename}')
        output=[]
        for analysis_id,group in old.groupby('analysis_id',sort=False):
            species,contrast,comparison,target=analysis_id.split('|')
            mem=membership[membership.species.eq(species)&membership.contrast.eq(contrast)]
            apa=B.add_bins(B.covariate_table(events,species,contrast))
            apa=apa[apa.eligible_units.ge(2)&apa.gene_id.isin(mem[mem.in_G1_nested].gene_id)]
            expr=g0[g0.species.eq(species)&g0.contrast.eq(contrast)]
            assert set(expr.gene_id)==set(mem.gene_id)
            bg=apa.copy() if comparison=='G2_vs_G1' else expr.copy()
            col='stratum_full' if comparison=='G2_vs_G1' else 'expression_stratum'
            ann=annotations[database][species]
            target_ids=set(mem[mem.in_G1_nested if comparison=='G1_vs_G0' else mem.in_G2].gene_id)
            if target!='all':
                target_ids &= set(mem[mem.recurrent_direction.eq('negative' if target=='proximal' else 'positive')].gene_id)
            original_bg_n=len(bg);original_target_n=len(target_ids)
            bg=bg[bg.gene_id.map(lambda g:bool(ann.get(g)))].copy()
            fg=target_ids&set(bg.gene_id)
            assert len(bg)==group.N.iloc[0] and len(fg)==group.n.iloc[0]
            bg['stratum']=bg[col];bg['target']=bg.gene_id.isin(fg)
            frozen.append(bg[['gene_id','stratum','target']].assign(database=database,analysis_id=analysis_id))
            pools={key:set(part.gene_id) for key,part in bg.groupby('stratum')}
            family=N.term_family(set(bg.gene_id),ann)
            assert set(family)==set(group.term_id)
            data=[]
            for row in group.itertuples():
                term=family[row.term_id]
                assert len(term)==row.K and len(term&fg)==row.hits
                strata=[(len(pool),len(pool&term),len(pool&fg)) for pool in pools.values()]
                p,mu,sd=exact(strata,row.hits)
                data.append(dict(analysis_id=analysis_id,term_id=row.term_id,p_exact=p,matched_mean=mu,matched_sd=sd,
                    naive_mean=row.n*row.K/row.N,excess=row.hits-mu,
                    matched_fold=row.hits/mu if mu else np.nan))
            extra=pd.DataFrame(data)
            merged=group.merge(extra,on=['analysis_id','term_id'],validate='one_to_one')
            merged['q_exact_family']=nullmod.bh(merged.p_exact)
            output.append(merged)
            summaries.append(dict(database=database,analysis_id=analysis_id,terms=len(group),
                background_unannotated=original_bg_n-len(bg),target_unannotated=original_target_n-len(fg),
                target_original=original_target_n,target_annotated=len(fg),background_annotated=len(bg),
                nominal_mc=int((group.p_matched<.05).sum()),nominal_exact=int((merged.p_exact<.05).sum()),
                significant_hyper=int((group.q_hyper_family<.05).sum()),significant_exact=int((merged.q_exact_family<.05).sum()),
                minimum_mc=group.p_matched.min(),minimum_exact=merged.p_exact.min(),minimum_exact_q=merged.q_exact_family.min(),
                min_rank_to_pass_at_mc_floor=int(np.ceil(len(group)/10001/.05)),
                median_expectation_ratio=(merged.matched_mean/merged.naive_mean).median(),
                original_hyper_hits_median_expectation_ratio=(merged.loc[merged.q_hyper_family<.05,'matched_mean']/merged.loc[merged.q_hyper_family<.05,'naive_mean']).median()))
            print(database,analysis_id,'exact q min',merged.q_exact_family.min(),flush=True)
        out=pd.concat(output)
        out['q_exact_grid']=nullmod.bh(out.p_exact)
        save(out,HERE/f'{database.lower()}_exact.tsv.gz')
    save(pd.concat(frozen),HERE/'functional_frozen_strata.tsv.gz')
    save(pd.DataFrame(summaries),HERE/'functional_summary.tsv')
    (HERE/'functional_input_sha256.json').write_text(json.dumps(INPUTS,indent=2)+'\n')


if __name__=='__main__': main()
