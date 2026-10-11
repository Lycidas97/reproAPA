"""Select new recurrent RBP candidates before querying unchanged GEX counts."""
import argparse

import numpy as np
import pandas as pd

from downstream_runtime import Context,Study
from downstream_regulatory import ROUTES

KEY=['species','study_id','route','factor']


def run(ctx):
    recurrent=ctx.read(ctx.analysis/'site_statistics/exact_site_cross_study.tsv.gz')
    recurrent=recurrent[recurrent.category.eq('RBP')&recurrent.recurrent_strict].copy();records=[]
    for row in recurrent.to_dict('records'):
        # A physical match may have both direct and transferred provenance.
        # Each evidence class keeps its own inventory; expression is queried once.
        shared=set(row['provenance_a'].split(';'))&set(row['provenance_b'].split(';'))
        for code,evidence in [('RBP','direct_species'),('RBP_transfer','transferred_motif')]:
            if code in shared:records.append(row|dict(evidence=evidence))
    sites=pd.DataFrame(records,columns=list(recurrent.columns)+['evidence'])
    sites=sites.drop_duplicates(['species','factor','evidence','gene_id','chrom','strand','start0','end0','contrast'])
    candidates=sites.groupby(['species','factor','evidence']).agg(target_genes=('gene_id','nunique'),site_contrast_links=('factor','size')).reset_index()
    candidates['display_selected']=False
    for sp in ROUTES:
        idx=candidates[candidates.species.eq(sp)&candidates.evidence.eq('direct_species')].sort_values(['target_genes','factor'],ascending=[False,True]).index
        candidates.loc[idx if sp=='mouse' else idx[:6],'display_selected']=True
    ctx.table(candidates,'candidates.tsv');ctx.table(sites,'candidate_sites.tsv.gz')
    num=ctx.functions('notes/2026-10-02-rbp-expression/analyze.py',['library_rows','state_rows','contrast_rows','site_rows'],{'KEY':KEY})
    libraries=[];mapping=[]
    for sp,routes in ROUTES.items():
        factors=sorted(candidates[candidates.species.eq(sp)].factor.unique())
        if not factors:continue
        for sid,route in routes:
            study=Study(ctx,sid);study.load_gex();meta=study.meta[study.meta.route.eq(route)&study.meta.purity.ge(.8)].copy()
            symbols=study.gex_symbols().str.upper();values=np.full((len(meta),len(factors)),-1,dtype=np.int64)
            for j,factor in enumerate(factors):
                idx=np.flatnonzero(symbols.eq(factor.upper()))
                if len(idx)==1:values[:,j]=study.gex_counts[:,idx[0]].toarray().ravel()[meta.index]
                mapping.append(dict(species=sp,study_id=sid,factor=factor,gene_ids=';'.join(study.gex_genes.gene_id.iloc[idx]),
                    status='unique_exact_casefold_symbol' if len(idx)==1 else 'missing_symbol' if not len(idx) else 'ambiguous_symbol'))
            total=study.gex_total[meta.index];meta=meta.reset_index(drop=True)
            libraries.append(num.library_rows(meta,values,total,factors,sid,sp))
    lib=pd.concat(libraries,ignore_index=True) if libraries else pd.DataFrame(columns=KEY+['library_key','state','time','n_metacells','n_cells','qualified','mapped','count','total_gex','pooled_CPM','mean_log2CPM1','mean_CPM','detected_metacells','detected_fraction'])
    if len(lib):
        states=num.state_rows(lib);contrasts,pairs=num.contrast_rows(lib,states)
    else:
        states=pd.DataFrame(columns=KEY+['state','n_libraries']);contrasts=pd.DataFrame(columns=KEY+['contrast','delta_log2CPM1']);pairs=pd.DataFrame(columns=KEY+['library_key','delta_log2CPM1'])
    if not len(pairs.columns):pairs=pd.DataFrame(columns=KEY+['library_key','delta_log2CPM1'])
    for frame,name in [(lib,'library_expression.tsv.gz'),(states,'state_expression.tsv'),(contrasts,'expression_contrasts.tsv'),(pairs,'paired_library_effects.tsv')]:
        frame=frame.merge(candidates,on=['species','factor'],validate='many_to_many');ctx.table(frame,name)
    retention=num.site_rows(sites) if len(sites) else pd.DataFrame(columns=['species','factor','evidence','contrast','route','threshold','n_target_genes','n_sites','weighted_site_change'])
    ctx.table(retention,'site_retention.tsv');ctx.table(mapping,'gene_mapping.tsv',columns=['species','study_id','factor','gene_ids','status'])
    cross=[]
    for key,z in contrasts.groupby(['species','factor','contrast']):
        assert len(z)==2;z=z.sort_values('study_id');a,b=z.iloc[0],z.iloc[1];da,db=a.delta_log2CPM1,b.delta_log2CPM1
        status='unavailable' if not np.isfinite(da+db) else 'both_increased' if min(da,db)>1e-12 else 'both_decreased' if max(da,db)<-1e-12 else 'zero_or_tied' if min(abs(da),abs(db))<=1e-12 else 'opposite'
        cross.append(dict(zip(['species','factor','contrast'],key),study_a=a.study_id,study_b=b.study_id,effect_a=da,effect_b=db,
            time_effect_a=a.time_equal_delta,time_effect_b=b.time_equal_delta,direction=status,
            expressed_in_both_studies=bool(a.expression_in_either_state and b.expression_in_either_state),
            sign_stable_time_weighting=bool(np.isfinite(da+db+a.time_equal_delta+b.time_equal_delta) and np.sign(da)==np.sign(a.time_equal_delta) and np.sign(db)==np.sign(b.time_equal_delta))))
    ctx.table(cross,'cross_study_summary.tsv',columns=['species','factor','contrast','study_a','study_b','effect_a','effect_b','direction'])
    ctx.finish(dict(new_candidates_selected_before_expression=True,direct_transfer_separate=True,all_gene_CPM_denominator=True,
        library_equal_and_time_equal=True,PM_same_physical_libraries=True,expression_flags_never_filter_candidates=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'rbp_expression',a.out))
