"""Rebuild shared-ortholog display scores on the unchanged cell embedding."""
import argparse
import numpy as np
import pandas as pd
from downstream_runtime import Context,Study
from downstream_moran import scores
from downstream_cross import CODES


def run(ctx):
    roster=ctx.read(ctx.analysis/'cross_overlap/new_shared_members.tsv.gz')
    roster=roster[roster.layer.eq('supported50')&roster.threshold.eq(.02)].copy()
    cells=ctx.read(ctx.resource('frozen_cell_embedding'))
    assert cells.cell_id.is_unique and len(cells)==551201
    cells['join_study']=cells.study_id.replace({'WOT_GSE115943':'GSE115943'})
    frames=[];coverage=[];gene_parts=[]
    for sid in dict.fromkeys(s for s,_ in CODES.values()):
        study=Study(ctx,sid);species=study.f.species.iloc[0];key=species+'_gene_id'
        selected=roster[key].drop_duplicates()
        meta,genes,value,_=scores(study,pd.DataFrame({'gene_id':selected}))
        n=np.isfinite(value).sum(0)
        mean=np.divide(np.nansum(value,axis=0),n,out=np.full(len(genes),np.nan),where=n>0)
        for contrast in ['PS','PM']:
            chosen=set(roster.loc[roster.contrast.eq(contrast),key]);take=genes.isin(chosen)
            measured=np.isfinite(value[:,take]).sum(1)
            centered=np.divide(np.nansum(value[:,take]-mean[take],axis=1),measured,
                out=np.full(len(meta),np.nan),where=measured>0)
            f=pd.DataFrame(dict(study_id=sid,species=species,contrast=contrast,metacell_id=meta.metacell_id,
                centered_DWUI=centered,qualified_genes=measured));frames.append(f)
            g=pd.DataFrame(dict(study_id=sid,contrast=contrast,gene_id=sorted(chosen)))
            g['has_supported_multiPAS_axis']=g.gene_id.isin(genes);gene_parts.append(g)
            subset=cells[cells.join_study.eq(sid)]
            joined=subset.merge(f.drop(columns=['study_id','species']),on='metacell_id',how='left',validate='many_to_one')
            assert len(joined)==len(subset) and joined.qualified_genes.notna().all()
            coverage.append(dict(study_id=sid,species=species,contrast=contrast,selected_genes=len(chosen),
                measurable_axis_genes=int(take.sum()),metacells=len(meta),scored_metacells=int(np.isfinite(centered).sum()),
                displayed_cells=len(joined),scored_cells=int(joined.centered_DWUI.notna().sum()),
                scored_cell_fraction=float(joined.centered_DWUI.notna().mean())))
        print('shared_projection',sid,flush=True)
    ctx.table(pd.concat(frames,ignore_index=True),'shared_metacell_scores.tsv.gz')
    ctx.table(pd.concat(gene_parts,ignore_index=True),'shared_gene_eligibility.tsv')
    ctx.table(coverage,'shared_projection_coverage.tsv')
    ctx.finish(dict(current_shared_rosters=True,frozen_cell_embedding=True,all_six_planned_studies=True,
        supported50_native_counts=True,gene_depth20=True,within_study_gene_centering=True,
        equal_measured_gene_average=True,missing_scores_not_zero=True),
        interpretation='Descriptive projection only; cell scores inherit metacell measurements and are not independent observations')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'shared_projection',a.out))
