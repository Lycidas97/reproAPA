"""Whole-genome native DWUI on every frozen metacell for Figure 1C."""
import argparse
import numpy as np
import pandas as pd

from downstream_runtime import Context,Study,STUDIES
from downstream_moran import scores


def run(ctx):
    frames=[];audit=[]
    for sid in STUDIES:
        study=Study(ctx,sid)
        meta,genes,value,primary=scores(study,pd.DataFrame({'gene_id':list(study.genes)}))
        n=np.isfinite(value).sum(axis=0)
        mean=np.divide(np.nansum(value,axis=0),n,out=np.full(len(genes),np.nan),where=n>0)
        covered=np.isfinite(value).sum(axis=1)
        centered=np.divide(np.nansum(value-mean,axis=1),covered,out=np.full(len(meta),np.nan),where=covered>0)
        raw=np.divide(np.nansum(value,axis=1),covered,out=np.full(len(meta),np.nan),where=covered>0)
        frames.append(pd.DataFrame(dict(metacell_id=meta.metacell_id,study_id=sid,species=study.f.species.iloc[0],
            centered_DWUI=centered,raw_DWUI=raw,qualified_genes=covered,primary_pas_UMI=primary)))
        audit.append(dict(study_id=sid,genes=len(genes),metacells=len(meta),finite_scores=int(np.isfinite(centered).sum())))
        print('GLOBAL_FIGURE_SCORE',sid,flush=True)
    frame=pd.concat(frames,ignore_index=True)
    assert len(frame)==11104 and frame.metacell_id.is_unique
    ctx.table(frame,'metacell_scores.tsv.gz');ctx.table(audit,'coverage.tsv')
    ctx.finish(dict(frozen_metacell_membership=True,full_primary_axis=True,supported50_native_counts=True,
        gene_depth20=True,within_study_gene_centering=True,equal_measured_gene_average=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'figure_scores',a.out))
