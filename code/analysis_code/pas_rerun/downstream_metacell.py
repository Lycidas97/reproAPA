"""Recompute every PAS coverage field on unchanged metacell membership."""
import argparse
import json

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_parallel import ordered_map
from downstream_runtime import Context


def coverage_study(ctx,rec):
    summaries=[]
    root=ctx.native/'studies'/rec['study_id'];f=ctx.read(root/'features.tsv.gz',low_memory=False);m=ctx.read(root/'metacells.tsv.gz')
    x=sparse.load_npz(ctx.register(root/'metacell_by_pas_counts.npz')).tocsr()
    assert x.shape==(len(m),len(f)) and np.issubdtype(x.dtype,np.integer) and (x.data>=0).all()
    assert int(x.sum())==rec['molecules'] and len(m)==rec['metacells'] and int(m.n_cells.sum())==rec['cells']
    assert np.array_equal(np.asarray(x.sum(0)).ravel(),f.current_target_umi)
    # The imported registry contains legacy PAS coverage columns; never carry those into a new result.
    columns=[c for c in m if 'pas' not in c.lower() and not c.startswith('genes_with_ge2')]
    out=m[columns].copy();out['total_pas_umi']=np.asarray(x.sum(1)).ravel();out['detected_pas']=np.asarray((x>0).sum(1)).ravel()
    primary=f.keep_for_apa_primary.astype(str).str.lower().eq('true').to_numpy()
    sensitivity=f.keep_for_sensitivity.astype(str).str.lower().eq('true').to_numpy()
    for label,mask in [('primary',primary),('sensitivity',sensitivity),('relaxed_only',sensitivity&~primary),
        ('supported50',primary&f.matched_source_count_50nt.ge(1).to_numpy())]:
        mat=x[:,mask];out[label+'_pas_umi']=np.asarray(mat.sum(1)).ravel();out[label+'_detected_pas']=np.asarray((mat>0).sum(1)).ravel()
    chosen=f[primary&f.primary_gene_id.notna().to_numpy()];cols=chosen.feature_index_1based.to_numpy()-1;genes=pd.Index(sorted(chosen.primary_gene_id.unique()))
    projection=sparse.csr_matrix((np.ones(len(chosen),np.int64),(np.arange(len(chosen)),genes.get_indexer(chosen.primary_gene_id))),shape=(len(chosen),len(genes)))
    counts=x[:,cols]@projection;detected=(x[:,cols]>0).astype(np.int64)@projection
    out['primary_detected_genes']=np.asarray((counts>0).sum(1)).ravel()
    out['genes_with_ge2_detected_primary_sites']=np.asarray((detected>=2).sum(1)).ravel()
    eligible=np.asarray(projection.sum(0)).ravel()>=2
    for threshold in [10,20,50]:out[f'genes_with_ge2_primary_atlas_sites_and_ge{threshold}_umi']=np.asarray((counts[:,eligible]>=threshold).sum(1)).ravel()
    out['study_id']=rec['study_id'];out['species']=rec['species'];ctx.table(out,rec['study_id']+'/metacell_pas_coverage.tsv.gz')
    summaries.extend(dict(study_id=rec['study_id'],species=rec['species'],library_key=lib,metacells=len(q),cells=int(q.n_cells.sum()),
        molecules=int(q.total_pas_umi.sum()),primary_molecules=int(q.primary_pas_umi.sum()),supported50_molecules=int(q.supported50_pas_umi.sum()),
        median_qualified_primary_genes=q.genes_with_ge2_primary_atlas_sites_and_ge20_umi.median()) for lib,q in out.groupby('library_key'))
    return summaries


def run(ctx):
    audit=json.loads(ctx.register(ctx.native/'frozen_native_membership_audit.json').read_text())
    assert audit['status']=='PASS'
    intake=json.loads(ctx.register(ctx.native/'intake.native.json').read_text());summaries=[]
    for values in ordered_map(ctx,coverage_study,intake['studies']):summaries.extend(values)
    ctx.table(summaries,'library_coverage.tsv')
    ctx.finish(dict(frozen_membership_audit_passed=True,integer_total_and_column_conservation=True,
        all_PAS_coverage_fields_recomputed=True,no_GEX_or_membership_changes=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'metacell',a.out))
