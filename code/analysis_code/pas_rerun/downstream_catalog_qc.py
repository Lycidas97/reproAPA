"""Uniform-catalog QC, with caller-specific filters kept separate."""
import argparse
import json

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_parallel import ordered_map
from downstream_runtime import Context, STUDIES
from downstream_pas_qc import truth, sequence_evidence, read_profiles


def qc_study(ctx,sid):
    import pysam
    intake=json.loads(ctx.register(ctx.native/'intake.native.json').read_text())
    registry=ctx.read(ctx.resource('native_libraries'))
    summaries=[];strata=[];reasons=[];chosen_libraries=[];availability=[]
    folder='WOT_GSE115943' if sid=='GSE115943' else sid
    root=ctx.native/'studies'/folder
    f=ctx.read(root/'features.tsv.gz',low_memory=False)
    meta=ctx.read(root/'metacells.tsv.gz')
    matrix=sparse.load_npz(ctx.register(root/'metacell_by_pas_counts.npz')).tocsr()
    np.testing.assert_array_equal(np.asarray(matrix.sum(0)).ravel(),f.current_target_umi)
    depth=pd.DataFrame(dict(library_key=meta.library_key,UMI=np.asarray(matrix.sum(1)).ravel())).groupby('library_key').UMI.sum()
    rank=depth.rename('UMI').reset_index().sort_values(['UMI','library_key'])
    libkey=rank.iloc[(len(rank)-1)//2].library_key
    lib=registry[registry.library_key.eq(libkey)];assert len(lib)==1
    library=next(lib.itertuples(index=False));sp=f.species.dropna().iloc[0]
    counts=np.asarray(matrix[meta.library_key.eq(libkey).to_numpy()].sum(0)).ravel()
    origins=pd.DataFrame(dict(atlas_column_0based=f.atlas_column_0based,
                             library_id=library.library_id,native_molecules=counts))
    common=intake.get('source_id')=='polyapipe_common'
    if common:
        origins=origins.assign(misprime=f.native_misprime,chrom=f.chrom,strand=f.strand,
            native_cut_site_1based=f.cut_site+f.strand.eq('-').astype(int),site_id=f.native_peak_id)
    chosen_libraries.append(dict(study_id=sid,species=sp,library_key=libkey,
        selection='lower median selected-atlas-cell full-catalog molecule count',selected_UMI=int(counts.sum())))
    layers=[('all',np.ones(len(f),bool)),('primary',truth(f.keep_for_apa_primary)),
            ('supported50',truth(f.keep_for_apa_primary)&f.matched_source_count_50nt.ge(1)),
            ('sensitivity',truth(f.keep_for_sensitivity)),('annotation_excluded',~truth(f.keep_for_apa_primary))]
    for layer,mask in layers:
        q=f[mask];r=dict(study_id=sid,species=sp,layer=layer,sites=len(q),UMI=int(q.current_target_umi.sum()),
                        zero_count_sites=int(q.current_target_umi.eq(0).sum()),genes=q.primary_gene_id.nunique())
        for window in [25,50,100]:
            hit=q[f'matched_source_count_{window}nt'].ge(1)
            r[f'supported{window}_sites']=int(hit.sum());r[f'supported{window}_UMI']=int(q.loc[hit,'current_target_umi'].sum())
        summaries.append(r)
        for col in ['confidence_tier','primary_region_class','primary_assignment_status',
                    'internal_priming_flag','multi_gene_ambiguous','measurement_exclusion']:
            for value,z in q.groupby(col,dropna=False):
                strata.append(dict(study_id=sid,species=sp,layer=layer,field=col,value=value,sites=len(z),UMI=int(z.current_target_umi.sum())))
    for col in ['tier_reason','failure_reason','neighbor_gene_exclusion_reason','measurement_exclusion']:
        for value,z in f.groupby(col,dropna=False):
            reasons.append(dict(study_id=sid,species=sp,field=col,reason=value,sites=len(z),UMI=int(z.current_target_umi.sum())))
    ctx.table(f[['atlas_pas_id','primary_gene_id','confidence_tier','primary_region_class','keep_for_apa_primary',
                 'current_target_umi','quantifiable','coordinate_supported','analysis_view','measurement_exclusion']+
                [f'matched_source_count_{w}nt' for w in [25,50,100]]],folder+'/PAS_QC.tsv.gz')
    fasta_path=ctx.resource(sp+'_genome');ctx.register(str(fasta_path)+'.fai')
    with pysam.FastaFile(str(fasta_path)) as fasta:
        sequence_evidence(ctx,folder,f,origins,library,sp,fasta,native_misprime_available=common)
    # Local read-end QC is not a new count matrix or a six-rule BAM recount.
    read_profiles(ctx,folder,f,origins,library)
    availability.append(dict(study_id=sid,caller=intake['caller'],analysis_route=intake.get('analysis_route','base'),
                              native_misprime='common_GFF_audit_flag' if common else 'not_applicable',
                              ASNR='not_applied' if common else 'parallel_primary_data_route'))
    return summaries,strata,reasons,chosen_libraries,availability


def run(ctx):
    import pysam
    intake=json.loads(ctx.register(ctx.native/'intake.native.json').read_text())
    if intake['caller']=='polyApipe' and intake.get('source_id')!='polyapipe_common':
        from downstream_pas_qc import run as native_run
        return native_run(ctx)
    registry=ctx.read(ctx.resource('native_libraries'))
    summaries=[];strata=[];reasons=[];chosen_libraries=[];availability=[]
    for result in ordered_map(ctx,qc_study,STUDIES):
        for target,values in zip([summaries,strata,reasons,chosen_libraries,availability],result):target.extend(values)
    for data,name in [(summaries,'catalog_summary.tsv'),(strata,'support_and_annotation_strata.tsv'),
                      (reasons,'filter_reason_audit.tsv'),(chosen_libraries,'sequence_read_library_selection.tsv'),
                      (availability,'caller_filter_semantics.tsv')]:
        ctx.table(data,name)
    ctx.finish(dict(all_seven_catalogues=True,site_and_UMI_denominators=True,new_support_and_eligibility=True,
                    caller_filters_not_fabricated=True,local_read_QC_only_no_matrix_recount=True,
                    sequence_and_matched_background_recomputed=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'pas_qc',a.out))
