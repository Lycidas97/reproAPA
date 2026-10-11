"""Build effect-independent strict terminal-UTR/native PAS anchors and read pinned sequences."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TABLES=ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
CONFIG=ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'
CODES={'H1':'GSE118258','H2':'GSE150311','H3':'GSE178325','H4':'GSE221739','W':'GSE115943','C':'GSE103221'}


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def save(frame,name):
    frame.to_csv(HERE/name,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)


def upstream_bounds(cut,strand,span=100):
    return (cut-span,cut) if strand=='+' else (cut,cut+span)


def native_order(cuts,strand):
    return sorted(set(map(int,cuts)),reverse=strand=='-')


def prepare():
    inputs=[TABLES/'recurrence_controls/recurrence_frozen_covariates.tsv.gz',TABLES/'tandem/tandem_utr_eligibility.tsv.gz',TABLES/'tandem/annotation_terminal_utrs.tsv.gz',TABLES/'cis/detected_to_annotation.tsv.gz',TABLES/'cross_overlap/ortholog_matched_features.tsv.gz',CONFIG,HERE/'PROTOCOL.zh.md',HERE/'protocol.json']
    # Response fields are deliberately excluded while defining sequence eligibility.
    f=pd.read_csv(inputs[0],sep='\t',usecols=['species','contrast','threshold','unit','gene_id','gene_symbol','series_pattern','can_negative','can_positive','PAS_number','baseline','mean_GEX','APA_coverage'],float_precision='round_trip')
    f=f[f.contrast.eq('P_M')&f.threshold.eq(.02)].copy()
    e=pd.read_csv(inputs[1],sep='\t');u=pd.read_csv(inputs[2],sep='\t');a=pd.read_csv(inputs[3],sep='\t')
    a['study_id']=a.code.map(CODES);assert a.study_id.notna().all()
    orth=pd.read_csv(inputs[4],sep='\t')[['human_gene_id','mouse_gene_id','human_gene_symbol','mouse_gene_symbol','homology_type','is_high_confidence']]
    assert len(orth)==16092 and orth.homology_type.eq('ortholog_one2one').all() and orth.is_high_confidence.eq(1).all()
    assert orth.human_gene_id.is_unique and orth.mouse_gene_id.is_unique
    f=f[(f.species.eq('human')&f.gene_id.isin(orth.human_gene_id))|(f.species.eq('mouse')&f.gene_id.isin(orth.mouse_gene_id))].copy()
    rows=[];queries=[]
    for (sp,gene),g in f.groupby(['species','gene_id'],sort=True):
        units=sorted(g.unit.unique());z=e[e.species.eq(sp)&e.study_id.isin(units)&e.gene_id.eq(gene)]
        assert len(z)==len(units)
        context=dict(species=sp,gene_id=gene,gene_symbol=g.gene_symbol.iloc[0],eligible_studies=len(units),measured_pattern=';'.join(g.sort_values('unit').unit+':'+g.sort_values('unit').series_pattern),proximal_capacity=int(g.can_negative.sum()),distal_capacity=int(g.can_positive.sum()),PAS_number=float(g.PAS_number.median()),baseline=float(g.baseline.median()),mean_GEX=float(g.mean_GEX.median()),APA_coverage=float(g.APA_coverage.median()),full_primary_axis_tandem=bool(z.tandem_utr_compatible.all()),status='pending')
        if len(units)<2:context['status']='fewer_than_two_studies';rows.append(context);continue
        if not context['full_primary_axis_tandem']:context['status']='non_tandem_full_primary_axis';rows.append(context);continue
        transcript_sets=[set(t.split(';')) for t in z.compatible_transcripts]
        common=set.intersection(*transcript_sets)
        if not common:context['status']='no_shared_full_axis_UTR_transcript';rows.append(context);continue
        sites=a[a.species.eq(sp)&a.gene_id.eq(gene)&a.study_id.isin(units)].copy()
        assert set(sites.study_id)==set(units)
        if sites.chrom.nunique()!=1 or sites.strand.nunique()!=1:context['status']='multichromosome_or_strand';rows.append(context);continue
        strand=sites.strand.iloc[0];chrom=sites.chrom.iloc[0];ordered=native_order(sites.native_cut,strand)
        if len(ordered)<2:context['status']='fewer_than_two_supported_native_PAS';rows.append(context);continue
        pc,dc=ordered[0],ordered[-1];bounds=[upstream_bounds(c,strand) for c in [pc,dc]]
        # Every contributor's full-primary-axis transcript is intersected before checking windows.
        fit=u[u.species.eq(sp)&u.gene_id.eq(gene)&u.transcript_id.isin(common)&u.strand.eq(strand)&u.chrom.astype(str).str.removeprefix('chr').eq(str(chrom).removeprefix('chr'))].copy()
        fit=fit[fit.utr_start.sub(1).le(min(x[0] for x in bounds))&fit.utr_end.ge(max(x[1] for x in bounds))]
        if fit.empty:context['status']='upstream100_crosses_strict_UTR_boundary';rows.append(context);continue
        fit=fit.sort_values('transcript_id');selected=fit.iloc[0]
        context.update(status='sequence_pending',chrom=chrom,strand=strand,proximal_cut=pc,distal_cut=dc,anchor_span=abs(dc-pc),supported_union_PAS=len(ordered),native_cut_list=';'.join(map(str,ordered)),compatible_transcripts=';'.join(sorted(common)),window_compatible_transcripts=';'.join(fit.transcript_id),chosen_transcript=selected.transcript_id,utr_start0=int(selected.utr_start)-1,utr_end0=int(selected.utr_end),upstream100_overlap_nt=max(0,100-abs(dc-pc)))
        rows.append(context)
        for anchor,cut in [('proximal',pc),('distal',dc)]:
            lo,hi=upstream_bounds(cut,strand);queries.append(dict(species=sp,gene_id=gene,anchor=anchor,chrom=chrom,strand=strand,cut_site=cut,start0=lo,end0=hi))
    inventory=pd.DataFrame(rows);query=pd.DataFrame(queries)
    save(inventory,'background_before_sequence.tsv.gz');save(query,'sequence_queries.tsv')
    (HERE/'input_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)):digest(p) for p in inputs},indent=2)+'\n')
    print(inventory.groupby(['species','status']).size().to_string(),flush=True)



if __name__=='__main__':prepare()
