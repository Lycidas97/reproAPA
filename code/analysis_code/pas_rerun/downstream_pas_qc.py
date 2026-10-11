"""Native catalogue support, sequence controls and read-end profiles for all studies."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_runtime import Context,STUDIES,seed
from downstream_recount import canon
from polyapipe_prepare import digest


def truth(x):return x.astype(str).str.lower().eq('true')


def sequence(fasta,chrom,cut,strand,refs=None):
    refs={canon(c):c for c in fasta.references} if refs is None else refs;chrom=refs.get(canon(chrom));center=int(cut)-int(strand=='+')
    if chrom is None or center<100 or center+101>fasta.get_reference_length(chrom):return None
    seq=fasta.fetch(chrom,center-100,center+101).upper()
    if strand=='-':seq=seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
    return seq if len(seq)==201 and set(seq)<=set('ACGT') else None


def region_class(model,point,strand):
    if not model.get('transcript'):return 'no_transcript'
    a,b=model['transcript'][0]
    if not a<=point<b:return 'outside_transcript'
    if any(a<=point<b for a,b in model.get('CDS',[])):return 'CDS'
    if any(a<=point<b for a,b in model.get('exon',[])):
        coding=model.get('CDS',[])
        if not coding:return 'noncoding_exon'
        lo=min(a for a,b in coding);hi=max(b for a,b in coding)
        if (strand=='+' and point>=hi) or (strand=='-' and point<lo):return 'three_prime_UTR'
        if (strand=='+' and point<lo) or (strand=='-' and point>=hi):return 'five_prime_UTR'
        return 'other_exon'
    return 'intron'


def near(cuts,point):
    j=np.searchsorted(cuts,point)
    return (j<len(cuts) and abs(cuts[j]-point)<=100) or (j>0 and abs(cuts[j-1]-point)<=100)


def sequence_evidence(ctx,folder,f,origins,library,sp,fasta,*,native_misprime_available=True):
    primary=f[truth(f.keep_for_apa_primary)].copy();curves=[];summaries=[];metrics=[]
    references_by_chrom={canon(c):c for c in fasta.references}
    def evaluate(frame,category,alignment):
        motif=np.zeros(196,int);adenine=np.zeros(201,int);n=0;missing=0;up=0;down=[];runs=[]
        for row in frame.itertuples():
            seq=sequence(fasta,row.chrom,row.cut_site,row.strand,references_by_chrom)
            if seq is None:missing+=1;continue
            n+=1;adenine+=np.frombuffer(seq.encode(),dtype='S1')==b'A'
            hits=np.array([seq[i:i+6] in ['AATAAA','ATTAAA'] for i in range(196)])
            motif+=hits;up+=hits[60:91].any();a=seq[101:121].count('A')/20;run=max([len(v) for v in re.findall('A+',seq[101:121])]+[0]);down.append(a);runs.append(run)
            if alignment=='representative':metrics.append(dict(category=category,atlas_pas_id=row.atlas_pas_id,cut_interbase=row.cut_site,
                upstream_signal_40_10=bool(hits[60:91].any()),downstream_A_fraction=a,downstream_max_A_run=run,GC=(seq.count('G')+seq.count('C'))/201))
        for feature,counts in [('canonical_signal_start',motif),('adenine',adenine)]:
            curves.extend(dict(category=category,alignment=alignment,feature=feature,position_nt=i-100,count=int(v),denominator=n,
                               frequency=v/n if n else np.nan) for i,v in enumerate(counts))
        summaries.append(dict(category=category,alignment=alignment,input_sites=len(frame),extracted_sites=n,missing_sequence=missing,
            upstream_signal_sites=int(up),upstream_signal_fraction=up/n if n else np.nan,mean_downstream_A_fraction=np.mean(down) if down else np.nan,
            median_downstream_A_run=np.median(runs) if runs else np.nan))
    for name,frame in [('supported50',primary[primary.matched_source_count_50nt.ge(1)]),('unsupported50',primary[primary.matched_source_count_50nt.eq(0)])]:
        for alignment,col in [('representative','cut_site'),('member_min','member_cut_site_min'),('member_max','member_cut_site_max')]:
            evaluate(frame.assign(cut_site=frame[col]),name,alignment)
    if native_misprime_available:
        rejected=origins[origins.library_id.eq(library.library_id)&truth(origins.misprime)].copy()
        rejected['cut_site']=rejected.native_cut_site_1based-rejected.strand.eq('-').astype(int)
        rejected['atlas_pas_id']=rejected.site_id
        rejected=rejected.drop_duplicates(['chrom','strand','cut_site'])
        evaluate(rejected,'native_misprime_flagged_in_median_library','representative')
    else:
        # This caller does not have polyApipe's misprime flag. Absence is not a
        # false flag or an observed zero-rate estimate.
        rejected=pd.DataFrame(columns=['atlas_pas_id','chrom','strand','cut_site'])
    excluded=f[~truth(f.keep_for_apa_primary)&f.atlas_column_0based.isin(origins.loc[origins.library_id.eq(library.library_id),'atlas_column_0based'])]
    evaluate(excluded,'annotation_excluded_in_median_library','representative')
    rng=np.random.default_rng(seed('sequence_background',folder));chosen=primary.iloc[np.sort(rng.choice(len(primary),min(1000,len(primary)),replace=False))]
    transcripts=set(chosen.primary_transcript_id.dropna());models=defaultdict(lambda:defaultdict(list))
    with gzip.open(ctx.resource(sp+'_annotation_gtf'),'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if fields[2] not in ['transcript','exon','CDS']:continue
            match=re.search(r'transcript_id "([^"]+)"',fields[8])
            if match and match[1] in transcripts:models[match[1]][fields[2]].append((int(fields[3])-1,int(fields[4])))
    allcuts={k:np.sort(g.cut_site.to_numpy()) for k,g in primary.groupby(['chrom','strand'])};references={}
    for key,g in primary.groupby(['chrom','strand']):
        values=g[[c for c in f if c.endswith('_nearest_cut')]].to_numpy().ravel();references[key]=np.unique(values[np.isfinite(values)&(values>=0)])
    controls=[];ledger=[]
    for row in chosen.itertuples():
        seq=sequence(fasta,row.chrom,row.cut_site,row.strand,references_by_chrom);model=models.get(row.primary_transcript_id,{})
        cls=region_class(model,row.cut_site-int(row.strand=='+'),row.strand);status='no_match';matched=None
        if seq and cls not in ['no_transcript','outside_transcript']:
            gc=(seq.count('G')+seq.count('C'))/201;start,end=model['transcript'][0]
            for point in rng.integers(start,end,size=100):
                cut=int(point)+int(row.strand=='+');key=(row.chrom,row.strand)
                if region_class(model,point,row.strand)!=cls or near(allcuts[key],cut) or near(references[key],cut):continue
                candidate=sequence(fasta,row.chrom,cut,row.strand,references_by_chrom)
                if candidate and abs((candidate.count('G')+candidate.count('C'))/201-gc)<=.05:
                    status='matched';matched=cut;controls.append(dict(atlas_pas_id=row.atlas_pas_id,chrom=row.chrom,strand=row.strand,cut_site=cut));break
        else:status='unavailable_'+cls if seq else 'unavailable_sequence'
        ledger.append(dict(atlas_pas_id=row.atlas_pas_id,transcript_id=row.primary_transcript_id,position_class=cls,cut_interbase=row.cut_site,control_cut_interbase=matched,status=status))
    evaluate(pd.DataFrame(controls,columns=['atlas_pas_id','chrom','strand','cut_site']),'matched_genomic_background','representative')
    for frame,name,columns in [(curves,'sequence_curves.tsv.gz',None),(summaries,'sequence_summary.tsv',None),
                               (metrics,'sequence_site_metrics.tsv.gz',['category','atlas_pas_id']),
                               (ledger,'sequence_background_matches.tsv.gz',['atlas_pas_id','status'])]:ctx.table(frame,folder+'/'+name,columns=columns)
    ctx.table(rejected,folder+'/native_misprime_candidates.tsv.gz')
    ctx.table([dict(native_misprime_available=native_misprime_available,
                   role='caller_specific_audit_only')],folder+'/caller_QC_availability.tsv')


def read_profiles(ctx,folder,f,origins,library):
    import pysam
    primary=f[truth(f.keep_for_apa_primary)].copy();primary['isolated']=False
    for _,q in primary.groupby(['chrom','strand']):
        q=q.sort_values('cut_site');d=np.diff(q.cut_site);primary.loc[q.index,'isolated']=(np.r_[np.inf,d]>400)&(np.r_[d,np.inf]>400)
    values=origins[origins.library_id.eq(library.library_id)].groupby('atlas_column_0based').native_molecules.sum()
    primary['library_UMI']=primary.atlas_column_0based.map(values).fillna(0)
    rng=np.random.default_rng(seed('read_profiles',folder));selections=[];ledger=[]
    for category,supported in [('supported50',True),('unsupported50',False)]:
        frame=primary[primary.matched_source_count_50nt.ge(1).eq(supported)];pool=frame[frame.isolated&frame.library_UMI.gt(0)].sort_values('atlas_pas_id')
        chosen=pool.iloc[np.sort(rng.choice(len(pool),min(1000,len(pool)),replace=False))].assign(category=category);selections.append(chosen)
        ledger.append(dict(category=category,catalog_sites=len(frame),isolated_sites=int(frame.isolated.sum()),eligible_sites=len(pool),requested=1000,selected=len(chosen)))
    selected=pd.concat(selections,ignore_index=True);ctx.table(selected,folder+'/frozen_read_profile_sites.tsv.gz');ctx.table(ledger,folder+'/read_profile_selection.tsv')
    binding=ctx.read(library.members_path);assert digest(library.members_path)==library.members_sha256
    valid=set(binding.barcode.str.split('-').str[0]);assert len(valid)==len(binding)
    bam_path=Path(library.bam_path);stat=bam_path.stat();assert stat.st_size==library.bam_bytes and stat.st_mtime_ns==library.bam_mtime_ns
    assert digest(ctx.register(library.bai_path))==library.bai_sha256
    curves=[];sites=[];edges=np.arange(-402.5,403,5);centers=(edges[1:]+edges[:-1])/2
    with pysam.AlignmentFile(str(bam_path),'rb') as bam:
        refs={canon(c):c for c in bam.references}
        for category in ['supported50','unsupported50']:
            q=selected[selected.category.eq(category)];raw=np.zeros(161,np.int64);normalized=np.zeros(161);zero=0;missing=0
            for i,row in enumerate(q.itertuples()):
                chrom=refs.get(canon(row.chrom));cut=int(row.cut_site)-int(row.strand=='+');ends=[];status='valid'
                if chrom is None or cut<400 or cut+400>=bam.get_reference_length(chrom):missing+=1;status='unavailable_window'
                else:
                    for read in bam.fetch(chrom,cut-400,cut+401):
                        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary or ('-' if read.is_reverse else '+')!=row.strand:continue
                        if not read.has_tag('CB') or not read.has_tag('UB') or str(read.get_tag('CB')).split('-')[0] not in valid:continue
                        endpoint=read.reference_start if row.strand=='-' else read.reference_end-1
                        relative=(endpoint-cut)*(1 if row.strand=='+' else -1)
                        if -400<=relative<=400:ends.append(relative)
                hist=np.histogram(ends,bins=edges)[0];raw+=hist
                if ends:normalized+=hist/len(ends)
                elif status=='valid':zero+=1;status='zero_accepted_ends'
                sites.append(dict(category=category,atlas_pas_id=row.atlas_pas_id,status=status,accepted_read_ends=len(ends),
                    upstream_300_to_20=sum(-300<=v<=20 for v in ends),near_cut_50=sum(abs(v)<=50 for v in ends)))
                if (i+1)%200==0:print('READ_PROFILES',folder,category,i+1,flush=True)
            curves.extend(dict(category=category,position_nt=p,raw_read_ends=int(n),normalized_site_sum=v,selected_sites=len(q),
                zero_windows=zero,unavailable_windows=missing,mean_site_fraction=v/len(q) if len(q) else np.nan) for p,n,v in zip(centers,raw,normalized))
    final=bam_path.stat();assert (stat.st_size,stat.st_mtime_ns)==(final.st_size,final.st_mtime_ns)
    ctx.table(curves,folder+'/read_profile_curves.tsv.gz');ctx.table(sites,folder+'/read_profile_site_counts.tsv.gz',columns=['category','atlas_pas_id','status','accepted_read_ends'])
    ctx.table([dict(library_key=library.library_key,bam_path=str(bam_path),bam_bytes=stat.st_size,bam_mtime_ns=stat.st_mtime_ns,accepted_barcodes=len(valid),
        counting='accepted reads, not UMI-deduplicated molecules',coordinate='read-end base: interbase cut minus 1 on plus, plus 0 on minus')],folder+'/read_profile_audit.tsv')


def run(ctx):
    import pysam
    assert (ctx.native/'annotation.complete.json').exists()
    manifest=ctx.read(ctx.resource('native_libraries'));summaries=[];strata=[];filter_rows=[];libraries=[]
    for sid in STUDIES:
        folder='WOT_GSE115943' if sid=='GSE115943' else sid;root=ctx.native/'studies'/folder
        f=ctx.read(root/'features.tsv.gz',low_memory=False);sp=f.species.iloc[0]
        counts=sparse.load_npz(ctx.register(root/'metacell_by_pas_counts.npz'));np.testing.assert_array_equal(np.asarray(counts.sum(0)).ravel(),f.current_target_umi)
        mapping=ctx.read(root/'native_peak_to_atlas.tsv.gz');audit=ctx.read(root/'native_library_audit.tsv')
        assert audit.integer_count_conservation.all() and int(audit.integer_molecules.sum())==int(counts.sum())
        chosen=audit.sort_values(['integer_molecules','library_key']).iloc[(len(audit)-1)//2]
        lib=manifest[manifest.library_key.eq(chosen.library_key)];assert len(lib)==1;library=next(lib.itertuples())
        libraries.append(dict(study_id=sid,species=sp,library_key=library.library_key,selection='lower median full native molecule depth; frozen before sequence/read inspection'))
        mis=mapping.assign(flag=truth(mapping.misprime)).groupby('atlas_column_0based').flag.agg(['any','all'])
        f['native_misprime_any']=f.atlas_column_0based.map(mis['any']);f['native_misprime_all']=f.atlas_column_0based.map(mis['all'])
        assert f.native_misprime_any.notna().all()
        for layer,mask in [('all',np.ones(len(f),bool)),('primary',truth(f.keep_for_apa_primary)),('supported50',truth(f.keep_for_apa_primary)&f.matched_source_count_50nt.ge(1)),
                           ('sensitivity',truth(f.keep_for_sensitivity)),('annotation_excluded',~truth(f.keep_for_apa_primary))]:
            part=f[mask];base=dict(study_id=sid,species=sp,layer=layer,sites=len(part),UMI=int(part.current_target_umi.sum()),
                zero_count_sites=int(part.current_target_umi.eq(0).sum()),genes=part.primary_gene_id.nunique())
            for window in [25,50,100]:
                hit=part[f'matched_source_count_{window}nt'].ge(1);base[f'supported{window}_sites']=int(hit.sum());base[f'supported{window}_UMI']=int(part.loc[hit,'current_target_umi'].sum())
            summaries.append(base)
            for field in ['confidence_tier','primary_region_class','primary_assignment_status','internal_priming_flag','multi_gene_ambiguous','native_misprime_any','native_misprime_all']:
                for value,q in part.groupby(field,dropna=False):strata.append(dict(study_id=sid,species=sp,layer=layer,field=field,value=value,sites=len(q),UMI=int(q.current_target_umi.sum())))
        for field in ['tier_reason','failure_reason','neighbor_gene_exclusion_reason']:
            for reason,q in f.groupby(field,dropna=False):filter_rows.append(dict(study_id=sid,species=sp,field=field,reason=reason,sites=len(q),UMI=int(q.current_target_umi.sum())))
        ctx.table(f[['atlas_pas_id','primary_gene_id','confidence_tier','primary_region_class','keep_for_apa_primary','current_target_umi','native_misprime_any','native_misprime_all']+
                    [f'matched_source_count_{w}nt' for w in [25,50,100]]],folder+'/PAS_QC.tsv.gz')
        fasta_path=ctx.resource(sp+'_genome');ctx.register(str(fasta_path)+'.fai')
        with pysam.FastaFile(str(fasta_path)) as fasta:sequence_evidence(ctx,folder,f,mapping,library,sp,fasta)
        print('SEQUENCE_COMPLETE',folder,flush=True);read_profiles(ctx,folder,f,mapping,library)
        print('PAS_QC_COMPLETE',folder,flush=True)
    ctx.table(summaries,'catalog_summary.tsv');ctx.table(strata,'support_and_annotation_strata.tsv');ctx.table(filter_rows,'filter_reason_audit.tsv')
    ctx.table(libraries,'sequence_read_library_selection.tsv')
    ctx.finish(dict(all_seven_catalogues=True,site_and_molecule_denominators_separate=True,native_integer_counts_conserved=True,
        median_library_frozen_before_read_scan=True,full_primary_sequence_profiles=True,sequence_background_reselected=True,
        read_profile_sites_reselected=True,zero_and_unavailable_read_windows_separate=True,native_misprime_is_not_ASNR=True),
        coordinate='0-based interbase input; curves centered on transcript retained endpoint nucleotide',
        background_reference_exclusion='nearest reference source cuts represented in this primary catalogue, not full reference truth set',
        rejected_category='native misprime flag and unified annotation exclusion are reported separately; no ASNR-filtered native cohort exists')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'pas_qc',a.out))
