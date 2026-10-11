"""Re-intersect new terminal features with frozen external structures and AGO2 tracks."""
import argparse
from collections import defaultdict
import csv
import gzip
import json
import re
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from downstream_runtime import Context

LONG='notes/2026-10-02-GSE282319-long-read/analyze.py'
AGO='notes/2026-10-02-rad23b-ago2-gse25310/check_overlap.py'
GENES={'Rad23b':'ENSMUSG00000028426','Hmgn2':'ENSMUSG00000003038','Eif4ebp2':'ENSMUSG00000020091'}
TRACKS=[('CLIPSeq_Ago2_Leung2011a_mm9','KO_merged','merged'),('CLIPSeq_Ago2_Leung2011b_mm9','WT_merged','merged'),
        ('HITSCLIP_Ago2Leung2011a_mm9','WT1A','WT1'),('HITSCLIP_Ago2Leung2011b_mm9','WT1B','WT1'),('HITSCLIP_Ago2Leung2011c_mm9','WT2','WT2')]


def numerics(ctx):
    return ctx.functions(LONG,['chains','mapped','lift','introns','exon_models'],dict(defaultdict=defaultdict,gzip=gzip,re=re))


def mapped_interval(math,forward,reverse,chrom,start,end):
    try:
        c,a,b,chain=math.lift(forward,reverse,chrom,int(start),int(end))
        assert b-a==end-start
        return dict(status='unique_gapfree_reciprocal',chrom=c,start0=a,end0=b,chain_id=chain)
    except AssertionError as error:
        return dict(status='unmappable_or_nonreciprocal',reason=str(error))


def retention(exons,strand,start,end):
    terminal=exons[-1] if strand=='+' else exons[0]
    overlap=sum(max(0,min(b,end)-max(a,start)) for a,b in exons)
    status='included' if terminal[0]<=start<end<=terminal[1] else 'excluded' if overlap==0 else 'partial_or_spliced'
    return status,overlap


def published_longread_rows(path,transcripts):
    """Select current transcript IDs from the complete published classification."""
    import openpyxl
    wanted=set(transcripts);records=[]
    workbook=openpyxl.load_workbook(path,read_only=True,data_only=True)
    try:
        rows=workbook['Classification-file'].iter_rows(values_only=True)
        header=list(next(rows));key=header.index('Transcript name ID')
        for row in rows:
            if row[key] in wanted:records.append(dict(zip(header,row)))
    finally:workbook.close()
    frame=pd.DataFrame(records,columns=header).set_index('Transcript name ID')
    assert frame.index.is_unique,'Duplicate published transcript classification'
    assert set(frame.index)==wanted,'Missing published transcript classification: '+str(sorted(wanted-set(frame.index)))
    return frame


def longread(ctx):
    math=numerics(ctx);forward=math.chains(ctx.resource('chain_mm39_mm10'));reverse=math.chains(ctx.resource('chain_mm10_mm39'))
    features=ctx.read(ctx.analysis/'examples/external_query_features.tsv');lifted=[];mapping=[]
    named=features[['gene','gene_id']].drop_duplicates();assert named.gene.is_unique
    genes=dict(GENES);genes.update(dict(named.itertuples(index=False,name=None)))
    for row in features.to_dict('records'):
        m=mapped_interval(math,forward,reverse,row['chrom'],row['start0'],row['end0'])
        mapping.append(dict(row,**{'mapping_'+k:v for k,v in m.items()}))
        if m['status']=='unique_gapfree_reciprocal':lifted.append(dict(row,mm10_chrom=m['chrom'],mm10_start0=m['start0'],mm10_end0=m['end0']))
    families={};family_audit=[];wanted=set()
    for route in ['CIP','WOT']:
        models=json.loads(ctx.register(ctx.analysis/f'intervals/annotation/mouse/{route}/models.json').read_text())
        for gene,gid in genes.items():
            if gid not in models:
                family_audit.append(dict(gene=gene,gene_id=gid,route=route,status='no_new_terminal_family'));continue
            family=models[gid]['family'];anchor=family['exon_anchor']-int(family['strand']=='+')
            m=mapped_interval(math,forward,reverse,family['chrom'],anchor,anchor+1)
            ids={t for model in family['models'] for t in model['transcript_ids'].split(';')};wanted|=ids
            family_audit.append(dict(gene=gene,gene_id=gid,route=route,status=m['status'],anchor_mm39=anchor,
                anchor_mm10=m.get('start0',np.nan),reference_transcripts=';'.join(sorted(ids))))
            if m['status']=='unique_gapfree_reciprocal':families[gene,route]=dict(family,anchor_mm10=m['start0'],mm10_chrom=m['chrom'],reference_ids=ids)
    references,ref_lines=math.exon_models(ctx.resource('mouse_ensembl_gtf'),wanted,'transcript_id')
    mapped_refs={};reference_audit=[]
    for tid in sorted(wanted):
        if tid not in references:
            reference_audit.append(dict(transcript_id=tid,status='missing_secondary_reference'));continue
        r=references[tid];parts=[mapped_interval(math,forward,reverse,r['chrom'],a,b) for a,b in r['exons']]
        good=all(p['status']=='unique_gapfree_reciprocal' for p in parts)
        reference_audit.append(dict(transcript_id=tid,gene_id=r['gene'],status='mapped' if good else 'exon_mapping_failed'))
        if good:mapped_refs[tid]=dict(r,exons_mm10=[(p['start0'],p['end0']) for p in parts])
    txs,tx_lines=math.exon_models(ctx.resource('longread_gtf'),set(genes),'gene_id')
    qcs={};qcids={t.split('.',1)[1] for t in txs}
    with gzip.open(ctx.resource('longread_classification'),'rt') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            if row['isoform'] in qcids:
                assert row['isoform'] not in qcs;qcs[row['isoform']]=row
    abundance={}
    with gzip.open(ctx.resource('longread_tpm'),'rt') as handle:
        reader=csv.DictReader(handle,delimiter='\t');samples=reader.fieldnames[1:];assert len(samples)==len(set(samples))
        for row in reader:
            if row['transcript_name_id'] in txs:
                assert row['transcript_name_id'] not in abundance
                abundance[row['transcript_name_id']]={s:float(row[s]) for s in samples}
    assert set(abundance)==set(txs)
    inventory=ctx.read(ctx.resource('longread_sample_inventory'));assert set(inventory.TPM_column)==set(samples) and inventory.GSM.is_unique
    published=published_longread_rows(ctx.resource('longread_raw_supplement2'),txs)
    ctx.resource('longread_raw_sample_metadata');ctx.resource('longread_raw_supplement1');ctx.resource('longread_raw_supplement2')
    txrows=[];measurements=[];compatibility=[];pairs=[]
    for tid,m in txs.items():
        gene=m['gene'];qc=qcs[tid.split('.',1)[1]];assert qc['associated_gene']==genes[gene] and int(qc['exons'])==len(m['exons'])
        pub=published.loc[tid];assert str(pub['Filtered out for low expression']).lower()=='false'
        assert str(pub['Supported either by junction or TSS']).lower()=='true' and int(pub['Number of exons'])==len(m['exons'])
        assert str(pub['Structural category'])==qc['structural_category'] and str(pub['All canonical'])==qc['all_canonical']
        terminal=m['exons'][-1] if m['strand']=='+' else m['exons'][0];boundary=terminal[0] if m['strand']=='+' else terminal[1]
        parts=[mapped_interval(math,reverse,forward,m['chrom'],a,b) for a,b in m['exons']]
        mapped_ok=all(p['status']=='unique_gapfree_reciprocal' for p in parts)
        qc_ok=qc['RTS_stage']=='FALSE' and qc['all_canonical']=='canonical' and qc['polyA_motif_found']=='TRUE'
        values=abundance[tid];m.update(terminal=terminal,boundary=boundary,qc_ok=qc_ok,mapped_ok=mapped_ok)
        txrows.append(dict(gene=gene,transcript_id=tid,chrom=m['chrom'],strand=m['strand'],exons_mm10=json.dumps(m['exons']),
            exons_mm39=json.dumps([(p['start0'],p['end0']) for p in parts]) if mapped_ok else '',reciprocal_exons=mapped_ok,
            terminal_start0_mm10=terminal[0],terminal_end0_mm10=terminal[1],n_exons=len(m['exons']),qc_ok=qc_ok,
            max_TPM=max(values.values()),n_samples_TPM_positive=sum(v>0 for v in values.values()),n_samples_TPM_ge1=sum(v>=1 for v in values.values()),
            **{k:qc[k] for k in ['structural_category','RTS_stage','all_canonical','polyA_motif','polyA_dist','polyA_motif_found']}))
        for route in ['CIP','WOT']:
            fam=families.get((gene,route));frame=bool(fam and boundary==fam['anchor_mm10'] and m['chrom']==fam['mm10_chrom'] and m['strand']==fam['strand'])
            same=[t for t,r in mapped_refs.items() if fam and t in fam['reference_ids'] and math.introns(r['exons_mm10'])==math.introns(m['exons'])]
            compatibility.append(dict(gene=gene,transcript_id=tid,route=route,terminal_anchor_exact=frame,reference_intron_chain_matches=';'.join(same),
                secondary_reference_models_complete=bool(fam and all(set(z['transcript_ids'].split(';'))&set(mapped_refs) for z in fam['models']))))
        for f in (z for z in lifted if z['gene_id']==genes[gene]):
            status,overlap=retention(m['exons'],m['strand'],f['mm10_start0'],f['mm10_end0'])
            measurements.append(dict(gene=gene,transcript_id=tid,feature_id=f['feature_id'],status=status,exonic_overlap_nt=overlap,feature_length_nt=f['mm10_end0']-f['mm10_start0']))
    for feature in lifted:
        gene=next(g for g,gid in genes.items() if gid==feature['gene_id']);lo,hi=feature['mm10_start0'],feature['mm10_end0']
        fams=[families.get((gene,r)) for r in ['CIP','WOT']]
        if any(f is None for f in fams) or len({(f['mm10_chrom'],f['strand'],f['anchor_mm10']) for f in fams})!=1:continue
        compatible={tid:m for tid,m in txs.items() if m['gene']==gene and m['qc_ok'] and m['mapped_ok'] and m['boundary']==fams[0]['anchor_mm10']
                    and m['chrom']==fams[0]['mm10_chrom'] and m['strand']==fams[0]['strand']}
        for tid,m in compatible.items():
            if retention(m['exons'],m['strand'],lo,hi)[0]!='included':continue
            for uid,n in compatible.items():
                status,overlap=retention(n['exons'],n['strand'],lo,hi)
                proximal=n['terminal'][1]<=lo+1 if n['strand']=='+' else n['terminal'][0]>=hi-1
                if math.introns(m['exons'])!=math.introns(n['exons']) or overlap==hi-lo or not proximal:continue
                row=dict(gene=gene,feature_id=feature['feature_id'],included_transcript=tid,proximal_transcript=uid,
                    same_complete_intron_chain=True,same_terminal_anchor_both_new_families=True,proximal_overlap_nt=overlap,strict_full_exclusion=overlap==0)
                for platform in ['all','R9','R10']:
                    ss=inventory.loc[inventory.platform.eq(platform) if platform!='all' else np.ones(len(inventory),bool),'TPM_column']
                    for threshold in [0,1]:row[f'both_TPM_{"positive" if threshold==0 else "ge1"}_samples_{platform}']=sum(
                        abundance[tid][s]>0 and abundance[uid][s]>0 if threshold==0 else abundance[tid][s]>=1 and abundance[uid][s]>=1 for s in ss)
                pairs.append(row)
    summary=[]
    for f in mapping:
        rs=[r for r in measurements if r['feature_id']==f['feature_id']];ps=[r for r in pairs if r['feature_id']==f['feature_id']]
        summary.append(dict(feature_id=f['feature_id'],gene=f['gene'],kind=f['kind'],mapping_status=f['mapping_status'],
            isoforms=len(rs),included=sum(r['status']=='included' for r in rs),excluded=sum(r['status']=='excluded' for r in rs),
            strict_pairs=sum(r['strict_full_exclusion'] for r in ps),boundary1nt_pairs=sum(not r['strict_full_exclusion'] for r in ps)))
    ctx.table(mapping,'coordinate_liftover.tsv',columns=['feature_id','mapping_status']);ctx.table(family_audit,'family_eligibility.tsv')
    ctx.table(reference_audit,'secondary_reference_audit.tsv',columns=['transcript_id','status']);ctx.table(txrows,'transcript_audit.tsv',columns=['gene','transcript_id','qc_ok','reciprocal_exons'])
    ctx.table(compatibility,'transcript_family_compatibility.tsv',columns=['gene','transcript_id','route','terminal_anchor_exact']);ctx.table(measurements,'feature_retention.tsv',columns=['gene','transcript_id','feature_id','status'])
    ctx.table(pairs,'compatible_isoform_pairs.tsv',columns=['gene','feature_id','included_transcript','proximal_transcript','strict_full_exclusion'])
    ctx.table(summary,'structural_summary.tsv',columns=['feature_id','gene','kind','mapping_status','isoforms','strict_pairs'])
    ctx.table([dict(transcript_id=tid,gene=txs[tid]['gene'],**values) for tid,values in abundance.items()],'candidate_TPM.tsv',columns=['transcript_id','gene']+samples);ctx.table(inventory,'sample_inventory.tsv')
    (ctx.out/'reference_exons.mm39.gtf').write_text(''.join(ref_lines));(ctx.out/'candidate_exons.mm10.gtf').write_text(''.join(tx_lines))
    ctx.table([dict(gene=g,status='new_feature_available' if features.gene_id.eq(gid).any() else 'no_new_recurrent_feature',
        new_features=int(features.gene_id.eq(gid).sum()),filtered_isoforms=sum(m['gene']==g for m in txs.values())) for g,gid in genes.items()],'historical_three_gene_tracking.tsv')
    ctx.finish(dict(new_features_reintersected=True,unique_gapfree_reciprocal_mapping_required=True,missing_family_and_reference_explicit=True,
        exact_intron_chain_and_anchor_required=True,published_QC_and_TPM_samples_joined=True,TPM_used_only_descriptively=True),
        raw_read_validation=False,differential_usage_test=False,biological_replicate_independence='not inferred from platforms or passages')


def checked_mm9_sequence(data,key):
    payload=json.loads(data)
    assert payload['genome']=='mm9' and payload['chrom']==key[0]
    assert int(payload['start'])==key[1] and int(payload['end'])==key[2]
    sequence=payload['dna'].upper();assert len(sequence)==key[2]-key[1]
    return sequence


def ago2(ctx):
    import pysam
    math=numerics(ctx);chains={s:math.chains(ctx.resource('chain_'+s)) for s in ['mm39_mm10','mm10_mm39','mm10_mm9','mm9_mm10']}
    helpers=ctx.functions(AGO,['read_bed','overlap_bp'])
    features=ctx.read(ctx.analysis/'examples/external_query_features.tsv');f=features[features.gene_id.eq(GENES['Rad23b'])]
    f=f[f.kind.eq('shared_region')|(f.kind.eq('site')&f.category.eq('miRNA'))]
    queries=[];mappings=[]
    for r in f.to_dict('records'):
        first=mapped_interval(math,chains['mm39_mm10'],chains['mm10_mm39'],r['chrom'],r['start0'],r['end0'])
        second=mapped_interval(math,chains['mm10_mm9'],chains['mm9_mm10'],first['chrom'],first['start0'],first['end0']) if first['status']=='unique_gapfree_reciprocal' else first
        mappings.append(dict(r,mm10=json.dumps(first),mm9=json.dumps(second),mapping_status=second['status']))
        if second['status']=='unique_gapfree_reciprocal':queries.append(dict(r,chrom=second['chrom'],start0=second['start0'],end0=second['end0']))
    inventory=[];summary=[];hits=[];mapped_hits=[]
    for prefix,label,unit in TRACKS:
        metadata=json.loads(ctx.resource('ago2_'+prefix+'_json').read_text())[0];bed=helpers.read_bed(ctx.resource('ago2_'+prefix+'_bed'))
        assert str(metadata['PMID'])=='21258322' and len(bed)==int(metadata['sites'])
        inventory.append(dict(dataset=prefix,label=label,biological_unit=unit,sites=len(bed),methods=metadata['methods'],description=metadata['description']))
        for q in queries:
            overlap=[(b,helpers.overlap_bp(q,b)) for b in bed if helpers.overlap_bp(q,b)>0]
            summary.append(dict(feature_id=q['feature_id'],dataset=prefix,label=label,biological_unit=unit,overlap_clusters=len(overlap),overlap_bp_sum=sum(n for _,n in overlap)))
            for b,n in overlap:
                hit=dict(feature_id=q['feature_id'],dataset=prefix,label=label,biological_unit=unit,**b,overlap_bp=n)
                hits.append(hit)
                m=mapped_interval(math,chains['mm9_mm10'],chains['mm10_mm9'],b['chrom'],b['start0'],b['end0'])
                out=mapped_interval(math,chains['mm10_mm39'],chains['mm39_mm10'],m['chrom'],m['start0'],m['end0']) if m['status']=='unique_gapfree_reciprocal' else m
                mapped_hits.append(dict(hit,**{'mm39_'+k:v for k,v in out.items()},archived_flanks_may_be_present=True))
    support=[]
    for query in queries:
        counts={r['label']:r['overlap_clusters'] for r in summary if r['feature_id']==query['feature_id']}
        support.append(dict(feature_id=query['feature_id'],WT1_any_technical_library=counts['WT1A']>0 or counts['WT1B']>0,
            WT2=counts['WT2']>0,WT_biological_units_supported=int(counts['WT1A']>0 or counts['WT1B']>0)+int(counts['WT2']>0)))
    sequences=[];original=f.set_index('feature_id');downloaded={};cache={};cache_root=None
    if 'ago2_mm9_sequence_cache' in ctx.config['resources']:
        manifest=ctx.resource('ago2_mm9_sequence_cache');cache_root=manifest.parent
        entries=json.loads(manifest.read_text())['requests']
        cache={(a['chrom'],a['start'],a['end']):a for a in entries}
        assert len(cache)==len(entries)
    fasta_path=ctx.resource('mouse_genome');ctx.register(str(fasta_path)+'.fai')
    with pysam.FastaFile(str(fasta_path)) as fasta:
        for q in queries:
            key=(q['chrom'],q['start0'],q['end0'])
            if key not in downloaded:
                url='https://api.genome.ucsc.edu/getData/sequence?'+urllib.parse.urlencode(dict(genome='mm9',chrom=key[0],start=key[1],end=key[2]))
                if cache_root is not None and key in cache:
                    item=cache[key];assert item['url']==url
                    source=ctx.register(cache_root/item['file']);assert ctx.inputs[str(source)]==item['sha256']
                    data=source.read_bytes()
                else:
                    with urllib.request.urlopen(url,timeout=90) as response:data=response.read()
                oldseq=checked_mm9_sequence(data,key)
                path=ctx.out/f'ucsc_mm9_sequence_{len(downloaded)}.json';path.write_bytes(data);ctx.register(path)
                downloaded[key]=(url,oldseq)
            url,oldseq=downloaded[key];r=original.loc[q['feature_id']]
            chrom=r.chrom if r.chrom in fasta.references else r.chrom.removeprefix('chr')
            seq=fasta.fetch(chrom,int(r.start0),int(r.end0)).upper();assert len(seq)==len(oldseq)==q['end0']-q['start0']
            oriented=seq if r.strand=='+' else seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
            sequences.append(dict(feature_id=q['feature_id'],mm39_positive_strand_sequence=seq,mm9_positive_strand_sequence=oldseq,
                transcript_oriented_sequence=oriented,sequence_identical=seq==oldseq,mm9_sequence_source=url))
    ctx.table(mappings,'query_mapping.tsv',columns=['feature_id','mapping_status']);ctx.table(queries,'queries.mm9.tsv',columns=list(features.columns))
    ctx.table(inventory,'archive_inventory.tsv');ctx.table(summary,'overlap_summary.tsv',columns=['feature_id','dataset','label','biological_unit','overlap_clusters','overlap_bp_sum'])
    ctx.table(hits,'overlapping_clusters.mm9.tsv',columns=['feature_id','dataset','chrom','start0','end0','strand','overlap_bp'])
    ctx.table(mapped_hits,'overlapping_clusters.mm39.tsv',columns=['feature_id','dataset','mm39_status','mm39_start0','mm39_end0'])
    ctx.table(support,'biological_unit_support.tsv',columns=['feature_id','WT1_any_technical_library','WT2','WT_biological_units_supported'])
    ctx.table(sequences,'sequence_identity.tsv',columns=['feature_id','mm39_positive_strand_sequence','mm9_positive_strand_sequence','sequence_identical'])
    ctx.table([dict(gene='Rad23b',new_queries=len(f),mapped_queries=len(queries),status='tested' if queries else 'no_current_mappable_recurrent_feature')],'query_eligibility.tsv')
    ctx.finish(dict(new_interval_and_seed_queries=True,no_fixed_415nt_interval=True,unique_reciprocal_mapping=True,
        same_strand_positive_overlap=True,WT1_technical_libraries_one_biological_unit=True,missing_features_explicit=True),
        sequence_identity_across_assemblies='new query sequence compared and mismatches retained',
        limitations=['Archived cluster flanks may be present; original read bounds are not recovered.',
                     'Merged WT and KO archives lack biological replicate membership.',
                     'No read reprocessing, new peak calling or WT versus KO differential test.'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--node',choices=['longread','ago2'],required=True)
    a=p.parse_args();ctx=Context(a.config,a.node,a.out);(longread if a.node=='longread' else ago2)(ctx)
