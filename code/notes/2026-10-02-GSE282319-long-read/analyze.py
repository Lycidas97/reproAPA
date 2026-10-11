"""Audit three pre-existing manuscript examples against public long-read models.

All analysis coordinates are zero-based half-open. Only gap-free, uniquely
mapped chain blocks with exact reciprocal mappings are accepted. TPM is
descriptive detection evidence, never a read count or a replicate-level test.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv
import gzip
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RAW = HERE / 'raw'
REFERENCE = Path('/Users/lycidas/work/科研/projects/2026-09-13_human-mouse-ortholog-mapping/sources/web/ensembl116/genes/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz')
GENES = {'Rad23b': 'ENSMUSG00000028426', 'Hmgn2': 'ENSMUSG00000003038',
         'Eif4ebp2': 'ENSMUSG00000020091'}


def write(name, rows):
    assert rows, name
    with (HERE / name).open('w') as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def chains(path):
    blocks = defaultdict(list)
    with gzip.open(path, 'rt') as f:
        for line in f:
            x = line.split()
            if not x:
                continue
            if x[0] == 'chain':
                _, score, tc, ts, strand, t, te, qc, qs, orient, q, qe, cid = x
                assert strand == '+'
                t, q, qs = int(t), int(q), int(qs)
                continue
            n = int(x[0])
            blocks[tc].append((t, t + n, qc, q, orient, qs, cid))
            if len(x) == 3:
                t += n + int(x[1])
                q += n + int(x[2])
    return blocks


def mapped(blocks, chrom, lo, hi):
    found = []
    for a, b, dest, q, orient, size, cid in blocks[chrom]:
        if a <= lo < hi <= b:
            u, v = q + lo - a, q + hi - a
            if orient == '-':
                u, v = size - v, size - u
            found.append((dest, u, v, orient, cid))
    assert len(found) == 1, ('Ambiguous or gapped interval', chrom, lo, hi, found)
    return found[0]


def lift(forward, reverse, chrom, lo, hi):
    dest, u, v, orient, cid = mapped(forward, chrom, lo, hi)
    back = mapped(reverse, dest, u, v)
    assert back[:3] == (chrom, lo, hi) and back[3] == orient
    assert orient == '+', 'These loci must retain genomic orientation'
    return dest, u, v, cid


def introns(exons):
    return tuple((a[1], b[0]) for a, b in zip(exons, exons[1:]))


def exon_models(path, wanted, attribute):
    models = {}
    kept = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#') or not any('"' + t + '"' in line for t in wanted):
                continue
            x = line.rstrip().split('\t')
            if x[2] != 'exon':
                continue
            a = dict(re.findall(r'(\w+) "([^"]+)"', x[8]))
            if a.get(attribute) not in wanted:
                continue
            tid = a['transcript_id']
            m = models.setdefault(tid, dict(gene=a['gene_id'], chrom='chr' + x[0].removeprefix('chr'),
                                            strand=x[6], exons=[]))
            assert (m['chrom'], m['strand']) == ('chr' + x[0].removeprefix('chr'), x[6])
            m['exons'].append((int(x[3]) - 1, int(x[4])))
            kept.append(line)
    for m in models.values():
        m['exons'].sort()
        assert all(a[1] <= b[0] for a, b in zip(m['exons'], m['exons'][1:]))
    return models, kept


def main():
    fw = chains(RAW / 'mm39ToMm10.over.chain.gz')
    rv = chains(RAW / 'mm10ToMm39.over.chain.gz')
    sources = [ROOT / f'notes/2026-10-01-mpro-biology/annotation/mouse/{r}/models.json' for r in ['CIP', 'WOT']]
    study_models = [json.loads(p.read_text()) for p in sources]
    families = {g: study_models[0][gid]['family'] for g, gid in GENES.items()}
    assert all({k: v for k, v in families[g].items() if k != 'pas_ids'} ==
               {k: v for k, v in study_models[1][gid]['family'].items() if k != 'pas_ids'}
               for g, gid in GENES.items())
    features = [dict(gene='Rad23b', label='Rad23b_shared415', chrom='chr4', strand='+', start0=55390761, end0=55391176),
                dict(gene='Rad23b', label='Rad23b_miR290302_seed', chrom='chr4', strand='+', start0=55390802, end0=55390809)]
    displayed = ROOT / 'notes/2026-09-29-core-sensitivity/figure5-coordinate-audit/displayed_site_sensitivity.tsv'
    selected = {}
    for r in csv.DictReader(displayed.open(), delimiter='\t'):
        selected[r['site_label']] = dict(gene=r['gene_symbol_native'], label=r['site_label'], chrom=r['chrom'],
                                         strand=r['strand'], start0=int(r['start0']), end0=int(r['end0']))
    features.extend(selected.values())
    recurrent = ROOT / 'table/Table.S19.figure5-union/recurrent_region_pairs.tsv'
    for r in csv.DictReader(recurrent.open(), delimiter='\t'):
        if r['gene_symbol_a'] in ['Hmgn2', 'Eif4ebp2']:
            key = (r['gene_symbol_a'], int(r['shared_start0']), int(r['shared_end0']))
            selected[key] = dict(gene=key[0], label=f'{key[0]}_shared_{key[1]}_{key[2]}', chrom=r['chrom'],
                                 strand=r['strand'], start0=key[1], end0=key[2])
    features.extend(v for k, v in selected.items() if isinstance(k, tuple))
    lift_rows = []
    for f in features:
        c, lo, hi, cid = lift(fw, rv, f['chrom'], f['start0'], f['end0'])
        f.update(mm10_chrom=c, mm10_start0=lo, mm10_end0=hi, chain_id=cid)
        lift_rows.append(dict(kind='feature', **f, reciprocal_exact=True))
    refids = set()
    for gene, fam in families.items():
        for m in fam['models']:
            refids.update(m['transcript_ids'].split(';'))
        # exon_anchor is a GTF 1-based boundary in the frozen source models.
        anchor = fam['exon_anchor'] - (1 if fam['strand'] == '+' else 0)
        c, lo, hi, cid = lift(fw, rv, fam['chrom'], anchor, anchor + 1)
        fam['anchor_mm10'] = lo
        lift_rows.append(dict(kind='anchor', gene=gene, label='terminal_splice_boundary', chrom=fam['chrom'],
                              strand=fam['strand'], start0=anchor, end0=anchor + 1, mm10_chrom=c,
                              mm10_start0=lo, mm10_end0=hi, chain_id=cid, reciprocal_exact=True))
    write('coordinate_liftover.tsv', lift_rows)
    refs, ref_lines = exon_models(REFERENCE, refids, 'transcript_id')
    # The secondary Ensembl 116 annotation need not retain every frozen 2024-A
    # transcript. Preserve missing IDs explicitly; each frozen model needs a
    # retained representative before any whole-chain comparison is permitted.
    missing_reference_ids = sorted(refids - set(refs))
    assert all(set(m['transcript_ids'].split(';')) & set(refs)
               for fam in families.values() for m in fam['models'])
    refmapped = {}
    for tid, m in refs.items():
        exons = [lift(fw, rv, m['chrom'], a, b)[1:3] for a, b in m['exons']]
        refmapped[tid] = dict(**m, mm10_exons=exons)
        gene = next(g for g, gid in GENES.items() if gid == m['gene'])
        term = exons[-1] if m['strand'] == '+' else exons[0]
        assert (term[0] if m['strand'] == '+' else term[1]) == families[gene]['anchor_mm10']
    (HERE / 'reference_exons.mm39.gtf').write_text(''.join(ref_lines))
    (HERE / 'reference_models.json').write_text(json.dumps(refmapped, indent=2) + '\n')
    txs, tx_lines = exon_models(RAW / 'GSE282319_Repro_filter_tst.gtf.gz', set(GENES), 'gene_id')
    (HERE / 'candidate_exons.mm10.gtf').write_text(''.join(tx_lines))
    qcs = {}
    gene_qcs = []
    qcids = {t.split('.', 1)[1] for t in txs}
    with gzip.open(RAW / 'GSE282319_Repro_unfiltered_tst_classification.txt.gz', 'rt') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['associated_gene'] in GENES.values() or row['isoform'] in qcids:
                gene_qcs.append(row)
                qcs[row['isoform']] = row
    write('unfiltered_candidate_classification.tsv', gene_qcs)
    abundance = {}
    with gzip.open(RAW / 'GSE282319_Isoquant.transcript_name_id.tpm.txt.gz', 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        samples = reader.fieldnames[1:]
        assert len(samples) == len(set(samples))
        for row in reader:
            if row['transcript_name_id'] in txs:
                assert row['transcript_name_id'] not in abundance
                abundance[row['transcript_name_id']] = {s: float(row[s]) for s in samples}
    assert set(abundance) == set(txs), ('Missing TPM rows', set(txs) - set(abundance))
    txrows, audit, pairrows = [], [], []
    for tid, m in txs.items():
        gene = m['gene']
        exons = m['exons']
        term = exons[-1] if m['strand'] == '+' else exons[0]
        boundary = term[0] if m['strand'] == '+' else term[1]
        frame = boundary == families[gene]['anchor_mm10']
        q = qcs.get(tid.split('.', 1)[1], {})
        assert q, ('Filtered ID absent from raw classification', tid)
        assert q['associated_gene'] == GENES[gene]
        assert int(q['exons']) == len(exons)
        mapped_exons = [lift(rv, fw, m['chrom'], a, b)[1:3] for a, b in exons]
        samechain = [t for t, r in refmapped.items() if r['gene'] == GENES[gene] and introns(r['mm10_exons']) == introns(exons)]
        m.update(terminal=term, frame=frame, qc=q, samechain=samechain)
        values = abundance[tid]
        row = dict(gene=gene, transcript_id=tid, chrom=m['chrom'], strand=m['strand'], n_exons=len(exons),
                   exons_mm10=json.dumps(exons), exons_mm39=json.dumps(mapped_exons),
                   terminal_start0_mm10=term[0], terminal_end0_mm10=term[1], terminal_anchor_exact=frame,
                   reference_intron_chain_matches=';'.join(samechain),
                   structural_category=q['structural_category'], associated_transcript=q['associated_transcript'],
                   diff_to_TSS=q['diff_to_TSS'], diff_to_TTS=q['diff_to_TTS'], RTS_stage=q['RTS_stage'],
                   all_canonical=q['all_canonical'], polyA_motif=q['polyA_motif'], polyA_dist=q['polyA_dist'],
                   polyA_motif_found=q['polyA_motif_found'], perc_A_downstream_TTS=q['perc_A_downstream_TTS'],
                   n_samples_TPM_positive=sum(v > 0 for v in values.values()),
                   n_samples_TPM_ge1=sum(v >= 1 for v in values.values()), max_TPM=max(values.values()), sum_TPM=sum(values.values()))
        txrows.append(row)
        for feature in [f for f in features if f['gene'] == gene]:
            lo, hi = feature['mm10_start0'], feature['mm10_end0']
            overlap = sum(max(0, min(b, hi) - max(a, lo)) for a, b in exons)
            status = 'included' if term[0] <= lo and hi <= term[1] else ('excluded' if overlap == 0 else 'partial_or_spliced')
            audit.append(dict(gene=gene, transcript_id=tid, feature=feature['label'], status=status,
                              exonic_overlap_nt=overlap, feature_length_nt=hi - lo, terminal_anchor_exact=frame,
                              reference_intron_chain_matches=';'.join(samechain), n_samples_TPM_positive=row['n_samples_TPM_positive'],
                              n_samples_TPM_ge1=row['n_samples_TPM_ge1'], max_TPM=row['max_TPM']))
    write('transcript_audit.tsv', txrows)
    write('feature_retention.tsv', audit)
    write('candidate_TPM.tsv', [dict(gene=txs[t]['gene'], transcript_id=t, **abundance[t]) for t in sorted(txs)])
    qc_ok = lambda m: m['qc']['RTS_stage'] == 'FALSE' and m['qc']['all_canonical'] == 'canonical' and m['qc']['polyA_motif_found'] == 'TRUE'
    for feature in features:
        lo, hi = feature['mm10_start0'], feature['mm10_end0']
        eligible = [(t, m) for t, m in txs.items() if m['gene'] == feature['gene'] and m['frame'] and qc_ok(m)]
        for t, m in eligible:
            if not (m['terminal'][0] <= lo and hi <= m['terminal'][1]):
                continue
            for u, n in eligible:
                if introns(m['exons']) != introns(n['exons']):
                    continue
                overlap = sum(max(0, min(b, hi) - max(a, lo)) for a, b in n['exons'])
                if overlap == hi - lo:
                    continue
                proximal = n['terminal'][1] <= lo + 1 if n['strand'] == '+' else n['terminal'][0] >= hi - 1
                if not proximal:
                    continue
                pairrows.append(dict(gene=feature['gene'], feature=feature['label'], included_transcript=t,
                                     proximal_transcript=u, same_complete_intron_chain=True, same_terminal_anchor=True,
                                     proximal_overlap_nt=overlap, strict_full_exclusion=overlap == 0,
                                     inclusion_TPM_positive_samples=sum(v > 0 for v in abundance[t].values()),
                                     proximal_TPM_positive_samples=sum(v > 0 for v in abundance[u].values()),
                                     both_TPM_positive_samples=sum(abundance[t][s] > 0 and abundance[u][s] > 0 for s in samples),
                                     same_chain_as_reference=bool(m['samechain']), included_sum_TPM=sum(abundance[t].values()),
                                     proximal_sum_TPM=sum(abundance[u].values())))
    write('compatible_isoform_pairs.tsv', pairrows)
    # Keep the entire GEO sample metadata; never infer R9/R10 biological independence.
    metadata, current = {}, None
    with gzip.open(RAW / 'GSE282319_family.soft.gz', 'rt') as f:
        for line in f:
            if line.startswith('^SAMPLE = '):
                current = line.strip().split(' = ')[1]
                metadata[current] = defaultdict(list)
            elif current and line.startswith('!Sample_'):
                k, _, v = line.strip().partition(' = ')
                metadata[current][k.removeprefix('!Sample_')].append(v)
    (HERE / 'sample_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    titles = {v['title'][0]: (gsm, v) for gsm, v in metadata.items()}
    aliases = {'d18_OKS_ctrl_pc10': 'd18_OKS_ctrl_PCS_R10', 'd18_2iL_OKMS_ctrl_pc10': 'd18_2iL_OKMS_ctrl_PCS_R10',
               'd18_2iL_OKS_ctrl_pc10': 'd18_2iL_OKS_ctrl_PCS_R10'}
    meta_rows = []
    for s in samples:
        title = aliases.get(s, s)
        if title not in titles and title + '_PCS_R10' in titles:
            title += '_PCS_R10'
        assert title in titles, ('Unmatched TPM sample', s)
        gsm, meta = titles[title]
        chars = dict(x.split(': ', 1) for x in meta['characteristics_ch1'])
        meta_rows.append(dict(TPM_column=s, GSM=gsm, title=title, platform='R9' if 'R9' in title else 'R10',
                              cell_line=chars['cell line'], cell_type=chars['cell type'], genotype=chars['genotype'],
                              treatment=chars['treatment'], biological_replicate_independence='not established from title/platform',
                              SRA=';'.join(meta.get('relation', []))))
    write('sample_inventory.tsv', meta_rows)
    assert len(set(r['GSM'] for r in meta_rows)) == len(samples) == len(metadata)
    summaries = []
    for f in features:
        rows = [r for r in audit if r['feature'] == f['label']]
        pairs = [r for r in pairrows if r['feature'] == f['label']]
        summaries.append(dict(gene=f['gene'], feature=f['label'], mm10_start0=f['mm10_start0'], mm10_end0=f['mm10_end0'],
                              n_filtered_gene_isoforms=len(rows), n_anchor_included=sum(r['terminal_anchor_exact'] and r['status'] == 'included' for r in rows),
                              n_anchor_excluded=sum(r['terminal_anchor_exact'] and r['status'] == 'excluded' for r in rows),
                              n_same_chain_strict_pairs=sum(r['strict_full_exclusion'] for r in pairs),
                              n_same_chain_boundary1nt_pairs=sum(not r['strict_full_exclusion'] for r in pairs),
                              max_same_sample_co_detection=max((r['both_TPM_positive_samples'] for r in pairs), default=0)))
    write('structural_summary.tsv', summaries)
    files = sources + [displayed, recurrent, REFERENCE] + list(RAW.iterdir())
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.is_file()}
    (HERE / 'input_sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
    verification = dict(status='PASS', selection='Three pre-existing manuscript examples; no outcome-based locus selection',
                        filtered_isoforms=dict(Counter(m['gene'] for m in txs.values())), unique_TPM_samples=len(samples),
                        all_filtered_candidate_ids_join_raw_QC=True, all_TPM_columns_join_unique_GSM=True,
                        all_reference_and_long_read_exons_gapfree_unique_reciprocal_lift=True,
                        reference_terminal_boundaries_match_frozen_models=True,
                        secondary_reference='Ensembl 116 mm39; frozen source geometry remains primary',
                        secondary_reference_missing_frozen_transcript_ids=missing_reference_ids,
                        coordinate_system='0-based half-open; frozen exon_anchor converted from GTF 1-based boundary',
                        primary_pairs='Exact terminal boundary and identical complete intron chain; canonical junctions, RTS false, polyA motif true',
                        boundary_sensitivity='Record <=1 nt overlap separately from strict full exclusion',
                        raw_read_validation=False, differential_usage_test=False, manuscript_modified=False)
    (HERE / 'verification.json').write_text(json.dumps(verification, indent=2) + '\n')
    print(json.dumps(verification, indent=2))


if __name__ == '__main__':
    main()
