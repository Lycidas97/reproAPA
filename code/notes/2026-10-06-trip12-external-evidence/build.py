"""Add a Trip12 evidence overlay; never mutate the sealed ResultBundle."""
import ast
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUNDLE = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main'
ARCHIVE = ROOT / 'notes/2026-10-02-rad23b-ago2-gse25310/inputs'
TRACKS = [('CLIPSeq_Ago2_Leung2011a_mm9', 'KO_merged', 'merged'),
          ('CLIPSeq_Ago2_Leung2011b_mm9', 'WT_merged', 'merged'),
          ('HITSCLIP_Ago2Leung2011a_mm9', 'WT1A', 'WT1'),
          ('HITSCLIP_Ago2Leung2011b_mm9', 'WT1B', 'WT1'),
          ('HITSCLIP_Ago2Leung2011c_mm9', 'WT2', 'WT2')]
INPUTS = {}
MANIFEST = {r['path']: r for r in json.loads((BUNDLE / 'table-manifest.json').read_text())['tables']}


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def checked(path, expected=None):
    digest = sha(path)
    if expected is not None:
        assert digest == expected, str(path)
    INPUTS[str(path)] = digest
    return path


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader.fieldnames), list(reader)


def write(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with (gzip.open(path, 'wt') if path.suffix == '.gz' else path.open('w')) as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def trip(row):
    return any(str(v).upper() == 'TRIP12' or '|Trip12|' in str(v) for v in row.values())


def table(relative):
    name = 'tables/' + relative
    path = checked(BUNDLE / name, MANIFEST[name]['sha256'])
    assert path.stat().st_size == MANIFEST[name]['bytes']
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        columns = reader.fieldnames
        assert columns == MANIFEST[name]['columns']
        selected = [r for r in reader if trip(r)]
    write(HERE / 'tables' / relative, selected, columns)
    return selected


def functions(relative, wanted):
    path = checked(ROOT / relative)
    module = ast.parse(path.read_text())
    selected = [x for x in module.body if isinstance(x, ast.FunctionDef) and x.name in wanted]
    assert len(selected) == len(wanted)
    scope = dict(defaultdict=defaultdict, gzip=gzip, re=re)
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), 'exec'), scope)
    return scope


def main():
    records = json.loads(checked(BUNDLE / 'result-records.json').read_text())
    checked(BUNDLE / 'table-manifest.json')
    for group in ['candidate_screen', 'examples', 'geometry', 'longread']:
        for path in sorted((BUNDLE / 'tables' / group).glob('*.tsv*')):
            table(str(path.relative_to(BUNDLE / 'tables')))
    # GTF structures are not gene-labelled TSVs; retain exact reference IDs.
    families = table('longread/family_eligibility.tsv')
    references = {t for r in families for t in r['reference_transcripts'].split(';')}
    for filename, tokens in [('candidate_exons.mm10.gtf', {'Trip12'}),
                             ('reference_exons.mm39.gtf', references)]:
        name = 'tables/longread/' + filename
        path = checked(BUNDLE / name, MANIFEST[name]['sha256'])
        (HERE / name).write_text(''.join(line for line in path.read_text().splitlines(True)
                                      if any('"' + t + '"' in line for t in tokens)))
    long = functions('notes/2026-10-02-GSE282319-long-read/analyze.py', ['chains', 'mapped', 'lift', 'introns', 'exon_models'])
    bed = functions('notes/2026-10-02-rad23b-ago2-gse25310/check_overlap.py', ['read_bed', 'overlap_bp'])
    audit = json.loads(checked(BUNDLE / 'tables/ago2/audit.json', MANIFEST['tables/ago2/audit.json']['sha256']).read_text())
    expected = {k.split('__')[-1]: v for k, v in audit['inputs'].items() if '/resources/' in k}
    chains = {k: long['chains'](checked(ARCHIVE / filename, expected[filename])) for k, filename in
              [('39_10', 'mm39ToMm10.over.chain.gz'), ('10_39', 'mm10ToMm39.over.chain.gz'),
               ('10_9', 'mm10ToMm9.over.chain.gz'), ('9_10', 'mm9ToMm10.over.chain.gz')]}

    def lift(chrom, start, end):
        c, a, b, first = long['lift'](chains['39_10'], chains['10_39'], chrom, start, end)
        d, u, v, second = long['lift'](chains['10_9'], chains['9_10'], c, a, b)
        assert end - start == b - a == v - u
        return dict(mm10_chrom=c, mm10_start0=a, mm10_end0=b, mm9_chrom=d,
                    mm9_start0=u, mm9_end0=v, mm39_mm10_chain=first, mm10_mm9_chain=second,
                    mapping_status='unique_gapfree_reciprocal')

    features = table('examples/external_query_features.tsv')
    queries = [r for r in features if r['kind'] == 'shared_region' or r['category'] == 'miRNA']
    mapping = [dict(r, **lift(r['chrom'], int(r['start0']), int(r['end0']))) for r in queries]
    write(HERE / 'tables/ago2/query_mapping.tsv', mapping, list(mapping[0]))
    frozen = json.loads(checked(HERE / 'frozen-mm39-sequence.json').read_text())
    union = lift(frozen['chrom'], frozen['start0'], frozen['end0'])
    (HERE / 'mm9-sequence-request.json').write_text(json.dumps(dict(genome='mm9', chrom=union['mm9_chrom'],
                                                               start=union['mm9_start0'], end=union['mm9_end0']), indent=2) + '\n')
    # First invocation exposes an exact query; download its sequence then rerun.
    sequence_path = HERE / 'ucsc-mm9-sequence.json'
    if not sequence_path.exists():
        print('Prepared Trip12 tables and sequence request; AGO2 acceptance awaits sequence download.')
        return
    sequence = json.loads(checked(sequence_path).read_text())
    assert (sequence['genome'], sequence['chrom'], sequence['start'], sequence['end']) == (
        'mm9', union['mm9_chrom'], union['mm9_start0'], union['mm9_end0'])
    assert len(sequence['dna']) == len(frozen['dna']) == frozen['end0'] - frozen['start0']
    identity = []
    for r in mapping:
        s = frozen['dna'][int(r['start0'])-frozen['start0']:int(r['end0'])-frozen['start0']].upper()
        t = sequence['dna'][r['mm9_start0']-sequence['start']:r['mm9_end0']-sequence['start']].upper()
        assert len(s) == len(t) == int(r['end0']) - int(r['start0'])
        identity.append(dict(feature_id=r['feature_id'], sequence_identical=s == t,
                             mm39_positive_strand_sequence=s, mm9_positive_strand_sequence=t,
                             transcript_oriented_sequence=s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]))
    write(HERE / 'tables/ago2/sequence_identity.tsv', identity, list(identity[0]))
    inventory, overlaps, hits = [], [], []
    for prefix, label, unit in TRACKS:
        metadata = json.loads(checked(ARCHIVE / (prefix + '.json'), expected[prefix + '.json']).read_text())[0]
        peaks = bed['read_bed'](checked(ARCHIVE / (prefix + '.bed'), expected[prefix + '.bed']))
        assert str(metadata['PMID']) == '21258322' and len(peaks) == int(metadata['sites'])
        inventory.append(dict(dataset=prefix, label=label, biological_unit=unit, sites=len(peaks)))
        for r in mapping:
            q = dict(chrom=r['mm9_chrom'], start0=r['mm9_start0'], end0=r['mm9_end0'], strand=r['strand'])
            hs = [(p, bed['overlap_bp'](q, p)) for p in peaks if bed['overlap_bp'](q, p) > 0]
            overlaps.append(dict(feature_id=r['feature_id'], label=label, biological_unit=unit,
                                 overlap_clusters=len(hs), overlap_bp_sum=sum(n for _, n in hs)))
            hits.extend(dict(feature_id=r['feature_id'], label=label, biological_unit=unit, **p, overlap_bp=n) for p, n in hs)
    write(HERE / 'tables/ago2/archive_inventory.tsv', inventory, list(inventory[0]))
    write(HERE / 'tables/ago2/overlap_summary.tsv', overlaps, list(overlaps[0]))
    write(HERE / 'tables/ago2/overlapping_clusters.mm9.tsv', hits,
          ['feature_id', 'label', 'biological_unit', 'chrom', 'start0', 'end0', 'name', 'score', 'strand', 'overlap_bp'])
    mapped_hits = []
    for hit in hits:
        c, a, b, _ = long['lift'](chains['9_10'], chains['10_9'], hit['chrom'], hit['start0'], hit['end0'])
        d, u, v, _ = long['lift'](chains['10_39'], chains['39_10'], c, a, b)
        assert v-u == hit['end0']-hit['start0']
        mapped_hits.append(dict(hit, mm39_chrom=d, mm39_start0=u, mm39_end0=v,
                                mm39_mapping_status='unique_gapfree_reciprocal', archived_flanks_may_be_present=True))
    write(HERE / 'tables/ago2/overlapping_clusters.mm39.tsv', mapped_hits,
          ['feature_id', 'label', 'biological_unit', 'chrom', 'start0', 'end0', 'name', 'score', 'strand',
           'overlap_bp', 'mm39_chrom', 'mm39_start0', 'mm39_end0', 'mm39_mapping_status', 'archived_flanks_may_be_present'])
    support = []
    for r in mapping:
        counts = {x['label']: x['overlap_clusters'] for x in overlaps if x['feature_id'] == r['feature_id']}
        wt1 = counts['WT1A'] > 0 or counts['WT1B'] > 0
        wt2 = counts['WT2'] > 0
        support.append(dict(feature_id=r['feature_id'], WT1_any_technical_library=wt1, WT2=wt2,
                            WT_biological_units_supported=int(wt1)+int(wt2)))
    write(HERE / 'tables/ago2/biological_unit_support.tsv', support, list(support[0]))

    # Validate geometry against native site usages, with original contributors.
    geometry = table('geometry/axis-decomposition.tsv')
    native = table('geometry/native-sites.tsv.gz')
    for r in geometry:
        z = [x for x in native if x['route'] == r['route']]
        for state in ['P_usage', 'M_usage']:
            assert abs(sum(float(x[state]) for x in z) - 1) < 1e-12
        reconstructed = sum(float(x['native_contribution']) for x in z)
        decomposed = sum(float(r[k]) for k in ['within_family_DWUI', 'family_mass_DWUI', 'nonterminal_DWUI'])
        assert abs(reconstructed-float(r['native_delta_DWUI'])) < 1e-12
        assert abs(decomposed-float(r['native_delta_DWUI'])) < 1e-12

    transcripts = {r['transcript_id']: r for r in table('longread/transcript_audit.tsv')}
    tpm = {r['transcript_id']: r for r in table('longread/candidate_TPM.tsv')}
    long_audit = json.loads(checked(BUNDLE / 'tables/longread/audit.json', MANIFEST['tables/longread/audit.json']['sha256']).read_text())
    (HERE / 'longread-parent-audit.json').write_text(json.dumps(long_audit, indent=2) + '\n')
    raw_hashes = {k.split('__')[-1]: v for k, v in long_audit['inputs'].items() if '/resources/' in k}
    raw_root = ROOT / 'notes/2026-10-02-GSE282319-long-read/raw'
    gtf = checked(raw_root / 'GSE282319_Repro_filter_tst.gtf.gz', raw_hashes['GSE282319_Repro_filter_tst.gtf.gz'])
    raw_models, _ = long['exon_models'](gtf, {'Trip12'}, 'gene_id')
    assert set(raw_models) == set(transcripts) and len(transcripts) == 76
    for tid, model in raw_models.items():
        assert model['exons'] == [tuple(x) for x in json.loads(transcripts[tid]['exons_mm10'])]
    abundance = checked(raw_root / 'GSE282319_Isoquant.transcript_name_id.tpm.txt.gz', raw_hashes['GSE282319_Isoquant.transcript_name_id.tpm.txt.gz'])
    raw_tpm_ids = set()
    with gzip.open(abundance, 'rt') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            tid = r['transcript_name_id']
            if tid not in tpm:
                continue
            assert tid not in raw_tpm_ids
            raw_tpm_ids.add(tid)
            assert all(abs(float(v)-float(tpm[tid][k])) < 1e-12 for k, v in r.items() if k != 'transcript_name_id')
    assert raw_tpm_ids == set(tpm)
    pair = ['Trip12.G1389.16', 'Trip12.G1389.149']
    models = [json.loads(transcripts[t]['exons_mm39']) for t in pair]
    assert long['introns'](models[0]) == long['introns'](models[1])
    assert all(transcripts[t]['qc_ok'] == 'True' and transcripts[t]['reciprocal_exons'] == 'True' for t in pair)
    columns = [k for k in tpm[pair[0]] if k not in ['gene', 'transcript_id']]
    assert sum(float(tpm[pair[0]][k]) > 0 and float(tpm[pair[1]][k]) > 0 for k in columns) == 0
    boundary = []
    summaries = {r['feature_id']: r for r in table('longread/structural_summary.tsv')}
    for r in features:
        if r['kind'] != 'shared_region':
            continue
        a, b = int(r['start0']), int(r['end0'])
        overlaps_nt = [sum(max(0, min(v, b)-max(u, a)) for u, v in exons) for exons in models]
        assert overlaps_nt[0] == b-a
        boundary.append(dict(feature_id=r['feature_id'], length_nt=b-a, included_transcript=pair[0],
                             proximal_transcript=pair[1], included_overlap_nt=overlaps_nt[0],
                             proximal_overlap_nt=overlaps_nt[1], strict_full_exclusion=overlaps_nt[1] == 0,
                             published_strict_pairs=int(summaries[r['feature_id']]['strict_pairs']),
                             both_TPM_positive_samples=0, same_complete_intron_chain=True,
                             interpretation='structural_compatibility_only; partial_overlap_is_not_strict_exclusion'))
    write(HERE / 'tables/longread/boundary_pair_audit.tsv', boundary, list(boundary[0]))

    # Separate the complete inventory from the three displayed evidence tiers.
    all_sites = table('site_statistics/site_factor_inventory.tsv.gz')
    def reference_rows(relative):
        name = 'tables/' + relative
        path = checked(BUNDLE / name, MANIFEST[name]['sha256'])
        columns, rows = read(path)
        write(HERE / name, rows, columns)
        return rows
    seeds = reference_rows('external_mirna/external_seed_expression.tsv')
    expression = reference_rows('rbp_expression/state_expression.tsv')
    psc = {r['seed'] for r in seeds if r['species'] == 'mouse' and r['expression_group'] == 'PSC'}
    p_routes = defaultdict(set)
    for r in expression:
        if r['species'] == 'mouse' and r['state'] == 'P' and r['route'] in ['CIP', 'WOT'] and r['expressed_any_2CPM'] == 'True':
            p_routes[r['factor']].add(r['route'])
    expressed = {f for f, routes in p_routes.items() if routes == {'CIP', 'WOT'}}
    regions = [r for r in features if r['kind'] == 'shared_region']
    tier_rows = []
    for r in all_sites:
        if r['contrast'] != 'MP' or r['route'] not in ['CIP', 'WOT'] or r['strict_changed'] != 'True':
            continue
        iv = [x for x in regions if int(r['start0']) >= int(x['start0']) and int(r['end0']) <= int(x['end0'])]
        if not iv:
            continue
        assert len(iv) == 1
        if r['category'] == 'miRNA':
            show = r['factor'] in psc
            tier = 'PSC_miRNA' if show else 'other_miRNA'
        elif r['provenance'] == 'RBP':
            show, tier = True, 'direct_RBP_motif'
        else:
            show = r['provenance'] == 'RBP_transfer' and any(n in expressed for n in r['names'].split(';'))
            tier = 'transferred_RBP_P_expressed' if show else 'transferred_RBP_not_P_expressed_both_routes'
        tier_rows.append(dict(r, shared_interval_feature_id=iv[0]['feature_id'], evidence_tier=tier,
                              displayed_in_F5F=show, evidence_scope='predicted_motif_or_seed; not_experimental_Trip12_binding'))
    write(HERE / 'tables/figure5/site-display-evidence.tsv', tier_rows, list(tier_rows[0]))
    tier_summary = []
    for route in ['CIP', 'WOT']:
        for tier in sorted({r['evidence_tier'] for r in tier_rows}):
            z = [r for r in tier_rows if r['route'] == route and r['evidence_tier'] == tier]
            tier_summary.append(dict(route=route, evidence_tier=tier, annotations=len(z),
                                     distinct_coordinates=len({(r['start0'], r['end0']) for r in z})))
        assert {r['evidence_tier']: r['distinct_coordinates'] for r in tier_summary if r['route'] == route and r['evidence_tier'] in
                ['PSC_miRNA', 'direct_RBP_motif', 'transferred_RBP_P_expressed']} == {
                    'PSC_miRNA': 16, 'direct_RBP_motif': 2, 'transferred_RBP_P_expressed': 28}
    write(HERE / 'tables/figure5/site-display-summary.tsv', tier_summary, list(tier_summary[0]))
    table('panel_statistics/inclusion-paired-means.tsv')
    result = dict(source=records['source'], query_features=len(mapping), shared_regions=2,
                  sequence_identical_features=sum(r['sequence_identical'] for r in identity),
                  AGO2_queries_with_WT_unit_support=sum(r['WT_biological_units_supported'] > 0 for r in support),
                  AGO2_hit_rows=len(hits), geometry_residual_verified=True,
                  longread_raw_models_and_TPM_verified=True,
                  longread_boundary_pair=boundary, site_display_summary=tier_summary, inputs=INPUTS,
                  limitations=['Archived AGO2 absence is not evidence of biological absence.',
                               'Long-read TPM describes detection; no differential usage test.',
                               'Native and interval-library contributor sets differ.',
                               'No independent biological replication is inferred from input callers.'])
    (HERE / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    drawings = [
        ('S26A', 'Paired interval inclusion', ['examples/inclusion_paired_libraries.tsv.gz', 'examples/inclusion_state_summary.tsv',
                                             'panel_statistics/inclusion-paired-means.tsv'], 'Two shared intervals; keep 20/50 UMI and paired/all-state means separate'),
        ('S26B', 'Native DWUI and MPRO', ['examples/example_endpoint_effects.tsv.gz', 'geometry/native-contributors.tsv.gz'],
         'supported50_p08; P_M; same_library; CIP and WOT'),
        ('S26C', 'Threshold and leave-one-library influence', ['examples/inclusion_effect_sensitivity.tsv'],
         'Observed ranges and leave-one influence are not confidence intervals'),
        ('S26D', 'External long-read structures', ['longread/transcript_audit.tsv', 'longread/coordinate_liftover.tsv',
                                                  'longread/compatible_isoform_pairs.tsv', 'longread/candidate_TPM.tsv',
                                                  'longread/boundary_pair_audit.tsv'],
         '706 nt: strict pair; 842 nt: 2 nt boundary retention; zero jointly detected samples'),
        ('S26E', 'Archived AGO2 queries', ['ago2/query_mapping.tsv', 'ago2/sequence_identity.tsv', 'ago2/archive_inventory.tsv',
                                          'ago2/overlap_summary.tsv', 'ago2/overlapping_clusters.mm39.tsv', 'ago2/biological_unit_support.tsv'],
         'WT1A/B is one unit; WT2 another; merged tracks have no recovered replicate membership'),
        ('S26F', 'Complete native PAS and terminal family', ['geometry/native-sites.tsv.gz', 'geometry/axis-decomposition.tsv',
                                                            'examples/annotated_example_PAS.tsv'],
         'Full normalization span differs from participating PAS span in WOT'),
        ('S26G', 'Exact native DWUI decomposition', ['geometry/axis-decomposition.tsv', 'geometry/native-sites.tsv.gz'],
         'Three components sum to native DWUI; anchored algebra is not a molecular mechanism'),
        ('S26H', 'All sites and displayed evidence tiers', ['figure5/site-display-evidence.tsv', 'figure5/site-display-summary.tsv',
                                                         'external_mirna/external_seed_expression.tsv', 'rbp_expression/state_expression.tsv'],
         'Different coordinates and annotations are distinct denominators; all are predicted sites'),
    ]
    map_rows, components = [], []
    for panel, title, tables, guard in drawings:
        artifacts = []
        for name in tables:
            path = HERE / 'tables' / name
            columns, rows = read(path)
            assert rows, name
            artifacts.append(dict(path='tables/' + name, sha256=sha(path), bytes=path.stat().st_size, columns=columns))
            map_rows.append(dict(panel=panel, title=title, path='tables/' + name, sha256=sha(path), constraint=guard))
        components.append(dict(id=panel + '.Trip12', panel=panel, title=title, tables=artifacts, constraint=guard))
    write(HERE / 'drawing-map.tsv', map_rows, list(map_rows[0]))
    (HERE / 'components.json').write_text(json.dumps(dict(schema='example_evidence_overlay.v1',
        source=records['source'], parent_result_records_sha256=sha(BUNDLE / 'result-records.json'),
        parent_table_manifest_sha256=sha(BUNDLE / 'table-manifest.json'), components=components,
        status='data_verified; candidate_render_pending; production_unchanged'), indent=2) + '\n')
    output = [dict(path=str(p.relative_to(HERE)), bytes=p.stat().st_size, sha256=sha(p))
              for p in sorted(HERE.rglob('*')) if p.is_file() and p.name not in ['file-manifest.json', 'acceptance.json']]
    (HERE / 'file-manifest.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['inputs', 'source', 'limitations', 'site_display_summary']}, indent=2))


if __name__ == '__main__':
    main()
