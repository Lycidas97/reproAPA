"""Amendments to the single-file Supplementary Tables after the 15:30 main-text check (2026-10-06).

See MAIN-TEXT-CHECK-1530.zh.md. Each change below names the gap it closes. Values are copied from contract-r2
bundle files (SHA-256 checked against table-manifest.json) or from the Figure 2 panel values; the only computed
content is labelled 'derived' in the workbook: row counts, a subtraction, CDS length / 3 and a coordinate
difference, each from data shown in the same sheet.
"""
import hashlib, json, re
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BUNDLE = ROOT / 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle'
F2_VALUES = 'notes/2026-10-06-main-figures-polyapipe/candidates/run3/F2/panel-values.json'

PRIMARY_S16C = 'shared_p08|requalified|0.02|0.02'   # the only analysis reproducing Figure S16c (checked 15:40)
PDLIM7 = {'ENSG00000196923', 'ENSMUSG00000021493'}
PDLIM7_TX = ['ENST00000355572', 'ENST00000355841', 'ENSMUST00000069968', 'ENSMUST00000046246']

# item 5: titles that promised content the sheets did not hold
TITLES = {
    'S2d': 'Common one-to-one ortholog background (1,306 genes)',
    'S4a': 'Contrast design, series and source libraries',
}

# sheet blocks added from bundle files: (sheet, path) -> optional (description, filter) ; item numbers in comments
ADD = {
    ('S1a', 'frozen/marker_by_state.tsv'): None,                                   # 4 marker evidence
    ('S1a', 'state_grouping/marker_by_state.tsv.gz'): None,                         # 4
    ('S4a', 'global/planned_route_contrasts.tsv'): None,                            # 5 design
    ('S4c', 'intermediate/process/same_series_process_classes.tsv.gz'): (           # 1
        f"analysis == '{PRIMARY_S16C}' (primary analysis of Figure S11c; the 26 sensitivity analyses are in the "
        "source-data package)",
        lambda d: d[d['analysis'] == PRIMARY_S16C]),
    ('S5a', 'direction_difference/original_pair_direction_sensitivity.tsv'): None,  # 2 the 32 pair tests
    ('S5d', 'same_pas/v2_primary_annotated_usage.tsv.gz'): (                        # 3 usage per site and state
        'PDLIM7/Pdlim7 rows only (gene_id ENSG00000196923, ENSMUSG00000021493)',
        lambda d: d[d['gene_id'].isin(PDLIM7)]),
    ('S5d', 'same_pas/v2_primary_same_axis_gene_effects.tsv.gz'): (                 # 3 mass_min 0.918 / 0.886
        'PDLIM7/Pdlim7 rows at the primary 25-nt tolerance',
        lambda d: d[(d['human_gene_id'] == 'ENSG00000196923') & (d['tolerance'] == '25')]),
    ('S5d', 'supplement_statistics/same-PAS-example-exon-CDS.tsv.gz'): (            # 3 isoform structure
        'the four transcripts ending at the matched PASs (' + ', '.join(PDLIM7_TX) + ')',
        lambda d: d[d['transcript_id'].isin(PDLIM7_TX)]),
}

# filters on blocks that were already sheet content
FILTER = {
    ('S7b', 'examples/inclusion_effect_sensitivity.tsv'): (                         # 6
        "gene == TRIP12 (the other example genes are in the source-data package)",
        lambda d: d[d['gene'].str.upper() == 'TRIP12']),
}

# derived columns appended to an existing block
ADD_COLS = {
    ('S6c', 'weighted_retention/retention_summary.tsv'): [                          # 7 the 107 records
        ('genes_without_complete_sites', lambda d: (pd.to_numeric(d['genes']) - pd.to_numeric(d['genes_with_complete_sites'])).astype(int).astype(str),
         'genes_without_complete_sites = genes − genes_with_complete_sites (derived)'),
    ],
}


def _tsv(path):
    p = BUNDLE / 'tables' / path
    return pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False, na_filter=False,
                       compression='gzip' if path.endswith('.gz') else None)


def s4c_counts():
    d = _tsv('intermediate/process/same_series_process_classes.tsv.gz')
    d = d[d['analysis'] == PRIMARY_S16C]
    order = ['endpoint_aligned', 'reverse', 'within_threshold', 'direction_difference',
             'mixed_or_threshold_undetermined', 'endpoint_not_retained', 'unavailable']
    rows = []
    for sp, st in [('human', 'R'), ('mouse', 'R'), ('mouse', 'E')]:
        g = d[(d['species'] == sp) & (d['state'] == st)]['pair_class'].value_counts()
        if set(g.index) - set(order):
            raise ValueError('unexpected pair_class')
        rows.append({'species': sp, 'intermediate_state': st, **{c: str(int(g.get(c, 0))) for c in order},
                     'total': str(int(g.sum()))})
    return pd.DataFrame(rows)


def f2e_correlations():
    e = json.load(open(ROOT / F2_VALUES))['E']
    lab = {'heldout': 'centres refitted on the remaining routes', 'alone': 'full-fit centres'}
    rows = []
    for k, v in e.items():
        route, cls, kind = (k.split('|') + [''])[:3]
        if kind in lab:
            rows.append({'heldout_route': route, 'class': cls, 'held_out_route_means_compared_with': lab[kind],
                         'r': repr(float(v['r'])), 'RMSE': repr(float(v['RMSE'])), 'n_genes': str(int(v['n']))})
    return pd.DataFrame(rows)


def pdlim7_isoforms():
    tm = _tsv('geometry/transcript-models.tsv.gz')
    tm = tm[tm['transcript_id'].isin(PDLIM7_TX)][['species', 'gene_id', 'transcript_id', 'chrom', 'strand', 'start0', 'end0']]
    e = _tsv('supplement_statistics/same-PAS-example-exon-CDS.tsv.gz')
    c = e[(e['kind'] == 'CDS') & e['transcript_id'].isin(PDLIM7_TX)].copy()
    c['L'] = c['end0'].astype(int) - c['start0'].astype(int)
    L = c.groupby('transcript_id')['L'].sum()
    tm = tm.rename(columns={'start0': 'transcript_start0', 'end0': 'transcript_end0'})
    tm['CDS_nt'] = tm['transcript_id'].map(L).astype(int).astype(str)
    tm['protein_aa'] = (tm['CDS_nt'].astype(int) // 3).astype(str)
    if any(int(x) % 3 for x in tm['CDS_nt']):
        raise ValueError('CDS length not a multiple of 3')
    tm['order'] = tm['transcript_id'].map({t: i for i, t in enumerate(PDLIM7_TX)})
    return tm.sort_values('order').drop(columns='order').reset_index(drop=True)


def pdlim7_separation():
    s = _tsv('supplement_statistics/same-PAS-example-sites.tsv.gz')
    rows = []
    for sp, col in [('human', 'human_id'), ('mouse', 'mouse_id')]:
        parts = [x.split('|') for x in s[col]]
        strand = {p[3] for p in parts}
        if len(strand) != 1:
            raise ValueError('strand')
        strand = strand.pop()
        pos = sorted(int(p[4]) for p in parts)
        up, down = (pos[1], pos[0]) if strand == '-' else (pos[0], pos[1])
        rows.append({'species': sp, 'gene_id': parts[0][1], 'chrom': parts[0][2], 'strand': strand,
                     'upstream_PAS': str(up), 'downstream_PAS': str(down), 'separation_nt': str(abs(up - down))})
    return pd.DataFrame(rows)


# derived blocks: sheet -> [(name, function, description, input paths relative to ROOT)]
B_ = 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/'
DERIVED = {
    'S4c': [('intermediate_state_class_counts', s4c_counts,                         # 1
             'Derived: number of genes per pair_class in the same_series_process_classes block (primary analysis); '
             'these are the counts drawn in Figure S11c.',
             [B_ + 'intermediate/process/same_series_process_classes.tsv.gz'])],
    'S3b': [('held-out_route_mean_correlations', f2e_correlations,                  # 7
             'Copied from the Figure 2e panel values: correlation (r) and RMSE between held-out-route means of '
             'full-fit class members and either the centres refitted without that route or the full-fit centres.',
             [F2_VALUES])],
    'S5d': [('alternative_terminal_exon_isoforms', pdlim7_isoforms,                 # 3
             'Derived from the exon/CDS block: CDS_nt is the summed CDS length (stop codon excluded) and '
             'protein_aa = CDS_nt / 3; transcript coordinates are copied from Ensembl release 116 models.',
             [B_ + 'geometry/transcript-models.tsv.gz', B_ + 'supplement_statistics/same-PAS-example-exon-CDS.tsv.gz']),
            ('matched_PAS_separation', pdlim7_separation,
             'Derived from the matched-site block: genomic distance between the two matched PASs '
             '(upstream = proximal on the minus strand).',
             [B_ + 'supplement_statistics/same-PAS-example-sites.tsv.gz'])],
}


def manifest_sha():
    return {t['path'].replace('tables/', '', 1): t['sha256'] for t in json.load(open(BUNDLE / 'table-manifest.json'))['tables']}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
