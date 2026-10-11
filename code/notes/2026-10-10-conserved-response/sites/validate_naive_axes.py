"""Check actual hN vectors against the frozen complete Liu primary feature axes."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'


def main():
    source = ROOT / 'notes/2026-10-10-marker-reference-sensitivity/inputs.npz'
    with source.open('rb') as stream:
        source_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    assert source_hash == '6652eb94bc70ad1e361992d2991a4e24df9c7d0fe83e694db155ce156bbaccb6'
    with np.load(source, allow_pickle=False) as archive:
        records = json.loads(str(archive['GSE150311_features_json']))
    columns = ['atlas_pas_id', 'primary_gene_id', 'chrom', 'strand', 'cut_site', 'keep_for_apa_primary', 'matched_source_count_50nt']
    features = pd.DataFrame([{key: row[key] for key in columns} for row in records])
    del records
    primary = features[features.keep_for_apa_primary.astype(str).str.lower().eq('true') & features.primary_gene_id.notna()]
    e = pd.read_csv(TABLES / 'tandem/tandem_utr_eligibility.tsv.gz', sep='\t')
    e = e[e.study_id.eq('GSE150311')].set_index('gene_id')
    u = pd.read_csv(TABLES / 'tandem/annotation_terminal_utrs.tsv.gz', sep='\t')
    u = u[u.species.eq('human')]
    member = pd.read_csv(HERE.parent / 'function/comparison_membership_features.tsv.gz', sep='\t', low_memory=False)
    genes = set(member[member.set_name.isin(['hN79', 'dual46'])].human_gene_id)
    pairs = pd.read_csv(TABLES / 'naive_state/library_pair_effects.tsv.gz', sep='\t')
    pairs = pairs[pairs.version.eq('v2_primary') & pairs.setting.eq('supported50_p08') & pairs.node_a.eq('N') & pairs.node_b.eq('P') &
                  pairs.library_a.str.endswith(('6', '8')) & pairs.library_b.str.endswith('7') & pairs.gene_id.isin(genes)]
    rows = []; max_delta_error = 0.
    for gene, group in pairs.groupby('gene_id', sort=True):
        assert len(group) == 2 and group.library_a.nunique() == 2
        full = primary[primary.primary_gene_id.eq(gene)]
        supported = full[full.matched_source_count_50nt.ge(1)]
        expected = set(supported.atlas_pas_id)
        assert all(set(json.loads(value)) == expected for value in group.PAS_ids), gene
        assert full.chrom.nunique() == 1 and full.strand.nunique() == 1
        span = float(full.cut_site.max() - full.cut_site.min())
        assert gene in e.index and span == e.loc[gene, 'native_axis_span_nt'], gene
        pair_flag = bool(group.tandem_utr_compatible.all())
        full_flag = bool(e.loc[gene, 'tandem_utr_compatible'])
        assert pair_flag == full_flag, gene
        lookup = supported.set_index('atlas_pas_id')
        for pair in group.itertuples():
            ids = json.loads(pair.PAS_ids)
            signed = lookup.loc[ids, 'cut_site'].to_numpy(float) * np.where(lookup.loc[ids, 'strand'].eq('+'), 1, -1)
            delta = np.asarray(json.loads(pair.usage_a)) - np.asarray(json.loads(pair.usage_b))
            error = abs(float(delta @ ((signed - signed.min()) / span)) - pair.delta_DWUI)
            assert error < 1e-12, (gene, error)
            max_delta_error = max(max_delta_error, error)
        fits = []
        if full_flag:
            candidates = set(e.loc[gene, 'compatible_transcripts'].split(';'))
            transcript = u[u.gene_id.eq(gene) & u.transcript_id.isin(candidates)]
            retained = full.cut_site.to_numpy() - (1 if full.strand.iloc[0] == '+' else 0)
            for row in transcript.itertuples():
                if str(row.chrom).removeprefix('chr') == str(full.chrom.iloc[0]).removeprefix('chr') and row.strand == full.strand.iloc[0] and retained.min() >= row.utr_start - 1 and retained.max() < row.utr_end:
                    fits.append(row.transcript_id)
            assert set(fits) == candidates, gene
        rows.append(dict(gene_id=gene, actual_condition_pairs=len(group), full_primary_PAS=len(full), supported50_PAS=len(supported),
                         actual_pair_axes_exact=True, actual_full_span_exact=True, pair_tandem=pair_flag,
                         full_axis_UTR_compatible=full_flag, all_full_axis_common_transcripts_verified=';'.join(sorted(fits)),
                         actual_axis_verified=True))
    result = pd.DataFrame(rows)
    assert set(result.gene_id) == genes
    result.to_csv(HERE / 'naive_actual_axis_validation.tsv', sep='\t', index=False)
    report = dict(status='PASS', genes=len(result), actual_condition_pairs=len(pairs), strict_UTR_genes=int(result.full_axis_UTR_compatible.sum()),
                  input_npz_sha256=source_hash, only_features_array_read=True, native_cuts='positive retained cut-1; negative retained cut',
                  max_actual_delta_replay_error=max_delta_error,
                  proof='Both actual supported50 hN/P axes equal the supported subset of the frozen full primary axis; common terminal UTR contains every full-axis retained PAS base.')
    (HERE / 'naive_actual_axis_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
