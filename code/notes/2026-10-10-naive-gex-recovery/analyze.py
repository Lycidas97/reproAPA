"""Recover genuine naive-state GEX from unchanged frozen all-gene MC counts."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.12g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def csr_dense_rows(archive, rows):
    """Decode only requested CSR rows, retaining integer counts without scipy."""
    assert archive['format'].item() in [b'csr', 'csr']
    data, indices, indptr = (archive[k] for k in ['data', 'indices', 'indptr'])
    out = np.zeros((len(rows), int(archive['shape'][1])), dtype=np.int64)
    for j, row in enumerate(rows):
        start, stop = indptr[row:row + 2]
        # add.at conserves counts even if a CSR archive has duplicate columns.
        np.add.at(out[j], indices[start:stop], data[start:stop])
    return out


def expression_rows(meta, counts, totals, genes):
    assert len(meta) == len(counts) == len(totals) and (totals > 0).all()
    pieces = []
    keys = ['library_key', 'sample_accession', 'node', 'condition', 'experimental_time']
    for key, q in meta.groupby(keys, sort=True, dropna=False):
        ix = q.index.to_numpy(); c = counts[ix]
        cpm = c / totals[ix, None] * 1e6
        record = dict(zip(keys, key))
        record.update(state=record.pop('node'), time=record.pop('experimental_time'),
                      species='human', study_id='GSE150311', route='GSE150311',
                      n_metacells=len(ix), n_cells=int(q.n_cells.sum()),
                      qualified=len(ix) >= 3 and q.n_cells.sum() >= 30,
                      mapped=True, total_gex=int(totals[ix].sum()))
        f = genes.copy()
        for column, value in record.items():
            f[column] = value
        f['count'] = c.sum(0)
        f['pooled_CPM'] = c.sum(0) / totals[ix].sum() * 1e6
        f['mean_log2CPM1'] = np.log2(cpm + 1).mean(0)
        f['mean_CPM'] = cpm.mean(0)
        f['detected_metacells'] = (c > 0).sum(0)
        f['detected_fraction'] = (c > 0).mean(0)
        pieces.append(f)
    return pd.concat(pieces, ignore_index=True)


def contrast_rows(libraries, plans):
    q = libraries[libraries.qualified].set_index(['gene_id', 'sample_accession', 'state'])
    rows = []
    genes = sorted(libraries.gene_id.unique())
    for contrast, pairs in plans.items():
        parts = []; observed_pairs = []
        for sample_a, state_a, sample_b, state_b in pairs:
            a = q.xs((sample_a, state_a), level=['sample_accession', 'state']) if (sample_a, state_a) in q.index.droplevel('gene_id') else pd.DataFrame()
            b = q.xs((sample_b, state_b), level=['sample_accession', 'state']) if (sample_b, state_b) in q.index.droplevel('gene_id') else pd.DataFrame()
            if a.empty or b.empty:
                continue
            a = a.reindex(genes); b = b.reindex(genes)
            assert a.mean_log2CPM1.notna().all() and b.mean_log2CPM1.notna().all()
            observed_pairs.append((sample_a, state_a, sample_b, state_b))
            parts.append(pd.DataFrame(dict(gene_id=genes,
                target_GEX=a.mean_log2CPM1.to_numpy(), reference_GEX=b.mean_log2CPM1.to_numpy(),
                target_mean_CPM=a.mean_CPM.to_numpy(), reference_mean_CPM=b.mean_CPM.to_numpy(),
                target_pooled_CPM=a.pooled_CPM.to_numpy(), reference_pooled_CPM=b.pooled_CPM.to_numpy())))
        if not parts:
            continue
        f = pd.concat(parts).groupby('gene_id', sort=True).mean().reset_index()
        f['contrast'] = contrast
        f['n_pairs'] = len(parts)
        f['planned_pairs'] = len(pairs)
        f['complete_contributors'] = len(parts) == len(pairs)
        f['delta_log2CPM1'] = f.target_GEX - f.reference_GEX
        f['mean_GEX'] = (f.target_GEX + f.reference_GEX) / 2
        f['mean_CPM'] = (f.target_mean_CPM + f.reference_mean_CPM) / 2
        f['target_libraries'] = ';'.join(sorted(set(a for a, _, _, _ in observed_pairs)))
        f['reference_libraries'] = ';'.join(sorted(set(b for _, _, b, _ in observed_pairs)))
        f['planned_target_libraries'] = ';'.join(sorted(set(a for a, _, _, _ in pairs)))
        f['planned_reference_libraries'] = ';'.join(sorted(set(b for _, _, b, _ in pairs)))
        rows.append(f)
    return pd.concat(rows, ignore_index=True)


def main():
    protocol = json.loads((HERE / 'protocol.json').read_text())
    sources = json.loads((HERE / 'local_input_sha256.json').read_text())
    for path in [Path(__file__), HERE / 'protocol.json', TABLES / 'naive_state/library_pair_effects.tsv.gz']:
        sources[str(path)] = digest(path)
    for path in sorted((HERE / 'raw').glob('*')):
        sources[str(path)] = digest(path)
    remote_hashes = {}
    for line in (HERE / 'remote_input_sha256.txt').read_text().splitlines():
        value, remote = line.split(maxsplit=1)
        local = HERE / 'raw' / Path(remote).name
        assert digest(local) == value
        remote_hashes[remote] = value
    archive = np.load(HERE / 'raw/metacell_by_gene_counts.npz', allow_pickle=False)
    gmeta = pd.read_csv(HERE / 'raw/metacells.tsv.gz', sep='\t')
    genes = pd.read_csv(HERE / 'raw/genes.tsv.gz', sep='\t')
    assert genes.gene_id.is_unique and gmeta.metacell_id.is_unique
    assert np.array_equal(archive['shape'], [len(gmeta), len(genes)])
    assert np.issubdtype(archive['data'].dtype, np.integer) and (archive['data'] >= 0).all()
    full_rows = csr_dense_rows(archive, np.arange(len(gmeta)))
    assert np.array_equal(full_rows.sum(1), gmeta.total_gex.to_numpy())
    meta = pd.read_csv(TABLES / 'naive_state/frozen_selected_metadata.tsv.gz', sep='\t')
    ix = pd.Index(gmeta.metacell_id).get_indexer(meta.metacell_id)
    assert (ix >= 0).all()
    assert np.array_equal(ix, meta.matrix_row.to_numpy())
    assert np.array_equal(gmeta.n_cells.iloc[ix], meta.n_cells)
    assert np.array_equal(gmeta.library_key.iloc[ix], meta.library_key)
    assert np.array_equal(gmeta.condition.iloc[ix], meta.condition)
    full_selected = meta.copy()
    qualified = meta.purity.ge(.8).to_numpy()
    meta = meta[qualified].reset_index(drop=True)
    counts = full_rows[ix][qualified]
    total = gmeta.total_gex.to_numpy()[ix][qualified]
    libraries = expression_rows(meta, counts, total, genes)
    save(libraries, 'all_gene_library_expression.tsv.gz')
    inventory = libraries.drop_duplicates(['library_key', 'state'])[
        ['sample_accession', 'library_key', 'state', 'condition', 'time', 'n_metacells', 'n_cells', 'qualified', 'total_gex']]
    save(inventory, 'qualified_inventory.tsv')
    contrasts = contrast_rows(libraries, protocol['fixed_naive_comparisons'])
    contrasts = contrasts.merge(genes, on='gene_id', validate='many_to_one')
    save(contrasts, 'all_gene_contrast_expression.tsv.gz')
    requested = protocol['requested_factors']
    mappings = []
    for factor in requested:
        g = genes[genes.gene_symbol.str.upper().eq(factor)]
        mappings.append(dict(factor=factor, status='unique_exact_casefold_symbol' if len(g) == 1 else 'absent' if not len(g) else 'ambiguous',
                             gene_ids=';'.join(g.gene_id)))
    mapping = pd.DataFrame(mappings)
    assert mapping.status.eq('unique_exact_casefold_symbol').all()
    mapping['gene_id'] = mapping.gene_ids
    save(mapping, 'factor_gene_mapping.tsv')
    factor_libraries = libraries.merge(mapping[['factor', 'gene_id']], on='gene_id', validate='many_to_one')
    save(factor_libraries, 'factor_library_expression.tsv.gz')
    save(contrasts.merge(mapping[['factor', 'gene_id']], on='gene_id', validate='many_to_one'), 'factor_contrast_expression.tsv')
    factor_summary = []
    for factor, z in factor_libraries.groupby('factor', sort=True):
        z = z.set_index(['sample_accession', 'state'])
        n6 = z.loc[('GSM4546546', 'N')]; n8 = z.loc[('GSM4546548', 'N')]
        p7 = z.loc[('GSM4546547', 'P')]; p8 = z.loc[('GSM4546548', 'P')]
        factor_summary.append(dict(factor=factor, gene_id=n6.gene_id,
            N6_P7_delta_log2CPM1=n6.mean_log2CPM1-p7.mean_log2CPM1,
            N8_P7_delta_log2CPM1=n8.mean_log2CPM1-p7.mean_log2CPM1,
            N8_P8_delta_log2CPM1=n8.mean_log2CPM1-p8.mean_log2CPM1,
            hN_P7_delta_log2CPM1=(n6.mean_log2CPM1+n8.mean_log2CPM1)/2-p7.mean_log2CPM1,
            hN_P7_delta_pooled_CPM=(n6.pooled_CPM+n8.pooled_CPM)/2-p7.pooled_CPM,
            hN_P7_delta_mean_CPM=(n6.mean_CPM+n8.mean_CPM)/2-p7.mean_CPM))
    save(pd.DataFrame(factor_summary), 'factor_naive_summary.tsv')
    covariates = genes.copy()
    for label, sample, state in [('N6', 'GSM4546546', 'N'), ('N8', 'GSM4546548', 'N'), ('P7', 'GSM4546547', 'P')]:
        z = libraries[libraries.sample_accession.eq(sample) & libraries.state.eq(state)].set_index('gene_id')
        assert len(z) == len(genes) and z.qualified.all()
        for column in ['mean_log2CPM1', 'mean_CPM', 'pooled_CPM']:
            covariates[label + '_' + column] = covariates.gene_id.map(z[column])
    covariates['mean_GEX'] = (covariates.N6_mean_log2CPM1 + covariates.N8_mean_log2CPM1 + 2 * covariates.P7_mean_log2CPM1) / 4
    covariates['mean_CPM'] = (covariates.N6_mean_CPM + covariates.N8_mean_CPM + 2 * covariates.P7_mean_CPM) / 4
    pairs = pd.read_csv(TABLES / 'naive_state/library_pair_effects.tsv.gz', sep='\t')
    pairs = pairs[pairs.setting.eq('supported50_p08') & pairs.node_a.eq('N') & pairs.node_b.eq('P') & pairs.library_b.str.endswith('GSM4546547')]
    member_sets = [set(pairs[pairs.library_a.str.endswith(sample)].gene_id) for sample in ['GSM4546546', 'GSM4546548']]
    common = set.intersection(*member_sets)
    assert len(common) == 2575
    covariates['naive_APA_2575_member'] = covariates.gene_id.isin(common)
    assert common.issubset(set(covariates.gene_id))
    save(covariates, 'naive_gex_covariates.tsv.gz')
    for path, value in sources.items():
        assert digest(path) == value, path
    output_hashes = {str(p): digest(p) for p in sorted(HERE.glob('*.tsv*'))}
    audit = dict(status='PASS', source_hashes=sources, remote_source_hashes=remote_hashes, output_sha256=output_hashes,
        frozen_metacells_total=len(gmeta), frozen_selected_metacells=len(full_selected), purity_qualified_selected_metacells=len(meta),
        genes=len(genes), all_gene_count_sum=int(full_rows.sum()), qualified_naive_APA_genes=len(common),
        library_state_groups=len(inventory), target_factor_unique_mapping=len(mapping),
        normalization=protocol['normalization'], gene_specific_APA_contributor_filter=False,
        symbol_repair='not needed for requested factors: all 17 exact unique original GEX symbols; ID placeholders retained for other genes',
        checks=dict(raw_GEX_archive_matches_remote=True, raw_GEX_matrix_integer=True, all_gene_GEX_denominator_conserved=True,
                    frozen_native_state_identity=True, cell_numbers_identical=True, library_and_condition_identical=True,
                    N_GEX_observed_directly=True, no_PAS_recount=True),
        scope='Descriptive fixed libraries/conditions; no independent biological replication or protein-activity inference')
    (HERE / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({k: audit[k] for k in ['status', 'genes', 'frozen_metacells_total', 'purity_qualified_selected_metacells', 'qualified_naive_APA_genes', 'target_factor_unique_mapping']}, indent=2))


if __name__ == '__main__':
    main()
