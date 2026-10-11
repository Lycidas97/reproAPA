"""Exact within-library aggregation of all v7 PAS, with source-block readback."""
import argparse
import gc
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from human_metacell_pas_utils import aggregate
from prepare_human_metacell_pas import OUT, sha256, dump, FIELDS, COUNT_FEATURE_COLUMNS


def gene_matrices(counts, features, genes):
    """Use all eligible sites as the gene denominator, including zero-use sites."""
    primary = features.keep_for_apa_primary.to_numpy(bool)
    assert features.loc[primary, 'primary_gene_id'].notna().all()
    columns = pd.Index(genes.primary_gene_id).get_indexer(features.loc[primary, 'primary_gene_id'])
    assert (columns >= 0).all()
    mapping = sparse.csr_matrix((np.ones(primary.sum(), dtype=np.int64),
        (np.flatnonzero(primary), columns)), shape=(len(features), len(genes)))
    sums = (counts @ mapping).tocsr()
    present = counts.copy()
    present.data = np.ones(len(present.data), dtype=np.int64)
    detected = (present @ mapping).tocsr()
    assert np.array_equal(np.asarray(sums.sum(1)).ravel(), np.asarray(counts[:, primary].sum(1)).ravel())
    return sums, detected, mapping


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--library', required=True)
    args = p.parse_args()
    start = time.monotonic()
    manifest = pd.read_csv(OUT / 'input_manifest.tsv', sep='\t').set_index('library_key')
    record = manifest.loc[args.library]
    out = OUT / 'libraries' / args.library
    source = Path(record.source_projection)
    study = OUT / 'studies' / record.pas_dataset_id
    assert sha256(out / 'cell_to_pas_and_metacell.tsv.gz') == record.selected_binding_sha256
    assert sha256(source / 'cells.tsv.gz') == record.source_cells_sha256
    assert sha256(source / 'projection_audit.json') == record.source_audit_sha256
    feature_contract = pd.read_csv(OUT / 'study_feature_summary.tsv', sep='\t').set_index('dataset_id').loc[record.pas_dataset_id]
    assert sha256(study / 'features.tsv.gz') == feature_contract.compact_features_sha256
    assert sha256(source / 'features.tsv.gz') == record.source_features_sha256
    binding = pd.read_csv(out / 'cell_to_pas_and_metacell.tsv.gz', sep='\t', compression='gzip', dtype={f: str for f in FIELDS})
    metadata = pd.read_csv(out / 'metacell_metadata.tsv.gz', sep='\t', compression='gzip', dtype={f: str for f in FIELDS}).set_index('metacell_id')
    barcodes = pd.read_csv(source / 'cells.tsv.gz', sep='\t', compression='gzip')
    features = pd.read_csv(study / 'features.tsv.gz', sep='\t', compression='gzip', usecols=COUNT_FEATURE_COLUMNS, low_memory=False)
    genes = pd.read_csv(study / 'primary_genes.tsv.gz', sep='\t', compression='gzip')
    status = pd.read_csv(source / 'feature_by_library_status.tsv.gz', sep='\t', compression='gzip')
    assert status.atlas_pas_id.is_unique and set(status.atlas_pas_id) == set(features.atlas_pas_id)
    assert status.library_id.eq(record.library_id).all() and status.dataset_id.eq(record.pas_dataset_id).all()
    status = status.set_index('atlas_pas_id').loc[features.atlas_pas_id].reset_index()
    assert set(status.terminal_status).issubset({'observed', 'not_observed'})
    assert status.loc[status.terminal_status.eq('not_observed'), 'projected_umi'].eq(0).all()
    path = source / 'feature_by_cell.mtx.gz'
    before = path.stat()
    matrix_hash = sha256(path)
    with threadpool_limits(limits=2):
        full = mmread(path, spmatrix=True).T.tocsr()
    full.sum_duplicates(); full.eliminate_zeros()
    assert full.shape == (int(record.source_pas_cells), int(record.atlas_features))
    assert np.issubdtype(full.dtype, np.integer) and (full.data > 0).all()
    assert full.nnz == int(record.source_projected_nonzero)
    assert int(full.sum()) == int(record.source_projected_umi)
    all_totals = np.asarray(full.sum(0)).ravel()
    assert np.array_equal(all_totals, status.projected_umi.to_numpy())
    index = pd.Index(barcodes.original_barcode).get_indexer(binding.barcode)
    assert (index >= 0).all() and np.array_equal(index + 1, binding.pas_column_index_1based)
    assert np.array_equal(barcodes.cell_key.iloc[index], binding.pas_cell_key)
    matrix = ad.AnnData(X=full, obs=pd.DataFrame(index=pd.Index(barcodes.original_barcode, name='barcode')),
        var=pd.DataFrame(index=pd.Index(features.atlas_pas_id, name='atlas_pas_id')))
    result = aggregate(matrix, binding)
    y = result.X.tocsr(); y.eliminate_zeros()
    assert set(result.obs_names) == set(metadata.index)
    assert np.array_equal(result.obs.n_cells, metadata.loc[result.obs_names, 'n_cells'])
    selected = full[index]
    selected_totals = np.asarray(selected.sum(0)).ravel()
    assert np.array_equal(np.asarray(y.sum(0)).ravel(), selected_totals)
    cell_coverage = binding[['cell_id', 'barcode', 'metacell_id', 'library_key', 'local_state_id', 'time_interpretation']].copy()
    cell_coverage['pas_umi'] = np.asarray(selected.sum(1)).ravel()
    cell_coverage['detected_pas'] = np.diff(selected.indptr)
    cell_coverage['pas_match_status'] = np.where(cell_coverage.pas_umi.eq(0), 'matched_zero_pas_umi', 'matched_with_counts')
    cell_coverage.to_csv(out / 'cell_pas_coverage.tsv.gz', sep='\t', compression='gzip', index=False)
    del selected, matrix
    gc.collect()
    sparse.save_npz(out / 'metacell_by_pas_counts.npz', y)
    pd.DataFrame({'metacell_id': result.obs_names}).to_csv(out / 'matrix_rows.tsv.gz', sep='\t', compression='gzip', index=False)
    coverage = metadata.loc[result.obs_names].copy()
    coverage['pas_dataset_id'] = record.pas_dataset_id
    coverage['total_pas_umi'] = np.asarray(y.sum(1)).ravel()
    coverage['detected_pas'] = np.diff(y.indptr)
    primary = features.keep_for_apa_primary.to_numpy(bool)
    sensitivity = features.keep_for_sensitivity.to_numpy(bool)
    for name, mask in [('primary', primary), ('sensitivity', sensitivity), ('relaxed_only', sensitivity & ~primary)]:
        part = y[:, mask]
        coverage[name + '_pas_umi'] = np.asarray(part.sum(1)).ravel()
        coverage[name + '_detected_pas'] = np.diff(part.indptr)
    gene_y, gene_sites, mapping = gene_matrices(y, features, genes)
    sparse.save_npz(out / 'metacell_by_primary_gene_counts.npz', gene_y)
    sparse.save_npz(out / 'metacell_by_primary_gene_detected_sites.npz', gene_sites)
    coverage['primary_detected_genes'] = np.diff(gene_y.indptr)
    coverage['genes_with_ge2_detected_primary_sites'] = np.asarray((gene_sites >= 2).sum(1)).ravel()
    eligible = genes.n_primary_atlas_pas.ge(2).to_numpy()
    for threshold in [10, 20, 50]:
        coverage[f'genes_with_ge2_primary_atlas_sites_and_ge{threshold}_umi'] = np.asarray((gene_y[:, eligible] >= threshold).sum(1)).ravel()
    coverage['time_interpretation'] = binding.time_interpretation.iloc[0]
    coverage.to_csv(out / 'metacell_pas_coverage.tsv.gz', sep='\t', compression='gzip')
    gene_library = genes.copy()
    gene_library['n_primary_pas_observed_in_library'] = np.asarray(mapping.T @ status.terminal_status.eq('observed').to_numpy(np.int64)).ravel()
    gene_library['primary_umi_selected_cells'] = np.asarray(gene_y.sum(0)).ravel()
    gene_library['metacells_with_primary_umi'] = np.asarray((gene_y > 0).sum(0)).ravel()
    gene_library.to_csv(out / 'primary_gene_library_coverage.tsv.gz', sep='\t', compression='gzip', index=False)
    ledger = pd.DataFrame({'feature_index_1based': features.feature_index_1based,
        'terminal_status': status.terminal_status, 'n_mapped_sample_pas': status.n_mapped_sample_pas,
        'source_library_umi': all_totals, 'selected_cell_umi': selected_totals,
        'metacells_nonzero': np.asarray((y > 0).sum(0)).ravel()})
    ledger.to_csv(out / 'feature_library_coverage.tsv.gz', sep='\t', compression='gzip', index=False)
    strata = features[['confidence_tier', 'primary_region_class']].copy()
    strata['selected_cell_umi'] = selected_totals
    strata.groupby(['confidence_tier', 'primary_region_class'], dropna=False).agg(
        atlas_pas=('selected_cell_umi', 'size'), detected_pas=('selected_cell_umi', lambda z: int(z.gt(0).sum())),
        selected_cell_umi=('selected_cell_umi', 'sum')).reset_index().to_csv(out / 'coverage_by_annotation.tsv', sep='\t', index=False)
    # Reopen the delivered matrix; compare independent source cell-column sums.
    check = sparse.load_npz(out / 'metacell_by_pas_counts.npz').tocsr()
    assert check.shape == y.shape and (check != y).nnz == 0
    assert np.array_equal(np.asarray(check.sum(0)).ravel(), selected_totals)
    sampled = sorted({str(result.obs_names[0]), str(result.obs_names[-1]),
                      str(coverage.n_cells.idxmin()), str(coverage.n_cells.idxmax())})
    block_checks = []
    for mc in sampled:
        mask = binding.metacell_id.eq(mc).to_numpy()
        expected = np.asarray(full[index[mask]].sum(0)).ravel()
        got = check.getrow(result.obs_names.get_loc(mc)).toarray().ravel()
        assert np.array_equal(expected, got), f'Incorrect membership aggregation: {mc}'
        block_checks.append(dict(metacell_id=mc, cells=int(mask.sum()), umi=int(expected.sum()),
            all_pas_exact=True, source_vector_sha256=__import__('hashlib').sha256(expected.astype('<i8').tobytes()).hexdigest()))
    dump(out / 'direct_source_block_readback.json', block_checks)
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    summary = dict(status='PAS_LIBRARY_AGGREGATION_AND_SOURCE_READBACK_PASS', library_key=args.library,
        library_id=record.library_id, dataset=record.dataset, pas_dataset_id=record.pas_dataset_id,
        cells=len(binding), metacells=y.shape[0], atlas_features=y.shape[1], nonzero_metacell_pas=y.nnz,
        source_cells=int(record.source_pas_cells), source_projected_umi=int(record.source_projected_umi),
        selected_pas_umi=int(y.sum()), primary_pas_umi=int(coverage.primary_pas_umi.sum()),
        sensitivity_pas_umi=int(coverage.sensitivity_pas_umi.sum()),
        zero_pas_cells=int(cell_coverage.pas_umi.eq(0).sum()), zero_pas_metacells=int(coverage.total_pas_umi.eq(0).sum()),
        all_selected_barcodes_matched=True, all_features_retained=True,
        per_pas_count_conservation=True, source_library_per_pas_counts_match_v7=True,
        source_matrix_sha256=matrix_hash, source_matrix_unchanged_during_run=True,
        direct_source_metacell_blocks=len(block_checks), output_matrix_sha256=sha256(out / 'metacell_by_pas_counts.npz'),
        elapsed_seconds=round(time.monotonic() - start, 2))
    dump(out / 'summary.json', summary)
    (out / 'complete.ok').write_text('PAS_V7_LIBRARY_COMPLETE\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
