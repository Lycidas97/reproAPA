"""Bind human membership only to passed PAS releases; record unavailable studies."""
import gzip
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from build_human_metacells import P, OUT as GEX, sha256 as file_sha256, dump, FIELDS

OUT = P / os.environ.get('HUMAN_METACELL_PAS_OUTPUT', 'results_metacell_pas_v1')
ATLAS = Path('/root/nfsdata/apadetect/data/reprogramming_differentiation_apa/analyses/apadetect_reproapa_study_atlas_v7_merge_first_6aabc22_2026-09-05')
NEW_ATLAS = ATLAS.parent / 'apadetect_gse221739_study_atlas_v1_6aabc22_2026-09-09'
DATASET_ALIASES = {'gse178325_complete_32lib': 'gse178325_human_chemical_pluripotency_reprogramming'}
COUNT_FEATURE_COLUMNS = ['atlas_pas_id', 'feature_index_1based', 'keep_for_apa_primary',
                         'keep_for_sensitivity', 'primary_gene_id', 'confidence_tier', 'primary_region_class']



def sha256(path, decompressed=False):
    if not decompressed:
        return file_sha256(path)
    import hashlib
    h = hashlib.sha256()
    with gzip.open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    assert not (OUT / 'preflight_summary.json').exists(), 'PAS binding already frozen'
    manifest = pd.read_csv(GEX / 'prepared/libraries.tsv', sep='\t')
    study_scope = os.environ.get('HUMAN_METACELL_PAS_STUDY')
    if study_scope:
        manifest = manifest[manifest.study_id.eq(study_scope)].copy()
        assert len(manifest), f'No frozen libraries for requested study: {study_scope}'
    cells = pd.read_csv(GEX / 'release/cell_to_metacell.tsv.gz', sep='\t', compression='gzip', dtype=str).fillna('')
    nodes = pd.read_csv(GEX / 'release/metacell_metadata.tsv.gz', sep='\t', compression='gzip', dtype={f: str for f in FIELDS})
    releases, publications = {}, []
    for root in [ATLAS, NEW_ATLAS]:
        marker = root / 'release/release_complete.json'
        if not marker.exists():
            continue
        release = json.loads(marker.read_text())
        assert release['status'] == 'passed' and release['all_projection_count_conservation_passed']
        releases[str(root)] = release
        frame = pd.read_csv(root / 'release/library_projection_summary.tsv', sep='\t')
        frame['atlas_root'] = str(root)
        publications.append(frame)
    pub = pd.concat(publications, ignore_index=True)
    assert pub.library_id.is_unique
    bound, availability, selected = [], [], []
    for r in manifest.itertuples():
        candidates = pub[pub.library_id.eq(r.library_id)]
        if not len(candidates):
            # GSE118258's human GEX registry uses short day IDs; published PAS uses GEO accessions.
            candidates = pub[pub.library_id.eq(f'{r.dataset_id}_{r.sample_accession}')]
        assert len(candidates) <= 1
        if not len(candidates):
            availability.append(dict(library_key=r.library_key, study_id=r.study_id, cells=r.primary_cells,
                metacells=int(nodes.library_key.eq(r.library_key).sum()), status='UPSTREAM_PAS_RELEASE_PENDING',
                source_atlas=str(NEW_ATLAS) if r.study_id == 'GSE221739' else '', pas_library_id=''))
            continue
        row = candidates.iloc[0]
        assert row.dataset_id == DATASET_ALIASES.get(r.dataset_id, r.dataset_id) and row.status == 'passed'
        assert row.barcode_order_passed and row.count_conservation_passed
        source = Path(row.atlas_root) / 'studies' / row.dataset_id / 'libraries' / row.library_id
        barcodes = pd.read_csv(source / 'cells.tsv.gz', sep='\t', compression='gzip')
        assert barcodes.original_barcode.is_unique and barcodes.cell_key.is_unique
        assert barcodes.library_id.eq(row.library_id).all()
        assert np.array_equal(barcodes.cell_index_1based, np.arange(1, len(barcodes)+1))
        with h5py.File(r.input_path, 'r') as h:
            raw_barcodes = h['matrix/barcodes'][:].astype(str)
        # This independently validates aliases using the complete called-cell axis, not just an accession string.
        assert np.array_equal(raw_barcodes, barcodes.original_barcode), f'GEX/PAS called-cell axis mismatch: {r.library_key}'
        audit = json.loads((source / 'projection_audit.json').read_text())
        assert audit['status'] == 'passed' and audit['n_cells'] == len(barcodes) == row.n_cells
        selected.append((r, row, source, barcodes))
        availability.append(dict(library_key=r.library_key, study_id=r.study_id, cells=r.primary_cells,
            metacells=int(nodes.library_key.eq(r.library_key).sum()), status='PAS_RELEASE_AVAILABLE',
            source_atlas=row.atlas_root, pas_library_id=row.library_id))
    feature_rows, fingerprints = [], {}
    for dataset in sorted({row.dataset_id for _, row, _, _ in selected}):
        row = next(row for _, row, _, _ in selected if row.dataset_id == dataset)
        source = Path(row.atlas_root) / 'studies' / dataset
        dest = OUT / 'studies' / dataset
        dest.mkdir(parents=True, exist_ok=True)
        f = pd.read_csv(source / 'matrix/features.tsv.gz', sep='\t', compression='gzip')
        annotation = pd.read_csv(source / 'atlas/pas_atlas.tsv', sep='\t', low_memory=False)
        assert f.atlas_pas_id.is_unique and annotation.atlas_pas_id.is_unique
        assert set(f.atlas_pas_id) == set(annotation.atlas_pas_id)
        annotation = annotation.set_index('atlas_pas_id').loc[f.atlas_pas_id].reset_index()
        assert np.array_equal(f.feature_index_1based, annotation.atlas_index_1based)
        features = f.merge(annotation, on='atlas_pas_id', validate='one_to_one', sort=False)
        for col in ['keep_for_pas_quant', 'keep_for_apa_primary', 'keep_for_sensitivity']:
            assert features[col].dtype == bool
        assert (~features.keep_for_apa_primary | features.keep_for_sensitivity).all()
        features.to_csv(dest / 'features.tsv.gz', sep='\t', compression='gzip', index=False)
        primary = features.keep_for_apa_primary
        assert features.loc[primary, 'primary_gene_id'].notna().all()
        genes = features[primary].groupby('primary_gene_id', sort=True).agg(
            gene_symbol=('primary_gene_symbol', 'first'), n_primary_atlas_pas=('atlas_pas_id', 'size'))
        genes.insert(0, 'gene_index_1based', np.arange(1, len(genes)+1))
        genes.reset_index().to_csv(dest / 'primary_genes.tsv.gz', sep='\t', compression='gzip', index=False)
        fingerprints[dataset] = sha256(source / 'matrix/features.tsv.gz', decompressed=True)
        feature_rows.append(dict(dataset_id=dataset, atlas_features=len(features), primary_pas=int(primary.sum()),
            sensitivity_pas=int(features.keep_for_sensitivity.sum()), primary_genes=len(genes),
            genes_with_multiple_primary_atlas_pas=int(genes.n_primary_atlas_pas.ge(2).sum()),
            compact_features_sha256=sha256(dest / 'features.tsv.gz'), source_atlas=str(source),
            source_annotation_sha256=sha256(source / 'atlas/pas_atlas.tsv')))
        print('PAS_FEATURE_AXIS_READY', dataset, len(features), int(primary.sum()), flush=True)
    bindings = []
    for r, pubrow, source, barcodes in selected:
        dest = OUT / 'libraries' / r.library_key
        dest.mkdir(parents=True, exist_ok=True)
        assert sha256(source / 'features.tsv.gz', True) == fingerprints[pubrow.dataset_id]
        current = cells[cells.library_key.eq(r.library_key)].copy()
        index = pd.Index(barcodes.original_barcode).get_indexer(current.barcode)
        assert (index >= 0).all() and len(current) == r.primary_cells
        current['local_state_id'] = current.original_state_id
        current['pas_library_id'] = pubrow.library_id
        current['pas_dataset_id'] = pubrow.dataset_id
        current['pas_column_index_1based'] = index + 1
        current['pas_cell_key'] = barcodes.cell_key.iloc[index].to_numpy()
        current['time_interpretation'] = current.time_kind
        current.to_csv(dest / 'cell_to_pas_and_metacell.tsv.gz', sep='\t', compression='gzip', index=False)
        expected_nodes = nodes[nodes.library_key.eq(r.library_key)].sort_values('metacell_id')
        expected_nodes.to_csv(dest / 'metacell_metadata.tsv.gz', sep='\t', compression='gzip', index=False)
        original = GEX / 'libraries' / r.library_key / 'cell_to_metacell.tsv.gz'
        bound.append(dict(library_key=r.library_key, library_id=pubrow.library_id, original_gex_library_id=r.library_id,
            dataset=r.dataset_id, study_id=r.study_id, pas_dataset_id=pubrow.dataset_id,
            cells=len(current), metacells=len(expected_nodes), source_pas_cells=len(barcodes),
            source_cells_outside_primary=len(barcodes)-len(current), atlas_features=int(pubrow.n_atlas_pas),
            source_projected_umi=int(pubrow.projected_umi), source_projected_nonzero=int(pubrow.projected_nonzero),
            source_projection=str(source), source_matrix_bytes=(source / 'feature_by_cell.mtx.gz').stat().st_size,
            source_features_sha256=sha256(source / 'features.tsv.gz'), source_cells_sha256=sha256(source / 'cells.tsv.gz'),
            source_audit_sha256=sha256(source / 'projection_audit.json'), frozen_assignments_sha256=sha256(original),
            selected_binding_sha256=sha256(dest / 'cell_to_pas_and_metacell.tsv.gz'), original_assignments=str(original),
            all_selected_barcodes_matched=True, complete_gex_pas_called_barcode_order_exact=True))
        bindings.append(current)
        print('PAS_BOUND', r.library_key, len(current), flush=True)
    pd.DataFrame(bound).to_csv(OUT / 'input_manifest.tsv', sep='\t', index=False)
    pd.DataFrame(availability).to_csv(OUT / 'availability.tsv', sep='\t', index=False)
    pd.DataFrame(feature_rows).to_csv(OUT / 'study_feature_summary.tsv', sep='\t', index=False)
    pd.concat(bindings, ignore_index=True).to_csv(OUT / 'all_cell_to_pas_and_metacell.tsv.gz', sep='\t', compression='gzip', index=False)
    dump(OUT / 'preflight_summary.json', dict(status='AVAILABLE_PAS_BOUND_PENDING_EXPLICIT',
        available_libraries=len(bound), pending_libraries=len(manifest)-len(bound),
        available_cells=sum(r['cells'] for r in bound), available_metacells=sum(r['metacells'] for r in bound),
        releases=releases, dataset_registry_aliases=DATASET_ALIASES, missing_selected_barcodes=0, membership_changed=False,
        gex_release_sha256=sha256(GEX / 'release/summary.json')))


if __name__ == '__main__':
    main()
