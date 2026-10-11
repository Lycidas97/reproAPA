"""Bind frozen GEX membership to the released v7 PAS axes without filtering."""
import hashlib
import gzip
import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd


ATLAS = Path('/root/nfsdata/apadetect/data/reprogramming_differentiation_apa/analyses/apadetect_reproapa_study_atlas_v7_merge_first_6aabc22_2026-09-05')
OUT = Path('results_metacell_pas_v7_v1')
BASE = Path('results_metacell_on_working_baseline_v1')
MANIFEST = Path('results_full_mouse_projection_v1/full_library_manifest.tsv')


def sha256(path, decompressed=False):
    h = hashlib.sha256()
    opener = gzip.open if decompressed else open
    with opener(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def main():
    OUT.mkdir(exist_ok=True)
    (OUT / 'provenance').mkdir(exist_ok=True)
    release = json.loads((ATLAS / 'release/release_complete.json').read_text())
    assert release['status'] == 'passed' and release['all_projection_count_conservation_passed']
    assert sha256(ATLAS / 'release/library_projection_summary.tsv') == '52ece636206115d353cd5c610cd797bd73671b24db7e218d190f551f76e7bd4e'
    published = pd.read_csv(ATLAS / 'release/library_projection_summary.tsv', sep='\t').set_index('library_id')
    manifest = pd.read_csv(MANIFEST, sep='\t')
    cells = pd.read_csv(BASE / 'cells_with_metacell_on_baseline.tsv.gz', sep='\t')
    nodes = pd.read_csv(BASE / 'metacell_coordinates.tsv.gz', sep='\t')
    assert manifest.library_key.is_unique and published.index.is_unique
    assert cells.cell_id.is_unique and nodes.metacell_id.is_unique
    assert (len(manifest), len(cells), len(nodes)) == (102, 336087, 6772)
    studies = sorted(published.loc[manifest.library_id, 'dataset_id'].unique())
    study_rows = []
    feature_fingerprints = {}
    annotation_cols = ['atlas_pas_id', 'atlas_index_1based', 'species', 'genome_build', 'chrom',
        'peak_start', 'peak_end', 'cut_site', 'strand', 'primary_gene_id', 'primary_gene_symbol',
        'primary_gene_biotype', 'primary_transcript_id', 'primary_transcript_biotype',
        'primary_region_class', 'primary_assignment_status', 'primary_assignment_confidence',
        'multi_gene_ambiguous', 'transcript_model_conflict', 'candidate_gene_ids',
        'primary_model_transcript_end', 'primary_model_3p_distance', 'primary_strand_relation',
        'atlas_status', 'confidence_tier', 'keep_for_pas_quant', 'keep_for_apa_primary',
        'keep_for_sensitivity', 'internal_priming_flag', 'sequence_qc_status',
        'external_support_50bp', 'external_support_100bp', 'n_member_pas', 'n_libraries',
        'member_cut_site_min', 'member_cut_site_max', 'cluster_width_bp', 'tier_reason', 'failure_reason']
    for study in studies:
        source = ATLAS / 'studies' / study
        dest = OUT / 'studies' / study
        dest.mkdir(parents=True, exist_ok=True)
        features = pd.read_csv(source / 'matrix/features.tsv.gz', sep='\t')
        anno_path = source / 'atlas/pas_atlas.tsv'
        annotation = pd.read_csv(anno_path, sep='\t', usecols=annotation_cols, low_memory=False)
        assert features.atlas_pas_id.is_unique and annotation.atlas_pas_id.is_unique
        assert set(features.atlas_pas_id) == set(annotation.atlas_pas_id)
        annotation = annotation.set_index('atlas_pas_id').loc[features.atlas_pas_id].reset_index()
        assert np.array_equal(features.feature_index_1based, np.arange(1, len(features) + 1))
        assert np.array_equal(annotation.atlas_index_1based, features.feature_index_1based)
        for col in ['keep_for_pas_quant', 'keep_for_apa_primary', 'keep_for_sensitivity']:
            assert annotation[col].dtype == bool
        assert (~annotation.keep_for_apa_primary | annotation.keep_for_sensitivity).all()
        features = features.merge(annotation, on='atlas_pas_id', validate='one_to_one', sort=False)
        features.to_csv(dest / 'features.tsv.gz', sep='\t', index=False)
        primary = features.keep_for_apa_primary & features.primary_gene_id.notna()
        genes = features[primary].groupby('primary_gene_id', sort=True).agg(
            gene_symbol=('primary_gene_symbol', 'first'), n_primary_atlas_pas=('atlas_pas_id', 'size'))
        genes.insert(0, 'gene_index_1based', np.arange(1, len(genes) + 1))
        genes.reset_index().to_csv(dest / 'primary_genes.tsv.gz', sep='\t', index=False)
        feature_fingerprints[study] = sha256(source / 'matrix/features.tsv.gz', decompressed=True)
        sr = dict(dataset_id=study, atlas_features=len(features), primary_pas=int(primary.sum()),
            sensitivity_pas=int(features.keep_for_sensitivity.sum()), primary_genes=len(genes),
            genes_with_multiple_primary_atlas_pas=int(genes.n_primary_atlas_pas.ge(2).sum()),
            source_atlas=str(anno_path), source_annotation_sha256=sha256(anno_path),
            source_features_sha256=sha256(source / 'matrix/features.tsv.gz'),
            source_features_content_sha256=feature_fingerprints[study], compact_features_sha256=sha256(dest / 'features.tsv.gz'))
        study_rows.append(sr)
        print(json.dumps(sr), flush=True)
    rows, bindings = [], []
    for r in manifest.itertuples():
        pub = published.loc[r.library_id]
        assert pub.status == 'passed' and pub.barcode_order_passed and pub.count_conservation_passed
        source = ATLAS / 'studies' / pub.dataset_id / 'libraries' / r.library_id
        dest = OUT / 'libraries' / r.library_key
        dest.mkdir(parents=True, exist_ok=True)
        required = ['feature_by_cell.mtx.gz', 'features.tsv.gz', 'cells.tsv.gz',
                    'projection_audit.json', 'projection_complete.json', 'feature_by_library_status.tsv.gz']
        assert all((source / f).is_file() for f in required), str(source)
        assert sha256(source / 'features.tsv.gz', decompressed=True) == feature_fingerprints[pub.dataset_id]
        audit = json.loads((source / 'projection_audit.json').read_text())
        assert audit['status'] == 'passed' and audit['library_id'] == r.library_id
        assert audit['dataset_id'] == pub.dataset_id
        original = Path('results_metacell_pilot_gsm3195672_v1/target_50' if r.library_key == 'GSM3195672'
            else f'results_mouse_metacell_full_v1/library/{r.library_key}') / 'cell_to_metacell.tsv.gz'
        member = pd.read_csv(original, sep='\t', usecols=['cell_id', 'barcode', 'metacell_id', 'library_id', 'library_key', 'local_state_id'])
        current = cells[cells.library_key.eq(r.library_key)].copy()
        assert len(current) == len(member) == int(r.n_cells)
        assert member.cell_id.is_unique and member.barcode.is_unique
        assert set(member.cell_id) == set(current.cell_id)
        member = member.set_index('cell_id').loc[current.cell_id]
        assert np.array_equal(member.metacell_id, current.metacell_id)
        assert np.array_equal(member.local_state_id, current.local_state_id)
        assert member.library_id.eq(r.library_id).all()
        barcodes = pd.read_csv(source / 'cells.tsv.gz', sep='\t')
        assert barcodes.original_barcode.is_unique and barcodes.cell_key.is_unique
        assert barcodes.library_id.eq(r.library_id).all()
        assert np.array_equal(barcodes.cell_index_1based, np.arange(1, len(barcodes) + 1))
        index = pd.Index(barcodes.original_barcode).get_indexer(member.barcode)
        assert (index >= 0).all(), f'Missing selected barcodes: {r.library_key}'
        assert len(barcodes) == int(pub.n_cells) == audit['n_cells']
        current['barcode'] = member.barcode.to_numpy()
        current['library_id'] = r.library_id
        current['pas_dataset_id'] = pub.dataset_id
        current['pas_column_index_1based'] = index + 1
        current['pas_cell_key'] = barcodes.cell_key.iloc[index].to_numpy()
        current['time_interpretation'] = ('D3_capture_pool_reference_identity_unresolved'
            if r.dataset.startswith('gse171127_d3') else 'observed_sampling_day')
        current.to_csv(dest / 'cell_to_pas_and_metacell.tsv.gz', sep='\t', index=False)
        expected_nodes = nodes[nodes.library_key.eq(r.library_key)].sort_values('metacell_id')
        assert current.metacell_id.value_counts().reindex(expected_nodes.metacell_id).eq(expected_nodes.n_cells.to_numpy()).all()
        expected_nodes.to_csv(dest / 'metacell_metadata.tsv.gz', sep='\t', index=False)
        row = dict(library_key=r.library_key, library_id=r.library_id, dataset=r.dataset,
            study_id=r.study_id, pas_dataset_id=pub.dataset_id, cells=len(current),
            metacells=len(expected_nodes), source_pas_cells=len(barcodes),
            source_cells_outside_primary=len(barcodes) - len(current), atlas_features=int(pub.n_atlas_pas),
            source_projected_umi=int(pub.projected_umi), source_projected_nonzero=int(pub.projected_nonzero),
            source_projection=str(source), source_matrix_bytes=(source / 'feature_by_cell.mtx.gz').stat().st_size,
            source_features_sha256=sha256(source / 'features.tsv.gz'),
            source_cells_sha256=sha256(source / 'cells.tsv.gz'), source_audit_sha256=sha256(source / 'projection_audit.json'),
            frozen_assignments_sha256=sha256(original), selected_binding_sha256=sha256(dest / 'cell_to_pas_and_metacell.tsv.gz'),
            original_assignments=str(original), all_selected_barcodes_matched=True)
        rows.append(row); bindings.append(current)
        print(f'BOUND {len(rows)}/102 {r.library_key}: {len(current)} cells', flush=True)
    pd.DataFrame(rows).to_csv(OUT / 'input_manifest.tsv', sep='\t', index=False)
    pd.DataFrame(study_rows).to_csv(OUT / 'study_feature_summary.tsv', sep='\t', index=False)
    pd.concat(bindings, ignore_index=True).to_csv(OUT / 'all_cell_to_pas_and_metacell.tsv.gz', sep='\t', index=False)
    assert sum(r['cells'] for r in rows) == 336087 and sum(r['metacells'] for r in rows) == 6772
    provenance = dict(checked_at_utc=datetime.now(timezone.utc).isoformat(), atlas_root=str(ATLAS),
        atlas_release=release, source_manifest_sha256=sha256(MANIFEST),
        frozen_cells_sha256=sha256(BASE / 'cells_with_metacell_on_baseline.tsv.gz'),
        frozen_metacells_sha256=sha256(BASE / 'metacell_coordinates.tsv.gz'),
        annotation_policy='All coordinate-atlas PAS retained; v7 primary and sensitivity masks are annotations, not destructive filters',
        statistical_units='Frozen within-library metacells; not assumed independent biological replicates',
        cross_study_pas_axis_merged=False, d3_reference_identity_resolved=False)
    dump(OUT / 'provenance/input_contract.json', provenance)
    dump(OUT / 'preflight_summary.json', dict(status='ALL_SELECTED_BARCODES_BOUND', libraries=102,
        cells=336087, metacells=6772, studies=len(studies), source_cells_outside_primary=sum(r['source_cells_outside_primary'] for r in rows),
        missing_selected_cells=0, pas_matrices_read=False))


if __name__ == '__main__':
    main()
