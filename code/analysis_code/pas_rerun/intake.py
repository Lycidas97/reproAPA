#!/usr/bin/env python3
"""Validate existing PAS tables and sparse counts without changing GEX membership."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def table(path, columns):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        require(set(columns) <= set(reader.fieldnames or []), f'{path}: missing columns {set(columns)-set(reader.fieldnames or [])}')
        rows = list(reader)
    require(bool(rows), f'{path}: empty table')
    return rows


def unique(rows, key, context):
    values = [r[key] for r in rows]
    require(all(values) and len(set(values)) == len(values), f'{context}: blank/duplicate {key}')
    return values


def boolean(value):
    require(value in {'True', 'False', 'true', 'false', '0', '1'}, f'Invalid Boolean: {value!r}')
    return value in {'True', 'true', '1'}


def artifact(item, base):
    require(isinstance(item, dict) and item.get('path') and item.get('sha256'), 'Unbound input path/sha256')
    p = (base / item['path']).resolve()
    require(p.is_file(), f'Missing input: {p}')
    require(digest(p) == item['sha256'], f'Input hash mismatch: {p}')
    return p


def load_matrix(path, orientation, n_cells, n_pas):
    import numpy as np
    from scipy import sparse
    from scipy.io import mmread
    require(orientation in {'cell_by_pas', 'pas_by_cell'}, 'Matrix orientation must be explicit')
    x = sparse.load_npz(path) if str(path).endswith('.npz') else mmread(path, spmatrix=True)
    require(sparse.issparse(x), f'Counts must be sparse: {path}')
    require(np.issubdtype(x.dtype, np.integer) and np.all(x.data >= 0), f'Non-integer/negative PAS counts: {path}')
    if orientation == 'pas_by_cell':
        x = x.T
    require(x.shape == (n_cells, n_pas), f'Matrix/axis dimension mismatch: {path}: {x.shape} != {(n_cells,n_pas)}')
    # Duplicate sparse coordinates are storage details, not extra molecules.
    # Sum them exactly as the existing metacell aggregation does.
    x = sparse.csr_matrix(x, dtype=np.int64)
    x.sum_duplicates()
    x.eliminate_zeros()
    require(np.all(x.data >= 0), f'Integer overflow after coalescing: {path}')
    return x


def validate(manifest_path, *, scope_path=None, full=True):
    p = Path(manifest_path).resolve()
    m = json.loads(p.read_text())
    require(m.get('schema_version') == 1, 'Unsupported intake schema')
    require(m.get('preserve_expression_framework') is True, 'This rerun requires the frozen expression framework')
    require(m.get('pas_release_id') and m.get('caller') and m.get('quantification'), 'Missing PAS release/caller/quantification provenance')
    require(m.get('cut_coordinate_system') in {'0-based-genomic-base', '0-based-interbase'}, 'cut_site must explicitly declare genomic-base or interbase coordinates (not inferred)')
    for key in ['version', 'command', 'deduplication_unit']:
        require(m['quantification'].get(key), f'Missing quantification.{key}')
    scope = table(scope_path or HERE/'scope.tsv', ['study', 'library_id', 'species', 'expression_atlas_cells'])
    expected = {r['library_id']: r for r in scope}
    libs = m.get('libraries', [])
    unique(libs, 'library_id', 'intake libraries')
    require(set(r['library_id'] for r in libs) == set(expected), 'Library scope differs from frozen expression atlas; omissions are not allowed')
    if scope_path is None:
        frozen = table(HERE/'frozen-membership.tsv', ['library_id', 'role', 'sha256'])
        hashes = {(r['library_id'], r['role']): r['sha256'] for r in frozen}
        for lib in libs:
            for role in ['membership', 'metacells']:
                require(lib[role].get('sha256') == hashes[(lib['library_id'], role)], f'Frozen {role} changed: {lib["library_id"]}')
    studies = m.get('studies', [])
    unique(studies, 'study_id', 'studies')
    require(set(r['study_id'] for r in studies) == {r['study'] for r in scope}, 'Study scope differs from frozen atlas')
    catalog, files = {}, {}

    def get(item):
        path = artifact(item, p.parent)
        files[str(path)] = item['sha256']
        return path

    required_resources = {'gex', 'state_annotations', 'embeddings', 'library_design', 'contrast_design', 'genome_fasta', 'gene_annotation', 'external_pas', 'orthologs', 'sequence_resources', 'recount_inputs', 'external_validation', 'legacy_rosters'}
    resources = m.get('resources', {})
    require(required_resources <= set(resources), f'Missing resource roles: {sorted(required_resources-set(resources))}')
    for role, items in resources.items():
        require(isinstance(items, list) and items, f'Unbound resource role: {role}')
        for item in items:
            get(item)
    for study in studies:
        require(study.get('assembly') in {'GRCh38', 'GRCm39'}, f'Undeclared/nonproject assembly: {study["study_id"]}')
        f = table(get(study['features']), ['atlas_pas_id', 'feature_index_1based', 'chrom', 'cut_site', 'strand', 'primary_gene_id', 'keep_for_apa_primary'])
        unique(f, 'atlas_pas_id', study['study_id'])
        require([int(r['feature_index_1based']) for r in f] == list(range(1,len(f)+1)), 'PAS feature order must match 1-based matrix order')
        for row in f:
            require(row['strand'] in {'+', '-'} and row['chrom'] and int(row['cut_site']) >= 0, f'Invalid PAS coordinate: {row}')
            if boolean(row['keep_for_apa_primary']):
                require(row['primary_gene_id'] not in {'', 'NA', 'nan'}, f'Primary PAS lacks unique gene: {row["atlas_pas_id"]}')
        catalog[study['study_id']] = f
    summaries, global_cells, global_mcs = [], set(), set()
    for lib in libs:
        lid = lib['library_id']; reference = expected[lid]
        require(lib['study_id'] == reference['study'], f'Library/study mismatch: {lid}')
        feature = catalog[lib['study_id']]
        study = next(s for s in studies if s['study_id'] == lib['study_id'])
        require(study['assembly'] == ('GRCh38' if reference['species'] == 'human' else 'GRCm39'), f'Assembly/species mismatch: {lid}')
        paths = {k: get(lib[k]) for k in ['counts', 'cells', 'feature_status', 'membership', 'metacells']}
        cells = table(paths['cells'], ['cell_index_1based', 'original_barcode', 'cell_key'])
        barcodes = unique(cells, 'original_barcode', lid)
        unique(cells, 'cell_key', lid)
        require([int(r['cell_index_1based']) for r in cells] == list(range(1,len(cells)+1)), f'Cell row order mismatch: {lid}')
        members = table(paths['membership'], ['cell_id', 'barcode', 'metacell_id', 'library_key'])
        ids = unique(members, 'cell_id', lid); unique(members, 'barcode', lid)
        require(not global_cells.intersection(ids), f'Cell IDs reused across libraries: {lid}'); global_cells.update(ids)
        require(len(members) == int(reference['expression_atlas_cells']), f'Frozen cell count changed: {lid}')
        require({r['barcode'] for r in members} <= set(barcodes), f'Frozen cells missing from PAS matrix: {lid}; zero cells must retain rows')
        require({r['library_key'] for r in members} == {lib['library_key']}, f'Membership library key mismatch: {lid}')
        mc = table(paths['metacells'], ['metacell_id', 'n_cells'])
        mids = unique(mc, 'metacell_id', lid)
        require(not global_mcs.intersection(mids), f'Metacell IDs reused across libraries: {lid}'); global_mcs.update(mids)
        sizes = {}
        for r in members:
            sizes[r['metacell_id']] = sizes.get(r['metacell_id'], 0)+1
        require(sizes == {r['metacell_id']: int(r['n_cells']) for r in mc}, f'Frozen metacell membership/size mismatch: {lid}')
        status = table(paths['feature_status'], ['atlas_pas_id', 'terminal_status', 'projected_umi'])
        unique(status, 'atlas_pas_id', lid)
        require(set(r['atlas_pas_id'] for r in status) == set(r['atlas_pas_id'] for r in feature), f'Call-status PAS axis mismatch: {lid}')
        lookup = {r['atlas_pas_id']: r for r in status}
        totals = [int(lookup[r['atlas_pas_id']]['projected_umi']) for r in feature]
        for r in status:
            require(r['terminal_status'] in {'observed', 'not_observed'} and int(r['projected_umi']) >= 0, f'Invalid call status: {lid}')
            require(r['terminal_status'] != 'not_observed' or int(r['projected_umi']) == 0, f'Projected absent PAS has nonzero UMI: {lid}')
        if full:
            x = load_matrix(paths['counts'], lib['orientation'], len(cells), len(feature))
            require(x.sum(axis=0).A1.tolist() == totals, f'PAS count conservation failed: {lid}')
        summaries.append(dict(library_id=lid, study_id=lib['study_id'], cells=len(cells), frozen_cells=len(members), metacells=len(mc), pas=len(feature), umi=sum(totals)))
    return dict(status='validated' if full else 'metadata_only', pas_release_id=m['pas_release_id'], files=files, libraries=summaries, frozen_cells=len(global_cells), frozen_metacells=len(global_mcs))
