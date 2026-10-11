"""Apply native ASNR/junction rules to polyApipe GFF intervals, before library union.

This is postfiltering of frozen native counts, without UMI redistribution or rerunning
discovery. Original reproAPA junction evidence is reused for the matching library.
"""
from pathlib import Path
from collections import defaultdict
import ast
import bisect
import gzip
import hashlib
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = ROOT / 'notes/2026-10-02-tool-validation-coverage-audit'
PRIOR = AUDIT / 'junction-audit'
ATLAS = ROOT / 'notes/2026-10-02-scape-polyapipe-comparison/atlas/study_comparison/GSE118258'
CORE = ROOT / 'analysis_code/reproducibility_20260930/github/reproapa-analysis/modules/05_pas_catalog/code/apadetect_core'
sys.path.insert(0, str(CORE))
from junction_filter import junction_overlap_mask

# Load just the original ASNR functions; importing the calling CLI requires pysam.
source = ast.parse((CORE / 'apadetect.py').read_text())
names = {'canonical_chromosome', 'load_asnr_intervals', 'asnr_overlap_mask'}
functions = [node for node in source.body if isinstance(node, ast.FunctionDef) and node.name in names]
assert len(functions) == len(names)
exec(compile(ast.Module(body=functions, type_ignores=[]), str(CORE / 'apadetect.py'), 'exec'))
COORD = ['chrom', 'strand', 'cut_site']
CLASSES = np.array(['neither', 'ASNR_only', 'junction_only', 'both'])
inputs = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path, **kwargs):
    inputs[str(path.relative_to(ROOT))] = sha(path)
    return pd.read_csv(path, sep='\t', low_memory=False, **kwargs)


def save(frame, filename):
    frame.to_csv(HERE / filename, sep='\t', index=False, na_rep='NA',
                 compression={'method': 'gzip', 'mtime': 0} if filename.endswith('.gz') else None)


def nearest(query, target):
    """Both inputs must already use the same interbase coordinate convention."""
    result = np.full(len(query), np.nan)
    groups = {key: np.sort(group.cut_site.unique()) for key, group in target.groupby(['chrom', 'strand'])}
    for key, group in query.groupby(['chrom', 'strand']):
        ref = groups.get(key)
        if ref is None or not len(ref):
            continue
        points = group.cut_site.to_numpy()
        j = np.searchsorted(ref, points)
        result[group.index] = np.minimum(abs(points - ref[np.clip(j, 0, len(ref)-1)]),
                                         abs(points - ref[np.clip(j-1, 0, len(ref)-1)]))
    return result


def project_matrix(matrix, meta, flags, sites):
    """Remove a peak only in the library where it fails; retain other libraries."""
    blocks, order = [], []
    for library, group in flags.groupby('library_key', sort=False):
        rows = np.flatnonzero(meta.library_key.eq(library))
        assert not group.union_column.duplicated().any()
        retained = group[group.keep_both].union_column.to_numpy()
        keep = np.zeros(len(sites), dtype=bool)
        keep[retained] = True
        block = matrix[rows].copy()
        assert int(block.sum()) == int(group.integer_molecules.sum())
        block.data *= keep[block.indices]
        block.eliminate_zeros()
        assert int(block.sum()) == int(group.loc[group.keep_both, 'integer_molecules'].sum())
        blocks.append(block)
        order.extend(rows)
    assert sorted(order) == list(range(len(meta)))
    retained_columns = np.sort(flags.loc[flags.keep_both, 'union_column'].unique())
    joined = sparse.vstack(blocks, format='csr')[np.argsort(order)][:, retained_columns]
    return joined, retained_columns


def intent_checks():
    # One shared coordinate failing in library A must not erase library B's count.
    toy = sparse.csr_matrix([[10, 20], [30, 40]], dtype=np.int64)
    meta = pd.DataFrame(dict(library_key=['A', 'B']))
    flags = pd.DataFrame(dict(library_key=['A', 'A', 'B', 'B'], union_column=[0, 1, 0, 1],
                               integer_molecules=[10, 20, 30, 40], keep_both=[False, True, True, False]))
    output, columns = project_matrix(toy, meta, flags, pd.DataFrame(index=range(2)))
    assert columns.tolist() == [0, 1] and output.toarray().tolist() == [[0, 20], [30, 0]]
    # GFF [101,350] is BED [100,350), preserving both terminal bases.
    assert 350 - (101 - 1) == 250
    intervals = pd.DataFrame(dict(start=[100, 100], end=[350, 350], strand=['+', '-']))
    assert junction_overlap_mask(intervals, {'+': [(90, 20)], '-': [(0, 20)]}, 10, 100).all()
    assert not junction_overlap_mask(intervals, {'+': [(450, 20)], '-': [(360, 20)]}, 10, 100).any()


intent_checks()
previous = json.loads((PRIOR / 'verification.json').read_text())
assert previous['status'] == 'verified'
trace_path = PRIOR / 'poly-native-nearest-candidate.tsv.gz'
assert sha(trace_path) == previous['output_sha256'][trace_path.name]
trace = read(trace_path)
assert len(trace) == 395311 and int(trace.integer_molecules.sum()) == 572251017
assert not trace.duplicated(['library_id', 'native_column']).any()
snapshot = json.loads((HERE / 'original/source-manifest.json').read_text())
for record in snapshot:
    path = HERE / 'original' / record['archive_path']
    assert path.stat().st_size == record['size'] and sha(path) == record['sha256']
asnr_bytes = gzip.decompress((AUDIT / 'frozen-human-asnr.bed.gz').read_bytes())
assert hashlib.sha256(asnr_bytes).hexdigest() == '475f2472fc87345463fbeea44e56439ba81068a0adf6e8607cf894a49ab4b0f1'
asnr_path = HERE / 'original/frozen-human-asnr.bed'
asnr_path.write_bytes(asnr_bytes)
inputs[str(asnr_path.relative_to(ROOT))] = sha(asnr_path)

parts, audits = [], []
for library_id, group in trace.groupby('library_id', sort=True):
    accession = library_id.rsplit('_', 1)[-1]
    gff_path = HERE / 'original' / accession / 'joint_polyA_peaks.gff'
    gff = read(gff_path, comment='#', header=None,
               names=['chrom', 'source', 'type', 'start1', 'end', 'score', 'strand', 'phase', 'attributes'])
    gff['site_id'] = gff.attributes.str.extract(r'peak="([^"]+)"')[0]
    gff['peak_depth_gff'] = gff.attributes.str.extract(r'peakdepth="(\d+)"')[0].astype(int)
    gff['start'] = gff.start1 - 1
    gff['cut_site'] = np.where(gff.strand.eq('+'), gff.end, gff.start1)
    assert gff.site_id.is_unique and gff.start.ge(0).all() and gff.end.gt(gff.start).all()
    combined = group.merge(gff[COORD + ['site_id', 'start', 'end', 'peak_depth_gff']],
                          on=COORD + ['site_id'], validate='one_to_one')
    assert len(combined) == len(group) == len(gff)
    assert combined.peak_depth.eq(combined.peak_depth_gff).all()
    combined['postfilter_asnr'] = asnr_overlap_mask(combined, asnr_path)
    junction_path = PRIOR / 'original' / accession / 'junctions.bed'
    assert sha(junction_path) == previous['input_sha256'][str(junction_path.relative_to(ROOT))]
    junctions = read(junction_path, header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    combined['postfilter_junction'] = False
    for chrom, candidates in combined.groupby('chrom'):
        current = junctions[junctions.chrom.eq(chrom)]
        points = {strand: sorted(zip(part.start.astype(int), part.score.astype(int)))
                  for strand, part in current.groupby('strand')}
        combined.loc[candidates.index, 'postfilter_junction'] = junction_overlap_mask(candidates, points, 10, 100)
    combined['postfilter_class'] = CLASSES[combined.postfilter_asnr.astype(int) + 2 * combined.postfilter_junction.astype(int)]
    combined['keep_both'] = ~(combined.postfilter_asnr | combined.postfilter_junction)
    combined['proxy50_removed'] = combined.distance_nt.le(50) & combined.nearest_filter_class.isin(['ASNR_only', 'junction_only', 'both'])
    combined['proxy50_evaluable'] = combined.distance_nt.le(50) & combined.nearest_candidates.eq(1)
    # Independent direct interval tests on random and high-count peaks.
    merged_asnr = load_asnr_intervals(asnr_path)
    sample = pd.concat([combined.sample(100, random_state=20261003), combined.nlargest(100, 'integer_molecules')]).drop_duplicates('site_id')
    for row in sample.itertuples():
        starts, ends = merged_asnr.get((canonical_chromosome(row.chrom), row.strand), ([], []))
        a = np.asarray(starts); b = np.asarray(ends)
        assert row.postfilter_asnr == bool(((a < row.end) & (b > row.start)).any())
        js = junctions[junctions.chrom.eq(row.chrom) & junctions.strand.eq(row.strand)].start
        lower, upper = ((row.start - 10, row.end + 100) if row.strand == '+' else (row.start - 100, row.end + 10))
        assert row.postfilter_junction == bool(((js >= lower) & (js < upper)).any())
    parts.append(combined)
    audits.append(dict(library_id=library_id, library_key=group.library_key.iloc[0], peaks=len(combined),
                       peak_width_min=int((gff.end-gff.start).min()), peak_width_max=int((gff.end-gff.start).max()),
                       retained_peaks=int(combined.keep_both.sum()),
                       retained_molecules=int(combined.loc[combined.keep_both, 'integer_molecules'].sum()),
                       independent_interval_checks=len(sample)))
    print(json.dumps(audits[-1]), flush=True)

flags = pd.concat(parts, ignore_index=True)
save(flags, 'native-peak-filter-audit.tsv.gz')
counts = flags.groupby(['library_key', 'postfilter_class']).agg(native_peak_records=('site_id', 'size'),
    integer_molecules=('integer_molecules', 'sum'), native_peak_depth=('peak_depth', 'sum')).reset_index()
save(counts, 'filter-classes-by-library.tsv')
save(counts.groupby('postfilter_class', as_index=False)[['native_peak_records', 'integer_molecules', 'native_peak_depth']].sum(), 'filter-classes-pooled.tsv')
primary = read(AUDIT / 'repro-primary-support.tsv.gz')
scenarios = {'before': np.ones(len(flags), bool), 'ASNR_only': ~flags.postfilter_asnr,
             'junction_only': ~flags.postfilter_junction, 'both_native_GFF': flags.keep_both,
             'both_proxy50': ~flags.proxy50_removed}
summary, distances, library_summaries = [], [], []
old_repro_anchors = set(primary.anchor_id.dropna())
for label, mask in scenarios.items():
    kept = flags.loc[mask]
    union = kept.groupby(COORD, as_index=False).agg(integer_molecules=('integer_molecules', 'sum'),
        native_libraries=('library_key', 'nunique'), native_peak_records=('site_id', 'size'),
        native_peak_depth=('peak_depth', 'sum'), old_anchor_id=('anchor_id', 'first'),
        nearest_repro_primary_distance_nt=('nearest_repro_primary_distance_nt', 'first'))
    missing = ~union.nearest_repro_primary_distance_nt.le(50)
    interbase = union[COORD].copy()
    interbase['cut_site'] -= interbase.strand.eq('-').astype(int)
    reverse_distances = nearest(primary, interbase)
    poly_anchors = set(union.old_anchor_id.dropna())
    shared = poly_anchors & old_repro_anchors
    record = dict(scenario=label, native_peak_records=len(kept), union_sites=len(union),
        counted_union_sites=int(union.integer_molecules.gt(0).sum()), integer_molecules=int(kept.integer_molecules.sum()),
        native_peak_depth=int(kept.peak_depth.sum()), sites_without_repro_primary50=int(missing.sum()),
        molecules_without_repro_primary50=int(union.loc[missing, 'integer_molecules'].sum()),
        repro_primary_sites_with_poly50=int((reverse_distances <= 50).sum()),
        repro_primary_molecules_with_poly50=int(primary.loc[reverse_distances <= 50, 'molecules'].sum()),
        old_strict_annotation_ids=len(poly_anchors), old_strict_shared_annotation_ids=len(shared),
        old_strict_poly_only_annotation_ids=len(poly_anchors-shared),
        old_strict_annotation_molecules=int(union.loc[union.old_anchor_id.notna(), 'integer_molecules'].sum()),
        old_strict_shared_annotation_molecules=int(union.loc[union.old_anchor_id.isin(shared), 'integer_molecules'].sum()))
    summary.append(record)
    for library, part in kept.groupby('library_key'):
        library_summaries.append(dict(scenario=label, library_key=library, native_peak_records=len(part),
                                      integer_molecules=int(part.integer_molecules.sum())))
    for tolerance in [25, 50, 100, 250, 500]:
        near = union.nearest_repro_primary_distance_nt.le(tolerance)
        distances.append(dict(scenario=label, tolerance_nt=tolerance, union_sites=len(union),
                              integer_molecules=int(union.integer_molecules.sum()), near_sites=int(near.sum()),
                              near_molecules=int(union.loc[near, 'integer_molecules'].sum())))
    if label in ['both_native_GFF', 'both_proxy50']:
        save(union, label + '-coordinate-union.tsv.gz')
save(pd.DataFrame(summary), 'before-after-summary.tsv')
save(pd.DataFrame(distances), 'distance-sensitivity.tsv')
save(pd.DataFrame(library_summaries), 'before-after-by-library.tsv')

sites = read(ATLAS / 'complete_coordinate_union.tsv.gz')
meta = read(ATLAS / 'metacell_index.tsv.gz')
matrix_path = ATLAS / 'metacell_by_union_site_counts.npz'
inputs[str(matrix_path.relative_to(ROOT))] = sha(matrix_path)
matrix = sparse.load_npz(matrix_path).tocsr()
assert matrix.shape == (583, 273066)
assert np.array_equal(matrix.sum(0).A1, sites.integer_molecules)
filtered, columns = project_matrix(matrix, meta, flags, sites)
assert np.issubdtype(filtered.dtype, np.integer) and (filtered.data >= 0).all()
assert int(filtered.sum()) == int(flags.loc[flags.keep_both, 'integer_molecules'].sum())
output_sites = sites.iloc[columns][COORD + ['site_id', 'union_column']].copy().reset_index(drop=True)
output_sites = output_sites.rename(columns={'union_column': 'original_union_column'})
output_sites['matrix_column'] = np.arange(len(output_sites))
output_sites['integer_molecules'] = filtered.sum(0).A1
expected = flags[flags.keep_both].groupby('union_column').integer_molecules.sum().reindex(columns)
assert np.array_equal(output_sites.integer_molecules, expected)
sparse.save_npz(HERE / 'postfilter-metacell-by-site-counts.npz', filtered)
save(output_sites, 'postfilter-matrix-sites.tsv.gz')
save(meta, 'postfilter-metacell-index.tsv.gz')
libraries = pd.Index(sorted(flags.library_key.unique()))
retained = flags[flags.keep_both]
available = sparse.csr_matrix((np.ones(len(retained), dtype=np.int8),
    (libraries.get_indexer(retained.library_key), pd.Index(columns).get_indexer(retained.union_column))),
    shape=(len(libraries), len(columns)))
assert available.data.min() == available.data.max() == 1
assert int(available.sum()) == len(retained) and (available.sum(0).A1 > 0).all()
nonzero = filtered.tocoo()
library_row = libraries.get_indexer(meta.library_key)
assert available.toarray()[library_row[nonzero.row], nonzero.col].all()
sparse.save_npz(HERE / 'postfilter-library-site-availability.npz', available)
save(pd.DataFrame(dict(library_row=np.arange(len(libraries)), library_key=libraries)), 'postfilter-library-index.tsv')
for name in ['analyze.py', 'original/source-manifest.json']:
    inputs[str((HERE / name).relative_to(ROOT))] = sha(HERE / name)
for name in ['apadetect.py', 'junction_filter.py']:
    inputs[str((CORE / name).relative_to(ROOT))] = sha(CORE / name)
verification = dict(status='verified', libraries=audits, native_peak_records=len(flags),
    rules=dict(ASNR='exact native same-strand interval overlap', junction_left=10, junction_right=100,
               junction_source='matching-library original retained junctions.bed', interval='native GFF [start1,end] converted to BED [start1-1,end)'),
    matrix_shape=list(filtered.shape), availability_shape=list(available.shape), retained_integer_molecules=int(filtered.sum()),
    matrix_integer_and_count_conservation=True, cross_library_filtering_intent_test=True,
    input_sha256=inputs, output_sha256={p.name: sha(p) for p in HERE.iterdir() if p.suffix in ['.gz', '.npz', '.tsv']},
    limitations=['Postfiltering original counts, without rerunning discovery, featureCounts or UMI assignment.',
                 'PolyApipe counting intervals differ from reproAPA derivative candidates; this is a specified harmonization sensitivity.',
                 'Original historical junction evidence is reused, not extracted anew from frozen QC cells.',
                 'Old strict annotation-ID summaries retain the known unextended-UTR restriction.',
                 'RNA/polyA aligned-block coverage and APA effects are not recomputed in this postfilter analysis.',
                 'Matrix zeros may mean a peak failed filtering or was not called in that library; consult the availability mask before APA effects.',
                 'Proxy50 retains unmatched/ambiguous peaks as unevaluated, not as verified filter passes.'])
(HERE / 'verification.json').write_text(json.dumps(verification, indent=2) + '\n')
print(pd.DataFrame(summary).to_string(index=False))
