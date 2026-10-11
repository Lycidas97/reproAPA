"""Read-only input extraction. No marker classification or APA effect is computed."""
from pathlib import Path
import gzip
import hashlib
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse

P_SHARED = ['DPPA2', 'DPPA4', 'SALL4', 'TDGF1']
P_LEGACY = {
    'mouse': ['DPPA4', 'ESRRB', 'ZFP42', 'DPPA5A', 'TDGF1'],
    'human': ['DPPA4', 'DNMT3B', 'ESRG', 'L1TD1', 'TDGF1']
}
ECM = ['COL1A1', 'COL3A1', 'DCN', 'LUM', 'FN1']
STUDIES = {
    'GSE103221': ('mouse', ['CIP']),
    'GSE115943': ('mouse', ['WOT']),
    'GSE118258': ('human', ['GSE118258']),
    'GSE150311': ('human', ['GSE150311']),
    'GSE178325': ('human', ['178325_0618', '178325_0809']),
    'GSE221739': ('human', ['GSE221739'])
}
arrays = {}
provenance = {}
mapping = []
audit = []


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path):
    provenance[str(path)] = sha(path)
    return pd.read_csv(path, sep='\t', low_memory=False)


def records(frame):
    return np.array(frame.to_json(orient='records', double_precision=15))


annotation = {}
for species in ['human', 'mouse']:
    entry = CONFIG['resources'][species + '_annotation_gtf']
    path = Path(entry['path'])
    provenance[str(path)] = sha(path)
    assert provenance[str(path)] == entry['sha256']
    names = {}
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2] == 'gene':
                attributes = dict(re.findall(r'(\w+) "([^\"]*)"', fields[8]))
                if 'gene_name' in attributes:
                    names[attributes['gene_id'].split('.')[0]] = attributes['gene_name']
    annotation[species] = names

states_entry = CONFIG['resources']['frozen_metacell_states']
frozen_states = read(states_entry['path'])
assert provenance[states_entry['path']] == states_entry['sha256']
assert frozen_states.metacell_id.is_unique
arrays['frozen_states_json'] = records(frozen_states)
ortholog_entry = CONFIG['resources']['ortholog']
orthologs = read(ortholog_entry['path'])
assert provenance[ortholog_entry['path']] == ortholog_entry['sha256']
arrays['orthologs_json'] = records(orthologs)

for sid, (species, wanted_routes) in STUDIES.items():
    folder = 'WOT_GSE115943' if sid == 'GSE115943' else sid
    root = Path(CONFIG['native_run']) / 'studies' / folder
    gex_root = Path(CONFIG['frozen_gex']) / folder
    meta = read(root / 'metacells.tsv.gz')
    meta['matrix_row'] = np.arange(len(meta))
    selected = read(root / 'selected_metacells.tsv.gz')
    assert meta.metacell_id.is_unique and selected.metacell_id.is_unique
    assert np.array_equal(meta.metacell_id.iloc[selected.matrix_row], selected.metacell_id)
    assert selected.groupby('library_key').route.nunique().max() == 1
    by_library = selected.drop_duplicates('library_key').set_index('library_key').route
    meta['route'] = meta.library_key.map(by_library)
    # Route is library metadata, never inferred from marker score or APA.
    if sid == 'GSE103221':
        expected = meta.dataset.map({'gse103221_cip': 'CIP', 'gse103221_osk': 'OSK'})
        assert expected.notna().all()
        assert meta.route.dropna().eq(expected[meta.route.notna()]).all()
        meta['route'] = expected
    elif sid == 'GSE115943':
        assert meta.dataset.eq('wot_common_to_2i').all()
        meta['route'] = 'WOT'
    elif sid != 'GSE178325':
        assert meta.route.dropna().eq(sid).all()
        meta['route'] = sid
    assert meta.route.notna().all(), (sid, meta.loc[meta.route.isna(), 'library_key'].unique())
    meta['in_planned_route'] = meta.route.isin(wanted_routes)
    meta['original_selected'] = meta.metacell_id.isin(selected.metacell_id)
    selected_by_id = selected.set_index('metacell_id')
    meta['original_node'] = meta.metacell_id.map(selected_by_id.node)
    meta['original_purity'] = meta.metacell_id.map(selected_by_id.purity)
    arrays[sid + '_selected_json'] = records(selected)
    arrays[sid + '_sampling_json'] = records(read(root / 'sampling_design.tsv'))
    arrays[sid + '_endpoint_json'] = records(read(root / 'endpoint_library_design.tsv'))
    gm = read(gex_root / 'metacells.tsv.gz')
    genes = read(gex_root / 'genes.tsv.gz')
    assert gm.metacell_id.is_unique
    ix = pd.Index(gm.metacell_id).get_indexer(meta.metacell_id)
    assert (ix >= 0).all()
    assert np.array_equal(gm.n_cells.iloc[ix], meta.n_cells)
    assert np.array_equal(gm.library_key.iloc[ix], meta.library_key)
    ids = genes.gene_id.astype(str).str.split('.').str[0]
    assert ids.is_unique
    symbols = genes.gene_symbol.fillna('').astype(str).copy()
    replacement = ids.map(annotation[species])
    placeholders = symbols.eq('') | symbols.str.split('.').str[0].eq(ids)
    symbols.loc[placeholders & replacement.notna()] = replacement[placeholders & replacement.notna()]
    markers = list(dict.fromkeys(P_SHARED + P_LEGACY[species] + ECM))
    positions = []
    for marker in markers:
        matches = np.flatnonzero(symbols.str.upper().eq(marker))
        mapping.append(dict(study_id=sid, species=species, marker=marker,
                            status='unique' if len(matches) == 1 else 'missing' if not len(matches) else 'ambiguous',
                            gene_id=';'.join(ids.iloc[matches]), gene_symbol=';'.join(symbols.iloc[matches])))
        assert len(matches) == 1, (sid, marker, len(matches))
        positions.append(matches[0])
    xp = gex_root / 'metacell_by_gene_counts.npz'
    provenance[str(xp)] = sha(xp)
    x = sparse.load_npz(xp).tocsr()[ix]
    totals = gm.total_gex.to_numpy()[ix]
    assert np.issubdtype(x.dtype, np.integer) and (x.data >= 0).all()
    assert np.array_equal(np.asarray(x.sum(1)).ravel(), totals)
    arrays[sid + '_marker_counts'] = x[:, positions].toarray()
    arrays[sid + '_markers'] = np.array(markers)
    arrays[sid + '_marker_gene_ids'] = ids.iloc[positions].to_numpy(str)
    arrays[sid + '_total_gex'] = totals
    meta['total_gex'] = totals
    del x
    features = read(root / 'features.tsv.gz')
    features = features[features.keep_for_apa_primary.astype(str).str.lower().eq('true') & features.primary_gene_id.notna()].copy()
    assert features.feature_index_1based.is_unique
    signed = features.cut_site * np.where(features.strand.eq('+'), 1, -1)
    shape = features.groupby('primary_gene_id').agg(n=('cut_site', 'nunique'), c=('chrom', 'nunique'), s=('strand', 'nunique'))
    features['position'] = (signed - signed.groupby(features.primary_gene_id).transform('min')) / (
        signed.groupby(features.primary_gene_id).transform('max') - signed.groupby(features.primary_gene_id).transform('min'))
    features['coordinate_eligible'] = features.primary_gene_id.isin(shape.query('n>=2 and c==1 and s==1').index)
    features = features.sort_values(['primary_gene_id', 'feature_index_1based']).reset_index(drop=True)
    pp = root / 'metacell_by_pas_counts.npz'
    provenance[str(pp)] = sha(pp)
    pas = sparse.load_npz(pp).tocsr()
    assert pas.shape[0] == len(meta)
    assert np.issubdtype(pas.dtype, np.integer) and (pas.data >= 0).all()
    meta['total_pas'] = np.asarray(pas.sum(1)).ravel()
    primary = pas[:, features.feature_index_1based.to_numpy() - 1].tocsc()
    primary.sort_indices()
    arrays[sid + '_pas_data'] = primary.data
    arrays[sid + '_pas_indices'] = primary.indices
    arrays[sid + '_pas_indptr'] = primary.indptr
    arrays[sid + '_pas_shape'] = np.array(primary.shape)
    arrays[sid + '_features_json'] = records(features)
    arrays[sid + '_meta_json'] = records(meta)
    audit.append(dict(study_id=sid, species=species, all_native_mc=len(meta),
                      old_selected_mc=len(selected), planned_route_mc=int(meta.in_planned_route.sum()),
                      planned_route_libraries=int(meta.loc[meta.in_planned_route, 'library_key'].nunique()),
                      unselected_planned_mc=int((meta.in_planned_route & ~meta.original_selected).sum()),
                      primary_pas=len(features), primary_sparse_nonzero=int(primary.nnz),
                      all_gex_joined=True, source_integer_counts=True))
    print(json.dumps(audit[-1]), file=sys.stderr, flush=True)
    del pas, primary

for path, expected in provenance.items():
    assert sha(path) == expected, ('input_changed_during_read', path)
arrays['mapping_json'] = np.array(json.dumps(mapping))
arrays['audit_json'] = np.array(json.dumps(audit))
arrays['provenance_json'] = np.array(json.dumps(provenance, sort_keys=True))
arrays['config_json'] = np.array(json.dumps(CONFIG, sort_keys=True))
arrays['panels_json'] = np.array(json.dumps({'P_shared': P_SHARED, 'P_legacy': P_LEGACY, 'ECM': ECM}, sort_keys=True))
np.savez_compressed(sys.stdout.buffer, **arrays)
print('READ_ONLY_INPUT_EXTRACTION_COMPLETE', file=sys.stderr, flush=True)
