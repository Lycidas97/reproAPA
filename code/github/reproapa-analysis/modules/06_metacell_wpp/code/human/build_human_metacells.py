"""Mouse-matched library-local SEACells and exact full-gene count aggregation."""
import argparse
import hashlib
import importlib.metadata
import json
import math
import random
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import SEACells

P = Path(__file__).resolve().parents[1]
CONFIG = json.loads((P / 'metacell_config.json').read_text())
OUT = P / CONFIG['output']
FIELDS = ['original_state_id', 'coarse_state', 'experimental_time', 'time_kind',
          'condition', 'donor', 'pool', 'capture', 'series', 'stage']


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def raw_counts(path, cells):
    """Exact library/barcode identity; retain every native Gene Expression feature."""
    assert cells.cell_id.is_unique and cells.barcode.is_unique
    with h5py.File(path, 'r') as f:
        g = f['matrix']
        barcodes = pd.Index(g['barcodes'][:].astype(str))
        assert barcodes.is_unique
        rows = barcodes.get_indexer(cells.barcode)
        assert (rows >= 0).all(), 'Missing exact raw GEX barcode'
        keep = g['features/feature_type'][:].astype(str) == 'Gene Expression'
        features = pd.DataFrame({'gene_id': g['features/id'][:].astype(str)[keep],
                                 'gene_symbol': g['features/name'][:].astype(str)[keep]})
        assert features.gene_id.is_unique
        x = sparse.csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]),
                              shape=tuple(g['shape'][:])).T.tocsr()[rows][:, keep]
    assert np.isfinite(x.data).all() and (x.data >= 0).all()
    assert np.equal(x.data, np.floor(x.data)).all(), 'Noninteger raw counts'
    x = x.astype(np.int64)
    x.sum_duplicates(); x.eliminate_zeros()
    return x, features


def prepare():
    dest = OUT / 'prepared'
    dest.mkdir(parents=True, exist_ok=True)
    assert not (dest / 'contract.json').exists(), 'Prepared contract already frozen'
    files = ['results_v1/prepared/primary_cells.tsv.gz', 'results_v1/prepared/libraries.tsv',
             'results_v1/prepared/source_provenance.json', 'results_v1/model/hvg_features.tsv',
             'results_temporal_v1/layout/time_coordinates.npy',
             'results_temporal_v1/layout/time_coordinates.tsv.gz',
             'data/state-flow-v1/local_to_coarse_state.tsv', 'figure_style.json', 'metacell_config.json']
    cells = pd.read_csv(P / files[0], sep='\t', dtype=str).fillna('')
    libs = pd.read_csv(P / files[1], sep='\t', dtype=str).fillna('')
    coords = pd.read_csv(P / files[5], sep='\t', usecols=['cell_id'])
    xy = np.load(P / files[4])
    states = pd.read_csv(P / files[6], sep='\t').set_index('state_id', verify_integrity=True)
    assert cells.cell_id.is_unique and len(cells) == CONFIG['expected_cells']
    assert libs.library_key.is_unique and len(libs) == CONFIG['expected_libraries']
    assert cells.original_state_id.nunique() == len(states) == CONFIG['expected_states']
    assert np.array_equal(cells.cell_id, coords.cell_id) and xy.shape == (len(cells), 2)
    cells['UMAP1'], cells['UMAP2'] = xy.T
    cells['coarse_state'] = cells.original_state_id.map(states.coarse_state)
    cells['coarse_label'] = cells.original_state_id.map(states.coarse_label)
    assert cells.coarse_state.notna().all()
    assert cells.groupby('library_key').size().reindex(libs.library_key).eq(libs.primary_cells.astype(int).to_numpy()).all()
    features = pd.read_csv(P / files[3], sep='\t')
    assert features.gene_id.is_unique and len(features) == CONFIG['n_hvg']
    libs['requested_metacells'] = np.maximum(2, np.ceil(libs.primary_cells.astype(int) / CONFIG['target_cells'])).astype(int)
    cells.to_csv(dest / 'cells.tsv.gz', sep='\t', compression='gzip', index=False)
    libs.to_csv(dest / 'libraries.tsv', sep='\t', index=False)
    features.to_csv(dest / 'hvg_features.tsv', sep='\t', index=False)
    for key, group in cells.groupby('library_key', sort=True):
        d = OUT / 'libraries' / key
        d.mkdir(parents=True, exist_ok=True)
        group.to_csv(d / 'input_cells.tsv.gz', sep='\t', compression='gzip', index=False)
    dump(dest / 'contract.json', {'status': 'INPUTS_FROZEN', 'source_hashes': {f: sha256(P / f) for f in files},
        'config': CONFIG, 'requested_metacells': int(libs.requested_metacells.sum()),
        'source_method': 'mouse scripts/run_seacells_library.py, accepted target_50',
        'construction': 'frozen human HVG counts; within-library normalization/log1p/PCA/SEACells',
        'adaptations': 'human feature identities and metadata; no unused library UMAP computation',
        'time_or_states_used_for_construction': False, 'cross_library_membership': False})
    print('PREPARED', len(cells), len(libs), int(libs.requested_metacells.sum()), flush=True)


def run_library(key):
    start = time.monotonic()
    dest = OUT / 'libraries' / key
    assert not (dest / 'complete.ok').exists(), 'Refusing to overwrite completed membership'
    record = pd.read_csv(OUT / 'prepared/libraries.tsv', sep='\t').set_index('library_key').loc[key]
    cells = pd.read_csv(dest / 'input_cells.tsv.gz', sep='\t', compression='gzip', dtype=str).fillna('')
    assert cells.library_key.eq(key).all() and len(cells) == int(record.primary_cells)
    cells[['UMAP1', 'UMAP2']] = cells[['UMAP1', 'UMAP2']].astype(float)
    source = Path(record.input_path)
    st = source.stat()
    raw_prov = json.loads((P / 'results_v1/prepared/source_provenance.json').read_text())
    frozen = next(r for r in raw_prov['raw_matrix_stats'] if r['path'] == str(source))
    assert st.st_size == frozen['bytes'] and st.st_mtime_ns == frozen['mtime_ns']
    source_hash = sha256(source)
    x, features = raw_counts(source, cells)
    hvg = pd.read_csv(OUT / 'prepared/hvg_features.tsv', sep='\t')
    cols = pd.Index(features.gene_id).get_indexer(hvg.gene_id)
    assert (cols >= 0).all()
    work = ad.AnnData(X=x[:, cols].astype(np.float32),
                     obs=pd.DataFrame(index=pd.Index(cells.cell_id, name='cell_id')))
    assert (np.asarray(work.X.sum(1)).ravel() > 0).all()
    random.seed(CONFIG['seed']); np.random.seed(CONFIG['seed'])
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    n_pcs = min(CONFIG['n_pcs'], work.n_obs - 1, work.n_vars - 1)
    sc.pp.pca(work, n_comps=n_pcs, random_state=CONFIG['seed'])
    model = SEACells.core.SEACells(work, build_kernel_on='X_pca',
        n_SEACells=int(record.requested_metacells), use_gpu=False, verbose=True,
        n_waypoint_eigs=min(CONFIG['n_waypoint_eigs'], n_pcs - 1), n_neighbors=CONFIG['n_neighbors'],
        convergence_epsilon=CONFIG['convergence_epsilon'], use_sparse=True)
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=CONFIG['min_iter'], max_iter=CONFIG['max_iter'])
    rss = np.asarray(model.RSS_iters, dtype=float)
    pd.DataFrame({'iteration': np.arange(len(rss)), 'rss': rss}).to_csv(dest / 'convergence.tsv', sep='\t', index=False)
    assert np.isfinite(rss).all() and len(rss) >= CONFIG['min_iter'] + 1
    final_converged = bool(abs(rss[-1] - rss[-2]) < model.convergence_threshold)
    assert final_converged, 'Final RSS change exceeds the convergence criterion'
    assignments = model.get_hard_assignments().reindex(work.obs_names)['SEACell']
    assert assignments.notna().all()
    labels = sorted(assignments.astype(str).unique())
    ids = {s: f'{key}.t{CONFIG["target_cells"]:03d}.MC{i:04d}' for i, s in enumerate(labels, 1)}
    cells['metacell_id'] = assignments.astype(str).map(ids).to_numpy()
    cells.to_csv(dest / 'cell_to_metacell.tsv.gz', sep='\t', compression='gzip', index=False)
    nodes, compositions = [], []
    for mc, group in cells.groupby('metacell_id', sort=True):
        xy = group[['UMAP1', 'UMAP2']].to_numpy(float)
        center = xy.mean(0)
        row = {'metacell_id': mc, 'study_id': record.study_id, 'library_key': key,
               'library_id': record.library_id, 'sample_accession': record.sample_accession,
               'n_cells': len(group), 'centroid_UMAP1': center[0], 'centroid_UMAP2': center[1],
               'umap_radius_q90': float(np.quantile(np.linalg.norm(xy - center, axis=1), .9))}
        for field in FIELDS:
            counts = group[field].replace('', 'missing').value_counts()
            maximum = int(counts.max())
            ties = sorted(counts[counts.eq(maximum)].index)
            row[field] = ties[0]
            row[field + '_purity'] = maximum / len(group)
            row[field + '_tie'] = len(ties) > 1
            for label, count in counts.items():
                compositions.append({'metacell_id': mc, 'field': field, 'label': label,
                                     'n_cells': int(count), 'fraction': int(count) / len(group)})
        eligible = group[group.original_state_id.eq(row['original_state_id'])].sort_values('cell_id')
        distances = np.linalg.norm(eligible[['UMAP1', 'UMAP2']].to_numpy(float) - center, axis=1)
        rep = eligible.iloc[int(np.argmin(distances))]
        row.update(representative_cell_id=rep.cell_id, UMAP1=rep.UMAP1, UMAP2=rep.UMAP2,
                   coarse_purity_ge08=row['coarse_state_purity'] >= .8 and not row['coarse_state_tie'],
                   coarse_purity_ge09=row['coarse_state_purity'] >= .9 and not row['coarse_state_tie'])
        nodes.append(row)
    metadata = pd.DataFrame(nodes).set_index('metacell_id')
    composition = pd.DataFrame(compositions)
    assert np.allclose(composition.groupby(['metacell_id', 'field']).fraction.sum(), 1)
    codes = metadata.index.get_indexer(cells.metacell_id)
    g = sparse.csr_matrix((np.ones(len(cells), dtype=np.int64), (codes, np.arange(len(cells)))),
                         shape=(len(metadata), len(cells)))
    y = (g @ x).tocsr()
    assert np.array_equal(np.asarray(x.sum(0)), np.asarray(y.sum(0)))
    metadata['gex_umi'] = np.asarray(y.sum(1)).ravel()
    metadata['detected_genes'] = np.asarray((y > 0).sum(1)).ravel()
    data = ad.AnnData(X=y, obs=metadata.copy(), var=features.set_index('gene_id'))
    data.uns['count_scale'] = 'sum of raw nonnegative integer UMI counts over exact member cells'
    data.write_h5ad(dest / 'all_gene_counts.h5ad', compression='lzf')
    metadata.to_csv(dest / 'metacell_metadata.tsv.gz', sep='\t', compression='gzip')
    composition.to_csv(dest / 'state_composition.tsv.gz', sep='\t', compression='gzip', index=False)
    reread = ad.read_h5ad(dest / 'all_gene_counts.h5ad')
    assert (reread.X != y).nnz == 0 and reread.obs_names.equals(metadata.index)
    checked = []
    for mc in sorted({metadata.index[0], metadata.index[-1], metadata.n_cells.idxmin(), metadata.n_cells.idxmax()}):
        expected = np.asarray(x[cells.metacell_id.eq(mc).to_numpy()].sum(0)).ravel()
        assert np.array_equal(reread.X.getrow(metadata.index.get_loc(mc)).toarray().ravel(), expected)
        checked.append(mc)
    assert (source.stat().st_size, source.stat().st_mtime_ns) == (st.st_size, st.st_mtime_ns)
    summary = dict(status='GEX_MEMBERSHIP_COUNTS_PASS', library_key=key, study_id=record.study_id,
        cells=len(cells), requested_metacells=int(record.requested_metacells), metacells=len(metadata),
        genes=len(features), gex_umi=int(y.sum()), zero_gex_metacells=int(metadata.gex_umi.eq(0).sum()),
        cells_min=int(metadata.n_cells.min()), cells_median=float(metadata.n_cells.median()), cells_max=int(metadata.n_cells.max()),
        iterations=len(rss)-1, final_rss_delta=float(abs(rss[-1]-rss[-2])),
        convergence_threshold=float(model.convergence_threshold), converged=final_converged,
        exact_cell_assignment=True, per_gene_counts_conserved=True, direct_source_readback_metacells=checked,
        source_matrix=str(source), source_sha256=source_hash, source_unchanged=True,
        construction_script_sha256=sha256(__file__), membership_sha256=sha256(dest / 'cell_to_metacell.tsv.gz'),
        h5ad_sha256=sha256(dest / 'all_gene_counts.h5ad'), elapsed_seconds=round(time.monotonic()-start, 2),
        versions={pkg: importlib.metadata.version(pkg) for pkg in ['SEACells', 'scanpy', 'anndata', 'numpy', 'scipy', 'palantir']})
    dump(dest / 'summary.json', summary)
    (dest / 'complete.ok').write_text('GEX_MEMBERSHIP_COUNTS_PASS\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--library')
    args = parser.parse_args()
    prepare() if args.prepare else run_library(args.library)
