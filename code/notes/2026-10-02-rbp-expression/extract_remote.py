"""Read original GEX matrices; CONFIG is supplied over stdin by collect.py."""
from pathlib import Path
import hashlib
import io
import json
import sys
import anndata as ad
import numpy as np
import pandas as pd

BASE = Path('/root/nfsdata/apadetect/analysis_projects')
M = BASE / 'mouse_time_aware_integration_pilot_v1'
H = BASE / 'human_time_aware_integration_v1'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


arrays, provenance, mappings = {}, {}, []
human = None
for study, spec in CONFIG['studies'].items():
    species, folder, route = spec
    path = M / 'results_apa_majority_shared_v1/profiles' / folder / 'selected_metacells.tsv.gz'
    provenance[str(path)] = sha(path)
    assert provenance[str(path)] == CONFIG['metadata_sha256'][study], 'Local/remote state inventory differs'
    meta = pd.read_csv(path, sep='\t')
    meta = meta[meta.route.eq(route) & meta.node.isin(list('SREPM')) & meta.purity.ge(.8)].reset_index(drop=True)
    genes = CONFIG['factors'][species]
    counts = np.full((len(meta), len(genes)), -1, dtype=np.int64)
    total = np.zeros(len(meta), dtype=np.int64)
    if species == 'human' and human is None:
        path = H / 'results_metacell_v1/release/human_metacell_all_gene_counts.h5ad'
        human = ad.read_h5ad(path, backed='r')
        provenance[str(path)] = sha(path)
    for library, ix in meta.groupby('library_key', sort=False).groups.items():
        if species == 'mouse':
            path = M / 'results_mouse_metacell_gex_completion_v1' / library / 'all_gene_counts.h5ad'
            a = ad.read_h5ad(path, backed='r')
            provenance[str(path)] = sha(path)
        else:
            a = human
        gene_ids = pd.Index(a.var.gene_ids.astype(str) if 'gene_ids' in a.var else a.var_names.astype(str)).str.replace(r'\.[0-9]+$', '', regex=True)
        col = next((c for c in ['gene_symbols', 'gene_symbol', 'symbol', 'gene_name'] if c in a.var), None)
        symbols = a.var[col].astype(str).to_numpy() if col else a.var_names.astype(str).to_numpy()
        symbols_upper = np.array([s.upper() for s in symbols])
        ri = a.obs_names.get_indexer(meta.loc[ix, 'metacell_id'])
        assert (ri >= 0).all()
        assert np.array_equal(a.obs.iloc[ri].n_cells, meta.loc[ix, 'n_cells'])
        x = a.X[ri].tocsr()
        assert (x.data >= 0).all() and (x.data == np.floor(x.data)).all()
        total[ix] = np.asarray(x.sum(1)).ravel().astype(np.int64)
        for j, factor in enumerate(genes):
            pos = np.flatnonzero(symbols_upper == factor.upper())
            ids = sorted(set(gene_ids[pos]))
            status = 'unique_gene_id' if len(ids) == 1 else 'absent' if not ids else 'ambiguous_gene_ids'
            mappings.append(dict(species=species, study_id=study, library_key=library, factor=factor,
                status=status, gene_ids=';'.join(ids), symbols=';'.join(sorted(set(symbols[pos]))), n_columns=len(pos)))
            if status == 'unique_gene_id':
                counts[ix, j] = np.asarray(x[:, pos].sum(1)).ravel().astype(np.int64)
        if species == 'mouse':
            a.file.close()
    assert (total > 0).all()
    arrays[study + '_meta_json'] = np.array(meta.to_json(orient='records'))
    arrays[study + '_factors'] = np.array(genes)
    arrays[study + '_counts'] = counts
    arrays[study + '_total_gex'] = total
    print(f'{study}: {len(meta)} metacells, {meta.library_key.nunique()} libraries, {len(genes)} RBP factors', file=sys.stderr, flush=True)
if human is not None:
    human.file.close()
arrays['provenance_json'] = np.array(json.dumps(provenance))
arrays['mapping_json'] = np.array(json.dumps(mappings))
buffer = io.BytesIO()
np.savez_compressed(buffer, **arrays)
sys.stdout.buffer.write(buffer.getvalue())
