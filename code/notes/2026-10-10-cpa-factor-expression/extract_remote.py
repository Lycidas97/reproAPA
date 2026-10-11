"""Read five pinned GEX matrices only; emit the small requested-factor subset."""
from pathlib import Path
import gzip
import hashlib
import io
import json
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    provenance[str(path)] = sha(path)
    return pd.read_csv(path, sep='\t')


arrays, provenance, mapping = {}, {}, []
native = Path(CONFIG['native_run']) / 'studies'
gex = Path(CONFIG['frozen_gex'])
annotations = {}
for species in ['human', 'mouse']:
    spec = CONFIG['resources'][species + '_annotation_gtf']
    path = Path(spec['path']); provenance[str(path)] = sha(path)
    assert provenance[str(path)] == spec['sha256']
    names = {}
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2] != 'gene':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            if 'gene_name' in attrs:
                gene = attrs['gene_id'].split('.')[0]
                assert gene not in names or names[gene] == attrs['gene_name']
                names[gene] = attrs['gene_name']
    annotations[species] = names

for study, spec in PROTOCOL['studies'].items():
    if study == 'GSE150311':
        continue  # Independently collected by the all-gene naive-expression task.
    species, routes = spec['species'], spec['routes']
    folder = 'WOT_GSE115943' if study == 'GSE115943' else study
    root = native / folder
    selected = read(root / 'selected_metacells.tsv.gz')
    full = read(root / 'metacells.tsv.gz')
    columns = [c for c in ['condition', 'sample_accession', 'time_kind', 'original_state_id'] if c in full]
    selected = selected.merge(full[['metacell_id'] + columns], on='metacell_id', validate='one_to_one')
    meta = selected[selected.route.isin(routes) & selected.node.isin(PROTOCOL['state_filter']['nodes']) & selected.purity.ge(.8)].reset_index(drop=True)
    gm = read(gex / folder / 'metacells.tsv.gz')
    genes = read(gex / folder / 'genes.tsv.gz')
    indices = pd.Index(gm.metacell_id).get_indexer(meta.metacell_id)
    assert (indices >= 0).all() and pd.Index(gm.metacell_id).is_unique
    path = gex / folder / 'metacell_by_gene_counts.npz'
    provenance[str(path)] = sha(path)
    x = sparse.load_npz(path).tocsr()[indices]
    assert (x.data >= 0).all() and (x.data == np.floor(x.data)).all()
    totals = gm.total_gex.to_numpy()[indices].astype(np.int64)
    assert (totals > 0).all() and np.array_equal(np.asarray(x.sum(1)).ravel(), totals)
    symbols = genes.gene_symbol.fillna('').astype(str).copy()
    ids = genes.gene_id.astype(str).str.split('.').str[0]
    placeholder = symbols.eq('') | symbols.str.split('.').str[0].eq(ids)
    replacement = ids.map(annotations[species])
    symbols.loc[placeholder & replacement.notna()] = replacement[placeholder & replacement.notna()]
    factors = PROTOCOL['primary_factors'] + PROTOCOL['secondary_shared_factors']
    values = np.full((len(meta), len(factors)), -1, dtype=np.int64)
    for j, factor in enumerate(factors):
        positions = np.flatnonzero(symbols.str.upper().eq(factor))
        status = 'unique_exact_casefold_symbol' if len(positions) == 1 else 'missing_symbol' if not len(positions) else 'ambiguous_symbol'
        mapping.append(dict(species=species, study_id=study, factor=factor, gene_ids=';'.join(ids.iloc[positions]), symbols=';'.join(symbols.iloc[positions]), status=status))
        if len(positions) == 1:
            values[:, j] = x[:, positions[0]].toarray().ravel().astype(np.int64)
    design = read(root / 'sampling_design.tsv')
    arrays[study + '_meta_json'] = np.array(meta.to_json(orient='records'))
    arrays[study + '_design_json'] = np.array(design[design.route.isin(routes)].to_json(orient='records'))
    arrays[study + '_factors'] = np.array(factors)
    arrays[study + '_counts'] = values
    arrays[study + '_total_gex'] = totals
    print(f'{study}: {len(meta)} selected metacells, {meta.library_key.nunique()} libraries, {len(factors)} factors', file=sys.stderr, flush=True)
arrays['provenance_json'] = np.array(json.dumps(provenance, sort_keys=True))
arrays['mapping_json'] = np.array(json.dumps(mapping))
buffer = io.BytesIO()
np.savez_compressed(buffer, **arrays)
sys.stdout.buffer.write(buffer.getvalue())
