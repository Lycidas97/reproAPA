#!/usr/bin/env python3
"""First migrated numerical adapter: aggregate fresh PAS counts on frozen metacells.

Runs inside run.py as the metacell node. Does not filter PASs or reuse old
PAS-column indexes, coverage flags, qualified genes, or supported-site matches.
"""
import json
import os
from pathlib import Path

import numpy as np
from scipy import sparse

from intake import artifact, digest, load_matrix, require, table


def aggregate_counts(x, cells, members, metacells):
    barcodes = {r['original_barcode']: i for i, r in enumerate(cells)}
    mc_index = {r['metacell_id']: i for i, r in enumerate(metacells)}
    indices = np.array([barcodes[r['barcode']] for r in members])
    groups = np.array([mc_index[r['metacell_id']] for r in members])
    require(len(set(indices)) == len(indices), 'Repeated selected cells')
    mapping = sparse.csr_matrix((np.ones(len(members), dtype=np.int64), (groups, indices)), shape=(len(metacells),len(cells)))
    result = (mapping @ x).tocsr()
    require(np.array_equal(result.sum(0).A1, x[indices].sum(0).A1), 'Selected-cell/PAS count conservation failed')
    require(np.array_equal(mapping.sum(1).A1, [int(r['n_cells']) for r in metacells]), 'Frozen membership sizes changed')
    return result


def main():
    manifest = Path(os.environ['PAS_INPUT_MANIFEST'])
    m = json.loads(manifest.read_text()); out = Path(os.environ['PAS_NODE_DIR'])
    artifacts, summary = [], []
    def record(path, role):
        artifacts.append({'path':str(path.relative_to(out)), 'sha256':digest(path), 'role':role})
    for lib in m['libraries']:
        folder = out/lib['library_id']
        require(folder.resolve().is_relative_to(out.resolve()), 'Invalid library folder')
        folder.mkdir()
        study = next(s for s in m['studies'] if s['study_id']==lib['study_id'])
        paths = {key: artifact(lib[key], manifest.parent) for key in ['counts','cells','membership','metacells']}
        features = table(artifact(study['features'],manifest.parent), ['atlas_pas_id'])
        cells = table(paths['cells'],['original_barcode'])
        members = table(paths['membership'],['barcode','metacell_id'])
        mc = table(paths['metacells'],['metacell_id','n_cells'])
        x = load_matrix(paths['counts'],lib['orientation'],len(cells),len(features))
        y = aggregate_counts(x,cells,members,mc)
        target = folder/'metacell_by_pas_counts.npz'; sparse.save_npz(target,y)
        require((sparse.load_npz(target)!=y).nnz==0, 'Matrix readback failed')
        record(target,'source_data')
        rows = folder/'matrix_rows.tsv'
        rows.write_text('metacell_id\n'+''.join(r['metacell_id']+'\n' for r in mc));record(rows,'source_data')
        cols = folder/'matrix_columns.tsv'
        cols.write_text('atlas_pas_id\n'+''.join(r['atlas_pas_id']+'\n' for r in features));record(cols,'source_data')
        summary.append(dict(library_id=lib['library_id'],selected_cells=len(members),metacells=y.shape[0],pas=y.shape[1],umi=int(y.sum()),zero_pas_metacells=int(np.sum(y.sum(1).A1==0))))
    audit=out/'aggregation-audit.json';audit.write_text(json.dumps(summary,indent=2)+'\n');record(audit,'audit')
    result=dict(node=os.environ['PAS_NODE_ID'],run_id=os.environ['PAS_RUN_ID'],input_fingerprint=os.environ['PAS_INPUT_FINGERPRINT'],status='complete',artifacts=artifacts,checks={'frozen_membership':True,'selected_count_conservation':True,'matrix_readback':True})
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
