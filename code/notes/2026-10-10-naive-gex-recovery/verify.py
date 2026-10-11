"""Independent sparse-column and existing RBP-function validation."""
from pathlib import Path
import ast
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    tables = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
    archive = np.load(HERE / 'raw/metacell_by_gene_counts.npz')
    genes = pd.read_csv(HERE / 'raw/genes.tsv.gz', sep='\t')
    gmeta = pd.read_csv(HERE / 'raw/metacells.tsv.gz', sep='\t')
    meta = pd.read_csv(tables / 'naive_state/frozen_selected_metadata.tsv.gz', sep='\t')
    meta = meta[meta.purity.ge(.8)].reset_index(drop=True)
    ix = pd.Index(gmeta.metacell_id).get_indexer(meta.metacell_id)
    factors = json.loads((HERE / 'protocol.json').read_text())['requested_factors']
    row_index = np.repeat(np.arange(len(gmeta)), np.diff(archive['indptr']))
    total = np.bincount(row_index, weights=archive['data'], minlength=len(gmeta)).astype(np.int64)
    assert np.array_equal(total, gmeta.total_gex)
    counts = np.full((len(gmeta), len(factors)), -1, np.int64)
    for j, factor in enumerate(factors):
        columns = genes.index[genes.gene_symbol.str.upper().eq(factor)].to_numpy()
        assert len(columns) == 1
        selected = archive['indices'] == columns[0]
        counts[:, j] = np.bincount(row_index[selected], weights=archive['data'][selected], minlength=len(gmeta)).astype(np.int64)
    legacy = ROOT / 'notes/2026-10-02-rbp-expression/analyze.py'
    frozen_audit_path = tables / 'rbp_expression/audit.json'
    frozen_audit = json.loads(frozen_audit_path.read_text())
    expected = {Path(p).name: h for p, h in frozen_audit['inputs'].items() if '/frozen_gex/GSE150311/' in p}
    assert len(expected) == 3
    assert all(digest(HERE / 'raw' / name) == value for name, value in expected.items())
    tree = ast.parse(legacy.read_text())
    functions = [f for f in tree.body if isinstance(f, ast.FunctionDef) and f.name == 'library_rows']
    assert len(functions) == 1
    ns = dict(np=np, pd=pd)
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(legacy), 'exec'), ns)
    baseline = ns['library_rows'](meta, counts[ix], total[ix], factors, 'GSE150311', 'human')
    current = pd.read_csv(HERE / 'factor_library_expression.tsv.gz', sep='\t')
    key = ['factor', 'library_key', 'state']
    b = baseline.set_index(key).sort_index(); c = current.set_index(key).sort_index()
    assert b.index.equals(c.index)
    numeric = ['n_metacells', 'n_cells', 'count', 'total_gex', 'pooled_CPM', 'mean_log2CPM1', 'mean_CPM', 'detected_metacells', 'detected_fraction']
    comparisons = []
    for column in numeric:
        error = (b[column] - c[column]).abs().max()
        assert np.allclose(b[column], c[column], rtol=1e-10, atol=1e-9), (column, error)
        comparisons.append(dict(comparison='independent_sparse_columns_original_rbp_function', field=column, n=len(b), max_absolute_error=float(error)))
    assert np.array_equal(b.qualified, c.qualified) and np.array_equal(b.mapped, c.mapped)
    existing_path = tables / 'rbp_expression/library_expression.tsv.gz'
    existing = pd.read_csv(existing_path, sep='\t')
    existing = existing[existing.study_id.eq('GSE150311') & existing.factor.isin(factors)].drop_duplicates(key)
    merged = current.merge(existing, on=key, suffixes=('_new', '_frozen'), validate='one_to_one')
    assert len(merged) == len(existing) and len(merged) > 0
    for column in numeric:
        error = (merged[column + '_new'] - merged[column + '_frozen']).abs().max()
        assert np.allclose(merged[column + '_new'], merged[column + '_frozen'], rtol=1e-9, atol=1e-8), (column, error)
        comparisons.append(dict(comparison='frozen_current_rbp_output', field=column, n=len(merged), max_absolute_error=float(error)))
    pd.DataFrame(comparisons).to_csv(HERE / 'verification.tsv', sep='\t', index=False)
    checks = dict(status='PASS', independent_sparse_column_reconstruction=True,
        original_library_rows_function_reproduced=True, legacy_function_sha256=digest(legacy),
        existing_output_sha256=digest(existing_path), target_factor_library_groups=len(b),
        frozen_prior_audit_sha256=digest(frozen_audit_path), raw_source_hashes_match_prior_frozen_audit=True,
        frozen_RBP_matching_groups=len(merged), frozen_RBP_matching_factors=sorted(merged.factor.unique()),
        checks=len(comparisons), tests_run=3, tests_skipped=0,
        all_result_hashes={str(p): digest(p) for p in sorted(HERE.glob('*.tsv*'))},
        code_sha256=digest(Path(__file__)))
    (HERE / 'validation.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps({k: checks[k] for k in ['status', 'target_factor_library_groups', 'frozen_RBP_matching_groups', 'frozen_RBP_matching_factors', 'checks']}, indent=2))


if __name__ == '__main__':
    main()
