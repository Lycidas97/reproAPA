"""Verify correspondence and vector invariants without response-set labels."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'


def main():
    paths = [TABLES / f for f in [
        'cis/reciprocal_PAS_candidates.tsv.gz',
        'cis/reciprocal_one_to_one.tsv.gz',
        'cis/annotation_catalog.tsv.gz',
        'cis/site_sequence_features.tsv.gz',
        'endpoints/route_effects.tsv.gz',
    ]]
    candidate, actual, catalog, sequence = [pd.read_csv(p, sep='\t') for p in paths[:4]]
    order = candidate.sort_values(
        ['reciprocal25', 'reciprocal50', 'forward_distance', 'backward_distance', 'human_id', 'mouse_id'],
        ascending=[False, False, True, True, True, True],
    )
    used_h, used_m, chosen = set(), set(), []
    for row in order.itertuples():
        if row.human_id in used_h or row.mouse_id in used_m:
            continue
        chosen.append(row.Index)
        used_h.add(row.human_id)
        used_m.add(row.mouse_id)
    expected = candidate.loc[chosen].reset_index(drop=True)
    pd.testing.assert_frame_equal(expected, actual, check_exact=True)
    assert actual.human_id.is_unique and actual.mouse_id.is_unique
    assert candidate.reciprocal25.eq(candidate[['forward_distance', 'backward_distance']].max(axis=1).le(25)).all()
    assert candidate.reciprocal50.eq(candidate[['forward_distance', 'backward_distance']].max(axis=1).le(50)).all()
    ids = catalog.set_index('id')
    assert catalog.id.is_unique and sequence.id.is_unique
    assert set(catalog.id) == set(sequence.id)
    for sp in ['human', 'mouse']:
        assert ids.loc[actual[sp + '_id'], 'gene_id'].to_numpy().tolist() == actual[sp + '_gene_id'].tolist()
        assert ids.loc[actual[sp + '_id'], 'species'].eq(sp).all()
    assert sequence.phastcons.dropna().between(0, 1).all()
    effects = pd.read_csv(paths[4], sep='\t', usecols=['PAS_ids', 'target_usage', 'reference_usage', 'delta_usage'])
    max_error = 0.
    for row in effects.itertuples():
        pas = json.loads(row.PAS_ids)
        target, reference, delta = [np.asarray(json.loads(getattr(row, name)), float) for name in ['target_usage', 'reference_usage', 'delta_usage']]
        assert len(pas) == len(target) == len(reference) == len(delta)
        assert len(set(pas)) == len(pas)
        assert np.isfinite(target).all() and np.isfinite(reference).all()
        assert target.min() >= -1e-12 and reference.min() >= -1e-12
        assert abs(target.sum() - 1) < 1e-12 and abs(reference.sum() - 1) < 1e-12
        error = float(np.max(np.abs(target - reference - delta)))
        assert error < 1e-12
        max_error = max(max_error, error)
    hashes = {}
    for p in paths:
        with p.open('rb') as stream:
            hashes[str(p.relative_to(ROOT))] = hashlib.file_digest(stream, 'sha256').hexdigest()
    report = dict(
        status='PASS', response_sets_read=False,
        candidate_rows=len(candidate), assigned_rows=len(actual),
        assignment='frozen_greedy_reproduced_exactly',
        catalog_rows=len(catalog), endpoint_vector_rows=len(effects),
        max_delta_vector_error=max_error, input_sha256=hashes,
        limitation='Conservation features share one mm10 track; this verifies frozen tables, not current cross-species event support.',
    )
    (HERE / 'frozen_resource_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'input_sha256'}), flush=True)


if __name__ == '__main__':
    main()
