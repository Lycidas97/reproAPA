"""Apply the predeclared BH families across all agent outputs."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCES = [
    ('background/tests.tsv', 'p_value'),
    ('cross-species/results/raw-pvalues.tsv', 'raw_p'),
    ('bridges/raw-tests.tsv', 'raw_p'),
    ('cycle/cycle_score_tests.tsv', 'p'),
    ('cycle/polyapipe/dual_alignment_tests.tsv', 'p'),
    ('cycle/apadetect_p50/dual_alignment_tests.tsv', 'p'),
    ('cycle/polyapipe/pm_route_mean_tests.tsv', 'p'),
    ('cycle/apadetect_p50/pm_route_mean_tests.tsv', 'p'),
]


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda: f.read(2**20), b''):
            h.update(x)
    return h.hexdigest()


def bh(p):
    order = np.argsort(p)
    q = np.empty(len(p))
    q[order] = np.minimum(1., np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return q


def main():
    frames, provenance = [], {}
    for relative, column in SOURCES:
        path = HERE / relative
        f = pd.read_csv(path, sep='\t')
        f['source_table'] = relative
        f['source_row'] = np.arange(1, len(f) + 1)
        f['raw_p'] = pd.to_numeric(f[column])
        valid = f.raw_p.notna()
        assert np.isfinite(f.loc[valid, 'raw_p']).all()
        assert f.loc[valid, 'raw_p'].between(0., 1.).all()
        assert f.family.notna().all()
        frames.append(f)
        provenance[relative] = digest(path)
    all_tests = pd.concat(frames, ignore_index=True)
    all_tests['q_BH_family'] = np.nan
    all_tests['family_valid_tests'] = 0
    all_tests['status_for_correction'] = np.where(all_tests.raw_p.notna(), 'estimable', 'unestimable_explicitly_reported')
    families = []
    for family, part in all_tests.groupby('family'):
        ix = part.index[part.raw_p.notna()]
        if len(ix):
            all_tests.loc[ix, 'q_BH_family'] = bh(all_tests.loc[ix, 'raw_p'].to_numpy())
        all_tests.loc[part.index, 'family_valid_tests'] = len(ix)
        families.append(dict(family=family, reported_rows=len(part), estimable_tests=len(ix),
                             unestimable_rows=len(part)-len(ix),
                             rejected_BH_05=int(all_tests.loc[ix, 'q_BH_family'].lt(.05).sum())))
    assert all_tests[(all_tests.family.eq('A3-cross-species'))].shape[0] == 4
    assert all_tests[(all_tests.family.eq('A5-interval-bridge'))].shape[0] == 4
    assert all_tests[(all_tests.family.eq('A6-shape-endpoint'))].shape[0] == 4
    scores = pd.read_csv(HERE/'cycle/cycle_score_library_states.tsv', sep='\t')
    assert set(scores.study_id) == {'GSE103221','GSE115943','GSE171127',
                                   'GSE118258','GSE150311','GSE178325','GSE221739'}
    all_tests['rejected_BH_05'] = all_tests.q_BH_family.lt(.05)
    all_tests.to_csv(HERE/'corrected-tests.tsv', sep='\t', index=False, float_format='%.17g')
    pd.DataFrame(families).to_csv(HERE/'test-families.tsv', sep='\t', index=False)
    assert all(digest(HERE/p) == h for p, h in provenance.items())
    (HERE/'correction-receipt.json').write_text(json.dumps(dict(
        status='complete', procedure='BH within predeclared task families, both PAS included',
        unestimable='NA retained with status; excluded from valid-test family size',
        input_sha256=provenance, code_sha256=digest(Path(__file__)),
        reported_rows=len(all_tests), estimable_tests=int(all_tests.raw_p.notna().sum()),
        all_seven_study_cycle_scores=True), indent=2)+'\n')
    print(pd.DataFrame(families).to_string(index=False))


if __name__ == '__main__':
    main()
