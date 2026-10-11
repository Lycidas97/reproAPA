"""Independent cycle arithmetic used by the integration audit."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPECTED_PROTOCOL = '85677ecc53b5bca9f6dbc5d59dc1789274e3c130733c23a57b1f62bcb488051a'


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    assert digest(HERE / 'PROTOCOL.zh.md') == EXPECTED_PROTOCOL
    source = ROOT / 'notes/2026-10-10-proliferation-apa/inputs.npz'
    archive = np.load(source, allow_pickle=False)
    results, errors = [], []
    for key in archive.files:
        if not key.endswith('_meta_json'):
            continue
        study = key[:-10]
        meta = pd.DataFrame(json.loads(str(archive[key])))
        counts = archive[study + '_marker_counts']
        total = archive[study + '_total_gex']
        controls = archive[study + '_control_counts']
        assert counts.shape == (len(meta), 18)
        assert np.issubdtype(counts.dtype, np.integer) and (counts >= 0).all()
        assert (total > 0).all() and controls.shape == (len(meta), 360)
        values = np.log1p(counts * 10000. / total[:, None])
        background = np.log1p(controls * 10000. / total[:, None]).reshape(len(meta), 18, 20).mean(2)
        cycle = (values[:, :10].mean(1) + values[:, 10:].mean(1)) / 2
        adjusted = ((values - background)[:, :10].mean(1) + (values - background)[:, 10:].mean(1)) / 2
        for name, value in [('Cycle', cycle), ('Adjusted_Cycle', adjusted)]:
            error = float(np.abs(meta[name] - value).max())
            assert error < 1e-12
            errors.append(dict(study_id=study, score=name, n_metacells=len(meta), max_abs_error=error))
            meta[name] = value
        fields = ['metacell_id', 'study_id', 'route', 'library_key', 'node', 'purity', 'n_cells', 'Cycle', 'Adjusted_Cycle']
        results.append(meta[fields])
    scores = pd.concat(results, ignore_index=True)
    assert scores.metacell_id.is_unique and scores.purity.ge(.8).all()
    production = pd.read_csv(ROOT / 'notes/2026-10-10-proliferation-apa/metacell_cycle_scores.tsv.gz', sep='\t')
    joined = scores.merge(production, on='metacell_id', suffixes=('_raw', '_frozen'), validate='one_to_one', how='outer', indicator=True)
    assert joined._merge.eq('both').all()
    for name in ['study_id', 'route', 'library_key', 'node', 'n_cells']:
        assert joined[name + '_raw'].eq(joined[name + '_frozen']).all()
    for name in ['Cycle', 'Adjusted_Cycle']:
        assert np.allclose(joined[name + '_raw'], joined[name + '_frozen'], rtol=0, atol=1e-12)
    scores.to_csv(HERE / 'independent_mc_scores.tsv.gz', sep='\t', index=False, float_format='%.17g', compression={'method': 'gzip', 'mtime': 0})
    result = dict(status='PASS', protocol_sha256=EXPECTED_PROTOCOL, source_sha256=digest(source), n_metacells=len(scores),
                  n_studies=len(results), checks=errors, score_denominator='all_gene_GEX_CP10K', tests_skipped=0)
    (HERE / 'input_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'checks'}, indent=2))


if __name__ == '__main__':
    main()
