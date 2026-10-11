"""Independent consistency checks on exported temporal summaries."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(name):
    return pd.read_csv(HERE / name, sep='\t')


def main():
    checks = {}
    manifest = json.loads((HERE / 'manifest.json').read_text())
    for path, expected in manifest['inputs'].items():
        assert sha(ROOT / path) == expected
    assert sha(HERE / 'analyze.py') == manifest['analysis_code_sha256']
    assert sha(HERE / 'protocol.json') == manifest['protocol_sha256']
    checks['all_input_code_protocol_hashes_unchanged'] = True
    keys = ['species', 'study_id', 'route', 'metric', 'subset']
    curves = read('curve_summary.tsv')
    libraries = read('curve_library_values.tsv')
    dates = read('curve_date_values.tsv')
    members = read('curve_fixed_gene_membership.tsv.gz')
    endpoints = read('curve_gene_endpoints.tsv.gz')
    summary_index = curves.set_index(keys).sort_index()
    library_index = libraries.set_index(keys).sort_index()
    member_index = members.set_index(keys).sort_index()
    for key, z in dates.groupby(keys):
        summary = summary_index.loc[key]
        selected = library_index.loc[key].reset_index()
        assert z.n_genes.eq(summary.n_genes).all()
        expected = selected.groupby(['time', 'sample_unit']).value.mean().groupby('time').mean()
        actual = z.set_index('time').value
        np.testing.assert_allclose(actual.sort_index(), expected.sort_index(), rtol=0, atol=2e-12)
        assert member_index.loc[key].gene_id.nunique() == summary.n_genes
        assert selected.library_key.nunique() == summary.n_libraries
    checks['date_means_independently_reconstructed_from_library_outputs'] = True
    for key, z in endpoints.groupby(keys):
        summary = summary_index.loc[key]
        assert summary.status == 'measured'
        assert len(z) == summary.n_genes
        assert abs(z.delta_DWUI.mean() - summary.last_minus_first_mean) < 2e-12
        np.testing.assert_allclose(z.last_value - z.first_value, z.delta_DWUI, rtol=0, atol=2e-12)
        if key[-1] == 'full_primary_tandem_UTR':
            assert abs(z.delta_retained_nt.mean() - summary.last_minus_first_mean_nt) < 2e-9
            assert abs(z.delta_retained_nt.median() - summary.last_minus_first_median_nt) < 2e-9
    assert not endpoints.route.isin(['GSE118258', '178325_0618', '178325_0809']).any()
    checks['endpoint_differences_and_UTR_length_summaries_verified'] = True
    checks['single_endpoint_routes_have_no_temporal_difference_rows'] = True
    components = read('PM_curve_component_summary.tsv')
    np.testing.assert_allclose(components.PM_last_minus_first, components.P_last_minus_first - components.M_last_minus_first, rtol=0, atol=2e-12)
    checks['same_PM_curve_contributor_temporal_identity_verified'] = True
    profiles = read('stage_gene_profiles.tsv.gz')
    np.testing.assert_allclose(profiles.R_S, profiles.R - profiles.S, rtol=0, atol=2e-12)
    np.testing.assert_allclose(profiles.E_R, profiles.E - profiles.R, rtol=0, atol=2e-12)
    np.testing.assert_allclose(profiles.P_E, profiles.P - profiles.E, rtol=0, atol=2e-12)
    np.testing.assert_allclose(profiles.R_S + profiles.E_R + profiles.P_E, profiles.P - profiles.S, rtol=0, atol=3e-12)
    checks['adjacent_state_step_identity_verified'] = True
    tests = read('date_level_tests.tsv')
    assert len(tests) == 2 and tests.route.eq('WOT').all() and tests.subset.eq('all_fixed').all()
    draws = np.load(HERE / 'date_slope_null_draws.npz')
    for row in tests.itertuples(index=False):
        sample = draws[row.metric]
        expected = (1 + int((np.abs(sample) >= abs(row.slope) - 1e-14).sum())) / (len(sample) + 1)
        assert abs(expected - row.p_value) < 1e-12
        z = dates[dates.route.eq('WOT') & dates.subset.eq('all_fixed') & dates.metric.eq(row.metric)]
        x = z.day.to_numpy()
        y = z.value.to_numpy()
        expected_slope = ((x - x.mean()) * (y - y.mean())).sum() / ((x - x.mean()) ** 2).sum()
        assert abs(expected_slope - row.slope) < 2e-12
    checks['exact_two_date_test_family_and_null_tails_verified'] = True
    assert 'Ran 5 tests' in (HERE / 'tests.log').read_text() and '\nOK\n' in (HERE / 'tests.log').read_text()
    checks['intent_tests_passed'] = 5
    checks['intent_tests_skipped'] = 0
    output = dict(status='passed', verification=checks,
                  output_sha256={p.name: sha(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name != 'validation.json'})
    (HERE / 'validation.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
