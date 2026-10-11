"""Independently check reported summaries and matched-reference aggregates."""
from pathlib import Path
from collections import Counter
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / '2026-10-10-marker-reference-sensitivity'
HASHES = {}


def read(name, **kwargs):
    path = SOURCE / name
    HASHES[name] = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kwargs)


def equal(actual, expected, label):
    if pd.isna(expected):
        assert pd.isna(actual), (label, actual, expected)
    else:
        assert np.isclose(actual, expected, rtol=0, atol=1e-12), (label, actual, expected)


def check_summaries(effects, summaries, keys, gradient=False):
    indexed = {key if isinstance(key, tuple) else (key,): part
               for key, part in effects.groupby(keys, dropna=False, sort=False)}
    assert not summaries.duplicated(keys + ['threshold']).any()
    for row in summaries.to_dict('records'):
        key = tuple(row[k] for k in keys)
        q = indexed.get(key, effects.iloc[:0])
        equal(row['N'], len(q), ('N', key))
        assert row['status'] == ('descriptive' if len(q) else 'not_estimable'), key
        threshold = row['threshold']
        if gradient:
            for layer in ['low', 'middle', 'high', 'whole']:
                delta = q['delta_' + layer].to_numpy(float)
                equal(row['mean_' + layer], delta.mean() if len(q) else np.nan, key)
                equal(row['n_proximal_' + layer], np.sum(delta <= -threshold), key)
                equal(row['n_distal_' + layer], np.sum(delta >= threshold), key)
                equal(row['fraction_proximal_' + layer], np.mean(delta <= -threshold) if len(q) else np.nan, key)
            change = q.difference_high_minus_low.to_numpy(float)
            equal(row['mean_difference_high_minus_low'], change.mean() if len(q) else np.nan, key)
            equal(row['fraction_proximal_high_minus_low'], row['fraction_proximal_high'] - row['fraction_proximal_low'], key)
            for label in ['increasing', 'decreasing', 'nonmonotonic']:
                equal(row['n_' + label], q.monotonicity.eq(label).sum(), key)
        else:
            delta = q.delta_DWUI.to_numpy(float)
            prox = delta <= -threshold
            distal = delta >= threshold
            equal(row['n_proximal'], prox.sum(), key)
            equal(row['n_distal'], distal.sum(), key)
            equal(row['n_within_threshold'], (~prox & ~distal).sum(), key)
            equal(row['fraction_proximal'], prox.mean() if len(q) else np.nan, key)
            equal(row['fraction_distal'], distal.mean() if len(q) else np.nan, key)
            equal(row['mean_delta'], delta.mean() if len(q) else np.nan, key)
            equal(row['median_delta'], np.median(delta) if len(q) else np.nan, key)
            equal(row['mean_proximal_conditional'], delta[prox].mean() if prox.any() else np.nan, key)
            libraries = set()
            if 'contributor_libraries' in q:
                for item in q.contributor_libraries:
                    libraries.update(json.loads(item))
            else:
                libraries.update(q.library_key)
            equal(row['n_contributor_libraries'], len(libraries), key)
    assert set(summaries.route) == {'CIP', 'WOT', 'GSE118258', 'GSE150311', '178325_0618', '178325_0809', 'GSE221739'}
    return len(summaries)


def check_layout(summaries, planned, backgrounds, library=False):
    scopes = planned[['study_id', 'route'] + (['library_key'] if library else [])].drop_duplicates()
    settings = [('supported50', v) for v in ['core_third_strict', 'core_quartile_strict',
                 'core_third_relaxed', 'native_third_strict']] + [('expanded_primary', 'core_third_strict')]
    columns = ['study_id', 'route'] + (['library_key'] if library else []) + ['setting', 'version', 'background', 'contrast', 'threshold']
    actual = set(summaries[columns].itertuples(index=False, name=None))
    required = {scope + (setting, version, background, contrast, threshold)
                for scope in scopes.itertuples(index=False, name=None)
                for setting, version in settings for background, contrasts in backgrounds.items()
                for contrast in contrasts for threshold in [.02, .05]}
    assert required.issubset(actual), ('missing_summary_scopes', sorted(required - actual)[:10])


def check_aggregates(pairs, routes):
    keys = ['study_id', 'route', 'version', 'setting', 'gene_id']
    assert not pairs.duplicated(keys + ['contrast', 'library_key']).any()
    assert pairs.time_a.eq(pairs.time_b).all() and pairs.sample_unit_a.eq(pairs.sample_unit_b).all()
    metadata_columns = ['study_id', 'route', 'library_key', 'time_a', 'sample_unit_a']
    metadata = pairs[metadata_columns].drop_duplicates()
    assert not metadata.duplicated(metadata_columns[:3]).any()
    sampling = {(r.study_id, r.route, r.library_key): (r.time_a, r.sample_unit_a)
                for r in metadata.itertuples(index=False)}
    indexed = {}
    for r in pairs.itertuples(index=False):
        key = tuple(getattr(r, k) for k in keys)
        indexed.setdefault(key, {}).setdefault(r.contrast, {})[r.library_key] = r.delta_DWUI
    intersections = {'common_four': {'D00', 'D01', 'D10', 'D11'},
                     'common_M_bridge': {'D00', 'D01'}, 'common_P_bridge': {'D00', 'D10'}}
    max_error = 0.
    weight_cache = {}
    expected_cache = {}
    trace_cache = {}
    for row in routes.itertuples(index=False):
        key = tuple(getattr(row, k) for k in keys)
        allpairs = indexed[key]
        background = row.background
        if background == 'full':
            expected = set(allpairs[row.contrast])
        else:
            token = (key, background)
            if token not in expected_cache:
                expected_cache[token] = set.intersection(*[set(allpairs.get(c, {})) for c in intersections[background]])
            expected = expected_cache[token]
        contributors = json.loads(row.contributor_libraries)
        assert len(contributors) == len(set(contributors)) and set(contributors) == expected, (key, background)
        # Same-library comparisons reduce the frozen hierarchy to equal dates,
        # equal biological units within dates, then equal captures per unit.
        signature = (row.study_id, row.route, tuple(contributors))
        if signature not in weight_cache:
            units = [sampling[(row.study_id, row.route, lib)] for lib in contributors]
            dates = {t for t, u in units}
            per_date = {t: len({u for tt, u in units if tt == t}) for t in dates}
            captures = Counter(units)
            weight_cache[signature] = np.array([1 / (len(dates) * per_date[t] * captures[(t, u)]) for t, u in units])
        weights = weight_cache[signature]
        stored = np.array(json.loads(row.weights))
        assert np.allclose(weights, stored, rtol=0, atol=1e-12), (key, background, 'weights')
        estimate = float(weights @ np.array([allpairs[row.contrast][lib] for lib in contributors]))
        equal(row.delta_DWUI, estimate, (key, background, row.contrast))
        equal(row.n_libraries, len(contributors), key)
        max_error = max(max_error, abs(estimate - row.delta_DWUI))
        if background in intersections:
            token = (key, background)
            trace = trace_cache.setdefault(token, [set(), row.contributor_libraries, row.weights])
            assert row.contrast not in trace[0] and row.contributor_libraries == trace[1] and row.weights == trace[2], token
            trace[0].add(row.contrast)
    for (key, background), trace in trace_cache.items():
        assert trace[0] == intersections[background], (key, background)
    return dict(rows=len(routes), max_abs_error=max_error)


def main():
    pairkeys = ['study_id', 'route', 'library_key', 'setting', 'version', 'contrast', 'gene_id']
    pairs = read('pair_effects.tsv.gz', usecols=['study_id', 'route', 'version', 'setting', 'gene_id',
                 'library_key', 'contrast', 'delta_DWUI', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b', 'common_four'])
    effect_fields = ['study_id', 'route', 'gene_id', 'contrast', 'n_libraries', 'contributor_libraries',
                     'weights', 'delta_DWUI', 'background', 'version', 'setting']
    routes = read('route_effects.tsv.gz', usecols=effect_fields)
    bridge = read('bridge_route_effects.tsv.gz', usecols=effect_fields)
    planned = read('library_classification.tsv', usecols=['study_id', 'route', 'library_key'])
    routekeys = ['study_id', 'route', 'version', 'setting', 'background', 'contrast']
    checks = {}
    route_summary = read('route_summary.tsv')
    bridge_summary = read('bridge_route_summary.tsv')
    ordinary = {'full': ['D00', 'D01', 'D10', 'D11'], 'common_four': ['D00', 'D01', 'D10', 'D11']}
    bridges = {'common_M_bridge': ['D00', 'D01'], 'common_P_bridge': ['D00', 'D10']}
    check_layout(route_summary, planned, ordinary)
    check_layout(bridge_summary, planned, bridges)
    checks['route_summary_rows'] = check_summaries(routes, route_summary, routekeys)
    checks['bridge_route_summary_rows'] = check_summaries(bridge, bridge_summary, routekeys)
    library_effects = pd.concat([pairs.assign(background='full'), pairs[pairs.common_four].assign(background='common_four')])
    library_summary = read('library_summary.tsv')
    check_layout(library_summary, planned, ordinary, library=True)
    checks['library_summary_rows'] = check_summaries(library_effects, library_summary, routekeys + ['library_key'])
    bridge_weights = read('bridge_pair_weights.tsv.gz', usecols=pairkeys + ['background'])
    bridge_library = bridge_weights.merge(pairs, on=['study_id', 'route', 'version', 'setting', 'gene_id', 'library_key', 'contrast'], validate='many_to_one')
    bridge_library_summary = read('bridge_library_summary.tsv')
    check_layout(bridge_library_summary, planned, bridges, library=True)
    checks['bridge_library_summary_rows'] = check_summaries(bridge_library, bridge_library_summary, routekeys + ['library_key'])
    checks['matched_aggregates'] = check_aggregates(pairs, pd.concat([routes, bridge]))
    gradient = read('gradient_route_effects.tsv.gz')
    gpair = read('gradient_pair_effects.tsv.gz')
    gkeys = ['study_id', 'route', 'panel']
    checks['gradient_route_summary_rows'] = check_summaries(gradient, read('gradient_route_summary.tsv'), gkeys, gradient=True)
    checks['gradient_library_summary_rows'] = check_summaries(gpair, read('gradient_library_summary.tsv'), gkeys + ['library_key'], gradient=True)
    for data in [gradient, gpair]:
        assert np.allclose(data.delta_high - data.delta_low, data.difference_high_minus_low, rtol=0, atol=1e-12)
        assert np.allclose(data.difference_high_minus_low, data.direct_M_low_minus_high, rtol=0, atol=1e-12)
    lost = read('../2026-10-10-marker-reference-review/shared_MC_pairs_lost.tsv')
    assert len(lost) == 15 and not lost.duplicated(pairkeys).any()
    marked = pairs.merge(lost.assign(_remove=True), on=pairkeys, how='left', validate='one_to_one')
    assert marked._remove.fillna(False).sum() == 15
    unique_pairs = marked[marked._remove.isna()].drop(columns='_remove').copy()
    jointkeys = ['study_id', 'route', 'library_key', 'version', 'setting', 'gene_id']
    unique_pairs['common_four'] = unique_pairs.groupby(jointkeys).contrast.transform('nunique').eq(4)
    unique_routes = read('unique_support_route_effects.tsv.gz', usecols=effect_fields)
    unique_bridge = read('unique_support_bridge_route_effects.tsv.gz', usecols=effect_fields)
    usummary = read('unique_support_route_summary.tsv')
    ubsummary = read('unique_support_bridge_route_summary.tsv')
    check_layout(usummary, planned, ordinary)
    check_layout(ubsummary, planned, bridges)
    checks['unique_route_summary_rows'] = check_summaries(unique_routes, usummary, routekeys)
    checks['unique_bridge_route_summary_rows'] = check_summaries(unique_bridge, ubsummary, routekeys)
    unique_library = pd.concat([unique_pairs.assign(background='full'),
                                unique_pairs[unique_pairs.common_four].assign(background='common_four')])
    ulibsummary = read('unique_support_library_summary.tsv')
    check_layout(ulibsummary, planned, ordinary, library=True)
    checks['unique_library_summary_rows'] = check_summaries(unique_library, ulibsummary, routekeys + ['library_key'])
    uweights = read('unique_support_pair_weights.tsv.gz', usecols=pairkeys + ['background'])
    ubridge_library = uweights[uweights.background.isin(bridges)].merge(unique_pairs, on=pairkeys, validate='many_to_one')
    ublibsummary = read('unique_support_bridge_library_summary.tsv')
    check_layout(ublibsummary, planned, bridges, library=True)
    checks['unique_bridge_library_summary_rows'] = check_summaries(ubridge_library, ublibsummary, routekeys + ['library_key'])
    checks['unique_matched_aggregates'] = check_aggregates(unique_pairs, pd.concat([unique_routes, unique_bridge]))
    checks['unique_removed_pairs'] = 15
    result = dict(status='PASS', checks=checks, inputs_sha256=HASHES)
    (HERE / 'result_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status='PASS', checks=checks), indent=2))


if __name__ == '__main__':
    main()
