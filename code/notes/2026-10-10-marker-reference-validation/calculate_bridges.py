"""Independently calculate both matched bridges from sealed library deltas.

Linearity permits scalar reconstruction without importing production code or
reconstructing PAS profiles a second time. All prespecified settings are kept.
"""
from collections import Counter
from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / '2026-10-10-marker-reference-sensitivity'
REVIEW = HERE.parent / '2026-10-10-marker-reference-review'
KEYS = ['study_id', 'route', 'version', 'setting', 'gene_id']
PAIRKEYS = ['study_id', 'route', 'library_key', 'setting', 'version', 'contrast', 'gene_id']
BRIDGES = {'common_M_bridge': ['D00', 'D01'], 'common_P_bridge': ['D00', 'D10']}


def main():
    pairs = pd.read_csv(SOURCE / 'pair_effects.tsv.gz', sep='\t', float_precision='round_trip',
        usecols=PAIRKEYS + ['delta_DWUI', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b'])
    assert not pairs.duplicated(PAIRKEYS).any()
    assert pairs.time_a.eq(pairs.time_b).all() and pairs.sample_unit_a.eq(pairs.sample_unit_b).all()
    lost = pd.read_csv(REVIEW / 'shared_MC_pairs_lost.tsv', sep='\t')
    lostkeys = set(lost[PAIRKEYS].itertuples(index=False, name=None))
    assert len(lostkeys) == 15
    metadata, groups, qualified_lost = {}, {}, set()
    for row in pairs.itertuples(index=False):
        library = (row.study_id, row.route, row.library_key)
        value = (row.time_a, row.sample_unit_a)
        assert library not in metadata or metadata[library] == value
        metadata[library] = value
        key = tuple(getattr(row, k) for k in KEYS)
        token = tuple(getattr(row, k) for k in PAIRKEYS)
        remove = token in lostkeys
        if remove:
            qualified_lost.add(token)
        groups.setdefault(key, {}).setdefault(row.contrast, {})[row.library_key] = (row.delta_DWUI, remove)
    assert qualified_lost == lostkeys
    cache, effects = {}, []
    for key, contrasts in groups.items():
        for rule in ['frozen_group_contribution', 'unique_MC_support']:
            for background, required in BRIDGES.items():
                subsets = [{lib for lib, (delta, remove) in contrasts.get(c, {}).items()
                            if rule == 'frozen_group_contribution' or not remove} for c in required]
                libraries = sorted(set.intersection(*subsets))
                if not libraries:
                    continue
                signature = (key[0], key[1], tuple(libraries))
                if signature not in cache:
                    units = [metadata[(key[0], key[1], lib)] for lib in libraries]
                    dates = {d for d, u in units}
                    counts = Counter(units)
                    per_date = {d: len({u for dd, u in units if d == dd}) for d in dates}
                    cache[signature] = np.array([1 / (len(dates) * per_date[d] * counts[(d, u)]) for d, u in units])
                weights = cache[signature]
                for c in required:
                    estimate = weights @ np.array([contrasts[c][lib][0] for lib in libraries])
                    effects.append(dict(zip(KEYS, key), rule=rule, background=background, contrast=c,
                        delta_DWUI=float(estimate), contributor_libraries=json.dumps(libraries, separators=(',', ':')),
                        weights=json.dumps(weights.tolist(), separators=(',', ':')), n_libraries=len(libraries)))
    frame = pd.DataFrame(effects)
    frame.to_csv(HERE / 'independent_bridge_route_effects.tsv.gz', sep='\t', index=False,
                 float_format='%.17g', compression={'method': 'gzip', 'mtime': 0})
    plans = pd.read_csv(SOURCE / 'library_classification.tsv', sep='\t')[['study_id', 'route']].drop_duplicates()
    indexed = {key: q for key, q in frame.groupby(['rule', 'background', 'study_id', 'route', 'version', 'setting', 'contrast'])}
    settings = [('supported50', v) for v in ['core_third_strict', 'core_quartile_strict', 'core_third_relaxed', 'native_third_strict']] + [('expanded_primary', 'core_third_strict')]
    summaries = []
    for sid, route in plans.itertuples(index=False, name=None):
        for rule in ['frozen_group_contribution', 'unique_MC_support']:
            for setting, version in settings:
                for background, contrasts in BRIDGES.items():
                    for c in contrasts:
                        q = indexed.get((rule, background, sid, route, version, setting, c), frame.iloc[:0])
                        delta = q.delta_DWUI.to_numpy()
                        libraries = set()
                        for item in q.contributor_libraries:
                            libraries.update(json.loads(item))
                        for threshold in [.02, .05]:
                            summaries.append(dict(study_id=sid, route=route, rule=rule, setting=setting, version=version,
                                background=background, contrast=c, threshold=threshold, N=len(q),
                                n_proximal=int(np.sum(delta <= -threshold)), n_distal=int(np.sum(delta >= threshold)),
                                fraction_proximal=float(np.mean(delta <= -threshold)) if len(q) else np.nan,
                                fraction_distal=float(np.mean(delta >= threshold)) if len(q) else np.nan,
                                mean_delta=float(delta.mean()) if len(q) else np.nan,
                                n_contributor_libraries=len(libraries)))
    out = pd.DataFrame(summaries)
    out.to_csv(HERE / 'independent_bridge_route_summary.tsv', sep='\t', index=False, float_format='%.17g', na_rep='NA')
    print(out[out.rule.eq('unique_MC_support') & out.setting.eq('supported50') & out.version.eq('core_third_strict') & out.background.eq('common_M_bridge') & out.threshold.eq(.02)].to_string(index=False))


if __name__ == '__main__':
    main()
