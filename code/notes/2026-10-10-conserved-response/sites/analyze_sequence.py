"""Eight frozen exploratory conditional sequence-score comparisons."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

from analyze_sites import HERE, ROOT, read, save

B = 20000


def bins(values):
    result = pd.Series(-1, index=values.index, dtype=int)
    valid = values.notna()
    if valid.any():
        cuts = np.quantile(values[valid], [1/3, 2/3], method='linear')
        result.loc[valid] = np.searchsorted(cuts, values[valid], side='left')
    return result


def bh(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    output = np.empty(len(values))
    output[order] = np.minimum(1, np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1])
    return output


def conditional_draws(frame, column, random, n=B):
    scores = frame.score.to_numpy(float)
    labels = frame.case.to_numpy(bool)
    groups = list(frame.groupby(column, sort=True, dropna=False).indices.values())
    cases = int(labels.sum())
    result = np.zeros(n)
    expected, variance = 0., 0.
    margins = []
    for indexes in groups:
        size, count = len(indexes), int(labels[indexes].sum())
        s = scores[indexes]
        shuffled = random.permuted(np.broadcast_to(labels[indexes], (n, size)), axis=1)
        assert np.all(shuffled.sum(axis=1) == count)
        result += shuffled @ s
        expected += count / size * float(s.sum())
        if size > 1:
            variance += count * (size - count) / (size * (size - 1)) * float(np.sum((s - s.mean()) ** 2))
        margins.append(dict(stratum=frame.iloc[indexes[0]][column], genes=size, cases=count,
                            controls=size - count, mobile=0 < count < size, score_min=float(s.min()), score_max=float(s.max())))
    if cases:
        return result / cases, expected / cases, variance / cases ** 2, pd.DataFrame(margins)
    return np.full(n, np.nan), np.nan, 0., pd.DataFrame(margins)


def main():
    contexts = read(HERE / 'site_event_contexts.tsv.gz')
    contexts = contexts[contexts.source_population.isin(['main48', 'hN79']) & contexts.gene_threshold.eq(.02) & contexts.sequence_candidate]
    windows = read(HERE / 'fixed_window_metadata.tsv').set_index('id')
    membership = read(HERE.parent / 'function/comparison_membership_features.tsv.gz', low_memory=False)
    features, failures = [], []
    for row in contexts.itertuples():
        pairs = json.loads(row.sequence_driver_pairs)
        flat = sorted(set(key for pair in pairs for key in pair))
        available = all(key in windows.index for key in flat)
        if not available:
            failures.append(dict(set_name=row.set_name, tolerance=row.tolerance, human_gene_id=row.human_gene_id, contrast=row.contrast, status='missing_fixed_window_metadata'))
            continue
        selected = windows.loc[flat]
        if not selected.finite_fraction.ge(.8).all() or selected[['gc_fraction', 'phastcons_replayed']].isna().any().any():
            failures.append(dict(set_name=row.set_name, tolerance=row.tolerance, human_gene_id=row.human_gene_id, contrast=row.contrast, status='low_finite_coverage_or_unknown_GC_score'))
            continue
        score = np.mean([windows.loc[list(pair)].phastcons_replayed.mean() for pair in pairs])
        gc = np.mean([windows.loc[list(pair)].gc_fraction.mean() for pair in pairs])
        features.append(dict(set_name=row.set_name, tolerance=row.tolerance, human_gene_id=row.human_gene_id, mouse_gene_id=row.mouse_gene_id,
                             case=row.parent_case, contrast=row.contrast, score=float(score), GC=float(gc),
                             finite_fraction=float(selected.finite_fraction.min()), matched_PAS=row.matched_PAS,
                             human_span=row.human_native_span, mouse_span=row.mouse_native_span,
                             matched_fraction=min(row.human_matched_fraction, row.mouse_matched_fraction)))
    feature = pd.DataFrame(features)
    save(feature, 'sequence_contrast_features.tsv')
    save(pd.DataFrame(failures, columns=['set_name', 'tolerance', 'human_gene_id', 'contrast', 'status']), 'sequence_feature_failures.tsv')
    outputs, summaries, all_margins, draws = [], [], [], {}
    for set_name, population in [('main', 'main48'), ('hN', 'hN79')]:
        for tolerance in [25, 50]:
            selected = feature[feature.set_name.eq(set_name) & feature.tolerance.eq(tolerance)]
            frame = selected.groupby(['human_gene_id', 'mouse_gene_id', 'case'], as_index=False).agg(
                score=('score', 'mean'), GC=('GC', 'mean'), finite_fraction=('finite_fraction', 'min'),
                matched_PAS=('matched_PAS', 'median'), human_span=('human_span', 'median'), mouse_span=('mouse_span', 'median'),
                matched_fraction=('matched_fraction', 'median'), contrasts=('contrast', lambda values: ';'.join(sorted(values))))
            parent = membership[membership.set_name.eq(population) & membership.threshold.eq(.02)][['human_gene_id', 'coarse_stratum']]
            frame = frame.merge(parent, on='human_gene_id', validate='one_to_one')
            frame['GC_bin'] = bins(frame.GC)
            frame['coverage_bin'] = np.where(frame.finite_fraction.eq(1), 'complete', '0.8_to_lt1')
            frame['matched_PAS_bin'] = np.where(frame.matched_PAS.le(2), '2', np.where(frame.matched_PAS.le(4), '3_to_4', 'ge5'))
            rank = pd.concat([frame.human_span.rank(method='average', pct=True), frame.mouse_span.rank(method='average', pct=True)], axis=1)
            frame['min_span_rank'] = rank.min(axis=1)
            frame['span_bin'] = bins(frame.min_span_rank)
            frame['fraction_bin'] = bins(frame.matched_fraction)
            frame['sequence_coarse'] = frame[['coarse_stratum', 'GC_bin', 'coverage_bin', 'matched_PAS_bin']].astype(str).agg('|'.join, axis=1)
            frame['sequence_fine'] = frame[['sequence_coarse', 'span_bin', 'fraction_bin']].astype(str).agg('|'.join, axis=1)
            frame['set_name'], frame['tolerance'] = set_name, tolerance
            outputs.append(frame)
            for model in ['coarse', 'fine']:
                key = f'{set_name}|{tolerance}|{model}'
                seed = int.from_bytes(hashlib.sha256(('20261010|sequence|' + key).encode()).digest()[:8], 'little')
                null, expected, variance, margins = conditional_draws(frame, 'sequence_' + model, np.random.default_rng(seed))
                cases, controls = int(frame.case.sum()), int((~frame.case).sum())
                observed = float(frame.loc[frame.case, 'score'].mean()) if cases else np.nan
                mobile_cases = int(margins.loc[margins.mobile, 'cases'].sum())
                valid = cases >= 3 and variance > 1e-24
                tail = int(np.sum(null >= observed - 1e-12)) if cases else B
                p = (tail + 1) / (B + 1) if valid else 1.
                summaries.append(dict(set_name=set_name, tolerance=tolerance, model=model, genes=len(frame), cases=cases, controls=controls,
                                      observed_mean=observed, conditional_expected=expected, null_low=float(np.nanquantile(null, .025)) if cases else np.nan,
                                      null_high=float(np.nanquantile(null, .975)) if cases else np.nan, null_sd=float(np.sqrt(variance)),
                                      Monte_Carlo_mean=float(np.nanmean(null)) if cases else np.nan, mobile_cases=mobile_cases,
                                      mobile_case_fraction=mobile_cases / cases if cases else np.nan, strata=len(margins),
                                      mobile_strata=int(margins.mobile.sum()), seed=seed, permutations=B, upper_tail_count=tail,
                                      p=p, estimable=valid, status='estimable' if valid else 'NE_case_count_lt3' if cases < 3 else 'NE_zero_conditional_variance'))
                margins['set_name'], margins['tolerance'], margins['model'] = set_name, tolerance, model
                all_margins.append(margins); draws[key] = null
    result = pd.DataFrame(summaries)
    assert len(result) == 8
    result['q8'] = bh(result.p)
    save(pd.concat(outputs, ignore_index=True), 'sequence_gene_features_and_strata.tsv')
    save(pd.concat(all_margins, ignore_index=True), 'sequence_stratum_margins.tsv')
    save(result, 'sequence_tests.tsv')
    np.savez_compressed(HERE / 'sequence_null_draws.npz', **draws)
    print(result.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
