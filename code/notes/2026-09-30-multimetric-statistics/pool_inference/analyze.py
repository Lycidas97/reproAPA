"""Frozen pair, independent day-by-pool P--M metric sign inference."""
from pathlib import Path
import hashlib
import json
import sys
import platform

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import extract as frozen

METRICS = ['DWUI', 'RWUI', 'PDUI', 'PPUI', 'MPRO']
FOLDER = 'gse221739_human_pluripotency_reprogramming'
INPUTS = {}


def digest(path):
    return hashlib.file_digest(Path(path).open('rb'), 'sha256').hexdigest()


def read(path, **kwargs):
    path = Path(path)
    INPUTS[str(path)] = dict(sha256=digest(path), bytes=path.stat().st_size)
    return pd.read_csv(path, sep='\t', low_memory=False, **kwargs)


def save(data, name):
    data.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def dump(value, name):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(q, 1)
    return result


def sign_statistics(delta, weights, linear_weights):
    n = len(delta)
    patterns = np.arange(1 << n)
    signs = 2 * ((patterns[:, None] >> np.arange(n)) & 1) - 1
    mean_delta = np.asarray(weights) @ delta
    transformed = signs @ (delta * np.asarray(weights)[:, None])
    effects = np.r_[mean_delta @ linear_weights, np.ptp(mean_delta)]
    observed = np.r_[np.abs(effects[:4]), effects[4]]
    null = np.column_stack((np.abs(transformed @ linear_weights), np.ptp(transformed, axis=1)))
    extreme = (null >= observed - 1e-12).sum(0)
    p = extreme / len(patterns)
    assert np.all(p >= 2 / len(patterns))
    assert np.allclose(delta.sum(1), 0, atol=1e-12)
    assert np.allclose(transformed.sum(1), 0, atol=1e-12)
    return effects, p, mean_delta, null, extreme


def freeze_design():
    protocol = HERE / 'PROTOCOL.zh.md'
    INPUTS[str(protocol)] = dict(sha256=digest(protocol), bytes=protocol.stat().st_size)
    pairs = read(HERE.parent / 'pair_effects.tsv.gz')
    pairs = pairs[pairs.study_id.eq('GSE221739') & pairs.contrast.eq('P_M')
                  & pairs.design.eq('same_library')].copy()
    assert not pairs.duplicated(['gene_id', 'library_key_a', 'library_key_b']).any()
    assert pairs.library_key_a.eq(pairs.library_key_b).all()
    assert pairs.sample_unit_a.eq(pairs.sample_unit_b).all() and pairs.time_a.eq(pairs.time_b).all()
    assert pairs.replicate_status_a.eq('pooled_chip_sample').all()
    assert pairs.replicate_status_b.eq('pooled_chip_sample').all()
    design = read(HERE.parents[2] / 'notes/2026-09-26-endpoint-counting/inputs/sampling_design.tsv')
    design = design[design.study_id.eq('GSE221739')].copy()
    aliases = design.groupby('library_key').agg(sample_unit=('sample_unit', 'first'),
        time=('time', 'first'), replicate_status=('replicate_status', 'first'),
        states=('state', lambda x: ';'.join(sorted(set(x)))), source_url=('source_url', 'first')).reset_index()
    assert design.groupby('library_key').sample_unit.nunique().max() == 1
    assert aliases.groupby('sample_unit').time.nunique().max() == 1
    aliases['in_frozen_PM_pair_set'] = aliases.library_key.isin(pairs.library_key_a)
    for suffix in ['a', 'b']:
        link = pairs.merge(aliases, left_on='library_key_' + suffix, right_on='library_key', validate='many_to_one')
        assert link['sample_unit_' + suffix].eq(link.sample_unit).all()
        assert link['time_' + suffix].eq(link.time).all()
    # Capture aggregation occurs before counting experimental units.
    units = pairs[['gene_id', 'time_a', 'sample_unit_a']].drop_duplicates()
    duplicate_dates = units.groupby(['gene_id', 'time_a']).sample_unit_a.nunique().rename('pools').reset_index()
    coverage = pairs.groupby('gene_id').agg(gene_symbol=('gene_symbol', 'first'),
        n_original_capture_pairs=('library_key_a', 'size'), n_independent_pools=('sample_unit_a', 'nunique'),
        n_dates=('time_a', 'nunique')).reset_index()
    duplicates = duplicate_dates[duplicate_dates.pools.ge(2)].groupby('gene_id').time_a.nunique()
    coverage['n_duplicate_pool_dates'] = coverage.gene_id.map(duplicates).fillna(0).astype(int)
    coverage['tested'] = coverage.n_independent_pools.ge(4) & coverage.n_duplicate_pool_dates.ge(2)
    coverage['exclusion_reason'] = np.where(coverage.tested, 'none',
        np.where(coverage.n_independent_pools.lt(4), 'fewer_than_4_independent_pools', 'fewer_than_2_duplicate_pool_dates'))
    save(pairs, 'frozen_original_pairs.tsv.gz')
    save(aliases, 'pool_capture_aliases.tsv')
    save(coverage, 'coverage_and_exclusions.tsv')
    dump(dict(stage='FROZEN_BEFORE_NEW_METRIC_TESTS', original_pairs=len(pairs),
              original_genes=len(coverage), testable_genes=int(coverage.tested.sum()),
              capture_libraries=pairs.library_key_a.nunique(), independent_pools=pairs.sample_unit_a.nunique(),
              dates=sorted(pairs.time_a.unique()), protocol_sha256=digest(protocol), inputs=INPUTS), 'freeze.json')
    return pairs, aliases, coverage


def main():
    pairs, aliases, coverage = freeze_design()
    selected = pairs[pairs.gene_id.isin(coverage.loc[coverage.tested, 'gene_id'])].copy()
    vectors, states, features, groups, axes, state_error = frozen.study_vectors(FOLDER, selected)
    save(states, 'reconstructed_capture_state_profiles.tsv.gz')
    save(features[['feature_index_1based', 'atlas_pas_id', 'primary_gene_id', 'cut_site', 'strand',
                   'position', 'rank', 'n', 'WPP', 'RWUI', 'PDUI', 'PPUI']], 'frozen_PAS_axes.tsv.gz')
    test_rows, pool_rows, reconstruction = [], [], []
    for gene, w in selected.groupby('gene_id', sort=True):
        ix = groups[gene]
        linear = features[['WPP', 'RWUI', 'PDUI', 'PPUI']].iloc[ix].to_numpy()
        pool_deltas, pool_meta = [], []
        for (date, pool), captures in w.groupby(['time_a', 'sample_unit_a'], sort=True):
            differences = []
            for row in captures.itertuples():
                a = vectors[(row.route, row.library_key_a, 'P', gene)]
                b = vectors[(row.route, row.library_key_b, 'M', gene)]
                joint = int(((a[1] + b[1] >= 10) & (a[2] + b[2] >= 2)).sum())
                assert joint == row.n_joint_PAS and joint >= 2
                d = a[0] - b[0]
                assert abs(d @ linear[:, 0] - row.APA_delta) < 1e-8
                assert abs(d @ linear[:, 1] - row.delta_RWUI) < 1e-8
                differences.append(d)
                reconstruction.append(dict(gene_id=gene, library_key=row.library_key_a,
                    DWUI_error=float(abs(d @ linear[:, 0] - row.APA_delta)),
                    RWUI_error=float(abs(d @ linear[:, 1] - row.delta_RWUI)), n_joint_PAS=joint))
            delta = np.mean(differences, axis=0)
            pool_deltas.append(delta)
            pool_meta.append(dict(gene_id=gene, gene_symbol=w.gene_symbol.iloc[0], time=date, sample_unit=pool,
                                  n_capture_pairs=len(captures), captures=';'.join(sorted(captures.library_key_a)),
                                  delta_usage=json.dumps(delta.tolist(), separators=(',', ':')),
                                  PAS_ids=json.dumps(features.atlas_pas_id.iloc[ix].tolist(), separators=(',', ':'))))
        meta = pd.DataFrame(pool_meta)
        weights = 1 / meta.time.nunique() / meta.groupby('time').sample_unit.transform('size').to_numpy()
        assert abs(weights.sum() - 1) < 1e-12
        n = len(meta)
        assert meta.sample_unit.nunique() == n and 4 <= n <= 6
        delta = np.stack(pool_deltas)
        effects, p, mean_delta, null, extreme = sign_statistics(delta, weights, linear)
        meta['weight'] = weights
        meta['PAS_number'] = len(ix)
        for j, metric in enumerate(METRICS[:4]):
            meta['delta_' + metric] = delta @ linear[:, j]
        meta['MPRO_per_pool_descriptive'] = np.ptp(delta, axis=1)
        pool_rows.extend(meta.to_dict('records'))
        for j, metric in enumerate(METRICS):
            test_rows.append(dict(gene_id=gene, gene_symbol=w.gene_symbol.iloc[0], metric=metric,
                estimate=float(effects[j]), p=float(p[j]), extreme_patterns=int(extreme[j]),
                full_sign_patterns=1 << n, minimum_two_sided_p=2 / (1 << n),
                n_independent_pools=n, n_dates=meta.time.nunique(),
                n_duplicate_pool_dates=int(meta.groupby('time').sample_unit.nunique().ge(2).sum()),
                n_capture_pairs=len(w), PAS_number=len(ix),
                dates=';'.join(sorted(meta.time.unique())), pools=';'.join(meta.sample_unit),
                null_statistic_min=float(null[:, j].min()), null_statistic_max=float(null[:, j].max()),
                balanced_delta_usage=json.dumps(mean_delta.tolist(), separators=(',', ':'))))
        if len(test_rows) % 1000 == 0:
            print('tested genes', len(test_rows) // 5, flush=True)
    result = pd.DataFrame(test_rows)
    assert len(result) == 5 * int(coverage.tested.sum())
    result['q_BH_metric_family'] = result.groupby('metric').p.transform(lambda x: bh(x.to_numpy()))
    result['q_BH_5metric_family'] = bh(result.p.to_numpy())
    calls = []
    for row in result.itertuples():
        thresholds = [.10, .20] if row.metric == 'MPRO' else [.02, .05]
        for threshold in thresholds:
            effect = row.estimate if row.metric == 'MPRO' else abs(row.estimate)
            reached = effect >= threshold
            calls.append(dict(gene_id=row.gene_id, metric=row.metric, threshold=threshold,
                effect_exceeds_threshold=reached, raw_p05=row.p < .05,
                metric_family_q05=row.q_BH_metric_family < .05,
                combined_family_q05=row.q_BH_5metric_family < .05,
                effect_and_raw_p05=reached and row.p < .05,
                effect_and_metric_family_q05=reached and row.q_BH_metric_family < .05,
                effect_and_combined_family_q05=reached and row.q_BH_5metric_family < .05))
    save(result, 'gene_metric_tests.tsv.gz')
    save(pd.DataFrame(calls), 'effect_threshold_test_calls.tsv.gz')
    save(pd.DataFrame(pool_rows), 'independent_pool_vectors.tsv.gz')
    save(pd.DataFrame(reconstruction), 'pair_reconstruction.tsv.gz')
    summary = result.groupby('metric').agg(genes=('gene_id', 'size'), raw_p05=('p', lambda x: int(x.lt(.05).sum())),
        metric_family_q05=('q_BH_metric_family', lambda x: int(x.lt(.05).sum())),
        combined_family_q05=('q_BH_5metric_family', lambda x: int(x.lt(.05).sum())),
        minimum_p=('p', 'min'), minimum_metric_q=('q_BH_metric_family', 'min'),
        minimum_combined_q=('q_BH_5metric_family', 'min')).reset_index()
    save(summary, 'summary.tsv')
    theoretical = result[result.metric.eq('DWUI')]
    dump(dict(status='COMPLETE', original_genes=len(coverage), tested_genes=int(coverage.tested.sum()),
              tested_gene_metric_rows=len(result), primary_family_size=len(theoretical), combined_family_size=len(result),
              minimum_theoretical_metric_BH=float(bh(theoretical.minimum_two_sided_p.to_numpy()).min()),
              input_sha256=INPUTS, extraction_inputs=frozen.INPUTS,
              protocol_sha256=digest(HERE / 'PROTOCOL.zh.md'), script_sha256=digest(Path(__file__)),
              state_reconstruction_max_error=float(state_error),
              pool_size_counts={str(k): int(v) for k, v in theoretical.n_independent_pools.value_counts().items()},
              versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__),
              inference_scope='same BJ cell line independent day-by-pool observed-state association',
              assumptions=['Independent day-by-pool PAS difference vectors symmetric about zero under null',
                           'Frozen dates fixed and equal-weighted; pools equal-weighted within dates',
                           'Capture aliases within one pool are not independent experimental replicates']), 'audit.json')
    print(summary.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
