"""Fixed-bin, class-conditional ortholog direction correspondence.

Classes and their input hashes must be frozen before execution. This analysis
does not infer a conserved PAS, isoform, trans regulator, or causal transition.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCAN = ROOT / 'notes/2026-10-10-ortholog-threshold-scan'
OLD = ROOT / 'notes/2026-10-10-recurrence-utr-ortholog'
spec = importlib.util.spec_from_file_location('fixed_threshold_scan', SCAN / 'analyze.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
A = S.A
LABELS = A.LABELS
KEY = ['human_gene_id', 'mouse_gene_id']
THRESHOLDS = [.01, .015, .02, .025, .03, .035, .04, .045, .05]
MODELS = ['measurement_PAS_capacity', 'full_assay']
CLASSES = ['both_tandem', 'both_ALE']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(frame, name):
    compression = {'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA',
                 float_format='%.17g', compression=compression)


def independently_replay_labels(events, members):
    """A third opposing study remains allowed, as in the frozen recurrence rule."""
    for threshold in THRESHOLDS:
        e = events[events.threshold.eq(threshold)]
        m = members[members.threshold.eq(threshold)]
        assert not e.duplicated(['species', 'unit', 'gene_id']).any()
        assert np.array_equal(e.event.to_numpy(), np.where(e.effect_max.le(-threshold),
            -1, np.where(e.effect_min.ge(threshold), 1, 0)))
        for species in ['human', 'mouse']:
            ee = e[e.species.eq(species)]
            counts = ee.groupby('gene_id').event.apply(
                lambda x: (int(x.eq(-1).sum()), int(x.eq(1).sum())))
            ids = m[species + '_gene_id']
            negative = np.array([counts[g][0] >= 2 for g in ids])
            positive = np.array([counts[g][1] >= 2 for g in ids])
            assert not (negative & positive).any()
            expected = np.where(negative, 'proximal', np.where(positive, 'distal', 'nonrecurrent'))
            assert np.array_equal(expected, m[species + '_direction'].to_numpy())


def correction(frame):
    f = frame.copy()
    primary = f.human_direction.eq('proximal') & f.mouse_direction.eq('proximal')
    for column in ['q_primary4', 'q_sensitivity4', 'q_scan36',
                   'q_secondary16_per_threshold', 'q_secondary_scan144']:
        f[column] = np.nan
    for threshold, column in [(.02, 'q_primary4'), (.05, 'q_sensitivity4')]:
        idx = f.index[primary & f.threshold.eq(threshold)]
        assert len(idx) == 4
        f.loc[idx, column] = A.bh(f.loc[idx, 'p_enrichment'])
    idx = f.index[primary]
    assert len(idx) == 36
    f.loc[idx, 'q_scan36'] = A.bh(f.loc[idx, 'p_enrichment'])
    directional = f.human_direction.isin(LABELS[:2]) & f.mouse_direction.isin(LABELS[:2])
    for _, idx in f[directional].groupby('threshold').groups.items():
        assert len(idx) == 16
        f.loc[idx, 'q_secondary16_per_threshold'] = A.bh(f.loc[idx, 'p_enrichment'])
    idx = f.index[directional]
    assert len(idx) == 144
    f.loc[idx, 'q_secondary_scan144'] = A.bh(f.loc[idx, 'p_enrichment'])
    return f


def select_features(features, pairs, species, threshold):
    out = features[features.species.eq(species) & features.threshold.eq(threshold)]
    assert out.gene_id.is_unique
    out = out.set_index('gene_id').reindex(pairs[species + '_gene_id'].to_list())
    out.index.name = 'gene_id'
    out = out.reset_index()
    assert out.notna().all().all()
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--protocol', type=Path, required=True)
    p.add_argument('--classification', type=Path, required=True)
    p.add_argument('--strict-eligibility', type=Path)
    args = p.parse_args()
    protocol_path = args.protocol.resolve()
    classification_path = args.classification.resolve()
    protocol = json.loads(protocol_path.read_text())
    assert protocol['thresholds'] == THRESHOLDS
    assert protocol['models'] == MODELS
    assert protocol['classes'] == CLASSES
    assert protocol['permutations'] == 100000
    assert protocol['classification_sha256'] == digest(classification_path)
    paths = [protocol_path, classification_path, SCAN / 'ortholog_assay_features.tsv.gz',
        SCAN / 'measured_ortholog_members.tsv.gz', SCAN / 'study_events.tsv.gz',
        SCAN / 'analyze.py', OLD / 'analyze.py',
        ROOT / 'notes/2026-10-10-directional-recurrence/analyze.py',
        ROOT / 'notes/2026-10-02-null-model-controls/recurrence.py',
        ROOT / 'notes/2026-10-02-null-model-controls/common.py',
        ROOT / 'analysis_code/pas_rerun/downstream_runtime.py']
    if 'parent_protocol_path' in protocol:
        parent = ROOT / protocol['parent_protocol_path']
        assert protocol['parent_protocol_sha256'] == digest(parent)
        paths.append(parent)
    if args.strict_eligibility:
        strict_path = args.strict_eligibility.resolve()
        assert protocol['strict_eligibility_sha256'] == digest(strict_path)
        paths.append(strict_path)
    hashes = {str(x.relative_to(ROOT)): digest(x) for x in paths}
    classes = pd.read_csv(classification_path, sep='\t')
    required = KEY + ['human_event_class', 'mouse_event_class', 'pair_class']
    assert set(required).issubset(classes.columns)
    assert len(classes) == 751 and not classes.duplicated(KEY).any()
    features = pd.read_csv(paths[2], sep='\t', float_precision='round_trip')
    features = features[features.contrast.eq('P_M')].copy()
    members = pd.read_csv(paths[3], sep='\t', float_precision='round_trip')
    members = members[members.contrast.eq('P_M')].copy()
    events = pd.read_csv(paths[4], sep='\t', float_precision='round_trip')
    events = events[events.contrast.eq('P_M')].copy()
    assert set(members.threshold) == set(THRESHOLDS)
    independently_replay_labels(events, members)
    canonical = members[members.threshold.eq(.02)][KEY].reset_index(drop=True)
    assert set(map(tuple, canonical.to_numpy())) == set(map(tuple, classes[KEY].to_numpy()))
    for threshold in THRESHOLDS:
        assert members[members.threshold.eq(threshold)][KEY].reset_index(drop=True).equals(canonical)
    # Eligibility does not alter bins: these bins were fixed in all 751 pairs.
    bins = ['PAS_number_bin', 'mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']
    for species in ['human', 'mouse']:
        fixed = select_features(features, canonical, species, .02)
        for threshold in THRESHOLDS:
            current = select_features(features, canonical, species, threshold)
            assert current[bins].equals(fixed[bins])
    joined = members.merge(classes[required], on=KEY, validate='many_to_one')
    results, diagnostics, mobility, supports, archives = [], [], [], [], {}
    save(joined, 'classified_ortholog_members.tsv.gz')
    for event_class in CLASSES:
        for threshold in THRESHOLDS:
            pairs = joined[joined.pair_class.eq(event_class) & joined.threshold.eq(threshold)].copy()
            pairs = pairs.reset_index(drop=True)
            assert len(pairs) > 0
            context = dict(contrast='P_M', pair_class=event_class,
                threshold=threshold, measured_orthologs=len(pairs))
            feats = {sp: select_features(features, pairs, sp, threshold) for sp in ['human', 'mouse']}
            labels = {sp: pairs[sp + '_direction'].map(dict(zip(LABELS, range(3)))).to_numpy(np.int8)
                      for sp in ['human', 'mouse']}
            h, m = labels['human'], labels['mouse']
            for model in MODELS:
                cells = {sp: A.partitions(feats[sp], model) for sp in ['human', 'mouse']}
                prob = {sp: A.class_probabilities(labels[sp], cells[sp]) for sp in ['human', 'mouse']}
                analytic = prob['human'].T @ prob['mouse']
                run_seed = A.seed('event_stratified_direction', 'P_M', event_class, threshold, model)
                null = A.cross_draws(h, m, cells['human'], cells['mouse'],
                    np.random.default_rng(run_seed), n=protocol['permutations'])
                archives['|'.join(map(str, (event_class, threshold, model)))] = null
                assert np.all(null.sum((1, 2)) == len(pairs))
                for a in range(3):
                    assert np.all(null.sum(2)[:, a] == int((h == a).sum()))
                    assert np.all(null.sum(1)[:, a] == int((m == a).sum()))
                for sp in ['human', 'mouse']:
                    diag, summary = S.feature_diagnostics(feats[sp], labels[sp], cells[sp],
                        context | dict(model=model, species=sp))
                    diagnostics.append(diag)
                    mobility.append(summary)
                for a, hd in enumerate(LABELS):
                    for b, md in enumerate(LABELS):
                        values = null[:, a, b]
                        row = context | dict(model=model, human_direction=hd,
                            mouse_direction=md, seed=run_seed,
                            **A.null_row(int(((h == a) & (m == b)).sum()), values, analytic[a, b]))
                        row.update(S.tail_precision(row['tail_count'], len(values)))
                        row['directional_test'] = a < 2 and b < 2
                        row['null_support_points'] = int(np.unique(values).size)
                        if row['directional_test']:
                            support, counts = np.unique(values, return_counts=True)
                            supports.extend(context | dict(model=model, human_direction=hd,
                                mouse_direction=md, overlap_count=int(v), draws=int(n),
                                probability=float(n / len(values)))
                                for v, n in zip(support, counts))
                        if not row['directional_test']:
                            row['p_enrichment'] = np.nan
                        results.append(row)
                print(event_class, threshold, model, 'fixed background', len(pairs), flush=True)
    result = correction(pd.DataFrame(results))
    save(result, 'class_direction_3x3.tsv')
    save(result[result.directional_test], 'directional_tests.tsv')
    save(result[result.human_direction.eq('proximal') & result.mouse_direction.eq('proximal')],
         'primary_proximal_overlap.tsv')
    save(pd.concat(diagnostics, ignore_index=True), 'stratum_diagnostics.tsv.gz')
    save(pd.DataFrame(mobility), 'mobility_summary.tsv')
    save(pd.DataFrame(supports), 'null_support.tsv')
    np.savez_compressed(HERE / 'null_draws.npz', **archives)
    # A stricter all-sites qualification is ancillary, with no new p-value family.
    if args.strict_eligibility:
        strict = pd.read_csv(strict_path, sep='\t')
        assert not strict.duplicated(KEY).any()
        assert 'both_full_axis_tandem' in strict
        assert strict.both_full_axis_tandem.notna().all()
        eligible = strict.loc[strict.both_full_axis_tandem & strict.pair_class.eq('both_tandem'), KEY]
        q = joined.merge(eligible, on=KEY, validate='many_to_one')
        rows, strict_mobility = [], []
        for threshold in THRESHOLDS:
            pairs = q[q.threshold.eq(threshold)].reset_index(drop=True)
            assert pairs.pair_class.eq('both_tandem').all()
            for model in MODELS:
                fs = {sp: select_features(features, pairs, sp, threshold) for sp in ['human', 'mouse']}
                ls = {sp: pairs[sp + '_direction'].map(dict(zip(LABELS, range(3)))).to_numpy(np.int8)
                      for sp in ['human', 'mouse']}
                cs = {sp: A.partitions(fs[sp], model) for sp in ['human', 'mouse']}
                ps = {sp: A.class_probabilities(ls[sp], cs[sp]) for sp in ['human', 'mouse']}
                expected = ps['human'].T @ ps['mouse']
                for sp in ['human', 'mouse']:
                    _, summary = S.feature_diagnostics(fs[sp], ls[sp], cs[sp],
                        dict(contrast='P_M', pair_class='both_tandem_full_axis',
                             threshold=threshold, model=model, species=sp))
                    strict_mobility.append(summary)
                for a, hd in enumerate(LABELS):
                    for b, md in enumerate(LABELS):
                        observed = int(((ls['human'] == a) & (ls['mouse'] == b)).sum())
                        rows.append(dict(threshold=threshold, model=model,
                            measured_orthologs=len(pairs), human_direction=hd, mouse_direction=md,
                            observed=observed, analytic_expected=expected[a, b],
                            excess=observed - expected[a, b],
                            fold=observed / expected[a, b] if expected[a, b] else np.nan))
        save(pd.DataFrame(rows), 'full_axis_tandem_descriptive.tsv')
        save(q, 'full_axis_tandem_members.tsv.gz')
        save(pd.DataFrame(strict_mobility), 'full_axis_tandem_mobility.tsv')
    assert all(digest(ROOT / path) == value for path, value in hashes.items())
    outputs = ['classified_ortholog_members.tsv.gz', 'class_direction_3x3.tsv',
        'directional_tests.tsv', 'primary_proximal_overlap.tsv', 'stratum_diagnostics.tsv.gz',
        'mobility_summary.tsv', 'null_draws.npz', 'null_support.tsv']
    if args.strict_eligibility:
        outputs += ['full_axis_tandem_descriptive.tsv', 'full_axis_tandem_members.tsv.gz',
                    'full_axis_tandem_mobility.tsv']
    validation = dict(status='PASS', fixed_pair_inventory=751,
        class_backgrounds=classes.pair_class.value_counts().to_dict(),
        fixed_bins_preserved=True, recurrence_labels_independently_replayed=True,
        permutations=100000, null_contexts=len(archives),
        exact_global_and_stratum_direction_margins=True, analytic_means_MC_checked=True,
        primary_family=4, sensitivity_family=4, scan_family=36,
        secondary_family_per_threshold=16, secondary_scan_family=144,
        descriptive_non_directional_cells_have_no_p_values=True,
        inference_scope='Class-conditional gene correspondence; not conserved PAS or causal mechanism.')
    (HERE / 'analysis_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    outputs.append('analysis_validation.json')
    manifest = dict(status='PASS', input_sha256=hashes,
        output_sha256={n: digest(HERE / n) for n in outputs},
        script_sha256=digest(Path(__file__)), python_executable=sys.executable,
        python=sys.version, numpy=np.__version__, pandas=pd.__version__)
    (HERE / 'analysis_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
