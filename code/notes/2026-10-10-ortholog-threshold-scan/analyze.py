"""Threshold scan on the existing fixed ortholog inventory and class null."""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'notes/2026-10-10-recurrence-utr-ortholog'
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
SPEC = importlib.util.spec_from_file_location('frozen_ortholog_direction', SOURCE / 'analyze.py')
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)
KEY = ['species', 'contrast', 'unit', 'gene_id']
LABELS = A.LABELS


def save(frame, name):
    compression = {'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA',
                 float_format='%.17g', compression=compression)


def threshold_events(base, threshold):
    """The min/max summarize every evaluable series; the mean is not a vote."""
    f = base.copy()
    f['threshold'] = threshold
    positive = f.effect_min.ge(threshold)
    negative = f.effect_max.le(-threshold)
    assert not (positive & negative).any()
    f['event'] = np.where(negative, -1, np.where(positive, 1, 0)).astype(np.int8)
    # This tolerance belongs only to the reviewed capacity guard, not the event.
    f['can_negative'] = f.capacity_negative.ge(threshold - 1e-9)
    f['can_positive'] = f.capacity_positive.ge(threshold - 1e-9)
    assert f.loc[f.event.eq(-1), 'can_negative'].all()
    assert f.loc[f.event.eq(1), 'can_positive'].all()
    return f


def apply_correction(frame):
    """The scan family includes all eight hypotheses, without threshold selection."""
    frame = frame.copy()
    frame.loc[~frame.tested, 'p_enrichment'] = np.nan
    for column in ['q_eight_per_threshold_model', 'q_scan_72_per_model',
                   'q_scan_144_across_models']:
        frame[column] = np.nan
    for _, ii in frame[frame.tested].groupby(['threshold', 'model']).groups.items():
        assert len(ii) == 8
        frame.loc[ii, 'q_eight_per_threshold_model'] = A.bh(frame.loc[ii, 'p_enrichment'])
    for _, ii in frame[frame.tested].groupby('model').groups.items():
        assert len(ii) == 72
        frame.loc[ii, 'q_scan_72_per_model'] = A.bh(frame.loc[ii, 'p_enrichment'])
    ii = frame.index[frame.tested]
    assert len(ii) == 144
    frame.loc[ii, 'q_scan_144_across_models'] = A.bh(frame.loc[ii, 'p_enrichment'])
    return frame


def tail_precision(tail_count, draws):
    """Sampling uncertainty of the null upper-tail probability, including zero tails."""
    rate = tail_count / draws
    z = 1.959963984540054
    denominator = 1 + z*z/draws
    center = (rate + z*z/(2*draws))/denominator
    radius = z*np.sqrt(rate*(1-rate)/draws + z*z/(4*draws*draws))/denominator
    return dict(p_tail_MC_SE=float(np.sqrt(rate*(1-rate)/draws)),
        p_tail_Wilson95_lo=float(max(0, center-radius)),
        p_tail_Wilson95_hi=float(min(1, center+radius)), permutation_p_floor=1/(draws+1))


def feature_diagnostics(features, labels, cells, context):
    rows = []
    for j, ii in enumerate(cells):
        part = features.iloc[ii]
        counts = np.bincount(labels[ii], minlength=3)
        row = context | dict(stratum=j, n=len(ii), proximal=int(counts[0]),
            distal=int(counts[1]), nonrecurrent=int(counts[2]), singleton=len(ii) == 1,
            homogeneous=np.count_nonzero(counts) == 1,
            expected_changed_fraction=1-float((counts/len(ii)) @ (counts/len(ii))),
            gene_ids=';'.join(part.gene_id))
        for column in ['measured_pattern', 'proximal_capacity', 'distal_capacity',
                       'PAS_number_bin', 'mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']:
            row[column + '_levels'] = ';'.join(map(str, sorted(part[column].unique())))
        for column in ['PAS_number', 'mean_GEX', 'APA_coverage', 'baseline']:
            row[column + '_min'] = float(part[column].min())
            row[column + '_max'] = float(part[column].max())
        rows.append(row)
    out = pd.DataFrame(rows)
    total = int(out.n.sum())
    summary = context | dict(measured_genes=total, n_strata=len(out),
        n_singleton_strata=int(out.singleton.sum()), largest_stratum=int(out.n.max()),
        singleton_gene_fraction=float(out.loc[out.singleton, 'n'].sum()/total),
        frozen_label_fraction=float(out.loc[out.homogeneous, 'n'].sum()/total),
        expected_changed_fraction=float((out.n * out.expected_changed_fraction).sum()/total))
    return out, summary


def paired_roster(pairs, part, labels, threshold):
    out = pairs.copy()
    for sp in ['human', 'mouse']:
        sf = part[part.species.eq(sp)]
        roster = A.observed_members(sf).set_index('gene_id')
        ids = out[sp + '_gene_id']
        out[sp + '_direction'] = np.asarray(LABELS)[labels[sp]]
        for direction, sign in [('proximal', -1), ('distal', 1)]:
            support = sf.assign(hit=sf.event.eq(sign)).groupby('gene_id').hit.sum()
            out[sp + '_' + direction + '_support'] = support.reindex(ids).to_numpy(int)
        measured = sf.groupby('gene_id').unit.nunique()
        out[sp + '_measured_studies'] = measured.reindex(ids).to_numpy(int)
        out[sp + '_opposing_study_allowed'] = np.minimum(
            out[sp + '_proximal_support'], out[sp + '_distal_support']).gt(0)
        expected = roster.direction.reindex(ids).fillna('nonrecurrent').to_numpy()
        assert np.array_equal(out[sp + '_direction'], expected)
    return out.assign(threshold=threshold)


def endpoint_replay(frame, measured, features, archives):
    """Validate the two previously reported points, including every null draw."""
    keys = ['contrast', 'threshold', 'model', 'human_direction', 'mouse_direction']
    old = pd.read_csv(SOURCE / 'ortholog_3x3.tsv', sep='\t', float_precision='round_trip')
    new = frame[frame.threshold.isin([.02, .05])]
    aa, bb = old.set_index(keys).sort_index(), new.set_index(keys).sort_index()
    assert aa.index.equals(bb.index)
    for column in old.columns.difference(keys + ['q_joint_16_per_model']):
        if pd.api.types.is_numeric_dtype(aa[column]):
            np.testing.assert_allclose(aa[column], bb[column], atol=0, rtol=0, equal_nan=True)
        else:
            assert aa[column].equals(bb[column]), column
    old_members = pd.read_csv(SOURCE / 'measured_ortholog_members.tsv.gz', sep='\t')
    pair_keys = ['contrast', 'threshold', 'human_gene_id', 'mouse_gene_id']
    aa = old_members.set_index(pair_keys).sort_index()
    bb = measured[measured.threshold.isin([.02, .05])].set_index(pair_keys).sort_index()
    endpoint_pairs = len(bb)
    assert aa.index.equals(bb.index)
    for column in aa:
        assert aa[column].equals(bb[column]), column
    old_features = pd.read_csv(SOURCE / 'ortholog_assay_features.tsv.gz', sep='\t',
                               float_precision='round_trip')
    feature_keys = ['contrast', 'threshold', 'species', 'gene_id']
    aa = old_features.set_index(feature_keys).sort_index()
    bb = features[features.threshold.isin([.02, .05])].set_index(feature_keys).sort_index()
    assert aa.index.equals(bb.index)
    for column in aa:
        if pd.api.types.is_numeric_dtype(aa[column]):
            np.testing.assert_allclose(aa[column], bb[column], atol=0, rtol=0)
        else:
            assert aa[column].equals(bb[column]), column
    counts = 0
    with np.load(SOURCE / 'null_draws.npz') as old_null:
        for key, values in archives.items():
            if key.split('|')[2] in ['0.02', '0.05']:
                np.testing.assert_array_equal(values[:10000], old_null[key])
                counts += 10000
    return dict(endpoint_rows=len(new), endpoint_pairs=endpoint_pairs, exact_draws=counts,
                cells_features_members_p_q8_exact=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--protocol', required=True, type=Path)
    args = parser.parse_args()
    args.protocol = args.protocol.resolve()
    protocol = json.loads(args.protocol.read_text())
    thresholds = protocol['thresholds']
    models = protocol['models']
    assert thresholds == [.01, .015, .02, .025, .03, .035, .04, .045, .05]
    assert models == ['measurement_PAS_capacity', 'full_assay']
    assert protocol['permutations'] == 100000
    paths = [args.protocol, SOURCE / 'protocol.json', SOURCE / 'analyze.py',
        SOURCE / 'ortholog_3x3.tsv', SOURCE / 'measured_ortholog_members.tsv.gz',
        SOURCE / 'ortholog_assay_features.tsv.gz', SOURCE / 'null_draws.npz',
        TABLES / 'recurrence_controls/recurrence_frozen_covariates.tsv.gz',
        TABLES / 'cross_overlap/ortholog_matched_features.tsv.gz',
        ROOT / 'notes/2026-10-10-directional-recurrence/analyze.py',
        ROOT / 'notes/2026-10-02-null-model-controls/recurrence.py',
        ROOT / 'notes/2026-10-02-null-model-controls/common.py',
        ROOT / 'analysis_code/pas_rerun/downstream_runtime.py']
    if 'unified_protocol_path' in protocol:
        paths.append(ROOT / protocol['unified_protocol_path'])
    hashes = {str(p.relative_to(ROOT)): A.digest(p) for p in paths}
    frozen = pd.read_csv(paths[7], sep='\t', float_precision='round_trip')
    assert set(frozen.threshold) == {.02, .05} and frozen.metric.eq('WPP').all()
    base = frozen[frozen.threshold.eq(.02)].sort_values(KEY).reset_index(drop=True)
    strict = frozen[frozen.threshold.eq(.05)].sort_values(KEY).reset_index(drop=True)
    fixed = [c for c in frozen if c not in ['threshold', 'event', 'can_negative',
                                           'can_positive', 'stratum']]
    assert base[fixed].equals(strict[fixed])
    assert not base.duplicated(KEY).any()
    for threshold in [.02, .05]:
        ref = frozen[frozen.threshold.eq(threshold)].sort_values(KEY).reset_index(drop=True)
        current = threshold_events(base, threshold)
        assert current[['event', 'can_negative', 'can_positive']].equals(
            ref[['event', 'can_negative', 'can_positive']].astype(
                {'event': 'int8', 'can_negative': 'bool', 'can_positive': 'bool'}))
    orth = pd.read_csv(paths[8], sep='\t')[['human_gene_id', 'mouse_gene_id',
        'human_gene_symbol', 'mouse_gene_symbol', 'is_high_confidence', 'homology_type']]
    assert len(orth) == 16092 and orth.is_high_confidence.eq(1).all()
    assert orth.homology_type.eq('ortholog_one2one').all()
    assert orth.human_gene_id.is_unique and orth.mouse_gene_id.is_unique
    paired, diagnostics, features, cross, mobility, events, archives = [], [], [], [], [], [], {}
    anchors = []
    expected_background = {'P_M': 751, 'P_S': 1656}
    inventories = {}
    for contrast, fixed_part in base.groupby('contrast', sort=True):
        pairs = A.common_pairs(fixed_part, orth)
        assert len(pairs) == expected_background[contrast]
        inventories[contrast] = pairs.copy()
        # The inventory and quantile bins are fixed before any threshold labels.
        fixed_features = {sp: A.gene_features(fixed_part[fixed_part.species.eq(sp)],
            pairs[sp + '_gene_id'].tolist()) for sp in ['human', 'mouse']}
        for threshold in thresholds:
            part = threshold_events(fixed_part, threshold)
            part = part.sort_values(KEY).reset_index(drop=True)
            events.append(part.drop(columns=['stratum']))
            current_pairs = A.common_pairs(part, orth)
            assert pairs.equals(current_pairs)
            context = dict(contrast=contrast, threshold=threshold, measured_orthologs=len(pairs))
            labels, feats = {}, {}
            for sp in ['human', 'mouse']:
                ids = pairs[sp + '_gene_id'].tolist()
                sf = part[part.species.eq(sp)]
                feats[sp] = A.gene_features(sf, ids)
                constant = [c for c in fixed_features[sp] if c not in
                            ['proximal_capacity', 'distal_capacity']]
                assert fixed_features[sp][constant].equals(feats[sp][constant])
                labels[sp] = A.class_labels(sf, ids)
                features.append(feats[sp].assign(species=sp, contrast=contrast, threshold=threshold))
            paired.append(paired_roster(pairs, part, labels, threshold).assign(contrast=contrast))
            h, m = labels['human'], labels['mouse']
            for model in models:
                cells = {sp: A.partitions(feats[sp], model) for sp in ['human', 'mouse']}
                prob = {sp: A.class_probabilities(labels[sp], cells[sp]) for sp in ['human', 'mouse']}
                analytic = prob['human'].T @ prob['mouse']
                run_seed = A.seed('ortholog_direction', contrast, threshold, model)
                null = A.cross_draws(h, m, cells['human'], cells['mouse'],
                                    np.random.default_rng(run_seed), n=protocol['permutations'])
                archives['|'.join(map(str, ('ortholog', contrast, threshold, model)))] = null
                for sp in ['human', 'mouse']:
                    diag, summary = feature_diagnostics(feats[sp], labels[sp], cells[sp],
                        context | dict(model=model, species=sp))
                    diagnostics.append(diag)
                    mobility.append(summary)
                for a, hd in enumerate(LABELS):
                    for b, md in enumerate(LABELS):
                        values = null[:, a, b]
                        row = context | dict(model=model, human_direction=hd, mouse_direction=md,
                            seed=run_seed, **A.null_row(int(((h == a) & (m == b)).sum()),
                                                      values, analytic[a, b]))
                        row['tested'] = a < 2 and b < 2
                        row.update(tail_precision(row['tail_count'], len(values)))
                        cross.append(row)
                        if threshold in [.02, .05]:
                            anchor = context | dict(model=model, human_direction=hd,
                                mouse_direction=md, seed=run_seed,
                                **A.null_row(int(((h == a) & (m == b)).sum()),
                                             values[:10000], analytic[a, b]), tested=a < 2 and b < 2)
                            anchors.append(anchor)
                print(contrast, threshold, model, 'fixed pairs', len(pairs), flush=True)
    result = apply_correction(pd.DataFrame(cross))
    anchor_result = pd.DataFrame(anchors)
    anchor_result.loc[~anchor_result.tested, 'p_enrichment'] = np.nan
    anchor_result['q_eight_per_threshold_model'] = np.nan
    for _, ii in anchor_result[anchor_result.tested].groupby(['threshold', 'model']).groups.items():
        assert len(ii) == 8
        anchor_result.loc[ii, 'q_eight_per_threshold_model'] = A.bh(anchor_result.loc[ii, 'p_enrichment'])
    members = pd.concat(paired, ignore_index=True)
    ffeatures = pd.concat(features, ignore_index=True)
    replay = endpoint_replay(anchor_result, members, ffeatures, archives)
    summary = pd.DataFrame(mobility)
    save(result, 'ortholog_3x3_threshold_scan.tsv')
    save(result[result.tested], 'directional_tests.tsv')
    save(anchor_result, 'anchor_10000_replay.tsv')
    save(members, 'measured_ortholog_members.tsv.gz')
    save(ffeatures, 'ortholog_assay_features.tsv.gz')
    save(pd.concat(diagnostics, ignore_index=True), 'ortholog_stratum_diagnostics.tsv.gz')
    save(summary, 'mobility_summary.tsv')
    save(pd.concat(events, ignore_index=True), 'study_events.tsv.gz')
    save(pd.concat([p.assign(contrast=c) for c, p in inventories.items()], ignore_index=True),
         'fixed_ortholog_background.tsv')
    np.savez_compressed(HERE / 'null_draws.npz', **archives)
    assert all(A.digest(ROOT / p) == value for p, value in hashes.items()), 'Input changed'
    (HERE / 'analysis_validation.json').write_text(json.dumps(dict(status='PASS',
        fixed_backgrounds=expected_background, threshold_count=len(thresholds),
        contexts=len(archives), draws_per_context=protocol['permutations'],
        tests=len(result[result.tested]), scan_family_per_model=72, scan_family_all_models=144,
        study_rows=len(base), input_fixed_fields_identical=True,
        labels_recomputed_from_all_series_extrema=True, bins_constant_across_thresholds=True,
        missing_observations_excluded=True, class_margins_exact=True,
        analytic_means_within_six_MC_SE=True, endpoint_replay=replay), indent=2) + '\n')
    output_names = ['ortholog_3x3_threshold_scan.tsv', 'directional_tests.tsv',
        'measured_ortholog_members.tsv.gz', 'ortholog_assay_features.tsv.gz',
        'ortholog_stratum_diagnostics.tsv.gz', 'mobility_summary.tsv', 'study_events.tsv.gz',
        'fixed_ortholog_background.tsv', 'null_draws.npz', 'analysis_validation.json']
    output_names.append('anchor_10000_replay.tsv')
    (HERE / 'analysis_manifest.json').write_text(json.dumps(dict(status='PASS',
        input_sha256=hashes, output_sha256={n: A.digest(HERE / n) for n in output_names},
        script_sha256=A.digest(Path(__file__)), python_executable=sys.executable,
        python=sys.version, numpy=np.__version__, pandas=pd.__version__,
        inference_scope='Conditional gene correspondence; not independent species replication.'),
        indent=2) + '\n')


if __name__ == '__main__':
    main()
