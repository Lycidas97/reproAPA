"""Independent class/label/conditional-null QA; no production analysis imports."""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = HERE.parent
OLD = ROOT / 'notes/2026-10-10-ortholog-threshold-scan'
THRESHOLDS = np.array([.01, .015, .02, .025, .03, .035, .04, .045, .05])
LABELS = ['proximal', 'distal', 'nonrecurrent']
MODELS = ['measurement_PAS_capacity', 'full_assay']
PAIRKEY = ['human_gene_id', 'mouse_gene_id']
checks = []


def verify(condition, reason):
    if not bool(condition):
        raise AssertionError(reason)
    checks.append(reason)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def canonical(values):
    a = np.asarray(values, float)
    distances = np.abs(a[:, None] - THRESHOLDS[None, :])
    verify((distances.min(1) < 1e-12).all(), 'Parameter thresholds map to the frozen nine values')
    return THRESHOLDS[distances.argmin(1)]


def boolean(values):
    text = values.astype(str).str.lower()
    verify(text.isin(['true', 'false']).all(), 'Boolean inputs contain no missing or arbitrary truthy strings')
    return text.eq('true')


def bh(p):
    p = np.asarray(p, float)
    ii = np.argsort(p, kind='stable')
    adjusted = p[ii] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty_like(p)
    out[ii] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return out


def save(df, name):
    df.to_csv(HERE / name, sep='\t', index=False, float_format='%.17g')


def species_group(features, model):
    columns = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
    if model == 'full_assay':
        columns += ['mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']
    # Independent Python tuple collection rather than production pandas grouping.
    index = {}
    for i, values in enumerate(features[columns].itertuples(index=False, name=None)):
        index.setdefault(values, []).append(i)
    return [np.asarray(ii, int) for ii in index.values()]


def conditional_prob(labels, groups):
    output = np.empty((len(labels), 3), float)
    for ii in groups:
        output[ii] = [(labels[ii] == j).mean() for j in range(3)]
    return output


def independent_primary_mc(primary_states, analytic):
    """Shuffle via independently ranked uniforms, not numpy.random.permuted."""
    records = []
    for context_index, ((label, model), (labs, groups)) in enumerate(primary_states.items()):
        rng = np.random.default_rng(2026101100 + context_index)
        counts = np.empty(100000, np.int32)
        n = len(labs['human'])
        for start in range(0, len(counts), 500):
            size = min(500, len(counts) - start)
            permuted = {}
            for species in ['human', 'mouse']:
                permuted[species] = np.tile(labs[species], (size, 1))
                for ii in groups[species]:
                    if np.unique(labs[species][ii]).size <= 1:
                        continue
                    random_rank = rng.random((size, len(ii))).argsort(axis=1)
                    permuted[species][:, ii] = labs[species][ii][random_rank]
            counts[start:start + size] = ((permuted['human'] == 0) & (permuted['mouse'] == 0)).sum(axis=1)
        selected = analytic[(analytic.threshold == .02) & analytic.pair_class.eq(label) & analytic.model.eq(model)
                            & analytic.human_direction.eq('proximal') & analytic.mouse_direction.eq('proximal')]
        verify(len(selected) == 1, f'{label}/{model} unique primary analytic row')
        row = selected.iloc[0]
        tail = int((counts >= row.observed).sum())
        rate = tail / len(counts)
        sd = float(counts.std(ddof=1))
        verify(abs(counts.mean() - row.analytic_expected) <= 6 * sd / math.sqrt(len(counts)) + 1e-10,
               f'{label}/{model} independent null mean agrees with analytic conditional expectation')
        records.append(dict(pair_class=label, model=model, threshold=.02, observed=int(row.observed),
            analytic_expected=float(row.analytic_expected), permutations=len(counts), independent_seed=2026101100 + context_index,
            independent_tail=tail, independent_p=(tail + 1) / (len(counts) + 1),
            independent_MC_SE=math.sqrt(rate * (1 - rate) / len(counts)),
            independent_null_mean=float(counts.mean()), independent_null_sd=sd,
            independent_support_n=len(np.unique(counts))))
    result = pd.DataFrame(records)
    save(result, 'independent_primary_mc.tsv')
    return result


def main():
    protocol_path = BASE / 'PROTOCOL.zh.md'
    freeze = json.loads((BASE / 'protocol_freeze.json').read_text())
    verify(digest(protocol_path) == freeze['protocol_sha256'], 'Protocol content matches frozen SHA256')
    paths = {name: ROOT / name for name in freeze['input_sha256']}
    initial_hashes = {name: digest(path) for name, path in paths.items()}
    verify(initial_hashes == freeze['input_sha256'], 'All frozen inputs match protocol SHA256 before QA')
    events = read(ROOT / 'notes/2026-10-10-shared-proximal-mechanisms/event_inventory.tsv.gz')
    members = read(OLD / 'measured_ortholog_members.tsv.gz')
    members = members[members.contrast.eq('P_M')].copy()
    members['threshold'] = canonical(members.threshold)
    features = read(OLD / 'ortholog_assay_features.tsv.gz')
    features = features[features.contrast.eq('P_M')].copy()
    features['threshold'] = canonical(features.threshold)
    studies = read(OLD / 'study_events.tsv.gz')
    studies = studies[studies.contrast.eq('P_M')].copy()
    studies['threshold'] = canonical(studies.threshold)
    background = read(OLD / 'fixed_ortholog_background.tsv')
    background = background[background.contrast.eq('P_M')].copy().sort_values(PAIRKEY).reset_index(drop=True)
    verify(len(background) == 751 and not background.duplicated(PAIRKEY).any(), 'Background contains all 751 unique measured ortholog pairs')
    verify(len(events) == 3467 and not events.duplicated(['species', 'study_id', 'route', 'gene_id']).any(), 'All 3467 old route events are unique and retained')
    verify(not events.driver_class.isna().any(), 'Unknown class is explicit, never missing or silently imputed')
    events['full_axis'] = boolean(events.full_primary_axis_tandem)
    original = read(BASE / 'classification/event_class_pairs.tsv').set_index(PAIRKEY).sort_index()
    rebuilt = background.copy()
    species_rows = []
    for species in ['human', 'mouse']:
        ids = set(background[species + '_gene_id'])
        verify(set(events.loc[events.species.eq(species), 'gene_id']) == ids, f'{species} event universe is exactly the 751 background genes')
        current = studies[studies.species.eq(species) & studies.threshold.eq(.02) & studies.gene_id.isin(ids)]
        expected_routes = set()
        for r in current.itertuples():
            for route in r.series_pattern.split(';'):
                expected_routes.add((r.unit, route, r.gene_id))
        observed_routes = set(events.loc[events.species.eq(species), ['study_id', 'route', 'gene_id']].itertuples(index=False, name=None))
        verify(expected_routes == observed_routes, f'{species} classification uses exactly frozen contributing study/series routes')
        classes, full_axes = {}, {}
        for gene, part in events[events.species.eq(species)].groupby('gene_id'):
            unique = set(part.driver_class)
            if unique == {'strict_tandem_UTR'}:
                label = 'tandem'
            elif unique == {'alternative_terminal_exon_compatible'}:
                label = 'ALE'
            else:
                label = next(iter(unique)) if len(unique) == 1 else 'mixed'
            classes[gene] = label
            full_axes[gene] = bool(part.full_axis.all())
            species_rows.append(dict(species=species, gene_id=gene, event_class=label,
                n_routes=len(part), n_studies=part.study_id.nunique(), full_axis=full_axes[gene]))
        rebuilt[species + '_event_class'] = rebuilt[species + '_gene_id'].map(classes)
        rebuilt[species + '_full_axis'] = rebuilt[species + '_gene_id'].map(full_axes)
        for threshold in THRESHOLDS:
            s = studies[studies.species.eq(species) & studies.threshold.eq(threshold) & studies.gene_id.isin(ids)]
            # Re-evaluate exact effect boundary; parameter key canonicalization does not modify effects.
            s = s.assign(expected_event=np.where(s.effect_max.le(-threshold), -1,
                                      np.where(s.effect_min.ge(threshold), 1, 0)))
            verify(s.expected_event.eq(s.event).all(), f'{species} study extrema reproduce frozen {threshold:g} event labels')
            support = s.groupby('gene_id').expected_event.apply(lambda v: ((v == -1).sum(), (v == 1).sum()))
            calculated = support.map(lambda x: 'proximal' if x[0] >= 2 else ('distal' if x[1] >= 2 else 'nonrecurrent'))
            old = members[members.threshold.eq(threshold)].set_index(species + '_gene_id')
            verify(calculated.reindex(old.index).eq(old[species + '_direction']).all(), f'{species} >=2-study support reproduces every {threshold:g} recurrence label')
            for direction, index in [('proximal', 0), ('distal', 1)]:
                values = support.map(lambda x: x[index]).reindex(old.index).to_numpy(int)
                verify(np.array_equal(values, old[species + '_' + direction + '_support']), f'{species} all measured study {direction} support counts reproduce {threshold:g}')
    rebuilt['pair_class'] = 'other/mixed'
    rebuilt.loc[rebuilt.human_event_class.eq('tandem') & rebuilt.mouse_event_class.eq('tandem'), 'pair_class'] = 'both_tandem'
    rebuilt.loc[rebuilt.human_event_class.eq('ALE') & rebuilt.mouse_event_class.eq('ALE'), 'pair_class'] = 'both_ALE'
    rebuilt['strict_full_axis_sensitivity'] = rebuilt.pair_class.eq('both_tandem') & rebuilt.human_full_axis & rebuilt.mouse_full_axis
    idx = rebuilt.set_index(PAIRKEY).sort_index()
    verify(idx.index.equals(original.index), 'Production class inventory retains all background identities')
    for column in ['human_event_class', 'mouse_event_class', 'pair_class', 'strict_full_axis_sensitivity']:
        verify(idx[column].eq(original[column]).all(), f'Independent aggregation matches production {column}')
    save(rebuilt, 'independent_class_inventory.tsv')
    save(pd.DataFrame(species_rows), 'independent_species_inventory.tsv')
    class_counts = rebuilt.groupby('pair_class').size().to_dict()
    compositions = []
    for threshold in THRESHOLDS:
        subset = members[members.threshold.eq(threshold)].merge(rebuilt, on=PAIRKEY, validate='one_to_one', suffixes=('', '_class'))
        core = subset[subset.human_direction.eq('proximal') & subset.mouse_direction.eq('proximal')]
        for label in ['both_tandem', 'both_ALE', 'other/mixed']:
            compositions.append(dict(threshold=threshold, pair_class=label, background_n=class_counts[label], core_n=int(core.pair_class.eq(label).sum())))
        if threshold in [.02, .05]:
            verify(len(core) == (27 if threshold == .02 else 7), f'Frozen total shared core remains {len(core)} at {threshold:g}')
    save(pd.DataFrame(compositions), 'independent_core_counts.tsv')
    analytic, primary_states, mobility = [], {}, []
    for threshold in THRESHOLDS:
        old = members[members.threshold.eq(threshold)].set_index(PAIRKEY)
        feature = features[features.threshold.eq(threshold)]
        for label in ['both_tandem', 'both_ALE', 'strict_full_axis_tandem']:
            choose = rebuilt.strict_full_axis_sensitivity if label == 'strict_full_axis_tandem' else rebuilt.pair_class.eq(label)
            pairs = rebuilt[choose].copy().reset_index(drop=True)
            pair_index = pd.MultiIndex.from_frame(pairs[PAIRKEY])
            directions = old.reindex(pair_index)
            for model in MODELS:
                labs, groups, probabilities = {}, {}, {}
                for species in ['human', 'mouse']:
                    labs[species] = directions[species + '_direction'].map({name: i for i, name in enumerate(LABELS)}).to_numpy(int)
                    sf = feature[feature.species.eq(species)].set_index('gene_id').reindex(pairs[species + '_gene_id']).reset_index()
                    verify(not sf.isna().any().any(), f'{label}/{threshold:g}/{model}/{species} features retain complete frozen qualification')
                    groups[species] = species_group(sf, model)
                    probabilities[species] = conditional_prob(labs[species], groups[species])
                    sizes = [len(ii) for ii in groups[species]]
                    homogeneous = [ii for ii in groups[species] if np.unique(labs[species][ii]).size == 1]
                    changed = sum(len(ii) * (1 - np.square([(labs[species][ii] == j).mean() for j in range(3)]).sum())
                                  for ii in groups[species]) / len(pairs)
                    mobility.append(dict(pair_class=label, threshold=threshold, model=model, species=species,
                        measured_genes=len(pairs), n_strata=len(sizes), n_singleton_strata=sum(n == 1 for n in sizes),
                        largest_stratum=max(sizes), singleton_gene_fraction=sum(n for n in sizes if n == 1) / len(pairs),
                        frozen_label_fraction=sum(len(ii) for ii in homogeneous) / len(pairs), expected_changed_fraction=changed))
                expected = probabilities['human'].T @ probabilities['mouse']
                for hi, human in enumerate(LABELS):
                    for mi, mouse in enumerate(LABELS):
                        analytic.append(dict(threshold=threshold, pair_class=label, model=model,
                            human_direction=human, mouse_direction=mouse, background_n=len(pairs),
                            observed=int(((labs['human'] == hi) & (labs['mouse'] == mi)).sum()),
                            analytic_expected=float(expected[hi, mi])))
                if threshold == .02 and label != 'strict_full_axis_tandem':
                    primary_states[(label, model)] = (labs, groups)
    analytic = pd.DataFrame(analytic)
    save(analytic, 'independent_analytic_expectations.tsv')
    mobility = pd.DataFrame(mobility)
    save(mobility, 'independent_mobility.tsv')
    independent_mc = independent_primary_mc(primary_states, analytic)
    output = dict(status='CLASS_LABEL_EXPECTATION_STAGE_PASS', checks_passed=len(checks), checks_skipped=0,
        pair_class_counts=class_counts, strict_full_axis_n=int(rebuilt.strict_full_axis_sensitivity.sum()),
        frozen_input_sha256=initial_hashes, production_protocol_sha256=digest(protocol_path), python=sys.executable,
        production_inference_checked=False)
    (HERE / 'preflight_validation.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2), flush=True)
    # Optional full validation is enabled only when the production inference table exists.
    table = BASE / 'overlap/class_direction_3x3.tsv'
    if not table.exists():
        verify({name: digest(path) for name, path in paths.items()} == initial_hashes, 'Frozen inputs unchanged after preflight QA')
        return
    production = read(table)
    production['threshold'] = canonical(production.threshold)
    check_production(production, analytic, independent_mc, mobility, paths, initial_hashes, output)


def check_production(production, analytic, independent_mc, mobility, paths, initial_hashes, output):
    keys = ['pair_class', 'threshold', 'model', 'human_direction', 'mouse_direction']
    independent = analytic[analytic.pair_class.ne('strict_full_axis_tandem')].set_index(keys).sort_index()
    tested = production.set_index(keys).sort_index()
    verify(independent.index.equals(tested.index), 'All 324 production directional and nonrecurrent cells are present exactly once')
    verify(np.array_equal(independent.observed, tested.observed), 'All 324 observed cell counts independently reproduce')
    delta = np.abs(independent.analytic_expected.to_numpy() - tested.analytic_expected.to_numpy())
    verify(delta.max() < 1e-10, 'All 324 conditional analytic expectations independently reproduce')
    verify(np.array_equal(independent.background_n, tested.measured_orthologs), 'Every table denominator uses its class-specific background')
    directional = production.human_direction.isin(LABELS[:2]) & production.mouse_direction.isin(LABELS[:2])
    verify(production.loc[~directional, 'p_enrichment'].isna().all(), 'Nonrecurrent cells have no inferential p-values')
    verify(production.loc[directional, 'p_enrichment'].notna().all(), 'Every planned directional hypothesis retains a p-value')
    for row in production.itertuples():
        verify(abs(row.excess - (row.observed - row.analytic_expected)) < 1e-10, 'Excess is observed minus exact conditional expectation')
        expected_fold = row.observed / row.analytic_expected if row.analytic_expected else np.nan
        verify(np.isclose(expected_fold, row.fold, atol=1e-12, rtol=1e-12, equal_nan=True), 'Fold uses class-conditional expectation, not full 751 expectation')
    # Check archived production samples independently; no production function calls.
    archives = np.load(BASE / 'overlap/null_draws.npz')
    verify(len(archives.files) == 36, 'All 36 class/model/threshold null contexts are archived')
    for key in archives.files:
        label, threshold_text, model = key.split('|')
        threshold = float(threshold_text)
        frame = production[production.pair_class.eq(label) & production.threshold.eq(threshold) & production.model.eq(model)]
        observed_matrix = frame.pivot(index='human_direction', columns='mouse_direction', values='observed').reindex(index=LABELS, columns=LABELS).to_numpy()
        draws = archives[key]
        verify(draws.shape == (100000, 3, 3), f'{key} has all planned samples and nine cells')
        verify((draws.sum((1, 2)) == observed_matrix.sum()).all(), f'{key} null preserves class-specific universe')
        verify((draws.sum(2) == observed_matrix.sum(1)).all(), f'{key} null preserves human direction margins')
        verify((draws.sum(1) == observed_matrix.sum(0)).all(), f'{key} null preserves mouse direction margins')
        for row in frame.itertuples():
            hi, mi = LABELS.index(row.human_direction), LABELS.index(row.mouse_direction)
            d = draws[:, hi, mi]
            tail = int((d >= row.observed).sum())
            verify(abs(d.mean() - row.null_mean) < 1e-12 and abs(d.std(ddof=1) - row.null_sd) < 1e-12, 'Archived samples reproduce reported null mean and SD')
            verify(len(np.unique(d)) == row.null_support_points, 'Archived samples reproduce null support/degeneracy')
            verify(bool(d.std(ddof=1) > 0) == bool(row.inferentially_mobile), 'Nonmobile states are explicitly retained')
            if row.directional_test:
                verify(tail == row.tail_count, 'Every reported directional upper-tail count independently matches archive')
                verify(abs((tail + 1) / (len(d) + 1) - row.p_enrichment) < 1e-15, 'Every directional p uses inclusive upper tail and plus-one correction')
            if row.null_sd == 0 and row.directional_test:
                verify(row.p_enrichment == 1, 'Frozen null keeps planned family position at p=1')
    archives.close()
    primary = production.human_direction.eq('proximal') & production.mouse_direction.eq('proximal')
    for column, select, n in [
        ('q_primary4', primary & production.threshold.eq(.02), 4),
        ('q_sensitivity4', primary & production.threshold.eq(.05), 4),
        ('q_scan36', primary, 36),
        ('q_secondary_scan144', directional, 144),
    ]:
        verify(int(select.sum()) == n, f'{column} includes exactly all {n} planned hypotheses')
        recalculated = bh(production.loc[select, 'p_enrichment'])
        verify(np.allclose(recalculated, production.loc[select, column], atol=1e-15, rtol=1e-12), f'{column} independent BH computation matches')
        verify(production.loc[~select, column].isna().all(), f'{column} remains missing outside its specified family')
    for threshold in THRESHOLDS:
        select = directional & production.threshold.eq(threshold)
        verify(int(select.sum()) == 16, f'Secondary family contains all 16 hypotheses at {threshold:g}')
        verify(np.allclose(bh(production.loc[select, 'p_enrichment']), production.loc[select, 'q_secondary16_per_threshold'], atol=1e-15, rtol=1e-12),
               f'Secondary BH16 independently matches at {threshold:g}')
    comparison = independent_mc.merge(production[primary & production.threshold.eq(.02)][['pair_class', 'model', 'p_enrichment', 'tail_count', 'permutations']],
                                      on=['pair_class', 'model'], validate='one_to_one', suffixes=('', '_production'))
    for row in comparison.itertuples():
        prod_rate = row.tail_count / row.permutations_production
        independent_rate = row.independent_tail / row.permutations
        variance = prod_rate * (1 - prod_rate) / row.permutations_production + independent_rate * (1 - independent_rate) / row.permutations
        allowance = 6 * math.sqrt(variance) + 2 / min(row.permutations, row.permutations_production)
        verify(abs(row.p_enrichment - row.independent_p) <= allowance,
               f'{row.pair_class}/{row.model} independently generated primary tail agrees within six combined MC SE plus discretization')
    save(comparison, 'independent_primary_mc_comparison.tsv')
    strict = read(BASE / 'overlap/full_axis_tandem_descriptive.tsv')
    strict['threshold'] = canonical(strict.threshold)
    strict_keys = ['threshold', 'model', 'human_direction', 'mouse_direction']
    aa = analytic[analytic.pair_class.eq('strict_full_axis_tandem')].set_index(strict_keys).sort_index()
    bb = strict.set_index(strict_keys).sort_index()
    verify(aa.index.equals(bb.index) and len(aa) == 162, 'Strict full-axis table retains all descriptive cells')
    verify(np.array_equal(aa.observed, bb.observed) and np.allclose(aa.analytic_expected, bb.analytic_expected, atol=1e-10, rtol=0),
           'Strict full-axis counts and class-conditional expectations independently match')
    verify(not any('p_enrichment' == c or c.startswith('q_') for c in strict.columns), 'Strict full-axis sensitivity remains descriptive as frozen')
    mobility_prod = read(BASE / 'overlap/mobility_summary.tsv')
    strict_mobility = read(BASE / 'overlap/full_axis_tandem_mobility.tsv')
    strict_mobility['pair_class'] = 'strict_full_axis_tandem'
    mobility_prod = pd.concat([mobility_prod, strict_mobility], ignore_index=True)
    mobility_prod['threshold'] = canonical(mobility_prod.threshold)
    mobility_keys = ['pair_class', 'threshold', 'model', 'species']
    aa = mobility.set_index(mobility_keys).sort_index()
    bb = mobility_prod.set_index(mobility_keys).sort_index()
    verify(aa.index.equals(bb.index), 'Mobility diagnostics include every planned class/model/species/threshold context')
    for c in aa.columns:
        verify(np.allclose(aa[c], bb[c], atol=1e-12, rtol=1e-12), f'Independent mobility matches {c}')
    manifest = json.loads((BASE / 'overlap/analysis_manifest.json').read_text())
    for name, hash_value in manifest['input_sha256'].items():
        verify(digest(ROOT / name) == hash_value, f'Production manifest input SHA verifies: {name}')
    for name, hash_value in manifest['output_sha256'].items():
        verify(digest(BASE / 'overlap' / name) == hash_value, f'Production manifest output SHA verifies: {name}')
    verify({name: digest(path) for name, path in paths.items()} == initial_hashes, 'Frozen inputs unchanged after independent QA')
    primary_rows = production[primary & production.threshold.eq(.02)]
    result = output | dict(status='PASS', checks_passed=len(checks), checks_failed=0, checks_skipped=0,
        production_inference_checked=True, all_324_counts_exact=True,
        max_analytic_expected_absolute_error=float(delta.max()), primary_mc_contexts=4,
        independent_mc_draws_per_context=100000,
        production_main_results=primary_rows[['pair_class', 'model', 'observed', 'analytic_expected', 'p_enrichment', 'q_primary4', 'q_scan36']].to_dict(orient='records'),
        independent_script_sha256=digest(Path(__file__)))
    (HERE / 'independent_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    (HERE / 'checks.txt').write_text('\n'.join(checks) + '\n')
    print(json.dumps(dict(status='PASS', checks_passed=len(checks), checks_failed=0, checks_skipped=0,
                         max_expected_error=float(delta.max())), indent=2), flush=True)


if __name__ == '__main__':
    main()
