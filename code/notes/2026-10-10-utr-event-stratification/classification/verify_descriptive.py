"""Independent CSV/standard-library reconstruction of the descriptive contract."""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, median

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
checks = 0


def check(value, message):
    global checks
    checks += 1
    if not value:
        raise AssertionError(message)


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def boolean(value):
    check(value in ['True', 'False'], 'Missing/invalid boolean must not become false')
    return value == 'True'


def classify(values):
    unique = set(values)
    if unique == {'strict_tandem_UTR'}:
        return 'tandem'
    if unique == {'alternative_terminal_exon_compatible'}:
        return 'ALE'
    if len(unique) == 1:
        return unique.pop()
    return 'mixed'


def near(a, b):
    return math.isclose(float(a), float(b), abs_tol=1e-12, rel_tol=0)


def main():
    # Tests encode conservative event interpretation rather than merely replaying output.
    check(classify(['strict_tandem_UTR', 'alternative_terminal_exon_compatible']) == 'mixed',
          'One incompatible route must prevent a homogeneous species class')
    check(classify(['unknown']) == 'unknown', 'Unknown must not be counted as a tandem or ALE event')
    check(classify(['terminal_exon_nonUTR']) == 'terminal_exon_nonUTR',
          'Terminal exon location alone is insufficient to call a UTR event')
    check(classify(['strict_tandem_UTR'] * 2) == classify(['strict_tandem_UTR']),
          'Extra coverage must not change a homogeneous structure label')
    event_path = ROOT / 'notes/2026-10-10-shared-proximal-mechanisms/event_inventory.tsv.gz'
    member_path = ROOT / 'notes/2026-10-10-ortholog-threshold-scan/measured_ortholog_members.tsv.gz'
    frozen_path = ROOT / 'notes/2026-10-10-ortholog-threshold-scan/study_events.tsv.gz'
    events = read(event_path)
    pairs = read(HERE / 'event_class_pairs.tsv')
    members = [row for row in read(member_path) if row['contrast'] == 'P_M']
    check(len(events) == 3467 and len(pairs) == 751 and len(members) == 6759, 'Fixed inventory sizes')
    by_gene = defaultdict(list)
    by_study = defaultdict(list)
    for row in events:
        by_gene[(row['species'], row['gene_id'])].append(row)
        by_study[(row['species'], row['study_id'], row['gene_id'])].append(row)
    pair_lookup = {(row['human_gene_id'], row['mouse_gene_id']): row for row in pairs}
    gene_pairs = {}
    for row in pairs:
        for species in ['human', 'mouse']:
            gene_pairs[(species, row[species + '_gene_id'])] = row
            group = by_gene[(species, row[species + '_gene_id'])]
            label = classify([event['driver_class'] for event in group])
            check(label == row[species + '_event_class'], 'Classify every route, including opposing routes')
            full = all(event['full_primary_axis_tandem'] == 'True' for event in group)
            check(full == boolean(row[species + '_all_full_axis_tandem']), 'Full-axis gate includes every route')
            check(len(group) == int(row[species + '_n_routes']), 'Route coverage')
            per_study = defaultdict(list)
            for event in group:
                per_study[event['study_id']].append(float(event['delta_DWUI']))
            expected = mean([mean(values) for values in per_study.values()])
            check(near(expected, row[species + '_delta_DWUI_equal_study_mean']), 'Study-equal gene effect')
            check(near(median([float(event['delta_DWUI']) for event in group]),
                       row[species + '_delta_DWUI_route_median']), 'Route-median auxiliary effect')
        classes = (row['human_event_class'], row['mouse_event_class'])
        expected = 'both_tandem' if classes == ('tandem', 'tandem') else (
                   'both_ALE' if classes == ('ALE', 'ALE') else 'other/mixed')
        check(row['pair_class'] == expected, 'Both species must satisfy the homogeneous pair class')
        gate = classes == ('tandem', 'tandem') and all(row[sp + '_all_full_axis_tandem'] == 'True'
                                                   for sp in ['human', 'mouse'])
        check(gate == boolean(row['strict_full_axis_sensitivity']), 'Strict sensitivity is an intersection')
    frozen = [row for row in read(frozen_path) if row['contrast'] == 'P_M' and row['metric'] == 'WPP'
              and (row['species'], row['gene_id']) in gene_pairs]
    frozen_base = [row for row in frozen if float(row['threshold']) == .02]
    check(len(frozen_base) == len(by_study), 'No source study was silently omitted')
    for row in frozen_base:
        group = by_study[(row['species'], row['unit'], row['gene_id'])]
        check(set(row['series_pattern'].split(';')) == {event['route'] for event in group},
              'Frozen eligible series must exactly match driver contributors')
        check(len(group) == int(row['n_routes']), 'Frozen route number')
        check(near(mean([float(event['delta_DWUI']) for event in group]), row['effect_mean']),
              'Frozen series mean')
    votes = defaultdict(Counter)
    for row in frozen:
        threshold = float(row['threshold'])
        event = -1 if float(row['effect_max']) <= -threshold else (
                 1 if float(row['effect_min']) >= threshold else 0)
        check(event == int(row['event']), 'All series extrema must support the threshold direction')
        votes[(threshold, row['species'], row['gene_id'])][event] += 1
    marginal = Counter()
    cross = Counter()
    paired_marginal = Counter()
    core_composition = Counter()
    for row in members:
        threshold = float(row['threshold'])
        pair = pair_lookup[(row['human_gene_id'], row['mouse_gene_id'])]
        for species in ['human', 'mouse']:
            support = votes[(threshold, species, row[species + '_gene_id'])]
            check(not (support[-1] >= 2 and support[1] >= 2), 'Opposite recurrence cannot silently share a label')
            label = 'proximal' if support[-1] >= 2 else ('distal' if support[1] >= 2 else 'nonrecurrent')
            check(label == row[species + '_direction'], 'Gene direction requires two supporting studies')
            for name, sign in [('proximal', -1), ('distal', 1)]:
                check(support[sign] == int(row[species + '_' + name + '_support']), 'Frozen study support count')
        scopes = ['all'] + (['both_full_axis_tandem'] if pair['strict_full_axis_sensitivity'] == 'True' else [])
        for scope in scopes:
            cross[(scope, threshold, pair['pair_class'], row['human_direction'], row['mouse_direction'])] += 1
            for species in ['human', 'mouse']:
                marginal[(scope, threshold, species, pair[species + '_event_class'], row[species + '_direction'])] += 1
                paired_marginal[(scope, threshold, species, pair['pair_class'], row[species + '_direction'])] += 1
        if row['human_direction'] == 'proximal' and row['mouse_direction'] == 'proximal':
            core_composition[(threshold, pair['pair_class'])] += 1
    for row in read(HERE / 'recurrence_direction_marginals.tsv'):
        key = (row['scope'], float(row['threshold']), row['species'], row['event_class'], row['direction'])
        check(marginal[key] == int(row['n_genes']), 'All threshold/species/class marginals')
        check(near(int(row['n_genes']) / int(row['class_background']), row['proportion']), 'Marginal denominator')
    for row in read(HERE / 'paired_stratum_direction_marginals.tsv'):
        key = (row['scope'], float(row['threshold']), row['species'], row['pair_class'], row['direction'])
        check(paired_marginal[key] == int(row['n_genes']), 'Pair-class marginals')
    for row in read(HERE / 'recurrence_cross_direction.tsv'):
        key = (row['scope'], float(row['threshold']), row['pair_class'], row['human_direction'], row['mouse_direction'])
        check(cross[key] == int(row['n_pairs']), 'Every 3 by 3 direction cell')
        check(near(int(row['n_pairs']) / int(row['class_background']), row['proportion']), 'Cross-cell denominator')
    check(core_composition[(.02, 'both_tandem')] == 17 and core_composition[(.02, 'both_ALE')] == 4
          and core_composition[(.02, 'other/mixed')] == 6, 'Recover the complete 27 core without selection')
    check(core_composition[(.05, 'both_tandem')] == 6 and core_composition[(.05, 'both_ALE')] == 1
          and core_composition[(.05, 'other/mixed')] == 0, 'Recover the strict seven without selection')
    sign_lookup = {(row['species'], row['study_id'], row['route'], row['gene_id']): row
                   for row in read(HERE / 'route_sign_audit.tsv')}
    for row in events:
        gain, loss = int(row['gain_cut']), int(row['loss_cut'])
        is_proximal = gain < loss if row['strand'] == '+' else gain > loss
        check(is_proximal == (row['gain_is_proximal'] == 'True'), 'Native strand determines positional direction')
        value = float(row['delta_DWUI'])
        net = 0 if abs(value) <= 1e-12 else (1 if value > 0 else -1)
        expected = -1 if is_proximal else 1
        label = 'no_valid_gain_loss' if not (float(row['gain_delta']) > 0 and float(row['loss_delta']) < 0) else (
                'zero_DWUI' if net == 0 else ('concordant' if net == expected else 'discordant'))
        key = (row['species'], row['study_id'], row['route'], row['gene_id'])
        check(sign_lookup[key]['sign_audit'] == label, 'Largest driver pair cannot redefine whole-axis direction')
    # Extra ALE qualification is semantic coverage, not a redefinition of the primary classes.
    utr_path = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/tandem/annotation_terminal_utrs.tsv.gz'
    interval_lookup = defaultdict(list)
    for row in read(utr_path):
        interval_lookup[(row['species'], row['gene_id'], row['chrom'].removeprefix('chr'), row['strand'])].append(
            (int(row['utr_start']) - 1, int(row['utr_end']), row['transcript_id']))
    route_coverage = read(HERE / 'ALE_driver_terminal_UTR_coverage_routes.tsv')
    covered_genes = defaultdict(list)
    for row in route_coverage:
        key = (row['species'], row['study_id'], row['route'], row['gene_id'])
        source = next(event for event in by_gene[(row['species'], row['gene_id'])]
                      if event['study_id'] == row['study_id'] and event['route'] == row['route'])
        intervals = interval_lookup[(source['species'], source['gene_id'],
                                     source['chrom'].removeprefix('chr'), source['strand'])]
        endpoints = []
        for end in ['gain', 'loss']:
            cut = int(source[end + '_cut'])
            base = cut - 1 if source['strand'] == '+' else cut
            fitted = {tid for lo, hi, tid in intervals if lo <= base < hi}
            check(bool(fitted) == boolean(row[end + '_terminal_UTR_covered']), 'ALE endpoint UTR coverage')
            check(';'.join(sorted(fitted)) == row[end + '_UTR_transcripts'], 'All compatible endpoint UTR IDs')
            endpoints.append(bool(fitted))
        covered_genes[(source['species'], source['gene_id'])].append(all(endpoints))
    covered_pairs = read(HERE / 'both_ALE_terminal_UTR_coverage_pairs.tsv')
    check(len(covered_pairs) == 52, 'Extra qualification retains every ALE pair')
    n_qualified = 0
    for row in covered_pairs:
        flags = []
        for species in ['human', 'mouse']:
            expected = all(covered_genes[(species, row[species + '_gene_id'])])
            check(expected == boolean(row[species + '_all_routes_both_ends_UTR_covered']),
                  'Extra UTR qualification includes every route endpoint')
            flags.append(expected)
        qualified = all(flags)
        check(qualified == boolean(row['both_species_all_routes_both_ends_UTR_covered']), 'Both species qualified')
        n_qualified += qualified
    check(n_qualified == 15, 'UTR-covered ALE subset denominator')
    manifest = json.loads((HERE / 'classification_contract.json').read_text())
    for name, expected in manifest['input_sha256'].items():
        check(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, 'Inputs remained unchanged')
    check(hashlib.sha256((HERE / 'event_class_pairs.tsv').read_bytes()).hexdigest() ==
          '5ba39cd3a5439557c7c3edf604ec6c9ca03c88acb0249996f4748572c17842c3', 'Frozen pair-class table retained its SHA')
    result = dict(passed=True, checks=checks, skipped=0, source_events=3467,
                  ortholog_pairs=751, thresholds=9, biological_p_values_computed=False)
    (HERE / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
