"""Describe frozen P–M response directions after terminal-structure stratification."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EVENTS = ROOT / 'notes/2026-10-10-shared-proximal-mechanisms/event_inventory.tsv.gz'
MEMBERS = ROOT / 'notes/2026-10-10-ortholog-threshold-scan/measured_ortholog_members.tsv.gz'
BACKGROUND = ROOT / 'notes/2026-10-10-ortholog-threshold-scan/fixed_ortholog_background.tsv'
STUDY_EVENTS = ROOT / 'notes/2026-10-10-ortholog-threshold-scan/study_events.tsv.gz'
UTRS = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/tandem/annotation_terminal_utrs.tsv.gz'
KEYS = ['human_gene_id', 'mouse_gene_id']
DIRECTIONS = ['proximal', 'distal', 'nonrecurrent']


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g')


def strict_bool(values):
    mapping = {'true': True, 'false': False}
    strings = values.astype(str).str.lower()
    assert strings.isin(mapping).all()
    return strings.map(mapping).astype(bool)


def aggregate_class(values):
    labels = sorted(set(values))
    if labels == ['strict_tandem_UTR']:
        return 'tandem'
    if labels == ['alternative_terminal_exon_compatible']:
        return 'ALE'
    return labels[0] if len(labels) == 1 else 'mixed'


def audit_contributors(events, background):
    frozen = pd.read_csv(STUDY_EVENTS, sep='\t', float_precision='round_trip')
    frozen = frozen[frozen.contrast.eq('P_M') & frozen.threshold.eq(.02) & frozen.metric.eq('WPP')]
    frozen = frozen[(frozen.species.eq('human') & frozen.gene_id.isin(background.human_gene_id)) |
                    (frozen.species.eq('mouse') & frozen.gene_id.isin(background.mouse_gene_id))]
    assert not frozen.duplicated(['species', 'unit', 'gene_id']).any()
    lookup = {(r.species, r.unit, r.gene_id): r for r in frozen.itertuples()}
    actual = set()
    rows = []
    for key, group in events.groupby(['species', 'study_id', 'gene_id'], sort=True):
        assert key in lookup
        ref = lookup[key]
        expected_routes = set(ref.series_pattern.split(';'))
        observed_routes = set(group.route)
        assert expected_routes == observed_routes and ref.n_routes == len(group)
        actual.add(key)
        differences = [abs(group.delta_DWUI.mean() - ref.effect_mean),
                       abs(group.delta_DWUI.min() - ref.effect_min),
                       abs(group.delta_DWUI.max() - ref.effect_max)]
        assert max(differences) < 1e-12
        rows.append(dict(species=key[0], study_id=key[1], gene_id=key[2],
                         expected_series_pattern=ref.series_pattern,
                         observed_series_pattern=';'.join(sorted(observed_routes)),
                         n_routes=len(group), max_abs_DWUI_difference=max(differences),
                         contribution_sets_identical=True))
    assert actual == set(lookup)
    save(pd.DataFrame(rows), 'frozen_contributor_audit.tsv')


def audit_ale_utrs(events, pairs):
    ale = pairs[pairs.pair_class.eq('both_ALE')]
    utrs = pd.read_csv(UTRS, sep='\t', dtype={'chrom': str})
    by_gene = {(sp, gene): group for (sp, gene), group in utrs.groupby(['species', 'gene_id'])}
    route_rows = []
    for species in ['human', 'mouse']:
        selected = events[events.species.eq(species) & events.gene_id.isin(ale[species + '_gene_id'])]
        for row in selected.itertuples():
            annotations = by_gene.get((species, row.gene_id), utrs.iloc[:0])
            annotations = annotations[annotations.strand.eq(row.strand) &
                annotations.chrom.str.removeprefix('chr').eq(str(row.chrom).removeprefix('chr'))]
            result = dict(species=species, study_id=row.study_id, route=row.route,
                          gene_id=row.gene_id, gene_symbol=row.gene_symbol,
                          driver_class=row.driver_class)
            for end in ['gain', 'loss']:
                cut = getattr(row, end + '_cut')
                base = cut - 1 if row.strand == '+' else cut
                compatible = annotations[annotations.utr_start.sub(1).le(base) & annotations.utr_end.gt(base)]
                result[end + '_retained_base0'] = base
                result[end + '_terminal_UTR_covered'] = not compatible.empty
                result[end + '_UTR_transcripts'] = ';'.join(sorted(compatible.transcript_id.unique()))
            result['both_driver_ends_UTR_covered'] = result['gain_terminal_UTR_covered'] & result['loss_terminal_UTR_covered']
            route_rows.append(result)
    route_table = pd.DataFrame(route_rows)
    save(route_table, 'ALE_driver_terminal_UTR_coverage_routes.tsv')
    gene_rows = []
    for (species, gene), group in route_table.groupby(['species', 'gene_id'], sort=True):
        gene_rows.append(dict(species=species, gene_id=gene, gene_symbol=group.gene_symbol.iloc[0],
            n_routes=len(group), all_routes_gain_terminal_UTR_covered=bool(group.gain_terminal_UTR_covered.all()),
            all_routes_loss_terminal_UTR_covered=bool(group.loss_terminal_UTR_covered.all()),
            all_routes_both_driver_ends_UTR_covered=bool(group.both_driver_ends_UTR_covered.all()),
            n_routes_with_nonUTR_end=int((~group.both_driver_ends_UTR_covered).sum())))
    gene_table = pd.DataFrame(gene_rows)
    save(gene_table, 'ALE_driver_terminal_UTR_coverage_genes.tsv')
    result = ale[KEYS + ['human_gene_symbol', 'mouse_gene_symbol', 'shared_proximal_002', 'shared_proximal_005']].copy()
    for species in ['human', 'mouse']:
        flags = gene_table[gene_table.species.eq(species)][['gene_id', 'all_routes_both_driver_ends_UTR_covered']]
        flags = flags.rename(columns={'gene_id': species + '_gene_id',
                        'all_routes_both_driver_ends_UTR_covered': species + '_all_routes_both_ends_UTR_covered'})
        result = result.merge(flags, on=species + '_gene_id', validate='one_to_one')
    result['both_species_all_routes_both_ends_UTR_covered'] = (
        result.human_all_routes_both_ends_UTR_covered & result.mouse_all_routes_both_ends_UTR_covered)
    save(result, 'both_ALE_terminal_UTR_coverage_pairs.tsv')


def main():
    events = pd.read_csv(EVENTS, sep='\t')
    members = pd.read_csv(MEMBERS, sep='\t', float_precision='round_trip')
    background = pd.read_csv(BACKGROUND, sep='\t')
    background = background[background.contrast.eq('P_M')].copy()
    members = members[members.contrast.eq('P_M')].copy()
    assert len(background) == 751 and not background.duplicated(KEYS).any()
    assert len(events) == 3467 and set(events.species) == {'human', 'mouse'}
    assert not events.duplicated(['species', 'study_id', 'route', 'gene_id']).any()
    assert events.driver_class.notna().all()
    events['full_primary_axis_tandem'] = strict_bool(events.full_primary_axis_tandem)
    events['gain_is_proximal'] = strict_bool(events.gain_is_proximal)
    expected = np.where(events.gain_is_proximal, -1, 1)
    net_sign = np.where(events.delta_DWUI.abs().le(1e-12), 0, np.sign(events.delta_DWUI))
    valid_driver = events.gain_delta.gt(0) & events.loss_delta.lt(0)
    events['driver_pair_sign'] = expected
    events['net_DWUI_sign'] = net_sign.astype(int)
    events['sign_audit'] = np.where(~valid_driver, 'no_valid_gain_loss',
                          np.where(net_sign == 0, 'zero_DWUI',
                          np.where(net_sign == expected, 'concordant', 'discordant')))
    audit_contributors(events, background)
    species_records = []
    for (species, gene), group in events.groupby(['species', 'gene_id'], sort=True):
        study_means = group.groupby('study_id').delta_DWUI.mean()
        species_records.append({
            'species': species, 'gene_id': gene, 'gene_symbol': group.gene_symbol.iloc[0],
            'event_class': aggregate_class(group.driver_class),
            'driver_classes': ';'.join(sorted(set(group.driver_class))),
            'all_full_axis_tandem': bool(group.full_primary_axis_tandem.all()),
            'n_routes': len(group), 'n_studies': group.study_id.nunique(),
            'route_pattern': ';'.join(sorted(group.study_id + ':' + group.route)),
            'delta_DWUI_route_median': group.delta_DWUI.median(),
            'delta_DWUI_route_min': group.delta_DWUI.min(),
            'delta_DWUI_route_max': group.delta_DWUI.max(),
            'delta_DWUI_equal_study_mean': study_means.mean(),
            'n_discordant_driver_routes': int(group.sign_audit.eq('discordant').sum()),
            'n_zero_DWUI_routes': int(group.sign_audit.eq('zero_DWUI').sum()),
        })
    species_genes = pd.DataFrame(species_records)
    assert len(species_genes) == 1502
    pairs = background.drop(columns=['contrast']).copy()
    for species in ['human', 'mouse']:
        group = species_genes[species_genes.species.eq(species)].drop(columns='species')
        group = group.rename(columns={col: species + '_' + col for col in group.columns
                                     if col not in ['gene_id', 'gene_symbol']})
        group = group.rename(columns={'gene_id': species + '_gene_id'})
        # Ortholog symbols are authoritative; retained route symbols have a separate field.
        group = group.rename(columns={'gene_symbol': species + '_route_gene_symbol'})
        pairs = pairs.merge(group, on=species + '_gene_id', how='left', validate='one_to_one')
        assert pairs[species + '_event_class'].notna().all()
    pairs['pair_class'] = np.where(pairs.human_event_class.eq('tandem') & pairs.mouse_event_class.eq('tandem'),
                                  'both_tandem',
                         np.where(pairs.human_event_class.eq('ALE') & pairs.mouse_event_class.eq('ALE'),
                                  'both_ALE', 'other/mixed'))
    pairs['both_full_axis_tandem'] = pairs.human_all_full_axis_tandem & pairs.mouse_all_full_axis_tandem
    pairs['strict_full_axis_sensitivity'] = pairs.pair_class.eq('both_tandem') & pairs.both_full_axis_tandem
    thresholds = sorted(members.threshold.unique())
    assert len(thresholds) == 9
    for threshold, group in members.groupby('threshold'):
        assert len(group) == 751 and not group.duplicated(KEYS).any()
        assert set(map(tuple, group[KEYS].to_numpy())) == set(map(tuple, pairs[KEYS].to_numpy()))
        assert group.human_direction.isin(DIRECTIONS).all() and group.mouse_direction.isin(DIRECTIONS).all()
    detailed = members.merge(pairs, on=KEYS, validate='many_to_one', suffixes=('', '_background'))
    for threshold, count in [(0.02, 27), (0.05, 7)]:
        core = detailed[np.isclose(detailed.threshold, threshold) & detailed.human_direction.eq('proximal')
                        & detailed.mouse_direction.eq('proximal')]
        assert len(core) == count
        core_ids = set(core.human_gene_id)
        pairs['shared_proximal_' + str(threshold).replace('.', '')] = pairs.human_gene_id.isin(core_ids)
    # Save the common classification contract early; all subsequent statistics are deterministic.
    class_path = HERE / 'event_class_pairs.tsv'
    encoded_pairs = pairs.sort_values(KEYS).to_csv(sep='\t', index=False, na_rep='NA', float_format='%.17g').encode()
    if class_path.exists():
        assert class_path.read_bytes() == encoded_pairs, 'Frozen class table would change; stop before writing'
    else:
        class_path.write_bytes(encoded_pairs)
    audit_ale_utrs(events, pairs)
    save(species_genes, 'species_gene_classes.tsv')
    save(events, 'route_sign_audit.tsv')
    save(species_genes.groupby(['species', 'event_class']).size().rename('n_genes').reset_index(),
         'species_class_counts.tsv')
    save(pairs.groupby(['human_event_class', 'mouse_event_class']).size().rename('n_pairs').reset_index(),
         'species_class_matrix.tsv')
    save(pairs.groupby('pair_class').size().rename('n_pairs').reset_index(), 'pair_class_counts.tsv')
    marginals = []
    pair_marginals = []
    cross = []
    effects = []
    scopes = [('all', detailed), ('both_full_axis_tandem', detailed[detailed.strict_full_axis_sensitivity])]
    for scope, subset in scopes:
        for threshold, td in subset.groupby('threshold'):
            for species in ['human', 'mouse']:
                class_col, direction_col = species + '_event_class', species + '_direction'
                for class_, cd in td.groupby(class_col):
                    for direction in DIRECTIONS:
                        chosen = cd[cd[direction_col].eq(direction)]
                        marginals.append(dict(scope=scope, threshold=threshold, species=species,
                            event_class=class_, direction=direction, n_genes=len(chosen),
                            class_background=len(cd), proportion=len(chosen) / len(cd)))
                        effects.append(dict(scope=scope, threshold=threshold, species=species,
                            event_class=class_, direction=direction, n_genes=len(chosen),
                            gene_route_median_DWUI_median=chosen[species + '_delta_DWUI_route_median'].median(),
                            gene_equal_study_mean_DWUI_median=chosen[species + '_delta_DWUI_equal_study_mean'].median(),
                            gene_route_median_DWUI_q25=chosen[species + '_delta_DWUI_route_median'].quantile(.25),
                            gene_route_median_DWUI_q75=chosen[species + '_delta_DWUI_route_median'].quantile(.75)))
            for class_, cd in td.groupby('pair_class'):
                for species in ['human', 'mouse']:
                    for direction in DIRECTIONS:
                        chosen = cd[cd[species + '_direction'].eq(direction)]
                        pair_marginals.append(dict(scope=scope, threshold=threshold, species=species,
                            pair_class=class_, direction=direction, n_genes=len(chosen),
                            class_background=len(cd), proportion=len(chosen) / len(cd)))
                for human in DIRECTIONS:
                    for mouse in DIRECTIONS:
                        chosen = cd[cd.human_direction.eq(human) & cd.mouse_direction.eq(mouse)]
                        cross.append(dict(scope=scope, threshold=threshold, pair_class=class_,
                            human_direction=human, mouse_direction=mouse, n_pairs=len(chosen),
                            class_background=len(cd), proportion=len(chosen) / len(cd)))
    save(pd.DataFrame(marginals), 'recurrence_direction_marginals.tsv')
    save(pd.DataFrame(pair_marginals), 'paired_stratum_direction_marginals.tsv')
    save(pd.DataFrame(cross), 'recurrence_cross_direction.tsv')
    save(pd.DataFrame(effects), 'effect_summary_by_species_class_direction.tsv')
    save(detailed, 'classified_response_members_all_thresholds.tsv')
    cores = detailed[detailed.human_direction.eq('proximal') & detailed.mouse_direction.eq('proximal')]
    save(cores, 'shared_proximal_core_all_thresholds.tsv')
    composition = cores.groupby(['threshold', 'pair_class', 'human_event_class', 'mouse_event_class']).size()
    composition = composition.rename('n_pairs').reset_index()
    composition['threshold_core_n'] = composition.threshold.map(cores.groupby('threshold').size())
    composition['proportion_of_core'] = composition.n_pairs / composition.threshold_core_n
    save(composition, 'shared_proximal_core_composition.tsv')
    save(cores[np.isclose(cores.threshold, .02) | np.isclose(cores.threshold, .05)], 'core_002_005_members.tsv')
    for name in sorted(pairs.pair_class.unique()):
        save(pairs[pairs.pair_class.eq(name)], 'genes_' + name.replace('/', '_') + '.tsv')
    sign_counts = events.groupby(['species', 'driver_class', 'sign_audit']).size().rename('n_routes').reset_index()
    save(sign_counts, 'driver_sign_counts.tsv')
    sign_sets = []
    for threshold in [None, *thresholds]:
        for species in ['human', 'mouse']:
            selected = events[events.species.eq(species)]
            label = 'all_background'
            if threshold is not None:
                core_ids = set(cores[np.isclose(cores.threshold, threshold)][species + '_gene_id'])
                selected = selected[selected.gene_id.isin(core_ids)]
                label = 'shared_proximal_core'
            sign_sets.append(dict(set=label, threshold=threshold, species=species,
                n_genes=selected.gene_id.nunique(), n_routes=len(selected),
                n_concordant_routes=int(selected.sign_audit.eq('concordant').sum()),
                n_discordant_routes=int(selected.sign_audit.eq('discordant').sum()),
                n_genes_any_discordant=selected.loc[selected.sign_audit.eq('discordant'), 'gene_id'].nunique(),
                n_zero_DWUI_routes=int(selected.sign_audit.eq('zero_DWUI').sum()),
                n_no_valid_gain_loss_routes=int(selected.sign_audit.eq('no_valid_gain_loss').sum())))
    save(pd.DataFrame(sign_sets), 'driver_sign_summary.tsv')
    manifest = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [EVENTS, MEMBERS, BACKGROUND, STUDY_EVENTS, UTRS]},
        'rules_sha256': sha(HERE / 'RULES.zh.md'), 'script_sha256': sha(HERE / 'classify_and_describe.py'),
        'n_ortholog_pairs': len(pairs), 'n_species_genes': len(species_genes), 'n_route_events': len(events),
        'thresholds': thresholds, 'statistics_are_descriptive': True,
    }
    (HERE / 'classification_contract.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'pairs': len(pairs), 'classes': pairs.pair_class.value_counts().to_dict(),
                      'full_axis_pairs': int(pairs.strict_full_axis_sensitivity.sum()),
                      'n_discordant_routes': int(events.sign_audit.eq('discordant').sum())}))


if __name__ == '__main__':
    main()
