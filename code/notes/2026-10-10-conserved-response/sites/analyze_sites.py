"""Fixed response contrasts, strict UTR qualification and unchanged top drivers."""
from pathlib import Path
from collections import defaultdict
import argparse
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
SETTING = 'supported50_p08'
CODES = {'H1': 'GSE118258', 'H2': 'GSE150311', 'H3': 'GSE178325', 'H4': 'GSE221739', 'W': 'GSE115943', 'C': 'GSE103221'}
INPUTS = {}


def read(path, **kwargs):
    path = Path(path).resolve()
    with path.open('rb') as stream:
        INPUTS[str(path.relative_to(ROOT))] = hashlib.file_digest(stream, 'sha256').hexdigest()
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kwargs)


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def top_drivers(delta):
    """Pick top sites once; an unmapped maximum cannot be replaced by a runner-up."""
    keys = sorted(delta)
    if not keys:
        return None, None
    gain = max(keys, key=lambda key: delta[key])
    loss = min(keys, key=lambda key: delta[key])
    return gain, loss


def merge_native(ids, values, mapping):
    merged = defaultdict(float)
    for pas, value in zip(ids, values):
        merged[mapping.get(pas, 'UNMAPPED:' + pas)] += float(value)
    assert abs(sum(merged.values()) - 1) < 1e-12
    return dict(merged)


def average_profiles(profiles):
    """Equal series/conditions within study, then equal independent studies."""
    by_study = defaultdict(list)
    for profile in profiles:
        by_study[profile['study_id']].append(profile)
    ids = sorted(set().union(*(set(p['target']) | set(p['reference']) for p in profiles)))
    study_values = {}
    for study, group in sorted(by_study.items()):
        target = {pas: float(np.mean([p['target'].get(pas, 0.) for p in group])) for pas in ids}
        reference = {pas: float(np.mean([p['reference'].get(pas, 0.) for p in group])) for pas in ids}
        study_values[study] = dict(target=target, reference=reference,
                                  delta={pas: target[pas] - reference[pas] for pas in ids})
    delta = {pas: float(np.mean([value['delta'][pas] for value in study_values.values()])) for pas in ids}
    return delta, study_values


def supports_drivers(profiles, gain, loss, threshold=.02):
    """All measurable series of a study must support the exact same two sites."""
    study = defaultdict(list)
    for p in profiles:
        g = p['target'].get(gain, 0.) - p['reference'].get(gain, 0.)
        l = p['target'].get(loss, 0.) - p['reference'].get(loss, 0.)
        study[p['study_id']].append(g >= threshold and l <= -threshold)
    return sum(all(values) for values in study.values())


def mean_supports_drivers(profiles, gain, loss, threshold=.02):
    """Descriptive study averages; this is not the original all-series study vote."""
    delta, studies = average_profiles(profiles)
    return sum(value['delta'][gain] >= threshold and value['delta'][loss] <= -threshold for value in studies.values())


def matched_mass(profiles, matched_ids):
    return min(sum(p[node].get(pas, 0.) for pas in matched_ids)
               for p in profiles for node in ['target', 'reference'])


def strict_utr(species, gene, profiles, eligibility):
    studies = sorted({p['study_id'] for p in profiles})
    if any(not p.get('axis_verified', True) for p in profiles):
        return False, 'actual_naive_axis_not_verified', '', np.nan
    rows = [eligibility.get((species, study, gene)) for study in studies]
    if any(row is None for row in rows):
        return False, 'missing_frozen_full_axis_UTR_qualification', '', np.nan
    if any(not row.tandem_utr_compatible for row in rows):
        return False, 'non_tandem_full_primary_axis', '', float(np.median([row.native_axis_span_nt for row in rows]))
    common = set.intersection(*(set(row.compatible_transcripts.split(';')) for row in rows))
    if not common:
        return False, 'no_common_terminal_UTR_transcript', '', float(np.median([row.native_axis_span_nt for row in rows]))
    if any(not p.get('pair_tandem', True) for p in profiles):
        return False, 'actual_naive_condition_not_tandem_qualified', '', float(np.median([row.native_axis_span_nt for row in rows]))
    return True, 'eligible', ';'.join(sorted(common)), float(np.median([row.native_axis_span_nt for row in rows]))


def main_profiles(frame, contrast, frozen, route_effects, mapping):
    profiles = {}
    for sp in ['human', 'mouse']:
        gene = frame[sp + '_gene_id']
        measurements = frozen[(frozen.species == sp) & (frozen.gene_id == gene) & (frozen.contrast == contrast)]
        output = []
        for row in measurements.itertuples():
            needed = row.series_pattern.split(';')
            for series in needed:
                route = series.removeprefix('178325_') if row.unit == 'GSE178325' else series
                if row.unit == 'GSE178325':
                    route = 'Guan2022_' + route
                # Routes use original table labels; the pattern inventory is checked below.
                current = route_effects[(route_effects.species == sp) & (route_effects.study_id == row.unit) &
                                        (route_effects.gene_id == gene) & (route_effects.contrast == contrast)]
                if row.unit == 'GSE178325':
                    current = current[current.route.str.endswith(series.removeprefix('178325_'))]
                elif row.unit in ['GSE103221', 'GSE115943']:
                    current = current[current.route == series]
                assert len(current) == 1, (sp, gene, contrast, row.unit, series, len(current))
                item = current.iloc[0]
                ids = json.loads(item.PAS_ids)
                target, reference = json.loads(item.target_usage), json.loads(item.reference_usage)
                output.append(dict(study_id=row.unit, series=series,
                                   original_delta=float(item.delta_DWUI), pair_tandem=True,
                                   target=merge_native(ids, target, mapping), reference=merge_native(ids, reference, mapping)))
        assert len({p['study_id'] for p in output}) >= 2
        profiles[sp] = output
    return profiles


def naive_profiles(frame, pairs, frozen, route_effects, mapping, dual_vectors=None):
    gene = frame['human_gene_id']
    selected = pairs[pairs.gene_id.eq(gene) & pairs.node_a.eq('N') & pairs.node_b.eq('P') &
                     pairs.library_a.str.endswith(('6', '8')) & pairs.library_b.str.endswith('7')]
    assert len(selected) == 2 and selected.library_a.nunique() == 2
    human = []
    for row in selected.itertuples():
        ids = json.loads(row.PAS_ids)
        human.append(dict(study_id='GSE150311', series=row.library_a,
                          original_delta=float(row.delta_DWUI), pair_tandem=bool(row.tandem_utr_compatible),
                          axis_verified=bool(row.actual_axis_verified),
                          target=merge_native(ids, json.loads(row.usage_a), mapping),
                          reference=merge_native(ids, json.loads(row.usage_b), mapping)))
    mouse = []
    gene = frame['mouse_gene_id']
    if dual_vectors is not None:
        selected = dual_vectors[dual_vectors.species.eq('mouse') & dual_vectors.gene_id.eq(gene)]
        assert len(selected) == 2 and selected.study_id.nunique() == 2
        for row in selected.itertuples():
            ids = json.loads(row.PAS_ids)
            mouse.append(dict(study_id=row.study_id, series=row.route, original_delta=float(row.DWUI_PM), pair_tandem=True,
                              target=merge_native(ids, json.loads(row.P_usage), mapping),
                              reference=merge_native(ids, json.loads(row.M_usage), mapping)))
        return dict(human=human, mouse=mouse)
    measurements = frozen[(frozen.species == 'mouse') & (frozen.gene_id == gene) & (frozen.contrast == 'P_M')]
    for row in measurements.itertuples():
        current = route_effects[route_effects.species.eq('mouse') & route_effects.study_id.eq(row.unit) &
                                route_effects.gene_id.eq(gene) & route_effects.contrast.eq('P_M')]
        assert len(current) == 1
        item = current.iloc[0]; ids = json.loads(item.PAS_ids)
        mouse.append(dict(study_id=row.unit, series=item.route, original_delta=float(item.delta_DWUI), pair_tandem=True,
                          target=merge_native(ids, json.loads(item.target_usage), mapping),
                          reference=merge_native(ids, json.loads(item.reference_usage), mapping)))
    assert len({p['study_id'] for p in mouse}) == 2
    return dict(human=human, mouse=mouse)


def characterize(frame, set_name, contrast, role, recurrent_side, profiles, eligibility, catalog, pairs, tolerance):
    row = dict(set_name=set_name, contrast=contrast, tolerance=tolerance,
               human_gene_id=frame['human_gene_id'], mouse_gene_id=frame['mouse_gene_id'],
               human_gene_symbol=frame.get('human_gene_symbol', frame.get('h_gene_symbol', '')),
               mouse_gene_symbol=frame.get('mouse_gene_symbol', frame.get('m_gene_symbol', '')),
               parent_case=bool(frame['case']), context_role=role, recurrent_side=recurrent_side)
    delta_by_species, driver_ids, study_values = {}, {}, {}
    for sp in ['human', 'mouse']:
        group = profiles[sp]
        strict, status, transcripts, span = strict_utr(sp, row[sp + '_gene_id'], group, eligibility)
        delta, studies = average_profiles(group)
        gain, loss = top_drivers(delta)
        row.update({sp + '_strict_UTR': strict, sp + '_UTR_status': status,
                    sp + '_common_transcripts': transcripts, sp + '_native_span': span,
                    sp + '_studies': len(studies), sp + '_series': len(group),
                    sp + '_original_delta': float(np.mean([np.mean([p['original_delta'] for p in group if p['study_id'] == study]) for study in studies])),
                    sp + '_gain_id': gain, sp + '_loss_id': loss,
                    sp + '_gain_delta': delta[gain], sp + '_loss_delta': delta[loss],
                    sp + '_gain_support_studies': supports_drivers(group, gain, loss),
                    sp + '_gain_support_studies_mean': mean_supports_drivers(group, gain, loss),
                    sp + '_driver_threshold': gain != loss and delta[gain] >= .02 and delta[loss] <= -.02,
                    sp + '_unmapped_mass_max': max(sum(value for key, value in p[node].items() if key.startswith('UNMAPPED:')) for p in group for node in ['target', 'reference'])})
        delta_by_species[sp], driver_ids[sp], study_values[sp] = delta, (gain, loss), studies
    current = pairs[pairs.human_gene_id.eq(row['human_gene_id']) & pairs['reciprocal' + str(tolerance)]]
    native_union = {sp: set(delta_by_species[sp]) for sp in ['human', 'mouse']}
    current = current[current.human_id.isin(native_union['human']) & current.mouse_id.isin(native_union['mouse'])].copy()
    if len(current):
        h = catalog.loc[current.human_id]
        m = catalog.loc[current.mouse_id]
        current['h_signed_cut'] = h.cut_site.to_numpy() * np.where(h.strand.eq('+'), 1, -1)
        current['m_signed_cut'] = m.cut_site.to_numpy() * np.where(m.strand.eq('+'), 1, -1)
        current = current.sort_values('h_signed_cut')
    order_ok = len(current) >= 2 and np.all(np.diff(current.m_signed_cut) > 0) and np.all(np.diff(current.h_signed_cut) > 0)
    row.update(matched_PAS=len(current), order_eligible=bool(order_ok))
    h_to_m = current.set_index('human_id').mouse_id.to_dict()
    m_to_h = current.set_index('mouse_id').human_id.to_dict()
    mapped_drivers = []
    for sp, to_other in [('human', h_to_m), ('mouse', m_to_h)]:
        gain, loss = driver_ids[sp]
        row[sp + '_drivers_have_partners'] = gain in to_other and loss in to_other and gain != loss
        for kind, key in [('gain', gain), ('loss', loss)]:
            if key in to_other:
                other = to_other[key]
                mapped_drivers.append((key, other) if sp == 'human' else (other, key))
        matched_ids = list(current[sp + '_id'])
        row[sp + '_min_retained_mass'] = matched_mass(profiles[sp], matched_ids)
        row[sp + '_matched_fraction'] = len(current) / len(native_union[sp]) if native_union[sp] else np.nan
        row[sp + '_matched_rank_delta'] = np.nan
        row[sp + '_matched_distance_delta'] = np.nan
        if order_ok and row[sp + '_min_retained_mass'] >= .8:
            cuts = current['h_signed_cut' if sp == 'human' else 'm_signed_cut'].to_numpy(float)
            rank, distance = np.linspace(0, 1, len(cuts)), (cuts - cuts.min()) / np.ptp(cuts)
            by_study = defaultdict(list)
            for p in profiles[sp]:
                values = {}
                for node in ['target', 'reference']:
                    u = np.asarray([p[node].get(key, 0.) for key in matched_ids])
                    values[node] = u / u.sum()
                by_study[p['study_id']].append(values['target'] - values['reference'])
            delta = np.mean([np.mean(values, axis=0) for values in by_study.values()], axis=0)
            row[sp + '_matched_rank_delta'] = float(delta @ rank)
            row[sp + '_matched_distance_delta'] = float(delta @ distance)
    row['strict_UTR'] = row['human_strict_UTR'] and row['mouse_strict_UTR']
    row['mass_eligible'] = order_ok and min(row['human_min_retained_mass'], row['mouse_min_retained_mass']) >= .8
    row['driver_judgeable'] = row['human_drivers_have_partners'] and row['mouse_drivers_have_partners'] and row['human_driver_threshold'] and row['mouse_driver_threshold']
    hg, hl = driver_ids['human']; mg, ml = driver_ids['mouse']
    row['double_driver_correspondence'] = row['driver_judgeable'] and h_to_m.get(hg) == mg and h_to_m.get(hl) == ml
    row['strict_double_driver'] = row['strict_UTR'] and row['mass_eligible'] and row['double_driver_correspondence']
    row['two_study_exact_event_both'] = row['strict_double_driver'] and row['human_gain_support_studies'] >= 2 and row['mouse_gain_support_studies'] >= 2
    row['human_condition_consistency'] = row['strict_double_driver'] and set_name == 'hN' and all(
        p['target'].get(hg, 0.) - p['reference'].get(hg, 0.) >= .02 and
        p['target'].get(hl, 0.) - p['reference'].get(hl, 0.) <= -.02 for p in profiles['human'])
    if not row['strict_UTR']:
        reason = 'strict_UTR_inapplicable'
    elif not order_ok:
        reason = 'fewer_than_two_corresponding_sites_or_crossing_order'
    elif not row['mass_eligible']:
        reason = 'retained_mass_below_0.8'
    elif not row['driver_judgeable']:
        reason = 'maximum_driver_unmatched_or_below_usage_threshold'
    elif not row['double_driver_correspondence']:
        reason = 'maximum_driver_pairs_differ'
    else:
        reason = 'eligible_double_driver'
    row['event_status'] = reason
    # Controls only use the recurring side, without selecting or forcing the other side's drivers.
    driver_pairs = []
    if role == 'shared_recurrent' and row['strict_double_driver']:
        driver_pairs = [(hg, mg), (hl, ml)]
    elif role == 'single_species_recurrent' and not row['parent_case'] and row['strict_UTR'] and row['mass_eligible']:
        side = recurrent_side
        if row[side + '_driver_threshold'] and row[side + '_drivers_have_partners']:
            gain, loss = driver_ids[side]
            driver_pairs = [(gain, h_to_m[gain]), (loss, h_to_m[loss])] if side == 'human' else [(m_to_h[gain], gain), (m_to_h[loss], loss)]
    row['sequence_candidate'] = bool(driver_pairs)
    row['sequence_driver_pairs'] = json.dumps(sorted(set(driver_pairs)), separators=(',', ':'))
    traces = []
    for sp in ['human', 'mouse']:
        gain, loss = driver_ids[sp]
        for p in profiles[sp]:
            traces.append(dict(set_name=set_name, contrast=contrast, tolerance=tolerance, human_gene_id=row['human_gene_id'],
                               mouse_gene_id=row['mouse_gene_id'], species=sp, study_id=p['study_id'], series=p['series'],
                               original_delta=p['original_delta'], gain_id=gain, loss_id=loss,
                               gain_delta=p['target'].get(gain, 0.) - p['reference'].get(gain, 0.),
                               loss_delta=p['target'].get(loss, 0.) - p['reference'].get(loss, 0.),
                               target=json.dumps(p['target'], sort_keys=True, separators=(',', ':')),
                               reference=json.dumps(p['reference'], sort_keys=True, separators=(',', ':'))))
    return row, traces


def run(membership_path):
    config = json.loads((HERE.parent / 'protocol.json').read_text())
    assert config['protocol_sha256'] == hashlib.sha256((HERE.parent / 'PROTOCOL.zh.md').read_bytes()).hexdigest()
    members = read(membership_path, low_memory=False, dtype={'measurable_mask': str})
    main = members[members.set_name.eq('main48') & members.threshold.eq(.02)].copy()
    naive = members[members.set_name.eq('hN79') & members.threshold.eq(.02)].copy()
    for column in ['P_S_measurable', 'P_M_measurable']:
        main[column] = main[column].map(lambda value: value is True or value == 'True')
    for name in ['PROTOCOL.zh.md', 'protocol.json', 'SITE-CONTRAST-ADDENDUM.zh.md', 'SITE-STUDY-VOTE-CLARIFICATION.zh.md']:
        path = HERE.parent / name
        INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    assert len(main) == 587 and main.case.sum() == 48
    assert len(naive) == 290 and naive.case.sum() == 79
    f = read(TABLES / 'recurrence_controls/recurrence_frozen_covariates.tsv.gz')
    f = f[f.threshold.eq(.02)]
    e = read(TABLES / 'tandem/tandem_utr_eligibility.tsv.gz')
    eligibility = {(r.species, r.study_id, r.gene_id): r for r in e.itertuples()}
    mapping_frame = read(TABLES / 'cis/detected_to_annotation.tsv.gz')
    assert mapping_frame.native_pas_id.is_unique
    mapping = mapping_frame.set_index('native_pas_id').id.to_dict()
    catalog = read(TABLES / 'cis/annotation_catalog.tsv.gz').set_index('id')
    pairs = read(TABLES / 'cis/reciprocal_one_to_one.tsv.gz')
    effects = read(TABLES / 'endpoints/route_effects.tsv.gz')
    effects = effects[effects.setting.eq(SETTING) &
                      ((effects.contrast.eq('P_S') & effects.design.eq('true_source')) |
                       (effects.contrast.eq('P_M') & effects.design.eq('same_library')))]
    naive_pairs = read(TABLES / 'naive_state/library_pair_effects.tsv.gz')
    naive_pairs = naive_pairs[naive_pairs.version.eq('v2_primary') & naive_pairs.setting.eq(SETTING)]
    naive_validation = read(HERE / 'naive_actual_axis_validation.tsv')
    naive_pairs = naive_pairs.merge(naive_validation[['gene_id', 'actual_axis_verified']], on='gene_id', how='left', validate='many_to_one')
    naive_pairs['actual_axis_verified'] = naive_pairs.actual_axis_verified.fillna(False).astype(bool)
    dual_vectors = read(TABLES / 'matched_pms/matched_state_vectors.tsv.gz')
    dual_vectors = dual_vectors[dual_vectors.setting.eq(SETTING)]
    rows, traces = [], []
    cache = {}
    for (population, threshold), current_members in members.groupby(['set_name', 'threshold'], sort=True):
        set_name = 'main' if population in ['main48', 'list45'] else 'hN'
        if set_name == 'main':
            current_members = current_members.copy()
            for column in ['P_S_measurable', 'P_M_measurable']:
                current_members[column] = current_members[column].map(lambda value: value is True or value == 'True')
        for index, frame in enumerate(current_members.to_dict('records')):
            contexts = []
            if set_name == 'main':
                for contrast in ['P_S', 'P_M']:
                    if not frame[contrast + '_measurable']:
                        continue
                    h, m = frame[contrast + '_human_direction'], frame[contrast + '_mouse_direction']
                    role = 'shared_recurrent' if h == m and h != 'nonrecurrent' else 'single_species_recurrent' if (h == 'nonrecurrent') != (m == 'nonrecurrent') else 'other_measurable'
                    recurring = 'human' if h != 'nonrecurrent' else 'mouse'
                    key = ('main', frame['human_gene_id'], contrast)
                    if key not in cache:
                        cache[key] = main_profiles(frame, contrast, f, effects, mapping)
                    contexts.append((contrast, role, recurring, cache[key]))
            else:
                role = 'shared_recurrent' if frame['case'] else 'single_species_recurrent'
                recurring = 'human' if frame['case'] else frame['control_side']
                key = (population == 'dual46', frame['human_gene_id'], 'hN')
                if key not in cache:
                    cache[key] = naive_profiles(frame, naive_pairs, f, effects, mapping, dual_vectors if population == 'dual46' else None)
                contexts.append(('hN_P__mouse_P_M', role, recurring, cache[key]))
            for contrast, role, side, profiles in contexts:
                for tolerance in [25, 50]:
                    row, trace = characterize(frame, set_name, contrast, role, side, profiles, eligibility, catalog, pairs, tolerance)
                    row.update(source_population=population, gene_threshold=threshold)
                    for item in trace:
                        item.update(source_population=population, gene_threshold=threshold)
                    rows.append(row); traces.extend(trace)
            if index % 100 == 0:
                print(population, threshold, index, len(rows), flush=True)
    output = pd.DataFrame(rows)
    save(output, 'site_event_contexts.tsv.gz')
    save(pd.DataFrame(traces), 'site_usage_traces.tsv.gz')
    stages = ['strict_UTR', 'order_eligible', 'mass_eligible', 'driver_judgeable', 'double_driver_correspondence', 'two_study_exact_event_both']
    attrition = []
    group_columns = ['source_population', 'gene_threshold', 'set_name', 'tolerance', 'parent_case', 'context_role']
    for key, group in output.groupby(group_columns, sort=True):
        current = group.copy()
        for stage in ['all_contexts'] + stages:
            if stage != 'all_contexts':
                current = current[current[stage]]
            attrition.append(dict(zip(group_columns, key), stage=stage,
                                  contexts=len(current), unique_orthologs=current.human_gene_id.nunique()))
    save(pd.DataFrame(attrition), 'event_attrition.tsv')
    save(output[output.parent_case], 'case_site_event_contexts.tsv')
    primary = output[output.source_population.isin(['main48', 'hN79']) & output.gene_threshold.eq(.02)]
    query_ids = sorted(set().union(*(set(item) for values in primary[primary.sequence_candidate].sequence_driver_pairs.map(json.loads) for item in values)))
    sequence = read(TABLES / 'cis/site_sequence_features.tsv.gz')
    for sp in ['human', 'mouse']:
        query = read(TABLES / f'cis/liftover/{sp}_sequence_hg38_mm10.audit.tsv.gz')
        query = query[query.id.isin(query_ids)]
        if sp == 'human':
            all_query = query
        else:
            all_query = pd.concat([all_query, query], ignore_index=True)
    all_query = all_query.merge(sequence[['id', 'phastcons']], on='id', validate='one_to_one')
    save(all_query, 'fixed_window_queries.tsv')
    (HERE / 'site_input_hashes.json').write_text(json.dumps(INPUTS, indent=2) + '\n')
    print('complete', len(output), 'contexts; fixed-window IDs', len(all_query), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--membership', required=True)
    args = parser.parse_args()
    run(args.membership)
