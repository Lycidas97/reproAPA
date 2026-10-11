"""Frozen S-M reference contrasts and observed-S P-S measurements."""
from pathlib import Path
from collections import defaultdict
import hashlib
import importlib.util
import json
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / 'notes/2026-10-10-marker-reference-sensitivity'
PROTOCOL = HERE.parent / 'PROTOCOL.zh.md'
PROTOCOL_SHA = '85677ecc53b5bca9f6dbc5d59dc1789274e3c130733c23a57b1f62bcb488051a'
STUDIES = {'GSE103221': ('mouse', ['CIP']), 'GSE115943': ('mouse', ['WOT']),
           'GSE118258': ('human', ['GSE118258']), 'GSE150311': ('human', ['GSE150311']),
           'GSE178325': ('human', ['178325_0618', '178325_0809']),
           'GSE221739': ('human', ['GSE221739'])}
SCORES = ['Cycle', 'Adjusted_Cycle']
VALUE_COLS = ['delta_DWUI', 'delta_Cycle_exact', 'delta_Adjusted_Cycle_exact',
              'delta_Cycle_fullcohort', 'delta_Adjusted_Cycle_fullcohort']
DESIGNS = ['main', 'late_M', 'same_library', 'same_time_cross_library']
BACKGROUNDS = ['available_contributors', 'all_candidate_pairs_intersection']


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


spec = importlib.util.spec_from_file_location('frozen_marker_measurement', OLD / 'analyze.py')
measurement = importlib.util.module_from_spec(spec)
spec.loader.exec_module(measurement)


def save(frame, name):
    pd.DataFrame(frame).to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                             compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def encoded(value):
    return json.dumps(np.asarray(value).tolist(), separators=(',', ':'))


def cohorts(z, sid, meta, cycles):
    selected = measurement.frame(z, sid + '_selected_json')
    selected = selected[selected.route.isin(STUDIES[sid][1]) & selected.purity.ge(.8) & selected.node.isin(['S', 'M', 'P'])]
    assert np.array_equal(meta.metacell_id.iloc[selected.matrix_row], selected.metacell_id)
    sampling = measurement.frame(z, sid + '_sampling_json')
    out = {}
    for (route, lib, state), q in selected.groupby(['route', 'library_key', 'node'], sort=True):
        if len(q) < 3 or q.n_cells.sum() < 30:
            continue
        d = sampling[sampling.route.eq(route) & sampling.library_key.eq(lib) & sampling.state.eq(state)]
        assert len(d) == 1
        rec = d.iloc[0].to_dict()
        rr = q.matrix_row.to_numpy(int)
        exact = cycles.loc[q.metacell_id]
        assert np.array_equal(exact.n_cells.to_numpy(), meta.n_cells.iloc[rr].to_numpy())
        assert np.array_equal(exact.library_key.to_numpy(), meta.library_key.iloc[rr].to_numpy())
        rec.update(rows=rr, documented_source=state == 'S' and str(rec['time']) in ['D0', 'source', 'Stage_I_day_0'])
        for score in SCORES:
            rec[score] = float(exact[score].mean())
        out[(route, lib, state)] = rec
    return out


def comparisons(sid, groups, endpoint):
    sm, ps = [], []
    late = set(zip(endpoint.loc[endpoint.contrast.eq('P_M') & endpoint.design.eq('same_library'), 'route'],
                   endpoint.loc[endpoint.contrast.eq('P_M') & endpoint.design.eq('same_library'), 'library_key_b']))
    for route in STUDIES[sid][1]:
        ss = [info for (r, _, state), info in groups.items() if r == route and state == 'S']
        ms = [info for (r, _, state), info in groups.items() if r == route and state == 'M']
        source = [info for info in ss if info['documented_source']]
        main_refs = {info['library_key'] for info in source or ss}
        for m in ms:
            for s in ss:
                schemes = []
                if s['library_key'] in main_refs:
                    schemes.append('main')
                    if (route, m['library_key']) in late:
                        schemes.append('late_M')
                if s['library_key'] == m['library_key']:
                    schemes.append('same_library')
                elif s['time'] == m['time']:
                    schemes.append('same_time_cross_library')
                for design in schemes:
                    sm.append(pair_design(sid, route, m, s, design, 'M', 'S',
                                          'documented_source' if source and design in ['main', 'late_M'] else 'observed_S'))
        if sid in ['GSE150311', 'GSE221739']:
            pp = [info for (r, _, state), info in groups.items() if r == route and state == 'P']
            for p in pp:
                for s in ss:
                    ps.append(pair_design(sid, route, p, s, 'observed_S_all_P', 'P', 'S', 'observed_S'))
    return sm, ps


def pair_design(sid, route, a, b, design, target, reference, source_type):
    rec = dict(species=STUDIES[sid][0], study_id=sid, route=route, design=design,
               target_state=target, reference_state=reference, source_type=source_type,
               library_key_a=a['library_key'], library_key_b=b['library_key'],
               time_a=a['time'], time_b=b['time'], sample_unit_a=a['sample_unit'], sample_unit_b=b['sample_unit'],
               replicate_status_a=a['replicate_status'], replicate_status_b=b['replicate_status'])
    return rec


def profile(counts, rr, meta):
    result = measurement.state(counts, rr, meta)
    for score in SCORES:
        result[score] = float(meta[score].iloc[result['rows']].mean()) if result['n_mc'] else np.nan
    return result


def effects(z, sid, meta, groups, designs, wanted):
    records, eligibility = [], []
    full = measurement.frame(z, sid + '_features_json')
    full['column'] = np.arange(len(full))
    features = full[full.coordinate_eligible & full.matched_source_count_50nt.ge(1) & full.primary_gene_id.isin(wanted)]
    data, indices, indptr = (z[sid + suffix] for suffix in ['_pas_data', '_pas_indices', '_pas_indptr'])
    nrows = int(z[sid + '_pas_shape'][0])
    by_gene = {gene: part for gene, part in features.groupby('primary_gene_id', sort=True) if len(part) >= 2}
    for gene in sorted(wanted):
        part = by_gene.get(gene)
        profiles = {}
        if part is not None:
            counts = measurement.dense_gene(data, indices, indptr, nrows, part.column.to_numpy())
            for key, info in groups.items():
                profiles[key] = profile(counts, info['rows'], meta)
        for d in designs:
            ak = (d['route'], d['library_key_a'], d['target_state'])
            bk = (d['route'], d['library_key_b'], d['reference_state'])
            a, b = profiles.get(ak), profiles.get(bk)
            why = 'fewer_than_two_supported_coordinate_PAS' if part is None else measurement.reason(a, b)
            eligibility.append(d | dict(gene_id=gene, eligible=not why, reason=why,
                n_qualified_MC_a=a['n_mc'] if a else 0, n_qualified_MC_b=b['n_mc'] if b else 0,
                n_qualified_cells_a=a['n_cells'] if a else 0, n_qualified_cells_b=b['n_cells'] if b else 0))
            if why:
                continue
            # M-S was the target hierarchy; sign is reversed only after the calculation.
            sign = -1 if d['target_state'] == 'M' else 1
            joint = ((a['molecules'] + b['molecules']) >= 10) & ((a['nonzero_mc'] + b['nonzero_mc']) >= 2)
            rec = d | dict(gene_id=gene, gene_symbol=part.primary_gene_symbol.iloc[0],
                contrast='S_M' if sign == -1 else 'P_S', setting='supported50_p08',
                delta_DWUI=sign * float((a['usage'] - b['usage']) @ part.position.to_numpy()),
                PAS_ids=encoded(part.atlas_pas_id.to_numpy()), position=encoded(part.position.to_numpy()),
                usage_a=encoded(a['usage']), usage_b=encoded(b['usage']),
                pooled_molecules_a=encoded(a['molecules']), pooled_molecules_b=encoded(b['molecules']),
                nonzero_MC_a=encoded(a['nonzero_mc']), nonzero_MC_b=encoded(b['nonzero_mc']),
                joint_supported_PAS=int(joint.sum()),
                qualified_MC_ids_a=encoded(meta.metacell_id.iloc[a['rows']].to_numpy()),
                qualified_MC_ids_b=encoded(meta.metacell_id.iloc[b['rows']].to_numpy()),
                n_qualified_MC_a=a['n_mc'], n_qualified_MC_b=b['n_mc'],
                n_qualified_cells_a=a['n_cells'], n_qualified_cells_b=b['n_cells'],
                n_fullcohort_MC_a=len(groups[ak]['rows']), n_fullcohort_MC_b=len(groups[bk]['rows']))
            for score in SCORES:
                rec['delta_' + score + '_exact'] = sign * (a[score] - b[score])
                rec['delta_' + score + '_fullcohort'] = sign * (groups[ak][score] - groups[bk][score])
            records.append(rec)
    return pd.DataFrame(records), pd.DataFrame(eligibility)


def aggregate(pairs, by):
    rows = []
    for key, part in pairs.groupby(by, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        part = part.sort_values(['library_key_a', 'library_key_b']).reset_index(drop=True)
        weight = hierarchical_weights(part)
        rec = dict(zip(by, key), n_pairs=len(part), n_target_libraries=part.library_key_a.nunique(),
                   n_reference_libraries=part.library_key_b.nunique(),
                   contributor_pairs=json.dumps(list(zip(part.library_key_a, part.library_key_b)), separators=(',', ':')),
                   contributor_times=json.dumps(list(zip(part.time_a, part.time_b)), separators=(',', ':')),
                   contributor_sample_units=json.dumps(list(zip(part.sample_unit_a, part.sample_unit_b)), separators=(',', ':')),
                   weights=encoded(weight))
        for field in VALUE_COLS:
            rec[field] = float(weight @ part[field].to_numpy())
        rows.append(rec)
    return pd.DataFrame(rows)


def hierarchical_weights(part):
    """Same frozen target/reference hierarchy, without repeated DataFrame merges."""
    keys = list(part[['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b']].itertuples(index=False, name=None))
    captures = defaultdict(int)
    children = [defaultdict(set) for _ in range(4)]
    for key in keys:
        captures[key] += 1
        for level in range(4):
            children[level][key[:level]].add(key[:level + 1])
    weights = np.array([1 / captures[key] / np.prod([len(children[level][key[:level]]) for level in range(4)]) for key in keys])
    assert abs(weights.sum() - 1) < 1e-12
    return weights


def summary(row, current):
    rec = row | dict(N=len(current), status='descriptive' if len(current) else 'not_estimable',
                      reason='' if len(current) else 'no_genes_pass_fixed_background')
    if len(current):
        v = current.delta_DWUI
        rec.update(mean_DWUI=float(v.mean()), median_DWUI=float(v.median()),
                   proximal=int(v.le(-.02).sum()), distal=int(v.ge(.02).sum()),
                   below_threshold=int(v.between(-.02, .02, inclusive='neither').sum()),
                   n_gene_negative=int(v.lt(0).sum()), n_gene_positive=int(v.gt(0).sum()))
        for field in VALUE_COLS[1:]:
            rec['mean_' + field] = float(current[field].mean())
            rec['median_' + field] = float(current[field].median())
    else:
        for field in ['mean_DWUI', 'median_DWUI', 'proximal', 'distal', 'below_threshold', 'n_gene_negative', 'n_gene_positive']:
            rec[field] = np.nan
        for field in VALUE_COLS[1:]:
            rec['mean_' + field] = rec['median_' + field] = np.nan
    return rec


def sm_summaries(pairs, candidates, fixed_sets):
    routes, summaries, stability = [], [], []
    for sid, (species, route_names) in STUDIES.items():
        for route in route_names:
            for design in DESIGNS:
                cp = candidates[candidates.study_id.eq(sid) & candidates.route.eq(route) & candidates.design.eq(design)]
                p = pairs[pairs.study_id.eq(sid) & pairs.route.eq(route) & pairs.design.eq(design)]
                for name in ['core27', 'complete_PM_proximal']:
                    ids = set(fixed_sets.loc[fixed_sets.species.eq(species) & fixed_sets.gene_set.eq(name), 'gene_id'])
                    candidate_count = len(cp)
                    measurable = p[p.gene_id.isin(ids)]
                    counts = measurable.groupby('gene_id').size()
                    fixed = set(counts.index[counts.eq(candidate_count)]) if candidate_count else set()
                    for background in BACKGROUNDS:
                        current = measurable if background == 'available_contributors' else measurable[measurable.gene_id.isin(fixed)]
                        context = dict(species=species, study_id=sid, route=route, contrast='S_M', design=design,
                                       background=background, gene_set=name, n_fixed_roster=len(ids), n_candidate_pairs=candidate_count,
                                       source_type='documented_source' if route in ['CIP', 'WOT', 'GSE118258', '178325_0618', '178325_0809'] and design in ['main', 'late_M'] else 'observed_S')
                        if len(current):
                            result = aggregate(current, ['gene_id'])
                            for field, value in context.items():
                                result[field] = value
                            routes.append(result)
                        else:
                            result = pd.DataFrame()
                        summaries.append(summary(context, result))
                        if not candidate_count:
                            continue
                        for field, level in [('library_key_a', 'target_library'), ('time_a', 'target_time')]:
                            for label in sorted(cp[field].unique()):
                                part = current[current[field].eq(label)]
                                v = aggregate(part, ['gene_id']) if len(part) else pd.DataFrame()
                                stability.append(summary(context | dict(level=level, target_label=label), v))
    return pd.concat(routes, ignore_index=True), pd.DataFrame(summaries), pd.DataFrame(stability)


def main():
    assert digest(PROTOCOL) == PROTOCOL_SHA
    paths = [PROTOCOL, OLD / 'inputs.npz', OLD / 'analyze.py',
             ROOT / 'notes/2026-10-10-proliferation-apa/metacell_cycle_scores.tsv.gz',
             HERE / 'fixed_gene_sets.tsv', HERE / 'preflight_validation.json',
             HERE.parent / 'contrast-dose/true_source_PS_common_orthologs.tsv']
    hashes = {str(path): digest(path) for path in paths}
    fixed_sets = pd.read_csv(HERE / 'fixed_gene_sets.tsv', sep='\t')
    observed_panel = pd.read_csv(paths[-1], sep='\t')
    assert len(observed_panel) == 1499 and observed_panel.human_gene_id.is_unique
    if '--aggregate-only' in sys.argv:
        return recover_aggregation(hashes, fixed_sets)
    cycles = pd.read_csv(paths[3], sep='\t', float_precision='round_trip').set_index('metacell_id')
    assert cycles.index.is_unique
    z = np.load(OLD / 'inputs.npz', allow_pickle=False)
    all_sm, all_sme, all_ps, all_pse, sm_designs, ps_designs = [], [], [], [], [], []
    for sid, (species, routes) in STUDIES.items():
        meta = measurement.frame(z, sid + '_meta_json')
        for score in SCORES:
            meta[score] = meta.metacell_id.map(cycles[score])
        group = cohorts(z, sid, meta, cycles)
        sm, ps = comparisons(sid, group, measurement.frame(z, sid + '_endpoint_json'))
        sm_designs.extend(sm); ps_designs.extend(ps)
        wanted = set(fixed_sets.loc[fixed_sets.species.eq(species) & fixed_sets.gene_set.eq('complete_PM_proximal'), 'gene_id'])
        smp, sme = effects(z, sid, meta, group, sm, wanted)
        all_sm.append(smp); all_sme.append(sme)
        if ps:
            psp, pse = effects(z, sid, meta, group, ps, set(observed_panel.human_gene_id))
            all_ps.append(psp); all_pse.append(pse)
        print('COMPLETE', sid, 'SM_eligible_pairs', len(smp), 'observed_PS_pairs', len(psp) if ps else 0, flush=True)
    sm = pd.concat(all_sm, ignore_index=True); ps = pd.concat(all_ps, ignore_index=True)
    smc = pd.DataFrame(sm_designs); psc = pd.DataFrame(ps_designs)
    save(smc, 'SM_candidate_pairs.tsv'); save(psc, 'observed_PS_candidate_pairs.tsv')
    save(pd.concat(all_sme, ignore_index=True), 'SM_pair_eligibility.tsv.gz')
    save(pd.concat(all_pse, ignore_index=True), 'observed_PS_pair_eligibility.tsv.gz')
    save(sm, 'SM_pair_effects.tsv.gz'); save(ps, 'observed_PS_pair_effects.tsv.gz')
    finish_aggregation(sm, ps, smc, fixed_sets, hashes)


def finish_aggregation(sm, ps, smc, fixed_sets, hashes):
    # Deliver observed P-S first so the separate contrast analysis need not wait.
    psr = aggregate(ps, ['species', 'study_id', 'route', 'gene_id', 'design', 'source_type'])
    psr['contrast'] = 'P_S'; psr['setting'] = 'supported50_p08'
    save(psr, 'observed_PS_route_effects.tsv.gz')
    print('OBSERVED_PS_AGGREGATION_COMPLETE', len(psr), flush=True)
    rr, ss, ls = sm_summaries(sm, smc, fixed_sets)
    save(rr, 'SM_route_gene_effects.tsv.gz'); save(ss, 'SM_route_summary.tsv'); save(ls, 'SM_library_date_summary.tsv')
    assert all(digest(path) == expected for path, expected in hashes.items())
    validation = dict(status='COMPLETE_AWAITING_INDEPENDENT_QA', protocol_sha256=PROTOCOL_SHA, input_sha256=hashes,
        script_sha256=digest(Path(__file__)), SM_eligible_pair_rows=len(sm), observed_PS_eligible_pair_rows=len(ps),
        SM_route_gene_rows=len(rr), SM_summary_rows=len(ss), observed_PS_route_gene_rows=len(psr),
        no_gene_MC_or_species_p=True, support='Original selected S/M/P mutually exclusive; joint nonzero_MC are distinct.',
        primary_score='Same gene-qualified MC as APA; full pre-gene cohort raw/adjusted scores are fixed sensitivities.')
    (HERE / 'analysis_validation.json').write_text(json.dumps(validation, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(validation | {'input_sha256': 'see_file'}, ensure_ascii=False), flush=True)


def recover_aggregation(hashes, fixed_sets):
    sources = [HERE / name for name in ['SM_pair_effects.tsv.gz', 'observed_PS_pair_effects.tsv.gz',
               'SM_candidate_pairs.tsv', 'SM_pair_eligibility.tsv.gz', 'observed_PS_pair_eligibility.tsv.gz']]
    source_hashes = {str(path): digest(path) for path in sources}
    hashes.update(source_hashes)
    sm, ps, smc = [pd.read_csv(path, sep='\t', float_precision='round_trip') for path in sources[:3]]
    assert len(sm) == 48769 and len(ps) == 104903
    # Compare diverse contributor patterns against the unchanged original function.
    checks, maximum = 0, 0.
    for _, part in pd.concat([sm, ps], ignore_index=True).groupby(['study_id', 'design', 'gene_id'], sort=True):
        if checks >= 200:
            break
        error = float(np.max(np.abs(hierarchical_weights(part) - measurement.weights(part))))
        maximum = max(maximum, error); checks += 1
    assert maximum < 1e-14
    finish_aggregation(sm, ps, smc, fixed_sets, hashes)
    assert all(digest(path) == expected for path, expected in source_hashes.items())
    (HERE / 'aggregation_recovery_validation.json').write_text(json.dumps(dict(status='PASS',
        sealed_pair_and_qualification_sha256=source_hashes, unchanged_source_rows={'SM': len(sm), 'observed_PS': len(ps)},
        frozen_weight_function_comparisons=checks, max_weight_abs_error=maximum,
        reason='Stopped expensive repeated pandas aggregation; recovered from sealed effect vectors using equivalent hierarchy.'), indent=2) + '\n')


if __name__ == '__main__':
    main()
