"""Frozen marker reference sensitivity, with baseline acceptance before new effects."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
PROTOCOL_SHA = '1fc9ff8f14ee74302e032b52d4a6a68b6e89463b73f5bc457fa5ece266c75040'
ADDENDUM_SHA = 'ddbfe768ef6210ebda71cdf584e70dd36a0f5e9f3691ed48d5716445758b956a'
QA_ADDENDUM_SHA = 'f0b3be44713f72e8845c57384c978e5ef6a0f87f3679c22054df411fb1b2836b'
LOST_KEYS_SHA = 'f9fb32fd29a90af67ff19daf5ef57a29274196e50c8a55567997e3c6ffa0524d'
STUDIES = {'GSE103221': ('mouse', ['CIP']), 'GSE115943': ('mouse', ['WOT']),
           'GSE118258': ('human', ['GSE118258']), 'GSE150311': ('human', ['GSE150311']),
           'GSE178325': ('human', ['178325_0618', '178325_0809']),
           'GSE221739': ('human', ['GSE221739'])}
P_CORE = ['DPPA2', 'DPPA4', 'SALL4', 'TDGF1']
P_NATIVE = {'human': ['DPPA4', 'DNMT3B', 'ESRG', 'L1TD1', 'TDGF1'],
            'mouse': ['DPPA4', 'ESRRB', 'ZFP42', 'DPPA5A', 'TDGF1']}
ECM = ['COL1A1', 'COL3A1', 'DCN', 'LUM', 'FN1']
VERSIONS = {'core_third_strict': ('core', 1/3, True),
            'core_quartile_strict': ('core', 1/4, True),
            'core_third_relaxed': ('core', 1/3, False),
            'native_third_strict': ('native', 1/3, True)}
CONTRASTS = {'D00': ('original_P', 'original_M'), 'D01': ('original_P', 'Mstar'),
             'D10': ('Pstar', 'original_M'), 'D11': ('Pstar', 'Mstar')}
INPUTS = {}
RESOLUTIONS = []


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path, **kwargs):
    INPUTS[str(path)] = digest(path)
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kwargs)


def frame(z, key):
    return pd.read_json(io.StringIO(str(z[key])), orient='records')


def save(data, name):
    pd.DataFrame(data).to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                            compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def encode(value):
    return json.dumps(np.asarray(value).tolist(), separators=(',', ':'))


def weights(data):
    """Exact source date/unit then target date/unit/capture hierarchy of the main pipeline."""
    levels = [['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'],
              ['time_a', 'sample_unit_a', 'time_b'], ['time_a', 'sample_unit_a'], ['time_a'], []]
    q = data.reset_index(drop=True).copy()
    out = 1 / q.groupby(levels[0], dropna=False).time_a.transform('size').to_numpy(float)
    for child, parent in zip(levels[:-1], levels[1:]):
        distinct = q[child].drop_duplicates()
        if not parent:
            out /= len(distinct)
        else:
            n = distinct.groupby(parent, dropna=False).size().rename('_n').reset_index()
            out /= q[parent].merge(n, on=parent, how='left', validate='many_to_one')._n.to_numpy()
    assert abs(out.sum() - 1) < 1e-12
    return out


def state(counts, rows, meta):
    """Equal MC proportions after the frozen 20-UMI gene gate, never pooled UMI usage."""
    c = counts[rows]
    rows = np.asarray(rows)[c.sum(1) >= 20]
    c = counts[rows]
    return dict(rows=rows, n_mc=len(rows), n_cells=int(meta.n_cells.iloc[rows].sum()),
                usage=(c / c.sum(1)[:, None]).mean(0) if len(c) else np.zeros(counts.shape[1]),
                molecules=c.sum(0), nonzero_mc=(c > 0).sum(0))


def reason(a, b):
    if a is None or b is None:
        return 'cohort_absent_or_below_pre_gene_3MC_30cells'
    if min(a['n_mc'], b['n_mc']) < 3:
        return 'insufficient_gene_qualified_metacells'
    if min(a['n_cells'], b['n_cells']) < 30:
        return 'insufficient_gene_qualified_cells'
    joint = (a['molecules'] + b['molecules'] >= 10) & (a['nonzero_mc'] + b['nonzero_mc'] >= 2)
    return '' if joint.sum() >= 2 else 'fewer_than_two_joint_supported_PAS'


def dense_gene(data, indices, indptr, rows, columns):
    out = np.zeros((rows, len(columns)), dtype=data.dtype)
    for j, column in enumerate(columns):
        a, b = indptr[column:column + 2]
        out[indices[a:b], j] = data[a:b]
    return out


def programme(cp, positions):
    values = np.log1p(cp[:, positions])
    sd = values.std(0, ddof=0)
    variable = sd > 1e-12
    standardized = np.zeros_like(values)
    standardized[:, variable] = (values[:, variable] - values[:, variable].mean(0)) / sd[variable]
    return standardized.mean(1), int(variable.sum()), (cp[:, positions] >= 1).sum(1)


def classify(pscore, escore, pexpressed, eexpressed, fraction, strict):
    high_fraction = .75 if fraction == .25 else 2/3
    pl, ph = np.quantile(pscore, [fraction, high_fraction], method='linear')
    el, eh = np.quantile(escore, [fraction, high_fraction], method='linear')
    result = np.full(len(pscore), 'unassigned/mixed', dtype=object)
    if ph <= pl or eh <= el:
        return result, np.full(len(pscore), 'programme_quantiles_not_separated', dtype=object), (pl, ph, el, eh), False
    pre_p = (pscore >= ph) & (escore <= el)
    pre_m = (escore >= eh) & (pscore <= pl)
    p = pre_p & (pexpressed >= 2) & ((eexpressed <= 1) if strict else True)
    m = pre_m & (eexpressed >= 2) & ((pexpressed <= 1) if strict else True)
    assert not (p & m).any()
    result[p] = 'Pstar'; result[m] = 'Mstar'
    why = np.full(len(pscore), 'relative_quantiles_not_passed', dtype=object)
    why[pre_p | pre_m] = 'absolute_marker_support_not_passed'
    why[p | m] = 'assigned'
    return result, why, (pl, ph, el, eh), True


def original_groups(meta):
    groups = {}
    old = meta[meta.in_planned_route & meta.original_selected & meta.original_purity.ge(.8)]
    for (route, library, node), q in old[old.original_node.isin(['P', 'M'])].groupby(['route', 'library_key', 'original_node']):
        if len(q) >= 3 and q.n_cells.sum() >= 30:
            groups[(route, library, 'original_' + node)] = q.index.to_numpy()
    return groups


def planned_pairs(z, sid):
    design = frame(z, sid + '_endpoint_json')
    q = design[design.route.isin(STUDIES[sid][1]) & design.contrast.eq('P_M') & design.design.eq('same_library')].copy()
    assert q.library_key_a.eq(q.library_key_b).all()
    assert not q.duplicated(['route', 'library_key_a']).any()
    return q


def features_for(z, sid, setting):
    f = frame(z, sid + '_features_json')
    f['column'] = np.arange(len(f))
    q = f[f.coordinate_eligible].copy()
    if setting == 'supported50':
        q = q[q.matched_source_count_50nt.ge(1)]
    return f, q


def gene_matrices(z, sid, setting):
    full, features = features_for(z, sid, setting)
    data, indices, indptr = (z[sid + suffix] for suffix in ['_pas_data', '_pas_indices', '_pas_indptr'])
    nrows = int(z[sid + '_pas_shape'][0])
    for gene, part in features.groupby('primary_gene_id', sort=True):
        if len(part) >= 2:
            yield gene, part, dense_gene(data, indices, indptr, nrows, part.column.to_numpy())


def aggregate_baseline(pairs):
    records = []
    if not len(pairs):
        return pd.DataFrame()
    for (sid, route, gene), q in pairs.groupby(['study_id', 'route', 'gene_id'], sort=True):
        ww = weights(q)
        a = ww @ np.stack(q.target_usage); b = ww @ np.stack(q.reference_usage)
        records.append(dict(study_id=sid, route=route, gene_id=gene, n_pairs=len(q),
                            target_usage=a, reference_usage=b, delta_usage=a - b,
                            delta_DWUI=float((a - b) @ q.position.iloc[0])))
        pairs.loc[q.index, 'weight'] = ww
    return pd.DataFrame(records)


def baseline(z):
    expected_pairs = read(TABLES / 'state_vectors/pair_effects.primary.tsv.gz', usecols=[
        'study_id', 'route', 'gene_id', 'library_key_a', 'library_key_b', 'contrast', 'design', 'setting',
        'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b', 'weight', 'delta_DWUI', 'PAS_ids',
        'target_usage', 'reference_usage', 'delta_usage', 'metacell_ids_a', 'metacell_ids_b',
        'n_mc_a', 'n_mc_b', 'n_cells_a', 'n_cells_b', 'molecules_a', 'molecules_b'])
    expected_pairs = expected_pairs[expected_pairs.setting.eq('supported50_p08') & expected_pairs.contrast.eq('P_M') & expected_pairs.design.eq('same_library')]
    expected_routes = read(TABLES / 'state_vectors/route_effects.primary.tsv.gz', usecols=[
        'study_id', 'route', 'gene_id', 'contrast', 'design', 'setting', 'delta_DWUI',
        'target_usage', 'reference_usage', 'delta_usage', 'PAS_ids', 'n_pairs'])
    expected_routes = expected_routes[expected_routes.setting.eq('supported50_p08') & expected_routes.contrast.eq('P_M') & expected_routes.design.eq('same_library')]
    reconstructed = []
    study_rows = []
    for sid in STUDIES:
        meta = frame(z, sid + '_meta_json'); groups = original_groups(meta); design = planned_pairs(z, sid)
        for gene, sites, counts in gene_matrices(z, sid, 'supported50'):
            states = {key: state(counts, rows, meta) for key, rows in groups.items()}
            for d in design.to_dict('records'):
                a = states.get((d['route'], d['library_key_a'], 'original_P'))
                b = states.get((d['route'], d['library_key_b'], 'original_M'))
                if reason(a, b):
                    continue
                reconstructed.append(dict(d, study_id=sid, gene_id=gene,
                    PAS_ids=sites.atlas_pas_id.tolist(), position=sites.position.to_numpy(),
                    target_usage=a['usage'], reference_usage=b['usage'], delta_usage=a['usage'] - b['usage'],
                    delta_DWUI=float((a['usage'] - b['usage']) @ sites.position.to_numpy()),
                    metacell_ids_a=meta.metacell_id.iloc[a['rows']].tolist(),
                    metacell_ids_b=meta.metacell_id.iloc[b['rows']].tolist(),
                    n_mc_a=a['n_mc'], n_mc_b=b['n_mc'], n_cells_a=a['n_cells'], n_cells_b=b['n_cells'],
                    molecules_a=a['molecules'], molecules_b=b['molecules']))
        study_rows.append(dict(study_id=sid, n_reconstructed_pairs=sum(r['study_id'] == sid for r in reconstructed)))
        print('BASELINE', study_rows[-1], flush=True)
    pairs = pd.DataFrame(reconstructed); routes = aggregate_baseline(pairs)
    checks = []
    key_p = ['study_id', 'route', 'gene_id', 'library_key_a', 'library_key_b']
    key_r = ['study_id', 'route', 'gene_id']
    for label, actual, expected, keys in [('pair', pairs, expected_pairs, key_p), ('route', routes, expected_routes, key_r)]:
        a = actual.set_index(keys); e = expected.set_index(keys)
        assert a.index.is_unique and e.index.is_unique
        assert set(a.index) == set(e.index), (label, 'baseline_key_mismatch', len(a), len(e))
        e = e.loc[a.index]
        for key in a.index:
            ar, er = a.loc[key], e.loc[key]
            out = dict(zip(keys, key)); out['level'] = label
            for column in ['target_usage', 'reference_usage', 'delta_usage']:
                out[column + '_max_abs_difference'] = float(np.max(np.abs(ar[column] - np.array(json.loads(er[column])))))
            out['delta_DWUI_abs_difference'] = abs(ar.delta_DWUI - er.delta_DWUI)
            if label == 'pair':
                assert ar.PAS_ids == json.loads(er.PAS_ids), ('PAS_order', key)
                for column in ['metacell_ids_a', 'metacell_ids_b']:
                    assert ar[column] == json.loads(er[column]), (column, key)
                for column in ['molecules_a', 'molecules_b']:
                    assert np.array_equal(ar[column], json.loads(er[column])), (column, key)
                for column in ['n_mc_a', 'n_mc_b', 'n_cells_a', 'n_cells_b']:
                    assert ar[column] == er[column], (column, key)
                for column in ['time_a', 'time_b', 'sample_unit_a', 'sample_unit_b']:
                    assert str(ar[column]) == str(er[column]), (column, key)
                out['weight_abs_difference'] = abs(ar.weight - er.weight)
            else:
                assert ar.n_pairs == er.n_pairs, ('n_pairs', key)
            checks.append(out)
    results = pd.DataFrame(checks)
    differences = [c for c in results if c.endswith('_difference')]
    maximum = float(results[differences].max().max())
    assert maximum < 1e-12, ('baseline_numeric_mismatch', maximum)
    save(results[results.level.eq('pair')], 'baseline_pair_replay.tsv.gz')
    save(results[results.level.eq('route')], 'baseline_route_replay.tsv.gz')
    for name, data in [('baseline_pair_effects.tsv.gz', pairs), ('baseline_route_effects.tsv.gz', routes)]:
        out = data.copy()
        for column in ['PAS_ids', 'position', 'target_usage', 'reference_usage', 'delta_usage',
                       'metacell_ids_a', 'metacell_ids_b', 'molecules_a', 'molecules_b']:
            if column in out:
                out[column] = out[column].map(encode)
        save(out, name)
    acceptance = dict(status='PASS', protocol_sha256=PROTOCOL_SHA, input_sha256=digest(HERE / 'inputs.npz'),
                      n_pairs=len(pairs), n_routes=len(routes), max_abs_difference=maximum,
                      studies=study_rows, baseline_source_hashes=INPUTS.copy())
    (HERE / 'baseline_validation.json').write_text(json.dumps(acceptance, indent=2) + '\n')
    return pairs, routes


class TSVStream:
    """Stream the complete qualification matrices rather than retaining millions of rows."""
    def __init__(self, name, fields):
        self.raw = (HERE / name).open('wb')
        self.gzip = gzip.GzipFile(fileobj=self.raw, mode='wb', mtime=0)
        self.text = io.TextIOWrapper(self.gzip, encoding='utf-8', newline='')
        self.writer = csv.DictWriter(self.text, fields, delimiter='\t', extrasaction='ignore')
        self.writer.writeheader()
        self.n = 0

    def write(self, row):
        self.writer.writerow(row)
        self.n += 1

    def close(self):
        self.text.close(); self.raw.close()


def metadata(z, sid):
    meta = frame(z, sid + '_meta_json')
    states = frame(z, 'frozen_states_json').set_index('metacell_id')
    meta['frozen_state_group'] = meta.metacell_id.map(states.state_group).fillna('unknown')
    meta['frozen_state_purity'] = meta.metacell_id.map(states.state_purity)
    label = 'coarse_state' if 'coarse_state' in meta else 'working_comparison_group'
    meta['native_original_label'] = meta[label].fillna('unknown')
    meta['native_original_label_field'] = label
    meta['species'] = STUDIES[sid][0]; meta['study_id'] = sid
    return meta


def marker_mapping(z):
    mapping = frame(z, 'mapping_json')
    broad = read(ROOT / 'notes/2026-09-13-integration-validation/inputs/S6-integration-orthologs.tsv')
    strict = frame(z, 'orthologs_json')
    records = []
    for row in mapping.to_dict('records'):
        side = row['species']; gene = row['gene_id']
        rows = broad[broad[side + '_gene_id'].eq(gene)]
        one = strict[strict[side + '_gene_id'].eq(gene)]
        records.append(dict(row, broad_membership=bool(len(rows)), strict_16092_membership=bool(len(one)),
                            broad_n_pairs=len(rows), broad_documented_one_to_one=encode(rows.documented_one_to_one),
                            broad_orthology_confidence=encode(rows.orthology_confidence),
                            broad_partner_gene_ids=encode(rows[('mouse' if side == 'human' else 'human') + '_gene_id']),
                            strict_partner_gene_ids=encode(one[('mouse' if side == 'human' else 'human') + '_gene_id']),
                            panel_roles=encode([p for p, genes in [('core_P', P_CORE), ('native_P', P_NATIVE[side]), ('ECM', ECM)] if row['marker'] in genes])))
    save(records, 'marker_mapping.tsv')


def score_and_groups(z, sid, meta):
    markers = z[sid + '_markers'].tolist()
    counts = z[sid + '_marker_counts']
    totals = z[sid + '_total_gex']
    positions = {m: i for i, m in enumerate(markers)}
    groups = original_groups(meta)
    design = planned_pairs(z, sid)
    old_libraries = set(zip(design.route, design.library_key_a))
    groups = {key: rows for key, rows in groups.items() if key[:2] in old_libraries}
    candidate_groups = {}
    for route, library in old_libraries:
        q = meta[meta.route.eq(route) & meta.library_key.eq(library) & meta.original_selected & meta.original_purity.ge(.8)]
        for node in ['P', 'M']:
            candidate_groups[(route, library, 'original_' + node)] = q.index[q.original_node.eq(node)].to_numpy()
    classifications, libraries, gradients, variability, expression = [], [], [], [], []
    for (route, library), q in meta[meta.in_planned_route].groupby(['route', 'library_key'], sort=True):
        rr = q.index.to_numpy()
        valid_total = bool((totals[rr] > 0).all())
        cp = np.divide(counts[rr] * 1e4, totals[rr, None], out=np.full((len(rr), len(markers)), np.nan), where=totals[rr, None] > 0)
        for j, marker in enumerate(markers):
            log = np.log1p(cp[:, j])
            variability.append(dict(study_id=sid, route=route, library_key=library, marker=marker,
                                    n_candidate_mc=len(rr), mean_logCP10K=float(log.mean()),
                                    SD_logCP10K_ddof0=float(log.std(ddof=0)), variable=bool(log.std(ddof=0) > 1e-12)))
            expression.extend(dict(study_id=sid, route=route, library_key=library, metacell_id=meta.metacell_id.iloc[r],
                                   marker=marker, gene_id=z[sid + '_marker_gene_ids'][j], integer_count=int(counts[r, j]),
                                   total_gex=int(totals[r]), CP10K=float(cp[k, j]), logCP10K=float(log[k])) for k, r in enumerate(rr))
        scores = {}
        for panel in ['core', 'native']:
            pmarkers = P_CORE if panel == 'core' else P_NATIVE[STUDIES[sid][0]]
            ps, pv, pe = programme(cp, [positions[m] for m in pmarkers])
            es, ev, ee = programme(cp, [positions[m] for m in ECM])
            scores[panel] = (ps, pv, pe, es, ev, ee)
        for version, (panel, fraction, strict) in VERSIONS.items():
            ps, pv, pe, es, ev, ee = scores[panel]
            why = '' if valid_total else 'nonpositive_total_GEX'
            if pv < 2 or ev < 2:
                why = why or 'fewer_than_two_variable_markers_in_programme'
            if why:
                assignment = np.full(len(rr), 'unassigned/mixed', dtype=object)
                reasons = np.full(len(rr), why, dtype=object)
                cuts = (np.nan,) * 4; classifiable = False
            else:
                assignment, reasons, cuts, classifiable = classify(ps, es, pe, ee, fraction, strict)
                if not classifiable:
                    why = 'programme_quantiles_not_separated'
            pl, ph, el, eh = cuts
            base = dict(species=STUDIES[sid][0], study_id=sid, route=route, library_key=library,
                        version=version, panel=panel, classifiable=classifiable,
                        reason=why or 'classifiable', n_candidate_mc=len(rr), n_candidate_cells=int(q.n_cells.sum()),
                        P_variable_markers=pv, ECM_variable_markers=ev,
                        P_q_low=pl, P_q_high=ph, ECM_q_low=el, ECM_q_high=eh,
                        n_Pstar=int((assignment == 'Pstar').sum()), n_Mstar=int((assignment == 'Mstar').sum()),
                        cells_Pstar=int(q.n_cells.to_numpy()[assignment == 'Pstar'].sum()),
                        cells_Mstar=int(q.n_cells.to_numpy()[assignment == 'Mstar'].sum()),
                        n_unassigned=int((assignment == 'unassigned/mixed').sum()),
                        n_original_P=int(q.original_node.eq('P').sum()), n_original_M=int(q.original_node.eq('M').sum()),
                        n_not_in_original_selected=int((~q.original_selected).sum()))
            libraries.append(base)
            for k, r in enumerate(rr):
                classifications.append(dict(base, metacell_id=meta.metacell_id.iloc[r], matrix_row=int(r), n_cells=int(meta.n_cells.iloc[r]),
                    original_selected=bool(meta.original_selected.iloc[r]), original_node=meta.original_node.iloc[r],
                    original_purity=meta.original_purity.iloc[r], native_original_label=meta.native_original_label.iloc[r],
                    native_original_label_field=meta.native_original_label_field.iloc[r], frozen_state_group=meta.frozen_state_group.iloc[r],
                    P_score=float(ps[k]), ECM_score=float(es[k]), P_expressed_markers=int(pe[k]),
                    ECM_expressed_markers=int(ee[k]), classification=assignment[k], classification_reason=reasons[k]))
            for label in ['Pstar', 'Mstar']:
                take = rr[assignment == label]
                candidate_groups[(route, library, version + '_' + label)] = take
                if classifiable and len(take) >= 3 and meta.n_cells.iloc[take].sum() >= 30:
                    groups[(route, library, version + '_' + label)] = take
        # Gradient is defined inside original M, on the full-candidate programme scale.
        old_m = q[q.original_selected & q.original_purity.ge(.8) & q.original_node.eq('M')]
        local = pd.Index(rr).get_indexer(old_m.index)
        for panel in ['core', 'native']:
            ps, pv, pe, es, ev, ee = scores[panel]
            rec = dict(species=STUDIES[sid][0], study_id=sid, route=route, library_key=library, panel=panel,
                       n_original_M=len(old_m), cells_original_M=int(old_m.n_cells.sum()), P_variable_markers=pv,
                       q_low=np.nan, q_high=np.nan, status='not_estimable', reason='')
            why = 'original_same_library_design_absent' if (route, library) not in old_libraries else ''
            if not valid_total or pv < 2:
                why = why or 'P_programme_unavailable'
            if len(old_m) < 3:
                why = why or 'fewer_than_three_original_M'
            cuts = np.quantile(ps[local], [1/3, 2/3], method='linear') if len(local) and valid_total else [np.nan, np.nan]
            rec['q_low'], rec['q_high'] = cuts
            if not np.isfinite(cuts).all() or cuts[1] <= cuts[0]:
                why = why or 'M_programme_quantiles_not_separated'
            masks = {'low': ps[local] <= cuts[0], 'middle': (ps[local] > cuts[0]) & (ps[local] < cuts[1]), 'high': ps[local] >= cuts[1]}
            for layer, mask in masks.items():
                take = old_m.index.to_numpy()[mask]
                candidate_groups[(route, library, 'gradient_' + panel + '_' + layer)] = take
                rec['n_' + layer] = len(take); rec['cells_' + layer] = int(meta.n_cells.iloc[take].sum())
                if len(take) < 3 or meta.n_cells.iloc[take].sum() < 30:
                    why = why or 'M_layer_below_3MC_30cells'
            if (route, library, 'original_P') not in groups:
                why = why or 'original_P_below_3MC_30cells'
            if not why:
                rec['status'] = 'pre_gene_eligible'
                for layer, mask in masks.items():
                    groups[(route, library, 'gradient_' + panel + '_' + layer)] = old_m.index.to_numpy()[mask]
            rec['reason'] = why
            gradients.append(rec)
    return groups, candidate_groups, classifications, libraries, gradients, variability, expression


def sampling_pairs(z, sid, meta):
    sampling = frame(z, sid + '_sampling_json')
    columns = ['route', 'library_key', 'time', 'sample_unit', 'replicate_status']
    q = sampling[columns].drop_duplicates()
    assert not q.duplicated(['route', 'library_key']).any()
    planned = meta[meta.in_planned_route][['route', 'library_key']].drop_duplicates()
    q = planned.merge(q, how='left', on=['route', 'library_key'], validate='one_to_one')
    missing = q[q.time.isna() | q.sample_unit.isna()]
    if len(missing):
        source = ROOT / 'table/production-v15/upstream/S2-sequencing-libraries.tsv'
        titles = read(source, usecols=['study', 'sample_id', 'geo_sample_title'])
        rules = {
            ('GSE115943', 'GSM3195673'): ('GSM3195673', 'D6_Dox_C2', 'D6', 'bio_well_D6_Dox_C2', 'biological_independent', 'bio_well_D6_Dox_C1'),
            ('GSE221739', 'GSE221739::D9A_rep1_GSM6894036'): ('GSM6894036', 'D9, pool A, replicate 1', 'D9', 'GSE221739:D9:poolA', 'pooled_chip_sample', 'GSE221739:D9:poolA')
        }
        for index, row in missing.iterrows():
            key = (sid, row.library_key)
            assert key in rules, ('unresolved_sampling_metadata', key)
            sample, title, time, unit, status, neighbor_unit = rules[key]
            source_row = titles[titles.study.eq(sid) & titles.sample_id.eq(sample)]
            assert len(source_row) == 1 and source_row.geo_sample_title.iloc[0] == title
            neighbor = q[q.time.eq(time) & q.sample_unit.eq(neighbor_unit) & q.replicate_status.eq(status)]
            assert len(neighbor), ('frozen_neighbor_absent', key)
            assert meta.loc[meta.library_key.eq(row.library_key), 'experimental_time'].eq(time).all()
            q.loc[index, ['time', 'sample_unit', 'replicate_status']] = [time, unit, status]
            RESOLUTIONS.append(dict(study_id=sid, route=row.route, library_key=row.library_key, time=time,
                sample_unit=unit, replicate_status=status, source_sample_id=sample, source_title=title,
                resolution='original_title_and_frozen_neighbor_rule', frozen_neighbor_libraries=encode(neighbor.library_key),
                source_path=str(source), source_sha256=INPUTS[str(source)], adds_original_endpoint=False))
    records = []
    for row in q.to_dict('records'):
        rec = dict(study_id=sid, route=row['route'], library_key=row['library_key'])
        for side in ['a', 'b']:
            rec['library_key_' + side] = row['library_key']
            rec['time_' + side] = row['time']; rec['sample_unit_' + side] = row['sample_unit']
            rec['replicate_status_' + side] = row['replicate_status']
        records.append(rec)
    return records


def summaries(data, background, version, setting, level='route', candidate_libraries=None):
    records = []
    for sid, (species, routes) in STUDIES.items():
        for route in routes:
            available_libraries = [] if level == 'route' else candidate_libraries[(sid, route)] if candidate_libraries is not None else data[(data.study_id == sid) & (data.route == route)].library_key.unique()
            scopes = [('', data[(data.study_id == sid) & (data.route == route)])] if level == 'route' else [
                (lib, data[(data.study_id == sid) & (data.route == route) & (data.library_key == lib)])
                for lib in available_libraries]
            for lib, subset in scopes:
                for contrast in CONTRASTS:
                    x = subset[subset.contrast.eq(contrast)]
                    delta = x.delta_DWUI.to_numpy(float)
                    for threshold in [.02, .05]:
                        prox = delta <= -threshold; distal = delta >= threshold
                        records.append(dict(species=species, study_id=sid, route=route, library_key=lib,
                            version=version, setting=setting, background=background, level=level,
                            contrast=contrast, threshold=threshold, status='descriptive' if len(x) else 'not_estimable',
                            N=len(x), n_proximal=int(prox.sum()), n_distal=int(distal.sum()),
                            n_within_threshold=int((~prox & ~distal).sum()),
                            fraction_proximal=float(prox.mean()) if len(x) else np.nan,
                            fraction_distal=float(distal.mean()) if len(x) else np.nan,
                            mean_delta=float(delta.mean()) if len(x) else np.nan,
                            median_delta=float(np.median(delta)) if len(x) else np.nan,
                            mean_proximal_conditional=float(delta[prox].mean()) if prox.any() else np.nan,
                            n_contributor_libraries=len(set(lib for items in x.contributor_libraries for lib in json.loads(items))) if len(x) and 'contributor_libraries' in x else int(bool(len(x))),
                            max_libraries_per_gene=int(x.n_libraries.max()) if len(x) and 'n_libraries' in x else int(bool(len(x)))))
    return records


def aggregate_new(pairs, common=False):
    """Every gene's four common effects share contributor libraries and exact weights."""
    if not len(pairs):
        return pd.DataFrame(columns=['study_id', 'route', 'gene_id', 'contrast', 'delta_DWUI', 'n_libraries'])
    records, weight_cache = [], {}
    for (sid, route, gene), q in pairs.groupby(['study_id', 'route', 'gene_id'], sort=True):
        joint = q.groupby('library_key').contrast.nunique()
        libs = set(joint[joint.eq(4)].index)
        for contrast, part in q.groupby('contrast', sort=True):
            if common:
                part = part[part.library_key.isin(libs)]
            if not len(part):
                continue
            signature = tuple(part[['library_key', 'time_a', 'sample_unit_a', 'time_b', 'sample_unit_b']].itertuples(index=False, name=None))
            if signature not in weight_cache:
                weight_cache[signature] = weights(part)
            ww = weight_cache[signature]
            pairs.loc[part.index, 'common_weight' if common else 'full_weight'] = ww
            if common:
                pairs.loc[part.index, 'common_four'] = True
            a = ww @ np.stack(part.target_usage); b = ww @ np.stack(part.reference_usage)
            record = dict(study_id=sid, route=route, gene_id=gene, contrast=contrast, n_libraries=len(part),
                contributor_libraries=encode(part.library_key), weights=encode(ww),
                target_usage=encode(a), reference_usage=encode(b), PAS_ids=part.PAS_ids.iloc[0],
                delta_DWUI=float((a-b) @ part.position.iloc[0]), background='common_four' if common else 'full')
            records.append(record)
    return pd.DataFrame(records)


def bridge_identities(pairs, routes, version, setting):
    records = []
    for level, q, keys in [('library', pairs[pairs.common_four], ['study_id', 'route', 'library_key', 'gene_id']),
                            ('route', routes, ['study_id', 'route', 'gene_id'])]:
        if not len(q):
            continue
        for key, part in q.groupby(keys, sort=True):
            d = part.set_index('contrast').delta_DWUI
            assert set(d.index) == set(CONTRASTS)
            observed = d.D11 - d.D00; expected = (d.D01 - d.D00) + (d.D10 - d.D00)
            error = abs(observed - expected)
            assert error < 1e-12, ('bridge_identity', key, error)
            records.append(dict(zip(keys, key), level=level, version=version, setting=setting,
                                observed_change=float(observed), reference_change=float(d.D01-d.D00),
                                target_change=float(d.D10-d.D00), abs_error=float(error)))
    return records


def gradient_aggregate(pairs):
    records = []
    if not len(pairs):
        return pd.DataFrame(columns=['study_id', 'route', 'gene_id', 'panel', 'delta_low', 'delta_middle', 'delta_high', 'delta_whole', 'difference_high_minus_low'])
    for (sid, route, panel, gene), q in pairs.groupby(['study_id', 'route', 'panel', 'gene_id'], sort=True):
        ww = weights(q); pairs.loc[q.index, 'weight'] = ww
        p = ww @ np.stack(q.P_usage)
        m = {layer: ww @ np.stack(q['M_' + layer + '_usage']) for layer in ['low', 'middle', 'high', 'whole']}
        delta = {layer: float((p-m[layer]) @ q.position.iloc[0]) for layer in m}
        observed = delta['high'] - delta['low']
        expected = float((m['low'] - m['high']) @ q.position.iloc[0])
        assert abs(observed - expected) < 1e-12
        records.append(dict(study_id=sid, route=route, panel=panel, gene_id=gene, n_libraries=len(q),
                            contributor_libraries=encode(q.library_key), weights=encode(ww),
                            **{'delta_' + layer: value for layer, value in delta.items()},
                            difference_high_minus_low=observed, direct_M_low_minus_high=expected,
                            identity_abs_error=abs(observed-expected), monotonicity=monotonicity(delta)))
    return pd.DataFrame(records)


def monotonicity(delta):
    lo, mid, hi = [delta[layer] for layer in ['low', 'middle', 'high']]
    if max(lo, mid, hi) - min(lo, mid, hi) < 1e-12:
        return 'flat_within_1e-12'
    if lo <= mid <= hi:
        return 'increasing'
    if lo >= mid >= hi:
        return 'decreasing'
    return 'nonmonotonic'


def gradient_summary(data, level='route', candidate_libraries=None):
    records = []
    for sid, (species, routes) in STUDIES.items():
        for route in routes:
            libraries = [''] if level == 'route' else candidate_libraries[(sid, route)]
            for library in libraries:
                for panel in ['core', 'native']:
                    q = data[data.study_id.eq(sid) & data.route.eq(route) & data.panel.eq(panel)]
                    if level == 'library':
                        q = q[q.library_key.eq(library)]
                    for threshold in [.02, .05]:
                        rec = dict(species=species, study_id=sid, route=route, library_key=library, panel=panel,
                                   level=level, threshold=threshold, N=len(q), status='descriptive' if len(q) else 'not_estimable',
                                   mean_difference_high_minus_low=float(q.difference_high_minus_low.mean()) if len(q) else np.nan,
                                   n_increasing=int(q.monotonicity.eq('increasing').sum()) if len(q) else 0,
                                   n_decreasing=int(q.monotonicity.eq('decreasing').sum()) if len(q) else 0,
                                   n_nonmonotonic=int(q.monotonicity.eq('nonmonotonic').sum()) if len(q) else 0)
                        for layer in ['low', 'middle', 'high', 'whole']:
                            delta = q['delta_' + layer].to_numpy(float)
                            rec['mean_' + layer] = float(delta.mean()) if len(q) else np.nan
                            rec['n_proximal_' + layer] = int((delta <= -threshold).sum())
                            rec['n_distal_' + layer] = int((delta >= threshold).sum())
                            rec['fraction_proximal_' + layer] = float((delta <= -threshold).mean()) if len(q) else np.nan
                        rec['fraction_proximal_high_minus_low'] = rec['fraction_proximal_high'] - rec['fraction_proximal_low']
                        records.append(rec)
    return records


def new_effects(z):
    marker_mapping(z)
    metas, cohorts, candidates, sampling = {}, {}, {}, {}
    classifications, library_class, gradient_class, variability, expression = [], [], [], [], []
    candidate_libraries = {}
    for sid in STUDIES:
        meta = metadata(z, sid); metas[sid] = meta
        groups, allgroups, c, l, g, v, e = score_and_groups(z, sid, meta)
        cohorts[sid] = groups; candidates[sid] = allgroups
        classifications.extend(c); library_class.extend(l); gradient_class.extend(g)
        variability.extend(v); expression.extend(e)
        sampling[sid] = sampling_pairs(z, sid, meta)
        for route in STUDIES[sid][1]:
            candidate_libraries[(sid, route)] = sorted(meta.loc[meta.in_planned_route & meta.route.eq(route), 'library_key'].unique())
    save(pd.concat(metas.values(), ignore_index=True), 'metacell_metadata.tsv.gz')
    save(classifications, 'metacell_classification.tsv.gz')
    save(library_class, 'library_classification.tsv')
    save(gradient_class, 'gradient_library_eligibility.tsv')
    save(variability, 'marker_variability.tsv')
    save(expression, 'metacell_marker_expression.tsv.gz')
    save(RESOLUTIONS, 'sampling_metadata_resolution.tsv')
    conf = pd.DataFrame(classifications)
    conf['original_node_for_confusion'] = conf.original_node.fillna('not_in_old_selected')
    confusion = conf.groupby(['species', 'study_id', 'route', 'library_key', 'version', 'classification',
        'original_selected', 'original_node_for_confusion', 'native_original_label', 'frozen_state_group'], dropna=False).agg(
            n_mc=('metacell_id', 'size'), n_cells=('n_cells', 'sum')).reset_index()
    save(confusion, 'classification_confusion.tsv')
    print('CLASSIFICATION_COMPLETE', len(classifications), flush=True)
    profile_fields = ['species', 'study_id', 'route', 'library_key', 'setting', 'gene_id', 'gene_symbol', 'cohort', 'group_id',
        'status', 'reason', 'n_pre_gene_mc', 'n_pre_gene_cells', 'n_mc', 'n_cells', 'qualified_MC_ids',
        'qualified_MC_pas_depths', 'PAS_ids', 'position', 'pooled_molecules', 'nonzero_MC', 'usage', 'raw_counts_source']
    profile = TSVStream('gene_group_profiles.tsv.gz', profile_fields)
    eligibility = TSVStream('pair_eligibility.tsv.gz', ['study_id', 'route', 'library_key', 'setting', 'version', 'contrast',
        'gene_id', 'status', 'reason', 'n_mc_a', 'n_mc_b', 'n_cells_a', 'n_cells_b', 'n_joint_PAS', 'target_group', 'reference_group'])
    gradient_eligibility = TSVStream('gradient_gene_eligibility.tsv.gz', ['study_id', 'route', 'library_key', 'panel', 'gene_id',
        'status', 'reason', 'n_mc_P', 'n_mc_low', 'n_mc_middle', 'n_mc_high', 'n_cells_P', 'n_cells_low',
        'n_cells_middle', 'n_cells_high', 'n_joint_low', 'n_joint_middle', 'n_joint_high'])
    eligible_pairs, gradient_pairs, library_eligibility, background, axes, cohort_coverage = [], [], [], [], [], []
    for sid, (species, routes) in STUDIES.items():
        meta, groups, allgroups = metas[sid], cohorts[sid], candidates[sid]
        design = planned_pairs(z, sid); oldkeys = set(zip(design.route, design.library_key_a))
        for (route, library, cohort), rr in allgroups.items():
            cohort_coverage.append(dict(study_id=sid, route=route, library_key=library, cohort=cohort,
                n_mc=len(rr), n_cells=int(meta.n_cells.iloc[rr].sum()), pre_gene_eligible=(route, library, cohort) in groups,
                candidate_MC_ids=encode(meta.metacell_id.iloc[rr])))
        full, _ = features_for(z, sid, 'supported50')
        cols = ['atlas_pas_id', 'feature_index_1based', 'primary_gene_id', 'primary_gene_symbol', 'chrom', 'strand',
                'cut_site', 'position', 'coordinate_eligible', 'matched_source_count_50nt']
        small = full[cols].copy(); small['study_id'] = sid; axes.append(small)
        for setting, versions in [('supported50', list(VERSIONS)), ('expanded_primary', ['core_third_strict'])]:
            _, retained = features_for(z, sid, setting)
            retained_counts = retained.groupby('primary_gene_id').size().to_dict()
            for gene, part in full.groupby('primary_gene_id', sort=True):
                n_retained = retained_counts.get(gene, 0)
                background.append(dict(study_id=sid, setting=setting, gene_id=gene, n_full_primary_PAS=len(part),
                    n_retained_PAS=n_retained, coordinate_eligible=bool(part.coordinate_eligible.all()),
                    status='measurable_axis_requalification_required' if n_retained >= 2 else 'not_measurable',
                    reason='' if n_retained >= 2 else 'unusable_primary_axis' if not part.coordinate_eligible.all() else 'fewer_than_two_retained_PAS'))
            contexts = []
            for version in versions:
                for d in sampling[sid]:
                    for contrast, (target, reference) in CONTRASTS.items():
                        tg = version + '_' + target if target in ['Pstar', 'Mstar'] else target
                        rg = version + '_' + reference if reference in ['Pstar', 'Mstar'] else reference
                        a = groups.get((d['route'], d['library_key'], tg)); b = groups.get((d['route'], d['library_key'], rg))
                        allowed = contrast == 'D11' or (d['route'], d['library_key']) in oldkeys
                        why = '' if allowed and a is not None and b is not None else 'original_same_library_design_absent' if not allowed else 'cohort_absent_or_below_pre_gene_3MC_30cells'
                        row = dict(d, species=species, setting=setting, version=version, contrast=contrast,
                                   target_cohort=tg, reference_cohort=rg, status='pre_gene_eligible' if not why else 'not_estimable', reason=why)
                        for side, co in [('a', tg), ('b', rg)]:
                            rr = allgroups.get((d['route'], d['library_key'], co), np.array([], int))
                            row['pre_gene_n_mc_' + side] = len(rr); row['pre_gene_n_cells_' + side] = int(meta.n_cells.iloc[rr].sum())
                        library_eligibility.append(row)
                        if not why:
                            contexts.append(row)
            # Profile cohorts only needed by pre-qualified comparisons/gradient; qualification failures remain in separate matrices.
            wanted = {(d['route'], d['library_key'], co) for d in contexts for co in [d['target_cohort'], d['reference_cohort']]}
            gradcontexts = [q for q in gradient_class if q['study_id'] == sid and q['status'] == 'pre_gene_eligible'] if setting == 'supported50' else []
            for q in gradcontexts:
                wanted.update((q['route'], q['library_key'], 'gradient_' + q['panel'] + '_' + layer) for layer in ['low', 'middle', 'high'])
                wanted.update((q['route'], q['library_key'], 'original_' + node) for node in ['P', 'M'])
            for gene, sites, counts in gene_matrices(z, sid, setting):
                states = {key: state(counts, groups[key], meta) for key in sorted(wanted)}
                pasids = encode(sites.atlas_pas_id); position = sites.position.to_numpy()
                for key, s in states.items():
                    route, library, co = key
                    why = 'insufficient_gene_qualified_metacells' if s['n_mc'] < 3 else 'insufficient_gene_qualified_cells' if s['n_cells'] < 30 else ''
                    profile.write(dict(species=species, study_id=sid, route=route, library_key=library, setting=setting,
                        gene_id=gene, gene_symbol=sites.primary_gene_symbol.iloc[0], cohort=co,
                        group_id='|'.join([sid, route, library, setting, gene, co]), status='available' if not why else 'not_measurable', reason=why,
                        n_pre_gene_mc=len(groups[key]), n_pre_gene_cells=int(meta.n_cells.iloc[groups[key]].sum()),
                        n_mc=s['n_mc'], n_cells=s['n_cells'], qualified_MC_ids=encode(meta.metacell_id.iloc[s['rows']]),
                        qualified_MC_pas_depths=encode(counts[s['rows']].sum(1)), PAS_ids=pasids, position=encode(position),
                        pooled_molecules=encode(s['molecules']), nonzero_MC=encode(s['nonzero_mc']),
                        usage=encode(s['usage']) if s['n_mc'] else 'NA', raw_counts_source='inputs.npz'))
                for d in contexts:
                    a = states[(d['route'], d['library_key'], d['target_cohort'])]
                    b = states[(d['route'], d['library_key'], d['reference_cohort'])]
                    why = reason(a, b)
                    joint = (a['molecules'] + b['molecules'] >= 10) & (a['nonzero_mc'] + b['nonzero_mc'] >= 2)
                    eligibility.write(dict(d, gene_id=gene, status='not_measurable' if why else 'eligible', reason=why,
                        n_mc_a=a['n_mc'], n_mc_b=b['n_mc'], n_cells_a=a['n_cells'], n_cells_b=b['n_cells'], n_joint_PAS=int(joint.sum()),
                        target_group='|'.join([sid, d['route'], d['library_key'], setting, gene, d['target_cohort']]),
                        reference_group='|'.join([sid, d['route'], d['library_key'], setting, gene, d['reference_cohort']])))
                    if why:
                        continue
                    eligible_pairs.append(dict(d, gene_id=gene, gene_symbol=sites.primary_gene_symbol.iloc[0], status='descriptive', reason='',
                        PAS_ids=pasids, position=position, target_usage=a['usage'], reference_usage=b['usage'],
                        delta_DWUI=float((a['usage'] - b['usage']) @ position), n_joint_PAS=int(joint.sum()),
                        target_group='|'.join([sid, d['route'], d['library_key'], setting, gene, d['target_cohort']]),
                        reference_group='|'.join([sid, d['route'], d['library_key'], setting, gene, d['reference_cohort']]),
                        n_mc_a=a['n_mc'], n_mc_b=b['n_mc'], n_cells_a=a['n_cells'], n_cells_b=b['n_cells'],
                        full_weight=np.nan, common_weight=np.nan, common_four=False))
                for q in gradcontexts:
                    route, library, panel = q['route'], q['library_key'], q['panel']
                    p = states[(route, library, 'original_P')]
                    m = {layer: states[(route, library, 'gradient_' + panel + '_' + layer)] for layer in ['low', 'middle', 'high']}
                    m['whole'] = states[(route, library, 'original_M')]
                    reasons = {layer: reason(p, m[layer]) for layer in ['low', 'middle', 'high']}
                    joints = {layer: int(((p['molecules'] + m[layer]['molecules'] >= 10) & (p['nonzero_mc'] + m[layer]['nonzero_mc'] >= 2)).sum()) for layer in reasons}
                    why = ';'.join(layer + ':' + value for layer, value in reasons.items() if value)
                    gradient_eligibility.write(dict(study_id=sid, route=route, library_key=library, panel=panel, gene_id=gene,
                        status='not_measurable' if why else 'eligible', reason=why, n_mc_P=p['n_mc'], n_cells_P=p['n_cells'],
                        **{'n_mc_' + layer: m[layer]['n_mc'] for layer in reasons},
                        **{'n_cells_' + layer: m[layer]['n_cells'] for layer in reasons},
                        **{'n_joint_' + layer: value for layer, value in joints.items()}))
                    if why:
                        continue
                    d = next(row for row in sampling[sid] if row['route'] == route and row['library_key'] == library)
                    delta = {layer: float((p['usage']-m[layer]['usage']) @ position) for layer in m}
                    observed = delta['high'] - delta['low']; direct = float((m['low']['usage']-m['high']['usage']) @ position)
                    assert abs(observed-direct) < 1e-12
                    gradient_pairs.append(dict(d, panel=panel, gene_id=gene, PAS_ids=pasids, position=position,
                        status='descriptive', P_usage=p['usage'], **{'M_' + layer + '_usage': value['usage'] for layer, value in m.items()},
                        **{'delta_' + layer: value for layer, value in delta.items()}, difference_high_minus_low=observed,
                        direct_M_low_minus_high=direct, identity_abs_error=abs(observed-direct), monotonicity=monotonicity(delta),
                        P_group='|'.join([sid, route, library, setting, gene, 'original_P']),
                        **{'M_' + layer + '_group': '|'.join([sid, route, library, setting, gene, 'gradient_' + panel + '_' + layer]) for layer in reasons},
                        weight=np.nan))
            print('NEW_EFFECTS', sid, setting, 'profiles', profile.n, 'eligibility', eligibility.n,
                  'measurable_pairs', len(eligible_pairs), 'gradient_pairs', len(gradient_pairs), flush=True)
    profile.close(); eligibility.close(); gradient_eligibility.close()
    save(library_eligibility, 'library_apa_eligibility.tsv')
    save(cohort_coverage, 'cohort_coverage.tsv')
    save(background, 'gene_background.tsv.gz')
    save(pd.concat(axes, ignore_index=True), 'primary_axis_features.tsv.gz')
    finish_effects(pd.DataFrame(eligible_pairs), pd.DataFrame(gradient_pairs), z, candidate_libraries,
                   dict(classification=len(classifications), profiles=profile.n, eligibility=eligibility.n, gradient_eligibility=gradient_eligibility.n))


def finish_effects(pairs, gradient_pairs, z, candidate_libraries, counts):
    allroutes, route_summary, library_summary, identities = [], [], [], []
    for setting, versions in [('supported50', list(VERSIONS)), ('expanded_primary', ['core_third_strict'])]:
        for version in versions:
            print('AGGREGATE_BEGIN', setting, version, flush=True)
            q = pairs[pairs.setting.eq(setting) & pairs.version.eq(version)].copy()
            full = aggregate_new(q, common=False); common = aggregate_new(q, common=True)
            pairs.loc[q.index, ['full_weight', 'common_weight', 'common_four']] = q[['full_weight', 'common_weight', 'common_four']]
            for label, effect in [('full', full), ('common_four', common)]:
                effect['version'] = version; effect['setting'] = setting; allroutes.append(effect)
                route_summary.extend(summaries(effect, label, version, setting))
                chosen = q if label == 'full' else q[q.common_four]
                library_summary.extend(summaries(chosen, label, version, setting, 'library', candidate_libraries))
            identities.extend(bridge_identities(q, common, version, setting))
            print('AGGREGATE_COMPLETE', setting, version, len(full), len(common), flush=True)
    routes = pd.concat(allroutes, ignore_index=True)
    if not len(gradient_pairs):
        gradient_pairs = pd.DataFrame(columns=['study_id', 'route', 'library_key', 'panel', 'gene_id', 'delta_low', 'delta_middle', 'delta_high', 'delta_whole', 'difference_high_minus_low', 'monotonicity'])
    gradient_routes = gradient_aggregate(gradient_pairs)
    for name, data in [('pair_effects.tsv.gz', pairs), ('gradient_pair_effects.tsv.gz', gradient_pairs)]:
        out = data.copy()
        for col in ['position', 'target_usage', 'reference_usage', 'P_usage'] + ['M_' + layer + '_usage' for layer in ['low', 'middle', 'high', 'whole']]:
            if col in out:
                out[col] = out[col].map(encode)
        save(out, name)
    save(routes, 'route_effects.tsv.gz'); save(route_summary, 'route_summary.tsv')
    save(library_summary, 'library_summary.tsv')
    save(pd.DataFrame(identities, columns=['study_id', 'route', 'library_key', 'gene_id', 'level', 'version', 'setting',
                                          'observed_change', 'reference_change', 'target_change', 'abs_error']), 'bridge_identities.tsv.gz')
    save(gradient_routes, 'gradient_route_effects.tsv.gz')
    save(gradient_summary(gradient_routes), 'gradient_route_summary.tsv')
    save(gradient_summary(gradient_pairs, 'library', candidate_libraries), 'gradient_library_summary.tsv')
    gident = []
    for level, source in [('library', gradient_pairs), ('route', gradient_routes)]:
        for row in source.to_dict('records'):
            gident.append({k: row[k] for k in ['study_id', 'route', 'panel', 'gene_id', 'difference_high_minus_low', 'direct_M_low_minus_high', 'identity_abs_error']} | {'level': level, 'library_key': row.get('library_key', '')})
    save(pd.DataFrame(gident, columns=['study_id', 'route', 'library_key', 'panel', 'gene_id', 'level',
                                     'difference_high_minus_low', 'direct_M_low_minus_high', 'identity_abs_error']), 'gradient_identities.tsv.gz')
    audit = frame(z, 'audit_json'); save(audit, 'input_coverage.tsv')
    (HERE / 'remote_input_sha256.json').write_text(str(z['provenance_json']) + '\n')
    (HERE / 'analysis_validation.json').write_text(json.dumps(dict(status='COMPLETE_AWAITING_INDEPENDENT_QA',
        protocol_sha256=PROTOCOL_SHA, input_hashes=INPUTS, n_classification_rows=counts['classification'],
        n_group_profile_rows=counts['profiles'], n_gene_pair_qualification_rows=counts['eligibility'],
        n_gradient_gene_qualification_rows=counts['gradient_eligibility'], n_measurable_pairs=len(pairs),
        n_gradient_pairs=len(gradient_pairs), n_bridge_identities=len(identities),
        max_bridge_identity_abs_error=max((q['abs_error'] for q in identities), default=0.),
        max_gradient_identity_abs_error=max((q['identity_abs_error'] for q in gident), default=0.)), indent=2) + '\n')


def recover_staged(z):
    """Resume aggregation from sealed qualification/profile streams without recalculating PAS."""
    names = ['gene_group_profiles.tsv.gz', 'pair_eligibility.tsv.gz', 'gradient_gene_eligibility.tsv.gz']
    source_hashes = {name: digest(HERE / name) for name in names}
    parts, n_eligibility = [], 0
    for chunk in read(HERE / 'pair_eligibility.tsv.gz', chunksize=250000):
        n_eligibility += len(chunk)
        parts.append(chunk[chunk.status.eq('eligible')])
    eligible = pd.concat(parts, ignore_index=True)
    grad = read(HERE / 'gradient_gene_eligibility.tsv.gz')
    n_gradient = len(grad); grad = grad[grad.status.eq('eligible')]
    needed = set(eligible.target_group) | set(eligible.reference_group)
    for row in grad.itertuples(index=False):
        prefix = '|'.join([row.study_id, row.route, row.library_key, 'supported50', row.gene_id])
        needed.update(prefix + '|' + co for co in ['original_P', 'original_M'] +
                      ['gradient_' + row.panel + '_' + layer for layer in ['low', 'middle', 'high']])
    profiles, n_profiles = {}, 0
    columns = ['group_id', 'gene_symbol', 'status', 'n_mc', 'n_cells', 'PAS_ids', 'position', 'usage']
    for chunk in read(HERE / 'gene_group_profiles.tsv.gz', usecols=columns, chunksize=250000):
        n_profiles += len(chunk)
        for row in chunk[chunk.group_id.isin(needed)].to_dict('records'):
            assert row['group_id'] not in profiles and row['status'] == 'available'
            row['position'] = np.array(json.loads(row['position']))
            row['usage'] = np.array(json.loads(row['usage']))
            assert np.isfinite(row['usage']).all() and abs(row['usage'].sum()-1) < 1e-12
            profiles[row['group_id']] = row
    assert set(profiles) == needed, ('sealed_profiles_missing', len(needed-set(profiles)))
    context_cols = ['study_id', 'route', 'library_key', 'setting', 'version', 'contrast']
    context = read(HERE / 'library_apa_eligibility.tsv').fillna('')
    contexts = {tuple(row[k] for k in context_cols): row for row in context.to_dict('records')}
    assert len(contexts) == len(context)
    pairs = []
    for row in eligible.to_dict('records'):
        a, b = [profiles[row[side + '_group']] for side in ['target', 'reference']]
        assert a['PAS_ids'] == b['PAS_ids']
        np.testing.assert_array_equal(a['position'], b['position'])
        assert (a['n_mc'], b['n_mc'], a['n_cells'], b['n_cells']) == tuple(row[k] for k in ['n_mc_a', 'n_mc_b', 'n_cells_a', 'n_cells_b'])
        d = dict(contexts[tuple(row[k] for k in context_cols)], **row)
        d.update(gene_symbol=a['gene_symbol'], status='descriptive', reason='', PAS_ids=a['PAS_ids'],
            position=a['position'], target_usage=a['usage'], reference_usage=b['usage'],
            delta_DWUI=float((a['usage']-b['usage']) @ a['position']),
            full_weight=np.nan, common_weight=np.nan, common_four=False)
        pairs.append(d)
    pairframe = pd.DataFrame(pairs)
    reference = read(HERE / 'baseline_pair_effects.tsv.gz')
    key = ['study_id', 'route', 'gene_id', 'library_key']
    replay = pairframe[pairframe.version.eq('core_third_strict') & pairframe.setting.eq('supported50') & pairframe.contrast.eq('D00')]
    oldkeys = reference.assign(library_key=reference.library_key_a)
    check = replay[key + ['delta_DWUI']].merge(oldkeys[key + ['delta_DWUI']], on=key, suffixes=('_recovery', '_accepted'), validate='one_to_one', how='outer')
    assert len(check) == len(replay) == len(reference) and check.notna().all().all()
    baseline_error = float(abs(check.delta_DWUI_recovery-check.delta_DWUI_accepted).max())
    assert baseline_error < 1e-12
    sampling = {(sid, row['route'], row['library_key']): row for sid in STUDIES
                for row in sampling_pairs(z, sid, metadata(z, sid))}
    gradient_pairs = []
    for row in grad.to_dict('records'):
        prefix = '|'.join([row['study_id'], row['route'], row['library_key'], 'supported50', row['gene_id']])
        p = profiles[prefix + '|original_P']
        m = {layer: profiles[prefix + '|gradient_' + row['panel'] + '_' + layer] for layer in ['low', 'middle', 'high']}
        m['whole'] = profiles[prefix + '|original_M']
        for layer in m:
            assert p['PAS_ids'] == m[layer]['PAS_ids']
            np.testing.assert_array_equal(p['position'], m[layer]['position'])
        delta = {layer: float((p['usage']-value['usage']) @ p['position']) for layer, value in m.items()}
        observed = delta['high']-delta['low']; direct = float((m['low']['usage']-m['high']['usage']) @ p['position'])
        assert abs(observed-direct) < 1e-12
        gradient_pairs.append(dict(sampling[(row['study_id'], row['route'], row['library_key'])],
            panel=row['panel'], gene_id=row['gene_id'], PAS_ids=p['PAS_ids'], position=p['position'],
            status='descriptive', P_usage=p['usage'], **{'M_' + layer + '_usage': value['usage'] for layer, value in m.items()},
            **{'delta_' + layer: value for layer, value in delta.items()}, difference_high_minus_low=observed,
            direct_M_low_minus_high=direct, identity_abs_error=abs(observed-direct), monotonicity=monotonicity(delta),
            P_group=prefix+'|original_P', **{'M_' + layer + '_group': prefix+'|gradient_'+row['panel']+'_'+layer for layer in ['low', 'middle', 'high']}, weight=np.nan))
    meta = read(HERE / 'metacell_metadata.tsv.gz', usecols=['study_id', 'route', 'library_key', 'in_planned_route'])
    candidates = {(sid, route): sorted(meta.loc[meta.study_id.eq(sid) & meta.route.eq(route) & meta.in_planned_route, 'library_key'].unique())
                  for sid, (_, routes) in STUDIES.items() for route in routes}
    classification = read(HERE / 'metacell_classification.tsv.gz', usecols=['metacell_id'])
    for name, expected in source_hashes.items():
        assert digest(HERE / name) == expected, ('sealed_source_changed', name)
    (HERE / 'recovery_validation.json').write_text(json.dumps(dict(status='PASS', reason='route_summary_must_not_require_library_key',
        source_hashes=source_hashes, n_recovered_pairs=len(pairframe), n_recovered_gradient_pairs=len(gradient_pairs),
        n_baseline_pairs_rechecked=len(check), max_baseline_abs_error=baseline_error), indent=2)+'\n')
    print('STAGED_RECOVERY_COMPLETE', len(pairframe), len(gradient_pairs), flush=True)
    finish_effects(pairframe, pd.DataFrame(gradient_pairs), z, candidates,
                   dict(classification=len(classification), profiles=n_profiles, eligibility=n_eligibility, gradient_eligibility=n_gradient))


def aggregate_bridge(pairs, background, contrasts):
    """Each two-way diagnostic background is matched independently of the other bridge."""
    effect_rows, weight_rows, changes = [], [], []
    weight_cache = {}
    for (sid, route, gene), allrows in pairs.groupby(['study_id', 'route', 'gene_id'], sort=True):
        q = allrows[allrows.contrast.isin(contrasts)]
        counts = q.groupby('library_key').contrast.nunique()
        shared = set(counts[counts.eq(len(contrasts))].index)
        q = q[q.library_key.isin(shared)]
        if not len(q):
            continue
        by_contrast = {}
        for contrast in contrasts:
            part = q[q.contrast.eq(contrast)].sort_values('library_key')
            key = tuple(part[['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b']].itertuples(index=False, name=None))
            if key not in weight_cache:
                weight_cache[key] = weights(part)
            ww = weight_cache[key]
            a = ww @ np.stack(part.target_usage); b = ww @ np.stack(part.reference_usage)
            d = float((a-b) @ part.position.iloc[0])
            assert abs(d - ww @ part.delta_DWUI.to_numpy()) < 1e-12
            rec = dict(study_id=sid, route=route, gene_id=gene, contrast=contrast, background=background,
                n_libraries=len(part), contributor_libraries=encode(part.library_key), weights=encode(ww),
                target_usage=encode(a), reference_usage=encode(b), PAS_ids=part.PAS_ids.iloc[0], delta_DWUI=d)
            effect_rows.append(rec); by_contrast[contrast] = (a, b, d, part.position.iloc[0])
            for row, weight in zip(part.to_dict('records'), ww):
                weight_rows.append(dict(study_id=sid, route=route, gene_id=gene, library_key=row['library_key'],
                    contrast=contrast, background=background, shared_weight=float(weight), delta_DWUI=row['delta_DWUI'],
                    time_a=row['time_a'], time_b=row['time_b'], sample_unit_a=row['sample_unit_a'], sample_unit_b=row['sample_unit_b']))
        for level, library, rows in [('route', '', by_contrast)] + [
                ('library', library, {r['contrast']: (r['target_usage'], r['reference_usage'], r['delta_DWUI'], r['position'])
                    for r in part.to_dict('records')}) for library, part in q.groupby('library_key')]:
            a0, b0, d0, axis = rows['D00']; a1, b1, d1, _ = rows[contrasts[1]]
            if background == 'common_M_bridge':
                np.testing.assert_allclose(a0, a1, atol=1e-12, rtol=0)
                direct = float((b0-b1) @ axis)
            else:
                np.testing.assert_allclose(b0, b1, atol=1e-12, rtol=0)
                direct = float((a1-a0) @ axis)
            change = d1-d0; error = abs(change-direct)
            assert error < 1e-12
            changes.append(dict(study_id=sid, route=route, gene_id=gene, library_key=library, background=background,
                                level=level, changed_contrast=contrasts[1], delta_change=float(change), direct_change=direct, identity_abs_error=error))
    effect_columns = ['study_id', 'route', 'gene_id', 'contrast', 'background', 'n_libraries', 'contributor_libraries',
                      'weights', 'target_usage', 'reference_usage', 'PAS_ids', 'delta_DWUI']
    return pd.DataFrame(effect_rows, columns=effect_columns), weight_rows, changes


def addendum_effects():
    assert digest(HERE / 'ADDENDUM.zh.md') == ADDENDUM_SHA
    INPUTS[str(HERE / 'ADDENDUM.zh.md')] = ADDENDUM_SHA
    pairs = read(HERE / 'pair_effects.tsv.gz')
    for column in ['position', 'target_usage', 'reference_usage']:
        pairs[column] = pairs[column].map(lambda value: np.array(json.loads(value)))
    assert not pairs.duplicated(['study_id', 'route', 'gene_id', 'library_key', 'version', 'setting', 'contrast']).any()
    meta = read(HERE / 'metacell_metadata.tsv.gz', usecols=['study_id', 'route', 'library_key', 'in_planned_route'])
    candidates = {(sid, route): sorted(meta.loc[meta.study_id.eq(sid) & meta.route.eq(route) & meta.in_planned_route, 'library_key'].unique())
                  for sid, (_, routes) in STUDIES.items() for route in routes}
    allroutes, allweights, allchanges, route_summaries, library_summaries = [], [], [], [], []
    for setting, versions in [('supported50', list(VERSIONS)), ('expanded_primary', ['core_third_strict'])]:
        for version in versions:
            print('ADDENDUM_BEGIN', setting, version, flush=True)
            q = pairs[pairs.setting.eq(setting) & pairs.version.eq(version)]
            for bg, contrasts in [('common_M_bridge', ['D00', 'D01']), ('common_P_bridge', ['D00', 'D10'])]:
                effects, rows, changes = aggregate_bridge(q, bg, contrasts)
                effects['setting'] = setting; effects['version'] = version; allroutes.append(effects)
                w = pd.DataFrame(rows, columns=['study_id', 'route', 'gene_id', 'library_key', 'contrast', 'background',
                                               'shared_weight', 'delta_DWUI', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b'])
                w['setting'] = setting; w['version'] = version; allweights.append(w)
                for row in changes:
                    row.update(setting=setting, version=version)
                allchanges.extend(changes)
                route_summaries.extend(row for row in summaries(effects, bg, version, setting) if row['contrast'] in contrasts)
                library_summaries.extend(row for row in summaries(w, bg, version, setting, 'library', candidates) if row['contrast'] in contrasts)
            print('ADDENDUM_COMPLETE', setting, version, flush=True)
    save(pd.concat(allroutes, ignore_index=True), 'bridge_route_effects.tsv.gz')
    save(pd.concat(allweights, ignore_index=True), 'bridge_pair_weights.tsv.gz')
    save(route_summaries, 'bridge_route_summary.tsv'); save(library_summaries, 'bridge_library_summary.tsv')
    save(pd.DataFrame(allchanges, columns=['study_id', 'route', 'gene_id', 'library_key', 'background', 'level',
        'changed_contrast', 'delta_change', 'direct_change', 'identity_abs_error', 'setting', 'version']), 'bridge_changes.tsv.gz')
    (HERE / 'addendum_validation.json').write_text(json.dumps(dict(status='COMPLETE_AWAITING_INDEPENDENT_QA',
        protocol_sha256=PROTOCOL_SHA, addendum_sha256=ADDENDUM_SHA, input_hashes=INPUTS,
        n_bridge_gene_route_rows=sum(len(q) for q in allroutes), n_bridge_pair_weight_rows=sum(len(q) for q in allweights),
        n_change_identities=len(allchanges), max_identity_abs_error=max((q['identity_abs_error'] for q in allchanges), default=0.)), indent=2) + '\n')
    print('ADDENDUM_EFFECTS_COMPLETE', flush=True)


def remove_failed_support_pairs(pairs, lost):
    keys = ['study_id', 'route', 'library_key', 'setting', 'version', 'contrast', 'gene_id']
    assert not lost.duplicated(keys).any() and not pairs.duplicated(keys).any()
    marked = pairs.merge(lost.assign(_remove=True), on=keys, how='left', validate='one_to_one')
    assert marked._remove.eq(True).sum() == len(lost)
    kept = marked.loc[marked._remove.ne(True)].drop(columns='_remove').copy()
    kept[['full_weight', 'common_weight']] = np.nan
    kept['common_four'] = False
    return kept


def unique_support_effects():
    """Fixed QA sensitivity: discard exact failed unique-source gates and rematch every background."""
    assert digest(HERE / 'QA-ADDENDUM.zh.md') == QA_ADDENDUM_SHA
    INPUTS[str(HERE / 'QA-ADDENDUM.zh.md')] = QA_ADDENDUM_SHA
    source = ROOT / 'notes/2026-10-10-marker-reference-review/shared_MC_pairs_lost.tsv'
    assert digest(source) == LOST_KEYS_SHA
    lost = read(source)
    pairs = read(HERE / 'pair_effects.tsv.gz')
    assert len(lost) == 15 and lost.contrast.eq('D01').all()
    pairs = remove_failed_support_pairs(pairs, lost)
    for column in ['position', 'target_usage', 'reference_usage']:
        pairs[column] = pairs[column].map(lambda value: np.array(json.loads(value)))
    meta = read(HERE / 'metacell_metadata.tsv.gz', usecols=['study_id', 'route', 'library_key', 'in_planned_route'])
    candidates = {(sid, route): sorted(meta.loc[meta.study_id.eq(sid) & meta.route.eq(route) & meta.in_planned_route, 'library_key'].unique())
                  for sid, (_, routes) in STUDIES.items() for route in routes}
    allroutes, allbridges, allweights, rs, ls, brs, bls, four_ids, change_ids = [], [], [], [], [], [], [], [], []
    weight_cols = ['study_id', 'route', 'library_key', 'gene_id', 'version', 'setting', 'contrast',
                   'delta_DWUI', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b']
    for setting, versions in [('supported50', list(VERSIONS)), ('expanded_primary', ['core_third_strict'])]:
        for version in versions:
            print('UNIQUE_SUPPORT_BEGIN', setting, version, flush=True)
            q = pairs[pairs.setting.eq(setting) & pairs.version.eq(version)].copy()
            full, common = aggregate_new(q, False), aggregate_new(q, True)
            for background, effect, selected, weight_column in [('full', full, q, 'full_weight'),
                    ('common_four', common, q[q.common_four], 'common_weight')]:
                effect['version'] = version; effect['setting'] = setting; allroutes.append(effect)
                rs.extend(summaries(effect, background, version, setting))
                ls.extend(summaries(selected, background, version, setting, 'library', candidates))
                w = selected[weight_cols].copy()
                w['background'] = background; w['weight'] = selected[weight_column]
                assert np.isfinite(w.weight).all()
                allweights.append(w)
            four_ids.extend(bridge_identities(q, common, version, setting))
            for background, contrasts in [('common_M_bridge', ['D00', 'D01']), ('common_P_bridge', ['D00', 'D10'])]:
                effect, rows, changes = aggregate_bridge(q, background, contrasts)
                effect['version'] = version; effect['setting'] = setting; allbridges.append(effect)
                w = pd.DataFrame(rows, columns=['study_id', 'route', 'gene_id', 'library_key', 'contrast', 'background',
                    'shared_weight', 'delta_DWUI', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b']).rename(columns={'shared_weight': 'weight'})
                w['version'] = version; w['setting'] = setting; allweights.append(w)
                brs.extend(row for row in summaries(effect, background, version, setting) if row['contrast'] in contrasts)
                bls.extend(row for row in summaries(w, background, version, setting, 'library', candidates) if row['contrast'] in contrasts)
                for row in changes:
                    row.update(version=version, setting=setting)
                change_ids.extend(changes)
            print('UNIQUE_SUPPORT_COMPLETE', setting, version, flush=True)
    save(pd.concat(allroutes, ignore_index=True), 'unique_support_route_effects.tsv.gz')
    save(rs, 'unique_support_route_summary.tsv'); save(ls, 'unique_support_library_summary.tsv')
    save(pd.concat(allbridges, ignore_index=True), 'unique_support_bridge_route_effects.tsv.gz')
    save(brs, 'unique_support_bridge_route_summary.tsv'); save(bls, 'unique_support_bridge_library_summary.tsv')
    save(pd.concat(allweights, ignore_index=True), 'unique_support_pair_weights.tsv.gz')
    save(lost, 'unique_support_removed_pairs.tsv')
    save(pd.DataFrame(four_ids, columns=['study_id', 'route', 'library_key', 'gene_id', 'level', 'version', 'setting',
        'observed_change', 'reference_change', 'target_change', 'abs_error']), 'unique_support_four_identities.tsv.gz')
    save(pd.DataFrame(change_ids, columns=['study_id', 'route', 'gene_id', 'library_key', 'background', 'level',
        'changed_contrast', 'delta_change', 'direct_change', 'identity_abs_error', 'setting', 'version']), 'unique_support_bridge_changes.tsv.gz')
    diagnostic = ROOT / 'notes/2026-10-10-marker-reference-review/shared_MC_joint_support.tsv.gz'
    INPUTS[str(diagnostic)] = digest(diagnostic)
    (HERE / 'unique_support_validation.json').write_text(json.dumps(dict(status='COMPLETE_AWAITING_INDEPENDENT_QA',
        qa_addendum_sha256=QA_ADDENDUM_SHA, lost_keys_sha256=LOST_KEYS_SHA, input_hashes=INPUTS,
        n_removed_qualified_pair_rows=len(lost), n_distinct_removed_gene_library=len(lost[['study_id', 'gene_id', 'library_key']].drop_duplicates()),
        n_remaining_qualified_pair_rows=len(pairs), n_four_identities=len(four_ids), n_bridge_identities=len(change_ids),
        max_four_identity_abs_error=max((row['abs_error'] for row in four_ids), default=0.),
        max_bridge_identity_abs_error=max((row['identity_abs_error'] for row in change_ids), default=0.),
        gradient='unchanged_original_P_and_M_are_disjoint_no_copy'), indent=2)+'\n')
    print('UNIQUE_SUPPORT_ALL_COMPLETE', flush=True)


def accepted_baseline():
    acceptance = json.loads((HERE / 'baseline_validation.json').read_text())
    assert acceptance['status'] == 'PASS' and acceptance['protocol_sha256'] == PROTOCOL_SHA
    for path, expected in acceptance['baseline_source_hashes'].items():
        assert digest(path) == expected, ('accepted_baseline_input_changed', path)
        INPUTS[path] = expected
    return acceptance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--new-only', action='store_true')
    parser.add_argument('--addendum-only', action='store_true')
    parser.add_argument('--aggregate-only', action='store_true')
    parser.add_argument('--unique-support-only', action='store_true')
    args = parser.parse_args()
    assert digest(HERE / 'PROTOCOL.zh.md') == PROTOCOL_SHA
    INPUTS[str(HERE / 'PROTOCOL.zh.md')] = PROTOCOL_SHA
    INPUTS[str(HERE / 'inputs.npz')] = digest(HERE / 'inputs.npz')
    if args.addendum_only:
        accepted_baseline()
        addendum_effects()
        return
    if args.unique_support_only:
        accepted_baseline()
        unique_support_effects()
        return
    z = np.load(HERE / 'inputs.npz', allow_pickle=False)
    if args.aggregate_only:
        accepted_baseline()
        recover_staged(z)
        return
    if args.new_only:
        accepted_baseline()
    else:
        pairs, routes = baseline(z)
    if args.baseline_only:
        return
    new_effects(z)


if __name__ == '__main__':
    main()
