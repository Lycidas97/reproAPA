"""Frozen-feature, exact-metacell state support audit; no manuscript mutations."""
from pathlib import Path
import gc
import hashlib
import io
import json
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACK = ROOT / 'analysis_code/recurrence_correction_2026-09-25/sources/MOUSE'
FROZEN = PACK / 'outputs/figure3-recurrence-common-reference-2026-09-25'
STATES = PACK / 'analysis/figure3-apa-gex-state-analysis-2026-09-19'
M = Path('/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot')
H = Path('/Users/lycidas/work/科研/projects/2026-09-09_human-time-aware-integration')
SPECS = {
    'GSE103221': ('mouse', 'CIP', 'gse103221_pluripotency_reprogramming_core'),
    'GSE115943': ('mouse', 'WOT', 'gse115943_pluripotency_reprogramming_core'),
    'GSE118258': ('human', 'GSE118258', 'gse118258_human_pluripotency_reprogramming'),
    'GSE178325': ('human', '178325_0618', 'gse178325_human_chemical_pluripotency_reprogramming'),
}
SOURCE = {'D0', 'source', 'Stage_I_day_0'}
REQUIRED = {'mouse': ['R', 'E'], 'human': ['R']}
AXES = {'mouse': 'GSE103221__GSE115943', 'human': 'GSE118258__GSE178325'}
VARIANTS = ['baseline', 'support_q95', 'support_q90']
CANDIDATE = 'ENSMUSG00000020091'
INPUTS = {}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def record(path):
    INPUTS[str(path)] = digest(path)
    return path


def read(path, **kwargs):
    return pd.read_csv(record(path), sep='\t', **kwargs)


def write(frame, filename):
    frame.to_csv(HERE / filename, sep='\t', index=False, na_rep='NA', float_format='%.12g')


def protected():
    paths = [p for f in ('main', 'supplement') for p in (ROOT / f).rglob('*') if p.is_file()]
    paths += [ROOT / f for f in ('Fig/Figure.3.pdf', 'Fig/SuppFig.S10.pdf', 'output/pdf/main.pdf', 'output/pdf/supplement.pdf')]
    return {str(p.relative_to(ROOT)): digest(p) for p in paths}


def nearest_support(dist_cross, dist_self, labels, other_labels, libraries, k=15):
    """Cross-study neighbors compete over ALL labels; caliper is within-study."""
    n = len(labels)
    order = np.argsort(dist_cross, axis=1, kind='stable')[:, :k]
    out = pd.DataFrame(index=np.arange(n))
    for label in sorted(set(other_labels)):
        out['vote_' + label] = (other_labels[order] == label).mean(axis=1)
    out['same_label_vote'] = (other_labels[order] == labels[:, None]).mean(axis=1)
    out['predicted_label'] = out.filter(like='vote_').idxmax(axis=1).str.removeprefix('vote_')
    out['within_mean5_distance'] = np.nan
    out['cross_same_mean5_distance'] = np.nan
    out['calibration'] = 'insufficient_neighbors'
    out['caliper_q95'] = np.nan
    out['caliper_q90'] = np.nan
    for label in sorted(set(labels)):
        ix = np.flatnonzero(labels == label)
        oi = np.flatnonzero(other_labels == label)
        one_library = len(set(libraries[ix])) == 1
        within = dist_self[np.ix_(ix, ix)].copy()
        forbidden = np.eye(len(ix), dtype=bool) if one_library else libraries[ix, None] == libraries[ix][None, :]
        within[forbidden] = np.inf
        v = np.sort(within, axis=1)[:, :5].mean(1) if len(ix) >= 6 else np.full(len(ix), np.nan)
        v[~np.isfinite(v)] = np.nan
        out.loc[ix, 'within_mean5_distance'] = v
        if np.isfinite(v).sum() < 5 or len(oi) < 5 or dist_cross.shape[1] < k:
            continue
        out.loc[ix, 'calibration'] = 'same_library_only' if one_library else 'cross_library'
        out.loc[ix, 'caliper_q95'] = np.nanquantile(v, .95)
        out.loc[ix, 'caliper_q90'] = np.nanquantile(v, .90)
        out.loc[ix, 'cross_same_mean5_distance'] = np.sort(dist_cross[np.ix_(ix, oi)], axis=1)[:, :5].mean(1)
    for q in (95, 90):
        out[f'support_q{q}'] = (out.same_label_vote >= .6) & (out.cross_same_mean5_distance <= out[f'caliper_q{q}'])
    return out


def expression():
    archive = np.load(record(HERE / 'expression_inputs.npz'), allow_pickle=False)
    (HERE / 'remote_input_hashes.json').write_text(str(archive['provenance_json']) + '\n')
    metas, ranks, programs = {}, {}, []
    for study, (species, route, folder) in SPECS.items():
        local = read(M / 'data/apa-majority-shared-v1/profiles' / folder / 'selected_metacells.tsv.gz')
        raw = pd.DataFrame(json.loads(str(archive[study + '_meta_json'])))
        assert local.metacell_id.tolist() == raw.metacell_id.tolist()
        for key in ['n_cells', 'library_key', 'route', 'node', 'matrix_row']:
            assert np.array_equal(local[key], raw[key]), (study, key)
        assert np.allclose(local.purity, raw.purity, atol=1e-9)
        panel_root = M if species == 'mouse' else H
        panel = read(panel_root / 'data/state-matching-evidence-v1/inputs_snapshot/features.tsv')
        panel = panel[panel.panel.eq('new_nonhvg')]
        assert panel.gene_id.tolist() == [x['gene_id'] for x in json.loads(str(archive[species + '_features_json']))]
        write(panel, species + '_heldout_features.tsv')
        keep = local.route.eq(route) & local.purity.ge(.8) & (~local.node.eq('S') | local.experimental_time.isin(SOURCE))
        meta = local[keep].copy().reset_index(drop=True)
        source_dir = (PACK / 'data/metacell-pas-v7-v1/studies' if species == 'mouse' else H / 'data/results_metacell_pas_v2/studies') / folder
        full = read(source_dir / 'metacells.tsv.gz')
        assert np.array_equal(full.metacell_id.iloc[meta.matrix_row], meta.metacell_id)
        column = 'local_state_id' if species == 'mouse' else 'original_state_id'
        meta['local_state_id'] = full.iloc[meta.matrix_row][column].to_numpy()
        meta['gex_umi'] = archive[study + '_total_gex'][keep]
        meta['panel_detected_genes'] = (archive[study + '_counts'][keep] > 0).sum(1)
        ranks[study] = rankdata(archive[study + '_counts'][keep], axis=1, method='average')
        p = json.loads(str(archive[species + '_program_json']))
        symbols = archive[species + '_marker_symbols'].tolist()
        marker = archive[study + '_marker_logcp10k'][keep]
        for name, genes in p.items():
            scores = marker[:, [symbols.index(g.upper()) for g in genes.split()]].mean(1)
            programs.append(pd.DataFrame({'study_id': study, 'metacell_id': meta.metacell_id, 'program': name, 'score': scores}))
        metas[study] = meta
    output = []
    for species in ['mouse', 'human']:
        a, b = [s for s in SPECS if SPECS[s][0] == species]
        cross = cdist(ranks[a], ranks[b], metric='correlation')
        for study, other, distance in [(a, b, cross), (b, a, cross.T)]:
            meta = metas[study]
            support = nearest_support(distance, cdist(ranks[study], ranks[study], metric='correlation'),
                                      meta.node.to_numpy(), metas[other].node.to_numpy(), meta.library_key.to_numpy())
            output.append(pd.concat([meta, support], axis=1))
    data = pd.concat(output, ignore_index=True)
    data['baseline'] = True
    write(data, 'metacell_expression_support.tsv.gz')
    program = pd.concat(programs, ignore_index=True).merge(data[['metacell_id', 'species', 'node', 'library_key', 'local_state_id', 'support_q95', 'support_q90']], on='metacell_id', validate='many_to_one')
    write(program, 'metacell_programs.tsv.gz')
    summary = []
    for keys, part in data.groupby(['species', 'study_id', 'node'], sort=False):
        row = dict(zip(['species', 'study_id', 'state'], keys))
        row.update(n_mc=len(part), n_cells=int(part.n_cells.sum()), n_libraries=part.library_key.nunique(),
                   median_same_label_vote=part.same_label_vote.median(), calibration=';'.join(part.calibration.unique()),
                   median_within_distance=part.within_mean5_distance.median(), median_cross_distance=part.cross_same_mean5_distance.median(),
                   caliper_q95=part.caliper_q95.iloc[0], caliper_q90=part.caliper_q90.iloc[0])
        for variant in VARIANTS[1:]:
            selected = part[part[variant]]
            row[variant + '_mc'] = len(selected)
            row[variant + '_fraction'] = len(selected) / len(part)
            row[variant + '_cells'] = int(selected.n_cells.sum())
            row[variant + '_libraries'] = selected.library_key.nunique()
        summary.append(row)
    write(pd.DataFrame(summary), 'expression_support_summary.tsv')
    write(data.groupby(['species', 'study_id', 'node', 'local_state_id', 'library_key'], dropna=False).agg(n_mc=('metacell_id', 'size'), n_cells=('n_cells', 'sum'), support_q95_mc=('support_q95', 'sum'), support_q90_mc=('support_q90', 'sum'), mean_same_label_vote=('same_label_vote', 'mean')).reset_index(), 'local_state_library_support.tsv')
    write(program.groupby(['species', 'study_id', 'node', 'library_key', 'program']).score.mean().reset_index(), 'library_programs.tsv')
    return data


def frozen_cohorts():
    classes = read(FROZEN / 'A_same_series_process_classes.tsv.gz')
    primary = classes[classes.analysis.eq('primary_0.02_0.02')].copy()
    keep = primary.groupby(['species', 'gene_id']).apply(lambda x: x.pair_measurable.all() and x.endpoint_retained.all(), include_groups=False)
    keys = set(keep[keep].index)
    primary = primary[[(s, g) in keys for s, g in zip(primary.species, primary.gene_id)]]
    assert primary.groupby('species').gene_id.nunique().to_dict() == {'human': 165, 'mouse': 295}
    return primary


def reconstruct(data, primary):
    regions = read(STATES / 'shared/PAS_regions.tsv.gz')
    sampling = read(STATES / 'raw/library_sampling_units.tsv').set_index('library_key')
    rows, pas_rows, coverage = [], [], []
    genes = primary[['species', 'gene_id', 'gene_symbol']].drop_duplicates()
    for study, (species, route, folder) in SPECS.items():
        meta = data[data.study_id.eq(study) & data.node.isin(list('SREP'))].reset_index(drop=True)
        source = (PACK / 'data/metacell-pas-v7-v1/studies' if species == 'mouse' else H / 'data/results_metacell_pas_v2/studies') / folder
        features = read(source / 'features.tsv.gz', usecols=['feature_index_1based', 'atlas_pas_id'])
        index = features.set_index('atlas_pas_id').feature_index_1based - 1
        mapping = regions[regions.axis_id.eq(AXES[species]) & regions.study_id.eq(study) & regions.gene_id.isin(genes[genes.species.eq(species)].gene_id)]
        columns = np.sort(index.loc[mapping.atlas_pas_id.drop_duplicates()].to_numpy(int))
        raw_to_local = {int(v): i for i, v in enumerate(columns)}
        matrix = sparse.load_npz(record(source / 'metacell_by_pas_counts.npz')).tocsr()[meta.matrix_row.to_numpy()][:, columns]
        for gene, block in mapping.groupby('gene_id', sort=False):
            reg = sorted(block.region.unique())
            if len(reg) < 2:
                continue
            counts = np.column_stack([np.asarray(matrix[:, [raw_to_local[int(c)] for c in index.loc[block[block.region.eq(r)].atlas_pas_id]]].sum(1)).ravel() for r in reg])
            axes = np.array([block[block.region.eq(r)].axis.iloc[0] for r in reg])
            total = counts.sum(1)
            fractions = np.divide(counts, total[:, None], out=np.full_like(counts, np.nan, dtype=float), where=total[:, None] > 0)
            wpp = fractions @ axes
            eligible = total >= 20
            symbol = genes[genes.gene_id.eq(gene)].gene_symbol.iloc[0]
            for variant in VARIANTS:
                keep = eligible & meta[variant].to_numpy(bool)
                for (library, state), ix in meta[keep].groupby(['library_key', 'node'], sort=False).groups.items():
                    ix = np.asarray(list(ix))
                    if len(ix) < 3 or meta.iloc[ix].n_cells.sum() < 30:
                        continue
                    row = dict(variant=variant, species=species, study_id=study, route=route, gene_id=gene, gene_symbol=symbol,
                               library_key=library, state=state, time=str(meta.iloc[ix].experimental_time.iloc[0]),
                               sample_unit=str(sampling.loc[library, 'sampling_unit']), APA=float(wpp[ix].mean()),
                               n_mc=len(ix), n_cells=int(meta.iloc[ix].n_cells.sum()), pas_counts=int(total[ix].sum()))
                    rows.append(row)
                    if gene == CANDIDATE:
                        for j, region in enumerate(reg):
                            pas_rows.append(dict(row, region=region, axis=axes[j], fraction=float(fractions[ix, j].mean())))
                if gene == CANDIDATE:
                    for state, ix in meta.groupby('node').groups.items():
                        ix = np.asarray(list(ix))
                        coverage.append(dict(variant=variant, study_id=study, state=state,
                                             expression_pool_mc=len(ix), pas_measurable_mc=int(eligible[ix].sum()),
                                             selected_pas_measurable_mc=int(keep[ix].sum())))
        del matrix
        gc.collect()
        print('PAS rebuilt', study, flush=True)
    result = pd.DataFrame(rows)
    write(result, 'library_state_WPP.tsv.gz')
    write(pd.DataFrame(pas_rows), 'Eif4ebp2_library_PAS.tsv')
    write(pd.DataFrame(coverage), 'Eif4ebp2_metacell_coverage.tsv')
    original = read(STATES / 'shared/state_values.tsv.gz')
    original = original[original.axis_id.isin(AXES.values()) & original.setting.eq('shared_p08') & original.state.isin(list('SREP'))]
    original = original[[SPECS.get(s, ('', 'not', ''))[1] == r for s, r in zip(original.study_id, original.route)]]
    original = original[~original.state.eq('S') | original.time.isin(SOURCE)]
    original = original.merge(genes[['species', 'gene_id']], on=['species', 'gene_id'], validate='many_to_one')
    keys = ['species', 'study_id', 'route', 'gene_id', 'library_key', 'state', 'time']
    check = original.merge(result[result.variant.eq('baseline')], on=keys, how='outer', suffixes=('_frozen', '_rebuilt'), indicator=True, validate='one_to_one')
    audit = dict(rows=len(check), key_mismatches=int(check._merge.ne('both').sum()), max_abs_WPP_error=float((check.APA_frozen - check.APA_rebuilt).abs().max()),
                 mc_count_mismatches=int(check.n_mc_frozen.ne(check.n_mc_rebuilt).sum()))
    (HERE / 'baseline_validation.json').write_text(json.dumps(audit, indent=2) + '\n')
    assert audit['key_mismatches'] == audit['mc_count_mismatches'] == 0 and audit['max_abs_WPP_error'] < 1e-9, audit
    return result


def weighted_window(group, value='APA'):
    if group.empty:
        return np.nan
    assert not group.library_key.duplicated().any()
    assert np.isfinite(group[value]).all()
    return group.groupby(['time', 'sample_unit'])[value].mean().groupby(level=0).mean().mean()


def pair_class(a, b, sign):
    values = np.array([a, b]) * sign
    if not np.isfinite(values).all():
        return 'unavailable'
    cats = np.where(values >= .02, 'endpoint_aligned', np.where(values <= -.02, 'reverse', 'within_threshold'))
    if cats[0] == cats[1]:
        return cats[0]
    return 'direction_difference' if set(cats) == {'endpoint_aligned', 'reverse'} else 'mixed_or_threshold_undetermined'


def effects_and_outcomes(library, primary):
    rows = []
    for keys, part in library.groupby(['variant', 'species', 'study_id', 'route', 'gene_id', 'gene_symbol'], sort=False):
        s = weighted_window(part[part.state.eq('S') & part.time.isin(SOURCE)])
        for state, at_state in part.groupby('state', sort=False):
            value = weighted_window(at_state)
            rows.append(dict(zip(['variant', 'species', 'study_id', 'route', 'gene_id', 'gene_symbol'], keys), state=state,
                             source_WPP=s, target_WPP=value, delta_WPP=value-s, n_libraries=at_state.library_key.nunique(),
                             n_mc=int(at_state.n_mc.sum()), n_cells=int(at_state.n_cells.sum())))
    effects = pd.DataFrame(rows)
    write(effects, 'state_effects.tsv')
    write(effects[effects.gene_id.eq(CANDIDATE)], 'Eif4ebp2_state_effects.tsv')
    assert effects.groupby(['variant', 'species', 'study_id', 'gene_id']).source_WPP.nunique().max() == 1
    lookup = effects.set_index(['variant', 'study_id', 'gene_id', 'state']).delta_WPP.to_dict()
    output = []
    for p in primary.itertuples():
        for variant in VARIANTS:
            delta = [lookup.get((variant, st, p.gene_id, p.state), np.nan) for st in (p.study_1, p.study_2)]
            end = [lookup.get((variant, st, p.gene_id, 'P'), np.nan) for st in (p.study_1, p.study_2)]
            measurable = bool(np.isfinite(delta + end).all())
            passed = measurable and bool((np.array(end) * p.direction_sign >= .02).all())
            current = pair_class(*delta, p.direction_sign) if measurable else 'unavailable'
            outcome = 'unavailable' if not measurable else ('endpoint_no_longer_passes' if not passed else ('stable' if current == p.pair_class else 'category_changed'))
            output.append(dict(variant=variant, species=p.species, gene_id=p.gene_id, gene_symbol=p.gene_symbol, state=p.state,
                               primary_class=p.pair_class, primary_opposite=p.pair_class == 'direction_difference', current_class=current,
                               measurable=measurable, endpoint_pass=passed, outcome=outcome,
                               delta_1=delta[0], delta_2=delta[1], endpoint_delta_1=end[0], endpoint_delta_2=end[1]))
            if variant == 'baseline':
                assert outcome == 'stable', (p, current, delta, end)
    output = pd.DataFrame(output)
    write(output, 'fixed_gene_state_outcomes.tsv')
    gene = []
    for keys, part in output.groupby(['variant', 'species', 'gene_id', 'gene_symbol'], sort=False):
        measurable = bool(part.measurable.all())
        passed = bool(part.endpoint_pass.all())
        original_opposite = bool(part.primary_opposite.any())
        current_opposite = bool(part.current_class.eq('direction_difference').any())
        status = 'unavailable' if not measurable else ('endpoint_no_longer_passes' if not passed else ('stable' if original_opposite == current_opposite else 'category_changed'))
        gene.append(dict(zip(['variant', 'species', 'gene_id', 'gene_symbol'], keys), primary_opposite=original_opposite,
                         current_opposite=current_opposite, outcome=status, exact_configuration_unchanged=bool(part.primary_class.eq(part.current_class).all())))
    gene = pd.DataFrame(gene)
    write(gene, 'fixed_complete_gene_outcomes.tsv')
    summary = []
    for scope, frame in [('complete_gene', gene), ('gene_state', output)]:
        for (variant, species), block in frame.groupby(['variant', 'species']):
            for cohort, part in [('all_original', block), ('original_opposite', block[block.primary_opposite])]:
                counts = part.outcome.value_counts()
                summary.append(dict(scope=scope, variant=variant, species=species, cohort=cohort, fixed_n=len(part),
                                    **{name: int(counts.get(name, 0)) for name in ['stable', 'category_changed', 'endpoint_no_longer_passes', 'unavailable']}))
    write(pd.DataFrame(summary), 'fixed_cohort_summary.tsv')
    return effects


def main():
    before = protected()
    data = expression()
    print('Expression support frozen; proceeding to PAS.', flush=True)
    primary = frozen_cohorts()
    write(primary, 'frozen_original_cohorts.tsv')
    library = reconstruct(data, primary)
    effects_and_outcomes(library, primary)
    after = protected()
    assert before == after
    (HERE / 'verification.json').write_text(json.dumps({'protected_files_unchanged': before == after, 'protected_hashes': after, 'local_inputs_sha256': INPUTS}, indent=2) + '\n')
    print('Completed. Formal figures and manuscript unchanged.', flush=True)


if __name__ == '__main__':
    main()
