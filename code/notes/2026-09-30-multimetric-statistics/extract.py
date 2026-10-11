"""Frozen-contributor APA metric comparisons; historical WPP is DWUI externally."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
PROJECTS = Path('/Users/lycidas/work/科研/projects')
MOUSE = PROJECTS / '2026-09-02_mouse-time-aware-integration-pilot'
HUMAN = PROJECTS / '2026-09-09_human-time-aware-integration'
STATE = MOUSE / 'outputs/figure3-apa-gex-state-analysis-2026-09-19'
PROFILES = MOUSE / 'data/apa-majority-shared-v1/profiles'
BASE = ['species', 'study_id', 'route', 'gene_id', 'gene_symbol', 'contrast', 'design']
STATEKEY = ['route', 'library_key', 'state', 'gene_id']
METRICS = ['WPP', 'RWUI', 'PDUI', 'PPUI']
INPUTS = {}


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read(path, **kwargs):
    path = Path(path)
    INPUTS[str(path)] = {'sha256': digest(path), 'bytes': path.stat().st_size}
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, name):
    path = HERE / name
    frame.to_csv(path, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                 compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)


def pair_weights(frame):
    """Same equal source-unit/date then target-unit/date hierarchy as manuscript."""
    q = frame.copy()
    levels = [BASE + ['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'],
              BASE + ['time_a', 'sample_unit_a', 'time_b'],
              BASE + ['time_a', 'sample_unit_a'], BASE + ['time_a'], BASE]
    q['weight'] = 1 / q.groupby(levels[0], dropna=False).gene_id.transform('size')
    for child, parent in zip(levels[:-1], levels[1:]):
        distinct = q[child].drop_duplicates()
        n = distinct.groupby(parent, dropna=False).size().rename('_n').reset_index()
        q = q.merge(n, on=parent, validate='many_to_one')
        q['weight'] /= q.pop('_n')
    assert np.allclose(q.groupby(BASE).weight.sum(), 1, atol=1e-12)
    return q


def select_pairs():
    q = read(STATE / 'matched_library_pair_effects.tsv.gz', low_memory=False)
    q = q[q.setting.eq('supported50_p08') &
          (((q.contrast.eq('P_S')) & q.design.eq('true_source')) |
           ((q.contrast.eq('P_M')) & q.design.eq('same_library')))].copy()
    q = q[~q.route.str.contains('reference|perturbation')].copy()
    assert not q.duplicated(BASE + ['library_key_a', 'library_key_b']).any()
    return q


def raw_states(folder, pairs):
    wanted = []
    for suffix in ['a', 'b']:
        x = pairs[['route', 'library_key_' + suffix, 'gene_id', 'contrast']].copy()
        x['state'] = 'P' if suffix == 'a' else x.contrast.str[-1]
        x = x.rename(columns={'library_key_' + suffix: 'library_key'})
        wanted.append(x[STATEKEY])
    wanted = pd.concat(wanted).drop_duplicates()
    path = MOUSE / 'analysis/figure3-apa-gex-state-analysis-2026-09-19/raw' / folder / 'apa_gex_states.tsv.gz'
    INPUTS[str(path)] = {'sha256': digest(path), 'bytes': path.stat().st_size}
    columns = STATEKEY + ['setting', 'APA', 'GEX_CPM', 'GEX_all', 'n_mc', 'n_cells', 'pas_counts']
    pieces = []
    for part in pd.read_csv(path, sep='\t', usecols=columns, chunksize=200000):
        part = part[part.setting.eq('supported50_p08')]
        pieces.append(part.merge(wanted, on=STATEKEY, validate='one_to_one'))
    data = pd.concat(pieces, ignore_index=True)
    assert not data.duplicated(STATEKEY).any()
    assert len(data) == len(wanted), (folder, 'missing frozen states', len(data), len(wanted))
    return data


def study_vectors(folder, pairs):
    species = pairs.species.iloc[0]
    src = (MOUSE / 'data/metacell-pas-v7-v1/studies' if species == 'mouse'
           else HUMAN / 'data/results_metacell_pas_v2/studies') / folder
    support = (MOUSE / 'data/pas-annotation-support-v1/studies' if species == 'mouse'
               else MOUSE / 'data/apa-majority-shared-v1/human_annotation') / folder
    meta = read(PROFILES / folder / 'selected_metacells.tsv.gz').sort_values('matrix_row')
    orig = read(src / 'metacells.tsv.gz')
    assert np.array_equal(orig.metacell_id.iloc[meta.matrix_row], meta.metacell_id)
    features = read(src / 'features.tsv.gz')
    features = features[features.keep_for_apa_primary & features.primary_gene_id.notna()].copy()
    shape = features.groupby('primary_gene_id').agg(n=('cut_site', 'nunique'),
               c=('chrom', 'nunique'), s=('strand', 'nunique'))
    features = features[features.primary_gene_id.isin(shape.query('n>=2 and c==1 and s==1').index)].copy()
    x = features.cut_site * np.where(features.strand.eq('+'), 1, -1)
    features['position'] = (x - x.groupby(features.primary_gene_id).transform('min')) / (
        x.groupby(features.primary_gene_id).transform('max') - x.groupby(features.primary_gene_id).transform('min'))
    features = features.merge(read(support / 'primary_pas_reference_matches.tsv.gz',
              usecols=['atlas_pas_id', 'matched_source_count_50nt']), on='atlas_pas_id', validate='one_to_one')
    features = features[features.matched_source_count_50nt.ge(1)].copy()
    features = features[features.primary_gene_id.isin(pairs.gene_id.unique())].reset_index(drop=True)
    assert features.atlas_pas_id.is_unique and features.feature_index_1based.is_unique
    signed = features.cut_site * np.where(features.strand.eq('+'), 1, -1)
    features['rank'] = signed.groupby(features.primary_gene_id).rank(method='first')
    features['n'] = features.groupby('primary_gene_id').atlas_pas_id.transform('size')
    features['RWUI'] = features['rank'] / features['n']
    features['PDUI'] = (features['rank'] == features['n']).astype(float)
    features['PPUI'] = features['rank'].eq(1).astype(float)
    features['WPP'] = features.position
    raw = raw_states(folder, pairs)
    countpath = src / 'metacell_by_pas_counts.npz'
    INPUTS[str(countpath)] = {'sha256': digest(countpath), 'bytes': countpath.stat().st_size}
    counts = sparse.load_npz(countpath).tocsr()
    assert np.issubdtype(counts.data.dtype, np.integer) and (counts.data >= 0).all()
    counts = counts[meta.matrix_row.to_numpy()][:, features.feature_index_1based.to_numpy() - 1].tocsr()
    meta = meta.reset_index(drop=True)
    genes = pd.Index(features.primary_gene_id.unique())
    gi = genes.get_indexer(features.primary_gene_id)
    mapping = sparse.csr_matrix((np.ones(len(features)), (np.arange(len(features)), gi)),
                                shape=(len(features), len(genes)))
    groups = {g: np.asarray(ix) for g, ix in features.groupby('primary_gene_id', sort=False).groups.items()}
    vectors, state_rows, audits = {}, [], []
    for (route, lib, state), old in raw.groupby(['route', 'library_key', 'state'], sort=False):
        rows = np.flatnonzero(meta.route.eq(route) & meta.library_key.eq(lib) & meta.node.eq(state) & meta.purity.ge(.8))
        c = counts[rows].tocsr()
        totals = (c @ mapping).toarray()
        good = totals >= 20
        n_mc = good.sum(0)
        n_cells = (good * meta.iloc[rows].n_cells.to_numpy()[:, None]).sum(0)
        coo = c.tocoo()
        take = good[coo.row, gi[coo.col]]
        sums = np.bincount(coo.col[take], weights=coo.data[take] /
                           totals[coo.row[take], gi[coo.col[take]]], minlength=len(features))
        site_count = np.bincount(coo.col[take], weights=coo.data[take], minlength=len(features))
        site_mc = np.bincount(coo.col[take], minlength=len(features))
        means = np.divide(sums, n_mc[gi], out=np.zeros(len(features)), where=n_mc[gi] > 0)
        for row in old.itertuples():
            ix = groups[row.gene_id]
            gindex = genes.get_loc(row.gene_id)
            usage = means[ix]
            assert np.isclose(usage.sum(), 1, atol=1e-12)
            assert n_mc[gindex] == row.n_mc and n_cells[gindex] == row.n_cells
            assert site_count[ix].sum() == row.pas_counts
            estimated = float(usage @ features.WPP.iloc[ix].to_numpy())
            assert abs(estimated - row.APA) < 1e-8, (folder, route, lib, state, row.gene_id, estimated, row.APA)
            key = (route, lib, state, row.gene_id)
            vectors[key] = (usage, site_count[ix], site_mc[ix], row.GEX_CPM, row.pas_counts / row.n_mc)
            rec = dict(species=species, study_id=pairs.study_id.iloc[0], route=route,
                       library_key=lib, state=state, gene_id=row.gene_id,
                       n_mc=row.n_mc, n_cells=row.n_cells, pas_counts=row.pas_counts)
            for metric in METRICS:
                rec[metric] = float(usage @ features[metric].iloc[ix].to_numpy())
            state_rows.append(rec)
            audits.append(abs(estimated - row.APA))
        print('STATE', folder, route, lib, state, len(old), flush=True)
    axis = features[['primary_gene_id', 'atlas_pas_id', 'cut_site', 'strand', 'position', 'rank', 'n']].copy()
    axis = axis.assign(species=species, study_id=pairs.study_id.iloc[0])
    return vectors, pd.DataFrame(state_rows), features, groups, axis, max(audits)


def calculate(folder, pairs):
    vectors, states, features, groups, axis, state_error = study_vectors(folder, pairs)
    pairs = pair_weights(pairs)
    added = []
    for row in pairs.itertuples():
        a = vectors[(row.route, row.library_key_a, 'P', row.gene_id)]
        b = vectors[(row.route, row.library_key_b, row.contrast[-1], row.gene_id)]
        ix = groups[row.gene_id]
        assert int(((a[1] + b[1] >= 10) & (a[2] + b[2] >= 2)).sum()) == row.n_joint_PAS
        rec = {'mean_GEX': (a[3] + b[3]) / 2, 'APA_coverage': (a[4] + b[4]) / 2}
        for metric in METRICS:
            w = features[metric].iloc[ix].to_numpy()
            rec['target_' + metric] = a[0] @ w
            rec['reference_' + metric] = b[0] @ w
            rec['delta_' + metric] = (a[0] - b[0]) @ w
        added.append(rec)
    pairs = pd.concat([pairs.reset_index(drop=True), pd.DataFrame(added)], axis=1)
    error = float((pairs.delta_WPP - pairs.APA_delta).abs().max())
    assert error < 1e-8, (folder, 'pair recovery', error)
    routes = []
    for key, part in pairs.groupby(BASE, sort=False):
        rec = dict(zip(BASE, key))
        usage_a, usage_b = [], []
        for row in part.itertuples():
            usage_a.append(vectors[(row.route, row.library_key_a, 'P', row.gene_id)][0])
            usage_b.append(vectors[(row.route, row.library_key_b, row.contrast[-1], row.gene_id)][0])
        weights = part.weight.to_numpy()
        ua = weights @ np.stack(usage_a)
        ub = weights @ np.stack(usage_b)
        difference = ua - ub
        ix = groups[rec['gene_id']]
        for metric in METRICS:
            w = features[metric].iloc[ix].to_numpy()
            rec['target_' + metric] = ua @ w
            rec['reference_' + metric] = ub @ w
            rec['delta_' + metric] = difference @ w
        rec.update(MPRO=float(difference.max() - difference.min()),
                   PAS_number=len(ix), supported_span=float(features.cut_site.iloc[ix].max() - features.cut_site.iloc[ix].min()),
                   rwui_min_weight=float(features.RWUI.iloc[ix].min()),
                   rwui_max_weight=float(features.RWUI.iloc[ix].max()),
                   reference_concentration=float(ub @ ub), reference_min_usage=float(ub.min()),
                   reference_max_usage=float(ub.max()), n_pairs=len(part),
                   n_target_libraries=part.library_key_a.nunique(), n_source_libraries=part.library_key_b.nunique(),
                   n_target_dates=part.time_a.nunique(),
                   mean_GEX=float(weights @ part.mean_GEX), APA_coverage=float(weights @ part.APA_coverage),
                   target_usage=json.dumps(ua.tolist()), reference_usage=json.dumps(ub.tolist()),
                   delta_usage=json.dumps(difference.tolist()),
                   PAS_ids=json.dumps(features.atlas_pas_id.iloc[ix].tolist()))
        rec['delta_RWUI_range_scaled'] = rec['delta_RWUI'] / (1 - 1 / len(ix))
        routes.append(rec)
    return pd.DataFrame(routes), pairs, states, axis, dict(study_id=pairs.study_id.iloc[0],
        n_pairs=len(pairs), n_routes=len(routes), n_states=len(states), max_state_error=state_error,
        max_pair_error=error)


def main():
    frozen = select_pairs()
    route_parts, pair_parts, state_parts, axes, audits = [], [], [], [], []
    for meta_path in sorted(PROFILES.glob('*/selected_metacells.tsv.gz')):
        identity = pd.read_csv(meta_path, sep='\t', nrows=1)
        pairs = frozen[frozen.study_id.eq(identity.study_id.iloc[0])].copy()
        if pairs.empty:
            continue
        r, p, s, a, audit = calculate(meta_path.parent.name, pairs)
        route_parts.append(r); pair_parts.append(p); state_parts.append(s); axes.append(a); audits.append(audit)
        save(r, meta_path.parent.name + '_route_effects.tsv.gz')
        print('STUDY_COMPLETE', audit, flush=True)
    routes = pd.concat(route_parts, ignore_index=True)
    pairs = pd.concat(pair_parts, ignore_index=True)
    expected = read(STATE / 'matched_study_effects.tsv.gz')
    expected = expected[expected.setting.eq('supported50_p08') &
        (((expected.contrast.eq('P_S')) & expected.design.eq('true_source')) |
         ((expected.contrast.eq('P_M')) & expected.design.eq('same_library')))]
    expected = expected[~expected.route.str.contains('reference|perturbation')]
    compare = routes.merge(expected[BASE + ['APA_delta']], on=BASE, how='outer',
                           validate='one_to_one', indicator=True)
    assert compare._merge.eq('both').all(), compare._merge.value_counts().to_dict()
    route_error = float((compare.delta_WPP - compare.APA_delta).abs().max())
    assert route_error < 1e-8, route_error
    save(routes, 'route_effects.tsv.gz')
    save(pairs, 'pair_effects.tsv.gz')
    save(pd.concat(state_parts, ignore_index=True), 'library_state_metrics.tsv.gz')
    save(pd.concat(axes, ignore_index=True), 'frozen_PAS_axes.tsv.gz')
    INPUTS[str(Path(__file__))] = {'sha256': digest(__file__)}
    (HERE / 'input_manifest.json').write_text(json.dumps(INPUTS, indent=2) + '\n')
    (HERE / 'extraction_audit.json').write_text(json.dumps(dict(status='PASS',
        historical_name_mapping={'WPP': 'DWUI'}, studies=audits, pair_records=len(pairs),
        route_records=len(routes), maximum_route_error=route_error), indent=2) + '\n')
    print('COMPLETE', len(routes), len(pairs), route_error, flush=True)


if __name__ == '__main__':
    main()
