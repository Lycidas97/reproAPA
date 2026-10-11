"""Frozen-design APA arithmetic shared by the seven-study comparison."""
import numpy as np
import pandas as pd

GROUP = ['species', 'study_id', 'route', 'gene_id', 'gene_symbol', 'contrast', 'design']
PAIR = ['library_key_a', 'library_key_b', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b']


def hierarchical_weights(frame):
    q = frame.reset_index(drop=True).copy()
    q['_row'] = np.arange(len(q))
    levels = [GROUP + ['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'],
        GROUP + ['time_a', 'sample_unit_a', 'time_b'],
        GROUP + ['time_a', 'sample_unit_a'], GROUP + ['time_a'], GROUP]
    q['_weight'] = 1 / q.groupby(levels[0], dropna=False).gene_id.transform('size')
    for child, parent in zip(levels[:-1], levels[1:]):
        n = q[child].drop_duplicates().groupby(parent, dropna=False).size().rename('_n').reset_index()
        q = q.merge(n, on=parent, validate='many_to_one')
        q['_weight'] /= q.pop('_n')
    assert np.allclose(q.groupby(GROUP, dropna=False)._weight.sum(), 1, atol=1e-12)
    return pd.Series(q.sort_values('_row')._weight.to_numpy(), index=frame.index)


def state_counts(counts, meta, route, library, state, common_mc=None):
    assert counts.shape[0] == len(meta)
    selected = meta.route.eq(route) & meta.library_key.eq(library) & meta.node.eq(state) & meta.purity.ge(.8)
    if common_mc is not None:
        selected &= meta.metacell_id.isin(common_mc)
    rows = np.flatnonzero(selected)
    data = counts[rows].toarray()
    assert np.issubdtype(data.dtype, np.integer) and (data >= 0).all()
    good = data.sum(1) >= 20
    data = data[good]; rows = rows[good]
    usage = (data / data.sum(1)[:, None]).mean(0) if len(rows) else np.zeros(counts.shape[1])
    return dict(metacell_ids=meta.metacell_id.iloc[rows].to_numpy(),
        n_mc=len(rows), n_cells=int(meta.n_cells.iloc[rows].sum()), usage=usage,
        molecules=data.sum(0), nonzero_mc=(data > 0).sum(0))


def pair_effect(a, b, sites):
    if a['n_mc'] < 3 or b['n_mc'] < 3 or a['n_cells'] < 30 or b['n_cells'] < 30:
        return None, 'below_20_molecules_per_MC_or_3_MC_or_30_cells_per_state'
    joint = (a['molecules'] + b['molecules'] >= 10) & (a['nonzero_mc'] + b['nonzero_mc'] >= 2)
    signed = sites.cut_site.to_numpy() * np.where(sites.strand.eq('+'), 1, -1)
    if np.unique(signed[joint]).size < 2 or sites.chrom.nunique() != 1 or sites.strand.nunique() != 1 or np.ptp(signed) == 0:
        return None, 'fewer_than_two_joint_supported_distinct_sites_or_multiple_chromosomes_strands'
    delta = a['usage'] - b['usage']
    assert abs(delta.sum()) < 1e-12
    position = sites.primary_axis_position.to_numpy()
    return dict(n_mc_a=a['n_mc'], n_mc_b=b['n_mc'], n_cells_a=a['n_cells'], n_cells_b=b['n_cells'],
        pas_counts_a=int(a['molecules'].sum()), pas_counts_b=int(b['molecules'].sum()),
        n_sites=len(sites), n_joint_sites=int(joint.sum()),
        delta_primary_axis=float(delta @ position) if np.isfinite(position).all() else np.nan,
        delta_own_axis=float(delta @ ((signed - signed.min()) / np.ptp(signed))),
        delta_nt=float(delta @ signed), delta_usage=delta,
        MPRO=float(np.ptp(delta)), inside_primary_axis=bool(np.isfinite(position).all() and
            ((position >= 0) & (position <= 1)).all())), None


def matched_pair_views(tool, repro, sites, meta, design):
    roles = [(design['library_key_a'], 'P'), (design['library_key_b'], design['contrast'][-1])]
    t = [state_counts(tool, meta, design['route'], lib, state) for lib, state in roles]
    r = [state_counts(repro, meta, design['route'], lib, state) for lib, state in roles]
    common = [np.intersect1d(a['metacell_ids'], b['metacell_ids']) for a, b in zip(t, r)]
    result = {}
    for view in ['matched_sites_independent_metacells', 'matched_sites_common_metacells']:
        if view.endswith('common_metacells'):
            t = [state_counts(tool, meta, design['route'], lib, state, ids)
                for (lib, state), ids in zip(roles, common)]
            r = [state_counts(repro, meta, design['route'], lib, state, ids)
                for (lib, state), ids in zip(roles, common)]
        te, tr = pair_effect(*t, sites); re, rr = pair_effect(*r, sites)
        result[view] = dict(tool=te, repro=re, reason_tool=tr, reason_repro=rr,
            common_mc_a=len(common[0]), common_mc_b=len(common[1]))
    return result


def fixed_original_pairs(available, frozen):
    assert not available.duplicated(PAIR).any() and not frozen.duplicated(PAIR).any()
    if frozen.empty:
        return None, dict(reason='no_frozen_primary_gene_comparison', original_pairs=0)
    selected = frozen[PAIR + ['weight']].rename(columns={'weight': 'frozen_weight'}).merge(
        available, on=PAIR, how='left', indicator=True, validate='one_to_one')
    missing = int(selected._merge.ne('both').sum())
    audit = dict(original_pairs=len(frozen), available_pairs=len(available),
        missing_original_pairs=missing, extra_available_pairs=len(available) - len(frozen) + missing)
    if missing:
        return None, dict(audit, reason='incomplete_frozen_library_pair_coverage')
    assert abs(selected.frozen_weight.sum() - 1) < 1e-12
    return selected.drop(columns='_merge'), dict(audit, reason=None)


def collapse_pair_effects(frame, weights):
    weights = np.asarray(weights)
    assert len(frame) == len(weights) and abs(weights.sum() - 1) < 1e-12
    delta = weights @ np.stack(frame.delta_usage)
    result = dict(delta_usage=delta, MPRO=float(np.ptp(delta)))
    for name in ['delta_primary_axis', 'delta_own_axis', 'delta_nt']:
        result[name] = float(weights @ frame[name])
    if 'repro_delta_usage' in frame:
        reference = weights @ np.stack(frame.repro_delta_usage)
        result.update(repro_delta_usage=reference, repro_MPRO=float(np.ptp(reference)),
            usage_change_MAE=float(np.abs(delta - reference).mean()))
        for name in ['repro_delta_primary_axis', 'repro_delta_own_axis', 'repro_delta_nt']:
            result[name] = float(weights @ frame[name])
    return result
