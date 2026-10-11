"""Evaluate all frozen-design gene backgrounds with auditable contributors."""
import numpy as np
import pandas as pd
from atlas_effect_math import GROUP, PAIR, state_counts, pair_effect, hierarchical_weights
from atlas_effect_math import fixed_original_pairs, collapse_pair_effects

VIEWS = ['native_terminal_sites', 'matched_sites_independent_metacells', 'matched_sites_common_metacells']
METRICS = ['delta_primary_axis', 'delta_own_axis', 'delta_nt', 'MPRO']


def effect_class(values, threshold):
    values = np.asarray(values, dtype=float)
    return np.where(np.isfinite(values), np.where(values <= -threshold, -1,
        np.where(values >= threshold, 1, 0)), np.nan)


def serializable(value):
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [serializable(item) for item in value]
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def gene_effects(gene, symbol, native_sites, native_counts, shared_sites, tool_shared_counts,
        repro_shared_counts, meta, design, frozen, primary_catalogue_gene):
    assert meta.metacell_id.is_unique
    template = design.copy(); template['gene_id'] = gene; template['gene_symbol'] = symbol
    contexts = dict(native_terminal_sites=(native_sites, native_counts, None),
        matched_sites_independent_metacells=(shared_sites, tool_shared_counts, repro_shared_counts),
        matched_sites_common_metacells=(shared_sites, tool_shared_counts, repro_shared_counts))
    cache = {}

    def states(view, role):
        key = (view, role)
        if key in cache:
            return cache[key]
        route, library, state = role
        sites, tool, repro = contexts[view]
        if view == 'matched_sites_common_metacells':
            a, b = states('matched_sites_independent_metacells', role)
            common = np.intersect1d(a['metacell_ids'], b['metacell_ids'])
            value = (state_counts(tool, meta, route, library, state, common),
                state_counts(repro, meta, route, library, state, common))
        else:
            value = (state_counts(tool, meta, route, library, state),
                state_counts(repro, meta, route, library, state) if repro is not None else None)
        cache[key] = value
        return value

    pair_records, route_records = [], []
    for (route, contrast, comparison_design), planned in template.groupby(['route', 'contrast', 'design'], sort=False):
        old = frozen[frozen.route.eq(route) & frozen.contrast.eq(contrast) & frozen.design.eq(comparison_design)]
        identity = planned.iloc[0][GROUP].to_dict()
        if len(old):
            assert abs(old.weight.sum() - 1) < 1e-12 and old.gene_id.eq(gene).all()
            identity['gene_symbol'] = old.gene_symbol.iloc[0]
        main = float(old.delta_WPP @ old.weight) if len(old) else np.nan
        background = 'frozen_primary_measurable' if len(old) else (
            'primary_catalogue_no_frozen_main_effect' if primary_catalogue_gene else 'native_terminal_gene_only')
        for view in VIEWS:
            sites, tool_matrix, repro_matrix = contexts[view]
            available = []
            tool_qualified, repro_qualified = 0, 0
            for planned_pair in planned.to_dict('records'):
                roles = [(route, planned_pair['library_key_a'], 'P'),
                    (route, planned_pair['library_key_b'], contrast[-1])]
                (a, ra), (b, rb) = [states(view, role) for role in roles]
                effect, reason = pair_effect(a, b, sites)
                reference, reference_reason = pair_effect(ra, rb, sites) if ra is not None else (None, None)
                tool_qualified += effect is not None
                repro_qualified += reference is not None
                pair_identity = dict(planned_pair, **{key: identity[key] for key in ['gene_id', 'gene_symbol']}, view=view)
                pair_records.append(dict(pair_identity, site_ids=sites.site_id.to_numpy(), tool=effect,
                    repro_same_PAS=reference, reason_tool=reason, reason_repro=reference_reason,
                    tool_state_a=a, tool_state_b=b, repro_state_a=ra, repro_state_b=rb))
                if effect is None or (view != 'native_terminal_sites' and reference is None):
                    continue
                value = dict(pair_identity, **effect)
                if reference is not None:
                    for key in METRICS + ['delta_usage']:
                        value['repro_' + key] = reference[key]
                value['metacell_ids_a'] = a['metacell_ids']; value['metacell_ids_b'] = b['metacell_ids']
                available.append(value)
            frame = pd.DataFrame(available) if available else pd.DataFrame(columns=GROUP + PAIR)
            fixed, audit = fixed_original_pairs(frame, old)
            for contributors in ['all_available', 'fixed_original_pairs_complete']:
                part = frame if contributors == 'all_available' else fixed
                missing = 'no_jointly_qualified_library_pair' if frame.empty else None
                if contributors == 'fixed_original_pairs_complete':
                    missing = audit['reason']
                row = dict(identity, view=view, contributors=contributors, background_class=background,
                    original_endpoint_design_pairs=len(planned), tool_qualified_pairs=int(tool_qualified),
                    repro_same_PAS_qualified_pairs=int(repro_qualified), jointly_qualified_pairs=len(frame),
                    frozen_main_pairs=len(old), frozen_main_delta_DWUI=main,
                    missing_reason=missing, evaluated_pairs=0 if part is None else len(part),
                    missing_original_pairs=audit.get('missing_original_pairs', 0),
                    extra_available_pairs=audit.get('extra_available_pairs', len(frame)), site_ids=sites.site_id.to_numpy())
                if view != 'native_terminal_sites':
                    row.update({'repro_' + name: np.nan for name in METRICS})
                if part is not None and len(part):
                    weights = hierarchical_weights(part) if contributors == 'all_available' else part.frozen_weight
                    row.update(collapse_pair_effects(part, weights))
                    row['inside_primary_axis'] = bool(part.inside_primary_axis.all())
                    row['contributing_library_pairs'] = part[PAIR].to_dict('records')
                    row['pair_weights'] = np.asarray(weights)
                    for side in ['a', 'b']:
                        ids = np.unique(np.concatenate(part['metacell_ids_' + side].to_list()))
                        selected = meta.metacell_id.isin(ids)
                        row['unique_tool_metacells_' + side] = len(ids)
                        row['unique_tool_cells_' + side] = int(meta.n_cells[selected].sum())
                else:
                    row.update({name: np.nan for name in METRICS})
                route_records.append(row)
    return pair_records, route_records


def agreement_summary(routes):
    rows = []
    keys = ['species', 'study_id', 'route', 'contrast', 'design', 'view', 'contributors']
    for identity, part in routes.groupby(keys, dropna=False, sort=False):
        comparisons = [('frozen_main', 'frozen_main_delta_DWUI')]
        if part.view.iloc[0] != 'native_terminal_sites':
            comparisons.append(('repro_same_PAS_same_contributors', 'repro_delta_primary_axis'))
        for comparison, refcol in comparisons:
            a = pd.to_numeric(part.delta_primary_axis, errors='coerce')
            b = pd.to_numeric(part[refcol], errors='coerce')
            valid = a.notna() & b.notna(); x = a[valid].to_numpy(); y = b[valid].to_numpy()
            result = dict(zip(keys, identity), comparison=comparison, union_gene_background=len(part),
                frozen_main_gene_background=int(part.frozen_main_delta_DWUI.notna().sum()),
                tool_measurable_genes=int(a.notna().sum()), reference_measurable_genes=int(b.notna().sum()),
                jointly_measurable_genes=int(valid.sum()), missing_tool_genes=int(a.isna().sum()),
                tool_count_qualified_genes=int(part.tool_qualified_pairs.gt(0).sum()),
                repro_same_PAS_count_qualified_genes=int(part.repro_same_PAS_qualified_pairs.gt(0).sum()),
                missing_reference_genes=int(b.isna().sum()), direction_agree=int((np.sign(x) == np.sign(y)).sum()),
                mean_abs_effect_difference=float(np.abs(x-y).mean()) if len(x) else np.nan,
                pearson=float(np.corrcoef(x,y)[0,1]) if len(x)>2 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan,
                spearman=float(pd.Series(x).rank().corr(pd.Series(y).rank())) if len(x)>2 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan)
            for threshold in [.02, .05]:
                ca = effect_class(x, threshold); cb = effect_class(y, threshold)
                result.update({f'{label}_{threshold}': int(mask.sum()) for label, mask in [
                    ('class_agree',ca == cb), ('both_change',(ca != 0) & (cb != 0)),
                    ('same_direction_change',(ca == cb) & (ca != 0)),
                    ('only_tool_change',(ca != 0) & (cb == 0)), ('only_repro_change',(ca == 0) & (cb != 0))]})
            rows.append(result)
    return pd.DataFrame(rows)
