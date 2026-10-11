"""Run-owned route, trajectory and homologous-axis extensions."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from downstream_runtime import Context, Study, bh
from pipeline import digest

from types import SimpleNamespace
B=1500
SEED=20261004
NODE='axis_extensions'
ROUTES={'CIP':'C','WOT':'W','GSE118258':'H1','GSE150311':'H2','178325_0618':'H3','178325_0809':'H3b','GSE221739':'H4'}

def run(ctx):
    global CTX,HERE,BASE,CHECKS,read,write,e
    CTX=ctx;HERE=ctx.out;BASE=ctx.analysis;CHECKS={}
    read=lambda name,**kw:pd.read_csv(ctx.register(BASE/name),sep='\t',**kw)
    write=lambda frame,name:ctx.table(frame,name)
    other_axes()
    homolog_summary(*build_homolog(*homolog_inputs()))
    e=SimpleNamespace(HERE=HERE,weighted_rank_corr=weighted_rank_corr,write=write,digest=digest)
    same_cohort()
    ctx.finish(dict(new_candidate_backgrounds=True,synchronized_gene_bootstrap=True,paired_layer_decomposition=True),
               numerical_checks=CHECKS,bootstrap=B,seed=SEED)

def ci_mean(x, rng, draws=B):
    x = np.asarray(x, float)
    vals = np.array([rng.choice(x, len(x), replace=True).mean() for _ in range(draws)])
    return float(x.mean()), *np.quantile(vals, [.025, .975])

def weighted_rank_corr(x, y, w):
    """Midranks in the reweighted empirical distribution (ties preserved)."""
    def ranks(a):
        _, ix = np.unique(a, return_inverse=True)
        mass = np.bincount(ix, weights=w)
        return (np.cumsum(mass) - mass / 2)[ix]
    a, b = ranks(x), ranks(y)
    a = a - np.average(a, weights=w)
    b = b - np.average(b, weights=w)
    den = np.sqrt(np.sum(w * a * a) * np.sum(w * b * b))
    return np.sum(w * a * b) / den if den > 0 else np.nan

def direction(x, y, w):
    ok = (abs(x) >= .02) & (abs(y) >= .02)
    return np.average((x[ok] * y[ok]) > 0, weights=w[ok]) if w[ok].sum() > 0 else np.nan

def homolog_inputs():
    one = read('cis/reciprocal_one_to_one.tsv.gz')
    assert one.human_id.is_unique and one.mouse_id.is_unique
    cat = read('cis/annotation_catalog.tsv.gz').set_index('id')
    usage = read('cis/annotated_state_usage.tsv.gz')
    total = usage.groupby(['code', 'gene_id', 'node']).usage.sum()
    assert np.allclose(total, 1, atol=1e-12, rtol=0)
    node = read('cross_overlap/all_mapped_nodes.tsv.gz')
    node = node[node.layer.eq('supported50') & node.measurement_pass]
    w = node.pivot(index=['code', 'gene_id'], columns='node', values='position_index')
    native = (w.P - w.S).dropna().rename('native').reset_index()
    u = usage[usage.node.isin(['P', 'S'])].pivot(
        index=['code', 'gene_id', 'id'], columns='node', values='usage').dropna().reset_index()
    u = u.rename(columns={'P': 'target', 'S': 'reference'}).merge(native, on=['code', 'gene_id'], validate='many_to_one')
    u['design'] = 'observed_PS'
    frames = [u]
    mapping = read('cis/detected_to_annotation.tsv.gz')
    maps = {c: q.set_index('native_pas_id').id.to_dict() for c, q in mapping.groupby('code')}
    cols = ['route', 'gene_id', 'contrast', 'design', 'PAS_ids', 'target_usage', 'reference_usage', 'weight', 'delta_DWUI']
    pairs = read('state_vectors/pair_effects.primary.tsv.gz', usecols=cols)
    records = []
    for (route, gene, contrast), q in pairs.groupby(['route', 'gene_id', 'contrast'], sort=False):
        code = ROUTES[route]
        lookup = maps['H3' if code == 'H3b' else code]
        assert q.PAS_ids.nunique() == 1
        ids = json.loads(q.PAS_ids.iloc[0])
        if not all(x in lookup for x in ids):
            # The cis catalog is limited to the fixed one-to-one ortholog universe.
            assert not any(x in lookup for x in ids)
            continue
        weight = q.weight.to_numpy()
        assert abs(weight.sum() - 1) < 1e-12
        a = weight @ np.stack(q.target_usage.map(json.loads))
        b = weight @ np.stack(q.reference_usage.map(json.loads))
        z = pd.DataFrame({'id': [lookup[x] for x in ids], 'target': a, 'reference': b}).groupby('id').sum()
        assert np.allclose(z.sum(), 1, atol=1e-12, rtol=0)
        records.extend(dict(code=code, gene_id=gene, id=idx, target=r.target, reference=r.reference,
                            native=float(weight @ q.delta_DWUI), design='qualified_' + contrast)
                       for idx, r in z.iterrows())
    frames.append(pd.DataFrame(records))
    all_usage = pd.concat(frames, ignore_index=True)
    all_usage['delta'] = all_usage.target - all_usage.reference
    write(all_usage, 'homolog_annotated_effects.tsv.gz')
    route = read('global/qualified_route_effects.tsv.gz')
    route = route[route.setting.eq('supported50_p08')].copy()
    route['code'] = route.route.map(ROUTES)
    check = all_usage[all_usage.design.ne('observed_PS')].drop_duplicates(['design', 'code', 'gene_id'])
    check = check.assign(contrast=check.design.str.replace('qualified_', '', regex=False))
    check = check.merge(route, on=['code', 'gene_id', 'contrast'], validate='one_to_one')
    assert np.allclose(check.native, check['delta_y'], atol=1e-12, rtol=0)
    CHECKS['qualified_native_effects_replayed'] = len(check)
    return one, cat, all_usage

def build_homolog(one, cat, usage):
    sites, genes, coverage = [], [], []
    for tolerance in [25, 50]:
        o = one[one['reciprocal' + str(tolerance)]].copy()
        for side in ['human', 'mouse']:
            o[side + '_signed_cut'] = o[side + '_id'].map(cat.cut_site) * o[side + '_id'].map(cat.strand).map({'+': 1, '-': -1})
        for design, all_u in usage.groupby('design', sort=False):
            codes = set(all_u.code)
            for h in sorted(codes - {'C', 'W'}):
                for m in ['C', 'W']:
                    key = dict(tolerance=tolerance, design=design, pair=h + '_' + m)
                    a = all_u[all_u.code.eq(h)].drop(columns=['code', 'design', 'gene_id'])
                    b = all_u[all_u.code.eq(m)].drop(columns=['code', 'design', 'gene_id'])
                    z = o.merge(a, left_on='human_id', right_on='id', validate='one_to_one').drop(columns='id')
                    z = z.merge(b, left_on='mouse_id', right_on='id', suffixes=('_h', '_m'), validate='one_to_one').drop(columns='id')
                    for name, val in key.items():
                        z[name] = val
                    z['n_matched'] = z.groupby('human_gene_id').human_id.transform('size')
                    z['mass_min'] = z.groupby('human_gene_id')[['target_h', 'reference_h', 'target_m', 'reference_m']].transform('sum').min(axis=1)
                    sites.append(z)
                    kept = 0
                    for gid, q in z.groupby('human_gene_id', sort=False):
                        q = q.sort_values('human_signed_cut')
                        if len(q) < 2 or not np.all(np.diff(q.mouse_signed_cut) > 0):
                            continue
                        rank = np.linspace(0, 1, len(q))
                        rec = dict(key, human_gene_id=gid, mouse_gene_id=q.mouse_gene_id.iloc[0],
                                   n_matched=len(q), mass_min=q.mass_min.iloc[0])
                        for side in ['h', 'm']:
                            a, b = q['target_' + side].to_numpy(), q['reference_' + side].to_numpy()
                            coord = q['human_signed_cut' if side == 'h' else 'mouse_signed_cut'].to_numpy()
                            dist = (coord - coord.min()) / (coord.max() - coord.min())
                            norm = a / a.sum() - b / b.sum() if min(a.sum(), b.sum()) > 0 else np.full(len(q), np.nan)
                            rec['aligned_' + side] = norm @ rank
                            rec['distance_' + side] = norm @ dist
                            rec['native_' + side] = q['native_' + side].iloc[0]
                        genes.append(rec)
                        kept += 1
                    coverage.append(dict(key, matched_sites=len(z), genes_any=z.human_gene_id.nunique(),
                                         genes_two_ordered=kept, genes_mass80=z.loc[z.mass_min.ge(.8) & z.n_matched.ge(2), 'human_gene_id'].nunique()))
    write(pd.concat(sites, ignore_index=True), 'homolog_site_pairs.tsv.gz')
    write(genes, 'homolog_gene_axes.tsv.gz')
    write(coverage, 'homolog_coverage.tsv')
    CHECKS['homolog_unique_mapping_and_unit_usage'] = True
    return pd.concat(sites, ignore_index=True), pd.DataFrame(genes)

def homolog_summary(sites, genes):
    rng = np.random.default_rng(SEED)
    rows = []
    for key, g in genes.groupby(['tolerance', 'design']):
        s = sites[(sites.tolerance == key[0]) & (sites.design == key[1])]
        for subset in ['two_sites', 'mass80']:
            z = g[g.mass_min.ge(.8)] if subset == 'mass80' else g
            z = z.dropna(subset=['aligned_h', 'aligned_m'])
            ids = sorted(z.human_gene_id.unique())
            look = {v: i for i, v in enumerate(ids)}
            bundles = []
            for pair, a in z.groupby('pair'):
                a = a.sort_values('human_gene_id')
                q = s[s.pair.eq(pair) & s.human_gene_id.isin(a.human_gene_id)].copy()
                if len(a) < 10:
                    rows.append(dict(tolerance=key[0], design=key[1], subset=subset, pair=pair, n_genes=len(a), status='fewer_than_10_genes'))
                    continue
                ix = np.array([look[x] for x in a.human_gene_id])
                si = np.array([look[x] for x in q.human_gene_id])
                sw = 1 / q.n_matched.to_numpy()
                xys = [(a['native_h'].to_numpy(), a['native_m'].to_numpy()),
                       (a.aligned_h.to_numpy(), a.aligned_m.to_numpy()),
                       (a.distance_h.to_numpy(), a.distance_m.to_numpy())]
                sx, sy = q.delta_h.to_numpy(), q.delta_m.to_numpy()
                bundles.append((pair, a, ix, si, sw, xys, sx, sy))
            if not bundles:
                continue
            values = np.empty((B + 1, len(bundles), 6))
            for draw in range(B + 1):
                counts = np.ones(len(ids)) if draw == 0 else rng.multinomial(len(ids), np.ones(len(ids)) / len(ids))
                for j, (_, a, ix, si, sw, xys, sx, sy) in enumerate(bundles):
                    wg, ws = counts[ix], counts[si] * sw
                    native, aligned, distance = [weighted_rank_corr(x, y, wg) for x, y in xys]
                    values[draw, j] = [native, aligned, distance, aligned - native,
                                      weighted_rank_corr(sx, sy, ws), direction(xys[1][0], xys[1][1], wg)]
            labels = ['native_rho', 'aligned_rho', 'matched_distance_rho', 'aligned_minus_native_rho',
                      'site_gene_equal_rho', 'aligned_direction_both_abs_ge02']
            for j, (pair, a, *_) in enumerate(bundles):
                for k, metric in enumerate(labels):
                    lo, hi = np.nanquantile(values[1:, j, k], [.025, .975])
                    rows.append(dict(tolerance=key[0], design=key[1], subset=subset, pair=pair, metric=metric,
                                     n_genes=len(a), estimate=values[0, j, k], low=lo, high=hi, status='estimated', bootstrap=B))
            avg = np.nanmean(values, axis=1)
            for k, metric in enumerate(labels):
                lo, hi = np.nanquantile(avg[1:, k], [.025, .975])
                rows.append(dict(tolerance=key[0], design=key[1], subset=subset, pair='pair_mean', metric=metric,
                                 n_genes=len(ids), estimate=avg[0, k], low=lo, high=hi, status='estimated', bootstrap=B))
            print('homolog', key, subset, len(ids), 'genes', len(bundles), 'pairs', flush=True)
    write(rows, 'homolog_summary.tsv')
    CHECKS['bootstrap_unit'] = 'ortholog gene, synchronized across pairs; pair effects conditional on sampled studies'

def other_axes():
    rng = np.random.default_rng(SEED + 1)
    a = read('profiles/architecture_gene_metrics.tsv.gz')
    assert np.allclose(a[['contribution_SR', 'contribution_RE', 'contribution_EP']].sum(1), 1)
    out = []
    for key, q in a.groupby(['species', 'route', 'cluster']):
        for segment in ['SR', 'RE', 'EP']:
            mean, lo, hi = ci_mean(q['contribution_' + segment], rng)
            out.append(dict(zip(['species', 'route', 'cluster'], key), segment=segment, n_genes=len(q),
                            contribution=mean, low=lo, high=hi, dominant_fraction=q.dominant_stage.eq(segment).mean()))
    write(out, 'trajectory_stage_summary.tsv')
    x = a[a.species.eq('mouse')].pivot(index=['gene_id', 'cluster'], columns='route', values=['transient_fraction_T', 'dominant_stage'])
    out = []
    for cluster, q in x.groupby(level='cluster'):
        same = q['dominant_stage']['CIP'].eq(q['dominant_stage']['WOT'])
        delta = q['transient_fraction_T']['CIP'].astype(float) - q['transient_fraction_T']['WOT'].astype(float)
        for label, vals in [('dominant_stage_agreement', same), ('mean_T_CIP_minus_WOT', delta)]:
            mean, lo, hi = ci_mean(vals, rng)
            out.append(dict(cluster=cluster, measure=label, n_genes=len(q), estimate=mean, low=lo, high=hi))
    write(out, 'trajectory_mouse_paired.tsv')
    write(a.groupby(['species', 'route', 'cluster']).agg(n_genes=('gene_id', 'size'), median_T=('transient_fraction_T', 'median')).reset_index(), 'trajectory_T_medians.tsv')
    f = read('global/qualified_route_effects.tsv.gz')
    out = []
    for key, q in f.groupby(['species', 'route', 'contrast']):
        w = q.pivot(index='gene_id', columns='setting', values='delta')
        common = w[['supported50_p08', 'expanded_primary_p08']].dropna()
        expansion = w.loc[w.supported50_p08.isna(), 'expanded_primary_p08'].dropna()
        supported = common.supported50_p08
        expanded = common.expanded_primary_p08
        assert w.loc[w.expanded_primary_p08.isna(), 'supported50_p08'].dropna().empty
        mean, lo, hi = ci_mean(expanded - supported, rng)
        total = w.expanded_primary_p08.mean() - supported.mean()
        common_part = len(common) / w.expanded_primary_p08.notna().sum() * mean
        entry_part = len(expansion) / w.expanded_primary_p08.notna().sum() * (expansion.mean() - supported.mean()) if len(expansion) else 0.
        assert np.isclose(total, common_part + entry_part)
        sa = np.where(supported >= .02, 1, np.where(supported <= -.02, -1, 0))
        sb = np.where(expanded >= .02, 1, np.where(expanded <= -.02, -1, 0))
        out.append(dict(zip(['species', 'route', 'contrast'], key), common_genes=len(common), expanded_only=len(expansion),
                        supported_mean=supported.mean(), expanded_common_mean=expanded.mean(), expanded_all_mean=w.expanded_primary_p08.mean(),
                        paired_change=mean, low=lo, high=hi, total_change=total, common_reestimate_contribution=common_part,
                        added_background_contribution=entry_part, labels_changed_02=int((sa != sb).sum()),
                        opposite_strong_02=int((sa * sb < 0).sum())))
    write(out, 'layer_paired_decomposition.tsv')
    write(read('global/direction_stability.tsv'), 'route_direction_audit.tsv')
    read('global/leave_target_date_summary.tsv').to_csv(HERE / 'leave_date_summary.tsv', sep='\t', index=False)
    out = []
    for unit in ['well', 'pool']:
        d = read(unit + '_inference/gene_metric_tests.tsv.gz')
        for metric, q in d.groupby('metric'):
            out.append(dict(unit=unit, metric=metric, genes=len(q), significant=int(q.joint_q05.sum()),
                            min_units=q.n_units.min(), max_units=q.n_units.max(), min_p=q.p.min(), min_joint_q=q.q_BH_joint.min(),
                            median_abs_effect=q.effect.abs().median(), n_abs_ge02=int(q.effect.abs().ge(.02).sum())))
        floors = d[d.metric.eq('DWUI')].groupby('n_units').size().reset_index(name='genes')
        floors['exact_two_sided_sign_min_p'] = 2. / 2. ** floors.n_units
        write(floors, unit + '_resolution.tsv')
    write(out, 'replicate_inference_audit.tsv')
    f = read('apa_gex/four_class_nulls.tsv.gz')
    f['p_two_sided'] = np.minimum(1., 2 * np.minimum(f.p_enrich, f.p_deplete))
    q = f.groupby(['species', 'threshold', 'null_model', 'baseline_bins', 'cell']).agg(
        observed=('observed', 'first'), null_min=('null_mean', 'min'), null_max=('null_mean', 'max'),
        enrichment_min=('enrichment', 'min'), enrichment_max=('enrichment', 'max'), p_both=('p_two_sided', 'max')).reset_index()
    q['q_all_32_per_binning'] = q.groupby('baseline_bins').p_both.transform(bh)
    q['APA_relation'] = np.where(q.cell.isin(['A', 'B']), 'same', 'opposite')
    q['GEX_relation'] = np.where(q.cell.isin(['A', 'C']), 'same', 'opposite')
    write(q, 'apa_gex_four_class_exploratory.tsv')
    write(read('apa_gex/prespecified_eight_test_family.tsv'), 'apa_gex_primary_tests.tsv')
    d = read('intermediate/WOT_date_stratified_tests.tsv')
    assert np.allclose(bh(d.p_value), d.q_all_requested_contrasts)
    write(d.groupby('contrast').agg(requested=('gene_id', 'size'), estimable=('status', lambda s: sum(s == 'tested')),
                                    joint_q05=('q_all_requested_contrasts', lambda s: sum(s < .05)),
                                    min_joint_q=('q_all_requested_contrasts', 'min')).reset_index(), 'date_conditioned_summary.tsv')
    write(d[d.q_all_requested_contrasts.lt(.05)], 'date_conditioned_hits.tsv')
    for src, dst in [('external_pas/native_comparison_summary.tsv', 'external_pas_summary.tsv'),
                     ('external_pas/vector_summary.tsv', 'external_pas_assignment_tests.tsv'),
                     ('longread/structural_summary.tsv', 'longread_structure.tsv'),
                     ('longread/family_eligibility.tsv', 'longread_eligibility_original_status.tsv'),
                     ('ago2/biological_unit_support.tsv', 'ago2_support.tsv')]:
        write(read(src), dst)
    CHECKS['layer_change_exact_decomposition'] = True
    CHECKS['date_FDR_replayed'] = len(d)
    print('remaining seven directions summarized', flush=True)

def same_cohort():
    path = e.HERE / 'homolog_gene_axes.tsv.gz'
    f = pd.read_csv(CTX.register(path),sep='\t')
    rng = np.random.default_rng(20261006)
    rows = []
    for tol in [25, 50]:
        q = f[f.tolerance.eq(tol) & f.mass_min.ge(.8)]
        ps = q[q.design.eq('qualified_P_S')]
        pm = q[q.design.eq('qualified_P_M')]
        joined = ps.merge(pm, on=['pair', 'human_gene_id'], suffixes=('_ps', '_pm'), validate='one_to_one')
        ids = sorted(joined.human_gene_id.unique())
        lookup = {v: i for i, v in enumerate(ids)}
        groups = list(joined.groupby('pair'))
        assert joined.n_matched_ps.eq(joined.n_matched_pm).all()
        values = np.empty((1501, len(groups), 3, 3))
        for draw in range(1501):
            w = np.ones(len(ids)) if draw == 0 else rng.multinomial(len(ids), np.ones(len(ids)) / len(ids))
            for j, (_, g) in enumerate(groups):
                weight = w[[lookup[x] for x in g.human_gene_id]]
                for k, axis in enumerate(['native', 'aligned', 'distance']):
                    a = e.weighted_rank_corr(g[axis+'_h_ps'], g[axis+'_m_ps'], weight)
                    b = e.weighted_rank_corr(g[axis+'_h_pm'], g[axis+'_m_pm'], weight)
                    values[draw, j, k] = [a, b, b - a]
        for j, (pair, g) in enumerate(groups + [('pair_mean', joined)]):
            v = values[:, j] if j < len(groups) else values.mean(axis=1)
            for a, axis in enumerate(['native', 'aligned', 'distance']):
                for k, metric in enumerate(['P_S_rho', 'P_M_rho', 'P_M_minus_P_S_rho']):
                    lo, hi = np.quantile(v[1:, a, k], [.025, .975])
                    rows.append(dict(tolerance=tol, pair=pair, n_genes=g.human_gene_id.nunique(), axis=axis, metric=metric,
                                     estimate=v[0, a, k], low=lo, high=hi))
    e.write(rows, 'homolog_same_cohort_reference_comparison.tsv')
    (e.HERE / 'provenance_followup.json').write_text(json.dumps(dict(input_sha256=e.digest(path), script_sha256=e.digest(__file__), seed=20261006), indent=2))

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,NODE,a.out))
