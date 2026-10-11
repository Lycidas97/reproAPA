"""Conditional measurement and matching sensitivity on the frozen APA inventory."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import time

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / 'notes/2026-09-30-multimetric-statistics'
M = Path('/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot')
H = Path('/Users/lycidas/work/科研/projects/2026-09-09_human-time-aware-integration')
SEED = 20261002
BOOT = 500
SPLIT = 100
THIN = 20
BASE = ['species', 'study_id', 'route', 'gene_id', 'gene_symbol', 'contrast', 'design']
SPECS = {
    'GSE103221': ('mouse', 'gse103221_pluripotency_reprogramming_core'),
    'GSE115943': ('mouse', 'gse115943_pluripotency_reprogramming_core'),
    'GSE118258': ('human', 'gse118258_human_pluripotency_reprogramming'),
    'GSE150311': ('human', 'gse150311_human_pluripotency_reprogramming'),
    'GSE178325': ('human', 'gse178325_human_chemical_pluripotency_reprogramming'),
    'GSE221739': ('human', 'gse221739_human_pluripotency_reprogramming'),
}
INPUTS = {}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


legacy = module('followup_weighting', OLD / 'extract.py')
nullmod = module('followup_nulls', OLD / 'nulls/analyze.py')


def read(path, **kw):
    path = Path(path)
    with path.open('rb') as f:
        INPUTS[str(path)] = hashlib.file_digest(f, 'sha256').hexdigest()
    return pd.read_csv(path, sep='\t', **kw)


def save(frame, path):
    path = Path(path)
    path.parent.mkdir(exist_ok=True, parents=True)
    frame.to_csv(path, sep='\t', index=False, na_rep='NA', float_format='%.12g',
                 compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)


def rng_for(*keys):
    key = '|'.join(map(str, keys))
    return np.random.default_rng(SEED + int(hashlib.sha256(key.encode()).hexdigest()[:12], 16))


def load_study(study, axis='native'):
    species, folder = SPECS[study]
    data = (M / 'data/metacell-pas-v7-v1/studies' if species == 'mouse'
            else H / 'data/results_metacell_pas_v2/studies') / folder
    meta = read(M / 'data/apa-majority-shared-v1/profiles' / folder / 'selected_metacells.tsv.gz')
    full = read(data / 'metacells.tsv.gz')
    assert np.array_equal(full.metacell_id.iloc[meta.matrix_row], meta.metacell_id)
    axes = read(OLD / 'frozen_PAS_axes.tsv.gz')
    axes = axes[axes.study_id.eq(study)].reset_index(drop=True)
    feature = read(data / 'features.tsv.gz', usecols=['atlas_pas_id', 'feature_index_1based'])
    axes = axes.merge(feature, on='atlas_pas_id', validate='one_to_one')
    countpath = data / 'metacell_by_pas_counts.npz'
    with countpath.open('rb') as f:
        INPUTS[str(countpath)] = hashlib.file_digest(f, 'sha256').hexdigest()
    counts = sparse.load_npz(countpath).tocsr()[meta.matrix_row.to_numpy()]
    counts = counts[:, axes.feature_index_1based.to_numpy() - 1].tocsr()
    assert np.issubdtype(counts.data.dtype, np.integer) and (counts.data >= 0).all()
    if axis == 'common':
        mapping = read(HERE / 'common_axes.tsv.gz')
        mapping = mapping[mapping.study_id.eq(study)]
        selected = axes.reset_index().merge(mapping, on=['study_id', 'atlas_pas_id', 'primary_gene_id'], validate='one_to_one')
        regions = selected[['primary_gene_id', 'cluster_id', 'common_position', 'common_cut_site']].drop_duplicates().sort_values('cluster_id')
        ci = pd.Index(regions.cluster_id).get_indexer(selected.cluster_id)
        transform = sparse.csr_matrix((np.ones(len(selected), np.int64),
            (selected['index'].to_numpy(), ci)), shape=(len(axes), len(regions)))
        counts = (counts @ transform).tocsr()
        axes = regions.rename(columns={'cluster_id': 'atlas_pas_id', 'common_position': 'position', 'common_cut_site': 'cut_site'}).reset_index(drop=True)
        signed_order = axes.groupby('primary_gene_id').position.rank(method='first')
        axes['rank'] = signed_order
        axes['n'] = axes.groupby('primary_gene_id').position.transform('size')
    return meta.reset_index(drop=True), axes, counts


def build_axes():
    axes = read(OLD / 'frozen_PAS_axes.tsv.gz')
    locations = []
    for study in ['GSE103221', 'GSE115943', 'GSE118258', 'GSE178325']:
        species, folder = SPECS[study]
        base = (M / 'data/metacell-pas-v7-v1/studies' if species == 'mouse'
                else H / 'data/results_metacell_pas_v2/studies') / folder
        feature = read(base / 'features.tsv.gz', usecols=['atlas_pas_id', 'chrom', 'genome_build'])
        locations.append(feature.assign(study_id=study))
    axes = axes.merge(pd.concat(locations), on=['study_id', 'atlas_pas_id'], validate='one_to_one')
    rows = []
    for species, studies in [('mouse', ['GSE103221', 'GSE115943']), ('human', ['GSE118258', 'GSE178325'])]:
        for gene, group in axes[axes.study_id.isin(studies)].groupby('primary_gene_id'):
            if (group.study_id.nunique() != 2 or group.strand.nunique() != 1
                    or group.chrom.nunique() != 1 or group.genome_build.nunique() != 1):
                continue
            clusters = []
            for row in group.sort_values(['cut_site', 'study_id', 'atlas_pas_id']).to_dict('records'):
                if not clusters or row['cut_site'] - clusters[-1][0]['cut_site'] > 50:
                    clusters.append([])
                clusters[-1].append(row)
            clusters = [c for c in clusters if {r['study_id'] for r in c} == set(studies)]
            if len(clusters) < 2:
                continue
            coord = np.array([np.mean([r['cut_site'] for r in c]) for c in clusters])
            oriented = coord * (1 if group.strand.iloc[0] == '+' else -1)
            if np.ptp(oriented) == 0:
                continue
            positions = (oriented - oriented.min()) / np.ptp(oriented)
            for j, cluster in enumerate(clusters):
                for row in cluster:
                    rows.append(dict(species=species, study_id=row['study_id'], primary_gene_id=gene,
                        atlas_pas_id=row['atlas_pas_id'], cluster_id=f'{species}|{gene}|{j}',
                        common_cut_site=coord[j], common_position=positions[j], cluster_span=max(r['cut_site'] for r in cluster)-min(r['cut_site'] for r in cluster)))
    frame = pd.DataFrame(rows)
    assert frame.cluster_span.le(50).all()
    assert not frame.duplicated(['study_id', 'atlas_pas_id']).any()
    save(frame, HERE / 'common_axes.tsv.gz')
    print('COMMON AXES', frame.groupby('species').primary_gene_id.nunique().to_dict(), flush=True)


def cycle_covariates(study, meta):
    path = ROOT / 'notes/2026-09-26-state-comparability-audit/expression_inputs.npz'
    with np.load(path, allow_pickle=False) as a:
        if study + '_meta_json' not in a:
            return None
        original = pd.DataFrame(json.loads(str(a[study + '_meta_json'])))
        species = SPECS[study][0]
        symbols = a[species + '_marker_symbols'].tolist()
        programs = json.loads(str(a[species + '_program_json']))
        marker = a[study + '_marker_logcp10k']
        cov = original[['metacell_id', 'n_cells']].copy()
        cov['gex_per_cell'] = a[study + '_total_gex'] / original.n_cells
        cov['panel_complexity'] = (a[study + '_counts'] > 0).sum(1)
        for name in ['Cycle_S', 'Cycle_G2M']:
            cov[name] = marker[:, [symbols.index(g.upper()) for g in programs[name].split()]].mean(1)
    out = meta.merge(cov.drop(columns='n_cells'), on='metacell_id', validate='one_to_one')
    assert len(out) == len(meta)
    for route, part in out.groupby('route'):
        for name in ['Cycle_S', 'Cycle_G2M']:
            edges = np.unique(np.quantile(part[name], [1/3, 2/3]))
            out.loc[part.index, name + '_bin'] = np.searchsorted(edges, part[name], side='right')
    out['cycle_bin'] = (out.Cycle_S_bin.astype(int) * 3 + out.Cycle_G2M_bin.astype(int)).astype(int)
    return out


def hypergeom_thin(counts, depth, repeats, rng):
    """Exactly depth molecules, without replacement, for every eligible metacell."""
    counts = np.asarray(counts, np.int64)
    n = counts.sum(1)
    assert (n >= depth).all()
    out = np.zeros((repeats, len(counts), counts.shape[1]), np.int64)
    remaining = np.full((repeats, len(counts)), depth, np.int64)
    population = n.copy()
    for j in range(counts.shape[1] - 1):
        good = counts[:, j]
        take = rng.hypergeometric(np.broadcast_to(good, remaining.shape),
                                 np.broadcast_to(population - good, remaining.shape), remaining)
        out[:, :, j] = take
        remaining -= take
        population -= good
    out[:, :, -1] = remaining
    assert np.all(out.sum(2) == depth) and np.all(out <= counts[None])
    return out


def state_record(c, cells, rng, resample=True):
    keep = c.sum(1) >= 20
    c, cells = c[keep], cells[keep]
    if len(c) < 3 or cells.sum() < 30:
        return None
    p = c / c.sum(1)[:, None]
    mean = p.mean(0)
    if not resample:
        return dict(mean=mean, counts=c.sum(0), positive=(c > 0).sum(0),
                    n_mc=len(c), n_cells=int(cells.sum()))
    weights = rng.multinomial(len(p), np.full(len(p), 1 / len(p)), size=BOOT) / len(p)
    boot = weights @ p
    split = np.full((SPLIT, c.shape[1]), np.nan)
    if len(c) >= 6:
        order = rng.permuted(np.broadcast_to(np.arange(len(c)), (SPLIT, len(c))), axis=1)
        half = len(c) // 2
        a, b = order[:, :half], order[:, half:]
        good = (cells[a].sum(1) >= 30) & (cells[b].sum(1) >= 30)
        split[good] = p[a[good]].mean(1) - p[b[good]].mean(1)
    thinned = {}
    for depth in [20, 40]:
        ix = c.sum(1) >= depth
        if ix.sum() < 3 or cells[ix].sum() < 30:
            continue
        t = hypergeom_thin(c[ix], depth, THIN, rng)
        thinned[depth] = dict(mean=(t / depth).mean(1), counts=t.sum(1), positive=(t > 0).sum(1), n_mc=int(ix.sum()), n_cells=int(cells[ix].sum()))
    return dict(mean=mean, boot=boot, split=split, counts=c.sum(0), positive=(c > 0).sum(0),
                c=c, p=p, cells=cells, mask=keep, n_mc=len(c), n_cells=int(cells.sum()), thinned=thinned)


def joint_good(a, b):
    return int(((a['counts'] + b['counts'] >= 10) & (a['positive'] + b['positive'] >= 2)).sum()) >= 2


def cycle_pair(a, b, bins_a, bins_b):
    bins_a, bins_b = bins_a[a['mask']], bins_b[b['mask']]
    common = sorted(set(bins_a) & set(bins_b))
    if not common:
        return None
    ia, ib = np.isin(bins_a, common), np.isin(bins_b, common)
    if min(ia.sum(), ib.sum()) < 3 or min(a['cells'][ia].sum(), b['cells'][ib].sum()) < 30:
        return None
    sa = {'counts': a['c'][ia].sum(0), 'positive': (a['c'][ia] > 0).sum(0)}
    sb = {'counts': b['c'][ib].sum(0), 'positive': (b['c'][ib] > 0).sum(0)}
    if not joint_good(sa, sb):
        return None
    weights = np.array([min((bins_a == z).sum(), (bins_b == z).sum()) for z in common], float)
    weights /= weights.sum()
    ma = weights @ np.array([a['p'][bins_a == z].mean(0) for z in common])
    mb = weights @ np.array([b['p'][bins_b == z].mean(0) for z in common])
    return ma, mb, int(ia.sum()), int(ib.sum()), len(common)


def make_route(key, ua, ub, axis, original, pairs, variant, replicate=0):
    rec = dict(zip(BASE, key))
    d = ua - ub
    w = axis.position.to_numpy()
    rw = axis['rank'].to_numpy() / len(axis)
    rec.update(variant=variant, replicate=replicate, delta_WPP=float(d @ w), reference_WPP=float(ub @ w),
               delta_RWUI=float(d @ rw), reference_RWUI=float(ub @ rw), rwui_min_weight=float(rw.min()),
               rwui_max_weight=float(rw.max()), MPRO=float(np.ptp(d)), reference_concentration=float(ub @ ub),
               reference_min_usage=float(ub.min()), reference_max_usage=float(ub.max()), PAS_number=len(axis),
               mean_GEX=float(original.mean_GEX), APA_coverage=float(original.APA_coverage),
               n_pairs=len(pairs), series_pattern=rec['route'])
    assert np.isclose(ua.sum(), 1) and np.isclose(ub.sum(), 1)
    return rec


def measure(study, axis_kind):
    started = time.monotonic()
    out = HERE / axis_kind / study
    out.mkdir(exist_ok=True, parents=True)
    meta, axes, counts = load_study(study, axis_kind)
    cov = cycle_covariates(study, meta) if axis_kind == 'native' else None
    pairs = read(OLD / 'pair_effects.tsv.gz')
    pairs = pairs[pairs.study_id.eq(study)].copy()
    old_routes = read(OLD / 'route_effects.tsv.gz')
    old_routes = old_routes[old_routes.study_id.eq(study)].set_index(BASE)
    old_states = read(OLD / 'library_state_metrics.tsv.gz')
    old_states = old_states[old_states.study_id.eq(study)].set_index(['route', 'library_key', 'state', 'gene_id'])
    axis_groups = {g: np.asarray(ix) for g, ix in axes.groupby('primary_gene_id', sort=False).groups.items()}
    state_groups = {key: np.flatnonzero(meta.route.eq(key[0]) & meta.library_key.eq(key[1]) & meta.node.eq(key[2]) & meta.purity.ge(.8))
                    for key in {(r.route, r.library_key_a, 'P') for r in pairs.itertuples()} |
                    {(r.route, r.library_key_b, r.contrast[-1]) for r in pairs.itertuples()}}
    summaries, variants, pseudo, contrasts, coverage, state_noise = [], [], [], [], [], []
    max_error = 0.0
    for gn, (gene, gp) in enumerate(pairs.groupby('gene_id', sort=True)):
        if gene not in axis_groups:
            coverage.append(dict(study_id=study, gene_id=gene, status='no_common_axis'))
            continue
        cols = axis_groups[gene]
        ax = axes.iloc[cols]
        weights = ax.position.to_numpy()
        cg = counts[:, cols].toarray()
        state = {}
        needed = {(r.route, r.library_key_a, 'P') for r in gp.itertuples()} | {(r.route, r.library_key_b, r.contrast[-1]) for r in gp.itertuples()}
        for k in sorted(needed):
            ix = state_groups[k]
            value = state_record(cg[ix], meta.n_cells.to_numpy()[ix], rng_for(study, axis_kind, gene, *k), axis_kind == 'native')
            if value is None:
                continue
            state[k] = value
            if axis_kind == 'native':
                orig = old_states.loc[(*k, gene)]
                err = abs(value['mean'] @ weights - orig.WPP)
                max_error = max(max_error, err)
                assert err < 1e-8 and value['n_mc'] == orig.n_mc and value['n_cells'] == orig.n_cells
                x = value['split'] @ weights
                ok = np.isfinite(x)
                if ok.any():
                    state_noise.append(dict(study_id=study, gene_id=gene, route=k[0], library_key=k[1], state=k[2],
                        n_mc=value['n_mc'], split_estimable=int(ok.sum()), abs_dwui_p90=float(np.quantile(abs(x[ok]), .9)),
                        abs_dwui_p95=float(np.quantile(abs(x[ok]), .95)), dwui002_fraction=float((abs(x[ok]) >= .02).mean()),
                        mpro010_fraction=float((np.ptp(value['split'][ok], axis=1) >= .1).mean())))
        for key, part in gp.groupby(BASE, sort=False):
            valid = np.array([(r.route, r.library_key_a, 'P') in state and
                              (r.route, r.library_key_b, r.contrast[-1]) in state and
                              joint_good(state[(r.route, r.library_key_a, 'P')], state[(r.route, r.library_key_b, r.contrast[-1])])
                              for r in part.itertuples()])
            baseline = part[valid].copy()
            if baseline.empty:
                coverage.append(dict(study_id=study, gene_id=gene, route=key[2], contrast=key[-2], status='no_qualified_pairs'))
                continue
            baseline = legacy.pair_weights(baseline.drop(columns='weight'))
            pr = list(baseline.itertuples())
            aa = [state[(r.route, r.library_key_a, 'P')] for r in pr]
            bb = [state[(r.route, r.library_key_b, r.contrast[-1])] for r in pr]
            pw = baseline.weight.to_numpy()
            ua = sum(w * a['mean'] for w, a in zip(pw, aa))
            ub = sum(w * b['mean'] for w, b in zip(pw, bb))
            original = old_routes.loc[key]
            rec = make_route(key, ua, ub, ax, original, baseline, axis_kind)
            variants.append(rec)
            if axis_kind == 'native':
                err = abs(rec['delta_WPP'] - original.delta_WPP)
                max_error = max(max_error, err)
                assert err < 1e-8 and len(baseline) == len(part)
                boot_a = sum(w * a['boot'] for w, a in zip(pw, aa))
                boot_b = sum(w * b['boot'] for w, b in zip(pw, bb))
                d = (boot_a - boot_b) @ weights
                pm = np.ptp(boot_a - boot_b, axis=1)
                noise_pm95 = float(np.quantile(np.ptp((boot_a - boot_b) - (ua - ub), axis=1), .95))
                summaries.append(rec | dict(dwui_lo=np.quantile(d, .025), dwui_hi=np.quantile(d, .975),
                    same_sign_probability=float((d * rec['delta_WPP'] > 0).mean()),
                    original_direction_002_probability=float(((d * np.sign(rec['delta_WPP'])) >= .02).mean()),
                    mpro_lo=np.quantile(pm, .025), mpro_hi=np.quantile(pm, .975), mpro010_probability=float((pm >= .1).mean()),
                    mpro_noise_p95=noise_pm95, mpro_exceeds_noise95=bool(rec['MPRO'] > noise_pm95),
                    interval_excludes_zero=bool(np.quantile(d, .025) > 0 or np.quantile(d, .975) < 0)))
                # Full-depth calibrated pseudo error: half of two independent within-state half differences.
                pu = .5 * sum(w * (a['split'] + b['split']) for w, a, b in zip(pw, aa, bb))
                for j in range(SPLIT):
                    pseudo.append(dict(zip(BASE, key)) | dict(replicate=j, delta_WPP=float(pu[j] @ weights), MPRO=float(np.ptp(pu[j])),
                        null_kind='half_scaled_sum_of_within_state_differences'))
                for depth in [20, 40]:
                    weight_cache = {tuple(np.ones(len(baseline), bool)): (baseline, pw)}
                    for j in range(THIN):
                        keep = np.array([depth in a['thinned'] and depth in b['thinned'] and
                            int(((a['thinned'][depth]['counts'][j] + b['thinned'][depth]['counts'][j] >= 10) &
                                (a['thinned'][depth]['positive'][j] + b['thinned'][depth]['positive'][j] >= 2)).sum()) >= 2
                            for a, b in zip(aa, bb)])
                        if not keep.any():
                            continue
                        mask_key = tuple(keep)
                        if mask_key not in weight_cache:
                            sub = legacy.pair_weights(baseline.loc[keep].drop(columns='weight'))
                            weight_cache[mask_key] = (sub, sub.weight.to_numpy())
                        sub, ww = weight_cache[mask_key]
                        at = [a['thinned'][depth]['mean'][j] for a, keepa in zip(aa, keep) if keepa]
                        bt = [b['thinned'][depth]['mean'][j] for b, keepb in zip(bb, keep) if keepb]
                        variants.append(make_route(key, ww @ np.stack(at), ww @ np.stack(bt), ax, original, sub, f'depth{depth}', j))
                if cov is not None:
                    matched, flags = [], []
                    for r, a, b in zip(pr, aa, bb):
                        ka, kb = (r.route, r.library_key_a, 'P'), (r.route, r.library_key_b, r.contrast[-1])
                        x = cycle_pair(a, b, cov.cycle_bin.to_numpy()[state_groups[ka]], cov.cycle_bin.to_numpy()[state_groups[kb]])
                        flags.append(x is not None)
                        if x is not None:
                            matched.append(x)
                    if any(flags):
                        sub = legacy.pair_weights(baseline.loc[flags].drop(columns='weight'))
                        ww = sub.weight.to_numpy()
                        ma, mb = ww @ np.stack([x[0] for x in matched]), ww @ np.stack([x[1] for x in matched])
                        variants.append(make_route(key, ma, mb, ax, original, sub, 'cycle_matched'))
                # Associations are conditional descriptive diagnostics, not independent gene tests.
                source_depth = sum(w * b['counts'].sum() / b['n_mc'] for w, b in zip(pw, bb))
                target_depth = sum(w * a['counts'].sum() / a['n_mc'] for w, a in zip(pw, aa))
                contrasts.append(dict(zip(BASE, key)) | dict(delta_WPP=rec['delta_WPP'], delta_log_pas_depth=np.log2(target_depth/source_depth),
                    source_pas_depth=source_depth, target_pas_depth=target_depth))
            coverage.append(dict(study_id=study, gene_id=gene, route=key[2], contrast=key[-2], status='available',
                                 original_pairs=len(part), retained_pairs=len(baseline)))
        if gn % 500 == 0:
            print('MEASURE', study, axis_kind, gn, len(variants), round(time.monotonic()-started, 1), flush=True)
    save(pd.DataFrame(variants), out / 'route_variants.tsv.gz')
    save(pd.DataFrame(summaries), out / 'conditional_intervals.tsv.gz')
    save(pd.DataFrame(pseudo), out / 'pseudo_route_effects.tsv.gz')
    save(pd.DataFrame(contrasts), out / 'depth_associations.tsv.gz')
    save(pd.DataFrame(coverage), out / 'coverage.tsv.gz')
    save(pd.DataFrame(state_noise), out / 'same_state_split_summary.tsv.gz')
    audit = dict(status='PASS', study=study, axis=axis_kind, original_genes=pairs.gene_id.nunique(),
                 interval_records=len(summaries), route_records=len(variants), max_baseline_error=max_error,
                 seconds=time.monotonic()-started, input_sha256=INPUTS, bootstrap=BOOT, split_half=SPLIT,
                 thinning=THIN, inference_scope='within-library metacell measurement; not biological-replicate uncertainty')
    (out / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print('COMPLETE', study, axis_kind, len(variants), max_error, flush=True)


def collapse_pattern(routes):
    nullmod.THRESHOLDS = {'WPP': [.02, .05], 'MPRO': [.1, .2]}
    events = nullmod.collapse(routes)
    patterns = routes.groupby(['species', 'contrast', 'study_id', 'gene_id']).route.agg(lambda x: ';'.join(sorted(x)))
    keys = pd.MultiIndex.from_frame(events[['species', 'contrast', 'unit', 'gene_id']])
    events['series_pattern'] = patterns.reindex(keys).to_numpy()
    return events


def run_nulls():
    routes = pd.concat([pd.read_csv(p, sep='\t') for p in sorted((HERE / 'native').glob('*/route_variants.tsv.gz'))], ignore_index=True)
    commons = [pd.read_csv(p, sep='\t') for p in sorted((HERE / 'common').glob('*/route_variants.tsv.gz'))]
    if commons:
        routes = pd.concat([routes, *commons], ignore_index=True)
    intervals = pd.concat([pd.read_csv(p, sep='\t') for p in sorted((HERE / 'native').glob('*/conditional_intervals.tsv.gz'))])
    confident = routes[(routes.variant == 'native') & (routes.replicate == 0)].merge(
        intervals[BASE + ['interval_excludes_zero', 'mpro_exceeds_noise95']], on=BASE, validate='one_to_one')
    confident['variant'] = 'interval_supported'
    summaries, members, diagnostics, archives = [], [], [], {}
    for variant in ['native', 'interval_supported', 'depth20', 'depth40', 'cycle_matched', 'common']:
        f = confident if variant == 'interval_supported' else routes[(routes.variant == variant) & (routes.replicate == 0)].copy()
        if f.empty:
            continue
        events = collapse_pattern(f)
        if variant == 'interval_supported':
            keys = pd.MultiIndex.from_frame(events[['species', 'contrast', 'unit', 'gene_id']])
            for metric, column in [('WPP', 'interval_excludes_zero'), ('MPRO', 'mpro_exceeds_noise95')]:
                support = f.groupby(['species', 'contrast', 'study_id', 'gene_id'])[column].all()
                keep = support.reindex(keys).to_numpy()
                events.loc[events.metric.eq(metric) & ~keep, 'event'] = 0
        save(events, HERE / 'nulls' / f'{variant}_study_events.tsv.gz')
        for (metric, species, contrast, threshold), frame in events.groupby(['metric', 'species', 'contrast', 'threshold']):
            if frame.unit.nunique() < 2:
                diagnostics.append(dict(variant=variant, metric=metric, species=species, contrast=contrast, threshold=threshold, status='insufficient_studies'))
                continue
            genes = sorted(frame.gene_id.unique())
            gi = {g: i for i, g in enumerate(genes)}
            for pattern_matching in [False, True] if variant == 'native' else [True]:
                tag = f'{variant}|{metric}|{species}|{contrast}|{threshold}|series{int(pattern_matching)}'
                rng = rng_for(tag)
                op, om = np.zeros(len(genes), int), np.zeros(len(genes), int)
                dp, dm = np.zeros((10000, len(genes)), np.int8), np.zeros((10000, len(genes)), np.int8)
                diagparts = []
                for unit, sub in frame.groupby('unit'):
                    sf = nullmod.stratify(sub, metric, threshold)
                    if pattern_matching:
                        sf['stratum'] = sf.stratum + '|' + sf.series_pattern
                    draws, _, diag = nullmod.shuffle(sf, rng, 10000)
                    idx = np.array([gi[g] for g in sf.gene_id])
                    op[idx] += (sf.event == 1).to_numpy()
                    om[idx] += (sf.event == -1).to_numpy()
                    dp[:, idx] += draws == 1
                    dm[:, idx] += draws == -1
                    diagparts.append(diag.assign(unit=unit))
                observed_mask = np.maximum(op, om) >= 2
                null = (np.maximum(dp, dm) >= 2).sum(1)
                row = dict(variant=variant, metric=metric, species=species, contrast=contrast, threshold=threshold,
                           series_matching=pattern_matching, n_studies=frame.unit.nunique(), n_genes=len(genes)) | nullmod.summarize(observed_mask.sum(), null)
                row['excess_vs_null_mean'] = row['observed'] - row['null_mean']
                row['excess_vs_null_median'] = row['observed'] - row['null_median']
                # These are transformed null ranges, explicitly NOT effect CIs.
                row['excess_null_range_lo'] = row['observed'] - row['null_p975']
                row['excess_null_range_hi'] = row['observed'] - row['null_p025']
                diag = pd.concat(diagparts)
                row['expected_labels_changed'] = float(diag.expected_labels_changed.sum())
                row['singleton_fraction'] = float(diag.loc[diag.singleton, 'genes'].sum() / len(frame))
                row['frozen_label_fraction'] = float(diag.loc[~diag.mixed_labels, 'genes'].sum() / len(frame))
                summaries.append(row)
                archives[tag] = null
                for i in np.flatnonzero(observed_mask):
                    members.append(dict(variant=variant, metric=metric, species=species, contrast=contrast,
                        threshold=threshold, series_matching=pattern_matching, gene_id=genes[i], direction=1 if op[i] >= om[i] else -1))
                print('NULL', tag, row['observed'], round(row['null_mean'], 2), row['empirical_p'], flush=True)
    summary = pd.DataFrame(summaries)
    # Keep the eight planned comparisons even when a sensitivity loses a study.
    summary['q_within_variant_metric'] = summary.groupby(['variant', 'metric', 'series_matching']).empirical_p.transform(
        lambda x: nullmod.bh(np.r_[x.to_numpy(), np.ones(8-len(x))])[:len(x)])
    summary['q_all_followup_tests'] = nullmod.bh(np.r_[summary.empirical_p.to_numpy(), np.ones(len(diagnostics))])[:len(summary)]
    save(summary, HERE / 'nulls/recurrence_summary.tsv')
    save(pd.DataFrame(members), HERE / 'nulls/recurrence_members.tsv.gz')
    save(pd.DataFrame(diagnostics), HERE / 'nulls/unavailable_tests.tsv')
    np.savez_compressed(HERE / 'nulls/null_draws.npz', **archives)
    # Gene identities/directions and missing comparisons remain separate from threshold changes.
    base = pd.DataFrame(members)
    baseline = base[(base.variant == 'native') & base.series_matching]
    overlaps = []
    for key, f in base[base.series_matching].groupby(['variant', 'metric', 'species', 'contrast', 'threshold']):
        ref = baseline[(baseline.metric == key[1]) & (baseline.species == key[2]) & (baseline.contrast == key[3]) & (baseline.threshold == key[4])]
        a, b = set(zip(f.gene_id, f.direction)), set(zip(ref.gene_id, ref.direction))
        overlaps.append(dict(zip(['variant', 'metric', 'species', 'contrast', 'threshold'], key)) | dict(original=len(b), current=len(a), same_direction_intersection=len(a & b)))
    save(pd.DataFrame(overlaps), HERE / 'nulls/member_retention.tsv')
    associations = pd.concat([pd.read_csv(p, sep='\t') for p in (HERE / 'native').glob('*/depth_associations.tsv.gz')])
    ar = []
    for key, f in associations.groupby(['species', 'study_id', 'route', 'contrast']):
        ar.append(dict(zip(['species', 'study_id', 'route', 'contrast'], key)) | dict(n=len(f),
            spearman_dwui_depth=float(spearmanr(f.delta_WPP, f.delta_log_pas_depth).statistic)))
    save(pd.DataFrame(ar), HERE / 'nulls/depth_associations.tsv')
    print('NULLS COMPLETE', len(summary), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['axes', 'measure', 'nulls'])
    parser.add_argument('--study', choices=list(SPECS))
    parser.add_argument('--axis', default='native', choices=['native', 'common'])
    args = parser.parse_args()
    if args.command == 'axes':
        build_axes()
    elif args.command == 'measure':
        measure(args.study, args.axis)
    else:
        run_nulls()
