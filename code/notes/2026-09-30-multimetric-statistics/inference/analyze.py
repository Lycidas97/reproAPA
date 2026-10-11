"""Independent-well, fixed-date PAS-composition validation for frozen WOT states."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import platform

import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.stats import t, beta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MOUSE = Path('/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot')
STUDY = 'gse115943_pluripotency_reprogramming_core'
SOURCE = MOUSE / 'data/metacell-pas-v7-v1/studies' / STUDY
OLD = MOUSE / 'analysis/figure3-apa-gex-state-analysis-2026-09-19'
OLD_OUT = MOUSE / 'outputs/figure3-apa-gex-state-analysis-2026-09-19'
SEED, B = 20260930, 131071
INPUTS = {}


def register(path):
    path = Path(path)
    INPUTS[str(path)] = dict(sha256=hashlib.file_digest(path.open('rb'), 'sha256').hexdigest(),
                             bytes=path.stat().st_size)
    return path


def read(path, **kw):
    return pd.read_csv(register(path), sep='\t', low_memory=False, **kw)


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def dump(value, name):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    ans = np.empty(len(p))
    ans[order] = np.minimum(q, 1)
    return ans


SIGN_CACHE = {}


def signs(n):
    if n not in SIGN_CACHE:
        if n <= 18:
            ix = np.arange(1 << (n - 1), dtype=np.uint32)
            s = np.column_stack((np.ones(len(ix)),
                                 2 * ((ix[:, None] >> np.arange(n - 1)) & 1).astype(float) - 1))
            exact = True
        else:
            rng = np.random.default_rng(SEED + n)
            s = 2 * rng.integers(0, 2, size=(B, n), dtype=np.int8).astype(float) - 1
            s[:, 0] = 1
            exact = False
        SIGN_CACHE[n] = (s, exact)
    return SIGN_CACHE[n]


def sign_tests(vector, weights):
    """Omnibus max score plus columns supplied for scalar two-sided tests."""
    weighted = np.asarray(vector, float) * np.asarray(weights)[:, None]
    denom = np.sqrt(np.square(weighted).sum(0))
    standardized = np.divide(weighted, denom, out=np.zeros_like(weighted), where=denom > 1e-14)
    s, exact = signs(len(weights))
    observed = np.abs(standardized.sum(0))
    null = np.abs(s @ standardized)
    site_n = vector.shape[1] - 2
    stats = np.r_[observed[:site_n].max(), observed[-2:]]
    nulls = np.column_stack((null[:, :site_n].max(1), null[:, -2:]))
    extreme = (nulls >= stats - 1e-10).sum(0)
    p = extreme / len(s) if exact else (extreme + 1) / (len(s) + 1)
    return p, dict(exact_sign_enumeration=exact, sign_patterns=len(s),
                   omnibus_statistic=float(stats[0]), count_WPP_statistic=float(stats[1]),
                   equal_MC_WPP_statistic=float(stats[2]))


def nested_summary(values, dates):
    x = pd.DataFrame(dict(value=values, date=dates)).groupby('date').value.agg(['mean', 'var', 'count'])
    assert x['count'].min() >= 2
    k = len(x)
    components = x['var'].to_numpy() / x['count'].to_numpy() / k ** 2
    variance = components.sum()
    denom = (np.square(components) / (x['count'].to_numpy() - 1)).sum()
    df = variance ** 2 / denom if denom > 0 else np.nan
    se = np.sqrt(variance)
    estimate = float(x['mean'].mean())
    margin = t.ppf(.975, df) * se if se > 0 else 0.
    loo = (x['mean'].sum() - x['mean']) / (k - 1)
    return dict(estimate=estimate, SE=float(se), df=float(df), CI95_low=float(estimate - margin),
                CI95_high=float(estimate + margin), date_min=float(x['mean'].min()),
                date_max=float(x['mean'].max()), leave_date_min=float(loo.min()),
                leave_date_max=float(loo.max()))


def load_inputs():
    cols = ['feature_index_1based', 'atlas_pas_id', 'primary_gene_id', 'primary_gene_symbol',
            'keep_for_apa_primary', 'chrom', 'strand', 'cut_site']
    f = read(SOURCE / 'features.tsv.gz', usecols=cols)
    f = f[f.keep_for_apa_primary & f.primary_gene_id.notna()].copy()
    shape = f.groupby('primary_gene_id').agg(n=('cut_site', 'nunique'), c=('chrom', 'nunique'),
                                            s=('strand', 'nunique'))
    f = f[f.primary_gene_id.isin(shape.query('n>=2 and c==1 and s==1').index)].copy()
    z = f.cut_site * np.where(f.strand.eq('+'), 1, -1)
    lo, hi = z.groupby(f.primary_gene_id).transform('min'), z.groupby(f.primary_gene_id).transform('max')
    f['position'] = (z - lo) / (hi - lo)
    f = f.merge(read(MOUSE / 'data/pas-annotation-support-v1/studies' / STUDY /
                    'primary_pas_reference_matches.tsv.gz',
                    usecols=['atlas_pas_id', 'matched_source_count_50nt']),
                on='atlas_pas_id', validate='one_to_one')
    f = f[f.matched_source_count_50nt.ge(1)].sort_values('feature_index_1based').copy()
    # The cached state vectors preserve feature-index order within each gene.
    axes = {g: x.reset_index(drop=True) for g, x in f.groupby('primary_gene_id') if len(x) >= 2}
    parts = []
    for chunk in read(OLD / 'raw' / STUDY / 'apa_gex_states.tsv.gz', chunksize=50000):
        parts.append(chunk[chunk.setting.eq('supported50_p08') & chunk.state.isin(['P', 'M', 'S'])])
    states = pd.concat(parts, ignore_index=True)
    assert states.replicate_status.eq('biological_independent').all()
    assert not states.duplicated(['gene_id', 'library_key', 'state']).any()
    assert states.n_mc.ge(3).all() and states.n_cells.ge(30).all()
    design = read(ROOT / 'notes/2026-09-26-endpoint-counting/inputs/sampling_design.tsv')
    design = design[design.study_id.eq('GSE115943')].copy()
    detailed = read(MOUSE / 'analysis/tandem-apa-v1/library_replicate_design.tsv')
    detailed = detailed[detailed.library_key.isin(design.library_key)].copy()
    save(detailed, 'independent_well_provenance.tsv')
    states = states.merge(detailed[['library_key', 'condition', 'confidence']],
                          on='library_key', validate='many_to_one')
    for row in states.itertuples():
        counts, n_mc = np.asarray(json.loads(row.site_counts)), np.asarray(json.loads(row.site_mc))
        assert len(counts) == len(n_mc) == len(axes[row.gene_id])
        assert np.all(counts >= 0) and np.array_equal(counts, counts.astype(int))
        assert counts.sum() == row.pas_counts
    save(f, 'supported50_PAS_axes.tsv.gz')
    save(design, 'sampling_design.tsv')
    return states, axes, design


def extract_pairs(states, axes):
    keys = ['gene_id', 'gene_symbol', 'library_key', 'time', 'sample_unit', 'replicate_status',
            'condition', 'confidence']
    p = states[states.state.eq('P')].merge(states[states.state.eq('M')], on=keys,
                                          suffixes=('_P', '_M'), validate='one_to_one')
    records, vectors = [], {}
    for row in p.itertuples():
        a, b = np.asarray(json.loads(row.site_counts_P)), np.asarray(json.loads(row.site_counts_M))
        ma, mb = np.asarray(json.loads(row.site_mc_P)), np.asarray(json.loads(row.site_mc_M))
        n_joint = int(((a + b >= 10) & (ma + mb >= 2)).sum())
        if n_joint < 2:
            continue
        assert row.condition == '2i'
        ax = axes[row.gene_id]
        delta = a / a.sum() - b / b.sum()
        count_wpp = float(delta @ ax.position.to_numpy())
        equal_wpp = row.APA_P - row.APA_M
        vectors[row.gene_id, row.library_key] = np.r_[delta, count_wpp, equal_wpp]
        rec = {k: getattr(row, k) for k in keys}
        rec.update(n_joint_PAS=n_joint, n_all_supported50_PAS=len(a),
                   count_WPP_delta=count_wpp, equal_MC_WPP_delta=equal_wpp,
                   P_MC=row.n_mc_P, M_MC=row.n_mc_M, P_cells=row.n_cells_P, M_cells=row.n_cells_M,
                   P_UMI=int(a.sum()), M_UMI=int(b.sum()),
                   P_site_counts=row.site_counts_P, M_site_counts=row.site_counts_M,
                   P_site_mc=row.site_mc_P, M_site_mc=row.site_mc_M)
        records.append(rec)
    pairs = pd.DataFrame(records)
    counts = pairs.groupby(['gene_id', 'time']).library_key.transform('nunique')
    pairs['date_replicates'] = counts
    pairs['replicate_complete_date'] = counts >= 2
    save(pairs, 'qualified_well_pairs.tsv.gz')
    return pairs, vectors


def test_all(pairs, vectors, axes):
    results, coverage = [], []
    for gene, all_pairs in pairs.groupby('gene_id', sort=True):
        w = all_pairs[all_pairs.replicate_complete_date].sort_values(['time', 'library_key'])
        n, nd = len(w), w.time.nunique()
        eligible = n >= 8 and nd >= 4
        coverage.append(dict(gene_id=gene, gene_symbol=all_pairs.gene_symbol.iloc[0],
                             qualified_wells=len(all_pairs), qualified_dates=all_pairs.time.nunique(),
                             replicate_complete_wells=n, replicate_complete_dates=nd,
                             tested=eligible, exclusion_reason='none' if eligible else 'fewer_than_4_duplicate_dates'))
        if not eligible:
            continue
        data = np.stack([vectors[gene, x.library_key] for x in w.itertuples()])
        weights = 1 / nd / w.date_replicates.to_numpy()
        assert abs(weights.sum() - 1) < 1e-12
        p, stat = sign_tests(data, weights)
        rec = dict(species='mouse', study_id='GSE115943', route='WOT', contrast='P_M',
                   setting='supported50_p08', gene_id=gene, gene_symbol=w.gene_symbol.iloc[0],
                   n_wells=n, n_dates=nd, dates=';'.join(sorted(w.time.unique())),
                   libraries=';'.join(w.library_key), n_PAS=data.shape[1] - 2,
                   omnibus_p=float(p[0]), count_WPP_p=float(p[1]), equal_MC_WPP_p=float(p[2]), **stat)
        for metric, col in [('count_WPP', -2), ('equal_MC_WPP', -1)]:
            rec.update({metric + '_' + k: v for k, v in nested_summary(data[:, col], w.time).items()})
        usage_delta = (data[:, :-2] * weights[:, None]).sum(0)
        rec['MPRO_abs_date_balanced'] = float(usage_delta.max() - usage_delta.min())
        rec['L1_usage_change_date_balanced'] = float(np.abs(usage_delta).sum())
        rec['count_vs_equal_MC_WPP_direction_same'] = bool(rec['count_WPP_estimate'] * rec['equal_MC_WPP_estimate'] > 0)
        results.append(rec)
        if len(results) % 100 == 0:
            print('tested', len(results), flush=True)
    recorded = {row['gene_id'] for row in coverage}
    for gene, ax in axes.items():
        if gene not in recorded:
            coverage.append(dict(gene_id=gene, gene_symbol=ax.primary_gene_symbol.iloc[0],
                                 qualified_wells=0, qualified_dates=0, replicate_complete_wells=0,
                                 replicate_complete_dates=0, tested=False,
                                 exclusion_reason='no_cached_eligible_same_well_PM_pair'))
    out = pd.DataFrame(results)
    columns = ['omnibus_p', 'count_WPP_p', 'equal_MC_WPP_p']
    for col in columns:
        out[col[:-1] + 'q_BH_family'] = bh(out[col])
    combined = bh(out[columns].to_numpy().ravel()).reshape(len(out), len(columns))
    for j, col in enumerate(columns):
        out[col[:-1] + 'q_BH_combined'] = combined[:, j]
    save(pd.DataFrame(coverage), 'gene_coverage_and_exclusions.tsv')
    save(out, 'gene_tests.tsv.gz')
    return out


def audit_historical():
    summary = []
    for name in ['fixed_date_biological_inference.tsv.gz', 'window_biological_inference.tsv.gz']:
        d = read(OLD_OUT / name)
        w = d[d.study_id.eq('GSE115943') & d.contrast.isin(['P_M', 'P_S'])].copy()
        assert w.setting.eq('supported50_p08').all()
        w['interpretation'] = np.where(w.contrast.eq('P_M'), 'same_well_state_association',
                                       'D0_to_late_2i_endpoint_combined_contrast')
        w['pure_state_effect_adjusted_for_date_medium'] = False
        save(w, 'historical_' + name)
        for contrast, x in w.groupby('contrast'):
            summary.append(dict(source=name, contrast=contrast, tests=len(x), genes=x.gene_id.nunique(),
                                historical_APA_q05=int(x.APA_q.lt(.05).sum()),
                                note='historical family retained; no q inherited by new tests'))
    save(pd.DataFrame(summary), 'historical_inference_summary.tsv')
    return summary


def source_endpoint_tests(states, axes):
    """Independent source/target summary inference; source variance enters once."""
    records, coverage = [], []
    for gene, w in states.groupby('gene_id', sort=True):
        source = w[w.state.eq('S') & w.time.eq('D0')].sort_values('library_key')
        target = w[w.state.eq('P')].sort_values(['time', 'library_key'])
        reason = 'fewer_than_two_D0_sources'
        dates = []
        if len(source) == 2:
            sc = [np.asarray(json.loads(r.site_counts)) for r in source.itertuples()]
            sm = [np.asarray(json.loads(r.site_mc)) for r in source.itertuples()]
            for date, gg in target.groupby('time'):
                if len(gg) != 2:
                    continue
                good = True
                for row in gg.itertuples():
                    cc, cm = np.asarray(json.loads(row.site_counts)), np.asarray(json.loads(row.site_mc))
                    if any(int(((cc + c >= 10) & (cm + m >= 2)).sum()) < 2 for c, m in zip(sc, sm)):
                        good = False
                if good:
                    dates.append(gg)
            reason = 'none' if len(dates) >= 4 else 'fewer_than_four_jointly_supported_duplicate_target_dates'
        coverage.append(dict(gene_id=gene, gene_symbol=w.gene_symbol.iloc[0], D0_sources=len(source),
                             duplicate_target_dates=len(dates), tested=reason == 'none', exclusion_reason=reason))
        if reason != 'none':
            continue
        rec = dict(species='mouse', study_id='GSE115943', route='WOT', contrast='P_S',
                   setting='supported50_p08', gene_id=gene, gene_symbol=w.gene_symbol.iloc[0],
                   source_units=2, target_units=sum(map(len, dates)), n_dates=len(dates),
                   dates=';'.join(d.time.iloc[0] for d in dates),
                   interpretation='D0_to_late_2i_endpoint_combined_contrast', source_covariance='single_source_variance_component',
                   pure_state_effect_adjusted_for_date_medium=False)
        pos = axes[gene].position.to_numpy()
        def values(gg, metric):
            if metric == 'equal_MC_WPP':
                return gg.APA.to_numpy()
            c = np.stack([json.loads(x) for x in gg.site_counts])
            return c @ pos / c.sum(1)
        for metric in ['count_WPP', 'equal_MC_WPP']:
            s = values(source, metric)
            ts = [values(d, metric) for d in dates]
            means = np.r_[s.mean(), [v.mean() for v in ts]]
            components = np.r_[s.var(ddof=1) / len(s),
                               [v.var(ddof=1) / len(v) / len(dates) ** 2 for v in ts]]
            dfs = np.r_[len(s) - 1, [len(v) - 1 for v in ts]]
            variance = components.sum()
            denom = (components ** 2 / dfs).sum()
            df = variance ** 2 / denom if denom > 0 else np.nan
            se = np.sqrt(variance)
            est = float(means[1:].mean() - means[0])
            margin = float(t.ppf(.975, df) * se) if se > 0 else np.nan
            p = float(2 * t.sf(abs(est / se), df)) if se > 0 else 1.
            shifts = means[1:] - means[0]
            loo = (shifts.sum() - shifts) / (len(shifts) - 1)
            fields = dict(estimate=est, SE=float(se), df=float(df), CI95_low=est - margin,
                          CI95_high=est + margin, p=p, source_variance=float(s.var(ddof=1)),
                          inference_status='estimable' if se > 0 else 'degenerate_observed_between_well_variance',
                          source_zero_observed_variance=bool(s.var(ddof=1) == 0),
                          source_variance_fraction=float(components[0] / variance) if variance > 0 else np.nan,
                          date_min=float(shifts.min()), date_max=float(shifts.max()),
                          leave_date_min=float(loo.min()), leave_date_max=float(loo.max()))
            rec.update({metric + '_' + k: v for k, v in fields.items()})
        records.append(rec)
    recorded = {row['gene_id'] for row in coverage}
    for gene, ax in axes.items():
        if gene not in recorded:
            coverage.append(dict(gene_id=gene, gene_symbol=ax.primary_gene_symbol.iloc[0],
                                 D0_sources=0, duplicate_target_dates=0, tested=False,
                                 exclusion_reason='no_cached_eligible_P_or_S_state'))
    out = pd.DataFrame(records)
    columns = ['count_WPP_p', 'equal_MC_WPP_p']
    for col in columns:
        out[col[:-1] + 'q_BH_family'] = bh(out[col])
    combined = bh(out[columns].to_numpy().ravel()).reshape(len(out), len(columns))
    for j, col in enumerate(columns):
        out[col[:-1] + 'q_BH_combined'] = combined[:, j]
    save(out, 'source_endpoint_WPP_tests.tsv.gz')
    save(pd.DataFrame(coverage), 'source_endpoint_coverage_and_exclusions.tsv')
    return out


def design_audit(states):
    design = states[['sample_unit', 'library_key', 'time', 'state', 'condition']].drop_duplicates()
    design = design[design.state.isin(['P', 'M']) | (design.state.eq('S') & design.time.eq('D0'))]
    records = []
    for a, b in [('P', 'M'), ('P', 'S')]:
        x = design[design.state.isin([a, b])]
        base = pd.get_dummies(x[['time', 'condition']], drop_first=False, dtype=float).to_numpy()
        rank0 = np.linalg.matrix_rank(base)
        rank1 = np.linalg.matrix_rank(np.column_stack((base, x.state.eq(a).astype(float))))
        records.append(dict(contrast=a + '_' + b, observations=len(x), biological_units=x.sample_unit.nunique(),
                            base_date_medium_rank=int(rank0), state_augmented_rank=int(rank1),
                            adjusted_state_identifiable=bool(rank1 > rank0),
                            interpretation='within_well_state_effect' if a + b == 'PM' else 'endpoint_combined_state_date_medium_effect'))
    save(pd.DataFrame(records), 'design_estimability.tsv')
    return records


def main():
    for path in [Path(__file__), HERE / 'PROTOCOL.zh.md', OLD / 'extract_states.py', OLD / 'analyze.py',
                 OLD / 'window_inference.py']:
        register(path)
    states, axes, _ = load_inputs()
    print('states loaded', len(states), flush=True)
    pairs, vectors = extract_pairs(states, axes)
    print('qualified pairs', len(pairs), 'genes', pairs.gene_id.nunique(), flush=True)
    design = design_audit(states)
    history = audit_historical()
    out = test_all(pairs, vectors, axes)
    source = source_endpoint_tests(states, axes)
    dump(dict(status='COMPLETE', structural_supported50_multi_PAS_genes=len(axes), tested_genes=len(out),
              primary_omnibus_q05=int(out.omnibus_q_BH_family.lt(.05).sum()),
              omnibus_combined_q05=int(out.omnibus_q_BH_combined.lt(.05).sum()),
              count_WPP_combined_q05=int(out.count_WPP_q_BH_combined.lt(.05).sum()),
              equal_MC_WPP_combined_q05=int(out.equal_MC_WPP_q_BH_combined.lt(.05).sum()),
              source_endpoint_genes=len(source), source_count_WPP_combined_q05=int(source.count_WPP_q_BH_combined.lt(.05).sum()),
              source_equal_MC_WPP_combined_q05=int(source.equal_MC_WPP_q_BH_combined.lt(.05).sum()),
              source_zero_count_WPP_variance=int(source.count_WPP_source_zero_observed_variance.sum()),
              exact_genes=int(out.exact_sign_enumeration.sum()), randomized_genes=int((~out.exact_sign_enumeration).sum()),
              design=design, historical_summary=history, seed=SEED, random_sign_patterns=B,
              versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
              inputs=INPUTS, limitations=['Observed-state sign symmetry assumption; no randomized assignment',
                  'WOT independent well provenance moderate confidence; same cell line',
                  'No new DEXSeq fit; summary composition inference with well-level sign calibration',
                  'P_S is an endpoint combined contrast; pure adjusted state effect not identifiable']), 'audit.json')
    print('COMPLETE', len(out), flush=True)


if __name__ == '__main__':
    main()
