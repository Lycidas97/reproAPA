"""Within-observed-library conditional count tests, distinct from well inference."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import platform

import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.stats import beta, fisher_exact, random_table

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROJECTS = Path('/Users/lycidas/work/科研/projects')
MOUSE = PROJECTS / '2026-09-02_mouse-time-aware-integration-pilot'
HUMAN = PROJECTS / '2026-09-09_human-time-aware-integration'
PROFILES = MOUSE / 'data/apa-majority-shared-v1/profiles'
SEED, B = 20260930, 131071
FAMILY = ['species', 'study_id', 'route', 'setting', 'contrast', 'design',
          'library_key_a', 'library_key_b']
STATE = ['route', 'library_key', 'state', 'gene_id']
METRICS = ['DWUI', 'RWUI', 'PDUI', 'PPUI', 'MPRO']
INPUTS = {}


def read(path, **kwargs):
    path = Path(path)
    INPUTS[str(path)] = dict(sha256=hashlib.file_digest(path.open('rb'), 'sha256').hexdigest(),
                             bytes=path.stat().st_size)
    return pd.read_csv(path, sep='\t', **kwargs)


def save(frame, name):
    frame.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(q, 1)
    return result


def count_effects(a, b, axis):
    delta = a / a.sum() - b / b.sum()
    return dict(count_DWUI=float(delta @ axis.position.to_numpy()),
                count_RWUI=float(delta @ (axis['rank'] / axis.n).to_numpy()),
                count_PDUI=float(delta[axis['rank'].to_numpy().argmax()]),
                count_PPUI=float(delta[axis['rank'].to_numpy().argmin()]),
                count_MPRO=float(delta.max() - delta.min()))


def fisher(job):
    index, a, b, seed = job
    table = np.stack((a, b)).astype(np.int64)
    assert np.all(table >= 0) and np.all(table.sum(1) > 0)
    table = table[:, table.sum(0) > 0]
    assert table.shape[1] >= 2
    if table.shape[1] == 2:
        result = fisher_exact(table, alternative='two-sided')
        # For 2x2 SciPy's statistic is an odds ratio rather than table mass.
        return index, float(result.pvalue), 'Fisher_exact_2x2', 0, -1, seed, table.shape[1]
    # SciPy 1.18.1's Fisher MC wrapper incorrectly identifies a batch by
    # comparing its element count with the table's UMI total. Use its public
    # fixed-margin distribution directly; log mass also avoids underflow.
    distribution = random_table(table.sum(1), table.sum(0), seed=np.random.default_rng(seed))
    observed = float(distribution.logpmf(table))
    tolerance = np.finfo(float).eps * 100 * max(1., abs(observed))
    extreme = 0
    for start in range(0, B, 10000):
        simulated = distribution.rvs(size=min(10000, B - start))
        extreme += int((distribution.logpmf(simulated) <= observed + tolerance).sum())
    p = (extreme + 1) / (B + 1)
    assert 0 <= extreme <= B
    return index, float(p), 'Fisher_MC_2xn', B, extreme, seed, table.shape[1]


def load_states(pairs, axes):
    by_axis = {(s, g): a.reset_index(drop=True) for (s, g), a in
               axes.groupby(['study_id', 'primary_gene_id'], sort=False)}
    states, reconstruction = [], []
    for profile in sorted(PROFILES.glob('*/selected_metacells.tsv.gz')):
        identity = pd.read_csv(profile, sep='\t', nrows=1)
        subset = pairs[pairs.study_id.eq(identity.study_id.iloc[0])]
        if subset.empty:
            continue
        wanted = []
        for suffix in ['a', 'b']:
            w = subset[['route', 'library_key_' + suffix, 'gene_id', 'contrast']].copy()
            w['state'] = 'P' if suffix == 'a' else w.contrast.str[-1]
            wanted.append(w.rename(columns={'library_key_' + suffix: 'library_key'})[STATE])
        wanted = pd.concat(wanted).drop_duplicates()
        path = MOUSE / 'analysis/figure3-apa-gex-state-analysis-2026-09-19/raw' / profile.parent.name / 'apa_gex_states.tsv.gz'
        cols = STATE + ['species', 'study_id', 'setting', 'site_counts', 'site_mc', 'pas_counts', 'n_mc', 'n_cells']
        pieces = []
        for chunk in read(path, usecols=cols, chunksize=100000):
            pieces.append(chunk[chunk.setting.eq('supported50_p08')].merge(wanted, on=STATE, validate='one_to_one'))
        data = pd.concat(pieces, ignore_index=True)
        assert len(data) == len(wanted) and not data.duplicated(STATE).any()
        for row in data.itertuples():
            a = np.asarray(json.loads(row.site_counts))
            m = np.asarray(json.loads(row.site_mc))
            ax = by_axis[row.study_id, row.gene_id]
            assert len(a) == len(m) == len(ax)
            assert np.all(a >= 0) and np.array_equal(a, a.astype(np.int64))
            assert a.sum() == row.pas_counts and row.n_mc >= 3 and row.n_cells >= 30
        # Select by immutable record identity, independently of all test statistics.
        data['identity_hash'] = [hashlib.sha256(':'.join(str(getattr(r, c)) for c in STATE).encode()).hexdigest()
                                 for r in data.itertuples()]
        selected = data.sort_values('identity_hash').groupby('state', sort=True).head(3)
        meta = read(profile).sort_values('matrix_row')
        species = subset.species.iloc[0]
        source = (MOUSE / 'data/metacell-pas-v7-v1/studies' if species == 'mouse'
                  else HUMAN / 'data/results_metacell_pas_v2/studies') / profile.parent.name
        features = read(source / 'features.tsv.gz', usecols=['atlas_pas_id', 'feature_index_1based']).set_index('atlas_pas_id')
        original = read(source / 'metacells.tsv.gz', usecols=['metacell_id'])
        assert np.array_equal(original.metacell_id.iloc[meta.matrix_row], meta.metacell_id)
        matrix_path = source / 'metacell_by_pas_counts.npz'
        INPUTS[str(matrix_path)] = dict(sha256=hashlib.file_digest(matrix_path.open('rb'), 'sha256').hexdigest(),
                                       bytes=matrix_path.stat().st_size)
        matrix = sparse.load_npz(matrix_path).tocsr()
        assert np.issubdtype(matrix.data.dtype, np.integer)
        for row in selected.itertuples():
            ax = by_axis[row.study_id, row.gene_id]
            mm = meta[meta.route.eq(row.route) & meta.library_key.eq(row.library_key) &
                      meta.node.eq(row.state) & meta.purity.ge(.8)]
            ix = features.loc[ax.atlas_pas_id, 'feature_index_1based'].to_numpy() - 1
            observed = matrix[mm.matrix_row.to_numpy()][:, ix].toarray()
            good = observed.sum(1) >= 20
            observed, mm = observed[good], mm[good]
            assert np.array_equal(observed.sum(0), np.asarray(json.loads(row.site_counts)))
            assert np.array_equal((observed > 0).sum(0), np.asarray(json.loads(row.site_mc)))
            assert len(mm) == row.n_mc and mm.n_cells.sum() == row.n_cells
            reconstruction.append(dict(species=species, study_id=row.study_id, route=row.route,
                library_key=row.library_key, state=row.state, gene_id=row.gene_id,
                integer_count_exact=True, MC_support_exact=True, eligibility_exact=True))
        del matrix
        states.append(data.drop(columns='identity_hash'))
        print('COUNT_CACHE', profile.parent.name, len(data), flush=True)
    states = pd.concat(states, ignore_index=True)
    assert not states.duplicated(['study_id'] + STATE).any()
    save(states, 'integer_count_states.tsv.gz')
    save(pd.DataFrame(reconstruction), 'direct_matrix_checks.tsv')
    vectors = {(r.study_id, r.route, r.library_key, r.state, r.gene_id):
               (np.asarray(json.loads(r.site_counts), dtype=np.int64), np.asarray(json.loads(r.site_mc)))
               for r in states.itertuples()}
    return vectors, by_axis


def main():
    pairs = read(ROOT / 'pair_effects.tsv.gz', low_memory=False).reset_index(drop=True)
    assert len(pairs) == 103830 and not pairs.duplicated(FAMILY + ['gene_id']).any()
    axes = read(ROOT / 'frozen_PAS_axes.tsv.gz')
    vectors, axes = load_states(pairs, axes)
    jobs, effects = [], []
    for index, row in enumerate(pairs.itertuples()):
        a, ma = vectors[row.study_id, row.route, row.library_key_a, 'P', row.gene_id]
        b, mb = vectors[row.study_id, row.route, row.library_key_b, row.contrast[-1], row.gene_id]
        assert a.sum() == row.pas_counts_a and b.sum() == row.pas_counts_b
        assert ((a + b >= 10) & (ma + mb >= 2)).sum() == row.n_joint_PAS
        identity = json.dumps([getattr(row, c) for c in FAMILY] + [row.gene_id], separators=(',', ':'))
        seed = int.from_bytes(hashlib.sha256(identity.encode()).digest()[:8], 'little') ^ SEED
        jobs.append((index, a, b, seed))
        effects.append(count_effects(a, b, axes[row.study_id, row.gene_id]))
    out = pairs[FAMILY + ['gene_id', 'gene_symbol', 'time_a', 'time_b', 'sample_unit_a', 'sample_unit_b',
                        'replicate_status_a', 'replicate_status_b', 'pas_counts_a', 'pas_counts_b', 'n_joint_PAS']].copy()
    for metric in METRICS[:-1]:
        out['equal_MC_' + metric] = pairs['delta_' + ('WPP' if metric == 'DWUI' else metric)]
    out = pd.concat([out, pd.DataFrame(effects)], axis=1)
    completed = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        for result in executor.map(fisher, jobs, chunksize=40):
            completed.append(result)
            if len(completed) % 1000 == 0:
                print('FISHER', len(completed), '/', len(jobs), flush=True)
    tests = pd.DataFrame(completed, columns=['index', 'p', 'method', 'MC_draws', 'MC_extreme', 'MC_seed', 'positive_count_PAS'])
    assert np.array_equal(tests['index'], np.arange(len(out)))
    out = pd.concat([out, tests.drop(columns='index')], axis=1)
    mc = out.method.eq('Fisher_MC_2xn')
    alpha = .05 / int(mc.sum())
    out['p_MC_lower'], out['p_MC_upper'] = out.p.copy(), out.p.copy()
    k = out.loc[mc, 'MC_extreme'].to_numpy()
    out.loc[mc, 'p_MC_lower'] = np.where(k == 0, 0., beta.ppf(alpha / 2, k, B - k + 1))
    out.loc[mc, 'p_MC_upper'] = np.where(k == B, 1., beta.ppf(1 - alpha / 2, k + 1, B - k))
    # Intervals describe the underlying MC tail probability; include corrected point p.
    out['p_MC_lower'] = np.minimum(out.p_MC_lower, out.p)
    out['p_MC_upper'] = np.maximum(out.p_MC_upper, out.p)
    for col, name in [('p', 'q_BH'), ('p_MC_lower', 'q_BH_MC_lower'), ('p_MC_upper', 'q_BH_MC_upper')]:
        out[name] = out.groupby(FAMILY, dropna=False)[col].transform(bh)
    out['MC_classification_uncertain'] = out.q_BH_MC_lower.lt(.05) & out.q_BH_MC_upper.ge(.05)
    out['q_point_supported'] = out.q_BH.lt(.05)
    out['q_MC_upper_supported'] = out.q_BH_MC_upper.lt(.05)
    for metric in METRICS:
        for threshold in ([.1, .2] if metric == 'MPRO' else [.02, .05]):
            tag = str(threshold).replace('.', 'p')
            effect = out['count_' + metric].abs().ge(threshold)
            out[metric + '_' + tag + '_point_call'] = out.q_point_supported & effect
            out[metric + '_' + tag + '_MC_stable_call'] = out.q_MC_upper_supported & effect
    save(out, 'gene_pair_Fisher_tests.tsv.gz')
    summaries = []
    for key, group in out.groupby(['species', 'study_id', 'route', 'contrast'], sort=True):
        rec = dict(zip(['species', 'study_id', 'route', 'contrast'], key))
        rec.update(pair_gene_tests=len(group), unique_genes=group.gene_id.nunique(),
            pair_families=len(group[FAMILY].drop_duplicates()),
            exact_2x2=int(group.method.eq('Fisher_exact_2x2').sum()),
            MC_2xn=int(group.method.eq('Fisher_MC_2xn').sum()),
            q_point_supported=int(group.q_point_supported.sum()),
            q_MC_stable_supported=int(group.q_MC_upper_supported.sum()),
            MC_classification_uncertain=int(group.MC_classification_uncertain.sum()))
        for metric in METRICS:
            for threshold in ([.1, .2] if metric == 'MPRO' else [.02, .05]):
                tag = str(threshold).replace('.', 'p')
                for kind in ['point_call', 'MC_stable_call']:
                    rec[metric + '_' + tag + '_' + kind] = int(group[metric + '_' + tag + '_' + kind].sum())
        rec['count_vs_equal_MC_DWUI_sign_disagreements'] = int((group.count_DWUI * group.equal_MC_DWUI < 0).sum())
        summaries.append(rec)
    save(pd.DataFrame(summaries), 'summary_by_route.tsv')
    family_summary = out.groupby(FAMILY, dropna=False).agg(
        genes=('gene_id', 'size'), q_point_supported=('q_point_supported', 'sum'),
        q_MC_stable_supported=('q_MC_upper_supported', 'sum'),
        MC_classification_uncertain=('MC_classification_uncertain', 'sum')).reset_index()
    save(family_summary, 'FDR_families.tsv')
    audit = dict(status='PASS', tests=len(out), families=len(family_summary),
        MC_tests=int(mc.sum()), exact_2x2_tests=int((~mc).sum()), MC_draws_per_test=B,
        seed=SEED, MC_joint_interval_coverage=.95, MC_interval_alpha_per_test=alpha,
        point_q_supported=int(out.q_point_supported.sum()),
        MC_upper_q_supported=int(out.q_MC_upper_supported.sum()),
        MC_uncertain=int(out.MC_classification_uncertain.sum()),
        python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        inputs=INPUTS, scope='observed-library conditional integer-count composition; not independent biological-unit inference')
    (HERE / 'execution.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: v for k, v in audit.items() if k != 'inputs'}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
