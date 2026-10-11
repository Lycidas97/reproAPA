"""Fresh PMS/PEM shared-axis estimates, paired deltas and named conditional nulls."""
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

from downstream_axis_extensions import weighted_rank_corr
B=3000
N_NULL=5000
SEED=20261008
NODE='same_pas'
ROUTES={'CIP':'C','WOT':'W','GSE118258':'H1','GSE150311':'H2'}

def export_pem(ctx):
    import gc
    inventory=[];rows=[]
    for folder in ['GSE118258','GSE150311','GSE221739','GSE103221','WOT_GSE115943']:
        meta=ctx.read(ctx.native/'studies'/folder/'selected_metacells.tsv.gz')
        z=meta[meta.purity.ge(.8)].groupby(['route','library_key','node']).agg(n_mc=('metacell_id','size'),n_cells=('n_cells','sum')).reset_index()
        z['study_folder']=folder;z['pre_gene_eligible']=z.n_mc.ge(3)&z.n_cells.ge(30);inventory.append(z)
    inv=pd.concat(inventory,ignore_index=True);ctx.table(inv,'frozen_library_state_inventory.tsv')
    mapping=ctx.read(ctx.analysis/'cis/detected_to_annotation.tsv.gz')
    for sid in ['GSE150311','GSE103221']:
        study=Study(ctx,sid);q=inv[inv.study_folder.eq(sid)&inv.pre_gene_eligible]
        eligible=q.pivot(index=['route','library_key'],columns='node',values='n_mc').reindex(columns=['P','M','E']).dropna().index
        code='H2' if sid=='GSE150311' else 'C'
        for gene in sorted(set(mapping.loc[mapping.code.eq(code),'gene_id'])&set(study.genes)):
            axis=study.axis(gene)
            if len(axis)<2:continue
            data=study.gene_counts(axis)
            for route,lib in eligible:
                for state in ['P','M','E']:
                    s=study.state(data,route,lib,state)
                    if not s['endpoint_eligible']:continue
                    d=study.design[study.design.library_key.eq(lib)].iloc[0]
                    rows.append(dict(code=code,route=route,study_id=sid,gene_id=gene,library_key=lib,state=state,time=d.time,sample_unit=d.sample_unit,
                        n_mc=s['n_mc'],n_cells=s['n_cells'],usage=json.dumps(s['usage'].tolist()),site_counts=json.dumps(s['site_counts'].tolist()),
                        site_mc=json.dumps(s['site_mc'].tolist()),PAS_ids=json.dumps(axis.atlas_pas_id.tolist()),native=float(s['usage']@axis.position)))
        del study;gc.collect()
    ctx.table(rows,'PEM_library_usage.tsv',columns=['code','route','study_id','gene_id','library_key','state','time','sample_unit','n_mc','n_cells','usage','site_counts','site_mc','PAS_ids','native'])

def run(ctx):
    global CTX,HERE,read,write
    CTX=ctx;HERE=ctx.out;read=ctx.read;write=lambda frame,name:ctx.table(frame,name)
    export_pem(ctx)
    variant=ctx.config['reference_variant']
    summaries(build_axes(variant),variant)
    capacity(variant)
    ctx.finish(dict(fresh_PMS_and_PEM=True,paired_differences_direct=True,rank_only_null_q=True,
                    eight_test_family_per_tolerance=True,coarse_and_capacity_nulls_separate=True),
               bootstrap=3000,null_draws=5000,bootstrap_unit='ortholog gene, synchronized over fixed study pairs')

def usage_inputs(variant):
    cis = CTX.analysis / 'cis'
    mapping = read(cis / 'detected_to_annotation.tsv.gz')
    maps = {c: q.set_index('native_pas_id').id.to_dict() for c, q in mapping.groupby('code')}
    pms = read(CTX.analysis / 'matched_pms/matched_state_vectors.tsv.gz')
    pms = pms[pms.setting.eq('supported50_p08')]
    records = []
    for r in pms.itertuples():
        code = ROUTES.get(r.route)
        if not code:
            continue
        ids = json.loads(r.PAS_ids)
        lookup = maps[code]
        if not all(p in lookup for p in ids):
            assert not any(p in lookup for p in ids)
            continue
        frame = pd.DataFrame(dict(id=[lookup[p] for p in ids],
            P=json.loads(r.P_usage), M=json.loads(r.M_usage), S=json.loads(r.S_usage))).groupby('id').sum()
        assert np.allclose(frame.sum(), 1, atol=1e-12, rtol=0)
        records.extend(dict(code=code, gene_id=r.gene_id, id=idx, design='PMS',
            P=z.P, M=z.M, S=z.S, native_P=r.P_DWUI, native_M=r.M_DWUI, native_S=r.S_DWUI,
            n_libraries=r.n_target_libraries) for idx, z in frame.iterrows())
    controls = read(HERE / 'PEM_library_usage.tsv')
    eligible_records, library_rows = [], []
    for (code, gene), g in controls.groupby(['code', 'gene_id']):
        good = []
        for lib, q in g.groupby('library_key'):
            if set(q.state) != {'P', 'E', 'M'}:
                continue
            assert q.PAS_ids.nunique() == 1
            q = q.set_index('state')
            counts = {s: np.array(json.loads(q.loc[s, 'site_counts'])) for s in ['P', 'E', 'M']}
            mc = {s: np.array(json.loads(q.loc[s, 'site_mc'])) for s in ['P', 'E', 'M']}
            joint = {co: int(((counts[a] + counts[b] >= 10) & (mc[a] + mc[b] >= 2)).sum())
                     for co, a, b in [('PM', 'P', 'M'), ('EM', 'E', 'M')]}
            library_rows.append(dict(variant=variant, code=code, gene_id=gene, library_key=lib,
                PM_joint_PAS=joint['PM'], EM_joint_PAS=joint['EM'], passed=min(joint.values()) >= 2))
            if min(joint.values()) >= 2:
                good.append(q.reset_index())
        if not good:
            continue
        z = pd.concat(good).pivot(index=['library_key', 'time', 'sample_unit'], columns='state', values='usage')
        # Equal date, then sample unit, then capture; one weight vector for all three states.
        keys = z.reset_index()[['library_key', 'time', 'sample_unit']]
        nw = keys.groupby(['time', 'sample_unit']).library_key.transform('size').to_numpy(float)
        units = keys[['time', 'sample_unit']].drop_duplicates().groupby('time').size()
        weights = 1 / nw / keys.time.map(units).to_numpy() / keys.time.nunique()
        assert abs(weights.sum() - 1) < 1e-12
        ids = json.loads(g.PAS_ids.iloc[0])
        assert g.PAS_ids.nunique() == 1
        lookup = maps[code]
        assert all(p in lookup for p in ids)
        out = pd.DataFrame({'id': [lookup[p] for p in ids]})
        for state in ['P', 'E', 'M']:
            out[state] = weights @ np.stack(z[state].map(json.loads))
        out = out.groupby('id').sum()
        assert np.allclose(out.sum(), 1, atol=1e-12, rtol=0)
        native = pd.concat(good).pivot(index=['library_key', 'time', 'sample_unit'], columns='state', values='native').reindex(z.index)
        native_mean = {s: float(weights @ native[s]) for s in ['P', 'E', 'M']}
        eligible_records.extend(dict(code=code, gene_id=gene, id=idx, design='PEM', P=r.P, E=r.E, M=r.M,
            native_P=native_mean['P'], native_E=native_mean['E'], native_M=native_mean['M'],
            n_libraries=len(z)) for idx, r in out.iterrows())
    write(library_rows, variant + '_control_library_eligibility.tsv.gz')
    usage = pd.concat([pd.DataFrame(records), pd.DataFrame(eligible_records)], ignore_index=True)
    write(usage, variant + '_annotated_usage.tsv.gz')
    return usage, read(cis / 'reciprocal_one_to_one.tsv.gz'), read(cis / 'annotation_catalog.tsv.gz')

def build_axes(variant):
    usage, orth, cat = usage_inputs(variant)
    cat = cat.set_index('id')
    genes, sites = [], []
    for tolerance in [25, 50]:
        pairs = orth[orth['reciprocal' + str(tolerance)]].copy()
        for side in ['human', 'mouse']:
            pairs[side + '_cut'] = pairs[side + '_id'].map(cat.cut_site) * pairs[side + '_id'].map(cat.strand).map({'+': 1, '-': -1})
        for design, a in usage.groupby('design'):
            states = ['P', 'M', 'S'] if design == 'PMS' else ['P', 'E', 'M']
            for h in sorted(set(a.code) - {'C', 'W'}):
                for m in ['C', 'W']:
                    human = a[a.code.eq(h)].drop(columns=['gene_id', 'code', 'design'])
                    mouse = a[a.code.eq(m)].drop(columns=['gene_id', 'code', 'design'])
                    z = pairs.merge(human, left_on='human_id', right_on='id', validate='one_to_one').drop(columns='id')
                    z = z.merge(mouse, left_on='mouse_id', right_on='id', suffixes=('_h', '_m'), validate='one_to_one').drop(columns='id')
                    for gid, q in z.groupby('human_gene_id'):
                        q = q.sort_values('human_cut')
                        if len(q) < 2 or not np.all(np.diff(q.mouse_cut) > 0):
                            continue
                        masses = {s + '_' + side: q[s + '_' + side].sum() for s in states for side in ['h', 'm']}
                        mass = min(masses.values())
                        if mass <= 0:
                            continue
                        rec = dict(variant=variant, tolerance=tolerance, design=design, pair=h + '_' + m,
                            human_gene_id=gid, mouse_gene_id=q.mouse_gene_id.iloc[0], n_matched=len(q), mass_min=mass,
                            h_n_libraries=q.n_libraries_h.iloc[0], m_n_libraries=q.n_libraries_m.iloc[0])
                        rank = np.linspace(0, 1, len(q))
                        for side in ['h', 'm']:
                            coord = q['human_cut' if side == 'h' else 'mouse_cut'].to_numpy()
                            distance = (coord - coord.min()) / np.ptp(coord)
                            for state in states:
                                value = q[state + '_' + side].to_numpy()
                                norm = value / value.sum() if value.sum() > 0 else np.full(len(q), np.nan)
                                rec['rank_' + state + '_' + side] = norm @ rank
                                rec['distance_' + state + '_' + side] = norm @ distance
                                rec['native_' + state + '_' + side] = q['native_' + state + '_' + side].iloc[0]
                        for axis in ['rank', 'distance', 'native']:
                            for side in ['h', 'm']:
                                cs = [('PS', 'P', 'S'), ('MS', 'M', 'S'), ('PM', 'P', 'M')] if design == 'PMS' else [('PM', 'P', 'M'), ('EM', 'E', 'M')]
                                for co, x, y in cs:
                                    rec[axis + '_' + co + '_' + side] = rec[axis + '_' + x + '_' + side] - rec[axis + '_' + y + '_' + side]
                                if design == 'PMS':
                                    assert abs(rec[axis + '_PM_' + side] - rec[axis + '_PS_' + side] + rec[axis + '_MS_' + side]) < 1e-12
                        genes.append(rec)
                        qq = q[['human_gene_id', 'mouse_gene_id', 'human_id', 'mouse_id']].copy()
                        qq['variant'] = variant; qq['tolerance'] = tolerance; qq['design'] = design; qq['pair'] = h + '_' + m
                        sites.append(qq)
    g = pd.DataFrame(genes)
    write(g, variant + '_same_axis_gene_effects.tsv.gz')
    write(pd.concat(sites), variant + '_matched_physical_sites.tsv.gz')
    return g

def covariance_terms(x, y, w):
    # x/y: columns PS, MS. PM = PS - MS exactly. Population weighting is consistent throughout.
    xc = x - np.average(x, axis=0, weights=w)
    yc = y - np.average(y, axis=0, weights=w)
    cov = (xc * w[:, None]).T @ yc / w.sum()
    pmx, pmy = x[:, 0] - x[:, 1], y[:, 0] - y[:, 1]
    den = np.sqrt(np.average((pmx - np.average(pmx, weights=w))**2, weights=w) *
                  np.average((pmy - np.average(pmy, weights=w))**2, weights=w))
    terms = np.array([cov[0, 0], cov[1, 1], -cov[0, 1], -cov[1, 0]])
    pearson = terms.sum() / den if den > 0 else np.nan
    return np.r_[terms, den, terms / den if den > 0 else np.full(4, np.nan), pearson]

def contrast_metrics(x, y, w, design):
    r = np.array([weighted_rank_corr(x[:, k], y[:, k], w) for k in range(x.shape[1])])
    return np.r_[r, r[2] - r[0], r[2] - r[1]] if design == 'PMS' else np.r_[r, r[0] - r[1]]

def summaries(genes, variant):
    rng = np.random.default_rng(SEED)
    rows, covrows, tests, null_audits = [], [], [], []
    for tol in [25, 50]:
        for design in ['PMS', 'PEM']:
            z = genes[genes.tolerance.eq(tol) & genes.design.eq(design) & genes.mass_min.ge(.8)].copy()
            if z.empty:
                continue
            contrasts = ['PS', 'MS', 'PM'] if design == 'PMS' else ['PM', 'EM']
            labels = [co + '_rho' for co in contrasts] + (['PM_minus_PS_rho', 'PM_minus_MS_rho'] if design == 'PMS' else ['PM_minus_EM_rho'])
            ids = sorted(z.human_gene_id.unique()); lookup = {g: i for i, g in enumerate(ids)}
            bundles = []
            for pair, q in z.groupby('pair'):
                q = q.sort_values('human_gene_id')
                if len(q) < 10:
                    continue
                ix = np.array([lookup[g] for g in q.human_gene_id])
                bundles.append((pair, q, ix))
            if not bundles:
                continue
            for axis in ['rank', 'distance', 'native']:
                values = np.empty((B + 1, len(bundles), len(labels)))
                covvals = np.empty((B + 1, len(bundles), 10)) if design == 'PMS' else None
                for draw in range(B + 1):
                    counts = np.ones(len(ids)) if draw == 0 else rng.multinomial(len(ids), np.ones(len(ids)) / len(ids))
                    for j, (_, q, ix) in enumerate(bundles):
                        x = q[[axis + '_' + co + '_h' for co in contrasts]].to_numpy()
                        y = q[[axis + '_' + co + '_m' for co in contrasts]].to_numpy()
                        values[draw, j] = contrast_metrics(x, y, counts[ix], design)
                        if design == 'PMS':
                            covvals[draw, j] = covariance_terms(x[:, :2], y[:, :2], counts[ix])
                for j, (pair, q, _) in enumerate(bundles + [('pair_mean', z, None)]):
                    v = values[:, j] if j < len(bundles) else values.mean(axis=1)
                    for k, label in enumerate(labels):
                        low, high = np.nanquantile(v[1:, k], [.025, .975])
                        rows.append(dict(variant=variant, tolerance=tol, design=design, axis=axis, pair=pair,
                            n_genes=q.human_gene_id.nunique(), metric=label, estimate=v[0, k], low=low, high=high, bootstrap=B))
                    if design == 'PMS':
                        cv = covvals[:, j] if j < len(bundles) else covvals.mean(axis=1)
                        cn = ['cov_PS_PS', 'cov_MS_MS', 'minus_cov_PS_MS', 'minus_cov_MS_PS', 'PM_sd_product',
                              'normalized_PS_PS', 'normalized_MS_MS', 'normalized_minus_PS_MS', 'normalized_minus_MS_PS', 'PM_pearson']
                        for k, label in enumerate(cn):
                            low, high = np.nanquantile(cv[1:, k], [.025, .975])
                            covrows.append(dict(variant=variant, tolerance=tol, design=design, axis=axis, pair=pair,
                                n_genes=q.human_gene_id.nunique(), component=label, estimate=cv[0, k], low=low, high=high))
                if axis != 'rank':
                    continue
                # Keep each mouse ortholog's whole contrast vector and pair-availability pattern together.
                # Null is conditional on observed studies, number of matched sites and minimum retained mass.
                features = z.groupby('human_gene_id').agg(n_sites=('n_matched', 'median'), mass=('mass_min', 'min'))
                signatures = z.groupby('human_gene_id').pair.agg(lambda x: '|'.join(sorted(x)))
                strata = signatures + '|sites' + features.n_sites.ge(3).astype(int).astype(str) + '|mass' + features.mass.ge(.9).astype(int).astype(str)
                pools = [np.array([lookup[i] for i in q.index]) for _, q in strata.groupby(strata)]
                matrices = []
                for pair, q, ix in bundles:
                    x = q[['rank_' + co + '_h' for co in contrasts]].to_numpy()
                    y = q[['rank_' + co + '_m' for co in contrasts]].to_numpy()
                    yi = np.full((len(ids), len(contrasts)), np.nan); yi[ix] = y
                    matrices.append((x, yi, ix))
                null = np.empty((N_NULL, len(labels)))
                for draw in range(N_NULL):
                    perm = np.arange(len(ids))
                    for pool in pools:
                        perm[pool] = rng.permutation(pool)
                    estimates = []
                    for x, y, ix in matrices:
                        yy = y[perm[ix]]
                        assert np.isfinite(yy).all()
                        estimates.append(contrast_metrics(x, yy, np.ones(len(ix)), design))
                    null[draw] = np.mean(estimates, axis=0)
                observed = values[0].mean(axis=0)
                for k, label in enumerate(labels):
                    p = (1 + int(np.sum(null[:, k] >= observed[k] - 1e-12))) / (N_NULL + 1)
                    tests.append(dict(variant=variant, tolerance=tol, design=design, metric=label, estimate=observed[k],
                        null_mean=null[:, k].mean(), null_low=np.quantile(null[:, k], .025), null_high=np.quantile(null[:, k], .975),
                        p_one_sided=p, n_genes=len(ids), null_draws=N_NULL))
                np.savez_compressed(HERE / (variant + '_' + design + '_' + str(tol) + '_null.npz'), null=null, labels=np.array(labels), observed=observed)
                null_audits.append(dict(variant=variant, tolerance=tol, design=design, n_strata=len(pools),
                    n_singletons=sum(len(p) == 1 for p in pools), singleton_genes=sum(len(p) for p in pools if len(p) == 1),
                    n_genes=len(ids), seed=SEED))
    t = pd.DataFrame(tests)
    assert len(t) == 16, 'Same-PAS tests need all 8 planned estimable tests at each tolerance; inspect eligibility outputs'
    for tol, ix in t.groupby('tolerance').groups.items():
        t.loc[ix, 'q_BH_family'] = bh(t.loc[ix, 'p_one_sided'])
        t.loc[ix, 'family_size'] = len(ix)
    write(rows, variant + '_correlation_summary.tsv')
    write(covrows, variant + '_covariance_decomposition.tsv')
    write(t, variant + '_conditional_null_tests.tsv')
    write(null_audits, variant + '_null_strata_audit.tsv')

def capacity(variant):
    SEED=20261009;B=5000
    results, audits = [], []
    for variant in [variant]:
        genes = read(HERE / (variant + '_same_axis_gene_effects.tsv.gz'))
        for design in ['PMS', 'PEM']:
            z = genes[genes.tolerance.eq(25) & genes.mass_min.ge(.8) & genes.design.eq(design)]
            contrasts = ['PS', 'MS', 'PM'] if design == 'PMS' else ['PM', 'EM']
            labels = [x + '_rho' for x in contrasts] + (['PM_minus_PS_rho', 'PM_minus_MS_rho'] if design == 'PMS' else ['PM_minus_EM_rho'])
            ids = sorted(z.human_gene_id.unique()); index = {g: i for i, g in enumerate(ids)}
            q = z.groupby('human_gene_id').agg(n_sites=('n_matched', 'median'), mass=('mass_min', 'min'))
            signatures = z.groupby('human_gene_id').pair.agg(lambda x: '|'.join(sorted(x)))
            strata = signatures + '|sites' + q.n_sites.ge(3).astype(int).astype(str) + '|mass' + q.mass.ge(.9).astype(int).astype(str)
            # Fixed 0.5 split of human and mouse reference positions. All components
            # are exchanged jointly; do not select a different best stratum model.
            references = ['S', 'M'] if design == 'PMS' else ['M']
            for state in references:
                for side in ['h', 'm']:
                    baseline = z.groupby('human_gene_id')['rank_' + state + '_' + side].mean()
                    strata += '|' + state + side + baseline.ge(.5).astype(int).astype(str)
            pools = [np.array([index[g] for g in a.index]) for _, a in strata.groupby(strata)]
            bundles = []
            for pair, g in z.groupby('pair'):
                g = g.sort_values('human_gene_id')
                ix = np.array([index[x] for x in g.human_gene_id])
                x = g[['rank_' + co + '_h' for co in contrasts]].to_numpy()
                y = np.full((len(ids), len(contrasts)), np.nan)
                y[ix] = g[['rank_' + co + '_m' for co in contrasts]].to_numpy()
                bundles.append((x, y, ix))
            rng = np.random.default_rng(SEED)
            null = np.empty((B, len(labels)))
            observed = np.mean([contrast_metrics(x, y[ix], np.ones(len(ix)), design) for x, y, ix in bundles], axis=0)
            for b in range(B):
                perm = np.arange(len(ids))
                for pool in pools:
                    perm[pool] = rng.permutation(pool)
                vals = []
                for x, y, ix in bundles:
                    yy = y[perm[ix]]
                    assert np.isfinite(yy).all()
                    vals.append(contrast_metrics(x, yy, np.ones(len(ix)), design))
                null[b] = np.mean(vals, axis=0)
            for k, label in enumerate(labels):
                results.append(dict(variant=variant, design=design, metric=label, estimate=observed[k],
                    null_mean=null[:, k].mean(), p_one_sided=(1 + np.sum(null[:, k] >= observed[k] - 1e-12)) / (B + 1)))
            np.savez_compressed(HERE / f'{variant}_{design}_capacity_null.npz', null=null, observed=observed, labels=np.array(labels))
            audits.append(dict(variant=variant, design=design, n_genes=len(ids), n_strata=len(pools),
                singleton_genes=sum(len(x) for x in pools if len(x) == 1), movable_genes=sum(len(x) for x in pools if len(x) > 1)))
    t = pd.DataFrame(results)
    for variant, ix in t.groupby('variant').groups.items():
        assert len(ix) == 8
        t.loc[ix, 'q_BH_family8'] = bh(t.loc[ix, 'p_one_sided'])
    write(t, 'capacity_sensitivity_tests.tsv')
    write(audits, 'capacity_sensitivity_strata.tsv')

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,NODE,a.out))
