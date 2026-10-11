"""Reselect ortholog overlap from each run; retain the reviewed 20-test family."""
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

B=10000
NODE='shape_ortholog'
sha=digest

def run(ctx):
    global CTX,HERE,BASE,SOURCES,read,save
    CTX=ctx;HERE=ctx.out;BASE=ctx.analysis;SOURCES={}
    read=lambda name:ctx.read(BASE/name,low_memory=False)
    save=lambda frame,name:ctx.table(frame,name)
    calculate()
    ctx.finish(dict(fresh_profile_candidates=True,shared_joint_centers=True,twenty_test_family=True),
        selection='Candidate selection uses full data; correspondence conditional on selected ortholog pairs')

def permute_classes(frame, key, rng, b=B):
    h = frame.human_class.to_numpy(int)
    m = frame.mouse_class.to_numpy(int)
    draws = np.broadcast_to(m, (b, len(m))).copy()
    expected = np.zeros((3, 3))
    diagnostics = []
    for label, q in frame.groupby(key, sort=True):
        ii = q.index.to_numpy()
        hh, mm = h[ii], m[ii]
        draws[:, ii] = rng.permuted(np.broadcast_to(mm, (b, len(ii))), axis=1)
        for c in range(3):
            assert np.all((draws[:, ii] == c).sum(1) == (mm == c).sum())
        expected += np.outer(np.bincount(hh, minlength=3), np.bincount(mm, minlength=3)) / len(ii)
        diagnostics.append(dict(stratum=str(label), n=len(ii), mixed_mouse_classes=len(np.unique(mm)) > 1,
                                expected_changed=float(len(ii) * (1 - np.sum((np.bincount(mm, minlength=3) / len(ii)) ** 2)))))
    observed = np.zeros((3, 3), int)
    null = np.zeros((b, 3, 3), int)
    for a in range(3):
        for c in range(3):
            observed[a, c] = ((h == a) & (m == c)).sum()
            null[:, a, c] = ((h[None, :] == a) & (draws == c)).sum(1)
    assert np.all(null.sum((1, 2)) == len(frame))
    return observed, expected, null, diagnostics

def summaries(observed, expected, null, model):
    rows = []
    tests = [('diagonal', np.trace(observed), np.trace(expected), np.trace(null, axis1=1, axis2=2))]
    tests += [(f'C{a+1}_C{c+1}', observed[a, c], expected[a, c], null[:, a, c]) for a in range(3) for c in range(3)]
    for name, obs, exp, z in tests:
        p = (1 + np.sum(np.abs(z - exp) >= abs(obs - exp) - 1e-12)) / (len(z) + 1)
        rows.append(dict(model=model, test=name, observed=obs, expected=exp, excess=obs-exp,
                         null_mean=z.mean(), null_sd=z.std(ddof=1), null_p025=np.quantile(z, .025),
                         null_p975=np.quantile(z, .975), p_two_sided=p,
                         p_excess=(1 + (z >= obs).sum()) / (len(z) + 1),
                         p_depletion=(1 + (z <= obs).sum()) / (len(z) + 1)))
    return rows

def calculate():
    members = read('profiles/gene_membership_k3.tsv.gz')
    orth = read('cross_overlap/ortholog_matched_features.tsv.gz')
    nodes = read('profiles/all_route_node_positions.tsv.gz')
    fit_path = BASE / 'profiles/joint_fit.npz'
    fit = np.load(CTX.register(fit_path),allow_pickle=True)
    assert set(fit['species']) == {'human', 'mouse'}
    assert len(fit['species']) == len(members)
    assert all(np.array_equal(fit[k], members[col].to_numpy()) for k, col in [('species', 'species'), ('gene_ids', 'gene_id')])
    assert np.array_equal(fit['labels_k3'] + 1, members.cluster.str[1:].astype(int))
    SOURCES[str(fit_path)] = sha(fit_path)
    code_path = Path(__file__).with_name('downstream_profiles.py')
    SOURCES[str(code_path)] = sha(code_path)
    assert set(members.species) == {'human', 'mouse'}
    assert not members.duplicated(['species', 'gene_id']).any()
    assert orth.human_gene_id.is_unique and orth.mouse_gene_id.is_unique
    h = members[members.species.eq('human')].drop(columns='species').rename(columns={'gene_id':'human_gene_id','cluster':'human_cluster'})
    m = members[members.species.eq('mouse')].drop(columns='species').rename(columns={'gene_id':'mouse_gene_id','cluster':'mouse_cluster'})
    f = orth[['human_gene_id','mouse_gene_id','human_gene_symbol','mouse_gene_symbol']].merge(h, on='human_gene_id', validate='one_to_one').merge(m, on='mouse_gene_id', validate='one_to_one')
    f = f.sort_values('human_gene_id').reset_index(drop=True)
    for side in ['human','mouse']:
        f[side+'_class'] = f[side+'_cluster'].str[1:].astype(int)-1
        ix = pd.Index(fit['gene_ids']).get_indexer(f[side+'_gene_id'])
        assert (ix >= 0).all()
        for j, state in enumerate('SREP'):
            f[side+'_normalized_'+state] = fit['features'][ix, j]
    a = f[['human_normalized_'+x for x in 'SREP']].to_numpy()
    z = f[['mouse_normalized_'+x for x in 'SREP']].to_numpy()
    f['normalized_shape_distance'] = np.linalg.norm(a-z, axis=1)
    f['normalized_shape_cosine'] = np.sum(a*z, axis=1) / (np.linalg.norm(a,axis=1) * np.linalg.norm(z,axis=1))
    q = nodes[nodes.setting.eq('supported50_p08') & nodes.depth.eq(20) & nodes.route.isin(['CIP','WOT','GSE118258']) & nodes.node.isin(list('SREP'))].copy()
    for side, routes in [('human',['GSE118258']),('mouse',['CIP','WOT'])]:
        s = q[q.route.isin(routes)]
        complete = s.pivot(index=['route','gene_id'], columns='node', values='position_index').reindex(columns=list('SREP')).dropna()
        for route in routes:
            assert set(f[side+'_gene_id']) <= set(complete.loc[route].index)
        cov = s.groupby('gene_id')[['GEX_CPM','mean_pas_umi']].median()
        for name in ['GEX_CPM','mean_pas_umi']:
            f[side+'_'+name] = f[side+'_gene_id'].map(cov[name])
            assert f[side+'_'+name].notna().all()
            f[side+'_'+name+'_bin'] = np.minimum(np.ceil(f[side+'_'+name].rank(method='average', pct=True)*2).astype(int)-1,1)
    fields = [side+'_'+name+'_bin' for side in ['human','mouse'] for name in ['GEX_CPM','mean_pas_umi']]
    f['measurement_stratum'] = f[fields].astype(str).agg('|'.join,axis=1)
    f['unconditional_stratum'] = 'all_fixed_selected_orthologs'
    save(f,'ortholog_shape_members.tsv.gz')
    results=[]; diagnostic=[]; archive={}; cellrows=[]
    for j, (model,key) in enumerate([('selected_overlap_margins','unconditional_stratum'),('measurement_q2','measurement_stratum')]):
        observed, expected, null, diag = permute_classes(f,key,np.random.default_rng(2026100401+j))
        results += summaries(observed,expected,null,model)
        diagnostic += [dict(model=model,**d) for d in diag]
        archive[model] = null
        for a in range(3):
            for c in range(3):
                cellrows.append(dict(model=model,human_cluster=f'C{a+1}',mouse_cluster=f'C{c+1}',observed=observed[a,c],expected=expected[a,c]))
    results=pd.DataFrame(results); assert len(results)==20
    results['q_exploratory_20']=bh(results.p_two_sided)
    save(results,'shape_tests.tsv');save(cellrows,'shape_contingency_3x3.tsv');save(diagnostic,'null_mobility.tsv')
    np.savez_compressed(HERE/'null_draws.npz',**archive)
    checks=dict(gene_not_gene_route_identity=True,shared_joint_k3_centers_verified=True,selection_and_measurability_fixed=True,
                complete_profiles_in_all_fit_routes=True,each_stratum_class_margins_preserved=True,twenty_test_family_complete=True)
    (HERE/'PROVENANCE.json').write_text(json.dumps(dict(sources=SOURCES,checks=checks,n_selected_human=int(members.species.eq('human').sum()),n_selected_mouse=int(members.species.eq('mouse').sum()),n_common_selected_orthologs=len(f),permutations=B),indent=2))
    print(results[results.test.eq('diagonal')].to_string(index=False))

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,NODE,a.out))
