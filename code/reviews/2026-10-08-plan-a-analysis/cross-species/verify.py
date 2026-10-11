"""Independent input/estimate replay and exact conditional direction support."""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results'
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def direction_configurations(frame, side, hit_positions):
    labels = frame[side+'_label'].to_numpy(int)
    cols = [side+'_'+c for c in ['S_prox_capacity','M_prox_capacity','S_dist_capacity','M_dist_capacity','ref_q5','pas_q3']]
    selected = frame.loc[labels != 0]
    group_specs = []
    for _, part in selected.groupby(cols, dropna=False):
        ii = part.index.to_numpy()
        shared_bits = [j for j, pos in enumerate(hit_positions) if pos in set(ii)]
        if shared_bits:
            group_specs.append((len(ii), int((labels[ii] == -1).sum()), shared_bits))
    masks, probabilities = [], []
    for mask in itertools.product([0,1], repeat=len(hit_positions)):
        p = 1.
        for n, k, bits in group_specs:
            r, a = len(bits), sum(mask[j] for j in bits)
            if not 0 <= k-a <= n-r:
                p = 0.
                break
            p *= math.comb(n-r, k-a)/math.comb(n, k)
        if p:
            masks.append(mask)
            probabilities.append(p)
    probabilities = np.asarray(probabilities)
    assert abs(probabilities.sum()-1) < 1e-12
    return np.asarray(masks, bool), probabilities


def exact_direction_pmf(frame):
    hits = np.flatnonzero(frame.h_label.ne(0).to_numpy() & frame.m_label.ne(0).to_numpy())
    h, hp = direction_configurations(frame, 'h', hits)
    m, mp = direction_configurations(frame, 'm', hits)
    same = (h[:,None,:] == m[None,:,:]).sum(axis=2)
    prob = (hp[:,None]*mp[None,:]).ravel()
    pmf = np.bincount(same.ravel(), weights=prob, minlength=len(hits)+1)
    assert abs(pmf.sum()-1) < 1e-12
    return pmf, len(h), len(m)


def rate(k, n):
    return dict(power=k/n, power_ci_low=0 if k == 0 else beta.ppf(.025,k,n-k+1),
                power_ci_high=1 if k == n else beta.ppf(.975,k+1,n-k))


def main():
    global RESULTS
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=RESULTS)
    RESULTS = parser.parse_args().out
    audit = json.loads((RESULTS/'input-audit.json').read_text())
    for path, fingerprint in audit['inputs_sha256'].items():
        assert digest(path) == fingerprint
    assert digest(HERE/'analyze.py') == audit['code_sha256']
    assert audit['protocol_sha256'] in [digest(HERE/'PROTOCOL.initial.zh.md'),digest(HERE/'PROTOCOL.zh.md')]
    spec = importlib.util.spec_from_file_location('a3_frozen_analysis', HERE/'analyze.py')
    implementation = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(implementation)
    saved = json.loads((RESULTS/'result.json').read_text())
    for path, fingerprint in saved['outputs_sha256'].items():
        assert digest(RESULTS/path) == fingerprint
    tests = pd.read_csv(RESULTS/'raw-pvalues.tsv', sep='\t', float_precision='round_trip')
    draws = np.load(RESULTS/'simulation-draws.npz')
    powers = pd.read_csv(RESULTS/'power-curves.tsv', sep='\t', float_precision='round_trip')
    exact_rows, sensitivity, results = [], [], []
    for source in ['polyapipe','apadetect-p50']:
        f = pd.read_csv(RESULTS/(source+'-universe.tsv.gz'),sep='\t',float_precision='round_trip')
        base = ROOT/'notes/2026-10-06-dual-main-plot-data'/source/'supplement/tables'
        matched = pd.read_csv(base/'matched_pms/matched_state_vectors.tsv.gz',sep='\t')
        matched = matched[matched.setting.eq('supported50_p08')]
        candidate = pd.read_csv(base/'candidate_screen/dual_reference_candidates.tsv',sep='\t')
        orth = pd.read_csv(base/'cross_overlap/ortholog_matched_features.tsv.gz',sep='\t')
        human = set(matched.loc[matched.route.eq('GSE118258'),'gene_id'])
        mouse_rows = matched[matched.route.isin(['CIP','WOT'])]
        mouse_n = mouse_rows.groupby('gene_id').route.nunique()
        mouse = set(mouse_n[mouse_n.eq(2)].index)
        expected = orth[orth.human_gene_id.isin(human)&orth.mouse_gene_id.isin(mouse)]
        assert set(zip(f.human_gene_id,f.mouse_gene_id)) == set(zip(expected.human_gene_id,expected.mouse_gene_id))
        for species, side in [('human','h'),('mouse','m')]:
            d = candidate[candidate.species.eq(species)].set_index('gene_id').direction.map({'negative':-1,'positive':1})
            assert np.array_equal(f[species+'_gene_id'].map(d).fillna(0),f[side+'_label'])
            labels = f[side+'_label'].to_numpy(np.int8)
            cols = [side+'_'+c for c in ['S_prox_capacity','M_prox_capacity','S_dist_capacity','M_dist_capacity','ref_q5','pas_q3']]
            strata = [np.asarray(ii,int) for ii in f.groupby(cols).indices.values()]
            for conditional in [False,True]:
                use = [ii[labels[ii] != 0] for ii in strata] if conditional else strata
                use = [ii for ii in use if len(ii)]
                sample = implementation.draw_labels(labels,use,25,np.random.default_rng(9172))
                for ii in use:
                    for value in [-1,0,1]:
                        assert np.all((sample[:,ii] == value).sum(axis=1) == (labels[ii] == value).sum())
                if conditional:
                    assert np.array_equal(sample != 0,np.broadcast_to(labels != 0,sample.shape))
        hits = f.h_label.ne(0)&f.m_label.ne(0)
        observed = {'overlap_enrichment': int(hits.sum()),'direction_concordance':int((hits&f.h_label.eq(f.m_label)).sum())}
        for test, obs in observed.items():
            row = tests[tests.pas_source.eq(source)&tests.test.eq(test)].iloc[0]
            reference = draws[f'{source}|{test}|null']
            assert obs == row.observed and len(reference) == row.permutations
            p = (1+np.sum(reference >= obs))/(len(reference)+1)
            assert abs(p-row.raw_p)<1e-14
            assert abs(reference.mean()-row.null_mean)<1e-14
            for q in powers[powers.pas_source.eq(source)&powers.test.eq(test)].itertuples():
                odds = int(q.odds_ratio) if int(q.odds_ratio) == q.odds_ratio else q.odds_ratio
                sim = draws[f'{source}|{test}|OR{odds}']
                sorted_null = np.sort(reference)
                p = (1+len(reference)-np.searchsorted(sorted_null,sim,side='left'))/(len(reference)+1)
                assert np.sum(p <= q.alpha) == q.rejections
                assert len(sim) == q.simulations
        pmf, human_configs, mouse_configs = exact_direction_pmf(f)
        support = np.flatnonzero(pmf)
        t = np.arange(len(pmf))
        row = tests[tests.pas_source.eq(source)&tests.test.eq('direction_concordance')].iloc[0]
        assert abs(t@pmf-row.analytic_expected)<1e-12
        exact_rows.append(dict(pas_source=source, shared_members=int(hits.sum()),
            human_shared_direction_configurations=human_configs,mouse_shared_direction_configurations=mouse_configs,
            exact_null_expected=t@pmf,exact_support_min=int(support.min()),exact_support_max=int(support.max()),
            exact_minimum_attainable_upper_tail_p=pmf[support.max()],
            exact_upper_tail_at_observed=pmf[int(row.observed):].sum(),
            observed_same=int(row.observed)))
        pd.DataFrame(dict(same_direction=t,exact_probability=pmf)).to_csv(RESULTS/(source+'-direction-exact-pmf.tsv'),sep='\t',index=False,float_format='%.17g')
        reference = draws[f'{source}|direction_concordance|null']
        for odds in [1,1.5,2,3,4,6,8,12,16,32]:
            tilt = pmf*np.exp((t-t.max())*math.log(odds))
            tilt /= tilt.sum()
            seed = int.from_bytes(hashlib.sha256(f'20261008|{source}|direction_joint_tilt|{odds}'.encode()).digest()[:8],'little')
            rng = np.random.default_rng(seed)
            sim = rng.choice(t,size=3000,p=tilt)
            p = (1+len(reference)-np.searchsorted(np.sort(reference),sim,side='left'))/(len(reference)+1)
            for alpha in [.05,.0125]:
                sensitivity.append(dict(pas_source=source,test='direction_concordance',
                    generator='joint_exact_OR_to_T_tilt_postrun_sensitivity',odds_ratio=odds,alpha=alpha,
                    seed=seed,simulations=len(sim),simulated_mean=sim.mean(),analytic_tilt_expected=t@tilt,
                    expected_fraction=(t@tilt)/int(hits.sum()),rejections=int((p<=alpha).sum()),
                    **rate(int((p<=alpha).sum()),len(sim))))
        results.append(dict(pas_source=source,full_background_rebuilt=True,labels_replayed=True,
            all_p_and_power_replayed=True,exact_direction_expectation_verified=True,
            permutation_species_stratum_direction_margins_checked=True,
            direction_null_membership_fixed=True))
    pd.DataFrame(exact_rows).to_csv(RESULTS/'direction-exact-support.tsv',sep='\t',index=False,float_format='%.17g')
    pd.DataFrame(sensitivity).to_csv(RESULTS/'direction-power-joint-tilt-sensitivity.tsv',sep='\t',index=False,float_format='%.17g')
    report = dict(status='PASS',tests_skipped=0,checks=results,code_sha256=digest(Path(__file__)),
        updated_protocol_sha256=digest(HERE/'PROTOCOL.zh.md'),
        note='Original inference and predeclared power retained. Joint direction power is labeled postrun design sensitivity.')
    (RESULTS/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
