"""Independent region-score replay, matched exact null, and separate permutations."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]


def main():
    checks = []

    def check(condition, label):
        assert condition, label
        checks.append(label)

    freeze = json.loads((HERE / 'ER_proxy_freeze.json').read_text())
    for path, sha in freeze['sha256'].items():
        check(hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha, 'source hash')
    data = pd.read_csv(HERE / 'sources/GSE162971_RNAseq_counts_in_cUTRs_and_aUTRs.txt.gz', sep='\t')
    roster = pd.read_csv(HERE / 'ER_proxy_roster.tsv', sep='\t')
    check(len(roster) == 751 and roster.foreground_002.sum() == 27, 'fixed selection')
    duplicate = set(data[data.gene_symbol.duplicated(False)].gene_symbol)
    data = data[~data.gene_symbol.isin(duplicate)].set_index('gene_symbol')
    values, eligible = [], []
    for row in roster.itertuples():
        gene = row.mouse_gene_symbol
        if gene in duplicate:
            reason, score = 'duplicate_external_symbol', np.nan
        elif gene not in data.index:
            reason, score = 'not_in_external_table', np.nan
        else:
            ann = data.loc[gene]
            aM, aC = ann.C2C12P_Membrane_RNAseq_aUTR, ann.C2C12P_Cytosol_RNAseq_aUTR
            cM, cC = ann.C2C12P_Membrane_RNAseq_cUTR, ann.C2C12P_Cytosol_RNAseq_cUTR
            reason = 'cUTR_lt20' if cM + cC < 20 else 'aUTR_lt20' if aM + aC < 20 else 'eligible'
            score = np.log2(((aM + .5) * (cC + .5)) / ((aC + .5) * (cM + .5))) if reason == 'eligible' else np.nan
        check(row.reason == reason, 'coverage reason')
        check((np.isnan(score) and np.isnan(row.ER_region_proxy)) or np.isclose(score, row.ER_region_proxy, atol=2e-13), 'region ratio independent')
        values.append(score)
        eligible.append(reason == 'eligible')
    part = roster[np.array(eligible)].reset_index(drop=True)
    scores = np.array(values)[eligible]
    rank = rankdata(scores)
    target = part.foreground_002.to_numpy()
    n = int(target.sum())
    result = pd.read_csv(HERE / 'ER_proxy_results.tsv', sep='\t')
    fields = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
    features = pd.read_csv(PROJECT / 'notes/2026-10-10-recurrence-utr-ortholog/ortholog_assay_features.tsv.gz', sep='\t')
    features = features[features.contrast.eq('P_M') & features.threshold.eq(.02)]
    by_species = {s: features[features.species.eq(s)].set_index('gene_id') for s in ['human', 'mouse']}
    keys = [tuple(str(by_species[sp].loc[row[sp + '_gene_id'], f])
                  for sp in ['human', 'mouse'] for f in fields) for _, row in part.iterrows()]
    cuts_depth = np.unique(np.percentile(np.log2(part.external_cUTR_depth + 1), [100 / 3, 200 / 3]))
    cuts_ratio = np.unique(np.percentile(np.log2((part.external_aUTR_depth + .5) / (part.external_cUTR_depth + .5)), [100 / 3, 200 / 3]))
    depth_bins = np.digitize(np.log2(part.external_cUTR_depth + 1), cuts_depth)
    ratio_bins = np.digitize(np.log2((part.external_aUTR_depth + .5) / (part.external_cUTR_depth + .5)), cuts_ratio)
    model_keys = {'joint_COARSE': keys,
                  'joint_COARSE_external': [k + (int(d), int(b)) for k, d, b in zip(keys, depth_bins, ratio_bins)]}
    rng = np.random.default_rng(20261011)
    comparisons = []
    for model, model_labels in model_keys.items():
        groups = {}
        for i, key in enumerate(model_labels):
            groups.setdefault(key, []).append(i)
        dist, expectation, rank_variance = np.array([1.0]), 0.0, 0.0
        simulations = np.zeros(20000)
        for indices in groups.values():
            ii = np.array(indices)
            N, chosen, K = len(ii), int(target[ii].sum()), int((scores[ii] > 0).sum())
            dist = np.convolve(dist, hypergeom.pmf(np.arange(chosen + 1), N, K, chosen))
            expectation += chosen * float(rank[ii].mean())
            if N > 1:
                rank_variance += chosen * (N - chosen) / (N - 1) * float(rank[ii].var())
            if chosen == N:
                simulations += rank[ii].sum()
            elif chosen:
                for j in range(20000):
                    simulations[j] += rank[rng.permutation(ii)[:chosen]].sum()
        expectation /= n
        rank_variance /= n ** 2
        observed = float(rank[target].mean())
        p = (1 + int((simulations / n >= observed - 1e-12).sum())) / 20001
        binary = result[result.model.eq(model) & result.test.eq('E_gt0')].iloc[0]
        continuous = result[result.model.eq(model) & result.test.eq('mean_rank_E')].iloc[0]
        hits = int(((scores > 0) & target).sum())
        check(binary.hits == hits, 'binary observed')
        check(np.isclose(binary.p, dist[hits:].sum(), atol=2e-13), 'independent exact p')
        check(np.isclose(continuous.observed_mean_rank, observed, atol=2e-13), 'mean rank observed')
        check(abs(continuous.expected_mean_rank - expectation) <= 6 * np.sqrt(rank_variance / 20000), 'MC expectation accuracy')
        combined_se = np.sqrt((p * (1 - p) + continuous.p * (1 - continuous.p)) / 20000)
        check(abs(p - continuous.p) <= 6 * combined_se + 2 / 20001, 'independent random null agreement')
        check(abs(continuous.null_sd - np.sqrt(rank_variance)) <= 6 * np.sqrt(rank_variance) / np.sqrt(40000), 'MC variance accuracy')
        comparisons.append({'model': model, 'production_p': continuous.p, 'independent_p': p,
                            'analytic_mean_rank': expectation, 'analytic_rank_sd': np.sqrt(rank_variance),
                            'joint_Monte_Carlo_se': combined_se})
    p = result.p.fillna(1).to_numpy()
    order = np.argsort(p)
    q = np.empty(4)
    q[order] = np.minimum(1, np.minimum.accumulate((p[order] * 4 / np.arange(1, 5))[::-1])[::-1])
    check(np.allclose(q, result.q4, atol=2e-13), 'fixed-family q4')
    out = {'status': 'PASS', 'checks': len(checks), 'skipped': 0, 'seed': 20261011,
           'Monte_Carlo_comparisons': comparisons,
           'validator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE / 'ER_proxy_independent_validation.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
