"""Fixed-background comparison of external aUTR/cUTR ER-distribution proxy."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
OLD = PROJECT / 'notes/2026-10-10-rna-module-enrichment'


def main():
    freeze = json.loads((HERE / 'ER_proxy_freeze.json').read_text())
    for path, sha in freeze['sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha, path
    sys.path.insert(0, str(OLD))
    spec = importlib.util.spec_from_file_location('frozen_rna_go_er', OLD / 'analyze.py')
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    bg = pd.read_csv(OLD / 'ortholog_roster.tsv.gz', sep='\t')
    bg = bg[bg.contrast.eq('P_M')].reset_index(drop=True)
    assert len(bg) == 751 and bg.foreground_002.sum() == 27
    data = pd.read_csv(HERE / 'sources/GSE162971_RNAseq_counts_in_cUTRs_and_aUTRs.txt.gz', sep='\t')
    cols = ['C2C12P_Membrane_RNAseq_aUTR', 'C2C12P_Cytosol_RNAseq_aUTR',
            'C2C12P_Membrane_RNAseq_cUTR', 'C2C12P_Cytosol_RNAseq_cUTR']
    duplicate = set(data.loc[data.gene_symbol.duplicated(keep=False), 'gene_symbol'])
    data = data[~data.gene_symbol.isin(duplicate)]
    roster = bg.merge(data[['gene_symbol'] + cols], left_on='mouse_gene_symbol',
                      right_on='gene_symbol', how='left', validate='one_to_one')
    M_a, C_a, M_c, C_c = [roster[c] for c in cols]
    roster['external_cUTR_depth'] = M_c + C_c
    roster['external_aUTR_depth'] = M_a + C_a
    roster['reason'] = np.select([
        roster.mouse_gene_symbol.isin(duplicate), roster.gene_symbol.isna(),
        roster.external_cUTR_depth.lt(20), roster.external_aUTR_depth.lt(20)],
        ['duplicate_external_symbol', 'not_in_external_table', 'cUTR_lt20', 'aUTR_lt20'], default='eligible')
    roster['eligible'] = roster.reason.eq('eligible')
    roster['ER_region_proxy'] = np.where(roster.eligible,
        np.log2((M_a + .5) / (C_a + .5)) - np.log2((M_c + .5) / (C_c + .5)), np.nan)
    save_cols = ['human_gene_id', 'mouse_gene_id', 'human_gene_symbol', 'mouse_gene_symbol',
                 'foreground_002', 'foreground_005', 'eligible', 'reason'] + cols + [
                 'external_cUTR_depth', 'external_aUTR_depth', 'ER_region_proxy']
    roster[save_cols].to_csv(HERE / 'ER_proxy_roster.tsv', sep='\t', index=False)
    pairs = roster[roster.eligible].reset_index(drop=True)
    features = pd.read_csv(PROJECT / 'notes/2026-10-10-recurrence-utr-ortholog/ortholog_assay_features.tsv.gz', sep='\t')
    features = features[features.contrast.eq('P_M') & features.threshold.eq(.02)]
    _, coarse_labels = helpers.joint_codes(pairs, features, helpers.MODELS['joint_measurement_PAS_capacity'])
    external_log_depth = np.log2(pairs.external_cUTR_depth + 1)
    external_log_baseline = np.log2((pairs.external_aUTR_depth + .5) / (pairs.external_cUTR_depth + .5))
    boundaries = {}
    external_bins = []
    for name, values in [('cUTR_depth', external_log_depth), ('aUTR_cUTR_ratio', external_log_baseline)]:
        cuts = np.unique(np.quantile(values, [1 / 3, 2 / 3]))
        boundaries[name] = cuts.tolist()
        external_bins.append(np.searchsorted(cuts, values, side='right'))
    labels = {'joint_COARSE': coarse_labels,
              'joint_COARSE_external': np.array([c + '|' + str(d) + '|' + str(b)
                      for c, d, b in zip(coarse_labels, *external_bins)])}
    target = pairs.foreground_002.to_numpy(bool)
    E = pairs.ER_region_proxy.to_numpy()
    ranks = rankdata(E, method='average')
    results, margins, permutations = [], [], []
    rng = np.random.default_rng(20261010)
    for model, lab in labels.items():
        codes, _ = pd.factorize(lab, sort=True)
        _, _, binary = helpers.calculate(E > 0, target, codes)
        active = []
        for code in np.unique(codes):
            ii = np.flatnonzero(codes == code)
            n = int(target[ii].sum())
            margins.append({'model': model, 'stratum': int(code), 'label': lab[ii[0]], 'N': len(ii), 'n': n})
            if n:
                active.append((ii, n))
        mobile = sum(n for ii, n in active if n < len(ii))
        can_test = target.sum() >= 5 and mobile >= 3
        results.append({'model': model, 'test': 'E_gt0', 'status': 'tested' if can_test else 'coverage_no_test',
                        **binary, 'p': binary['p_enrichment'] if can_test else np.nan})
        obs = float(ranks[target].mean())
        sims = np.zeros(20000)
        for ii, n in active:
            if n == len(ii):
                sims += ranks[ii].sum()
            else:
                for j in range(len(sims)):
                    sims[j] += ranks[rng.choice(ii, n, replace=False)].sum()
        sims /= target.sum()
        sd = float(sims.std(ddof=1))
        p = (1 + int((sims >= obs - 1e-12).sum())) / (len(sims) + 1)
        status = 'tested' if can_test and sd > 1e-12 else 'coverage_or_frozen_no_test'
        results.append({'model': model, 'test': 'mean_rank_E', 'status': status,
                        'N': len(pairs), 'n': int(target.sum()), 'observed_mean_rank': obs,
                        'expected_mean_rank': float(sims.mean()), 'null_sd': sd, 'mobile_target': mobile,
                        'p': p if status == 'tested' else np.nan})
        permutations.extend({'model': model, 'iteration': j, 'mean_rank': float(value)} for j, value in enumerate(sims))
    result = pd.DataFrame(results)
    # Keep unavailable family members conservatively at p=1 for the fixed four-test family.
    result['q4'] = helpers.bh(result.p.fillna(1).to_numpy())
    result.to_csv(HERE / 'ER_proxy_results.tsv', sep='\t', index=False)
    pd.DataFrame(margins).to_csv(HERE / 'ER_proxy_strata.tsv', sep='\t', index=False)
    pd.DataFrame(permutations).to_csv(HERE / 'ER_proxy_permutations.tsv.gz', sep='\t', index=False)
    summary = {'N': len(pairs), 'n': int(target.sum()), 'strict_n': int(pairs.foreground_005.sum()),
               'background_E_median': float(np.median(E)), 'target_E_median': float(np.median(E[target])),
               'quantile_boundaries': boundaries, 'external_symbols': len(data),
               'duplicate_external_symbols': len(duplicate),
               'reasons': roster.reason.value_counts().to_dict(),
               'target_reasons': roster[roster.foreground_002].reason.value_counts().to_dict(),
               'protocol_sha256': hashlib.sha256((HERE / 'ER_PROXY_ADDENDUM.zh.md').read_bytes()).hexdigest(),
               'status': 'AWAITING_INDEPENDENT_QA'}
    (HERE / 'ER_proxy_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    print(result.to_string(index=False))


if __name__ == '__main__':
    main()
