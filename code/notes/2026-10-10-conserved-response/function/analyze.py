"""Frozen conditional functional analysis of shared gene-level APA responses."""
import argparse
import hashlib
import json
import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = HERE.parent
TABLES = ROOT / 'notes/2026-10-10-recurrence-utr-ortholog'
GO = ROOT / 'notes/2026-10-03-polyapipe-rerun/statistics_package_v1/resources'
THRESHOLDS = [.02, .05]
TRAITS = ['mean_GEX', 'APA_coverage', 'PAS_number']
CATEGORIES = {'GO:0003723': 'molecular_function', 'GO:0006397': 'biological_process',
              'GO:0006417': 'biological_process'}
PROTOCOL_HASH = '3884344b47e085974979f6f8efa364dd92c8b403593b8bccbb739f479cd65e3c'
INPUTS = {}


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def read(path, **kwargs):
    INPUTS[str(Path(path).relative_to(ROOT))] = digest(path)
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kwargs)


def save(frame, name):
    frame.to_csv(HERE/name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def bh(p):
    p = np.asarray(p, float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind='stable')
    result = np.empty(len(p))
    result[order] = np.minimum(1, np.minimum.accumulate(
        (p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1])
    return result


def bins(values):
    """Linear tertile cuts; preserve ties and a separate missing bin."""
    out = pd.Series(-1, index=values.index, dtype=int)
    valid = values.notna()
    if valid.any():
        cuts = np.quantile(values[valid], [1/3, 2/3], method='linear')
        out.loc[valid] = np.searchsorted(cuts, values[valid], side='left')
    return out


def attach_strata(frame, naive):
    frame = frame.copy()
    for trait in TRAITS:
        ranks = pd.concat([frame[sp+'_'+trait].rank(method='average', pct=True)
                           for sp in ['human', 'mouse']], axis=1)
        missing = ranks.isna().any(axis=1)
        for side, values in [('min', ranks.min(axis=1)), ('max', ranks.max(axis=1))]:
            frame[side+'_'+trait+'_rank'] = values.mask(missing)
            frame[side+'_'+trait+'_bin'] = bins(values.mask(missing))
    common = ['measurable_mask', 'direction_class']+[f'min_{c}_bin' for c in TRAITS]
    fine = common + ([f'max_{c}_bin' for c in TRAITS] if naive else
        [f'{sp}_{co}_measured_pattern' for sp in ['human', 'mouse'] for co in ['P_S', 'P_M']])
    for model, fields in [('coarse', common), ('fine', fine)]:
        frame[model+'_stratum'] = frame[fields].astype(str).agg('|'.join, axis=1)
    assert frame.human_gene_id.is_unique and frame.mouse_gene_id.is_unique
    return frame


def classify_same_contrast(frame):
    """Unmeasured contrasts are absent; cross-contrast recurrence is not control."""
    frame = frame.copy()
    frame['case_all_four'] = False
    frame['case_first_three'] = False
    direction_columns = []
    for co in ['P_S', 'P_M']:
        hd, md = frame[co+'_human_direction'], frame[co+'_mouse_direction']
        eligible = frame[co+'_measurable']
        same = eligible & hd.isin(['proximal', 'distal']) & hd.eq(md)
        frame['case_all_four'] |= same
        frame['case_first_three'] |= same if co == 'P_S' else same & hd.eq('proximal')
        for sp in ['human', 'mouse']:
            col = co+'_'+sp+'_direction'
            frame[col] = frame[col].where(eligible)
            direction_columns.append(col)
    for sp in ['human', 'mouse']:
        frame[sp+'_recurrent_any'] = frame[[co+'_'+sp+'_direction' for co in ['P_S', 'P_M']]].isin(
            ['proximal', 'distal']).any(axis=1)
    frame['control'] = frame.human_recurrent_any.ne(frame.mouse_recurrent_any)
    directions = frame[direction_columns]
    proximal, distal = directions.eq('proximal').any(axis=1), directions.eq('distal').any(axis=1)
    frame['direction_class'] = np.select([proximal & distal, proximal, distal],
                                         ['both', 'proximal', 'distal'], default='neither')
    frame['measurable_mask'] = frame.P_S_measurable.astype(int).astype(str)+frame.P_M_measurable.astype(int).astype(str)
    frame['control_side'] = np.select([frame.control & frame.human_recurrent_any,
                                     frame.control & frame.mouse_recurrent_any], ['human', 'mouse'], default='')
    assert not (frame.case_all_four & frame.control).any()
    return frame


def build_membership():
    labels = read(TABLES/'measured_ortholog_members.tsv.gz')
    features = read(TABLES/'ortholog_assay_features.tsv.gz')
    assert labels.is_high_confidence.eq(1).all() and labels.homology_type.eq('ortholog_one2one').all()
    members, full = [], []
    for threshold in THRESHOLDS:
        lab = labels[labels.threshold.eq(threshold)]
        frame = lab[['human_gene_id', 'mouse_gene_id', 'human_gene_symbol', 'mouse_gene_symbol']].drop_duplicates().set_index('human_gene_id')
        assert frame.index.is_unique and frame.mouse_gene_id.is_unique
        for co in ['P_S', 'P_M']:
            part = lab[lab.contrast.eq(co)].set_index('human_gene_id')
            frame[co+'_measurable'] = frame.index.isin(part.index)
            for sp in ['human', 'mouse']:
                frame[co+'_'+sp+'_direction'] = part[sp+'_direction'].reindex(frame.index)
        frame = classify_same_contrast(frame)
        for sp in ['human', 'mouse']:
            for trait in TRAITS+['baseline', 'proximal_capacity', 'distal_capacity']:
                values = []
                for co in ['P_S', 'P_M']:
                    f = features[features.threshold.eq(threshold)&features.species.eq(sp)&features.contrast.eq(co)].set_index('gene_id')
                    ids = frame.index.to_series() if sp == 'human' else frame.mouse_gene_id
                    values.append(ids.map(f[trait]).where(frame[co+'_measurable']))
                frame[sp+'_'+trait] = pd.concat(values, axis=1).median(axis=1)
            for co in ['P_S', 'P_M']:
                f = features[features.threshold.eq(threshold)&features.species.eq(sp)&features.contrast.eq(co)].set_index('gene_id')
                ids = frame.index.to_series() if sp == 'human' else frame.mouse_gene_id
                frame[sp+'_'+co+'_measured_pattern'] = ids.map(f.measured_pattern).where(frame[co+'_measurable']).fillna('NE')
        frame = frame.reset_index().sort_values('human_gene_id').reset_index(drop=True)
        full.append(frame.assign(threshold=threshold))
        for set_name, field in [('main48', 'case_all_four'), ('list45', 'case_first_three')]:
            part = frame[frame[field]|frame.control].copy()
            part['case'] = part[field]
            part['role'] = np.where(part.case, 'case', 'control')
            part = attach_strata(part, False)
            members.append(part.assign(set_name=set_name, threshold=threshold))
    joint = read(ROOT/'notes/2026-10-10-proximal-module-specificity/joint_membership_and_covariates.tsv.gz', low_memory=False)
    for threshold in THRESHOLDS:
        for set_name, target in [('hN79', 'P_M_recurrent'), ('dual46', 'P_M_PS_same_contributor_dual')]:
            part = joint[joint.block.eq('naive_direction')&joint.mouse_target.eq(target)&joint.threshold.eq(threshold)].copy()
            hp, hd = part.pooled_hN_P7_proximal, part.pooled_hN_P7_distal
            mp, md = part.proximal, part.distal
            part['case'] = hp & mp
            part['control'] = (hp & ~(mp|md)) | (mp & ~(hp|hd))
            part['control_side'] = np.select([part.control & hp, part.control & mp], ['human', 'mouse'], default='')
            part = part[part.case|part.control].copy()
            part['role'] = np.where(part.case, 'case', 'control')
            part['human_gene_symbol'] = part.pooled_hN_P7_h_gene_symbol
            part['mouse_gene_symbol'] = part.m_gene_symbol
            for sp, prefix in [('human', 'h_'), ('mouse', 'm_')]:
                for new, old in [('mean_GEX', 'GEX'), ('APA_coverage', 'coverage'), ('PAS_number', 'PAS'), ('baseline', 'baseline')]:
                    part[sp+'_'+new] = part[prefix+old]
                for direction in ['proximal', 'distal']:
                    prefix_capacity = 'pooled_hN_P7_' if sp == 'human' else 'm_'
                    part[sp+'_'+direction+'_capacity'] = part[prefix_capacity+'capacity_'+direction]
            part['measurable_mask'] = 'hN_PM' if set_name == 'hN79' else 'hN_PM_PS'
            part['direction_class'] = 'proximal'
            part = attach_strata(part, True)
            members.append(part.assign(set_name=set_name, threshold=threshold))
    frame = pd.concat(members, ignore_index=True, sort=False)
    assert len(frame[frame.set_name.eq('main48')&frame.threshold.eq(.02)]) == 587
    expected = {('main48', .02): 48, ('list45', .02): 45, ('hN79', .02): 79, ('dual46', .02): 46,
                ('main48', .05): 12, ('list45', .05): 12, ('hN79', .05): 16, ('dual46', .05): 8}
    for key, value in expected.items():
        part = frame[frame.set_name.eq(key[0])&frame.threshold.eq(key[1])]
        assert int(part.case.sum()) == value, key
        assert not part.duplicated(['human_gene_id', 'mouse_gene_id']).any()
    assert frame[[sp+'_'+trait for sp in ['human', 'mouse'] for trait in TRAITS]].notna().all().all()
    save(pd.concat(full, ignore_index=True), 'all_measured_ortholog_classification.tsv.gz')
    save(frame, 'comparison_membership_features.tsv.gz')
    counts = []
    for (name, threshold), part in frame.groupby(['set_name', 'threshold'], sort=True):
        for model in ['coarse', 'fine']:
            g = part.groupby(model+'_stratum', sort=True).case.agg(['count', 'sum'])
            mobile = g['sum'].gt(0)&g['sum'].lt(g['count'])
            counts.append(dict(set_name=name, threshold=threshold, model=model, N=len(part), cases=int(part.case.sum()),
                controls=int(part.control.sum()), strata=len(g), mobile_strata=int(mobile.sum()),
                mobile_cases=int(g.loc[mobile, 'sum'].sum()), mobile_case_fraction=g.loc[mobile, 'sum'].sum()/part.case.sum(),
                expected_changed_fraction=float((2*g['sum']*(1-g['sum']/g['count'])).sum()/len(part))))
    save(pd.DataFrame(counts), 'matching_feasibility.tsv')
    return frame


def ontology(path):
    terms, aliases = {}, {}
    for block in Path(path).read_text().split('\n[Term]\n')[1:]:
        fields = defaultdict(list)
        for line in block.split('\n[')[0].splitlines():
            if ': ' in line:
                key, value = line.split(': ', 1); fields[key].append(value)
        if not fields['id'] or fields['is_obsolete'] == ['true']:
            continue
        term = fields['id'][0]
        parents = [value.split()[0] for value in fields['is_a']]
        parents += [value.split()[1] for value in fields['relationship'] if value.startswith('part_of ')]
        terms[term] = dict(name=fields['name'][0], namespace=fields['namespace'][0], parents=parents)
        for alt in fields['alt_id']:
            aliases[alt] = term
    @lru_cache(None)
    def ancestors(term):
        if term not in terms:
            return frozenset()
        return frozenset([term]).union(*(ancestors(parent) for parent in terms[term]['parents']))
    return terms, aliases, ancestors


def annotations(frame):
    INPUTS[str((GO/'go-basic.obo').relative_to(ROOT))] = digest(GO/'go-basic.obo')
    terms, aliases, ancestors = ontology(GO/'go-basic.obo')
    descendants = []
    for root, namespace in CATEGORIES.items():
        assert terms[root]['namespace'] == namespace
        descendants += [dict(root_id=root, descendant_id=term, descendant_name=record['name'], namespace=namespace)
                        for term, record in terms.items() if root in ancestors(term)]
    save(pd.DataFrame(descendants), 'root_descendants.tsv')
    result, rows = {}, []
    for sp in ['human', 'mouse']:
        raw = read(GO/f'{sp}_go_ids.tsv', header=None, names=['gene_id', 'go_id'], dtype=str, keep_default_na=False)
        ids = set(frame[sp+'_gene_id']); direct = defaultdict(set)
        for gene, term in zip(raw.gene_id, raw.go_id):
            if gene in ids:
                term = aliases.get(term, term)
                if term in terms:
                    direct[gene].add(term)
        propagated = {gene: set().union(*(ancestors(term) for term in direct.get(gene, []))) for gene in ids}
        result[sp] = propagated
        for gene in sorted(ids):
            record = dict(species=sp, gene_id=gene, direct_terms=json.dumps(sorted(direct.get(gene, []))),
                          propagated_terms=json.dumps(sorted(propagated[gene])))
            for namespace in ['molecular_function', 'biological_process']:
                record[namespace+'_known'] = any(terms[t]['namespace'] == namespace for t in direct.get(gene, []))
            for root in CATEGORIES:
                record[root] = root in propagated[gene]
            rows.append(record)
    save(pd.DataFrame(rows), 'annotation_gene_membership.tsv.gz')
    return result, terms


def logchoose(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)


@lru_cache(maxsize=100000)
def hypergeom_component(N, K, n):
    lo, hi = max(0, n+K-N), min(K, n)
    logs = np.asarray([logchoose(K, k)+logchoose(N-K, n-k)-logchoose(N, n) for k in range(lo, hi+1)])
    probabilities = np.exp(logs-logs.max()); probabilities /= probabilities.sum()
    return lo, probabilities


def exact_distribution(parameters):
    offset, probabilities = 0, np.array([1.])
    for N, K, n in parameters:
        lo, component = hypergeom_component(int(N), int(K), int(n))
        offset += lo; probabilities = np.convolve(probabilities, component)
    assert abs(probabilities.sum()-1) < 1e-10
    return offset, probabilities/probabilities.sum()


def test_category(part, known, member, model):
    qualified = part[known].copy(); qualified['member'] = np.asarray(member)[np.asarray(known)]
    params, margins = [], []
    for stratum, g in qualified.groupby(model+'_stratum', sort=True):
        N, n, K = len(g), int(g.case.sum()), int(g.member.sum())
        params.append((N, K, n))
        margins.append(dict(stratum=stratum, N=N, K=K, n=n, hits=int((g.case&g.member).sum())))
    N, n, K = len(qualified), int(qualified.case.sum()), int(qualified.member.sum())
    hits = int((qualified.case&qualified.member).sum())
    offset, pmf = exact_distribution(params)
    x = np.arange(offset, offset+len(pmf)); expected = float(x@pmf); variance = float((x-expected)**2@pmf)
    upper_tail = min(1., float(pmf[x >= hits].sum()))
    status = 'TESTED' if n and K and K < N and variance > 1e-14 else 'NE_no_conditional_variation_or_empty'
    cdf = np.cumsum(pmf)
    interval = [int(x[min(np.searchsorted(cdf, q), len(x)-1)]) for q in [.025, .975]]
    mobile = [r for r in margins if 0 < r['n'] < r['N']]
    result = dict(N_before_annotation=len(part), cases_before_annotation=int(part.case.sum()), N=N, cases=n, controls=N-n,
        annotation_unknown_cases=int((part.case&~known).sum()), annotation_unknown_controls=int((~part.case&~known).sum()),
        K=K, hits=hits, expected=expected, matched_fold=hits/expected if expected else np.nan,
        excess=hits-expected, null_variance=variance, null_p025=interval[0], null_p975=interval[1],
        p=upper_tail if status == 'TESTED' else 1., status=status, strata=len(margins),
        mobile_cases=sum(r['n'] for r in mobile), mobile_case_fraction=sum(r['n'] for r in mobile)/n if n else np.nan,
        expected_changed_fraction=sum(2*r['n']*(1-r['n']/r['N']) for r in margins)/N if N else np.nan)
    return result, qualified, margins, offset, pmf


def namespace_masks(part, annotation, terms, root, mode):
    namespace = terms[root]['namespace']
    known, member = {}, {}
    for sp in ['human', 'mouse']:
        values = [annotation[sp][g] for g in part[sp+'_gene_id']]
        known[sp] = np.asarray([any(terms[t]['namespace'] == namespace for t in value) for value in values])
        member[sp] = np.asarray([root in value for value in values])
    if mode == 'AND':
        return known['human']&known['mouse'], member['human']&member['mouse']
    return known[mode], member[mode]


def summarize_balance(qualified, model, context):
    counts = qualified.groupby(model+'_stratum').case.agg(['count', 'sum'])
    weights = qualified[model+'_stratum'].map(counts['sum']/counts['count']).to_numpy()
    rows = []
    for sp in ['human', 'mouse']:
        for cov in TRAITS+['baseline', 'proximal_capacity', 'distal_capacity']:
            vals = qualified[sp+'_'+cov].to_numpy(float); good = np.isfinite(vals)
            cases, controls = good & qualified.case.to_numpy(), good & ~qualified.case.to_numpy()
            if not good.any() or not cases.any():
                continue
            obs = float(vals[cases].mean()); ctrl = float(vals[controls].mean()) if controls.any() else np.nan
            expected = float(np.average(vals[good], weights=weights[good])) if weights[good].sum() else np.nan
            sd = float(vals[good].std(ddof=1)) if good.sum() > 1 else 0.
            rows.append(context | dict(species=sp, covariate=cov, observed_case_mean=obs, original_control_mean=ctrl,
                conditional_expected_case_mean=expected, unmatched_case_control_SMD=(obs-ctrl)/sd if sd else np.nan,
                conditional_residual_SMD=(obs-expected)/sd if sd else np.nan, finite_genes=int(good.sum())))
    return rows


def run_functional(frame):
    annotation, terms = annotations(frame)
    results, memberships, allmargins, distributions, balance = [], [], [], [], []
    for (set_name, threshold), part in frame.groupby(['set_name', 'threshold'], sort=True):
        part = part.reset_index(drop=True)
        for mode in ['human', 'mouse', 'AND']:
            for root in CATEGORIES:
                known, member = namespace_masks(part, annotation, terms, root, mode)
                for model in ['coarse', 'fine']:
                    context = dict(set_name=set_name, threshold=threshold, annotation=mode, category=root, category_name=terms[root]['name'], model=model)
                    test_id = '|'.join(map(str, context.values()))
                    row, qualified, margins, offset, pmf = test_category(part, known, member, model)
                    primary = set_name in ['main48', 'hN79'] and threshold == .02 and mode == 'human' and model == 'coarse'
                    results.append(context | row | dict(test_id=test_id, primary=primary))
                    memberships.append(part[['human_gene_id', 'mouse_gene_id', 'case', model+'_stratum']].assign(
                        set_name=set_name, threshold=threshold, annotation=mode, category=root, model=model,
                        annotation_known=known, category_member=member))
                    allmargins += [dict(test_id=test_id, **m) for m in margins]
                    distributions += [dict(test_id=test_id, hit_count=offset+i, probability=float(p)) for i,p in enumerate(pmf)]
                    balance += summarize_balance(qualified, model, context)
    result = pd.DataFrame(results)
    assert len(result) == 144 and result.primary.sum() == 6
    result['q_primary6'] = np.nan; result.loc[result.primary, 'q_primary6'] = bh(result.loc[result.primary, 'p'])
    result['q_sensitivity138'] = np.nan; result.loc[~result.primary, 'q_sensitivity138'] = bh(result.loc[~result.primary, 'p'])
    result['q_joint144'] = bh(result.p)
    save(result, 'category_results.tsv')
    save(pd.concat(memberships, ignore_index=True), 'category_membership.tsv.gz')
    save(pd.DataFrame(allmargins), 'category_stratum_margins.tsv.gz')
    save(pd.DataFrame(distributions), 'category_exact_pmf.tsv.gz')
    save(pd.DataFrame(balance), 'covariate_balance.tsv.gz')
    exploratory = []
    for set_name in ['main48', 'hN79']:
        part = frame[frame.set_name.eq(set_name)&frame.threshold.eq(.02)].reset_index(drop=True)
        for term, record in sorted(terms.items()):
            if record['namespace'] not in ['biological_process', 'molecular_function']:
                continue
            known, member = namespace_masks(part, annotation, terms, term, 'human')
            K = int((known&member).sum())
            if not 10 <= K <= 500:
                continue
            row, _, _, _, _ = test_category(part, known, member, 'coarse')
            exploratory.append(dict(set_name=set_name, term_id=term, term_name=record['name'], namespace=record['namespace'], **row))
        print('exploratory finished', set_name, flush=True)
    explore = pd.DataFrame(exploratory)
    explore['q_set'] = explore.groupby('set_name').p.transform(bh)
    explore['q_joint_scan'] = bh(explore.p)
    save(explore, 'exploratory_GO.tsv.gz')
    print(result[result.primary][['set_name', 'category', 'N', 'cases', 'hits', 'expected', 'matched_fold', 'p', 'q_primary6', 'status']].to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--stage', choices=['membership', 'functional'], required=True)
    args = parser.parse_args()
    assert digest(BASE/'PROTOCOL.zh.md') == PROTOCOL_HASH
    INPUTS[str((BASE/'PROTOCOL.zh.md').relative_to(ROOT))] = PROTOCOL_HASH
    INPUTS[str((BASE/'protocol.json').relative_to(ROOT))] = digest(BASE/'protocol.json')
    protocol = json.loads((BASE/'protocol.json').read_text())
    if args.stage == 'membership':
        frame = build_membership()
        for name, sha in INPUTS.items():
            if name in protocol['inputs_sha256']:
                assert protocol['inputs_sha256'][name] == sha, name
        outputs = ['comparison_membership_features.tsv.gz', 'all_measured_ortholog_classification.tsv.gz', 'matching_feasibility.tsv']
        status = 'MEMBERSHIP_SEALED_BEFORE_ANNOTATION'
    else:
        frame = read(HERE/'comparison_membership_features.tsv.gz', low_memory=False)
        run_functional(frame)
        outputs = ['category_results.tsv', 'category_membership.tsv.gz', 'category_stratum_margins.tsv.gz',
                   'category_exact_pmf.tsv.gz', 'covariate_balance.tsv.gz', 'exploratory_GO.tsv.gz',
                   'annotation_gene_membership.tsv.gz', 'root_descendants.tsv']
        status = 'NUMERICAL_COMPLETE_AWAITING_INDEPENDENT_QA'
    (HERE/(args.stage+'_manifest.json')).write_text(json.dumps(dict(status=status, inputs_sha256=INPUTS,
        outputs_sha256={name: digest(HERE/name) for name in outputs}, code_sha256=digest(HERE/'analyze.py'),
        numpy=np.__version__, pandas=pd.__version__), indent=2)+'\n')


if __name__ == '__main__':
    main()
