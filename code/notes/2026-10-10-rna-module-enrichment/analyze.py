"""Prespecified RNA categories on fixed measured one-to-one ortholog inventories."""
import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from exact_null import bh, distribution, upper_tail, rejection_probability, detectable_effect

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = Path('/Users/lycidas/work/科研/projects/2026-09-26_reproapa-go-analysis/sources')
PREVIOUS = ROOT / 'notes/2026-10-10-recurrence-utr-ortholog'
CATEGORIES = {
    'RNA_binding': ('GO:0003723',),
    'RNA_processing': ('GO:0006396',),
    'RNA_splicing': ('GO:0008380',),
    'CPA_3prime_end': ('GO:0031124', 'GO:0180010', 'GO:0005847', 'GO:0005848', 'GO:0005849'),
    'mRNA_stability_regulation': ('GO:0043488',),
    'miRNA_processing_silencing': ('GO:0035195', 'GO:0035196'),
    'translation_regulation': ('GO:0006417',),
}
COARSE = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
FULL = COARSE + ['mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']
MODELS = {'joint_measurement_PAS_capacity': COARSE, 'joint_full_assay': FULL}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(frame, filename):
    compression = {'method': 'gzip', 'mtime': 0} if filename.endswith('.gz') else None
    frame.to_csv(HERE / filename, sep='\t', index=False, na_rep='NA',
                 float_format='%.17g', compression=compression)


def ontology(path):
    terms, aliases, obsolete = {}, {}, set()
    text = path.read_text()
    version = next(s.split(': ', 1)[1] for s in text.splitlines() if s.startswith('data-version: '))
    for block in text.split('\n[Term]\n')[1:]:
        record = defaultdict(list)
        for line in block.split('\n[')[0].splitlines():
            if ': ' in line:
                key, value = line.split(': ', 1)
                record[key].append(value)
        if not record['id']:
            continue
        term = record['id'][0]
        if record['is_obsolete'] == ['true']:
            obsolete.add(term)
            continue
        parents = [s.split()[0] for s in record['is_a']]
        parents += [s.split()[1] for s in record['relationship'] if s.startswith('part_of ')]
        terms[term] = dict(name=record['name'][0], namespace=record['namespace'][0], parents=parents)
        for alias in record['alt_id']:
            aliases[alias] = term

    @lru_cache(None)
    def closure(term):
        if term not in terms:
            return frozenset()
        return frozenset([term]).union(*(closure(parent) for parent in terms[term]['parents']))

    for roots in CATEGORIES.values():
        assert all(term in terms for term in roots), 'Prespecified category uses inactive GO term'
    return terms, aliases, obsolete, closure, version


def annotations(path, wanted, terms, aliases, obsolete, closure):
    direct, active, ignored = defaultdict(set), defaultdict(set), defaultdict(set)
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        if len(fields) < 2:
            continue
        gene, term = fields[:2]
        gene = gene.split('.')[0]
        if gene not in wanted or not term:
            continue
        term = aliases.get(term, term)
        if term not in terms:
            ignored[gene].add(term)
            continue
        direct[gene].add(term)
        active[gene].update(closure(term))
    return direct, active, ignored


def joint_codes(pairs, features, columns):
    parts = []
    for species in ['human', 'mouse']:
        part = features[features.species.eq(species)].set_index('gene_id')
        ids = pairs[species + '_gene_id']
        ordered = part.reindex(ids)[columns]
        assert len(ordered) == len(pairs) and not ordered.isna().any().any()
        parts.append(ordered.astype(str).agg('|'.join, axis=1).to_numpy())
    labels = np.array([h + '||' + m for h, m in zip(*parts)])
    return pd.factorize(labels, sort=True)[0], labels


def calculate(category_membership, foreground, codes):
    N = np.bincount(codes)
    K = np.bincount(codes, weights=category_membership.astype(int)).astype(int)
    n = np.bincount(codes, weights=foreground.astype(int)).astype(int)
    params = list(zip(N.tolist(), K.tolist(), n.tolist()))
    dist = distribution(params)
    hits = int((foreground & category_membership).sum())
    lo, probabilities = dist
    mean = float(np.arange(lo, lo + len(probabilities)) @ probabilities)
    variance = float(((np.arange(lo, lo + len(probabilities)) - mean) ** 2) @ probabilities)
    active = (n > 0) & (n < N)
    return params, dist, dict(
        N=int(N.sum()), n=int(n.sum()), K=int(K.sum()), hits=hits,
        expected=mean, fold=hits / mean if mean else np.nan, excess=hits - mean,
        null_sd=float(np.sqrt(variance)), p_enrichment=upper_tail(dist, hits),
        support_min=lo, support_max=lo + len(probabilities) - 1,
        inferentially_mobile=variance > 1e-15,
        strata=len(N), singleton_genes=int(N[N == 1].sum()),
        foreground_mobile=int(n[active].sum()), foreground_frozen=int(n[~active].sum()),
        expected_changed_foreground=float((n * (1 - n / N)).sum()),
        odds_ratio=float((hits * (int(N.sum()) - int(K.sum()) - int(n.sum()) + hits))
                         / ((int(n.sum()) - hits) * (int(K.sum()) - hits)))
        if (int(n.sum()) - hits) * (int(K.sum()) - hits) else np.nan,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--protocol-sha256', required=True)
    args = parser.parse_args()
    assert digest(args.protocol) == args.protocol_sha256, 'Unified protocol changed'
    input_paths = [args.protocol, PREVIOUS / 'measured_ortholog_members.tsv.gz',
                   PREVIOUS / 'ortholog_assay_features.tsv.gz', PREVIOUS / 'ortholog_3x3.tsv',
                   SOURCE / 'manifest.json', SOURCE / 'go-basic.obo',
                   SOURCE / 'human_go_ids.tsv', SOURCE / 'mouse_go_ids.tsv']
    hashes = {str(path): digest(path) for path in input_paths}
    source_manifest = json.loads((SOURCE / 'manifest.json').read_text())
    assert source_manifest['ensembl_release'] == 116
    for record in source_manifest['sources']:
        assert digest(SOURCE / record['file']) == record['sha256'], record['file']
    members = pd.read_csv(PREVIOUS / 'measured_ortholog_members.tsv.gz', sep='\t')
    features = pd.read_csv(PREVIOUS / 'ortholog_assay_features.tsv.gz', sep='\t', float_precision='round_trip')
    terms, aliases, obsolete, closure, ontology_version = ontology(SOURCE / 'go-basic.obo')
    mappings = {}
    for species in ['human', 'mouse']:
        wanted = set(members[species + '_gene_id'])
        mappings[species] = annotations(SOURCE / (species + '_go_ids.tsv'), wanted,
                                        terms, aliases, obsolete, closure)
    definitions = [dict(category=category, root_GO=term, term_name=terms[term]['name'],
                        namespace=terms[term]['namespace'])
                   for category, roots in CATEGORIES.items() for term in roots]
    save(pd.DataFrame(definitions), 'category_definitions.tsv')
    results, power, mde, cohorts, mapping_rows, evidence_rows, coverage, roster, strata = [], [], [], [], [], [], [], [], []
    for contrast, direction, expected_n in [('P_M', 'proximal', 27), ('P_S', 'distal', 14)]:
        pairs = members[members.contrast.eq(contrast) & members.threshold.eq(.02)].reset_index(drop=True)
        assert len(pairs) == {'P_M': 751, 'P_S': 1656}[contrast]
        assert pairs.human_gene_id.is_unique and pairs.mouse_gene_id.is_unique
        assert pairs.homology_type.eq('ortholog_one2one').all() and pairs.is_high_confidence.eq(1).all()
        foreground = (pairs.human_direction.eq(direction) & pairs.mouse_direction.eq(direction)).to_numpy()
        assert int(foreground.sum()) == expected_n
        strict = members[members.contrast.eq(contrast) & members.threshold.eq(.05)].reset_index(drop=True)
        assert list(zip(strict.human_gene_id, strict.mouse_gene_id)) == list(zip(pairs.human_gene_id, pairs.mouse_gene_id))
        strict_foreground = (strict.human_direction.eq(direction) & strict.mouse_direction.eq(direction)).to_numpy()
        assert int(strict_foreground.sum()) == {'P_M': 7, 'P_S': 4}[contrast]
        assert not (strict_foreground & ~foreground).any()
        context = dict(contrast=contrast, shared_direction=direction, threshold=.02)
        models = {'Fisher_unmatched': (np.zeros(len(pairs), int), np.repeat('all', len(pairs)))}
        feat = features[features.contrast.eq(contrast) & features.threshold.eq(.02)]
        for model, columns in MODELS.items():
            models[model] = joint_codes(pairs, feat, columns)
        for model, (codes, labels) in models.items():
            for code in np.unique(codes):
                ii = codes == code
                strata.append(context | dict(model=model, stratum=int(code), stratum_key=labels[ii][0],
                                             N=int(ii.sum()), n=int(foreground[ii].sum())))
        for i, pair in pairs.iterrows():
            base = dict(contrast=contrast, human_gene_id=pair.human_gene_id, mouse_gene_id=pair.mouse_gene_id,
                        human_gene_symbol=pair.human_gene_symbol, mouse_gene_symbol=pair.mouse_gene_symbol,
                        foreground_002=bool(foreground[i]), foreground_005=bool(strict_foreground[i]),
                        human_direction_002=pair.human_direction, mouse_direction_002=pair.mouse_direction,
                        human_direction_005=strict.iloc[i].human_direction,
                        mouse_direction_005=strict.iloc[i].mouse_direction)
            row = base.copy()
            for species in ['human', 'mouse']:
                direct, active, ignored = mappings[species]
                gene = pair[species + '_gene_id']
                row[species + '_any_active_GO'] = bool(active[gene])
                row[species + '_ignored_GO'] = ';'.join(sorted(ignored[gene]))
                for category, roots in CATEGORIES.items():
                    matches = [term for term in direct[gene] if set(roots).intersection(closure(term))]
                    row[species + '_' + category] = bool(matches)
                    mapping_rows.append(base | dict(annotation_species=species, category=category,
                                                   annotated=bool(matches), direct_GO=';'.join(sorted(matches))))
                    for term in sorted(matches):
                        evidence_rows.append(base | dict(annotation_species=species, category=category,
                                                         direct_GO=term, direct_GO_name=terms[term]['name'],
                                                         matched_root_GO=';'.join(sorted(set(roots).intersection(closure(term))))))
            roster.append(row)
        for species in ['human', 'mouse']:
            direct, active, ignored = mappings[species]
            ids = pairs[species + '_gene_id']
            annotation = np.array([bool(active[gene]) for gene in ids])
            for threshold, target in [(.02, foreground), (.05, strict_foreground)]:
                coverage.append(dict(contrast=contrast, annotation_species=species, threshold=threshold,
                                     N=len(pairs), n=int(target.sum()), any_active_GO=int(annotation.sum()),
                                     target_any_active_GO=int(annotation[target].sum()),
                                     unknown_GO=int((~annotation).sum()), target_unknown_GO=int((~annotation[target]).sum())))
            for category, roots in CATEGORIES.items():
                category_membership = np.array([bool(set(roots).intersection(active[gene])) for gene in ids])
                for threshold, target in [(.02, foreground), (.05, strict_foreground)]:
                    cohorts.append(dict(contrast=contrast, shared_direction=direction, threshold=threshold,
                                        annotation_species=species, category=category, N=len(pairs),
                                        n=int(target.sum()), K=int(category_membership.sum()),
                                        hits=int((category_membership & target).sum()),
                                        background_fraction=float(category_membership.mean()),
                                        target_fraction=float(category_membership[target].mean())))
                for model, (codes, labels) in models.items():
                    params, null, row = calculate(category_membership, foreground, codes)
                    ctx = context | dict(annotation_species=species, category=category, model=model)
                    results.append(ctx | row)
                    lo, probabilities = null
                    for j, probability in enumerate(probabilities):
                        mapping_rows.append(dict(contrast=contrast, annotation_species=species,
                                                 category=category, model=model, null_hits=lo + j,
                                                 null_probability=float(probability)))
                    for alpha_name, alpha in [('nominal', .05), ('Bonferroni_14', .05 / 14)]:
                        for odds in [1, 2, 3, 5]:
                            generated = distribution(params, odds)
                            glo, gp = generated
                            mean = float(np.arange(glo, glo + len(gp)) @ gp)
                            power.append(ctx | dict(alpha_name=alpha_name, alpha=alpha, conditional_odds=odds,
                                                    generated_expected_hits=mean,
                                                    rejection_probability=rejection_probability(null, generated, alpha)))
                        mde.append(ctx | dict(alpha_name=alpha_name, alpha=alpha, target_power=.8,
                                              maximum_searched_odds=1000000,
                                              **detectable_effect(params, alpha, maximum_odds=1000000)))
        print(contrast, 'complete', flush=True)
    result = pd.DataFrame(results)
    result['q_14_per_annotation_model'] = result.groupby(['annotation_species', 'model']).p_enrichment.transform(bh)
    result['q_56_joint_matched_family'] = np.nan
    matched = result.model.ne('Fisher_unmatched')
    assert int(matched.sum()) == 56
    result.loc[matched, 'q_56_joint_matched_family'] = bh(result.loc[matched, 'p_enrichment'])
    save(result, 'enrichment.tsv')
    save(pd.DataFrame(power), 'conditional_power.tsv')
    save(pd.DataFrame(mde), 'minimum_detectable_effect.tsv')
    save(pd.DataFrame(cohorts), 'threshold_descriptive_counts.tsv')
    combined = pd.DataFrame(mapping_rows)
    save(combined[combined.model.isna()].drop(columns=['model', 'null_hits', 'null_probability']), 'category_membership.tsv.gz')
    save(combined[combined.model.notna()][['contrast', 'annotation_species', 'category', 'model', 'null_hits', 'null_probability']],
         'exact_null_distributions.tsv.gz')
    save(pd.DataFrame(evidence_rows), 'annotation_evidence.tsv.gz')
    save(pd.DataFrame(coverage), 'annotation_coverage.tsv')
    save(pd.DataFrame(roster), 'ortholog_roster.tsv.gz')
    save(pd.DataFrame(roster)[lambda x: x.foreground_002], 'target_members.tsv')
    save(pd.DataFrame(strata), 'stratum_margins.tsv.gz')
    assert hashes == {str(path): digest(path) for path in input_paths}, 'Input changed during analysis'
    validation = dict(status='COMPLETE_AWAITING_INDEPENDENT_QA', protocol_sha256=args.protocol_sha256,
                      source_ontology_version=ontology_version, source_ensembl_release=116,
                      primary_family=14, matched_joint_family=56, tests_at_005=False,
                      input_sha256=hashes, output_sha256={p.name: digest(p) for p in HERE.iterdir() if p.suffix in {'.tsv', '.gz'}},
                      code_sha256={p.name: digest(p) for p in [HERE / 'analyze.py', HERE / 'exact_null.py', HERE / 'test_exact_null.py']})
    (HERE / 'analysis_validation.json').write_text(json.dumps(validation, indent=2) + '\n')


if __name__ == '__main__':
    main()
