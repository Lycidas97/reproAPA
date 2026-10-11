"""Exploratory, fixed-category functional test of the frozen P-M core."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
OLD = PROJECT / 'notes/2026-10-10-rna-module-enrichment'
PREVIOUS = PROJECT / 'notes/2026-10-10-recurrence-utr-ortholog'
SOURCE = Path('/Users/lycidas/work/科研/projects/2026-09-26_reproapa-go-analysis/sources')
CATEGORIES = {
    'ER_membrane': 'GO:0005789',
    'ER_to_Golgi_transport': 'GO:0006888',
    'vesicle_transport': 'GO:0016192',
    'cytoskeleton_organization': 'GO:0007010',
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(rows, filename):
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(HERE / filename, sep='\t', index=False)


def main():
    freeze = json.loads((HERE / 'freeze.json').read_text())
    for path, sha in freeze['sha256'].items():
        assert digest(path) == sha, path
    sys.path.insert(0, str(OLD))
    spec = importlib.util.spec_from_file_location('frozen_rna_go', OLD / 'analyze.py')
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    roster = pd.read_csv(OLD / 'ortholog_roster.tsv.gz', sep='\t')
    pairs = roster[roster.contrast.eq('P_M')].reset_index(drop=True)
    assert len(pairs) == 751 and pairs.human_gene_id.is_unique and pairs.mouse_gene_id.is_unique
    target = pairs.foreground_002.to_numpy(bool)
    strict = pairs.foreground_005.to_numpy(bool)
    assert target.sum() == 27 and strict.sum() == 7 and not (strict & ~target).any()
    raw = pd.read_csv(PREVIOUS / 'measured_ortholog_members.tsv.gz', sep='\t')
    raw = raw[raw.contrast.eq('P_M') & raw.threshold.eq(.02)].reset_index(drop=True)
    assert list(zip(raw.human_gene_id, raw.mouse_gene_id)) == list(zip(pairs.human_gene_id, pairs.mouse_gene_id))
    assert np.array_equal(target, (raw.human_direction.eq('proximal') & raw.mouse_direction.eq('proximal')).to_numpy())
    features = pd.read_csv(PREVIOUS / 'ortholog_assay_features.tsv.gz', sep='\t')
    features = features[features.contrast.eq('P_M') & features.threshold.eq(.02)]
    models = {name: helpers.joint_codes(pairs, features, cols)
              for name, cols in helpers.MODELS.items()}
    models['Fisher_unmatched'] = (np.zeros(len(pairs), int), np.repeat('all', len(pairs)))
    terms, aliases, obsolete, closure, version = helpers.ontology(SOURCE / 'go-basic.obo')
    for term in CATEGORIES.values():
        assert term in terms and term not in obsolete
    manifest = json.loads((SOURCE / 'manifest.json').read_text())
    assert manifest['ensembl_release'] == 116
    for source in manifest['sources']:
        assert digest(SOURCE / source['file']) == source['sha256']
    save([{'category': c, 'root_GO': g, 'name': terms[g]['name'], 'namespace': terms[g]['namespace']}
          for c, g in CATEGORIES.items()], 'category_definitions.tsv')
    results, memberships, evidence, coverage, nulls, strata = [], [], [], [], [], []
    for model, (codes, labels) in models.items():
        for code in np.unique(codes):
            ii = codes == code
            strata.append({'model': model, 'stratum': int(code), 'key': labels[ii][0],
                           'N': int(ii.sum()), 'n': int(target[ii].sum())})
    for species in ['human', 'mouse']:
        ids = pairs[species + '_gene_id'].to_list()
        direct, active, ignored = helpers.annotations(SOURCE / (species + '_go_ids.tsv'),
                                                      set(ids), terms, aliases, obsolete, closure)
        annotated = np.array([bool(active[g]) for g in ids])
        coverage.append({'annotation_species': species, 'N': len(ids), 'n': int(target.sum()),
                         'any_active_GO': int(annotated.sum()),
                         'target_any_active_GO': int(annotated[target].sum())})
        for category, root in CATEGORIES.items():
            member = np.array([root in active[g] for g in ids])
            for i, gene in enumerate(ids):
                matched = sorted(g for g in direct[gene] if root in closure(g))
                memberships.append({'annotation_species': species, 'category': category,
                                    'human_gene_id': pairs.iloc[i].human_gene_id,
                                    'mouse_gene_id': pairs.iloc[i].mouse_gene_id,
                                    'human_symbol': pairs.iloc[i].human_gene_symbol,
                                    'mouse_symbol': pairs.iloc[i].mouse_gene_symbol,
                                    'foreground_002': bool(target[i]), 'foreground_005': bool(strict[i]),
                                    'member': bool(member[i]), 'direct_GO': ';'.join(matched)})
                for term in matched:
                    evidence.append({'species': species, 'gene_id': gene, 'category': category,
                                     'root_GO': root, 'direct_GO': term, 'name': terms[term]['name']})
            for model, (codes, _) in models.items():
                _, dist, result = helpers.calculate(member, target, codes)
                context = {'annotation_species': species, 'category': category, 'model': model}
                results.append(context | result)
                lo, probs = dist
                nulls.extend(context | {'hits': lo + j, 'probability': float(p)} for j, p in enumerate(probs))
    result = pd.DataFrame(results)
    result['q4_per_annotation_model'] = result.groupby(['annotation_species', 'model']).p_enrichment.transform(helpers.bh)
    result['q16_matched_family'] = np.nan
    matched = result.model.ne('Fisher_unmatched')
    assert matched.sum() == 16
    result.loc[matched, 'q16_matched_family'] = helpers.bh(result.loc[matched, 'p_enrichment'])
    result['q8_unmatched_description'] = np.nan
    result.loc[~matched, 'q8_unmatched_description'] = helpers.bh(result.loc[~matched, 'p_enrichment'])
    save(result, 'enrichment.tsv')
    save(memberships, 'memberships.tsv.gz')
    membership = pd.DataFrame(memberships)
    save(membership[membership.foreground_002], 'target_memberships.tsv')
    save(evidence, 'annotation_evidence.tsv.gz')
    save(coverage, 'annotation_coverage.tsv')
    save(nulls, 'exact_null_distributions.tsv.gz')
    save(strata, 'stratum_margins.tsv')
    counts = []
    for (species, category), part in membership.groupby(['annotation_species', 'category']):
        for threshold, flag in [(.02, 'foreground_002'), (.05, 'foreground_005')]:
            counts.append({'annotation_species': species, 'category': category, 'threshold': threshold,
                           'N': len(part), 'K': int(part.member.sum()), 'n': int(part[flag].sum()),
                           'hits': int((part.member & part[flag]).sum())})
    save(counts, 'threshold_descriptive_counts.tsv')
    for path, sha in freeze['sha256'].items():
        assert digest(path) == sha, path
    (HERE / 'analysis_validation.json').write_text(json.dumps({
        'status': 'AWAITING_INDEPENDENT_QA', 'ontology_version': version,
        'protocol_sha256': digest(HERE / 'PROTOCOL.zh.md'),
        'inputs_sha256': freeze['sha256'],
        'outputs_sha256': {p.name: digest(p) for p in sorted(HERE.glob('*.tsv*'))},
        'code_sha256': digest(__file__),
    }, indent=2) + '\n')
    print(result[['annotation_species', 'category', 'model', 'hits', 'expected',
                  'p_enrichment', 'q16_matched_family', 'inferentially_mobile']].to_string(index=False))


if __name__ == '__main__':
    main()
