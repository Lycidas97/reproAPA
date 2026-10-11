"""Independent ontology traversal and SciPy conditional-hypergeometric replay."""
from pathlib import Path
from functools import lru_cache
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
SOURCE = Path('/Users/lycidas/work/科研/projects/2026-09-26_reproapa-go-analysis/sources')
checks = 0


def check(condition, label):
    global checks
    assert condition, label
    checks += 1


def close(a, b, label):
    check(np.isclose(a, b, rtol=2e-11, atol=2e-13), label)


def main():
    definitions = pd.read_csv(HERE / 'category_definitions.tsv', sep='\t')
    roots = dict(zip(definitions.category, definitions.root_GO))
    graph, aliases, obsolete = {}, {}, set()
    for block in (SOURCE / 'go-basic.obo').read_text().split('[Term]\n')[1:]:
        lines = block.split('\n[')[0].splitlines()
        ids = [line[4:] for line in lines if line.startswith('id: GO:')]
        if len(ids) != 1:
            continue
        term = ids[0]
        if 'is_obsolete: true' in lines:
            obsolete.add(term)
            continue
        parents = set()
        for line in lines:
            if line.startswith('is_a: GO:'):
                parents.add(line.split()[1])
            elif line.startswith('relationship: part_of GO:'):
                parents.add(line.split()[2])
            elif line.startswith('alt_id: GO:'):
                aliases[line.split()[1]] = term
        graph[term] = parents

    @lru_cache(None)
    def ancestors(term):
        term = aliases.get(term, term)
        if term not in graph or term in obsolete:
            return frozenset()
        seen, todo = set(), [term]
        while todo:
            t = todo.pop()
            if t in seen:
                continue
            seen.add(t)
            todo.extend(graph.get(t, set()))
        return frozenset(seen)

    roster = pd.read_csv(PROJECT / 'notes/2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz', sep='\t')
    pairs = roster[roster.contrast.eq('P_M')].reset_index(drop=True)
    fg = pairs.foreground_002.to_numpy(bool)
    strict = pairs.foreground_005.to_numpy(bool)
    check(len(pairs) == 751 and fg.sum() == 27 and strict.sum() == 7, 'fixed sizes')
    raw = pd.read_csv(PROJECT / 'notes/2026-10-10-recurrence-utr-ortholog/measured_ortholog_members.tsv.gz', sep='\t')
    for threshold, expected in [(.02, fg), (.05, strict)]:
        part = raw[raw.contrast.eq('P_M') & raw.threshold.eq(threshold)].reset_index(drop=True)
        check(part.human_gene_id.to_list() == pairs.human_gene_id.to_list(), 'raw identity')
        check(np.array_equal(expected, (part.human_direction.eq('proximal') & part.mouse_direction.eq('proximal')).to_numpy()), 'raw foreground')
    membership = pd.read_csv(HERE / 'memberships.tsv.gz', sep='\t')
    result = pd.read_csv(HERE / 'enrichment.tsv', sep='\t')
    features = pd.read_csv(PROJECT / 'notes/2026-10-10-recurrence-utr-ortholog/ortholog_assay_features.tsv.gz', sep='\t')
    features = features[features.contrast.eq('P_M') & features.threshold.eq(.02)]
    coarse = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
    models = {'joint_measurement_PAS_capacity': coarse,
              'joint_full_assay': coarse + ['mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin'],
              'Fisher_unmatched': []}
    all_p = []
    source_terminal_markers = 0
    for species in ['human', 'mouse']:
        wanted = set(pairs[species + '_gene_id'])
        ann = {gene: set() for gene in wanted}
        with (SOURCE / (species + '_go_ids.tsv')).open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                if line.strip() == '[success]':
                    source_terminal_markers += 1
                    continue
                gene, term = line.rstrip('\n').split('\t', 1)
                gene = gene.split('.')[0]
                if gene in ann:
                    ann[gene].update(ancestors(term.strip()))
        for category, root in roots.items():
            hits = np.array([root in ann[g] for g in pairs[species + '_gene_id']])
            reported = membership[membership.annotation_species.eq(species) & membership.category.eq(category)]
            check(reported.human_gene_id.to_list() == pairs.human_gene_id.to_list(), 'membership order')
            check(np.array_equal(hits, reported.member.to_numpy(bool)), 'independent GO membership')
            for model, columns in models.items():
                labels = []
                for i, pair in pairs.iterrows():
                    values = []
                    for sp in ['human', 'mouse']:
                        f = features[features.species.eq(sp)].set_index('gene_id')
                        values.extend(f.loc[pair[sp + '_gene_id'], columns].astype(str).to_list())
                    labels.append(tuple(values))
                groups = {}
                for i, label in enumerate(labels):
                    groups.setdefault(label, []).append(i)
                dist = np.array([1.0])
                for indices in groups.values():
                    N, K, n = len(indices), int(hits[indices].sum()), int(fg[indices].sum())
                    support = np.arange(n + 1)
                    dist = np.convolve(dist, hypergeom.pmf(support, N, K, n))
                observed = int((hits & fg).sum())
                mean = float(np.arange(len(dist)) @ dist)
                p = float(dist[observed:].sum())
                sd = float(np.sqrt(((np.arange(len(dist)) - mean) ** 2) @ dist))
                row = result[result.annotation_species.eq(species) & result.category.eq(category) & result.model.eq(model)].iloc[0]
                check(row.hits == observed and row.K == hits.sum(), 'observed/K')
                close(row.expected, mean, 'conditional expectation')
                close(row.p_enrichment, p, 'conditional upper tail')
                close(row.null_sd, sd, 'conditional null variance')
                if model == 'Fisher_unmatched':
                    N, K, n = len(fg), int(hits.sum()), int(fg.sum())
                    fisher = fisher_exact([[observed, n - observed], [K - observed, N - K - n + observed]], alternative='greater').pvalue
                    close(p, fisher, 'Fisher independent')
                all_p.append((species, category, model, p))

    def bh(values):
        ordered = sorted(enumerate(values), key=lambda x: x[1])
        output = np.ones(len(values))
        bound = 1.0
        for rank in range(len(values), 0, -1):
            index, value = ordered[rank - 1]
            bound = min(bound, value * len(values) / rank)
            output[index] = bound
        return output

    for mask, name in [(result.model.ne('Fisher_unmatched'), 'q16_matched_family'),
                       (result.model.eq('Fisher_unmatched'), 'q8_unmatched_description')]:
        check(np.allclose(result.loc[mask, name], bh(result.loc[mask, 'p_enrichment']), atol=2e-13), name)
    for _, rows in result.groupby(['annotation_species', 'model']):
        check(np.allclose(rows.q4_per_annotation_model, bh(rows.p_enrichment), atol=2e-13), 'q4')
    close(hypergeom.sf(1, 4, 2, 2), 1 / 6, 'toy exact')
    check(hypergeom.pmf(2, 4, 4, 2) == 1, 'degenerate all-hit')
    close(fisher_exact([[2, 0], [0, 2]], alternative='greater').pvalue, 1 / 6, 'toy Fisher')
    freeze = json.loads((HERE / 'freeze.json').read_text())
    for path, sha in freeze['sha256'].items():
        check(hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha, 'source hash ' + path)
    validation = {'status': 'PASS', 'checks': checks, 'skipped': 0,
                  'non_annotation_BioMart_success_markers': source_terminal_markers,
                  'independent_GO_graph_and_annotation': True,
                  'independent_SciPy_null': True,
                  'validator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE / 'independent_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
