"""Retain the upstream cohort; annotate legacy resource coverage without filtering."""
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NOTES = Path('/Users/lycidas/Documents/ChatGPT/reproAPA/notes')
CLAIMS = NOTES / '2026-09-14-integration-reliability/data/A05_claim_registry_v1/claim_index.tsv'
V4 = NOTES / '2026-09-13-figure5-site-change-v4/data'
REMOTE_CHECKED_SHA = '449b59569f9bfc0c4b6b187893883f481c6e484bbe40ff0efc3d5d77a3442ee3'


def read_tsv(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def build(claims, models, sites):
    indexed = defaultdict(list)
    keys = [(c['species'], c['gene_id'], c['contrast']) for c in claims]
    assert len(set(keys)) == len(keys), 'Duplicate upstream claims'
    for claim in claims:
        assert claim['species'] in {'human', 'mouse'}
        indexed[(claim['species'], claim['gene_id'])].append(claim)
    routes = defaultdict(set)
    by_gene = defaultdict(dict)
    for key, model in models.items():
        sp, route, gid = key.split('|')
        assert model['gene_id'] == gid and model['route'] == route
        routes[sp].add(route)
        by_gene[(sp, gid)][route] = model
    counts = Counter()
    for site in sites:
        key = (site['species'], site['route'], site['gene_id'])
        assert '|'.join(key) in models, 'Site without model'
        counts[(*key, site['category'])] += 1

    pools, coverage, extensions = [], [], []
    for (sp, gid), gene_claims in sorted(indexed.items()):
        symbols = {c['gene_symbol'] for c in gene_claims}
        assert len(symbols) == 1
        available = by_gene[(sp, gid)]
        pools.append(dict(species=sp, gene_id=gid, gene_symbol=next(iter(symbols)),
                          cohort='A05_fixed_candidate_pool',
                          contrasts=';'.join(sorted(c['contrast'] for c in gene_claims)),
                          legacy_model_routes=';'.join(sorted(available)),
                          legacy_model_count=len(available),
                          resource_status='legacy_model_available_not_revalidated' if available else 'needs_model_lookup_or_build',
                          current_site_effect_status='not_recomputed'))
        for c in sorted(gene_claims, key=lambda x: x['contrast']):
            contrast = c['contrast']
            assert len(contrast) == 2, 'Explicit state mapping required for this contrast'
            assert c['delta_definition'] == f'{contrast[1]} minus {contrast[0]}'
            for route in sorted(routes[sp]):
                model = available.get(route)
                if model is None:
                    status = 'needs_model_lookup_or_build'
                    states, missing = '', ';'.join(contrast)
                else:
                    states = ';'.join(model['states'])
                    missing = ';'.join(s for s in contrast if s not in model['states'])
                    status = 'needs_state_usage' if missing else 'legacy_states_present_not_revalidated'
                row = dict(species=sp, gene_id=gid, gene_symbol=c['gene_symbol'],
                           claim_id=c['claim_id'], contrast=contrast, delta_definition=c['delta_definition'],
                           route=route, legacy_states=states, missing_states=missing,
                           coverage_status=status, current_site_effect_status='not_recomputed')
                for cat in ['miRNA', 'RBP', 'RBP_transfer']:
                    n = counts[(sp, route, gid, cat)] if model is not None else None
                    row[f'legacy_{cat}_records'] = '' if n is None else n
                    row[f'{cat}_resource_status'] = ('unknown_no_model' if n is None else
                                                    'legacy_records_present' if n else 'no_legacy_category_records')
                coverage.append(row)
    for (sp, gid), available in sorted(by_gene.items()):
        if (sp, gid) not in indexed:
            extensions.append(dict(species=sp, gene_id=gid,
                                   gene_symbol=next(iter(available.values()))['gene_symbol'],
                                   legacy_model_routes=';'.join(sorted(available)),
                                   cohort='legacy_extension_outside_A05', in_primary_denominator=False))
    return pools, coverage, extensions


def main():
    assert hashlib.sha256(CLAIMS.read_bytes()).hexdigest() == REMOTE_CHECKED_SHA, 'Upstream changed: review version first'
    claims = read_tsv(CLAIMS)
    models = json.loads((V4 / 'genes.json').read_text())
    sites = read_tsv(V4 / 'site_registry.tsv.gz')
    pools, coverage, extensions = build(claims, models, sites)
    out = ROOT / 'data'
    out.mkdir(exist_ok=True)
    summary = {}
    for sp in ['mouse', 'human']:
        selected = [r for r in pools if r['species'] == sp]
        original = [r for r in claims if r['species'] == sp]
        assert {r['gene_id'] for r in selected} == {r['gene_id'] for r in original}
        assert len(selected) == int(original[0]['fixed_candidate_denominator'])
        write_tsv(out / f'{sp}_candidates.tsv', selected)
        write_tsv(out / f'{sp}_gene_contrasts.tsv', original)
        summary[sp] = dict(candidates=len(selected), gene_contrasts=len(original),
                           with_legacy_model=sum(r['legacy_model_count'] > 0 for r in selected),
                           without_legacy_model=sum(r['legacy_model_count'] == 0 for r in selected),
                           legacy_extension_genes=sum(r['species'] == sp for r in extensions))
    write_tsv(out / 'study_contrast_coverage.tsv', coverage)
    write_tsv(out / 'model_state_backlog.tsv', [r for r in coverage if r['coverage_status'].startswith('needs_')], list(coverage[0]))
    write_tsv(out / 'legacy_extension_candidates.tsv', extensions)
    manifest = dict(cohort_version='candidate-pool-v1', upstream='A05_claim_registry_v1',
                    same_site_requirement='none_for_admission; within_species_only_for_replication_annotation',
                    effect_threshold_for_admission=None, baseline_support_required=False,
                    species_files={sp: f'data/{sp}_candidates.tsv' for sp in summary},
                    contrast_files={sp: f'data/{sp}_gene_contrasts.tsv' for sp in summary},
                    reviewed_routes=sorted({r['route'] for r in coverage}),
                    upstream_remote_sha256_checked=REMOTE_CHECKED_SHA,
                    sources={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in [CLAIMS, V4 / 'genes.json', V4 / 'site_registry.tsv.gz']},
                    summary=summary, coverage_counts=dict(Counter(r['coverage_status'] for r in coverage)),
                    site_effects_recomputed=False)
    (ROOT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest['summary'], indent=2))


if __name__ == '__main__':
    main()
