"""Acceptance checks for source identity, strict boundaries and evidence scope."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows(relative):
    p = HERE / relative
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    checks = {}
    def require(name, value):
        checks[name] = bool(value)
        assert value, name
    for r in json.loads((HERE / 'file-manifest.json').read_text()):
        p = HERE / r['path']
        require('file:' + r['path'], p.stat().st_size == r['bytes'] and
                hashlib.sha256(p.read_bytes()).hexdigest() == r['sha256'])
    overlay = json.loads((HERE / 'components.json').read_text())
    require('polyApipe primary only', overlay['source']['run_id'] == 'polyapipe_common_recount_20261005_v2_primary'
            and overlay['source']['reference_variant'] == 'v2_primary')
    require('all S26 components', {r['panel'] for r in overlay['components']} == {'S26'+x for x in 'ABCDEFGH'})
    for component in overlay['components']:
        for r in component['tables']:
            require(component['id'] + ':' + r['path'], hashlib.sha256((HERE/r['path']).read_bytes()).hexdigest() == r['sha256'])
    mapping = rows('tables/ago2/query_mapping.tsv')
    require('current Trip12 minus-strand queries only', len(mapping) == 70 and
            all(r['gene_id'] == 'ENSMUSG00000026219' and r['strand'] == '-' and
                r['mapping_status'] == 'unique_gapfree_reciprocal' for r in mapping))
    identities = rows('tables/ago2/sequence_identity.tsv')
    require('all query IDs have sequence evidence', {r['feature_id'] for r in identities} == {r['feature_id'] for r in mapping})
    require('assembly sequence identity', all(r['sequence_identical'] == 'True' for r in identities))
    overlaps = rows('tables/ago2/overlap_summary.tsv')
    require('every query has all five tracks', len(overlaps) == 350 and
            len({(r['feature_id'], r['label']) for r in overlaps}) == 350)
    require('technical WT libraries do not create biological support',
            all(r['WT_biological_units_supported'] == '0' for r in rows('tables/ago2/biological_unit_support.tsv')))
    boundary = {int(r['length_nt']): r for r in rows('tables/longread/boundary_pair_audit.tsv')}
    require('706 nt strict; 842 nt boundary is not silently admitted',
            boundary[706]['strict_full_exclusion'] == 'True' and boundary[706]['proximal_overlap_nt'] == '0' and
            boundary[842]['strict_full_exclusion'] == 'False' and boundary[842]['proximal_overlap_nt'] == '2' and
            boundary[842]['published_strict_pairs'] == '0')
    require('structural compatibility does not become joint detection',
            all(r['both_TPM_positive_samples'] == '0' for r in boundary.values()))
    selected = rows('tables/candidate_screen/example_selection.tsv')
    require('original example rank is preserved', len(selected) == 1 and selected[0]['rank'] == '2')
    geometry = rows('tables/geometry/axis-decomposition.tsv')
    require('both native contributor sets retained',
            {r['route']: int(r['n_native_pairs']) for r in geometry} == {'CIP': 3, 'WOT': 5})
    require('exact native decomposition', all(abs(sum(float(r[k]) for k in
            ['within_family_DWUI', 'family_mass_DWUI', 'nonterminal_DWUI'])-float(r['native_delta_DWUI'])) < 1e-12 for r in geometry))
    require('two strict native low-shift observations', all(abs(float(r['native_delta_DWUI'])) < .01
            and float(r['native_MPRO']) >= .1 for r in geometry))
    pairs = rows('tables/examples/inclusion_paired_libraries.tsv.gz')
    for route in ['CIP', 'WOT']:
        for cutoff in ['20', '50']:
            for a, b in [('84700478', '84701184'), ('84701184', '84702026')]:
                z = [r for r in pairs if r['route'] == route and r['cutoff'] == cutoff and r['kind'] == 'interval'
                     and r['start0'] == a and r['end0'] == b]
                n = 3 if route == 'CIP' else 24 if cutoff == '20' else 20
                require(f'paired interval contributors {route}/{cutoff}/{a}', len(z) == n and len({r['library_key'] for r in z}) == n)
                require(f'interval strict effect {route}/{cutoff}/{a}', sum(float(r['delta']) for r in z)/n <= -.1)
    v = json.loads((HERE / 'verification.json').read_text())
    require('raw long-read model and abundance verification completed', v['longread_raw_models_and_TPM_verified'])
    out = dict(checks=len(checks), passed=sum(checks.values()), failures=[k for k,v in checks.items() if not v])
    (HERE/'acceptance.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out))


if __name__ == '__main__':
    main()
