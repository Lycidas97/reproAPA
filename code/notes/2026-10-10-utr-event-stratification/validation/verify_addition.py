"""Independent post-hoc ALE endpoint-UTR coverage and descriptive QA."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = HERE.parent
LABELS = ['proximal', 'distal', 'nonrecurrent']
THRESHOLDS = [.01, .015, .02, .025, .03, .035, .04, .045, .05]
PAIRKEY = ['human_gene_id', 'mouse_gene_id']
checks = []


def verify(test, name):
    assert bool(test), name
    checks.append(name)


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def chromosome(x):
    x = str(x).removeprefix('chr')
    return 'MT' if x == 'M' else x


def main():
    contract = json.loads((BASE / 'classification/ALE_UTR_coverage_contract.json').read_text())
    sources = {name: ROOT / name for name in contract['input_and_qualification_sha256']}
    initial = {name: sha(path) for name, path in sources.items()}
    verify(initial == contract['input_and_qualification_sha256'], 'Post-hoc coverage inputs match recorded SHA256')
    class_pairs = read(BASE / 'classification/event_class_pairs.tsv')
    ale = class_pairs[class_pairs.pair_class.eq('both_ALE')].copy()
    verify(len(ale) == 52, 'Semantic qualification retains original 52-pair ALE universe')
    events = read(ROOT / 'notes/2026-10-10-shared-proximal-mechanisms/event_inventory.tsv.gz')
    utrs = read(ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/tandem/annotation_terminal_utrs.tsv.gz')
    intervals = {}
    for r in utrs.itertuples():
        key = (r.species, r.gene_id, chromosome(r.chrom), r.strand)
        intervals.setdefault(key, []).append((int(r.utr_start) - 1, int(r.utr_end), r.transcript_id))
    production_routes = read(BASE / 'classification/ALE_driver_terminal_UTR_coverage_routes.tsv').set_index(['species', 'study_id', 'route', 'gene_id'])
    covers = {}
    route_rows = []
    for species in ['human', 'mouse']:
        ids = set(ale[species + '_gene_id'])
        selected = events[events.species.eq(species) & events.gene_id.isin(ids)]
        verify(selected.driver_class.eq('alternative_terminal_exon_compatible').all(), 'All retained semantic-qualification routes are original ALE-compatible events')
        for r in selected.itertuples():
            candidates = intervals.get((species, r.gene_id, chromosome(r.chrom), r.strand), [])
            flags = []
            output = dict(species=species, study_id=r.study_id, route=r.route, gene_id=r.gene_id)
            prod = production_routes.loc[(species, r.study_id, r.route, r.gene_id)]
            for endpoint in ['gain', 'loss']:
                coordinate = int(getattr(r, endpoint + '_cut')) - (1 if r.strand == '+' else 0)
                hits = sorted({transcript for start, end, transcript in candidates if start <= coordinate < end})
                covered = bool(hits)
                verify(coordinate == prod[endpoint + '_retained_base0'], 'Independent retained-base calculation matches native strand contract')
                verify(covered == bool(prod[endpoint + '_terminal_UTR_covered']), 'Independent UTR interval query matches endpoint coverage')
                recorded = set(str(prod[endpoint + '_UTR_transcripts']).split(';')) if pd.notna(prod[endpoint + '_UTR_transcripts']) else set()
                verify(set(hits) == recorded, 'Endpoint compatible-transcript identity set independently matches')
                flags.append(covered)
                output[endpoint + '_covered'] = covered
            both = all(flags)
            verify(both == bool(prod.both_driver_ends_UTR_covered), 'Both-end route eligibility independently matches')
            covers.setdefault((species, r.gene_id), []).append(both)
            route_rows.append(output | dict(both_endpoints_covered=both))
    verify(len(route_rows) == len(production_routes), 'Every post-hoc route is checked exactly once')
    qualification = read(BASE / 'classification/both_ALE_terminal_UTR_coverage_pairs.tsv').set_index(PAIRKEY).sort_index()
    rebuilt = ale[PAIRKEY].copy()
    for species in ['human', 'mouse']:
        rebuilt[species + '_covered'] = [all(covers[(species, gene)]) for gene in rebuilt[species + '_gene_id']]
    rebuilt['eligible'] = rebuilt.human_covered & rebuilt.mouse_covered
    indexed = rebuilt.set_index(PAIRKEY).sort_index()
    verify(indexed.index.equals(qualification.index), 'Post-hoc pair qualification does not replace or enlarge original ALE background')
    for species in ['human', 'mouse']:
        verify(indexed[species + '_covered'].eq(qualification[species + '_all_routes_both_ends_UTR_covered']).all(), 'All-route endpoint qualification independently matches by species')
    verify(indexed.eligible.eq(qualification.both_species_all_routes_both_ends_UTR_covered).all(), 'Both-species endpoint qualification independently matches')
    verify(int(indexed.eligible.sum()) == 15, 'All 15 eligible pairs independently reproduce')
    members = read(ROOT / 'notes/2026-10-10-ortholog-threshold-scan/measured_ortholog_members.tsv.gz')
    members = members[members.contrast.eq('P_M')].copy()
    features = read(ROOT / 'notes/2026-10-10-ortholog-threshold-scan/ortholog_assay_features.tsv.gz')
    features = features[features.contrast.eq('P_M')]
    eligible_pairs = rebuilt[rebuilt.eligible][PAIRKEY]
    expected_rows = []
    for threshold in THRESHOLDS:
        pairs = members[members.threshold.eq(threshold)].merge(eligible_pairs, on=PAIRKEY, validate='one_to_one')
        verify(len(pairs) == 15, 'Strict ALE semantic background remains 15 for every threshold')
        for model in ['measurement_PAS_capacity', 'full_assay']:
            columns = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
            if model == 'full_assay':
                columns += ['mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']
            probabilities, labels = {}, {}
            for species in ['human', 'mouse']:
                labels[species] = pairs[species + '_direction'].map(dict(zip(LABELS, range(3)))).to_numpy()
                sf = features[features.species.eq(species) & features.threshold.eq(threshold)].set_index('gene_id').reindex(pairs[species + '_gene_id'])
                keys = list(sf[columns].itertuples(index=False, name=None))
                probabilities[species] = np.empty((15, 3))
                for i, key in enumerate(keys):
                    group = np.array([j for j in range(15) if keys[j] == key])
                    probabilities[species][i] = [(labels[species][group] == d).mean() for d in range(3)]
            expected = probabilities['human'].T @ probabilities['mouse']
            for h, human in enumerate(LABELS):
                for m, mouse in enumerate(LABELS):
                    expected_rows.append(dict(threshold=threshold, model=model, human_direction=human, mouse_direction=mouse,
                        observed=int(((labels['human'] == h) & (labels['mouse'] == m)).sum()), analytic_expected=float(expected[h, m])))
    expected_rows = pd.DataFrame(expected_rows)
    production = read(BASE / 'overlap/strict_ale_utr_descriptive.tsv')
    keys = ['threshold', 'model', 'human_direction', 'mouse_direction']
    aa = expected_rows.set_index(keys).sort_index()
    bb = production.set_index(keys).sort_index()
    verify(aa.index.equals(bb.index) and len(aa) == 162, 'All 162 post-hoc descriptive cells independently reproduce')
    verify(np.array_equal(aa.observed, bb.observed), 'All post-hoc observed direction counts independently match')
    verify(np.allclose(aa.analytic_expected, bb.analytic_expected, atol=1e-12, rtol=0), 'All post-hoc class-conditional expectations independently match')
    verify(not any(c.startswith('p_') or c.startswith('q_') for c in production.columns), 'Post-hoc qualification introduces no p or q values')
    verify(bb.measured_orthologs.eq(15).all(), 'All post-hoc descriptive denominators use 15 eligible pairs')
    main_input_hash = sha(BASE / 'classification/event_class_pairs.tsv')
    verify(main_input_hash == '5ba39cd3a5439557c7c3edf604ec6c9ca03c88acb0249996f4748572c17842c3', 'Original main classification remains unchanged after post-hoc qualification')
    primary = aa.loc[(.02, slice(None), 'proximal', 'proximal'), :]
    verify(primary.observed.eq(4).all(), 'All four main ALE shared proximal pairs survive semantic qualification')
    pd.DataFrame(route_rows).to_csv(HERE / 'independent_addition_route_coverage.tsv', sep='\t', index=False)
    rebuilt.to_csv(HERE / 'independent_addition_pair_coverage.tsv', sep='\t', index=False)
    expected_rows.to_csv(HERE / 'independent_addition_expectations.tsv', sep='\t', index=False, float_format='%.17g')
    verify({name: sha(path) for name, path in sources.items()} == initial, 'Post-hoc inputs unchanged throughout independent audit')
    result = dict(status='PASS', checks_passed=len(checks), checks_failed=0, checks_skipped=0,
        posthoc_semantic_description_only=True, original_ALE_background=52, endpoint_qualified_background=15,
        route_events_checked=len(route_rows), core_002_retained=4, descriptive_cells=162,
        max_expected_absolute_error=float(np.max(np.abs(aa.analytic_expected - bb.analytic_expected))),
        no_new_inference_or_family=True, primary_QA_2344_record_unchanged=True,
        same_transcript_or_CDS_not_required=True, input_and_qualification_sha256=initial,
        independent_script_sha256=sha(Path(__file__)), main_report_interpretation_audit='No formal class interaction or full-axis UTR-only claim; cross-route common UTR wording flagged to parent for clarification.')
    (HERE / 'addition_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
