"""Independent table/coordinate checks; no imports of production analysis functions."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', low_memory=False)


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    for manifest in ['event_sequence_sealed_manifest.json', 'mirna_sealed_manifest.json']:
        seal = json.loads((HERE / manifest).read_text())
        for name, expected in seal['files'].items():
            assert digest(HERE / name) == expected, (manifest, name)
    inputs = {}
    for name in ['site_input_hashes.json', 'interval_input_hashes.json', 'chain_input_hashes.json']:
        for path, expected in json.loads((HERE / name).read_text()).items():
            assert digest(ROOT / path) == expected, path
            assert path not in inputs or inputs[path] == expected
            inputs[path] = expected
    source = ROOT / 'notes/2026-10-10-marker-reference-sensitivity/inputs.npz'
    assert digest(source) == json.loads((HERE / 'naive_actual_axis_validation.json').read_text())['input_npz_sha256']
    inputs[str(source.relative_to(ROOT))] = digest(source)
    config_path = ROOT / 'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'
    config = json.loads(config_path.read_text())
    for sp in ['human', 'mouse']:
        path = ROOT / f'notes/2026-10-05-apadetect-analysis-intake/completion-inputs/external/{sp}_mirna_annotations.tsv'
        assert digest(path) == config['resources'][sp + '_mirna_annotations']['sha256']
    contexts = read(HERE / 'site_event_contexts.tsv.gz')
    unique = ['source_population', 'gene_threshold', 'tolerance', 'human_gene_id', 'contrast']
    assert len(contexts) == 6888 and not contexts.duplicated(unique).any()
    selected = contexts[contexts.sequence_candidate]
    assert (selected.parent_case & selected.context_role.eq('shared_recurrent') |
            ~selected.parent_case & selected.context_role.eq('single_species_recurrent')).all()
    assert selected.strict_UTR.all() and selected.mass_eligible.all()
    assert selected[selected.parent_case].strict_double_driver.all()
    assert contexts[contexts.two_study_exact_event_both].strict_double_driver.all()
    assert contexts[contexts.two_study_exact_event_both].human_gain_support_studies.ge(2).all()
    assert contexts[contexts.two_study_exact_event_both].mouse_gain_support_studies.ge(2).all()
    original_query = read(HERE / 'before-method-repair/fixed_window_queries.tsv')
    current_query = read(HERE / 'fixed_window_queries.tsv')
    metadata = read(HERE / 'fixed_window_metadata.tsv')
    assert len(original_query) == len(metadata) == 116 and len(current_query) == 112
    assert set(current_query.id) <= set(original_query.id) == set(metadata.id)
    window = current_query.merge(metadata, on='id', validate='one_to_one')
    assert np.max(np.abs(window.phastcons - window.phastcons_replayed)) == 0
    assert window.finite_fraction.ge(.8).all() and window.ACGT_bases.eq(150).all()
    tests = read(HERE / 'sequence_tests.tsv')
    features = read(HERE / 'sequence_gene_features_and_strata.tsv')
    assert len(tests) == 8 and int((~tests.estimable).sum()) == 6
    with np.load(HERE / 'sequence_null_draws.npz', allow_pickle=False) as archive:
        for test in tests.itertuples():
            frame = features[features.set_name.eq(test.set_name) & features.tolerance.eq(test.tolerance)]
            column = 'sequence_' + test.model
            assert int(frame.case.sum()) == test.cases
            observed = frame.loc[frame.case, 'score'].mean()
            assert abs(observed - test.observed_mean) < 1e-14
            expected, variance = 0., 0.
            for _, block in frame.groupby(column, sort=True):
                score = block.score.to_numpy()
                n, k = len(block), int(block.case.sum())
                expected += k * score.mean()
                if n > 1:
                    variance += k * (n-k) / (n * (n-1)) * np.sum((score-score.mean()) ** 2)
            expected /= test.cases
            variance /= test.cases ** 2
            assert abs(expected - test.conditional_expected) < 1e-14
            assert abs(np.sqrt(variance) - test.null_sd) < 1e-14
            null = archive[f'{test.set_name}|{test.tolerance}|{test.model}']
            assert len(null) == 20000
            tail = int((null >= observed - 1e-12).sum())
            assert tail == test.upper_tail_count
            assert abs(test.p - ((tail+1)/20001 if test.estimable else 1)) < 1e-14
    p = tests.p.to_numpy(); order = np.argsort(p)
    q = np.empty(8); q[order] = np.minimum(1, np.minimum.accumulate((p[order] * 8 / np.arange(1,9))[::-1])[::-1])
    assert np.max(np.abs(q - tests.q8)) == 0
    intervals = read(HERE / 'strict_case_intervals.tsv')
    effects = read(HERE / 'strict_case_interval_inclusion.tsv')
    seeds = read(HERE / 'strict_case_predicted_seeds.tsv')
    query = read(HERE / 'seed_mapping_queries.tsv')
    pairs = read(HERE / 'reciprocal_predicted_seed_pairs.tsv')
    utr = read(BASE / 'tandem/annotation_terminal_utrs.tsv.gz')
    physical = {sp: read(BASE / f'external_mirna/{sp}_physical_seed_sites.tsv.gz') for sp in ['human', 'mouse']}
    keys = ['source_population', 'tolerance', 'contrast', 'human_gene_id', 'species']
    assert len(intervals) == 22 and len(seeds) == 1163 and len(query) == 70
    assert not intervals.duplicated(keys).any() and query.query_id.is_unique
    for row in intervals.itertuples():
        proximal, distal = json.loads(row.proximal_native_cuts), json.loads(row.distal_native_cuts)
        expected_core = (max(proximal), min(distal)) if row.strand == '+' else (max(distal), min(proximal))
        assert expected_core == (row.start0, row.end0) and row.start0 < row.end0
        tx = utr[utr.species.eq(row.species) & utr.gene_id.eq(row.gene_id) & utr.transcript_id.isin(row.common_transcripts.split(';'))]
        assert len(tx) == len(row.common_transcripts.split(';'))
        assert tx.utr_start.sub(1).le(row.start0).all() and tx.utr_end.ge(row.end0).all()
        current = effects[(effects[keys] == pd.Series({k:getattr(row,k) for k in keys})).all(axis=1)]
        means = current.groupby('study_id')[['target_inclusion','reference_inclusion','delta_inclusion']].mean().mean()
        for name in means.index:
            assert abs(means[name]-getattr(row,name)) < 1e-14
    assert np.max(np.abs(effects.target_inclusion-effects.reference_inclusion-effects.delta_inclusion)) < 1e-14
    for row in seeds.itertuples():
        assert row.start0 >= row.core_start0 and row.end0 <= row.core_end0 and row.end0-row.start0 in [7,8]
        original = physical[row.species]
        assert len(original[original.gene_id.eq(row.gene_id) & original.seed.eq(row.seed) & original.chrom.eq(row.chrom) &
                            original.strand.eq(row.strand) & original.start0.eq(row.start0) & original.end0.eq(row.end0)]) == 1
    maps = {r['query_id']:r for r in map(json.loads,(HERE/'seed_chain_maps.jsonl').read_text().splitlines())}
    qlookup = query.set_index('query_id')
    for pair in pairs.itertuples():
        h,m = maps[pair.human_query_id], maps[pair.mouse_query_id]
        assert h['status'] == m['status'] == 'unique'
        for source,target in [(h,qlookup.loc[pair.mouse_query_id]),(m,qlookup.loc[pair.human_query_id])]:
            hit=source['hits'][0]
            assert (hit['chrom'],hit['strand'],hit['start0'],hit['end0']) == (target.chrom,target.strand,target.start0,target.end0)
        assert qlookup.loc[pair.human_query_id,'seed'] == qlookup.loc[pair.mouse_query_id,'seed'] == pair.seed
    assert len(pairs) == 46 and pairs.drop_duplicates(['human_query_id','mouse_query_id']).shape[0] == 23
    report = dict(status='PASS', contexts=6888, fixed_current_windows=112, eight_sequence_tests_exact=True,
                  sequence_NE=6, sequence_estimable=2, interval_species_contexts=22, predicted_seed_context_rows=1163,
                  small_seed_queries=70, homologous_predicted_pairs_including_tolerance_duplicates=46,
                  unique_homologous_predicted_pairs=23, all_sealed_files_and_input_hashes_stable=True)
    (HERE/'final_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    (HERE/'final_input_hashes.json').write_text(json.dumps(inputs,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__ == '__main__':
    main()
