#!/usr/bin/env python3
"""Verify human input coverage and reject fabricated time or replicate structure."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    config = yaml.safe_load(Path(parser.parse_args().config).read_text())
    root = Path(config['output_root'])
    prep = root / 'prepared'
    read = lambda name: pd.read_csv(prep / name, sep='\t', dtype=str)
    primary, sensitivity = read('primary_cells.tsv.gz'), read('sensitivity_cells.tsv.gz')
    reference, libraries = read('reference_cells.tsv.gz'), read('libraries.tsv')
    technical, edges = read('technical_exclude_cells.tsv.gz'), read('temporal_edge_plan.tsv')
    checks = {}

    def check(name, value):
        checks[name] = bool(value)
        if not value:
            raise ValueError(name)

    check('all_qc_tiers_disjoint', pd.concat([primary.cell_id, sensitivity.cell_id, technical.cell_id]).is_unique)
    check('primary_count_and_unique_ids', len(primary) == 215114 and primary.cell_id.is_unique)
    check('sensitivity_count_and_unique_ids', len(sensitivity) == 42878 and sensitivity.cell_id.is_unique)
    check('technical_count', len(technical) == 79941)
    check('all_59_libraries_retained', libraries.library_key.is_unique and len(libraries) == 59 and primary.library_key.nunique() == 59)
    check('reference_primary_only', reference.cell_id.is_unique and set(reference.cell_id) <= set(primary.cell_id))
    check('every_library_in_reference', set(reference.library_key) == set(libraries.library_key))
    check('all_62_original_states_retained', primary.original_state_id.nunique() == 62)
    check('all_original_states_in_reference', set(reference.original_state_id) == set(primary.original_state_id))
    state_coverage = primary.groupby('original_state_id').size().rename('primary_cells').to_frame()
    state_coverage['reference_cells'] = reference.groupby('original_state_id').size()
    state_coverage.to_csv(prep / 'reference_state_coverage.tsv', sep='\t')
    for study, spec in config['studies'].items():
        scope = Path(config['clustering_root']) / spec['clustering_release']
        original = pd.read_csv(scope / 'frozen/cell_clusters.tsv.gz', sep='\t', dtype=str, usecols=['cell_id', 'cluster_id'])
        joined = primary[primary.study_id.eq(study)][['cell_id', 'cluster_id']].merge(original, on='cell_id', how='outer', suffixes=('_new', '_old'), validate='one_to_one')
        check(study + '_original_cluster_identity_unchanged', joined.cluster_id_new.eq(joined.cluster_id_old).all())
    pooled = libraries[libraries.study_id.eq('GSE150311')]
    check('pooled_and_culture_samples_have_no_fabricated_day', pooled.actual_day.isna().all())
    check('pooled_study_has_no_temporal_edges', not edges.route_id.eq('GSE150311').any())
    chemical = libraries[libraries.study_id.eq('GSE178325')]
    check('chemical_stage_day_not_global_day', chemical.actual_day.isna().all())
    check('ordinal_edges_have_no_elapsed_day', edges.loc[edges.relation.ne('actual_day'), 'elapsed_days'].isna().all())
    exact = edges[edges.relation.eq('actual_day')]
    check('numerical_edges_move_forward', pd.to_numeric(exact.elapsed_days).gt(0).all())
    new = libraries[libraries.study_id.eq('GSE221739')]
    check('new_study_has_no_invented_d0', pd.to_numeric(new.actual_day).min() == 3)
    source_counts = primary.groupby('reference_stratum').size()
    reference_counts = reference.groupby('reference_stratum').size()
    check('reference_strata_obey_fixed_cap', reference_counts.equals(source_counts.clip(upper=config['reference_cells_per_stratum'])))
    available = primary.groupby(['reference_stratum', 'library_key']).size()
    selected = reference.groupby(['reference_stratum', 'library_key']).size().reindex(available.index, fill_value=0)
    for stratum in source_counts.index:
        quota = selected.loc[stratum]
        remaining = quota[quota.lt(available.loc[stratum])]
        check('balanced_library_quota_' + stratum, remaining.empty or quota.max() - remaining.min() <= 1)
    lookup = libraries.set_index('library_key')
    reviews = []
    for row in edges.itertuples(index=False):
        keys = json.loads(row.source_libraries) + json.loads(row.target_libraries)
        check('same_study_' + row.route_id + '_' + row.time_from, lookup.loc[keys, 'study_id'].nunique() == 1)
        selected_endpoint = lookup.loc[keys, 'condition'].str.startswith('TRA-1-60_').any()
        reviews.append({'route_id': row.route_id, 'time_from': row.time_from, 'time_to': row.time_to,
                        'time_semantics': row.relation,
                        'mass_interpretation': 'conditional_on_sorted_capture_not_population_fractions' if selected_endpoint else 'conditional_sample_distribution',
                        'cells_and_split_captures_are_not_independent_biological_replicates': True,
                        'selection_strata_must_be_retained': bool(selected_endpoint), 'status': 'DESIGN_READY_NOT_SOLVED'})
    pd.DataFrame(reviews).to_csv(prep / 'temporal_design_review.tsv', sep='\t', index=False)
    validation = json.loads((root / 'reference/validation.json').read_text())
    check('raw_reference_reopened_and_validated', validation['status'] == 'PASS' and validation['per_cell_counts_conserved'])
    report = {'status': 'PASS', 'phase': 'INTEGRATION_PREPARATION_COMPLETE', 'checks': checks,
              'reference_cells': len(reference), 'primary_cells': len(primary), 'sensitivity_cells': len(sensitivity),
              'reference_matrix_ready': True, 'joint_model_trained': False, 'human_ot_solved': False, 'joint_layout_ready': False}
    (root / 'preparation_acceptance.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
