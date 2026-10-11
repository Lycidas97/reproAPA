#!/usr/bin/env python3
"""Supplementary-table plan for the polyApipe v2-primary manuscript (2026-10-06), consolidated to 8 workbooks.

Reads only: main/ text, contract r2 bundle manifest, supplementary-figure receipts (S1-S30), the S31 source pin,
main-figure render receipts and table/S1-S7 TSVs. Writes table-plan.tsv, table-sources.tsv, checks.json.
No statistics are recomputed.
"""
import csv, gzip, hashlib, json, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
NOTES = ROOT / 'notes'
BUNDLE = NOTES / '2026-10-06-supplement-polyapipe/contract-r2/bundle'
CAND = NOTES / '2026-10-06-supplement-polyapipe/candidates'
ALT = NOTES / '2026-10-06-alt-catalogue-suppfig'
MAINREC = {
    'F1': NOTES / '2026-10-06-main-figures-polyapipe/candidates/run3/F1/render-receipt.json',
    'F2': NOTES / '2026-10-06-main-figures-polyapipe/candidates/run3/F2/render-receipt.json',
    'F3': NOTES / '2026-10-06-main-figures-polyapipe/candidates/F3/render-receipt.json',
    'F4': NOTES / '2026-10-05-figure4-polyapipe-relayout/candidate-20261006-0536/render-receipt.json',
    'F5': NOTES / '2026-10-06-main-figures-polyapipe/candidates/run7/F5/render-receipt.json',
}
ORDER = ['D01', 'D02', 'D03', 'G1', 'G2', 'G3'] + [f'D{i:02d}' for i in range(4, 28)]
SNUM = {w: f'S{i+1}' for i, w in enumerate(ORDER)}
SNUM['A01'] = 'S31'

UP = [
    ('S1', '-', 'Reprocessed single-cell datasets', 'reuse', 'table/S1-datasets.tsv'),
    ('S2', '-', 'Sequencing libraries and Cell Ranger settings', 'reuse', 'table/S2-sequencing-libraries.tsv'),
    ('S3', '-', 'Expression-atlas libraries (161 retained of 210 aligned)', 'reuse', 'table/S3-expression-atlas-libraries.tsv'),
    ('S4', '-', 'Local expression partitions', 'reuse', 'table/S4-local-partitions.tsv'),
    ('S5', '-', 'Local-state annotations', 'reuse', 'table/S5-local-state-annotations.tsv'),
    ('S6', '-', 'Orthologs used for expression integration', 'reuse', 'table/S6-integration-orthologs.tsv'),
    ('S7', '-', 'Mouse annotation programmes', 'reuse', 'table/S7-mouse-annotation-programs.tsv'),
]

# number, evidence id, title, status, citing paragraphs, supp figures, receipt panel prefixes, core files, note
P = [
    ('S8', 'E01a', 'Expression states, marker evidence, partition stability and sampling inventory', 'revise',
     'R1P1', 'S1-S6', ['D01', 'D02', 'D03A', 'D03D', 'D03E', 'D03F', 'G1', 'G2', 'G3'],
     ['frozen/states.tsv', 'frozen/state_libraries.tsv', 'state_grouping/state_partition_summary.tsv',
      'state_grouping/state_library_inventory.tsv', 'global/state_time_inventory.tsv'],
     'Adds the G1-G3 marker/partition-stability tables (new Figures S4-S6); frozen GEX, no PAS values.'),
    ('S9', 'E02a', 'PAS catalogue, reference support and count QC', 'revise',
     'R1P2; Methods PAS identification; Methods robustness', 'S7', ['D04'],
     ['pas_qc/catalog_summary.tsv', 'pas_qc/support_and_annotation_strata.tsv', 'pas_qc/filter_reason_audit.tsv',
      'pas_qc/caller_filter_semantics.tsv'],
     'Carries 64.0-91.0% of primary PASs / 93.4-97.2% of UMIs within 50 nt. Caller-migration and call-availability audits cited by Methods are not in the current bundle.'),
    ('S10', 'E03a', 'Route means and distal excess across 16 measurement settings', 'update',
     'R1P3', 'S8a,b', ['D05A', 'D05B'],
     ['global/global_direction_summary.tsv', 'global/direction_stability.tsv', 'global/source_library_inventory.tsv'],
     'Primary values, 16-setting ranges and distal excess of Figure 1d.'),
    ('S11', 'E03b', 'Tandem-UTR restriction of route means', 'update',
     'R1P3', 'S8c', ['D05C'], ['tandem/tandem_utr_summary.tsv'],
     'Carries +0.027/+0.005/+0.015 vs -0.025/-0.054 (primary setting).'),
    ('S12', 'E01b', 'Documented-source design and the 1,306-ortholog common background', 'update',
     'R1P4', 'S3c', [], ['global/background_gene_membership.tsv.gz'],
     'Species x induction-class design; ortholog background behind the Figure 1d ticks.'),
    ('S13', 'E03c', 'Observed-state references and same-library reference triangle', 'update',
     'R1P5', 'S8d-f', ['D05D', 'D05F'],
     ['geometry/state-reference-summary.tsv', 'naive_state/condition_and_same_library_summary.tsv'],
     'R/E/P/M minus observed S per route (Figure 1d upper).'),
    ('S14', 'E04a', 'Profile-class resolution and candidate membership', 'update',
     'R2P1', 'S9a', ['D06A'],
     ['profiles/gene_membership_k3.tsv.gz', 'geometry/k3-assignment-margins.tsv', 'profiles/shape_centers.tsv',
      'profiles/resolution_comparison.tsv'],
     'Full lists of the 902 human and 796 mouse candidates with class (C1 progressive proximal, C2 reversal, C3 distal).'),
    ('S15', 'E04b', 'Route hold-out: shape recovery and assignment stability', 'update',
     'R2P2', 'S9d,e', ['D06C', 'D06D', 'D06E'], ['supplement_statistics/LOO-gene-assignments.tsv.gz'], ''),
    ('S16', 'E05a', 'Ortholog class correspondence (139 pairs) and composition nulls', 'update',
     'R2P3', 'S10a,b', ['D07A', 'D07B', 'D07C'], ['shape_ortholog/ortholog_shape_members.tsv.gz'], ''),
    ('S17', 'E17a+E17c', 'Functional tests of the profile classes: class ORA and ranked APA/GEX tests', 'merge',
     'R2P3', 'S17', ['D14'], [],
     'Merges v15 S17 (class ORA) and v15 S42 (route-specific APA-amplitude ranked tests); both are panels of the new Figure S17.'),
    ('S18', 'E05c', 'Adjacent-state contributions to path length', 'update',
     'R2P4', 'S10c,e', ['D07F'], ['profiles/architecture_route_summary.tsv'], ''),
    ('S19', 'E05b', 'Observed path geometry (D, L, T) and its descriptive limits', 'update',
     'R2P5', 'S9c; S10c,d', ['D06B', 'D07E'], ['profiles/trajectory_geometry.tsv.gz'],
     'Observed D/L and T only; no noise calibration.'),
    ('S20', 'E06a', 'Expression-graph organization (Moran) and APA-GEX profile agreement', 'update',
     'R2P6', 'S11', ['D08'], [], ''),
    ('S21', 'E07a', 'Within-species recurrence design: studies, series, references and estimators', 'update',
     'R3P1; Methods data collection, quantification, statistics', 'S12; S3c', ['D03C'],
     ['global/source_library_inventory.tsv'], 'Cohort/estimator definition table cited by Methods.'),
    ('S22', 'E07b', 'Recurrent genes, conditional nulls and three-study support', 'update',
     'R3P2', 'S12a,e,f', ['D09A', 'D09E', 'D09F'],
     ['recurrence/recurrence_summary.tsv', 'recurrence/recurrence_members.tsv.gz'],
     'Complete member lists: human/mouse P-S 483/604, P-M 268/420; 37 human P-M genes in all three studies.'),
    ('S23', 'E07c+D13', 'Direction composition of recurrent genes and intermediate-state alignment', 'revise',
     'R3P3', 'S12b; S16', ['D09B', 'D13'], ['intermediate/process/endpoint_eligibility.tsv'],
     'Adds the new intermediate-state analysis (41/293, 134/530, 233/527; reversals 15/29/5).'),
    ('S24', 'E08a', 'Measurement controls, culture-well inference and external PAS-seq', 'update',
     'R3P4; Methods robustness', 'S12c,d; S13', ['D09C', 'D09D', 'D10'],
     ['measurement/nulls/recurrence_summary.tsv', 'measurement/nulls/thinning_stability.tsv',
      'well_inference/gene_metric_tests.tsv.gz', 'well_inference/coverage_and_exclusions.tsv',
      'matched_pms/well_significance_overlap.tsv'],
     'Well family 5 x 614 (370 q<.05; 134/143 recurrent).'),
    ('S25', 'E09a', 'Common-contributor P/M/S comparisons and dual-reference candidates', 'update',
     'R3P5', 'S14; S21f,g', ['D11'],
     ['candidate_screen/dual_reference_candidates.tsv', 'candidate_screen/dual_reference_ladder.tsv',
      'matched_pms/matched_state_vectors.tsv.gz'],
     'Mouse 420 -> 395 -> 222; human 268 -> 250 -> 128; Pdlim7.'),
    ('S26', 'E10a', 'APA recurrence across opposing expression responses', 'update',
     'R3P6', 'S15', ['D12'], ['apa_gex/current_candidate_shared_state_relations.tsv.gz'],
     'Eight-test family; identity and direction now supported in both species.'),
    ('S27', 'E17b', 'Nested functional backgrounds (G0/G1/G2)', 'update',
     'R3P7', 'S18', ['D15'], ['functional/nested_universe_membership.tsv.gz'], ''),
    ('S28', 'E11b', 'Cross-species sharing, continuous correspondence and correction families', 'update',
     'R4P1; Methods orthologue mapping', 'S19b-f', ['D16B', 'D16C', 'D16D', 'D16E', 'D16F'],
     ['cross_sharing/sharing_summary.tsv', 'cross_sharing/sharing_members.tsv.gz',
      'cross_continuous/continuous_gene_block_summary.tsv'],
     'Renumbered from v15 S29 so that tables follow first citation (R4P1 before R4P2).'),
    ('S29', 'E11a', 'Shared proximal candidates and their conditional nulls', 'update',
     'R4P2', 'S19a; S20', ['D16A', 'D17'], ['cross_overlap/overlap_null_results.tsv'],
     'Renumbered from v15 S28. 23 P-S and 20 P-M candidates; three nulls.'),
    ('S30', 'E12a', 'Same-PAS cross-species correspondence', 'update',
     'R4P3', 'S21a,b,e', ['D18A', 'D18B', 'D18E'], ['same_pas/v2_primary_conditional_null_tests.tsv'],
     'v2 primary only (v3 removed).'),
    ('S31', 'E12b', 'Paired reference differences and the E-M control', 'update',
     'R4P4', 'S21c,d', ['D18C'], [], ''),
    ('S32', 'E12c', 'PDLIM7/Pdlim7 exon, CDS and PAS structure', 'update',
     'R4P5', 'S21f,g; S14e', ['D18H'],
     ['supplement_statistics/same-PAS-example-sites.tsv.gz', 'supplement_statistics/same-PAS-example-coordinate-check.tsv.gz'],
     ''),
    ('S33', 'E13a', 'Cis features and predictive-increment tests', 'update',
     'R4P6', 'S22; S23', ['D19', 'D20'], [],
     'Four primary tests in both prediction directions (min q = 0.069) + 32-setting family.'),
    ('S34', 'E14a', 'DWUI-MPRO union pool, interval models and inclusion', 'update',
     'R5P1; Methods data collection, quantification, statistics', 'S24a', ['D21A'],
     ['dwui_pool/candidates.tsv', 'union_pool/union_summary.tsv'],
     'Cohort/estimator definition table cited by Methods.'),
    ('S35', 'E14b', 'Selected-cohort interval agreement', 'update',
     'R5P2', 'S24c,d', ['D21C', 'D21D'], ['panel_statistics/observed-interval-agreement.tsv'], 'Descriptive only.'),
    ('S36', 'E14c', 'Weighted predicted-site retention', 'update',
     'R5P3', 'S24e,f', ['D21E'],
     ['weighted_retention/retention_summary.tsv', 'weighted_retention/interval_retention.tsv.gz'],
     'Loss/gain counts incl. the 107 exact-zero mouse RBP records.'),
    ('S37', 'E18a+D24', 'Matched sequence backgrounds, site density and external miRNA seed groups', 'revise',
     'R5P3; Methods site retention', 'S26; S27', ['D23', 'D24'], [],
     'Absorbs the external miRNA seed-group analysis (Figure S27, source of the Figure 5f PSC-miRNA row).'),
    ('S38', 'E18b', 'Candidate-RBP expression context and the Trip12 RBP filter', 'revise',
     'R5P3', 'S28', ['D25'], ['rbp_expression/expression_contrasts.tsv'],
     'Adds the 27-RBP filter behind Figure 5f (25 kept; Celf3, Rbm24 removed).'),
    ('S39', 'E15a', 'Exact-site partial retention and examples', 'update',
     'R5P4', 'S25', ['D22'],
     ['candidate_screen/partial_retention_relations.tsv.gz', 'site_statistics/cross_study_partial_retention.tsv'],
     '976/704 mouse and 235/231 human relationships; Eif4ebp2 (Figure 5e), Snx30 and SH3PXD2A examples.'),
    ('S40', 'E14d', 'All 51 low-net-shift mouse P-M genes: screen outcomes and axis decomposition', 'revise',
     'R5P5', 'S24b; S30', ['D21B', 'D27'],
     ['candidate_screen/low_shift_candidates.tsv', 'geometry/screen-outcomes.tsv', 'geometry/interval-eligibility.tsv'],
     'v15 title "All 64 ..." is obsolete: 51 -> 8 -> 4.'),
    ('S41', 'E16a', 'Trip12: interval inclusion, native axis, long-read models, AGO2 archive and site tiers', 'revise',
     'R5P6; Methods site retention', 'S29', ['D26'], ['trip12_evidence/ago2/overlap_summary.tsv'],
     'Replaces the v15 Rad23b content.'),
    ('S42', 'ALT', 'Main analyses repeated with alternative PAS calling', 'new',
     'Methods alternative PAS calling; Discussion', 'S31', ['A01'], [],
     'New; replaces the uncited v15 S42 (merged into S17). Gene-level percentages removed from the Discussion belong here.'),
]
# files that a panel straddling two tables reads, or that several figures read: explicit home (None = display config, not data)
OVERRIDE = {
    'frozen/palette.json': None, 'frozen/marker_selection.json': 'S8',
    'state_vectors/route_effects.primary.tsv.gz': 'S10',
    'profiles/all_resolution_memberships.tsv.gz': 'S14',
    'shape_ortholog/shape_contingency_3x3.tsv': 'S16',
    'profiles/architecture_route_summary.tsv': 'S18', 'profiles/trajectory_geometry.tsv.gz': 'S19',
    'site_statistics/region_consistency_statistics.tsv': 'S35',
    'candidate_screen/low_shift_attrition.tsv': 'S40', 'geometry/low-shift-cohort.tsv': 'S40',
    'geometry/axis-decomposition.tsv': 'S40',
    'candidate_screen/recurrent_region_pairs.tsv.gz': 'S35', 'site_statistics/site_factor_inventory.tsv.gz': 'S36',
}
MAIN = {
    'frozen/frozen_cell_embedding.tsv.gz': 'S8', 'frozen/S2A_state_geometry.tsv': 'S8',
    'axis_extensions/route_direction_audit.tsv': 'S10', 'figure_scores/metacell_scores.tsv.gz': 'S13',
    'geometry/state-reference-summary.tsv': 'S13', 'global/global_direction_summary.tsv': 'S10',
    'moran/k3_group_wpp_metacells.tsv.gz': 'S20', 'moran/moran_results.tsv': 'S20',
    'profiles/shape_centers.tsv': 'S14', 'profiles/leave_route_centroids.tsv': 'S15',
    'profiles/leave_route_projection.tsv': 'S15', 'profiles/route_shape_support.tsv': 'S15',
    'shape_ortholog/shape_tests.tsv': 'S16',
    'recurrence/recurrence_summary.tsv': 'S22', 'recurrence/study_event_covariates.tsv.gz': 'S22',
    'measurement/nulls/recurrence_summary.tsv': 'S24', 'recurrence_controls/recurrence_ladder.tsv': 'S22',
    'measurement/nulls/member_retention.tsv': 'S24', 'measurement/same-state/recurrence_summary.tsv': 'S24',
    'candidate_screen/dual_reference_ladder.tsv': 'S25', 'candidate_screen/dual_reference_candidates.tsv': 'S25',
    'matched_pms/matched_state_vectors.tsv.gz': 'S25', 'matched_pms/weighted_triplets.tsv.gz': 'S25',
    'candidate_screen/example_selection.tsv': 'S25',
    'apa_gex/prespecified_eight_test_family.tsv': 'S26', 'apa_gex/four_class_nulls.tsv.gz': 'S26',
    'cross_sharing/sharing_summary.tsv': 'S28', 'cross_continuous/continuous_gene_block_summary.tsv': 'S28',
    'cross_continuous/continuous_pair_summary.tsv': 'S28', 'cross_overlap/overlap_null_results.tsv': 'S29',
    'same_pas/v2_primary_correlation_summary.tsv': 'S30', 'same_pas/v2_primary_conditional_null_tests.tsv': 'S31',
    'same_pas/capacity_sensitivity_tests.tsv': 'S31', 'geometry/PAS-transcript-coordinate-check.tsv.gz': 'S32',
    'axis_extensions/homolog_site_pairs.tsv.gz': 'S30', 'cross_trace/member_endpoint_effects.tsv.gz': 'S32',
    'weighted_retention/retention_summary.tsv': 'S36', 'weighted_retention/gene_retention.tsv': 'S36',
    'candidate_screen/low_shift_candidates.tsv': 'S40', 'site_statistics/partial_retention_summary.tsv': 'S39',
    'site_statistics/cross_study_partial_retention.tsv': 'S39', 'site_statistics/exact_site_cross_study.tsv.gz': 'S39',
    'examples/inclusion_features.tsv': 'S39', 'examples/inclusion_library_states.tsv.gz': 'S39',
    'examples/inclusion_state_summary.tsv': 'S39',
    'examples/example_endpoint_effects.tsv.gz': 'S41', 'examples/inclusion_paired_libraries.tsv.gz': 'S41',
}
MAIN_PS = {
    'profile-': 'S14', 'transient-survival': 'S19', 'pooled-class-geometry': 'S19', 'signed-stage-means': 'S18',
    'stage-contribution': 'S18', 'recurrence-direction-counts': 'S23', 'well-denominators': 'S24',
    'example-cumulative-PAS-curves': 'S25', 'same-pas-rank-bin': 'S31', 'observed-interval-agreement': 'S35',
    'inclusion-paired-means': 'S39',
}


# ------------------------------------------------------------------ consolidation to 8 workbooks (user, 2026-10-06 07:06)
NEW = [
    ('S1', 'Datasets, libraries and expression states'),
    ('S2', 'PAS measurement and route-level PAS position'),
    ('S3', 'State-resolved APA profile classes'),
    ('S4', 'Within-species recurrence and its context'),
    ('S5', 'Cross-species comparisons and cis-feature models'),
    ('S6', 'Interval inclusion and predicted-site retention'),
    ('S7', 'Low-net-shift genes and the Trip12 example'),
    ('S8', 'Main analyses repeated with alternative PAS calling'),
]
# old evidence unit -> (new table, sheet letter, sheet title). Sheet order follows first citation in the text.
OLD2NEW = {
    'S8': ('S1', 'a', 'Expression states, marker evidence, partition stability and sampling'),
    'S1': ('S1', 'b', 'Datasets'), 'S2': ('S1', 'c', 'Sequencing libraries'),
    'S3': ('S1', 'd', 'Expression-atlas libraries'), 'S4': ('S1', 'e', 'Local expression partitions'),
    'S5': ('S1', 'f', 'Local-state annotations'), 'S6': ('S1', 'g', 'Integration orthologs'),
    'S7': ('S1', 'h', 'Mouse annotation programmes'),
    'S9': ('S2', 'a', 'PAS catalogue and reference support'),
    'S10': ('S2', 'b', 'Route means across 16 settings'), 'S11': ('S2', 'c', 'Tandem-UTR restriction'),
    'S12': ('S2', 'd', 'Design and 1,306-ortholog background'), 'S13': ('S2', 'e', 'Observed-state references'),
    'S14': ('S3', 'a', 'Resolution and candidate membership'), 'S15': ('S3', 'b', 'Route hold-out'),
    'S16': ('S3', 'c', 'Ortholog class correspondence'), 'S17': ('S3', 'd', 'Functional tests of the classes'),
    'S18': ('S3', 'e', 'Stage contributions and path geometry'), 'S19': ('S3', 'e', 'Stage contributions and path geometry'),
    'S20': ('S3', 'f', 'Moran structure and APA-GEX profile agreement'),
    'S21': ('S4', 'a', 'Design, series and estimators'), 'S22': ('S4', 'b', 'Recurrent genes and conditional nulls'),
    'S23': ('S4', 'c', 'Direction composition and intermediate states'),
    'S24': ('S4', 'd', 'Measurement controls, wells and external PAS-seq'),
    'S25': ('S4', 'e', 'Dual-reference candidates'), 'S26': ('S4', 'f', 'APA-GEX response identity and direction'),
    'S27': ('S4', 'g', 'Nested functional backgrounds'),
    'S28': ('S5', 'a', 'Sharing tests and correction families'), 'S29': ('S5', 'b', 'Shared proximal candidates'),
    'S30': ('S5', 'c', 'Same-PAS correspondence and reference differences'),
    'S31': ('S5', 'c', 'Same-PAS correspondence and reference differences'),
    'S32': ('S5', 'd', 'PDLIM7/Pdlim7 structure'), 'S33': ('S5', 'e', 'Cis features and increment tests'),
    'S34': ('S6', 'a', 'Union pool and interval models'), 'S35': ('S6', 'b', 'Interval agreement'),
    'S36': ('S6', 'c', 'Weighted predicted-site retention'),
    'S37': ('S6', 'd', 'Site density and external miRNA seed groups'), 'S38': ('S6', 'e', 'Candidate-RBP expression'),
    'S39': ('S6', 'f', 'Exact-site partial retention'),
    'S40': ('S7', 'a', 'The 51 low-net-shift genes'), 'S41': ('S7', 'b', 'Trip12'),
    'S42': ('S8', 'a', 'Alternative PAS calling'),
}
# small summary / test-family files promoted from figure source to sheet content (keyed by old unit)
EXTRA = {
    'S9': ['panel_statistics/primary-PAS-and-UMI-support.tsv'],
    'S15': ['profiles/leave_route_projection.tsv', 'profiles/leave_route_centroids.tsv'],
    'S16': ['shape_ortholog/shape_tests.tsv', 'shape_ortholog/shape_contingency_3x3.tsv'],
    'S17': ['profile_gex/enrichment_all_tests.tsv.gz', 'profile_gex/rank_enrichment_all_tests.tsv.gz'],
    'S20': ['moran/moran_results.tsv', 'profile_gex/APA_GEX_correlation_summary.tsv'],
    'S22': ['recurrence_controls/recurrence_ladder.tsv', 'recurrence/leave_one_study_out_nulls.tsv',
            'recurrence_controls/guan_series_unit_sensitivity.tsv'],
    'S23': ['panel_statistics/recurrence-direction-counts.tsv', 'intermediate/process/coverage_flow.tsv',
            'intermediate/process/fixed_cohort_outcomes.tsv'],
    'S24': ['measurement/nulls/member_retention.tsv', 'measurement/same-state/recurrence_summary.tsv',
            'external_pas/vector_summary.tsv'],
    'S26': ['apa_gex/prespecified_eight_test_family.tsv'],
    'S27': ['functional/functional_families.tsv', 'supplement_statistics/endpoint-GO-KEGG-display.tsv'],
    'S28': ['direction_difference/direction_contrast_tests.tsv', 'cross_continuous/continuous_pair_summary.tsv'],
    'S29': ['cross_overlap/new_shared_members.tsv.gz', 'cross_overlap/overlap_threshold_layer_sensitivity.tsv'],
    'S30': ['same_pas/v2_primary_correlation_summary.tsv', 'same_pas/v2_primary_same_axis_gene_effects.tsv.gz'],
    'S31': ['same_pas/capacity_sensitivity_tests.tsv'],
    'S33': ['cis/recomputed_models.tsv', 'cis/recomputed_coefficients.tsv', 'cis/recomputed_cv.tsv'],
    'S34': ['intervals/model_effects.tsv'],
    'S35': ['site_statistics/region_consistency_statistics.tsv', 'site_statistics/all_changed_region_pairs.tsv.gz'],
    'S37': ['motif_background/density_summary.tsv', 'motif_background/gene_density.tsv',
            'external_mirna/primary_tests.tsv', 'external_mirna/complete_test_family.tsv',
            'external_mirna/external_seed_expression.tsv'],
    'S38': ['rbp_expression/candidates.tsv', 'rbp_expression/state_expression.tsv'],
    'S39': ['site_statistics/partial_retention_summary.tsv'],
    'S40': ['candidate_screen/low_shift_attrition.tsv', 'geometry/axis-decomposition.tsv'],
    'S41': ['examples/inclusion_effect_sensitivity.tsv', 'trip12_evidence/ago2/biological_unit_support.tsv',
            'trip12_evidence/longread/boundary_pair_audit.tsv', 'trip12_evidence/figure5/site-display-evidence.tsv',
            'trip12_evidence/figure5/site-display-summary.tsv'],
}
# row filters to apply when a sheet is written (long-format sources)
FILTERS = {
    'global/background_gene_membership.tsv.gz': 'background = common_one_to_one_ortholog_all_settings (the 1,306-ortholog background)',
    'recurrence/recurrence_members.tsv.gz': 'one series_matching value (True/False rows are identical in count)',
    'apa_gex/current_candidate_shared_state_relations.tsv.gz': 'tested class only (same-direction APA, opposite GEX)',
    'same_pas/v2_primary_same_axis_gene_effects.tsv.gz': 'tolerance 25 nt',
}
EXCEL_MAX = 1048575
# the current text still uses v15 numbers, in which sharing is S29 and the proximal sets S28; plan units swap them
V15_TO_PLAN = {28: 29, 29: 28}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def nrows(p):
    op = gzip.open if str(p).endswith('.gz') else open
    try:
        with op(p, 'rt', errors='replace') as f:
            return sum(1 for _ in f) - 1
    except Exception:
        return ''


def main():
    man = {t['path']: t for t in json.load(open(BUNDLE / 'table-manifest.json'))['tables']}
    plan = {r[0]: r for r in P}
    src, issues = {}, []

    def add(tab, path, role, via, pin=None):
        key = (tab, path)
        if key in src:
            src[key]['via'] = ';'.join(sorted(set(src[key]['via'].split(';') + [via])))
            if src[key]['role'] == 'core_list' and role != 'core_list':
                src[key]['role'] = 'core_list'
            return
        src[key] = dict(role=role, via=via, pin=pin)

    # 1. supplementary-figure receipts
    unassigned = []
    for d in sorted(CAND.glob('*')):
        if d.name.startswith('_') or not d.is_dir():
            continue
        rec = json.load(open(next(d.glob('*.receipt.json'))))
        for t, v in rec['tables'].items():
            rel = t.replace('tables/', '', 1)
            via = SNUM[d.name] + ':' + ','.join(v['panels'])
            if rel in OVERRIDE:
                if OVERRIDE[rel]:
                    add(OVERRIDE[rel], rel, 'figure_source', via, v['sha256'])
                continue
            homes = sorted({n for n, r in plan.items() for pnl in v['panels'] if any(pnl.startswith(x) for x in r[6])})
            if len(homes) != 1:
                unassigned.append((d.name, rel, v['panels'], homes))
                continue
            add(homes[0], rel, 'figure_source', via, v['sha256'])
    # 2. S31 (alternative calling) pinned sources
    for k, v in json.load(open(ALT / 'sources.json'))['files'].items():
        add('S42', v['path'], 'figure_source_external', 'S31:' + k, v['sha256'])
    # 3. core lists (listed last so they keep role core_list)
    for n, r in plan.items():
        for c in list(r[7]) + EXTRA.get(n, []):
            if (n, c) in src:
                src[(n, c)]['role'] = 'core_list'
            else:
                add(n, c, 'core_list', 'plan')
    # 4. main-figure reads
    main_unassigned = []
    for fig, rp in MAINREC.items():
        for t, s in json.load(open(rp))['tables'].items():
            rel = t.replace('tables/', '', 1)
            tab = MAIN.get(rel)
            if tab is None and rel.startswith('panel_statistics/'):
                tab = next((v for k, v in MAIN_PS.items() if rel.split('/', 1)[1].startswith(k)), None)
            if tab is None and rel in OVERRIDE:
                tab = OVERRIDE[rel]
            if tab is None:
                tab = next((t2 for (t2, p2) in src if p2 == rel), None)
            if tab is None:
                main_unassigned.append((fig, rel))
                continue
            add(tab, rel, 'main_figure_source', fig, s if isinstance(s, str) else None)

    # 5. resolve paths, SHA, rows
    rows = []
    for (tab, path), r in src.items():
        if r['role'] == 'figure_source_external':
            p, mrec = ROOT / path, None
        else:
            p, mrec = BUNDLE / 'tables' / path, man.get('tables/' + path)
        ok = p.exists()
        cur = sha(p) if ok else ''
        exp = r['pin'] or (mrec or {}).get('sha256', '')
        if not ok:
            issues.append(f'missing {tab} {path}')
        elif exp and cur != exp:
            issues.append(f'sha-mismatch {tab} {path}')
        elif r['role'] != 'figure_source_external' and mrec is None:
            issues.append(f'not-in-manifest {tab} {path}')
        size = p.stat().st_size if ok else ''
        count = ''
        if ok and not path.endswith('.json') and (r['role'] == 'core_list' or size < 50e6):
            count = nrows(p)
        cols = (mrec or {}).get('columns') or []
        if ok and not cols and not path.endswith('.json'):
            op = gzip.open if path.endswith('.gz') else open
            with op(p, 'rt') as fh:
                cols = fh.readline().rstrip('\n').split('\t')
        rows.append(dict(table=tab, role=r['role'], via=r['via'], path=path, bytes=size, rows=count, sha256=cur,
                         columns=','.join(cols)))
    for n, e, t, s, f in UP:
        p = ROOT / f
        rows.append(dict(table=n, role='upstream', via='table/', path=f, bytes=p.stat().st_size, rows=nrows(p),
                         sha256=sha(p), columns=open(p).readline().rstrip('\n').replace('\t', ',')))
    # 6. consolidation: old unit -> new table / sheet; sheet content vs source-data archive
    for r in rows:
        nt, sh, title = OLD2NEW[r['table']]
        r['old_unit'] = r.pop('table')
        r['table'], r['sheet'], r['sheet_title'] = nt, nt + sh, title
        r['layer'] = 'sheet' if r['role'] in ('upstream', 'core_list', 'figure_source_external') else 'source_data'
        r['filter'] = FILTERS.get(r['path'], '') if r['layer'] == 'sheet' else ''
        if r['layer'] == 'sheet' and r['rows'] != '' and int(r['rows']) > EXCEL_MAX and not r['filter']:
            issues.append(f"sheet exceeds Excel rows {r['sheet']} {r['path']} {r['rows']}")
    rank = {'upstream': 0, 'core_list': 1, 'figure_source_external': 1, 'figure_source': 2, 'main_figure_source': 3}
    rows.sort(key=lambda r: (int(r['table'][1:]), r['sheet'], r['layer'] != 'sheet', rank[r['role']], r['path']))
    fields = ['table', 'sheet', 'sheet_title', 'layer', 'old_unit', 'role', 'via', 'path', 'bytes', 'rows', 'filter',
              'sha256', 'columns']
    rows = [{k: r[k] for k in fields} for r in rows]

    # 7. citations in the current main text and their replacements
    cites = []
    pat = re.compile(r'Tables?~((?:S\d+)(?:(?:--|,\s*|\s+and\s+)S?\d+)*)')
    for f in sorted((ROOT / 'main').rglob('*.tex')):
        para = ''
        for line in open(f):
            m = re.match(r'%\s*(R\dP\d+)', line)
            if m:
                para = m.group(1)
            for m in pat.finditer(line):
                g = m.group(1)
                if '--' in g:
                    a, b = [int(x) for x in re.findall(r'\d+', g)][:2]
                    nums = list(range(a, b + 1))
                else:
                    nums = [int(x) for x in re.findall(r'\d+', g)]
                nums = [V15_TO_PLAN.get(n, n) for n in nums]
                sheets = []
                for n in nums:
                    sh = OLD2NEW[f'S{n}'][0] + OLD2NEW[f'S{n}'][1]
                    if sh not in sheets:
                        sheets.append(sh)
                sec = f.relative_to(ROOT / 'main').as_posix()
                if sec.startswith('backbone/'):
                    new_text = 'rewrite by hand (S1-S8 overview)'
                elif len(sheets) == 1:
                    new_text = 'Table~' + sheets[0]
                else:
                    new_text = 'Tables~' + ', '.join(sheets[:-1]) + ' and ' + sheets[-1]
                cites.append(dict(file=sec, para=para, old=m.group(0), new=new_text, old_units=nums))
    drv = [l.split('{')[1].split('}')[0] for l in open(ROOT / 'main.tex') if l.startswith('\\input{')]
    order = []
    for d in drv:
        for c in cites:
            if c['file'] == d.replace('main/', '', 1) + '.tex' and not c['file'].startswith('backbone/'):
                for n in c['old_units']:
                    sh = ''.join(OLD2NEW[f'S{n}'][:2])
                    if sh not in order:
                        order.append(sh)
    tab_order = []
    for sh in order:
        if sh[:-1] not in tab_order:
            tab_order.append(sh[:-1])
    all_sheets = sorted({''.join(v[:2]) for v in OLD2NEW.values()}, key=lambda x: (int(x[1:-1]), x[-1]))
    uncited = [sh for sh in all_sheets if sh not in order]

    # 8. write
    with open(HERE / 'table-sources.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    status = {r[0]: (r[1], r[3], r[4], r[5], r[8]) for r in P}
    status.update({u[0]: (u[1], u[3], 'Methods', '', 'frozen GEX inventory; unchanged') for u in UP})
    with open(HERE / 'table-plan.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['table', 'table_title', 'sheet', 'sheet_title', 'v15_table', 'evidence_id', 'status', 'cited_by',
                    'supp_figures', 'sheet_files', 'sheet_rows', 'source_data_files', 'note'])
        for nt, ttl in NEW:
            sheets = []
            for old, (t2, sh, st) in OLD2NEW.items():
                if t2 == nt:
                    sheets.append((sh, st, old))
            for sh in sorted({x[0] for x in sheets}):
                olds = [x[2] for x in sheets if x[0] == sh]
                fs = [r for r in rows if r['sheet'] == nt + sh]
                sf = [r for r in fs if r['layer'] == 'sheet']
                v15 = ['S%d' % V15_TO_PLAN.get(int(o[1:]), int(o[1:])) for o in olds]
                w.writerow([nt, ttl, nt + sh, next(x[1] for x in sheets if x[0] == sh), '+'.join(v15),
                            '+'.join(status[o][0] for o in olds), '/'.join(sorted({status[o][1] for o in olds})),
                            '; '.join(status[o][2] for o in olds if status[o][2]),
                            '; '.join(status[o][3] for o in olds if status[o][3]),
                            len(sf), sum(int(r['rows'] or 0) for r in sf), len(fs) - len(sf),
                            ' '.join(status[o][4] for o in olds if status[o][4])])
    with open(HERE / 'citation-changes.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['file', 'paragraph', 'current', 'replace_with'])
        for c in cites:
            w.writerow([c['file'], c['para'], c['old'], c['new']])
    checks = dict(files=len(rows), sheet_files=sum(r['layer'] == 'sheet' for r in rows), issues=issues,
                  receipt_reads_unassigned=unassigned, main_reads_unassigned=main_unassigned,
                  first_citation_order_tables=tab_order, first_citation_order_sheets=order,
                  sheets_not_cited_in_text=uncited, citations=cites)
    json.dump(checks, open(HERE / 'checks.json', 'w'), indent=1)
    print('files', len(rows), '| sheet files', checks['sheet_files'], '| issues', len(issues),
          '| unassigned receipt/main reads', len(unassigned), len(main_unassigned))
    for x in issues[:30]:
        print('  I', x)
    print('table order:', tab_order)
    print('sheet order:', order)
    print('sheets not cited in text:', uncited)


if __name__ == '__main__':
    main()
