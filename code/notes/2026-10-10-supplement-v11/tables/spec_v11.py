"""v11 supplementary-table numbering: one place for the old->new sheet map, order, group titles, related figures and the
new (2026-10-10) worksheets. Imported by assemble_workbook.py, verify_tables.py, write_tex.py and ../legends/build_legends.py."""

# ------------------------------------------------------------------ old worksheet -> v11 worksheet
MAP = {**{f'S1{c}': f'S1{c}' for c in 'abcdefgh'}, **{f'S2{c}': f'S2{c}' for c in 'abcde'},
       'S3a': 'S8a', 'S3c': 'S8b', 'S3d': 'S8c', 'S3e': 'S8d',
       'S4a': 'S3a', 'S4b': 'S3b', 'S4c': 'S3c', 'S4d': 'S3d', 'S4f': 'S3e', 'S4g': 'S3f',
       **{f'S5{c}': f'S4{c}' for c in 'abcde'}, **{f'S6{c}': f'S5{c}' for c in 'abcdefg'},
       'S7a': 'S6b', 'S7b': 'S6c', 'S7c': 'S6d', 'S8a': 'S9a', 'S8b': 'S9b',
       'S9b': 'S7g', 'S9e': 'S6a', 'S9f': 'S10a'}
# old worksheets not carried into v11 (pruned 2026-10-11): no v11 claim, legend or text citation, and their figures were dropped
DROP_OLD = {'S3b': 'route hold-out of the profile classes (former Figure 2e, S6d,e)',
            'S3f': 'Moran structure and APA-GEX profile agreement (former S8)',
            'S4e': 'dual-reference candidates (former S12a-e)',
            'S9a': 'dual-reference conditional nulls', 'S9c': 'cross-species dual-reference comparison and power bounds',
            'S9d': 'profile-class entry into endpoint recurrence (former Figure 2d)'}
# blocks dropped from carried-over worksheets (old sheet -> old block numbers); kept blocks are renumbered 1..k
DROP_BLOCKS = {'S4c': [1, 3, 4, 5, 6],                         # intermediate-state analysis (former S11); keeps the direction counts of S6b
               'S3d': [1],                                     # class-wise ORA (former S9a-c); keeps the rank tests of S7d
               'S9b': [1, 2, 3, 6, 8, 22, 23, 24, 25]}         # dual-reference blocks; keeps cycle, common-support and P-M candidate blocks
TITLE_OVERRIDE = {'S3c': 'Direction composition of recurrent genes',
                  'S8c': 'Rank tests of GO and KEGG terms on APA and expression'}
OLD_TABLE = {'S1': 'S1', 'S2': 'S2', 'S3': 'S8', 'S4': 'S3', 'S5': 'S4', 'S6': 'S5', 'S7': 'S6', 'S8': 'S9'}   # whole-table refs (S9 split)

GROUPS = {'S1': 'Datasets, libraries and expression states', 'S2': 'PAS measurement and route-level PAS position',
          'S3': 'Within-species recurrence and its context', 'S4': 'Cross-species comparisons and cis-feature models',
          'S5': 'Interval inclusion and predicted-site retention', 'S6': 'Low-net-shift genes and interval examples',
          'S7': 'Candidate explanations of the shared proximal responses', 'S8': 'State-resolved APA profile classes',
          'S9': 'Main analyses repeated with alternative PAS calling', 'S10': 'Reference-state sensitivity'}

# related figures (v11 numbers); bare S-numbers are supplementary figures
REL = {'S1a': 'S1–S3', 'S1h': 'S2', 'S2a': 'S4', 'S2b': 'S5a,b', 'S2c': 'S5c', 'S2d': 'S5a', 'S2e': 'S5d–f', 'S2f': 'Figure 2e',
       'S3a': 'S3c; S6a–f', 'S3b': 'S6a,e,f', 'S3c': 'S6b', 'S3d': 'S6c,d,g–i', 'S3e': 'S8a–d', 'S3f': 'S7a–c,e',
       'S3g': 'Figure 2a', 'S3h': 'Figure 2c',
       'S4a': 'S9b–f', 'S4b': 'S9a', 'S4c': 'S11a–e', 'S4d': 'S11f,g', 'S4e': 'S19f–i', 'S4f': 'Figure 3a,b', 'S4g': 'Figure 3c; S10a–c',
       'S4h': 'S10d,e', 'S4i': 'S12',
       'S5a': 'S13a', 'S5b': 'S13c,d', 'S5c': 'S13e,f', 'S5d': 'S16a,b; S17a–f', 'S5e': 'S17g,h', 'S5f': 'S16c–f',
       'S6a': 'Figure 4b', 'S6b': 'S13b; S14', 'S6c': 'S15', 'S6d': 'Figure 4d; S21',
       'S7a': 'S18a,b,e', 'S7b': 'S18c', 'S7c': 'S18d', 'S7d': 'S19a', 'S7e': 'S19b', 'S7f': 'S19c–e',
       'S8a': 'S20a', 'S8b': 'S20a', 'S8c': 'S7d', 'S8d': 'S20b–d', 'S9a': 'S21', 'S9b': 'S21'}
# worksheets specified after the initial results (marked † in Contents); value = analysis descriptor kept from the source sheet
POSTHOC_OLD = {'S6a': 'Interval and gene-wide comparison', 'S7g': 'P–M and P/M/source-S cycle diagnostics',
               'S10a': 'Reference-state sensitivity of P–M comparisons'}

# ------------------------------------------------------------------ new worksheets (all specified after the initial results)
N = 'notes/2026-10-10-'
ALL = 'All rows; no statistical recomputation.'
NEW = {
 'S2f': dict(title='Route-mean composition by sign', blocks=[
   dict(t='Route means split into proximal and distal contributions', f='notes/2026-10-10-supplement-v11/tables/derived/route_mean_composition.tsv', q=None,
        rule='Values plotted in Figure 2e, computed from the frozen global tables by tables/derive_composition.py and checked against global_direction_summary (Table S2b.3).',
        note='P–S on the 1,306 one-to-one orthologs common to all measurement settings (supported50_p08), each gene averaged over sampled dates. distal_arm and proximal_arm are the sums of positive and of negative ΔDWUI divided by n_genes, so together they give mean_dDWUI. Responders have |ΔDWUI| ≥ 0.02. Species rows weight studies equally after averaging the two Guan2022 series, and responder percentiles are averaged across studies. The human_minus_mouse row splits the species difference into the two arms.')]),
 'S3g': dict(title='Direction-specific recurrence', blocks=[
   dict(t='Proximal and distal recurrent genes against the measurement-matched null', f=N+'directional-recurrence/directional_recurrence.tsv', q=None, rule=ALL,
        note='Recurrent genes counted separately by direction for each species, reference and threshold, against the L3 null used for recurrence, which permutes study-level events within strata of expression, APA coverage, PAS number, reference position and room to shift (10,000 permutations, one-sided). q_threshold_8: BH over the eight species–reference–direction tests at one threshold. q_joint_16: BH over both thresholds.'),
   dict(t='Distal share of the recurrent sets', f=N+'directional-recurrence/direction_composition.tsv', q=None, rule=ALL,
        note='Distal share, distal/(proximal + distal), against the same permutations by a two-sided equal-tailed test. q_threshold_4: BH over the four species–reference sets at one threshold. q_joint_8: BH over both thresholds.')]),
 'S3h': dict(title='Strict tandem-3′ UTR subset', blocks=[
   dict(t='Shortening and lengthening among recurrent tandem genes', f=N+'recurrence-utr-ortholog/tandem_directional_recurrence.tsv', q=None, rule=ALL,
        note='A recurrent gene entered the subset when its whole primary PAS axis lay in one annotated terminal 3′ UTR of a transcript compatible across all contributing studies. Proximal (shortening) and distal (lengthening) counts are tested against the L3 null restricted to the subset, which permutes study-level events within each study and remaining stratum (10,000 permutations, one-sided). q_eight_per_threshold: BH over the eight direction tests at one threshold. q_joint_16: BH over both thresholds.'),
   dict(t='Usage-weighted change in retained 3′ UTR length', f=N+'recurrence-utr-ortholog/tandem_recurrent_length.tsv', q=None, rule=ALL,
        note='Change in usage-weighted retained 3′ UTR length (nt) of the recurrent tandem genes per route, and as a gene-equal, study-equal summary. n_opposing counts genes whose own route shifted against the recurrent direction. Descriptive.')]),
 'S4f': dict(title='Direction pairing between species', blocks=[
   dict(t='Direction pairing at 0.02 and 0.05', f=N+'recurrence-utr-ortholog/ortholog_3x3.tsv', q=None, rule=ALL,
        note='Each one-to-one ortholog measurable in at least two studies per species is labelled proximal, distal or not recurrent in each species (P–S 1,656 and P–M 751 orthologs). Observed counts per cell of the 3 × 3 table against a joint null that permutes labels within each species, within strata of measurement and PAS capacity (measurement_PAS_capacity) or of assay features (full_assay), 10,000 permutations. q_eight_per_threshold_model: BH over the eight tested pairings per threshold and model. Rows with tested == False are cells outside the tested family.'),
   dict(t='Threshold scan from 0.010 to 0.050', f=N+'ortholog-threshold-scan/ortholog_3x3_threshold_scan.tsv', q=None, rule=ALL,
        note='The same pairings at nine thresholds, 100,000 permutations each. q_scan_72_per_model: BH over the 72 pairing–threshold tests of one model. q_scan_144_across_models: BH over both models. Monte Carlo standard errors and Wilson intervals of the tail probability are listed, and permutation_p_floor is the smallest attainable p.')]),
 'S4g': dict(title='Event classes of paired responses', blocks=[
   dict(t='Event classes of the P–M ortholog pairs', f=N+'utr-event-stratification/classification/species_class_matrix.tsv', q=None, rule=ALL,
        note='Each species\' gene is classed from the stored route events: tandem when, in every evaluable route, the switching PAS pair (the PASs with the largest usage gain and loss) lay in one annotated terminal 3′ UTR (a property of this pair, not the whole-axis strict tandem subset), ALE when every route\'s switching PAS pair was annotation-compatible with alternative terminal exons, and otherwise mixed or one of two minor classes. Counts over the 751-pair background.'),
   dict(t='Direction pairing within pairs classed alike', f=N+'utr-event-stratification/overlap/directional_tests.tsv', q=None, rule=ALL,
        note='Proximal and distal pairings restricted to pairs that are tandem in both species (both_tandem) or ALE in both (both_ALE), with both null models of Table S4f and nine thresholds. q_primary4: BH over two classes and two models at 0.02. q_scan36: BH over the scan. The sensitivity and secondary q columns name their own families.'),
   dict(t='Orthologs proximal in both species, all thresholds', f=N+'utr-event-stratification/classification/shared_proximal_core_all_thresholds.tsv', q=None, rule=ALL,
        note='P–M orthologs proximal in both species at each threshold (27 at 0.02), with per-species study support, event class, switching-PAS-pair classes and route-level ΔDWUI summaries.')]),
 'S4h': dict(title='Corresponding PAS switches and homologous seed matches', blocks=[
   dict(t='Attrition of corresponding PAS switches', f=N+'conserved-response/sites/event_attrition.tsv', q=None, rule=ALL,
        note='Contexts and unique orthologs left after each criterion for a corresponding PAS switch, per gene set and matching tolerance (25 nt, and 50 nt for sensitivity): reciprocal matching of human and mouse PASs after liftOver, at least 80% of usage on matched sites in every study and state, corresponding gaining and losing PASs, and support from two studies per species.'),
   dict(t='Changed intervals of the switches that met all criteria', f=N+'conserved-response/sites/strict_case_intervals.tsv', q=None, rule=ALL,
        note='Interval between the corresponding proximal and distal PASs in each species, with pooled inclusion in the target and reference states.'),
   dict(t='Interval inclusion per study and series', f=N+'conserved-response/sites/strict_case_interval_inclusion.tsv', q=None, rule=ALL,
        note='Inclusion of the same intervals in target and reference states of each study and series.'),
   dict(t='Predicted miRNA seed matches at homologous positions', f=N+'conserved-response/sites/reciprocal_predicted_seed_pairs.tsv', q=None, rule=ALL,
        note='Predicted seed matches in a changed interval whose whole seed maps gap-free, uniquely and reciprocally into the other species\' interval. status gives the outcome of each query. Predicted sites do not establish binding.')]),
 'S4i': dict(title='Naive-state specificity of the shared proximal set', blocks=[
   dict(t='Direct comparisons of naive and control enrichment', f=N+'proximal-module-specificity/direct_specificity_results.tsv', q=None, rule=ALL,
        note='For each mouse target set, the statistic is O_naive/E_naive − O_control/E_control on one common one-to-one ortholog background, one joint stratification and one set of joint permutation draws, so response sets of different size are compared on the same null. One-sided p. q_specificity14_within_model_threshold: BH over the 14 comparisons of one model and threshold. q_all_specificity_exploratory: BH over all rows. Human naive expression is a pseudobulk proxy recovered from the frozen expression matrices.')]),
 'S7a': dict(title='Cell-cycle scores and within-state PAS position', blocks=[
   dict(t='Cell-cycle score contrasts between states', f=N+'proliferation-apa/state_score_contrasts.tsv', q=None, rule=ALL,
        note='Cycle_S and Cycle_G2M are the mean ln(1 + counts per 10,000) of 10 S-phase and 8 G2/M markers, and Cycle is their equal-weight mean. Adjusted_* subtracts, for each marker, the mean of 20 expression-matched non-cycle genes. P–M is compared within the same physical libraries (paired, libraries weighted equally) and P–S or M–S against the documented source or observed source-like cells (unpaired). Descriptive, with library-bootstrap ranges.'),
   dict(t='Per-library score differences', f=N+'proliferation-apa/paired_library_score_differences.tsv', q=None, rule=ALL,
        note='Library-level differences behind the paired contrasts.'),
   dict(t='Within-state association between cell-cycle score and mean DWUI', f=N+'proliferation-apa/route_associations.tsv', q=None, rule=ALL,
        note='Partial correlation between cell-cycle score and mean DWUI over metacells, on a gene panel fixed within each library and state, controlling for panel expression, PAS depth, total expression and PAS counts per cell and cell number. Summarized per route, state scope and APA variant with library-bootstrap and leave-one-library ranges. Descriptive. Cell-cycle scores are transcriptional programs, not measured division.')]),
 'S7b': dict(title='Core CPA-factor expression', blocks=[
   dict(t='Expression of 11 core cleavage and polyadenylation factors', f=N+'cpa-factor-expression/expression_contrasts.tsv', q=None, rule=ALL,
        note='Difference in mean log2(counts per million + 1), with counts per million over all genes, for P–M in the same libraries, for P–S against the documented source and for the listed source-like and naive contrasts. Descriptive. The RNA-direction columns compare the sign of the factor change with the net ΔDWUI of the route. RNA levels do not establish factor activity.')]),
 'S7c': dict(title='Upstream UGUA signature', blocks=[
   dict(t='Upstream motif counts of proximal responders', f=N+'cfim-sequence-signature/sequence_signature_tests.tsv', q=None, rule=ALL,
        note='Distal minus proximal motif count in windows upstream of the cleavage sites for proximal responders against not-recurrent genes or against distal responders, on the one-to-one orthologs eligible in both species (paired_strict) and on each species\' strict tandem-3′ UTR background (species_strict). The null permutes direction labels within strata of measurement, PAS capacity and number, baseline DWUI, PAS spacing, GC content, expression, APA coverage and core hexamer. family primary_8: BH over the eight primary TGTA tests (−100 to −40 nt). exploratory_208: other windows, control motifs (TGAT, CGTA) and models. No CFIm activity is measured.')]),
 'S7d': dict(title='S–M contrast on the shared proximal core', blocks=[
   dict(t='Shared core in S relative to M', f=N+'cycle-source-attribution/sm-core/SM_route_summary.tsv', q=None, rule=ALL,
        note='Contrast S_M = S − M, a contrast without P cells, on the 27-gene shared core (core27) and on all P–M proximal genes (complete_PM_proximal). Negative mean_DWUI means more proximal in S. Cell-cycle score differences (S − M) are given on the same contributors. Designs: main, late M only, same library and same time across libraries. Descriptive. Rows with status not_estimable have no contributing library and are not zero effects.')]),
 'S7e': dict(title='Within-date cell-cycle increments in Schiebinger2019–2i', blocks=[
   dict(t='Fixed gene panels', f=N+'cycle-source-attribution/wot-dose/fixed_panel_descriptive_statistics.tsv', q=None, rule=ALL,
        note='Correlation across wells of the same date between the cell-cycle score increment and proximal strength for three fixed panels (native85, ortholog45 and mouse_proximal8, 11 dates), with leave-one-date-out ranges. Descriptive.'),
   dict(t='Core subset of the exact test', f=N+'cycle-source-attribution/wot-dose/core_formal_statistics_replay.tsv', q=None, rule=ALL,
        note='Exact two-sided test on the core subset (7 dates) for the cell-cycle and CPA predictors. q2_BH: BH over the two predictors.')]),
 'S7f': dict(title='Candidate mechanisms and functional categories', blocks=[
   dict(t='Candidate mechanisms for the shared proximal set', f=N+'shared-proximal-mechanisms/mechanism_summary.tsv', q=None, rule=ALL,
        note='ESCC miRNA seed sites in distal regions, external QKI-knockdown and EMT targets, membrane-system GO categories, U1 perturbation (GSE140543, pooled tracks without treatment-level p) and PAS anchor architecture, each against its stated expectation. family names the BH family. conclusion is the analysis label, in which not_supported means not enriched relative to expectation.'),
   dict(t='Functional categories of same-direction orthologs', f=N+'conserved-response/function/category_results.tsv', q=None, rule=ALL,
        note='GO categories among the listed ortholog sets against matched expectation with null ranges. q_primary6, q_sensitivity138 and q_joint144 name their BH families.')]),
 'S10b': dict(title='Shared proximal set under marker-defined M', blocks=[
   dict(t='Direction pairing with original and marker-defined M', f=N+'marker-ortholog-sharing/analysis/directional_tests.tsv', q=None, rule=ALL,
        note='P–M direction pairing of one-to-one orthologs on the 728 orthologs measurable under all three definitions, with P fixed: original M (original_matched), marker-defined M* in human only (human_Mstar, the primary analysis) and in both species (both_Mstar). Labels are permuted within strata of measurement and PAS capacity (measurement_PAS_capacity) or of assay features (full_assay), 100,000 permutations per definition, threshold and model; one-sided upper-tail p = (1 + tail)/100,001. q8: BH over the eight direction cells of one definition, threshold and model. q72: over the nine thresholds of one definition and model. q144: over both models of one definition. q288: over all unique tests.'),
   dict(t='Shared proximal members', f=N+'marker-ortholog-sharing/analysis/shared_proximal_primary_members.tsv', q=None, rule=ALL,
        note='Orthologs proximal in both species at |ΔDWUI| ≥ 0.02 under each definition, with per-species study support. original_core27 marks members of the original 27-pair shared set.'),
   dict(t='The original 27 shared pairs under each definition', f=N+'marker-ortholog-sharing/analysis/original_core27_flow.tsv', q=None, rule=ALL,
        note='Status of the original 27 shared proximal pairs for each definition and threshold. Pairs no longer measurable in at least two human studies (GNPNAT1, KLC1) are not evaluable, which is not the same as no response.')]),
}

# ------------------------------------------------------------------ blocks appended to carried-over worksheets (2026-10-11)
# Results of the frozen polyApipe v2-primary run (same bundle as the existing blocks) that the Supplementary Methods
# assign to the worksheet but that the frozen workbook did not carry. Existing rows stay byte-identical; the new
# blocks follow the last block, keep the worksheet's label, header and row styles, and add one note row.
M = 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/motif_controls/'
APPEND = {
 'S5d': [
   dict(name='grid_summary', f=M + 'grid_summary.tsv',
        note='Matching grid. Each changed interval was rematched to an equal-length stable window in the same gene under 12 settings that cross a GC-difference limit of 0.05, 0.10, 0.20 or none (GCinf) with a relative-center-distance limit of 0.20, 0.40 or none (POSinf). The dynamic cohort uses every interval matched at a setting, and fixed_baseline_windows keeps the intervals matched under the primary limits (GC 0.05, position 0.20). RBP_transfer counts hits of human-derived matrices assigned to one-to-one mouse orthologs (Table S5g). Gene-equal mean differences are tested by within-gene sign flips (exact with at most 20 nonzero genes, otherwise 20,000 draws), and lo and hi come from 5,000 gene-bootstrap draws. loo_min and loo_max are leave-one-gene-out means. Cohorts with fewer than ten genes are descriptive (p = 1). q_all_settings: BH within one cohort over the 96 direct positions (12 settings × eight species–contrast–category positions, resource_family direct96) or over the 24 transferred positions (transfer24). q_within_setting: BH over the positions of one setting.'),
   dict(name='background_ladder', f=M + 'background_ladder.tsv',
        note='Random-priority backgrounds. For each target interval, one window was drawn by random priority among the equal-length stable windows of the same gene (same_gene) or of a donor gene (cross_gene), either unconstrained (length), within GC difference 0.05 (GC), or within GC difference 0.05 and relative-center distance 0.20 (GC_position). Each donor is a unique gene outside the changed set whose stable sequence can match every target interval of the gene, assigned by maximum bipartite matching. six_background_common holds the genes with a donor and all six backgrounds (48 tests). same_gene_original_fixed repeats the three same-gene backgrounds on all target genes (24 tests). Tests and the descriptive rule follow S5d.6. q_all_backgrounds: BH over the tests of one cohort (48 or 24). q_eight: BH over the eight species–contrast–category positions of one background.'),
   dict(name='motif_power', f=M + 'motif_power.tsv',
        note='Detection power of the matched-window test (Figure S16b). For each cohort, 1,000 datasets give every gene its absolute observed difference with a random sign and add 0, 2, 5, 10 or 20 sites/kb to all genes (all) or to a fixed quarter of them (quarter). Each dataset is tested with 4,095 random sign flips at alpha 0.05 and 0.05/8. Cohorts with fewer than ten genes are not tested, so their rejection rate is zero by design. lo and hi are Clopper–Pearson 95% intervals of the rejection rate.'),
 ]}

ORDER = (['Contents'] + [f'S1{c}' for c in 'abcdefgh'] + [f'S2{c}' for c in 'abcdef'] + [f'S3{c}' for c in 'abcdefgh']
         + [f'S4{c}' for c in 'abcdefghi'] + [f'S5{c}' for c in 'abcdefg'] + [f'S6{c}' for c in 'abcd']
         + [f'S7{c}' for c in 'abcdefg'] + [f'S8{c}' for c in 'abcd'] + ['S9a', 'S9b', 'S10a', 'S10b', 'Sources'])
TAB = {1: 'FF23415E', 0: 'FF8FA9C4'}

def table_of(sheet): return sheet[:-1]
def caption(sheet):
    t = table_of(sheet); parts = [f'Part of Table {t}: {GROUPS[t]}']
    if sheet in NEW or sheet in POSTHOC_OLD: parts.append('Specified after the initial results')
    rel = REL.get(sheet, '')
    if rel:
        parts.append(('Related figures: ' if 'Figure' in rel else 'Related supplementary figures: ')+rel)
    elif sheet in POSTHOC_OLD: parts.append(POSTHOC_OLD[sheet])
    return '  ·  '.join(parts)

assert sorted(set(MAP.values()) | set(NEW)) == sorted(s for s in ORDER if s not in ('Contents', 'Sources')), 'order mismatch'
assert not set(DROP_OLD) & set(MAP) and set(DROP_BLOCKS) <= set(MAP) and set(TITLE_OVERRIDE) <= set(MAP.values())
assert set(REL) <= set(ORDER) and set(POSTHOC_OLD) <= set(MAP.values()), 'REL/POSTHOC keys outside ORDER'
assert set(APPEND) <= set(MAP.values()) and not set(APPEND) & set(NEW) and not {INV for INV, v in MAP.items() if v in APPEND} & set(DROP_BLOCKS), 'APPEND targets'
