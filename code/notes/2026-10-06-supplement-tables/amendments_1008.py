"""Amendments after the 2026-10-08 figure/table consistency review.

Source: reviews/2026-10-08-figures-tables-consistency/REPORT.md (items P11-01, P11-04, P11-05, P11-06, P12-01,
P12-02). Author decisions (2026-10-08): withdraw S6b.2 as a whole; withdraw the four hidden_remodeling_H_* columns
of S3e.1. No value is recomputed; blocks are removed, columns are dropped or explanatory notes are added.
"""

# P11-01: the legacy inference block (status=inference, primary_test, permutation p/q) contradicts the current
# descriptive reporting (Supplementary Methods, Figure 5a, Figure S16d). S6b.1 holds every value reported now.
EXCLUDE = {
    ('S6b', 'site_statistics/region_consistency_statistics.tsv'),
}

# P12-01: H = L - D is not defined or used in the manuscript.
DROP_COLS = {
    ('S3e', 'profiles/architecture_route_summary.tsv'): (
        ['hidden_remodeling_H_mean', 'hidden_remodeling_H_median', 'hidden_remodeling_H_q25', 'hidden_remodeling_H_q75'],
        'Columns shown: the four hidden_remodeling_H_* summaries of the source file (H = L − D) are omitted; '
        'they are not used in this study.'),
}

# P11-05: explanatory notes written under a block caption (values unchanged)
NOTES = {
    ('S7b', 'trip12_evidence/figure5/site-display-evidence.tsv'): [
        'displayed_in_F5F (source column name) is True for the sites in the upper three evidence tiers '
        '(PSC_miRNA, direct_RBP_motif, transferred_RBP_P_expressed), i.e. the sites that pass the rug rule of '
        'Figure 5f; Trip12 is drawn in Figure S20h, and Figure 5f shows Riok3.'],
}

# P11-04: titles and related figures that must match the sheet content
TITLES = {
    'S2d': 'Common one-to-one ortholog background (1,306 P–S and 208 P–M genes)',
    'S4a': 'Contrast design, series and source libraries',
}
SUPP_FIGURES = {
    'S2d': 'S5a',
}

# P11-06: path rule of the Sources sheet
SOURCES_NOTE = ('Paths beginning with table/ or notes/ are relative to the project root; all other paths are relative '
                'to the analysis bundle (tables/). The inputs named for derived blocks follow the same rule.')

# P12-02: contrast codes
CONTRAST_KEY = ('Contrast codes. Every contrast is target minus reference. Target-first codes name the target state '
                'first: P_S, PS and P-S all denote P − S, and P_M, PM and P-M denote P − M (likewise R_S, E_S, M_S, '
                'E_M, E_R and R_M, e.g. E_R = E − R). Reference-first codes, used mainly in the interval analyses of Tables S6–S8, are '
                'SP = P − S, MP = P − M, SR = R − S, RE = E − R, EP = P − E and SM = M − S; semicolons separate '
                'several codes in one cell. Whether a P–S contrast uses the documented source or observed S as its '
                'reference follows the analysis, as stated in the corresponding figure legend and Supplementary Methods. '
                'Source codes are kept unchanged for traceability.')
