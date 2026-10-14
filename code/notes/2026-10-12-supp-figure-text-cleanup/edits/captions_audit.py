TARGET = 'supplement/supplementary-figures/main.tex'
# 2026-10-12 figure-caption audit (reviews/2026-10-12-supp-figure-caption-audit/REPORT.zh.md, items F1-F14)
EDITS = [
 # F1 S4a: dotted line in the legend ('all') was not described
 (r"within $\pm$25 to $\pm$100 nt (bar). Guan2022 series are summed.",
  r"within $\pm$25 to $\pm$100 nt (bar), with the dotted line at 1 (all). Guan2022 series are summed."),
 # F2 S5: the shared d/e drawing note sat between d and e
 ("\n\nIn d and e, points mark gene means, dark ticks medians, thick bars interquartile ranges and thin bars 10--90\\% ranges.\n\n"
  r"\textbf{e}, The three edges of that library at the primary setting (hN $-$ S, P $-$ S and hN $-$ P; hN $-$ S = hN $-$ P + P $-$ S for each gene).",
  "\n\n"
  r"\textbf{e}, The three edges of that library at the primary setting (hN $-$ S, P $-$ S and hN $-$ P; hN $-$ S = hN $-$ P + P $-$ S for each gene)."
  " In d and e, points mark gene means, dark ticks medians, thick bars interquartile ranges and thin bars 10--90\\% ranges."),
 # F3 S6e: markers not tied to a
 ("observed recurrent genes at both thresholds against the null 95\\% range.",
  "observed recurrent genes at both thresholds (markers as in a) against the null 95\\% range."),
 # F4 S6h: RWUI, PDUI and PPUI were undefined
 ("tested at well level (of 420).",
  "tested at well level (of 420). RWUI is the rank-weighted usage index, and PDUI and PPUI are the usage proportions of the most distal and the most proximal supported PAS."),
 # F5 S6i: configuration names were undefined (Supplementary Methods, external PAS-seq check)
 ("with source-choice ranges and, under each configuration name, matched/measurable gene counts.",
  r"with source-choice ranges and, under each configuration name, matched/measurable gene counts. Primary requires 20 shared counts with three controls per candidate, Count $\geq$ 50 and Count $\geq$ 100 raise the count requirement, and 2 controls uses two controls per candidate."),
 # F6 S7d: unclear closing sentence
 ("GEX discoveries are a separate result.",
  "Rank-test discoveries for expression (GEX) are a separate result (Table~S8c)."),
 # F7 S9a: q16 and the bold threshold label were unexplained
 ("with q over the 16 settings. Every null",
  "with q over the 16 settings (q16). The bold 0.02 marks the primary threshold. Every null"),
 # F8 S10c: Hs and Mm in the figure
 ("Squares give the human (left) and mouse (right) direction,",
  "Squares give the human (Hs, left) and mouse (Mm, right) direction,"),
 # F9 S12: the figure names the two mouse sets, the caption did not
 ("Mouse sets: proximal P--M recurrent genes, and their subset that was also proximal in P--S with the same contributing studies.",
  "Mouse sets: proximal P--M recurrent genes (P--M recurrent proximal set), and their subset that was also proximal in P--S with the same contributing studies (dual-reference proximal set)."),
 # F10 S15h: the five rows did not match the described tiers
 ("The tiers distinguish miRNA seed matches, species-native RBP assay matrices (Direct RBP), human-derived matrices assigned to mouse RBPs (Transferred RBP), and their expression-filtered subset. "
  "Transfer requires a high-confidence one-to-one ortholog with a supporting CisBP-RNA family and RNA-binding-domain record. "
  "The P-expressed subset requires at least 2 CPM in P of both mouse studies. "
  "PSC (pluripotent stem cell) miRNA groups are defined from external expression in Figure~\\ref{fig:supp-S17}a. ",
  "The rows are seed matches of PSC-expressed miRNAs (PSC miRNA, PSC for pluripotent stem cell, defined from external expression in Figure~\\ref{fig:supp-S17}a), species-native RBP assay matrices (Direct RBP), human-derived matrices assigned to mouse RBPs expressed in P (Transferred RBP), and, in gray, seed matches of the other miRNAs (Other miRNA) and transferred-matrix sites whose RBP is not expressed in P. "
  "Transfer requires a high-confidence one-to-one ortholog with a supporting CisBP-RNA family and RNA-binding-domain record. "
  "Expression in P requires at least 2 CPM in P of both mouse studies (Figure~\\ref{fig:supp-S17}h). "),
 # F11 S17c,d: row prefixes now Hs/Mm (D24.py); define them
 ("in mouse P--S (8 genes) and P--M (24 genes), with medians.",
  "in mouse P--S (8 genes) and P--M (24 genes), with medians. Row labels in c and d give species and contrast (Hs, human; Mm, mouse)."),
 # F12 S18d and S19e: Hs and Mm
 ("for proximal responders against not-recurrent genes and against distal responders,",
  "for proximal responders against not-recurrent genes and against distal responders in human (Hs) and mouse (Mm),"),
 ("(rank of the proximal core PAS, rank difference and log anchor span) against assay-matched expectation,",
  "(rank of the proximal core PAS, rank difference and log anchor span) in human (Hs) and mouse (Mm) against assay-matched expectation,"),
 # F13 S19f: 'driver count' was a term used nowhere else
 ("In the schematic, one matched switching PAS gives a driver count of 1,",
  "In the schematic, the gain PAS is matched and the loss PAS is not (one matched switching PAS),"),
 # F14 S20a, c: fitting routes, class sizes and the arrows
 (r"\textbf{a}, Class means of $\Delta$DWUI relative to S per route, and gene-level heat maps. Other routes show projected classes.",
  r"\textbf{a}, Class means of $\Delta$DWUI relative to S per route, and gene-level heat maps in blocks of state and route (route colors above), with the number of genes per class at right. Xing2020, Guo2019--CIP and Schiebinger2019--2i are the fitting routes, and the other routes show projected classes."),
 (r"\textbf{c}, Mean adjacent-step $\Delta$DWUI, the fraction",
  r"\textbf{c}, Mean adjacent-step $\Delta$DWUI (color, with arrows for the sign), the fraction"),
]
