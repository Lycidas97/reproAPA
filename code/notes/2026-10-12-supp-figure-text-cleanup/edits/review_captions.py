TARGET = 'supplement/supplementary-figures/main.tex'
# 2026-10-12 review fixes (reviews/2026-10-12-supp-figure-revision-check, items 3 and 4):
# S9a: family-size subscript as in the display conventions (q upright, as elsewhere in the captions);
# S15h: split the nested parenthesis that carried the label, the abbreviation and the definition source.
EDITS = [
 ("with q over the 16 settings (q16).",
  "with q over the 16 settings (q$_{16}$)."),
 ("seed matches of PSC-expressed miRNAs (PSC miRNA, PSC for pluripotent stem cell, defined from external expression in Figure~\\ref{fig:supp-S17}a),",
  "seed matches of PSC-expressed miRNAs (PSC miRNA; PSC, pluripotent stem cell), defined from external expression in Figure~\\ref{fig:supp-S17}a,"),
]
