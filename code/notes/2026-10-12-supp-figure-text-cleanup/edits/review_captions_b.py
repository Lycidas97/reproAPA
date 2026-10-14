TARGET = 'supplement/supplementary-figures/main.tex'
# 2026-10-12 review fix, follow-up to review_captions.py: in the S15h row list the participial "defined from ..." read as a
# further list item, so the definition becomes its own sentence ahead of the transfer and expression rules.
EDITS = [
 ("(PSC miRNA; PSC, pluripotent stem cell), defined from external expression in Figure~\\ref{fig:supp-S17}a, species-native",
  "(PSC miRNA; PSC, pluripotent stem cell), species-native"),
 ("whose RBP is not expressed in P. Transfer requires",
  "whose RBP is not expressed in P. PSC-expressed miRNAs are defined from external expression (Figure~\\ref{fig:supp-S17}a). Transfer requires"),
]
