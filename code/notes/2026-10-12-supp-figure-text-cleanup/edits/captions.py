TARGET = 'supplement/supplementary-figures/main.tex'
EDITS = [
 # S5 a: right axis is now titled 'Expanded − supported'
 ("Right, the supported-to-expanded change in each route mean, split into",
  "Right, the supported-to-expanded change in each route mean (expanded $-$ supported), split into"),
 # S6 i: grey numbers under each configuration label (former in-figure note)
 ("with source-choice ranges and matched/measurable gene counts.",
  "with source-choice ranges and, under each configuration name, matched/measurable gene counts."),
 # S7 e: subplot titles reduced to GO / KEGG
 ("\\textbf{e}, Detection probability of the corrected G2 tests against an injected odds ratio for GO and KEGG, one subplot each",
  "\\textbf{e}, Detection probability of the corrected G2 tests (L3 full null) against an injected odds ratio for GO and KEGG, one subplot each"),
 # S11 b: per-pair gene counts no longer printed above the scatter plots
 ("for every qualified gene pair at 25 nt, for P--M and P--S. Values beyond",
  "for every qualified gene pair at 25 nt (Xing2020 $\\times$ Guo2019--CIP, 189 genes; Xing2020 $\\times$ Schiebinger2019--2i, 158 genes), for P--M and P--S. Values beyond"),
 # S19 f: worked-example values of the schematic (former in-figure text)
 ("The schematic has no measured coordinates.",
  "In the schematic, one matched switching PAS gives a driver count of 1, two of five human PASs are matched (PAS fraction 2/5), and A\\_g = 4/9 = 0.44. The schematic has no measured coordinates."),
 # S20 b: formula block removed from the schematic
 ("\\textbf{b}, Path length $L$, endpoint displacement $D$ and $T=1-D/L$, and the percentage of gene--route profiles with $T\\geq t$.",
  "\\textbf{b}, Path length $L=|R-S|+|E-R|+|P-E|$ and endpoint displacement $D=|P-S|$ (absolute DWUI differences), with $T=1-D/L$, and the percentage of gene--route profiles with $T\\geq t$."),
]
