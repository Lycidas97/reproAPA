r"""Caption edits for the 2026-10-12 main-figure text cleanup (main/results split source).
Each replacement must match exactly once. Run without arguments for a dry run; with --apply to write."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
R = ROOT/'main/results'
F2 = R/'gene-level-responses-in-both-directions-recur-within-each-species.tex'
F3 = R/'a-proximal-response-is-shared-between-species.tex'
F4 = R/'local-3-utr-remodeling-escapes-gene-level-summaries-and-changes-predicted-regulatory-content.tex'
EDITS = [
 # 2b: the in-figure subheading 'Same-state halves, median recurrent/observed' was removed; caption gives the order and format.
 (F2, r"pseudo-contrasts between same-state halves of S or of P, which share no biological change, against the median observed P--S recurrence on the same background (Supplementary Figure~S6c,d,g,h).",
      r"pseudo-contrasts between same-state halves of S and of P (in that order), which share no biological change, against the median observed P--S recurrence on the same background (recurrent/observed; Supplementary Figure~S6c,d,g,h).",
      '2b same-state halves: order and number format'),
 # 2c: the in-figure title was removed; the caption names the gene class.
 (F2, r"\textbf{c}, Recurrent P--S genes whose primary PASs lie in one annotated terminal $3'$ UTR:",
      r"\textbf{c}, Recurrent P--S genes whose primary PASs lie in one annotated terminal $3'$ UTR (strict tandem $3'$ UTR):",
      '2c gene class named'),
 # 2d: the figure now shows q at 0.02 only; the four q at 0.05 move here.
 (F2, r"with 95\% null ranges and BH $q$ over eight tests at right.",
      r"with 95\% null ranges and, at right, BH $q$ at 0.02 over eight tests. At 0.05, $q$ is 0.0002 and 0.022 for human identity and direction and 0.0002 and 0.0008 for mouse.",
      '2d q at 0.05'),
 # 3d: the sign glyph in the heatmap cells and the gray cells were not described.
 (F3, r"and native $\Delta$DWUI per study (below $-0.3$ drawn at the scale end).",
      r"and native $\Delta$DWUI per study for P--S and P--M (cells with the sign printed; below $-0.3$ drawn at the scale end; gray, not available).",
      '3d heatmap cells'),
 # 4c: the counts in the flow boxes (51, 8) were not in the caption.
 (F4, r"\textbf{c}, Genes with $|\Delta\mathrm{DWUI}|<0.02$",
      r"\textbf{c}, The 51 genes with $|\Delta\mathrm{DWUI}|<0.02$",
      '4c count'),
 (F4, r"in both studies (filled, strict recurrent interval;",
      r"in both studies (filled, the 8 with a strict recurrent interval;",
      '4c filled count'),
 # 4e: the ΔI column sits at the right of the decomposition bars.
 (F4, r"(diamond, total; 3 and 8 library pairs) and $\Delta I$ of its 970-nt interval",
      r"(diamond, total; 3 and 8 library pairs) and, at right, $\Delta I$ of its 970-nt interval",
      '4e ΔI placement'),
]
apply = '--apply' in sys.argv
texts = {p: p.read_text() for p in {e[0] for e in EDITS}}
for p, old, new, note in EDITS:
    n = texts[p].count(old)
    print(f"[{note}] {p.name}: matches={n}")
    assert n == 1, (note, n)
    texts[p] = texts[p].replace(old, new)
if apply:
    for p, t in texts.items(): p.write_text(t)
    print('written:', [p.name for p in texts])
