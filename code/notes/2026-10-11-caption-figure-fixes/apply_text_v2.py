"""Apply caption and citation revisions (2026-10-11), v2.
Continues apply_text.py (never run). Changes from v1, all to keep pace with the 15:51 long-sentence split and the
house style (no semicolons outside parentheses, sentences <= 40 words):
  E0  the Fig. 1c clause becomes its own sentence after the 'Because this axis...' sentence.
  E2  the 2b-lower insert is three sentences instead of one.
  E9  matched to the split text ('...S11). Only the alternative PAS calls...'); the E-M control is two sentences and
      gives q under both nulls (Table S4c.1 capacity, S4c.2 coarse, 25 nt): E-M rho 0.181, q8 0.041/0.024;
      P-M rho 0.077, q8 0.275/0.294.
  E10 Riok3 first (letters d then e), Fgfr1op2 rephrased into two sentences.
  E12 the human P-M count (14 of 20, Figure 4f) is its own sentence.
  Captions: 3e square key names the two columns; 'bars, s.e.' -> 'error bars, s.e.'.
Every old string must occur exactly once. Writes a JSON log to out/apply_text_v2.log.json."""
import sys, json, hashlib
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from captions import CAP
ROOT = HERE.parents[1]
R = 'main/results/'
F1 = R+'the-global-apa-direction-diverges-between-human-and-mouse-routes.tex'
F2 = R+'gene-level-responses-in-both-directions-recur-within-each-species.tex'
F3 = R+'the-divergence-lies-mainly-in-the-breadth-of-the-proximal-response.tex'
F4 = R+'a-proximal-response-is-shared-between-species.tex'
F5 = R+'local-3-utr-remodeling-escapes-gene-level-summaries-and-changes-predicted-regulatory-content.tex'
fig1 = r'Fig.~\ref{fig:reprogramming-overview}'; fig2 = r'Fig.~\ref{fig:apa-recurrence-context}'
fig3 = r'Fig.~\ref{fig:cross-species-apa}'; fig4 = r'Fig.~\ref{fig:apa-regulatory-sites}'

CAP = dict(CAP)
for old, new in [(r"Squares mark $q_8<0.05$ under a coarse and a capacity-conditioned null.",
                  r"Filled squares mark $q_8<0.05$ under the coarse (left) and capacity-conditioned (right) nulls."),
                 (r"per study pair (bars, s.e.).", r"per study pair (error bars, s.e.).")]:
    assert CAP[3].count(old) == 1, old
    CAP[3] = CAP[3].replace(old, new)

EDITS = [
 # Figure 1
 (F1, r"spanning its primary-catalog PASs (" + fig1 + r"c; Materials and Methods). Because this axis can include terminal exons of shorter isoforms, DWUI measures relative PAS position and tracks $3'$ UTR length only in tandem-$3'$ UTR genes.",
      r"spanning its primary-catalog PASs (Materials and Methods). Because this axis can include terminal exons of shorter isoforms, DWUI measures relative PAS position and tracks $3'$ UTR length only in tandem-$3'$ UTR genes. Within each study, gene-centered DWUI differed between expression states (" + fig1 + r"c)."),
 (F1, r"P was proximal or near zero relative to observed source-like cells (Supplementary Figure~S5f and Table~S2e).",
      r"P was proximal or near zero relative to observed source-like cells (" + fig1 + r"d, upper; Supplementary Figure~S5f and Table~S2e)."),
 # Figure 2 (rows reordered: old d -> c, old e -> d, old c -> e)
 (F2, r"(observed/null-mean ratios 1.34--1.84; " + fig2 + r"b).",
      r"(observed/null-mean ratios 1.34--1.84; " + fig2 + r"b). Most recurrent P--S genes kept their direction when bootstrap intervals had to exclude zero (386 of 483 human and 513 of 604 mouse genes). Pseudo-contrasts between same-state halves of S or of P, which share no biological change, gave medians of only 10 to 45 recurrent genes against 404 to 515 observed on the same backgrounds. In Schiebinger2019--2i wells, 134 of 143 tested recurrent P--M members reached well-level $q<0.05$, against 370 of all 614 tested genes (" + fig2 + r"b, lower)."),
 (F2, r"$q_8=0.00013$ each; " + fig2 + r"d).", r"$q_8=0.00013$ each; " + fig2 + r"c)."),
 (F2, r"while gene expression (GEX) changed in opposite directions (" + fig2 + r"e).",
      r"while gene expression (GEX) changed in opposite directions (" + fig2 + r"d)."),
 (F3, r"lay mainly on the proximal side (" + fig2 + r"c).", r"lay mainly on the proximal side (" + fig2 + r"e)."),
 (F3, r"with more proximal responders in mouse routes (" + fig2 + r"c).", r"with more proximal responders in mouse routes (" + fig2 + r"e)."),
 # Figure 3 (letters d <-> e)
 (F4, r"shared proximal shifts exceeded expectation at eight, shared distal shifts at four and reversals at none (" + fig3 + r"b).",
      r"shared proximal shifts in P relative to M exceeded expectation at eight, shared distal shifts from source to P at four and reversals in either comparison at none (" + fig3 + r"b)."),
 (F4, r"loses all three C-terminal LIM domains (" + fig3 + r"e;", r"loses all three C-terminal LIM domains (" + fig3 + r"d;"),
 (F4, r"difference 0.005, $-0.145$ to 0.144; " + fig3 + r"d; Supplementary Figure~S11). Only the alternative PAS calls separated the two (Supplementary Figure~S21c).",
      r"difference 0.005, $-0.145$ to 0.144; " + fig3 + r"e, upper; Supplementary Figure~S11). Only the alternative PAS calls separated the two (Supplementary Figure~S21c). In Liu2020 paired with Guo2019--CIP (116 genes), the two species also responded alike in E relative to M, another comparison without pluripotent cells ($\rho=0.181$, $q_8=0.041$ and 0.024 under the coarse and capacity-conditioned nulls). P relative to M did not correspond in this pair ($\rho=0.077$, $q_8\geq0.27$; " + fig3 + r"e, lower)."),
 # Figure 4 (Riok3 before Fgfr1op2 in the text)
 (F5, r"""Remodeling elsewhere in a gene could even move its overall PAS position opposite to the change at the $3'$ end. \textit{Fgfr1op2} counted as a distal recurrent gene in both mouse studies ($\Delta$DWUI $+0.033$ and $+0.024$) although its 970-nt terminal interval lost 23 and 31 pp of inclusion in P (""" + fig4 + r"""e).""",
      r"""Remodeling elsewhere in a gene could even move its overall PAS position opposite to the change at the $3'$ end. In \textit{Riok3}, an 889-nt interval lost about 15 pp of inclusion in both mouse studies, but lower P usage of upstream PASs made the gene-wide effect $+0.0247$ in Guo2019--CIP (""" + fig4 + r"""d; Supplementary Table~S6d). \textit{Fgfr1op2} showed this opposition in both mouse studies. It counted as a distal recurrent gene ($\Delta$DWUI $+0.033$ and $+0.024$) although its 970-nt terminal interval lost 23 and 31 pp of inclusion in P (""" + fig4 + r"""e)."""),
 (F5, r""" This split is exact arithmetic and agnostic to mechanism. Similarly, in \textit{Riok3}, an 889-nt interval lost about 15 pp of inclusion in both mouse studies, but lower P usage of upstream PASs made the gene-wide effect $+0.0247$ in Guo2019--CIP (""" + fig4 + r"""d).""",
      r""" This split is exact arithmetic and agnostic to mechanism."""),
 (F5, r"241 of 246 mouse P--M genes and more in 48 of 60 human observed-S P--S genes (net retention, the change in complete predicted-site content weighted by interval inclusion; Materials and Methods).",
      r"241 of 246 mouse P--M genes (net retention, the change in complete predicted-site content weighted by interval inclusion; Materials and Methods). Human transcripts retained more in 48 of 60 observed-S P--S genes but fewer in 14 of 20 P--M genes."),
 # Supplement pointers to relabelled panels
 ('supplement/supplementary-figures/main.tex', r"as 2 $\times$ 2 tiles as in Figure~2e:", r"as 2 $\times$ 2 tiles as in Figure~2d:"),
 ('supplement/supplementary-figures/main.tex', r"$\times$1.74, as in Figure~2e)", r"$\times$1.74, as in Figure~2d)"),
 ('supplement/supplementary-figures/main.tex', r"shown as tiles as in Figure~3d (filled", r"shown as tiles as in Figure~3e (filled"),
 ('supplement/supplementary-tables/main.tex', r"Route-mean composition by sign (Figure~2c).", r"Route-mean composition by sign (Figure~2e)."),
 ('supplement/supplementary-tables/main.tex', r"Strict tandem-3$'$ UTR subset (Figure~2d).", r"Strict tandem-3$'$ UTR subset (Figure~2c)."),
 ('supplement/supplementary-methods/cross-species-comparisons.tex', r"native $\Delta$DWUI of the gene (Figure~3e)", r"native $\Delta$DWUI of the gene (Figure~3d)"),
]

def main(dry):
    texts, before, log = {}, {}, []
    for f, old, new in EDITS:
        if f not in texts:
            texts[f] = before[f] = (ROOT/f).read_text()
        n = texts[f].count(old); assert n == 1, (f, n, old[:100])
        texts[f] = texts[f].replace(old, new); log.append(dict(file=f, old=old, new=new))
    for fig, f in [(1, F1), (2, F2), (3, F4), (4, F5)]:
        if f not in texts:
            texts[f] = before[f] = (ROOT/f).read_text()
        lines = texts[f].split('\n')
        idx = [i for i, l in enumerate(lines) if l.startswith(r'\mainfigure{Figure.%d.pdf}' % fig)]
        assert len(idx) == 1, (fig, idx)
        log.append(dict(file=f, caption=fig, old=lines[idx[0]], new=CAP[fig]))
        lines[idx[0]] = CAP[fig]; texts[f] = '\n'.join(lines)
    if dry:
        print('dry run ok:', len(EDITS), 'edits,', 4, 'captions,', len(texts), 'files'); return
    for f, s in texts.items():
        assert (ROOT/f).read_text() == before[f], ('file changed during run', f)
        (ROOT/f).write_text(s)
    rec = dict(edits=log, files={f: hashlib.sha256(s.encode()).hexdigest() for f, s in texts.items()})
    (HERE/'out'/'apply_text_v2.log.json').write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print('written:', *texts, sep='\n  ')

if __name__ == '__main__':
    main(dry='--dry-run' in sys.argv)
