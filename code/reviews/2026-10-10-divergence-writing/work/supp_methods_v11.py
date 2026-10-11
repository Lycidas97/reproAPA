"""Revise supplement/supplementary-methods/*.tex to v11 numbering in drafts/supp-methods-v11/ (production untouched).
Explicit edits use @..@ markers so the generic renumbering pass does not touch numbers that are already v11."""
import re, sys, importlib.util, shutil
from pathlib import Path
RUN = Path(__file__).resolve().parents[1]; ROOT = RUN.parents[1]
SRC = ROOT / 'revision-backups/20261011-before-v11-apply/supplement/supplementary-methods'; DST = RUN / 'drafts/supp-methods-v11'
sp = importlib.util.spec_from_file_location('s', ROOT / 'notes/2026-10-10-supplement-v11/tables/spec_v11.py')
spec = importlib.util.module_from_spec(sp); sp.loader.exec_module(spec)
FIG = {1:1,2:2,3:3,4:4,5:5,10:6,13:9,14:11,16:13,17:16,18:17,19:14,20:15,21:21}

def drop_sub(t, title):
    i = t.index('\\subsubsection{' + title + '}')
    j = t.find('\\subsubsection{', i + 10)
    return t[:i] + (t[j:] if j >= 0 else '')
def rep(t, a, b, f):
    assert t.count(a) == 1, (f, t.count(a), a[:80]); return t.replace(a, b)
def insert_after_sub(t, title, block):
    i = t.index('\\subsubsection{' + title + '}'); j = t.find('\\subsubsection{', i + 10)
    if j < 0: return t.rstrip('\n') + '\n\n' + block.strip() + '\n'
    return t[:j] + block.strip() + '\n\n' + t[j:]

# additions (already v11) -> protect numbers
add = (RUN / 'drafts/supp-methods-additions-plan-1.tex').read_text()
blocks = {}
for m in re.finditer(r'^% (s\d) -> [^\n]*\n(?:%[^\n]*\n)*(.*?)(?=^% s\d -> |\Z)', add, re.S | re.M):
    b = m.group(2).strip()
    b = re.sub(r'fig:supp-S(\d+)\}', r'fig:supp-@S\1@}', b)
    b = re.sub(r'(Tables?~)(S\d+[a-z]?)((?:(?:, | and |--)S\d+[a-z]?)*)', lambda m: m.group(1) + re.sub(r'(S\d+[a-z]?)', r'@\1@', m.group(2) + m.group(3)), b)
    blocks[m.group(1)] = b
assert sorted(blocks) == ['s1','s2','s3','s4','s5','s6'], blocks.keys()

E = {
'dwui-metacells-and-state-comparisons.tex': [
  ('The overview of Figure~1d (upper) and Figure~2 applied', 'The overview of Figure~@1@d (upper) and the state-resolved profiles (Figure~\\ref{fig:supp-@S20@}) applied'),
  ('(Figure~1d, upper; Figure~', '(Figure~@1@d, upper; Figure~'),
  ('(Figure~1d, lower; Figure~', '(Figure~@1@d, lower; Figure~')],
'state-resolved-profiles.tex': [
  ('The silhouette was highest at $k=2$.', 'The silhouette was highest at $k=2$ (Table~@S8a@).'),
  ('(Figure~\\ref{fig:supp-S6}a).', '(Figure~\\ref{fig:supp-@S20@}a).'),
  ('\\subsubsection{Path geometry and the fixed-class null}', '\\subsubsection{Path geometry}'),
  ('Route summaries are medians with interquartile ranges over genes.', 'Route summaries are medians with interquartile ranges over genes (Figure~\\ref{fig:supp-@S20@}b--d; Table~@S8d@).'),
  ('(139 pairs; Figure~\\ref{fig:supp-S7}a)', '(139 pairs; Table~@S8b@)'),
  ('\\subsubsection{Functional tests of the profile classes}', '\\subsubsection{Ranked-set tests of profile amplitude}'),
  ('(7,780 tests per layer and measure across the three routes; Figure~\\ref{fig:supp-S9}a--c,g)', '(7,780 tests per layer and measure across the three routes; Figure~\\ref{fig:supp-@S7@}d; Table~@S8c@)'),
  ('Separately, ranked-set tests asked whether', 'Ranked-set tests asked whether')],
'within-species-recurrence.tex': [
  ('\\subsubsection{Matched references and metric membership}', '\\subsubsection{Matched references}'),
  (' P--S recurrence was not required (Figure~\\ref{fig:supp-S12}a,d; Table~S4e). Metric membership (Figure~\\ref{fig:supp-S12}b) compared the DWUI 0.02 and MPRO 0.10 recurrent sets of each species and contrast.', ' P--S recurrence was not required.'),
  ('(Figure~\\ref{fig:supp-S12}f--i; Table~S4f)', '(Figure~\\ref{fig:supp-@S8@}a--d; Table~@S3e@)'),
  ('(Figure~\\ref{fig:supp-S9}d--f)', '(Figure~\\ref{fig:supp-@S7@}a--c)'),
  ('(Figure~\\ref{fig:supp-S9}h)', '(Figure~\\ref{fig:supp-@S7@}e)')],
'candidate-intervals-and-predicted-sites.tex': [
  ('cohorts of Figures~1--4.', 'cohorts of Figures~@1@--@3@.'),
  ('The union pool used for Figure~5 added', 'The union pool used for Figure~@4@ added'),
  ('(Figure~5a; Figure~\\ref{fig:supp-S16}c,d; Table~S6b)', '(Figure~@4@a; Figure~\\ref{fig:supp-@S13@}c,d; Table~@S5b@)'),
  ('(Figure~5b; Figure~\\ref{fig:supp-S16}e,f; Table~S6c)', '(Figure~@4@f; Figure~\\ref{fig:supp-@S13@}e,f; Table~@S5c@)'),
  ('The transferred-motif rug of Figure~5f', 'The transferred-motif rug of Figure~@4@d'),
  ('(8 genes; Figure~5c)', '(8 genes; Figure~@4@c)'),
  ('The Figure~5f example, \\textit{Riok3}', 'The Figure~@4@d example, \\textit{Riok3}')],
'cross-species-comparisons.tex': [
  ('under three nulls (Figure~4c; Table~S5b)', 'under three nulls (Figure~\\ref{fig:supp-@S9@}a; Table~@S4b@)'),
  ('Per-study native $\\Delta$DWUI of the gene (Figure~4e)', 'Per-study native $\\Delta$DWUI of the gene (Figure~@3@e)'),
  ('Coefficients (Figure~\\ref{fig:supp-S15}g) have', 'Coefficients have'),
  ('Driver PASs were selected entirely in the discovery species', 'The switching PAS pair was selected entirely in the discovery species'),
  ('largest gain and largest loss were the drivers (ties', 'largest gain and largest loss formed the switching PAS pair (ties'),
  ('the number of matched drivers (0--2)', 'the number of matched switching PASs (0--2)'),
  ('the mean phastCons score of the drivers and their', 'the mean phastCons score of the switching PASs and their')],
'alternative-pas-calling-workflow.tex': [
  ("The conditional reference nulls, dual-reference ortholog tests and local/net-\\linebreak[4]displacement classification were evaluated separately on the two PAS workflows' eligible backgrounds (Tables~S9a,c,e). Cycle-score tests were shared because the GEX input was fixed, and adjustment feasibility was checked against each workflow's metadata with the formal primary candidate sets fixed (Table~S9b). Failed common support did not yield adjusted APA estimates. Profile-class entry into recurrence used polyApipe only (Table~S9d).",
   "The local/net-displacement classification was evaluated separately on the two PAS workflows' eligible backgrounds (Table~@S6a@). Cycle-score tests were shared because the GEX input was fixed, and adjustment feasibility was checked against each workflow's metadata with the formal primary candidate sets fixed (Table~@S7g@). Failed common support did not yield adjusted APA estimates.")],
'display-conventions-and-reproducibility.tex': [
  (' The conditional reference and dual-reference overlap tests used 20,000 draws. The dual-reference direction distribution was exact.', ''),
  ('The conditional reference nulls, cycle-support diagnosis, profile-to-endpoint comparison, local/net-displacement classification and example decompositions, and dual-reference cross-species tests were also designed after inspecting the earlier results. These analyses were not preregistered.',
   'The cycle-support diagnosis, the local/net-displacement classification and the example decompositions were also designed after inspecting the earlier results, as were the analyses marked $\\dagger$ in the Supplementary Tables (main-text Methods). No analysis was preregistered.'),
  ('Supplementary Tables~S1--S9 are supplied', 'Supplementary Tables~@S1@--@S10@ are supplied'),
  ('The source-data archive dated 9 October 2026 supplies', 'The source-data archive dated 11 October 2026 supplies'),
  ('The companion \\texttt{Supplementary\\_Source\\_Code\\_S9.zip} supplies the frozen source tables, protocols, original implementations and corrected-export code for Table~S9.', 'The companion \\texttt{Supplementary\\_Source\\_Code\\_v11b.zip} supplies the current implementations and pinned generator inputs, with \\texttt{Supplementary\\_Source\\_Code\\_S9.zip} retained as its required historical companion.'),
  (' Figure~3c counts the original P--M members, the genes measurable with common contributors and the P--S/P--M aligned genes, in that order, and its human analysis is confined to Xing2020.', '')],
}
DROP = {'state-resolved-profiles.tex': ['Route hold-out', 'Organization on expression graphs', 'APA--expression profile agreement'],
        'within-species-recurrence.tex': ['Conditional reference nulls for dual-reference alignment', 'Profile-class entry into endpoint recurrence', 'Intermediate states'],
        'cross-species-comparisons.tex': ['Dual-reference ortholog overlap and direction']}
PARA_DROP = {'state-resolved-profiles.tex': [
  (r'\nRobustness \(Figure~\\ref\{fig:supp-S6\}b\) refitted .*?adjusted Rand indices \(ARI\)\.\n', '\n'),
  (r' A sensitivity raised the minimum \$L\$ .*?\(Figure~\\ref\{fig:supp-S7\}e\)\.', ''),
  (r'To ask how much path beyond the endpoint displacement arises .*?in each draw and refitted \$k=3\$\. ', '')]}
SENT_DROP = {'state-resolved-profiles.tex': [r' Class over-representation compared each class with this background .*?and across classes\.']}
INSERT = {'dwui-metacells-and-state-comparisons.tex': [('Reference designs and estimators', 's1')],
          'cross-species-comparisons.tex': [('Shared proximal candidates', 's2'), ('Cis-feature models', 's6')],
          'within-species-recurrence.tex': [('Cycle scores and common support', 's4'), ('Cycle scores and common support', 's3')],
          'candidate-intervals-and-predicted-sites.tex': [('\\textit{Trip12} external evidence', 's5')]}

TABREF = re.compile(r'(Tables?~)(S\d+[a-z]?(?:(?:, | and |--|,)S?\d*[a-z]?)*)')
def _key(t):
    m = re.match(r'S(\d+)([a-z]?)$', t); return (int(m.group(1)), m.group(2))
def retable(t, f):
    def one(m):
        body = m.group(2)
        toks = re.findall(r'S\d+[a-z]?', body)
        conv = []
        for x in toks:
            v = spec.MAP.get(x) or spec.OLD_TABLE.get(x)
            assert v, ('unmapped or dropped table', f, x, m.group(0)); conv.append(v)
        if '--' in body:
            assert len(conv) == 2 and _key(conv[0]) < _key(conv[1]), ('range', f, m.group(0)); return m.group(1) + '@' + conv[0] + '@--@' + conv[1] + '@'
        conv = sorted(conv, key=_key)
        tail = ['@%s@' % c for c in conv]
        return m.group(1) + (tail[0] if len(tail) == 1 else ', '.join(tail[:-1]) + ' and ' + tail[-1])
    return TABREF.sub(one, t)
def refig(t, f):
    def one(m):
        o = int(m.group(1)); assert o in FIG, ('unmapped figure', f, o); return 'fig:supp-@S%d@}' % FIG[o]
    return re.sub(r'fig:supp-S(\d+)\}', one, t)

DST.mkdir(parents=True, exist_ok=True)
for src in sorted(SRC.glob('*.tex')):
    f = src.name; t = src.read_text()
    for title in DROP.get(f, []): t = drop_sub(t, title)
    for pat, r in PARA_DROP.get(f, []):
        t, n = re.subn(pat, r, t, flags=re.S); assert n == 1, (f, pat[:50], n)
    for pat in SENT_DROP.get(f, []):
        t, n = re.subn(pat, '', t, flags=re.S); assert n == 1, (f, pat[:50], n)
    for a, b in E.get(f, []): t = rep(t, a, b, f)
    # main-figure refs must all have been handled explicitly
    left = re.findall(r'Figures?~[1-5][a-z]?', t)
    assert not left, (f, left)
    t = refig(t, f); t = retable(t, f)
    for title, key in INSERT.get(f, []): t = insert_after_sub(t, title, blocks[key])
    t = t.replace('@', '')
    (DST / f).write_text(t)
print('written', len(list(DST.glob('*.tex'))), 'files to', DST)
