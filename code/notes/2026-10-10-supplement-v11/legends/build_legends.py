"""Build the v11 supplementary-figure legends (candidate; supplement/ is never written).
Source: supplement/supplementary-figures/main.tex (production legends) and ../captions/new-figures.tex (new figures).
Output: supplementary-figures.v11.tex. Every substitution must match exactly once in its block (asserted)."""
import re
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
src = (ROOT/'revision-backups/20261011-before-v11-apply/supplement/supplementary-figures/main.tex').read_text()
newcap = (HERE.parent/'captions/new-figures.tex').read_text()

starts = [m.start() for m in re.finditer(r'\\begin\{figure\}\[H\]', src)]
head = src[:starts[0]]
blocks = {}
for k, s0 in enumerate(starts):
    s1 = starts[k+1] if k+1 < len(starts) else len(src)
    b = src[s0:s1]
    old = int(re.search(r'SuppFig\.S(\d+)\.pdf', b).group(1))
    blocks[old] = b
tail_m = re.search(r'\\clearpage\s*\Z', blocks[21])

def renum(b, old, new):
    b = b.replace(f'SuppFig.S{old}.pdf', f'SuppFig.@{new}@.pdf').replace(f'fig:supp-S{old}}}', f'fig:supp-@{new}@}}')
    return b

def sub(b, pairs, n):
    for a, c in pairs:
        assert b.count(a) == 1, (n, a[:60], b.count(a))
        b = b.replace(a, c)
    return b

COPY = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 10, 9: 13, 11: 14, 13: 16, 14: 19, 15: 20, 16: 17, 17: 18, 21: 21}
FIX = {
 6: [('drawn as in Figure~3a: upper bar', 'drawn as segmented bars: upper bar'),
     ('as kept/original bars as in Figure~3b.', 'as kept/original bars as in Figure~2b.')],
 9: [('At the primary setting (0.02, supported layer, Figure~4c) it has q = 0.051 within this grid (q = 0.032 over the two contrasts of Figure~4c)',
      'At the primary setting (0.02, supported layer) it has q = 0.051 within this grid (q = 0.032 over the two contrasts)'),
     ('the within-species L0--L3 of Figure~S10a', 'the within-species L0--L3 of Figure~S6a')],
 11: [('shown as tiles as in Figure~4d', 'shown as tiles as in Figure~3d'),
      ('Source data: Tables~S4e, S5c and S5d.', 'Source data: Tables~S5c and S5d.')],   # f,g use S5d (PDLIM7); S4e (dual-reference) not carried
 13: [('The first two steps are those of Figure~5c:', 'The first two steps are those of Figure~4c:'),
      ('A further step, not part of Figure~5c,', 'A further step, not part of Figure~4c,'),
      ('with gene-block bootstrap 95\\% intervals, as in Figure~5a.', 'with gene-block bootstrap 95\\% intervals, as in Figure~4a.'),
      ('as in Figure~5b (a cohort distinct from that of Figure~5a)', 'as in Figure~4f (a cohort distinct from that of Figure~4a)')],
 14: [('(the further screening step of Figure~S16b; Figure~5c shows only the first two steps)', '(the further screening step of Figure~S13b; Figure~4c shows only the first two steps)')],
 15: [('a low-net-shift gene of Figure~5c with', 'a low-net-shift gene of Figure~4c with'),
      ('follow the rug rule of Figure~5f.', 'follow the rug rule of Figure~4d.')],
 17: [("row of Figure~5f and is not", "row of Figure~4d and is not"),
      ("third rug row of Figure~5f,", "third rug row of Figure~4d,"),
      ("strict intervals (Figure~\\ref{fig:supp-S20})", "strict intervals (Figure~\\ref{fig:supp-@15@})")],
 21: [('(Figure~5c and its further step in Figure~S16b)', '(Figure~4c and its further step in Figure~S13b)')],
}
out = {}
for new, old in COPY.items():
    b = renum(blocks[old], old, new)
    b = sub(b, FIX.get(new, []), new)
    out[new] = b

# ---- S7 = old S9 d–h, S8 = old S12 f–i (relettered)
def split_block(old, new, title, keep, intro=''):
    b = blocks[old]
    cont = b[b.index('(continued).}'):]
    paras = re.findall(r'\\textbf\{([a-z])\}, (.*?)(?=\n\n)', cont, re.S)
    src_line = re.search(r'Source data: [^\n]*', cont).group(0)
    body = []
    for o, n in keep.items():
        t = dict(paras)[o]
        body.append(f'\\textbf{{{n}}}, {t}')
    return ('\\begin{figure}[H]\n\\begin{adjustwidth}{-\\extralength}{0cm}\n\\centering\n'
            f'\\includegraphics[width=150mm]{{SuppFig.@{new}@.pdf}}\n\\end{{adjustwidth}}\n'
            f'\\caption{{{title} Panel descriptions follow on the next page.\\label{{fig:supp-@{new}@}}}}\n\\end{{figure}}\n\\clearpage\n'
            '\\begingroup\\small\\setlength{\\parindent}{0pt}\\setlength{\\parskip}{5pt}\n'
            f'\\noindent\\textbf{{Figure~\\ref{{fig:supp-@{new}@}} (continued).}} {title}\n\n'
            + (intro+'\n\n' if intro else '') + '\n\n'.join(body) + f'\n\n{src_line}\n\n\\endgroup\n\\clearpage\n\n')

s7 = split_block(9, 7, 'Functional enrichment of the recurrence sets.', {'d': 'a', 'e': 'b', 'f': 'c', 'g': 'd', 'h': 'e'})
s7 = sub(s7, [('\\textbf{d}, Ranked APA effect of every tested term in the supported against the expanded layer (r = 0.89).',
               '\\textbf{d}, Ranked APA effect of every tested GO or KEGG term (the S--R--E--P range of APA in the term members against the other complete-profile genes) in the supported against the expanded layer (r = 0.89).')], 7)
out[7] = s7
s8 = split_block(12, 8, 'APA--GEX response identity and direction.', {'f': 'a', 'g': 'b', 'h': 'c', 'i': 'd'})
s8 = sub(s8, [('Source data: Tables~S4e, S4f and S5d.', 'Source data: Table~S4f.'),   # S4e and S5d served the dropped panels a-e
              ('as 2 $\\times$ 2 tiles as in Figure~3e:', 'as 2 $\\times$ 2 tiles as in Figure~2e:'),
              ('(clipped at $\\times$0.57 and $\\times$1.74, as in Figure~3e)', '(clipped at $\\times$0.57 and $\\times$1.74, as in Figure~2e)')], 8)
out[8] = s8

# ---- table numbers: production legends carry pre-v11 worksheet numbers; map them with spec_v11.MAP / OLD_TABLE.
# Runs before the new legends are added (those are written in v11 numbers). Lists are re-sorted; ranges must stay ordered.
import sys
sys.path.insert(0, str(HERE.parent / 'tables'))
import spec_v11 as spec
TABREF = re.compile(r'(Tables?~)(S\d+[a-z]?(?:(?:, | and |--)S\d+[a-z]?)*)')
def _key(t):
    m = re.match(r'S(\d+)([a-z]?)$', t); return (int(m.group(1)), m.group(2))
def retable(b, n):
    def one(m):
        conv = []
        for t in re.findall(r'S\d+[a-z]?', m.group(2)):
            v = spec.MAP.get(t) or spec.OLD_TABLE.get(t)
            assert v, ('unmapped table', n, t)
            conv.append(v)
        if '--' in m.group(2):
            assert len(conv) == 2 and _key(conv[0]) < _key(conv[1]), ('range', n, m.group(0), conv)
            return m.group(1) + '--'.join(conv)
        conv = sorted(conv, key=_key)
        return m.group(1) + (conv[0] if len(conv) == 1 else ', '.join(conv[:-1]) + ' and ' + conv[-1])
    return TABREF.sub(one, b)
TABLE_CHANGES = {}
for n in sorted(out):
    before = TABREF.findall(out[n])
    out[n] = retable(out[n], n)
    TABLE_CHANGES[n] = (before, TABREF.findall(out[n]))

# ---- new figures
TITLE = {10: 'Event classes of the shared P--M proximal response and corresponding PAS switches.',
         12: 'Specificity of the shared proximal response to naive pluripotency.',
         18: 'Cell-cycle and CPA-factor expression context of P, and the upstream UGUA test.',
         19: 'Explanations tested for the shared P--M proximal set, and cis-feature prediction of cross-species sharing.',
         20: 'State-resolved APA profile classes and partial endpoint retention.'}
parts = re.split(r'% -{10,} (S\d+)\n', newcap)
cap = {int(parts[i][1:]): parts[i+1].strip() for i in range(1, len(parts), 2)}
for n, t in TITLE.items():
    body = cap[n].replace(f'fig:supp-S{n}}}', f'fig:supp-@{n}@}}')
    assert body.startswith('\\noindent\\textbf{Figure~\\ref{fig:supp-@%d@} (continued).} %s' % (n, t)), (n, body[:120])
    out[n] = ('\\begin{figure}[H]\n\\begin{adjustwidth}{-\\extralength}{0cm}\n\\centering\n'
              f'\\includegraphics[width=150mm]{{SuppFig.@{n}@.pdf}}\n\\end{{adjustwidth}}\n'
              f'\\caption{{{t} Panel descriptions follow on the next page.\\label{{fig:supp-@{n}@}}}}\n\\end{{figure}}\n\\clearpage\n'
              '\\begingroup\\small\\setlength{\\parindent}{0pt}\\setlength{\\parskip}{5pt}\n' + body + '\n\n\\endgroup\n\\clearpage\n\n')
assert sorted(out) == list(range(1, 22)), sorted(out)
# ---- check against spec_v11.REL: every cited worksheet must list this figure (forward, asserted);
# worksheets whose related figures are not cited back are reported only.
def rel_figs(rel):
    f = set()
    for part in re.split(r';\s*', rel):
        part = part.strip()
        if part.startswith('Figure'): continue
        r = re.match(r'S(\d+)–S(\d+)$', part)
        f |= set(range(int(r.group(1)), int(r.group(2)) + 1)) if r else {int(re.match(r'S(\d+)', part).group(1))}
    return f
def expand(g):
    toks = re.findall(r'S\d+[a-z]?', g)
    if '--' in g:
        (a, la), (b, lb) = _key(toks[0]), _key(toks[1])
        assert a == b, ('cross-table range', g)
        return {f'S{a}{chr(c)}' for c in range(ord(la), ord(lb) + 1)}
    return set(toks)
cited = {n: set().union(*[expand(g) for _, g in TABREF.findall(out[n])]) for n in out}
bad = [(n, t, spec.REL.get(t)) for n in cited for t in cited[n] if n not in rel_figs(spec.REL.get(t, ''))]
assert not bad, ('cited worksheet not related to figure', bad)
for t, rel in spec.REL.items():
    for n in sorted(rel_figs(rel)):
        if t not in cited.get(n, set()): print('not-cited-back', t, '-> S%d' % n, '(REL %s)' % rel)
for n, (b, a) in TABLE_CHANGES.items():
    if b != a: print('table-renum S%d' % n, [x[0] + x[1] for x in b], '->', [x[0] + x[1] for x in a])
txt = head + ''.join(out[n] for n in range(1, 22))
txt = re.sub(r'@(\d+)@', r'S\1', txt)
(HERE/'supplementary-figures.v11.tex').write_text(txt)
# report remaining main-figure references for review
for m in re.finditer(r'Figures?~[1-5][a-z]?', txt):
    print('main-ref', txt[max(0, m.start()-60):m.end()+20].replace('\n', ' '))
print('labels', len(re.findall(r'\\label\{fig:supp-S\d+\}', txt)))
