#!/usr/bin/env python3
"""v11 table catalogue for the Supplement (candidate; supplement/ is never written).
Old worksheets keep the item text of supplement/supplementary-tables/main.tex with the figure parenthetical replaced by
spec_v11.REL (v11 figure numbers); new worksheets use spec_v11.NEW titles; a dagger marks POSTHOC worksheets."""
import re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = next(q for q in [HERE, *HERE.parents] if (q / 'notes').is_dir())
sys.path.insert(0, str(HERE)); import spec_v11 as spec
PROD = (ROOT / 'revision-backups/20261011-before-v11-apply/supplement/supplementary-tables/main.tex').read_text()
OUT = HERE / 'supplementary-tables.v11.tex'
INV = {v: k for k, v in spec.MAP.items()}
POSTHOC = set(spec.NEW) | set(spec.POSTHOC_OLD)
EXTRA = {'S10a': 'Marker-defined P and M references.'}   # old S9f postdates the production catalogue

def tex(s):
    for a, b in [('\\', '\\textbackslash{}'), ('&', '\\&'), ('%', '\\%'), ('_', '\\_'), ('#', '\\#'), ('–', '--'), ('—', '---'),
                 ('′', "$'$"), ('−', '$-$')]:
        s = s.replace(a, b)
    return s
def figtex(rel):
    parts = [p.strip().replace('–', '--') for p in rel.split(';') if p.strip()]
    main = ['Figure~' + p.removeprefix('Figure ').strip() for p in parts if p.startswith('Figure')]
    supp = [p for p in parts if not p.startswith('Figure')]
    out = main[:]
    if supp: out.append(('Figures~' if len(supp) > 1 or '--S' in supp[0] else 'Figure~') + '; '.join(supp))
    return '; '.join(out)

old_item = {m.group(1): m.group(2).strip() for m in re.finditer(r'^\\item\[\] \\textbf\{(S\d+[a-z])\}\\quad (.*)$', PROD, re.M)}
assert set(spec.MAP) - {'S9f'} <= set(old_item), set(spec.MAP) - set(old_item)
def strip_fig(t):
    m = re.search(r'\s?\(([^()]*)\)\.$', t)
    if not m or 'Figure' not in m.group(1): return t.rstrip('.')
    keep = m.group(1)[:m.group(1).index('Figure')].rstrip('; ')
    return t[:m.start()] + (f' ({keep})' if keep else '')
lines = [r'\section{Supplementary tables}\label{sec:supp-tables}', '',
         r'Supplementary Tables~S1--S10 are supplied as one workbook, \texttt{Supplementary\_Tables\_\allowbreak S1-S10.xlsx}, with one worksheet per lettered table (for example, Table~S4b), grouped by topic. Each worksheet holds one or more data blocks with complete rows, denominators, negative results and correction families. Worksheet headers name the related figures or, where none applies, the analysis. Worksheets marked \textsuperscript{\dag} hold analyses specified after the initial results. The Contents sheet indexes the worksheets, and the Sources sheet gives the source file, row filter and SHA-256 hash of every data block. In these files, columns labeled WPP denote DWUI. Every contrast code denotes target minus reference: P\_S, PS and P-S denote P$-$S, and the reference-first codes of the interval tables denote, for example, SP = P$-$S and MP = P$-$M. The Contents sheet lists all codes. Plotting exports that are not needed to reuse the results are provided in the source-data archive.', '']
sheets = [s for s in spec.ORDER if s not in ('Contents', 'Sources')]
for t in dict.fromkeys(spec.table_of(s) for s in sheets):
    lines += [rf'\paragraph{{Table~{t}. {tex(spec.GROUPS[t])}.}}', r'\begin{itemize}[leftmargin=*,itemsep=1pt,topsep=2pt]']
    for s in [x for x in sheets if spec.table_of(x) == t]:
        if s in spec.NEW: body = tex(spec.NEW[s]['title'])
        elif s in EXTRA: body = EXTRA[s].rstrip('.')
        elif s in spec.TITLE_OVERRIDE: body = tex(spec.TITLE_OVERRIDE[s])
        else: body = strip_fig(old_item[INV[s]])
        rel = spec.REL.get(s, '')
        if rel:
            body = body[:-1] + f'; {figtex(rel)})' if body.endswith(')') and 'Figure' not in body else body + f' ({figtex(rel)})'
        dag = r'\textsuperscript{\dag}' if s in POSTHOC else ''
        lines.append(rf'\item[] \textbf{{{s}}}{dag}\quad {body}.')
    lines += [r'\end{itemize}', '']
lines.append(r'The companion \texttt{Supplementary\_Source\_Code\_v11b.zip} supplies the current implementations and pinned generator inputs. Its manifest records file hashes and its required historical companion, \texttt{Supplementary\_Source\_Code\_S9.zip}.')
OUT.write_text('\n'.join(lines) + '\n')
txt = OUT.read_text()
assert all(txt.count(rf'\textbf{{{s}}}') == 1 for s in sheets)
assert txt.count('{') == txt.count('}'), 'unbalanced braces'
print('ok', OUT.name, len(sheets), 'items')
