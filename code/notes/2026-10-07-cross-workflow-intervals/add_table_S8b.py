"""Add Table S8b (polyApipe strict interval pairs under the alternative calls) to the supplementary workbook plan."""
import csv, hashlib
from pathlib import Path
import pandas as pd
T = Path('notes/2026-10-06-supplement-tables'); D = Path('notes/2026-10-07-cross-workflow-intervals')
FILES = [D/'interval-pair-classes.tsv', D/'polyapipe-strict-interval-pairs-under-alternative-calls.tsv']
rows = list(csv.DictReader(open(T/'table-sources.tsv'), delimiter='\t')); cols = list(rows[0].keys())
rows = [r for r in rows if r['sheet'] != 'S8b']
nrows = 0
for f in FILES:
    b = f.read_bytes(); d = pd.read_csv(f, sep='\t', dtype=str, keep_default_na=False, na_filter=False); nrows += len(d)
    rows.append(dict(table='S8', sheet='S8b', sheet_title='Strict interval pairs under alternative PAS calls', layer='sheet', old_unit='',
                     role='crossworkflow_intervals', via='R5P7;SR;SM', path=str(f), bytes=str(len(b)), rows=str(len(d)), filter='',
                     sha256=hashlib.sha256(b).hexdigest(), columns=','.join(d.columns)))
with open(T/'table-sources.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
plan = list(csv.DictReader(open(T/'table-plan.tsv'), delimiter='\t')); pc = list(plan[0].keys())
plan = [p for p in plan if p['sheet'] != 'S8b']
i = max(k for k, p in enumerate(plan) if p['sheet'] == 'S8a')
new = {c: '' for c in pc}
new.update(table='S8', table_title=plan[i]['table_title'], sheet='S8b', sheet_title='Strict interval pairs under alternative PAS calls',
           status='new', cited_by='R5P7; Discussion; SR; SM', supp_figures='S21', sheet_files=str(len(FILES)), sheet_rows=str(nrows),
           source_data_files='0', note='2026-10-07: each polyApipe strict recurrent interval pair located in the APADETECT p50 coordinate regions (region containing the 3\'-most base) and classed reproduced / below threshold / direction conflict / not measurable.')
plan.insert(i+1, new)
with open(T/'table-plan.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, pc, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(plan)
p = Path('supplement/supplementary-tables/main.tex'); t = p.read_text()
a = r"\item[] \textbf{S8a}\quad Alternative PAS calling (Figure~S21)."
b = a + "\n" + r"\item[] \textbf{S8b}\quad Strict interval pairs under alternative PAS calls."
if 'S8b' not in t:
    assert t.count(a) == 1; p.write_text(t.replace(a, b))
print('ok', nrows)
