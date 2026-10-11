"""Corrected proximal-overlap-tests.tsv for Supplementary Figure S21 (A01) and Table S8a (2026-10-07).
polyApipe rows: bundle-main overlap_null_results.tsv (A1 correction; fields absent from the server receipt -> NA).
Alternative-calling rows: notes/2026-10-07-cross-null-followup/ladder-results/apadetect-p50/original-nulls/overlap_null_results.tsv.
Schema and row order follow notes/2026-10-05-dual-main-concordance/proximal-overlap-tests.tsv."""
import csv, hashlib, json
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
old_p = ROOT/'notes/2026-10-05-dual-main-concordance/proximal-overlap-tests.tsv'
old = list(csv.DictReader(open(old_p), delimiter='\t')); cols = list(old[0].keys())
src = {'polyapipe': HERE/'bundle-main/tables/cross_overlap/overlap_null_results.tsv',
       'alt': ROOT/'notes/2026-10-07-cross-null-followup/ladder-results/apadetect-p50/original-nulls/overlap_null_results.tsv'}
labels = list(dict.fromkeys(r['source'] for r in old)); print('sources', labels)
tables = {labels[0]: list(csv.DictReader(open(src['polyapipe']), delimiter='\t')),
          labels[1]: list(csv.DictReader(open(src['alt']), delimiter='\t'))}
out = []
for o in old:
    n = [r for r in tables[o['source']] if r['contrast'] == o['contrast'] and r['model'] == o['model']]
    assert len(n) == 1 and int(n[0]['observed']) == int(o['observed']), o
    row = {c: (n[0].get(c, '') if c != 'source' else o['source']) for c in cols}
    out.append({c: ('NA' if v == '' else v) for c, v in row.items()})
p = HERE/'alt-s21/proximal-overlap-tests.tsv'
with open(p, 'w', newline='') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(out)
m = json.loads((ROOT/'notes/2026-10-06-alt-catalogue-suppfig/sources.json').read_text())
m['files']['proximal'] = dict(path=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
m['note'] += ' | 2026-10-07: proximal replaced by the A1 eligibility-corrected nulls (notes/2026-10-07-a1-figure-revision/make_s21_inputs.py).'
(HERE/'alt-s21/sources-a1.json').write_text(json.dumps(m, indent=1, ensure_ascii=False)+'\n')
for r in out: print(r['source'], r['contrast'], r['model'], r['observed'], r['null_mean'], r['null_q025'], r['null_q975'], r['q_within_model_two_contrasts'])
