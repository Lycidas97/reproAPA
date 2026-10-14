"""Figure 1 (2026-10-12 main-figure text cleanup).
1d lower, documented-source P–S axis: the tick labels '−0.04' and '−0.02' stood 2.3 pt apart. '−0.04' is dropped
(its tick mark stays); '−0.02', '0', '+0.02' keep the data range (−0.02 to +0.02) labelled."""
import json
from fixlib import Job
J = Job(1)
J.drop_word('−0.04', (300, 558, 330, 571), note='1d lower x tick label (crowded with −0.02)')
rec = J.run()
print(json.dumps({k: rec[k] for k in ('removed_kinds', 'words_removed', 'words_added', 'sha256')}, indent=1, ensure_ascii=False))
json.dump(rec, open('out/receipt.fig1.json', 'w'), indent=1, ensure_ascii=False)
