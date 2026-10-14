TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D25.py'
EDITS = [
 ("F.text(x0b+cwb*len(kept)/2, y+3.0, f'Kept: expressed ≥ 2 CPM in P of both studies ({len(kept)} RBPs, {t[t.keep][[\"start0\", \"end0\"]].drop_duplicates().shape[0]} positions)', 7.5, ha='center')",
  "F.text(x0b+cwb*len(kept)/2, y+3.0, 'Kept', 7.5, ha='center')   # 2026-10-12: keep rule and counts moved to the caption"),
]
