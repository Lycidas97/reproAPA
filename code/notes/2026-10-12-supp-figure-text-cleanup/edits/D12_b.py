TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D12.py'
EDITS = [
 ("F.xaxis(gx, yy-7.0, [-1, 0, 1, 2, 3], title=_gex_label)",
  "gb = F.xaxis(gx, yy-7.0, [-1, 0, 1, 2, 3])\nF.text(min((gx.p0+gx.p1)/2, x+w-0.3-textwidth(_gex_label, 8.0)/2), gb+1.6, _gex_label, 8.0, ha='center')   # axis title kept inside the panel"),
]
