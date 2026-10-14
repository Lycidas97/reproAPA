TARGET = 'notes/2026-10-06-supplement-polyapipe/tools/compose.py'
EDITS = [
 ("    F.text((px.p0+px.p1)/2, y+3.0, f'{db}, G2 tests (L3 full null)', 7.5, ha='center')\n",
  "    F.text((px.p0+px.p1)/2, y+3.0, db, 7.5, ha='center')   # 2026-10-12: facet label only; G2 tests / L3 full null in the caption\n"),
]
