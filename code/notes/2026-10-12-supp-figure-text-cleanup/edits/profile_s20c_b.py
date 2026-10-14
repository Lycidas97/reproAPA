TARGET = 'notes/2026-10-10-supplement-v11/figures/profile.py'
EDITS = [
 ("            (356.0, 187.0, 453.0, 198.3)]       # 2026-10-12 S20c: '|ΔDWUI| / L' and 'Median D/L' ran together (2.5 pt gap); redrawn below\n",
  "            (356.0, 187.0, 453.0, 198.3),       # 2026-10-12 S20c: '|ΔDWUI| / L' and 'Median D/L' ran together (2.5 pt gap); redrawn below\n"
  "            (412.0, 274.0, 417.8, 283.5)]       # 2026-10-12 S20c: dot-axis '0' touched the bar-axis '1' (2.1 pt gap); '.5' drawn at the dotted line instead\n"),
 ("F.text(PAGE_W-4.5, S20C_HDR_TOP, 'Median D/L', 8.0, ha='right')",
  "F.text(PAGE_W-4.5, S20C_HDR_TOP, 'Median D/L', 8.0, ha='right')\n"
  "F.text((414.9+444.4)/2, 274.7, '.5', 8.0, ha='center')     # dot column 0–1 spans 414.9–444.4 pt; the dotted line sits at .5"),
]
