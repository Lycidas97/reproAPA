TARGET = 'notes/2026-10-10-supplement-v11/figures/profile.py'
EDITS = [
 ("F2_CLEAR = [(121.0, 188.8, 215.0, 218.6)]       # 'L = |R−S| + |E−R| + |P−E|', 'D = |P−S|', 'T = 1 − D/L'",
  "F2_CLEAR = [(121.0, 188.8, 215.0, 218.6),       # 'L = |R−S| + |E−R| + |P−E|', 'D = |P−S|', 'T = 1 − D/L'\n"
  "            (356.0, 187.0, 453.0, 198.3)]       # 2026-10-12 S20c: '|ΔDWUI| / L' and 'Median D/L' ran together (2.5 pt gap); redrawn below\n"
  "S20C_HDR_TOP = float(os.environ.get('S20C_HDR_TOP', '189.1'))"),
 ("F.t.text_top(4.5+0.2, R7[1]+0.5+0.4+DY7, 'd', 11.0, font='Bol')",
  "F.t.text_top(4.5+0.2, R7[1]+0.5+0.4+DY7, 'd', 11.0, font='Bol')\n"
  "# S20c column headers (2026-10-12): bar block 358.2–408.3 pt, dot column 414.9–444.4 pt in the former Figure 2 coordinates.\n"
  "# The bar header loses its spaces and is centred on the bars; the dot-column header is right-aligned to the 4.5 pt margin.\n"
  "F.text((358.2+408.3)/2, S20C_HDR_TOP, '|ΔDWUI|/L', 8.0, ha='center')\n"
  "F.text(PAGE_W-4.5, S20C_HDR_TOP, 'Median D/L', 8.0, ha='right')"),
]
