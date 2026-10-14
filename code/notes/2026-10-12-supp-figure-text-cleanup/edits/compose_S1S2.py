TARGET = 'notes/2026-10-06-supplement-polyapipe/tools/compose.py'
EDITS = [
 ("sub = (HC-HDRC-2*5.0-13.0)/3", "sub = (HC-HDRC-2*5.0-3.0)/3"),
 ("    F.text(X0+0.4, (sy.p0+sy.p1)/2-3.4, prog, 7.0)\n", "    F.text(X0+9.4, (sy.p0+sy.p1)/2-3.4, prog, 7.0)\n"),
 ("unit_rules(yf, y+HC-13.0)\nF.text((CX0+RR)/2, y+HC-7.4, 'Mean log1p(CP10k) of each qualified library × state (≥ 20 cells); color = display family', 7.0, ha='center', color=MUTED)\n",
  "unit_rules(yf, y+HC-3.0)\n# 2026-10-12: the explanatory line under the panel moved to the caption; the quantity is a y-axis title\n"
  "F.vtext(X0+3.8, y+HDRC+(HC-HDRC-3.0)/2, 'log1p(CP10k)', 8.0)\n"),
 ("MARKERS_H = ('HB, HC = 62.0, 96.0', 'HB, HC = 62.0, 131.0')", "MARKERS_H = ('HB, HC = 62.0, 96.0', 'HB, HC = 62.0, 121.0')"),
 ("('H1, H2 = 136.0, 111.5', 'H1, H2 = 136.0, 97.5'), ('top, bot = y+27.0, y+80.0', 'top, bot = y+27.0, y+66.0')",
  "('H1, H2 = 136.0, 111.5', 'H1, H2 = 136.0, 88.5'), ('top, bot = y+18.0, y+71.0', 'top, bot = y+18.0, y+57.0')"),
]
