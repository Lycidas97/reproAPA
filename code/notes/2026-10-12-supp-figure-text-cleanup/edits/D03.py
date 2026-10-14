TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D03.py'
EDITS = [
 ("    F.text(x+w-1.0, y+1.0, f'{len(q):,} cells shown', 7.0, ha='right', color=MUTED)\n",
  "    # 2026-10-12: the displayed-cell count moved to the caption\n"),
 ("HE = 88.0", "HE = 97.0"),
 ("F.text((sx.p0+sx.p1)/2, y+3.0, 'Native program score', 7.5, ha='center')\n", ""),
 ("F.xaxis(sx, yy-4.0, [0, 1, 2, 3])", "F.xaxis(sx, yy-4.0, [0, 1, 2, 3], title='Native program score')   # 2026-10-12: former panel title"),
 ("F.text((fx.p0+fx.p1)/2, y+3.0, 'Held-out correlation excess', 7.5, ha='center')\nFv = {}; rp = 9.0; y0 = y+21.0",
  "Fv = {}; rp = 8.0; y0 = y+19.0   # 2026-10-12: panel title became the x-axis title"),
 ("F.xaxis(fx, y0+11*rp+6.0, [0, 0.02, 0.04], ['0', '.02', '.04'])", "F.xaxis(fx, y0+11*rp+6.0, [0, 0.02, 0.04], ['0', '.02', '.04'], title='Held-out correlation excess')"),
]
# follow-up: y0 19.0 -> 19.6 (letter clear zone), applied with sed
