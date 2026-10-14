TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D27.py'
EDITS = [
 ("H1, H2 = 153.0, 150.0", "H1, H2 = 153.0, 160.0   # 2026-10-12: panel titles moved to axis titles or the caption"),
 ("F.text((sx.p0+sx.p1)/2, y+3.0, 'PAS-axis span (kb), 51 genes', 7.5, ha='center')\n", ""),
 ("F.yaxis(sy, sx.p0, T, TL, 'Schiebinger2019–2i')\nF.xaxis(sx, sy.p0+1.0, T, TL, 'Guo2019–CIP')",
  "F.yaxis(sy, sx.p0, T, TL, 'PAS-axis span, Schiebinger2019–2i (kb)')\nF.xaxis(sx, sy.p0+1.0, T, TL, 'PAS-axis span, Guo2019–CIP (kb)')"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, f'Terminal-family shift, P − M (nt), {len(ts)} genes', 7.5, ha='center')\n", ""),
 ("F.yaxis(by, bx.p0, T, TL, 'Schiebinger2019–2i')\nF.xaxis(bx, by.p0+1.0, T, TL, 'Guo2019–CIP')",
  "F.yaxis(by, bx.p0, T, TL, 'Shift, Schiebinger2019–2i (nt)')\nF.xaxis(bx, by.p0+1.0, T, TL, 'Terminal-family shift, Guo2019–CIP (nt)')"),
 ("F.text((cx1.p0+cx1.p1)/2, y+3.0, 'Terminal-family shift (nt)', 7.5, ha='center')\nF.text((cx2.p0+cx2.p1)/2, y+3.0, 'Contribution to native ΔDWUI × 1,000', 7.5, ha='center')\n", ""),
 ("F.xaxis(cx1, yend+1.0, [-600, -300, 0], ['−600', '−300', '0'])\nF.xaxis(cx2, yend+1.0, [-30, 0, 30], ['−30', '0', '30'])",
  "F.xaxis(cx1, yend+1.0, [-600, -300, 0], ['−600', '−300', '0'], title='Terminal-family shift (nt)')\nF.xaxis(cx2, yend+1.0, [-30, 0, 30], ['−30', '0', '30'], title='Contribution to native ΔDWUI × 1,000')"),
]
