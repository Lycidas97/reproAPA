TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D26.py'
EDITS = [
 ("H1, H2, H3, H4 = 118.0, 100.0, 112.0, 100.0", "H1, H2, H3, H4 = 122.0, 100.0, 120.0, 100.0   # 2026-10-12: panel titles and notes moved to axis titles or the caption"),
 ("    F.text(x0+gw/2, y+12.0, f'n = {len(z)}', 6.8, ha='center', color=MUTED)\n", ""),
 ("F.text(x+24.0+2*gw, y+h-7.4, '1 = 842 nt, 2 = 706 nt; bold = mean', 7.0, ha='center', color=MUTED)\n", ""),
 ("    F.text((sx.p0+sx.p1)/2, top, lab, 7.5, ha='center')\n", ""),
 ("    F.xaxis(sx, top+36.0, ticks, tl)\nF.text(x+w/2, y+h-7.4, f\"pairs: {Bv['CIP|pairs']} | {Bv['WOT|pairs']}\", 6.8, ha='center', color=MUTED)",
  "    F.xaxis(sx, top+36.0, ticks, tl, title=lab)"),
 ("F.text((cx_.p0+cx_.p1)/2, y+3.0, 'Interval ΔI, P − M (pp)', 7.5, ha='center')\n", ""),
 ("F.xaxis(cx_, yy-4.0, [-30, -20, -10, 0], ['−30', '−20', '−10', '0'])", "F.xaxis(cx_, yy-4.0, [-30, -20, -10, 0], ['−30', '−20', '−10', '0'], title='Interval ΔI, P − M (pp)')"),
 ("'Distance from the terminal-exon start, chr1 (−) (kb)', 8.0, ha='center')\nF.text(x+0.4, y+h-7.4, 'Same complete intron chain; 0 samples detect both models', 6.8, color=MUTED)",
  "'Distance from the terminal-exon start (kb)', 8.0, ha='center')"),
 ("cx0 = x+58.0; cw = (x+w-1.0-cx0)/3; rh = 11.8; e0 = y+24.0", "cx0 = x+58.0; cw = (x+w-1.0-cx0)/3; rh = 14.0; e0 = y+14.0"),
 ("    F.text(cx0+cw*(j+0.5), y+13.0, lab, 7.0, ha='center')\nF.text(cx0+cw*1.5, y+3.0, 'AGO2 CLIP overlap', 7.5, ha='center')", "    F.text(cx0+cw*(j+0.5), y+3.0, lab, 7.0, ha='center')"),
 ("F.text(x+w/2+20.0, y+h-7.4, 'WT1 A/B = one biological unit', 6.8, ha='center', color=MUTED)\n", ""),
 ("F.text((gx_.p0+gx_.p1)/2, y+3.0, 'ΔDWUI × 1,000', 7.5, ha='center')\n", ""),
 ("F.xaxis(gx_, y+62.0, [-8, -4, 0, 4], ['−8', '−4', '0', '4'])\nfor k, (col, cc, lab) in enumerate(COMP):\n    F.key_row(x+6.0, y+78.0+k*11.4, [('block', cc, lab)], size=7.0)\nF.key_row(x+82.0, y+78.0,",
  "F.xaxis(gx_, y+62.0, [-8, -4, 0, 4], ['−8', '−4', '0', '4'], title='ΔDWUI × 1,000')\nfor k, (col, cc, lab) in enumerate(COMP):\n    F.key_row(x+6.0, y+86.0+k*11.4, [('block', cc, lab)], size=7.0)\nF.key_row(x+82.0, y+86.0,"),
 ("F.text(x+0.4, y+h-7.4, 'Upper three tiers = Figure 5f rug; all are predicted sites, not binding', 6.8, color=MUTED)\n", "# 2026-10-12: tier note (Figure 4d rug; predicted, not binding) moved to the caption\n"),
]
# follow-up (sed): H4 100 -> 91 (h note removed)
