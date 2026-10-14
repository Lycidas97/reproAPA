TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D10.py'
EDITS = [
 ("F.text((lx.p0+lx.p1)/2+6, y+3.0, f'{len(p)} genes, both q < .05: {int(both.sum())}', 7.5, ha='center')\n",
  "# 2026-10-12: gene counts (len(p), both.sum()) moved to the caption; they remain in values['a']\n"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, 'Joint q < .05, fraction', 7.5, ha='center')\n", ""),
 ("    yy += 19.8\nF.xaxis(bx, yy-8.0, [0, 0.5, 1], ['0', '.5', '1'])", "    yy += 18.4\nF.xaxis(bx, yy-8.0, [0, 0.5, 1], ['0', '.5', '1'], title='Fraction with joint q < .05')"),
 ("F.text((tx.p0+tx.p1)/2, y+3.0, 'TVD, whole library', 7.5, ha='center')\n", ""),
 ("F.xaxis(tx, yy-10.0, [0, 0.01, 0.02], ['0', '.01', '.02'])\nF.key_row(None, y+h-17.4, [('line', (NULLBAR, None, 'bold'), 'Source range')], gap=6.0, center=(tx.p0+tx.p1)/2)\nF.text(x+w/2, y+h-7.4, 'Matched / measurable genes', 7.0, ha='center', color=MUTED)",
  "F.xaxis(tx, yy-10.0, [0, 0.01, 0.02], ['0', '.01', '.02'], title='TVD')\nF.key_row(None, y+h-7.4, [('line', (NULLBAR, None, 'bold'), 'Source range')], gap=6.0, center=(tx.p0+tx.p1)/2)\n# 2026-10-12: the 'Matched / measurable genes' note moved to the caption"),
]
