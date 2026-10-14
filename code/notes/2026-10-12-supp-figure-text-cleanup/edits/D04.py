TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D04.py'
EDITS = [
 ("H1, H2, H3 = 108.0, 99.0, 112.0", "H1, H2, H3 = 115.0, 99.0, 101.5   # 2026-10-12: titles -> axis titles (a), title rows removed (c, d)"),
 ("""F.text((c1.p0+c1.p1)/2, y+3.0, 'Primary PAS (thousands)', 7.5, ha='center', color=MUTED)
F.text((c2.p0+c2.p1)/2, y+3.0, 'PAS within window', 7.5, ha='center', color=MUTED)
F.text((c3.p0+c3.p1)/2, y+3.0, 'UMI within window', 7.5, ha='center', color=MUTED)
""", ""),
 ("""F.xaxis(c1, ay, [0, 40, 80, 120])
F.xaxis(c2, ay, [0.6, 0.7, 0.8, 0.9, 1.0], ['.6', '.7', '.8', '.9', '1'])
F.xaxis(c3, ay, [0.88, 0.92, 0.96, 1.0], ['.88', '.92', '.96', '1'])""",
  """F.xaxis(c1, ay, [0, 40, 80, 120], title='Primary PAS (thousands)')
F.xaxis(c2, ay, [0.6, 0.7, 0.8, 0.9, 1.0], ['.6', '.7', '.8', '.9', '1'], title='Fraction of PAS within window')
F.xaxis(c3, ay, [0.88, 0.92, 0.96, 1.0], ['.88', '.92', '.96', '1'], title='Fraction of UMIs within window')"""),
 ("PY0, PY1, AXY = y+22.0, y+80.0, y+82.0", "PY0, PY1, AXY = y+13.0, y+71.0, y+73.0"),
 ("F.text((CX[0].p0+CX[1].p1)/2, y+3.0, 'Read ends around primary PAS', 7.5, ha='center', color=MUTED)\n", ""),
 ("    F.text((cx.p0+cx.p1)/2, y+12.0, lab, 7.5, ha='center')", "    F.text((cx.p0+cx.p1)/2, y+3.0, lab, 7.5, ha='center')"),
 ("F.vtext(x+4.4, (PY0+PY1)/2+4.2, 'Mean fraction per site', 8.0)", "F.vtext(x+4.4, (PY0+PY1)/2+8.0, 'Mean fraction per site', 8.0)"),
 ("F.text((CX[0].p0+CX[1].p1)/2, LEG3, 'Line colors: studies as in a', 7.5, ha='center', color=MUTED)\n", "# 2026-10-12: 'Line colors: studies as in a' moved to the caption\n"),
 ("F.text((DX[0].p0+DX[1].p1)/2, y+3.0, 'Genes with complete S–R–E–P profiles', 7.5, ha='center', color=MUTED)\n", ""),
 ("    F.text((dx.p0+dx.p1)/2, y+12.0, lab, 7.5, ha='center')", "    F.text((dx.p0+dx.p1)/2, y+3.0, lab, 7.5, ha='center')"),
 ("'PAS UMI per state (depth)'", "'Minimum PAS UMIs per metacell'"),   # production label since the 2026-10-07 PDF relabel
]
# follow-up (sed): vtext offset 8.0 -> 12.6 (letter clear zone)
