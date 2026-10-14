TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D19.py'
EDITS = [
 ("WA = 300.0\nx, y, w, h = F.panel('a', X0, Y1, WA, H1)\nF.text(x+74.0, y+3.0, 'Schematic gene', 7.5, ha='center')",
  "WA = 150.0   # 2026-10-12: feature definitions and worked-example values moved to the caption; schematic only\nx, y, w, h = F.panel('a', X0, Y1, WA, H1)"),
 ("""F.text(x+0.4, y+68.0, 'Driver = 1', 7.0, color=MUTED)
F.text(x+0.4, y+77.0, 'PAS fraction = 2 / 5', 7.0, color=MUTED)
F.text(x+0.4, y+86.0, 'A_g = 4 / 9 = .44', 7.0, color=MUTED)
F.key_row(x+62.0, y+68.0, [('line', (INK, None, 'bold'), 'matched PAS')])
F.key_row(x+62.0, y+77.0, [('line', (AXC, None, 'bold'), 'unmatched PAS')])
F.key_row(x+62.0, y+86.0, [('block', WIN, 'sequence window')])""",
  """F.key_row(None, y+h-27.4, [('line', (INK, None, 'bold'), 'matched PAS'), ('line', (AXC, None, 'bold'), 'unmatched PAS')], center=x+w/2)
F.key_row(None, y+h-17.4, [('block', WIN, 'sequence window')], center=x+w/2)"""),
 ("""tx = x+128.0; ty = y+3.0
for grp, items in DEFS:
    F.text(tx, ty, f'{grp} group', 7.5, font='Bol'); ty += 9.4
    for k, v in items:
        wk = F.text(tx+4.0, ty, k, 7.0, font='Bol' if False else 'Reg', color=INK)
        F.text(tx+4.0+46.0, ty, v, 7.0, color=MUTED); ty += 8.6
    ty += 1.6
""", "# DEFS (feature groups and definitions) are given in the caption, not drawn\n"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, 'Model genes', 7.5, ha='center')\n", ""),
 ("F.xaxis(bx, y+22.0+3*16.0+6.0, [0, 500, 1000], ['0', '500', '1,000'])", "F.xaxis(bx, y+22.0+3*16.0+6.0, [0, 500, 1000], ['0', '500', '1,000'], title='Model genes')"),
]
