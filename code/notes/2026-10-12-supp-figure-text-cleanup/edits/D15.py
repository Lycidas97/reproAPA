TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D15.py'
EDITS = [
 ("H1, H2, H3 = 79.0, 170.0, 140.0", "H1, H2, H3 = 86.0, 170.0, 140.0   # 2026-10-12: panel titles -> x-axis titles"),
 ("F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Annotated genes (GO)', 7.5, ha='center')\n", ""),
 ("F.xaxis(ax_, y+21.0+3*9.4+4.6, [0, 5000, 10000, 15000], ['0', '5k', '10k', '15k'])", "F.xaxis(ax_, y+21.0+3*9.4+4.6, [0, 5000, 10000, 15000], ['0', '5k', '10k', '15k'], title='Annotated genes (GO)')"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, 'q < .05 positions', 7.5, ha='center')\n", ""),
 ("F.xaxis(bx, y+21.0+2*13.6+6.4, [0, 8, 16])", "F.xaxis(bx, y+21.0+2*13.6+6.4, [0, 8, 16], title='Positions with q < .05')"),
 ("F.text((c1.p0+c1.p1)/2, y+3.0, 'Expected → observed genes', 7.5, ha='center')\nF.text((c2.p0+c2.p1)/2, y+3.0, 'Matched O/E', 7.5, ha='center')\n", ""),
 ("F.xaxis(c1, ab, [0, 80, 160, 240]); F.xaxis(c2, ab, [1, 1.5, 2])", "F.xaxis(c1, ab, [0, 80, 160, 240], title='Genes'); F.xaxis(c2, ab, [1, 1.5, 2], title='Matched O/E')"),
]
# follow-up (sed): H2 170 -> 181 (room for c axis titles)
