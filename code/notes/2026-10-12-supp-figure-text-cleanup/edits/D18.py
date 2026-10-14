TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D18.py'
EDITS = [
 ("H1, H2, H3 = 132.0, 106.0, 116.0", "H1, H2, H3 = 132.0, 108.0, 116.0   # 2026-10-12: panel titles moved to axis titles or the caption"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, 'Same-PAS designs', 7.5, ha='center')\n", ""),
 ("    F.rect(cx-17.6, y+3.6, 5.6, 5.6, fc=STUDY[hr][1]); F.text(cx-10.4, y+3.0, '×', 7.0, color=MUTED)\n    F.rect(cx-5.0, y+3.6, 5.6, 5.6, fc=STUDY[mr][1]); F.text(cx+2.4, y+3.0, f'{len(z)} genes', 7.0, color=MUTED)",
  "    F.study_pair(cx-9.1, y+6.4, hr, mr)   # gene count (len(z)) moved to the caption"),
 ("F.text((cx_.p0+cx_.p1)/2, y+3.0, 'Rank axis, 25 nt', 7.5, ha='center')\n", ""),
 ("F.text((dx.p0+dx.p1)/2, y+3.0, 'Paired correlation difference, 25 nt', 7.5, ha='center')\n", ""),
 ("Dv = {}; yy = y+23.0", "Dv = {}; yy = y+20.0"),
 ("F.xaxis(dx, yy-6.0, [-0.4, -0.2, 0, 0.2, 0.4], ['−.4', '−.2', '0', '.2', '.4'])", "F.xaxis(dx, yy-6.0, [-0.4, -0.2, 0, 0.2, 0.4], ['−.4', '−.2', '0', '.2', '.4'], title='Correlation difference')"),
 ("    F.text(x+16.0+textwidth(sym, 8.0, 'BolIta')+3.0, y+1.0, f'{chrom} (−), Ensembl 116', 7.0, color=MUTED)\n", "    # chromosome, strand and Ensembl release: caption\n"),
]
