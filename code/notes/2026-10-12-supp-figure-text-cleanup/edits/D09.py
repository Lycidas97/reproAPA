TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D09.py'
EDITS = [
 ("F.text((sx.p0+sx.p1)/2, y+3.0, 'Recurrent genes, observed / null mean', 7.5, ha='center')\nrh = 7.6;",
  "# 2026-10-12: panel titles moved to axis titles or the caption throughout this figure\nrh = 7.1;"),
 ("['.5', '1', '1.5', '2', '2.5', '3', '3.5', '4'])", "['.5', '1', '1.5', '2', '2.5', '3', '3.5', '4'], title='Observed / null mean')"),
 ("H_B = 96.0", "H_B = 92.0"),
 ("F.text(b0+bw_/2, y+3.0, 'Distal | proximal recurrent genes', 7.5, ha='center')\n", ""),
 ("F.text((lx.p0+lx.p1)/2, y+3.0, 'Recurrent genes, |ΔDWUI| ≥ .02 (log)', 7.5, ha='center')\n", ""),
 ("    yy += 15.0\nF.xaxis(lx, yy-4.0, [0, 1, 2, 3], ['1', '10', '100', '1,000'])", "    yy += 13.5\nF.xaxis(lx, yy-4.0, [0, 1, 2, 3], ['1', '10', '100', '1,000'], title='Recurrent genes')"),
 ("rh = 10.6; r0 = y+30.0; bwd = 12.0\nF.text(cx0+cw*2, y+3.0, 'Members kept, same direction', 7.5, ha='center')\n", "rh = 11.6; r0 = y+23.0; bwd = 12.0\n"),
 ("F.text(cx0+cw*j_+1.0, y+11.0, sp_, 7.0); F.text(cx0+cw*j_+1.0, y+18.6, ct_, 7.0, color=MUTED)", "F.text(cx0+cw*j_+1.0, y+3.0, sp_, 7.0); F.text(cx0+cw*j_+1.0, y+10.6, ct_, 7.0, color=MUTED)"),
 ("F.text((ex.p0+ex.p1)/2, y+3.0, 'Human P–M', 7.5, ha='center')\n", ""),
 ("    F.rect(x, yy-6.0, w, 17.6, fc=BAND) if st != 'GSE150311' else None", "    F.rect(x, yy-5.6, w, 16.0, fc=BAND) if st != 'GSE150311' else None"),
 ("    yy += 19.0\nF.xaxis(ex, yy-8.0, [0, 60, 120, 180])", "    yy += 17.0\nF.xaxis(ex, yy-8.0, [0, 60, 120, 180], title='Recurrent genes')"),
 ("F.text((fx.p0+fx.p1)/2+4, y+3.0, 'Guan series', 7.5, ha='center')\n", ""),
 ("F.xaxis(fx, yy+30.0, [0, 700, 1400], ['0', '700', '1,400'])\nF.text((fx.p0+fx.p1)/2, y+h-17.4, 'Human P–S', 7.0, ha='center', color=MUTED)",
  "F.xaxis(fx, yy+30.0, [0, 700, 1400], ['0', '700', '1,400'], title='Recurrent genes')"),
]
