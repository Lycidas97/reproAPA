TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D12.py'
EDITS = [
 ("F.text((sx.p0+sx.p1)/2, y+3.0, 'Tested class, observed / null mean', 7.5, ha='center')\n", "# 2026-10-12: panel titles moved to axis titles or the caption throughout this figure\n"),
 ("rh = 10.4; r0 = y+24.0", "rh = 10.4; r0 = y+21.0"),
 ("ab = F.xaxis(sx, r0+rh*7+4+rh/2+1.0, [1, 1.5, 2], ['1', '1.5', '2'])", "ab = F.xaxis(sx, r0+rh*7+4+rh/2+1.0, [1, 1.5, 2], ['1', '1.5', '2'], title='Observed / null mean')"),
 ("F.text(gx0+cw*2+sg/2, y+3.0, 'Comparable genes', 7.5, ha='center')\n", ""),
 ("    F.text(sx0+cw, y+12.6, f'{sp.capitalize()} {comp:,}', 7.0, ha='center')\n    F.text(sx0+cw/2, y+21.0, 'Same', 6.8, ha='center', color=MUTED); F.text(sx0+1.5*cw, y+21.0, 'Opp.', 6.8, ha='center', color=MUTED)\nF.text(gx0-3.0, y+21.0, 'ΔGEX', 6.8, ha='right', color=MUTED)",
  "    F.text(sx0+cw, y+3.0, f'{sp.capitalize()} {comp:,}', 7.0, ha='center')\n    F.text(sx0+cw/2, y+11.4, 'Same', 6.8, ha='center', color=MUTED); F.text(sx0+1.5*cw, y+11.4, 'Opp.', 6.8, ha='center', color=MUTED)\nF.text(gx0-3.0, y+11.4, 'ΔGEX', 6.8, ha='right', color=MUTED)"),
 ("    by0 = y+29.0+b_*(2*ch+6.0)", "    by0 = y+19.4+b_*(2*ch+6.0)"),
 ("F.text((ax_.p0+ax_.p1)/2, y+3.0, 'ΔDWUI', 7.5, ha='center')\n", ""),
 ("F.text(min((gx.p0+gx.p1)/2, x+w-0.3-textwidth(_gex_label, 7.5)/2), y+3.0, _gex_label, 7.5, ha='center')\nDv = {}; yy = y+22.0",
  "Dv = {}; yy = y+18.0"),
 ("F.xaxis(ax_, yy-7.0, [-0.3, -0.15, 0], ['−.3', '−.15', '0'])\nF.xaxis(gx, yy-7.0, [-1, 0, 1, 2, 3])",
  "F.xaxis(ax_, yy-7.0, [-0.3, -0.15, 0], ['−.3', '−.15', '0'], title='ΔDWUI')\nF.xaxis(gx, yy-7.0, [-1, 0, 1, 2, 3], title=_gex_label)"),
]
