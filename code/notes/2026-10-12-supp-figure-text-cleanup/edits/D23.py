TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D23.py'
EDITS = [
 ("H1, H2 = 132.0, 150.0", "H1, H2 = 141.0, 150.0   # 2026-10-12: panel titles moved to axis titles or the caption"),
 ("F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Changed − stable sites per kb', 7.5, ha='center')\n", ""),
 ("F.xaxis(ax_, yy-4.0, [-50, 0, 50, 100], ['−50', '0', '50', '100'])", "F.xaxis(ax_, yy-4.0, [-50, 0, 50, 100], ['−50', '0', '50', '100'], title='Changed − stable sites per kb')"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, 'Mouse P–M, α = .05', 7.5, ha='center')\n", ""),
]
