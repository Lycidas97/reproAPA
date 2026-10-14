TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D14.py'
EDITS = [
 ("F.text((ex.p0+ex.p1)/2, y+3.0, 'APA rank effect', 7.5, ha='center')\n", "# 2026-10-12: panel title folded into the axis titles\n"),
 ("F.yaxis(ey, ex.p0, [-0.4, 0, 0.4], ['−.4', '0', '.4'], 'Expanded layer')", "F.yaxis(ey, ex.p0, [-0.4, 0, 0.4], ['−.4', '0', '.4'], 'Rank effect, expanded')"),
 ("F.xaxis(ex, ey.p0+1.0, [-0.4, 0, 0.4], ['−.4', '0', '.4'], 'Supported layer')", "F.xaxis(ex, ey.p0+1.0, [-0.4, 0, 0.4], ['−.4', '0', '.4'], 'Rank effect, supported')"),
]
