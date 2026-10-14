TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D24.py'
EDITS = [
 ("F.text((ax_.p0+ax_.p1)/2, y+3.0, f'Mouse seeds ({len(m)})', 7.5, ha='center')\n", "# 2026-10-12: panel titles moved to axis titles or the caption throughout this figure\n"),
 ("F.text((bx.p0+bx.p1)/2, y+3.0, 'Matched / all lost bp', 7.5, ha='center')\n", ""),
 ("['0', '.05', '.1'], 'Fraction')", "['0', '.05', '.1'], 'Matched / all lost bp')"),
 ("F.text((cx_.p0+cx_.p1)/2, y+3.0, 'Lost − control sites/kb', 7.5, ha='center')\n", ""),
 ("F.text((dx.p0+dx.p1)/2, y+3.0, 'Density contrast', 7.5, ha='center')\n", ""),
 ("F.text((ex.p0+ex.p1)/2, y+3.0, 'Mouse somatic − PSC', 7.5, ha='center')\n", ""),
 ("F.xaxis(ex, yy-4.0, [-0.2, -0.1, 0, 0.1], ['−.2', '−.1', '0', '.1'])", "F.xaxis(ex, yy-4.0, [-0.2, -0.1, 0, 0.1], ['−.2', '−.1', '0', '.1'], title='Somatic − PSC')"),
 ("F.text((fx.p0+fx.p1)/2, y+3.0, 'Mouse P–S', 7.5, ha='center')\n", ""),
]
# follow-up (python replace): row pitch c, d 11.4 -> 11.0, e 8.0 -> 7.6 (room for axis titles without taller rows)
