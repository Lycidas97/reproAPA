TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D21.py'
EDITS = [
 ("H1, H2, H3 = 96.0, 176.0, 132.0", "H1, H2, H3 = 99.0, 176.0, 142.0   # 2026-10-12: panel titles moved to axis titles or the caption"),
 ("F.text((ax_.p0+ax_.p1)/2, y+3.0, 'Genes per cohort', 7.5, ha='center')\n", ""),
 ("F.xaxis(ax_, y+21.0+3*14.0+7.0, [0, 1000, 2000], ['0', '1,000', '2,000'])", "F.xaxis(ax_, y+21.0+3*14.0+7.0, [0, 1000, 2000], ['0', '1,000', '2,000'], title='Genes')"),
 ("F.text(x+w-1.0, y+3.0, 'Low-net-shift screen, mouse P–M', 7.5, ha='right')\n", ""),
 ("F.text(x+62.0, y+3.0, 'Same sign in both studies', 7.5, ha='center')\n", ""),
 ("F.text((ex.p0+ex.p1)/2, y+3.0, 'Weighted predicted-site change', 7.5, ha='center')\nF.text((cxw.p0+cxw.p1)/2, y+3.0, 'Genes: loss | gain', 7.5, ha='center')\n", ""),
 ("['−200', '−50', '−5', '0', '5', '50', '200'])\nF.xaxis(cxw, yy-2.0, [-200, 0, 200], ['200', '0', '200'])",
  "['−200', '−50', '−5', '0', '5', '50', '200'], title='Weighted predicted-site change')\nF.xaxis(cxw, yy-2.0, [-200, 0, 200], ['200', '0', '200'], title='Genes')"),
 ("F.text((fx.p0+fx.p1)/2, y+3.0, '50 − 20 UMI estimate', 7.5, ha='center')\n", ""),
]
