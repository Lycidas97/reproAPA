TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D16.py'
EDITS = [
 ("    F.xaxis(sx, y+20.0+7*8.6+3.0+5.4, [0, 35, 70], title='Shared proximal candidates')",
  "    ab_ = F.xaxis(sx, y+20.0+7*8.6+3.0+5.4, [0, 35, 70])\nF.text(x+48+pw+5, ab_+1.6, 'Shared proximal candidates', 8.0, ha='center')   # one axis title for both subplots"),
]
