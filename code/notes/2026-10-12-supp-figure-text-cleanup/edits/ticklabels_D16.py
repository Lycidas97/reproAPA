TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D16.py'
# 2026-10-12 audit: P-S '70' touched P-M '0' (S9a, 3.75 pt gap, read as '700'); the left subplot drops its right-end label
EDITS = [("    ab_ = F.xaxis(sx, y+20.0+7*8.6+3.0+5.4, [0, 35, 70])",
          "    ab_ = F.xaxis(sx, y+20.0+7*8.6+3.0+5.4, [0, 35, 70], ['0', '35', '' if j == 0 else '70'])")]
