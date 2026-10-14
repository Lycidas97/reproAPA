TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D16.py'
# 2026-10-12 balance: content did not refill the panel after title/note rows were removed (S9c, 12.5 pt empty at the bottom)
EDITS = [("    yy = y+21.0+i*14.4", "    yy = y+21.0+i*16.9"),
         ("F.line([qx(0.05)]*2, [y+15, y+21.0+4*14.4+5], lw='ref', color=AXC, ls='ref')", "F.line([qx(0.05)]*2, [y+15, y+21.0+4*16.9+5], lw='ref', color=AXC, ls='ref')"),
         ("F.xaxis(qx, y+21.0+4*14.4+7.0, [0, 0.5, 1]", "F.xaxis(qx, y+21.0+4*16.9+7.0, [0, 0.5, 1]")]
