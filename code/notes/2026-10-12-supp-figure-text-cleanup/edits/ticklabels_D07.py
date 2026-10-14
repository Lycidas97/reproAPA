TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D07.py'
# 2026-10-12 audit: the '1' of each step panel touched the '0' of the next panel (S20d, gaps 1.8-3.8 pt); right-end labels dropped
EDITS = [("for s_ in cx:\n    F.xaxis(s_, ax_y, [0, 0.5, 1], ['0', '.5', '1'])",
          "for s_ in cx:\n    F.xaxis(s_, ax_y, [0, 0.5, 1], ['0', '.5', ''])")]
