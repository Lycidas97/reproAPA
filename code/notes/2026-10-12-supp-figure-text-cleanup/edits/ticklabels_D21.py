TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D21.py'
# 2026-10-12 audit: '−200' and '−50' overlapped on the symmetric log axis of S13e (1.6 pt gap); the ±50 ticks are dropped
EDITS = [("F.xaxis(ex, yy-2.0, [fl(-200), fl(-50), fl(-5), 0, fl(5), fl(50), fl(200)], ['−200', '−50', '−5', '0', '5', '50', '200'], title='Weighted predicted-site change')",
          "F.xaxis(ex, yy-2.0, [fl(-200), fl(-5), 0, fl(5), fl(200)], ['−200', '−5', '0', '5', '200'], title='Weighted predicted-site change')")]
