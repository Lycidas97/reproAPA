TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D24.py'
# 2026-10-12: wider 'Hs/Mm' row labels overlapped the test names in d by up to 2.7 pt; test names and axis move 4 pt right
EDITS = [
 ("dx = Lin(-1.0, 1.1, x+76.0, x+w-26.0)", "dx = Lin(-1.0, 1.1, x+80.0, x+w-26.0)"),
 ("        F.text(x+72.0, yy-3.4, TL[t_], 6.8, ha='right', color=MUTED)", "        F.text(x+76.0, yy-3.4, TL[t_], 6.8, ha='right', color=MUTED)"),
]
