TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D20.py'
EDITS = [
 ("F.text((lx.p0+lx.p1)/2, y+3.0, 'LRT, All cis vs M1', 7.5, ha='center')\n", "# 2026-10-12: panel title (LRT, All cis vs M1) moved to the caption\n"),
 ("F.text(cx0+4*cw, b0+4*rh+2.0, 'τ; outlined = primary test', 7.0, ha='center', color=MUTED)\nF.cbar_key(cx0+4*cw, y+h-7.4, 'q, 32 tests', '1', '≤ .1', lambda f_: seq(f_), wbar=40)",
  "F.text(cx0-2.0, y+13.4, 'τ', 7.0, ha='right', color=MUTED)   # outlined-cell and 32-test notes: caption\nF.cbar_key(cx0+4*cw, y+h-7.4, 'q', '1', '≤ .1', lambda f_: seq(f_), wbar=40)"),
]
