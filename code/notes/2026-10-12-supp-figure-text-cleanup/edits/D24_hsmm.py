TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D24.py'
# 2026-10-12 caption audit: row prefixes 'H '/'M ' -> 'Hs '/'Mm ' (M collided with the M state; Hs/Mm as in S10c, S18d, S19e)
EDITS = [
 ("        if k == 0: F.text(x+0.4, yy-3.4, CL[c].replace('Mouse ', 'M '), 7.0)",
  "        if k == 0: F.text(x+0.4, yy-3.4, CL[c].replace('Mouse ', 'Mm '), 7.0)"),
 ("        if k == 0: F.text(x+0.4, yy+rp/2-3.4, CL[c].replace('Human ', 'H ').replace('Mouse ', 'M '), 7.0)",
  "        if k == 0: F.text(x+0.4, yy+rp/2-3.4, CL[c].replace('Human ', 'Hs ').replace('Mouse ', 'Mm '), 7.0)"),
]
