TARGET = 'notes/2026-10-06-supplement-polyapipe/figures/D22.py'
# 2026-10-12 balance: content did not refill the panel after title/note rows were removed (S16d, about 18 pt between the axis title and the key after the row height grew)
EDITS = [("        ry = yy+k*13.0", "        ry = yy+k*14.0"), ("    yy += 28.0\nF.xaxis(bx, yy-12.0", "    yy += 35.0\nF.xaxis(bx, yy-12.0")]
