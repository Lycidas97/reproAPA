"""S19 tightening plan, in the coordinates of the 10-12 production S19 (before the 10-13 density edit).
Each region (x0, y0, x1, y1) moves by dy; regions are listed bottom-up within each column for shift order."""
SRC = 'revision-backups/20261013-before-figure-density-fixes/Fig/SuppFig.S19.pdf'
W = 453.54
U = -3.2            # e block up: d-to-e gap back to the 11.3 pt panel gap
EP = 1.0            # e row pitch 7.0 -> 8.0 pt
S3 = U+5*EP         # row 3 (f, g) follows the new bottom of row 2: +1.8
GP = 3.0            # g bar pitch 15 -> 12 pt
S4 = S3-3*GP        # row 4 (h, i) follows the new bottom of row 3: -7.2
HP = 2.0            # h row pitch 13.5 -> 11.5 pt
def regions():
    R = []
    # row 4, bottom-up
    R.append(((0.0, 471.3, 272.0, 504.9), S4-3*HP, 'h axis, title and key'))
    for k in (3, 2, 1, 0):
        R.append(((0.0, 419.8+13.5*k, 272.0, 431.2+13.5*k), S4-HP*k, f'h row {k}'))
    R.append(((0.0, 398.2, 272.0, 419.7), S4, 'h letter and header'))
    R.append(((272.0, 486.0, W, 504.9), S4-3*HP, 'i key'))
    R.append(((272.0, 398.2, W, 485.9), S4, 'i letter, headers and heat map'))
    # row 3, bottom-up
    R.append(((162.0, 364.9, W, 398.1), S3-3*GP, 'g axis, title and key'))
    for k in (3, 2, 1, 0):
        R.append(((162.0, 309.9+15*k, W, 320.1+15*k), S3-GP*k, f'g bar {k}'))
    R.append(((162.0, 289.9, W, 309.8), S3, 'g letter'))
    R.append(((0.0, 289.9, 162.0, 398.1), S3, 'f'))
    # row 2 right column (e), bottom-up
    R.append(((296.5, 264.45, W, 289.85), U+5*EP, 'e axis, ticks and title'))
    for i in range(5, -1, -1):
        top = 221.6+7.0*i
        R.append(((296.5, top-0.02, W, top+7.02), U+EP*i, f'e row {i}'))
    R.append(((296.5, 197.5, W, 221.55), U, 'e letter and q head'))
    return R
