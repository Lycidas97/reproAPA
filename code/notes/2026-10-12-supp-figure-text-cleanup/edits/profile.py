TARGET = 'notes/2026-10-10-supplement-v11/figures/profile.py'
EDITS = [
 ("""S7 = ROOT/'revision-backups/20261011-before-v11-apply/Fig/SuppFig.S7.pdf'; H7 = 384.4
R2 = (0.0, 0.0, PAGE_W, 322.0)                  # former Figure 2 rows a–c incl. their keys
R7 = (0.0, 125.8, PAGE_W, 273.6)                # former S7 c
DY7 = 322.0+GAP-126.3
H = R7[3]+DY7+4.5
dst = new_page(H)
place(dst, H, crop(F2, H2, R2), R2, 0.0, 0.0)
place(dst, H, crop(S7, H7, R7, extra_clear=[(4.0, 126.0, 13.5, 138.5)]), R7, 0.0, R7[1]+DY7)
F = supdraw.Fig(H, name('profile'))
F.t.text_top(4.5+0.2, 126.3+0.4+DY7, 'd', 11.0, font='Bol')""",
  """import os
R2 = (0.0, 0.0, PAGE_W, 322.0)                  # former Figure 2 rows a–c incl. their keys
# 2026-10-12 text cleanup: the formula block of former Figure 2b (L, D, T definitions) is removed (now in the caption), and
# d is cropped from a composed page holding only former S7 c (consolidated G, D07 C) rebuilt without its column title (V11_S20_D).
F2_CLEAR = [(121.0, 188.8, 215.0, 218.6)]       # 'L = |R−S| + |E−R| + |P−E|', 'D = |P−S|', 'T = 1 − D/L'
if os.environ.get('V11_S20_D'):
    S7 = Path(os.environ['V11_S20_D'])
    with pikepdf.open(S7) as _p: H7 = float(_p.pages[0].MediaBox[3])
    R7 = (0.0, 4.0, PAGE_W, H7-4.5); L7 = (4.0, 4.0, 13.5, 16.5)
else:
    S7 = ROOT/'revision-backups/20261011-before-v11-apply/Fig/SuppFig.S7.pdf'; H7 = 384.4
    R7 = (0.0, 125.8, PAGE_W, 273.6); L7 = (4.0, 126.0, 13.5, 138.5)
DY7 = 322.0+GAP-(R7[1]+0.5)
H = R7[3]+DY7+4.5
dst = new_page(H)
place(dst, H, crop(F2, H2, R2, extra_clear=F2_CLEAR), R2, 0.0, 0.0)
place(dst, H, crop(S7, H7, R7, extra_clear=[L7]), R7, 0.0, R7[1]+DY7)
F = supdraw.Fig(H, name('profile'))
F.t.text_top(4.5+0.2, R7[1]+0.5+0.4+DY7, 'd', 11.0, font='Bol')"""),
 ("dict(source='Fig/Figure.2.pdf', source_sha=sha(F2), region=R2, dy=0.0, letters='a, b, c unchanged'),",
  "dict(source='Fig/Figure.2.pdf', source_sha=sha(F2), region=R2, dy=0.0, letters='a, b, c unchanged', cleared=F2_CLEAR),"),
 ("dict(source='Fig/SuppFig.S7.pdf', source_sha=sha(S7), region=R7, dy=DY7, letters={'c': 'd'})],",
  "dict(source=str(S7), source_sha=sha(S7), region=R7, dy=DY7, letters={'c': 'd'})],"),
]
