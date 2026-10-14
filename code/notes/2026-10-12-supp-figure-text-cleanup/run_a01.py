"""Render Supplementary Figure S21 (A01) as production does (notes/2026-10-07-a1-figure-revision/run_a01.py: corrected
proximal-overlap table via alt-s21/sources-a1.json), into <run>/out/a01/<tag>. A01.py itself carries the 2026-10-12 text edits."""
import sys
from pathlib import Path
RUN = Path(__file__).resolve().parent; ROOT = RUN.parents[1]
A01 = ROOT/'notes/2026-10-06-alt-catalogue-suppfig/A01.py'
SRC = ROOT/'notes/2026-10-07-a1-figure-revision/alt-s21/sources-a1.json'
out = RUN/'out'/'a01'/(sys.argv[1] if len(sys.argv) > 1 else 'current')
src = A01.read_text()
for old, new in [("E = PinnedTables(BASE/'sources.json')", f"E = PinnedTables(Path({str(SRC)!r}))"),
                 ("OUT = BASE/'candidates'/(sys.argv[1] if len(sys.argv) > 1 else 'run1')", f"OUT = Path({str(out)!r})")]:
    assert src.count(old) == 1, old; src = src.replace(old, new)
sys.argv = [str(A01)]
exec(compile(src, str(A01), 'exec'), {'__name__': '__main__', '__file__': str(A01)})
