"""Re-render Supplementary Figure S21 (A01, printed S21) with the corrected proximal-overlap table into alt-s21/out.
The A01 script and its pinned sources.json are untouched; only the manifest path and output folder are swapped."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
A01 = HERE.parent/'2026-10-06-alt-catalogue-suppfig/A01.py'
src = A01.read_text()
for old, new in [("E = PinnedTables(BASE/'sources.json')", f"E = PinnedTables(Path({str(HERE/'alt-s21/sources-a1.json')!r}))"),
                 ("OUT = BASE/'candidates'/(sys.argv[1] if len(sys.argv) > 1 else 'run1')", f"OUT = Path({str(HERE/'alt-s21/out')!r})")]:
    assert src.count(old) == 1, old; src = src.replace(old, new)
sys.argv = [str(A01)]
exec(compile(src, str(A01), 'exec'), {'__name__': '__main__', '__file__': str(A01)})
