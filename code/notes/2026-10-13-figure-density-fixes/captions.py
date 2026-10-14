r"""Figure 4c now has a logarithmic MPRO axis; the caption says so. Must match exactly once; --apply writes."""
import sys
from pathlib import Path
F = Path(__file__).resolve().parents[2]/'main/results/local-3-utr-remodeling-escapes-gene-level-summaries-and-changes-predicted-regulatory-content.tex'
OLD = r"(filled, the 8 with a strict recurrent interval; shading, $|\Delta\mathrm{DWUI}|<0.01$ in both)."
NEW = r"(filled, the 8 with a strict recurrent interval; shading, $|\Delta\mathrm{DWUI}|<0.01$ in both; logarithmic MPRO axis)."
s = F.read_text(); n = s.count(OLD); print('matches', n); assert n == 1
if '--apply' in sys.argv: F.write_text(s.replace(OLD, NEW)); print('written')
