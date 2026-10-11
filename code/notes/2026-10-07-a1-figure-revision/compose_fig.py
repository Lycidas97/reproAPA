"""Recompose consolidated supplementary figures (plan letters given on the command line) into final-revision/<letter>/
from the A1-corrected supplement bundle (v2). The sealed contract-r2 bundle and consolidated/ are untouched."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'2026-10-06-supplement-polyapipe/tools'))
import supdraw
NEW = HERE/'bundle-supplement-v2'
supdraw.BUNDLE = NEW; supdraw.Data.__init__.__defaults__ = (NEW,)
import compose
compose.OUT = HERE/'final-revision'
for cid in sys.argv[1:]: compose.compose(cid)
