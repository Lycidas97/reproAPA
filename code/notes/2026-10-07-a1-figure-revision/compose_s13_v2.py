"""Recompose Supplementary Figure S13 (plan letter M) from the A1-corrected supplement bundle into s13-final/.
The sealed contract-r2 bundle and consolidated/M are untouched."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent/'2026-10-06-supplement-polyapipe/tools'
sys.path.insert(0, str(TOOLS))
import supdraw
NEW = HERE/'bundle-supplement-v2'
supdraw.BUNDLE = NEW
supdraw.Data.__init__.__defaults__ = (NEW,)
import compose
compose.OUT = HERE/'s13-final-v2'
compose.compose('M')
