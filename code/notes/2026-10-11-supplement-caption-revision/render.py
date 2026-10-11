"""Rebuild S1/S2 with native state IDs using the existing consolidated renderer."""
import importlib.util
from pathlib import Path

RUN = Path(__file__).resolve().parent
ROOT = RUN.parents[1]
path = ROOT / 'notes/2026-10-06-supplement-polyapipe/tools/compose.py'
spec = importlib.util.spec_from_file_location('compose', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.OUT = RUN / 'rendered'
for key in ('A', 'B'):
    module.compose(key)
