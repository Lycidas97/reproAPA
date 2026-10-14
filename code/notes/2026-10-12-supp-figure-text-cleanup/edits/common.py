TARGET = 'notes/2026-10-10-supplement-v11/figures/common.py'
EDITS = [
 ("OUT = FIGDIR.parent/'out'            # defined after the star import (supdraw has its own HERE/OUT)",
  "import os as _os\n# defined after the star import (supdraw has its own HERE/OUT); V11_OUT redirects a rebuild (2026-10-12 text cleanup)\nOUT = Path(_os.environ['V11_OUT']) if _os.environ.get('V11_OUT') else FIGDIR.parent/'out'"),
]
