TARGET = 'notes/2026-10-12-supp-figure-text-cleanup/build_compose.py'
# 2026-10-12: key M (v11 S9) must read the A1-corrected bundle v2 like production (notes/2026-10-08-caption-review-fixes/run_compose.py);
# the default contract-r2 bundle holds the 10-06 cross-overlap tables (S9a q16, S9f power)
EDITS = [("spec = importlib.util.spec_from_file_location('compose', ROOT/'notes/2026-10-06-supplement-polyapipe/tools/compose.py')",
          "KEYS = sys.argv[1:]\n"
          "if 'M' in KEYS:\n"
          "    assert KEYS == ['M'], 'run M alone (A1-corrected bundle v2, as production)'\n"
          "    sys.path.insert(0, str(ROOT/'notes/2026-10-06-supplement-polyapipe/tools'))\n"
          "    import supdraw\n"
          "    BUNDLE_V2 = ROOT/'notes/2026-10-07-a1-figure-revision/bundle-supplement-v2'\n"
          "    supdraw.BUNDLE = BUNDLE_V2; supdraw.Data.__init__.__defaults__ = (BUNDLE_V2,)\n"
          "spec = importlib.util.spec_from_file_location('compose', ROOT/'notes/2026-10-06-supplement-polyapipe/tools/compose.py')")]
