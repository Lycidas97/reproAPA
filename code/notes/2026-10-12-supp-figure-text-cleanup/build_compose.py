"""Compose the v11 supplementary figures that come from the consolidated per-figure scripts, under their v11 numbers.
Output: <run>/out/compose/<key>/SuppFig.S<n>.pdf (+ png, receipt, values). Fig/ is never written here.
Keys: A B C D E J M N P T S Q R (old consolidated plan letters, unchanged bands) and
      V7 = former S9 d-h (D15 a,b,c; D14 e; D15 e)   -- replaces the 2026-10-10 vector split of old S9
      V8 = former S12 f-i (D12 a-d)                   -- replaces the vector split of old S12
      X19 = former S15 a-d (D19 a,b; D20 a,b), cropped into S19 f-i by mech.py
      X20 = former S7 c (D07 c), cropped into S20 d by profile.py
Usage: python3 build_compose.py [key ...]"""
import importlib.util, sys
from pathlib import Path
RUN = Path(__file__).resolve().parent; ROOT = RUN.parents[1]
KEYS = sys.argv[1:]
if 'M' in KEYS:
    assert KEYS == ['M'], 'run M alone (A1-corrected bundle v2, as production)'
    sys.path.insert(0, str(ROOT/'notes/2026-10-06-supplement-polyapipe/tools'))
    import supdraw
    BUNDLE_V2 = ROOT/'notes/2026-10-07-a1-figure-revision/bundle-supplement-v2'
    supdraw.BUNDLE = BUNDLE_V2; supdraw.Data.__init__.__defaults__ = (BUNDLE_V2,)
spec = importlib.util.spec_from_file_location('compose', ROOT/'notes/2026-10-06-supplement-polyapipe/tools/compose.py')
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
C.OUT = RUN/'out'/'compose'
S = C.SPEC
S['V7'] = dict(title='functional enrichment (v11 S7)', bands=[[('D15', 'ab')], [('D15', 'c')], [('D14', 'e'), ('D15', 'e')]],
               patches={'D14': S['I']['patches']['D14'], 'D15': S['I']['patches']['D15']})
S['V8'] = dict(title='APA x GEX (v11 S8)', bands=[[('D12', 'abcd')]], patches=dict(S['K'].get('patches', {})))
S['X19'] = dict(title='former S15 a-d', bands=[[('D19', 'ab')], [('D20', 'ab')]], patches=dict(S['O'].get('patches', {})))
S['X20'] = dict(title='former S7 c', bands=[[('D07', 'C')]], patches=dict(S['G'].get('patches', {})))
V11 = dict(A=1, B=2, C=3, D=4, E=5, J=6, V7=7, V8=8, M=9, N=11, P=13, T=14, S=15, Q=16, R=17, X19='19-fi', X20='20-d')
C.PNUM.clear(); C.PNUM.update(V11)
if __name__ == '__main__':
    for k in (sys.argv[1:] or [k for k in V11]):
        C.compose(k)
