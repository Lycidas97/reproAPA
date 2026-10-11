"""Corrected Table S5b sources after the A1 eligibility correction (2026-10-07).

- cross-overlap-five-null-results.tsv : polyapipe-null-results.tsv verbatim (all five nulls, legacy and corrected
  expectations, expected events on ineligible genes); replaces overlap_null_results.tsv and the 10-06
  expression-coverage-null-comparison.tsv.
- overlap-threshold-layer-sensitivity.tsv : polyapipe-grid-results.tsv verbatim.
- expression-coverage-{label-mobility,balance,calibration-power}.tsv : P-S rows kept from 10-06 (the corrected P-S
  strata, expectations and draws are identical to 10-06: expected and null mean unchanged, no ineligible events);
  P-M rows from local-extended-diagnostics (local rerun of the corrected infer(); P-M expectations, null means and
  p equal the server receipt to full precision). Candidate movable fractions are taken from the server receipt.
- missingness and tercile edges do not depend on the strata and are reused unchanged.
"""
import csv, hashlib, shutil
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
A1 = ROOT/'notes/2026-10-07-cross-overlap-null-a1'; OLD = ROOT/'notes/2026-10-06-cross-overlap-expression-coverage/supplement-table-S5b'
LOC = HERE/'local-extended-diagnostics'; OUT = HERE/'table-S5b'
shutil.copyfile(A1/'polyapipe-null-results.tsv', OUT/'cross-overlap-five-null-results.tsv')
shutil.copyfile(A1/'polyapipe-grid-results.tsv', OUT/'overlap-threshold-layer-sensitivity.tsv')
srv = pd.read_csv(A1/'polyapipe-null-results.tsv', sep='\t', dtype=str)
loc = pd.read_csv(LOC/'expression-coverage-null-results.tsv', sep='\t', dtype=str)
for m in ['Q5_expression_coverage', 'Range_expression_coverage']:   # P-M identity check, and P-S unchanged
    s = srv[(srv.contrast == 'PM') & (srv.model == m)].iloc[0]; l = loc[(loc.contrast == 'PM') & (loc.model == m)].iloc[0]
    for a, b in [('analytic_expected', 'analytic_expected'), ('null_mean', 'null_mean'), ('p', 'p'), ('null_sd', 'null_sd')]:
        assert abs(float(s[a])-float(l[b])) < 1e-9, (m, a)
    s = srv[(srv.contrast == 'PS') & (srv.model == m)].iloc[0]
    assert abs(float(s.analytic_expected)-float(s.legacy_analytic_expected)) < 1e-12 and float(s.h_legacy_expected_ineligible_events) == 0
def merge(name, pm):
    old = pd.read_csv(OLD/name, sep='\t', dtype=str, keep_default_na=False)
    new = pd.concat([old[old.contrast == 'PS'], pm[old.columns]], ignore_index=True)
    new.to_csv(OUT/name, sep='\t', index=False)
mob = []
for m in ['Q5_expression_coverage', 'Range_expression_coverage']:
    s = srv[(srv.contrast == 'PM') & (srv.model == m)].iloc[0]; l = loc[(loc.contrast == 'PM') & (loc.model == m)].iloc[0]
    mob.append(dict(contrast='PM', model=m, fixed_overlap=l.fixed_overlap, observed=l.observed, null_unique=l.null_unique,
                    null_sd=l.null_sd, h_candidate_movable_fraction=s.h_candidate_movable_fraction,
                    h_expected_events_below_two_actual_studies='0', m_candidate_movable_fraction=s.m_candidate_movable_fraction,
                    m_expected_events_below_two_actual_studies='0'))
merge('expression-coverage-label-mobility.tsv', pd.DataFrame(mob))
for name in ['expression-coverage-balance.tsv', 'expression-coverage-calibration-power.tsv']:
    l = pd.read_csv(LOC/name.replace('expression-coverage-', ''), sep='\t', dtype=str, keep_default_na=False)
    merge(name, l[l.contrast == 'PM'])
for name in ['expression-coverage-missingness.tsv', 'expression-coverage-tercile-edges.tsv']:
    shutil.copyfile(OLD/name, OUT/name)
for p in sorted(OUT.iterdir()):
    print(p.name, sum(1 for _ in open(p))-1, hashlib.sha256(p.read_bytes()).hexdigest()[:16])
