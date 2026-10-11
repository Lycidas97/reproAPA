"""Inspect eligibility and schemas without inspecting new QKI/EMT target hits."""
from pathlib import Path
import datetime
import hashlib
import json
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE / 'source/xlrd-2.0.2-py2.py3-none-any.whl'))
import xlrd

def read_sheet(name):
    sheet = xlrd.open_workbook(HERE / 'source' / name).sheet_by_index(0)
    return pd.DataFrame([sheet.row_values(i) for i in range(2, sheet.nrows)],
                        columns=sheet.row_values(1))

def main():
    roster = pd.read_csv(ROOT / 'notes/2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz', sep='\t')
    roster = roster[roster.contrast.eq('P_M')].copy()
    assert len(roster) == 751 and roster.foreground_002.sum() == 27 and roster.foreground_005.sum() == 7
    assert roster.human_gene_id.is_unique and roster.mouse_gene_id.is_unique
    checks = {}
    for name, label in [('KRNB_A_2294222_SM0800.xls', 'EMT'), ('KRNB_A_2294222_SM0790.xls', 'QKI_KD')]:
        frame = read_sheet(name)
        duplicate = frame.name.duplicated(keep=False)
        finite = np.isfinite(pd.to_numeric(frame.effect, errors='coerce')) & np.isfinite(pd.to_numeric(frame.fdr_zero, errors='coerce'))
        valid = frame[~duplicate & finite]
        covered = roster.human_gene_id.isin(valid.name)
        checks[label] = dict(rows=len(frame), unique_gene_ids=frame.name.nunique(),
            duplicate_gene_ids_excluded=int(frame.loc[duplicate, 'name'].nunique()),
            finite_unambiguous_genes=len(valid), background_eligible=int(covered.sum()),
            core_eligible=int((roster.foreground_002 & covered).sum()),
            strict7_eligible=int((roster.foreground_005 & covered).sum()), columns=list(frame.columns),
            sha256=hashlib.sha256((HERE / 'source' / name).read_bytes()).hexdigest())
    event_path = HERE.parent / 'event_inventory.tsv.gz'
    events = pd.read_csv(event_path, sep='\t')
    for field, value, label in [('driver_class', 'strict_tandem_UTR', 'all_route_drivers_strict_tandem'),
                                ('full_primary_axis_tandem', True, 'all_route_axes_strict_tandem')]:
        available = events.groupby(['species', 'gene_id'])[field].apply(lambda x: bool(x.eq(value).all()))
        flags = np.array([bool(available.get(('human', r.human_gene_id), False) and available.get(('mouse', r.mouse_gene_id), False))
                          for r in roster.itertuples()])
        checks[label] = dict(background_eligible=int(flags.sum()), core_eligible=int((roster.foreground_002 & flags).sum()))
    checks.update(timestamp_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  no_target_hit_or_effect_sign_summary_before_protocol=True,
                  event_inventory_sha256=hashlib.sha256(event_path.read_bytes()).hexdigest())
    (HERE / 'preflight.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(checks, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
