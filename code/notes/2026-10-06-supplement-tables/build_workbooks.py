#!/usr/bin/env python3
"""Write Supplementary Tables S1-S8 as xlsx workbooks from table-plan.tsv / table-sources.tsv (2026-10-06).

One workbook per table. First tab README (sheet list, then one row per data tab with source path, SHA-256,
row filter, rows and columns). Then one tab per `layer == sheet` source file, named '<sheet>.<i> <stem>'.
Values are copied, not recomputed; SHA-256 of every source is checked against table-sources.tsv first.
"""
import csv, datetime, gzip, hashlib, json, re, sys
from pathlib import Path
import pandas as pd
import xlsxwriter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BUNDLE = ROOT / 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables'
OUT = HERE / 'workbooks'
EXCEL_MAX = 1048575
MAXSTR = 32767

# row filters (path -> (description, function)); descriptions are written into README
FILTERS = {
    'global/background_gene_membership.tsv.gz': (
        "background == 'common_one_to_one_ortholog_all_settings' (one-to-one ortholog background: 1,306 genes per P-S route and setting, 208 per P-M)",
        lambda d: d[d['background'] == 'common_one_to_one_ortholog_all_settings']),
    'recurrence/recurrence_members.tsv.gz': (
        "series_matching == False (the True rows are identical apart from this flag)",
        lambda d: d[d['series_matching'] == 'False']),
    'apa_gex/current_candidate_shared_state_relations.tsv.gz': (
        "status == 'measured' (relationship evaluable in both studies of the axis; unmeasured rows are in the source-data package)",
        lambda d: d[d['status'] == 'measured']),
    'same_pas/v2_primary_same_axis_gene_effects.tsv.gz': (
        "tolerance == 25 (primary 25-nt matching; the 50-nt rows are in the source-data package)",
        lambda d: d[d['tolerance'] == '25']),
}

GENERAL = ("Values are copied from the polyApipe v2-primary analysis outputs (contract r2 bundle) without "
           "recomputation. Each data tab holds one source file; the README rows give its path, SHA-256 and any "
           "row filter. Unfiltered sources and figure-plotting data are in the Supplementary Source Data package.")


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def src_path(r):
    if r['role'] == 'core_list':
        return BUNDLE / r['path']
    return ROOT / r['path']


LEADZERO = re.compile(r'^[+-]?0\d')


def typed(d):
    """dtype=str frame -> list of per-column kinds ('num' or 'str')."""
    kinds = []
    for c in d.columns:
        v = d[c]
        ne = v[v != '']
        kind = 'str'
        if len(ne) and not ne.str.match(LEADZERO).any():
            try:
                pd.to_numeric(ne)
                kind = 'num'
            except (ValueError, TypeError):
                pass
        kinds.append(kind)
    return kinds


# tab-name prefixes for sources whose file names collide within a sheet
LABEL = {
    'measurement/nulls/recurrence_summary.tsv': 'nulls',
    'measurement/same-state/recurrence_summary.tsv': 'same-state',
    'notes/2026-10-05-dual-main-concordance/inputs/polyapipe_common/global/global_direction_summary.tsv': 'polyApipe',
    'notes/2026-10-05-dual-main-concordance/inputs/polyapipe_common/profiles/gene_membership_k3.tsv.gz': 'polyApipe',
    'notes/2026-10-05-dual-main-concordance/inputs/reject_original_point_terminal_p50/global/global_direction_summary.tsv': 'alt',
    'notes/2026-10-05-dual-main-concordance/inputs/reject_original_point_terminal_p50/profiles/gene_membership_k3.tsv.gz': 'alt',
}


def tabname(sheet, i, n, path):
    stem = Path(path).name
    for ext in ('.gz', '.tsv', '.json'):
        stem = stem[:-len(ext)] if stem.endswith(ext) else stem
    stem = re.sub(r'^S\d+-', '', stem)          # upstream table/S1-datasets.tsv -> datasets
    if path in LABEL:
        stem = LABEL[path] + ' ' + stem
    stem = re.sub(r'[\[\]:*?/\\]', '_', stem)
    head = sheet if n == 1 else f'{sheet}.{i}'
    return f'{head} {stem}'[:31]


def main():
    plan = list(csv.DictReader(open(HERE / 'table-plan.tsv'), delimiter='\t'))
    srcs = [r for r in csv.DictReader(open(HERE / 'table-sources.tsv'), delimiter='\t') if r['layer'] == 'sheet']
    OUT.mkdir(exist_ok=True)
    log, issues = [], []
    today = datetime.date.today().isoformat()
    for t in dict.fromkeys(r['table'] for r in plan):
        tplan = [r for r in plan if r['table'] == t]
        title = tplan[0]['table_title']
        wb = xlsxwriter.Workbook(OUT / f'Table_{t}.xlsx', {'strings_to_urls': False, 'strings_to_formulas': False,
                                                         'strings_to_numbers': False, 'nan_inf_to_errors': True,
                                                         'constant_memory': True})
        b = wb.add_format({'bold': True}); h1 = wb.add_format({'bold': True, 'font_size': 13})
        hdr = wb.add_format({'bold': True, 'bottom': 1, 'bg_color': '#EEEEEE'})
        wrap = wb.add_format({'text_wrap': True, 'valign': 'top'})
        rd = wb.add_worksheet('README')
        # README is written after the data tabs are known; collect rows first
        tabs = []
        for s in tplan:
            files = [r for r in srcs if r['sheet'] == s['sheet']]
            for i, r in enumerate(files, 1):
                p = src_path(r)
                cur = sha(p)
                if cur != r['sha256']:
                    issues.append(f"sha-mismatch {r['sheet']} {r['path']}")
                d = pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False, na_filter=False,
                                compression='gzip' if r['path'].endswith('.gz') else None)
                n_src = len(d)
                fdesc = ''
                if r['path'] in FILTERS:
                    fdesc, fn = FILTERS[r['path']]
                    d = fn(d)
                if len(d) > EXCEL_MAX:
                    issues.append(f"too many rows {r['sheet']} {r['path']} {len(d)}")
                    continue
                name = tabname(s['sheet'], i, len(files), r['path'])
                ws = wb.add_worksheet(name)
                kinds = typed(d)
                for j, c in enumerate(d.columns):
                    ws.write_string(0, j, c, hdr)
                    w = max([len(c)] + [len(x) for x in d[c].head(200)])
                    ws.set_column(j, j, min(max(w + 2, 8), 50))
                vals = d.to_numpy()
                trunc = 0
                for ri, row in enumerate(vals, 1):
                    for j, v in enumerate(row):
                        if v == '':
                            continue
                        if kinds[j] == 'num':
                            x = float(v)
                            if x != x or x in (float('inf'), float('-inf')):
                                ws.write_string(ri, j, v)
                            else:
                                ws.write_number(ri, j, x)
                        else:
                            if len(v) > MAXSTR:
                                v = v[:MAXSTR]; trunc += 1
                            ws.write_string(ri, j, v)
                ws.freeze_panes(1, 0)
                ws.autofilter(0, 0, len(d), max(len(d.columns) - 1, 0))
                if trunc:
                    issues.append(f'truncated {trunc} cells {name}')
                tabs.append(dict(tab=name, sheet=s['sheet'], path=r['path'], role=r['role'], v15=r['old_unit'],
                                 rows_src=n_src, rows_tab=len(d), filter=fdesc, sha256=cur,
                                 columns=', '.join(d.columns)))
                log.append(dict(table=t, **tabs[-1]))
                print(t, name, n_src, len(d), flush=True)
        # README
        rd.set_column(0, 0, 26); rd.set_column(1, 1, 46); rd.set_column(2, 9, 18)
        rd.write(0, 0, f'Table {t}. {title}', h1)
        rd.merge_range(1, 0, 1, 7, GENERAL, wrap); rd.set_row(1, 45)
        rd.write(2, 0, f'Generated {today} by notes/2026-10-06-supplement-tables/build_workbooks.py')
        row = 4
        rd.write(row, 0, 'Sheets', b); row += 1
        head = ['Sheet', 'Title', 'Supplementary figures', 'Cited in', 'Data tabs', 'Previous table (v15)', 'Note']
        for j, x in enumerate(head):
            rd.write(row, j, x, hdr)
        for s in tplan:
            row += 1
            ntab = sum(1 for x in tabs if x['sheet'] == s['sheet'])
            cited = s['cited_by'] if not s['cited_by'].startswith('Methods') or t != 'S8' else s['cited_by']
            for j, x in enumerate([s['sheet'], s['sheet_title'], s['supp_figures'], s['cited_by'], ntab,
                                   s['v15_table'], s['note']]):
                rd.write(row, j, x, wrap if j in (1, 6) else None)
        row += 2
        rd.write(row, 0, 'Data tabs', b); row += 1
        head = ['Tab', 'Source file', 'Rows in tab', 'Rows in source', 'Row filter', 'SHA-256 of source', 'Columns']
        for j, x in enumerate(head):
            rd.write(row, j, x, hdr)
        for x in tabs:
            row += 1
            for j, v in enumerate([x['tab'], x['path'], x['rows_tab'], x['rows_src'], x['filter'], x['sha256'],
                                   x['columns']]):
                if isinstance(v, int):
                    rd.write_number(row, j, v)
                else:
                    rd.write_string(row, j, v)
        rd.activate()
        wb.close()
    with open(HERE / 'workbook-tabs.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(log[0]), delimiter='\t')
        w.writeheader(); w.writerows(log)
    json.dump({'tabs': len(log), 'issues': issues}, open(HERE / 'workbook-build.json', 'w'), indent=1)
    print('issues:', issues)
    return 1 if issues else 0


if __name__ == '__main__':
    sys.exit(main())
