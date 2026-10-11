#!/usr/bin/env python3
"""Independent checks of out/Supplementary_Tables_S1-S10.v11.xlsx against its source workbook and source TSVs. Writes out/...verify.json."""
import csv, html, json, math, re, sys, zipfile
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = next(q for q in [HERE, *HERE.parents] if (q / 'notes').is_dir())
sys.path.insert(0, str(HERE)); import spec_v11 as spec
SRC = ROOT / 'notes/2026-10-10-table-s9f/Supplementary_Tables_S1-S9.xlsx'
OUT = HERE.parent / 'out' / 'Supplementary_Tables_S1-S10.v11.xlsx'
problems, checks = [], {}
def need(ok, msg):
    if not ok: problems.append(msg)
A = re.compile(r'<c r="A(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', re.S)
CELL = re.compile(r'<c r="([A-Z]+)(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', re.S)
ROW = re.compile(r'<row r="(\d+)"[^>]*?(?:/>|>.*?</row>)', re.S)
LAB = re.compile(r'^(S\d+[a-z])(\.(\d+))?\s{2}')
COUNT = re.compile(r'(\d[\d,]*) rows')
def book(path):
    z = zipfile.ZipFile(path)
    wb = z.read('xl/workbook.xml').decode(); rels = z.read('xl/_rels/workbook.xml.rels').decode()
    R = {a['Id']: a['Target'] for a in (dict(re.findall(r'(\w+)="([^"]*)"', m)) for m in re.findall(r'<Relationship\b([^>]*?)/?>', rels))}
    sheets = [dict(re.findall(r'([\w:]+)="([^"]*)"', m)) for m in re.findall(r'<sheet\b([^>]*?)/>', wb)]
    path = {s['name']: 'xl/' + R[s['r:id']].lstrip('/').removeprefix('xl/') for s in sheets}
    SI = re.findall(r'<si>(.*?)</si>', z.read('xl/sharedStrings.xml').decode(), re.S)
    return z, [s['name'] for s in sheets], path, SI
zo, oorder, opath, OSI = book(SRC); zn, norder, npath, NSI = book(OUT)
def val(SI, attrs, inner):
    inner = inner or ''; t = re.search(r'\bt="(\w+)"', attrs); t = t.group(1) if t else 'n'
    if t == 's': return html.unescape(''.join(re.findall(r'<t[^>]*>(.*?)</t>', SI[int(re.search(r'<v>(\d+)</v>', inner).group(1))], re.S)))
    if t == 'inlineStr': return html.unescape(''.join(re.findall(r'<t[^>]*>(.*?)</t>', inner, re.S)))
    v = re.search(r'<v>(.*?)</v>', inner, re.S)
    if not v: return None
    return html.unescape(v.group(1)) if t in ('str', 'e') else float(v.group(1))
def rows_of(x): return [(int(m.group(1)), m.group(0)) for m in ROW.finditer(x.split('<sheetData>', 1)[1].split('</sheetData>', 1)[0])]
def labels_in(x, SI, name):
    out = {}
    for m in A.finditer(x):
        if int(m.group(1)) >= 5 and re.search(r'\bt="(s|inlineStr|str)"', m.group(2)):
            v = val(SI, m.group(2), m.group(3)); mm = LAB.match(v or '')
            if mm and mm.group(1) == name: out[int(m.group(1))] = v
    return out
NORM = lambda rx: re.sub(r' r="([A-Z]*)\d+"', r' r="\1#"', rx)
def strip(rx, rows): return A.sub(lambda m: '' if int(m.group(1)) in rows else m.group(0), rx)
def relink(x): return re.sub(r'location="(?:\'|&apos;)([^\'&]+)(?:\'|&apos;)!', lambda m: f'location="\'{spec.MAP.get(m.group(1), m.group(1))}\'!', x)
def tables_of(z, p):
    rp = p.replace('worksheets/', 'worksheets/_rels/') + '.rels'
    if rp not in z.namelist(): return []
    return [z.read('xl/tables/' + Path(t).name).decode() for t in re.findall(r'Target="([^"]+)"', z.read(rp).decode()) if 'table' in t]

# 1. order and dropped sheets
need(norder == spec.ORDER, 'sheet order differs from spec_v11.ORDER'); checks['sheets'] = len(norder)
need(set(oorder) == set(spec.MAP) | set(spec.DROP_OLD) | {'Contents', 'Sources'}, 'source sheets not fully accounted for')

# 2. old worksheets: kept rows identical apart from row numbers and the rewritten cells (A1, A2, block labels)
kept_rows = 0
for old, new in spec.MAP.items():
    xo = zo.read(opath[old]).decode(); xn = zn.read(npath[new]).decode()
    lo = labels_in(xo, OSI, old); ln = labels_in(xn, NSI, new)
    idx = {r: int(LAB.match(v).group(3) or 1) for r, v in lo.items()}
    starts = sorted(lo); drop = set(spec.DROP_BLOCKS.get(old, []))
    spans = [(a, starts[k + 1] - 1 if k + 1 < len(starts) else 10 ** 9) for k, a in enumerate(starts) if idx[a] in drop]
    app = spec.APPEND.get(new, [])
    cut = sorted(ln)[len(lo) - len(drop)] if app else 10 ** 9            # first appended label row
    ro = [NORM(strip(rx, {1, 2} | set(lo))) for r, rx in rows_of(xo) if not any(a <= r <= b for a, b in spans)]
    rn = [NORM(strip(rx, {1, 2} | set(ln))) for r, rx in rows_of(xn) if r < cut]
    need(ro == rn, f'{old}->{new}: kept rows differ'); kept_rows += len(rn)
    need(len(ln) == len(lo) - len(drop) + len(app), f'{new}: label count')
    single = len(ln) == 1
    need(all((LAB.match(v).group(2) is None) == single for v in ln.values()), f'{new}: label suffix form')
    if not drop:
        so, sn = xo.split('<sheetData>')[0] + xo.split('</sheetData>')[1], xn.split('<sheetData>')[0] + xn.split('</sheetData>')[1]
        so = relink(so)
        if app:                                                          # appended blocks change only these three elements
            def outside(x, k, appended=False):
                x = re.sub(r'<dimension ref="[^"]*"/>', '', x)
                x = re.sub(r'<tableParts count="\d+">((?:<tablePart r:id="rId\d+"/>)*)</tableParts>',
                           lambda m: ''.join(re.findall(r'<tablePart r:id="rId\d+"/>', m.group(1))[:k]), x)
                return re.sub(r'<col min="(\d+)" max="\d+" width="18" customWidth="1"/></cols>', '</cols>', x) if appended else x
            k = len(re.findall(r'<tablePart ', so)); need(outside(so, k) == outside(sn, k, True), f'{new}: parts outside sheetData differ')
        else: need(so == sn, f'{new}: parts outside sheetData differ')
    t1 = {r: val(NSI, a, i) for r, a, i in ((int(m.group(1)), m.group(2), m.group(3)) for m in A.finditer(xn)) if r in (1, 2)}
    need(t1[1] == f'Table {new}. ' + spec.TITLE_OVERRIDE.get(new, t1[1].split('. ', 1)[1]), f'{new}: title')
    need(t1[2] == spec.caption(new), f'{new}: caption')
checks['old_sheets_rows_identical'] = kept_rows

# 3. Excel tables (all sheets): unique names, header cells equal the column names; conditional formats point at table headers
names, cf_n = [], 0
for s in spec.ORDER:
    x = zn.read(npath[s]).decode(); rows = dict(rows_of(x)); heads = set()
    for tx in tables_of(zn, npath[s]):
        nm = re.search(r'displayName="([^"]+)"', tx).group(1); names.append(nm)
        need(re.fullmatch(rf'T_{s}_\d+', nm) is not None or (s, nm) == ('Sources', 'T_Sources'), f'table {nm} on {s}')
        c1, h, c2, last = re.search(r'\sref="([A-Z]+)(\d+):([A-Z]+)(\d+)"', tx).groups(); h = int(h); heads.add(h)
        cols = [html.unescape(c) for c in re.findall(r'<tableColumn [^>]*name="([^"]*)"', tx)]
        cells = {c: val(NSI, a, i) for c, _, a, i in CELL.findall(rows.get(h, ''))}
        need([cells.get(chr(65 + k) if k < 26 else 'A' + chr(39 + k)) for k in range(len(cols))] == cols or len(cols) > 26, f'{nm}: header cells')
        need(int(last) in rows, f'{nm}: last row {last} missing')
    for el in re.findall(r'<conditionalFormatting\b.*?</conditionalFormatting>', x, re.S):
        s0 = int(re.search(r'sqref="[A-Z]+(\d+)', el).group(1)); H = int(re.search(r'ROW\(\)-(\d+)', el).group(1)); cf_n += 1
        need(H == s0 - 1 and H in heads, f'{s}: conditional format {s0}/{H}')
    for loc in re.findall(r'location="\'([^\']+)\'!', x): need(loc in spec.ORDER, f'{s}: link to {loc}')
need(len(names) == len(set(names)), 'duplicate table names'); checks['excel_tables'] = len(names); checks['conditional_formats'] = cf_n

# 4. new worksheets: every cell equals its source TSV value
cells_checked = 0
for key, sh in spec.NEW.items():
    x = zn.read(npath[key]).decode(); grid = {(c, int(r)): (a, i) for c, r, a, i in CELL.findall(x)}
    labs = sorted(labels_in(x, NSI, key)); need(len(labs) == len(sh['blocks']), f'{key}: labels')
    for b, lr in zip(sh['blocks'], labs):
        with open(ROOT / b['f'], newline='', encoding='utf-8') as fh: data = list(csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE))
        h = lr + 3
        for k, rec in enumerate(data):
            for j, src in enumerate(rec):
                L, n = '', j + 1
                while n: n, rr = divmod(n - 1, 26); L = chr(65 + rr) + L
                a, i = grid.get((L, h + k), ('', None)); got = val(NSI, a, i); cells_checked += 1
                ok = got == src if k == 0 else got == 'NA' if src in ('', 'NA', 'nan', 'NaN') else (got == float(src) if isinstance(got, float) else got == src)
                if not ok: problems.append(f'{key} {L}{h+k}: {got!r} != {src!r}'); break
checks['new_cells_equal_to_source'] = cells_checked

# 4b. appended blocks of carried-over worksheets: every cell equals its source TSV value; empty source fields stay empty
app_cells = 0
for new, blocks in spec.APPEND.items():
    x = zn.read(npath[new]).decode(); grid = {(c, int(r)): (a, i) for c, r, a, i in CELL.findall(x)}
    labs = sorted(labels_in(x, NSI, new))[-len(blocks):]
    for b, lr in zip(blocks, labs):
        v = val(NSI, *grid[('A', lr)]); need(v.endswith(f"  {b['name']}    {sum(1 for _ in open(ROOT / b['f'])) - 1:,} rows"), f'{new}: label {v}')
        need(val(NSI, *grid[('A', lr + 1)]) == b['note'], f'{new}: note row {lr + 1}')
        with open(ROOT / b['f'], newline='', encoding='utf-8') as fh: data = list(csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE))
        h = lr + 2
        for k, rec in enumerate(data):
            for j, src in enumerate(rec):
                L, n = '', j + 1
                while n: n, rr = divmod(n - 1, 26); L = chr(65 + rr) + L
                a, i = grid.get((L, h + k), ('', None)); got = val(NSI, a, i); app_cells += 1
                ok = got == src if k == 0 else got is None if src == '' else (got == float(src) if isinstance(got, float) else got == src)
                if not ok: problems.append(f'{new} {L}{h+k}: {got!r} != {src!r}'); break
        need((L, h + len(data) - 1) in grid and (('A', h + len(data)) not in grid), f'{new}: block {lr} extent')
checks['appended_cells_equal_to_source'] = app_cells

# 5. Contents and Sources against the worksheet labels
ct = zn.read(npath['Contents']).decode(); crow = dict(rows_of(ct))
links = re.findall(r'<hyperlink ref="A(\d+)" location="\'([^\']+)\'!A1"', ct)
need([s for _, s in links] == spec.ORDER[1:-1], 'Contents links/order')
all_labels = []
for r, s in links:
    x = zn.read(npath[s]).decode(); labs = labels_in(x, NSI, s).values()
    ids = [LAB.match(v).group(0).strip() for v in labs]; all_labels += ids
    c = {cc: val(NSI, a, i) for cc, _, a, i in CELL.findall(crow[int(r)])}
    need(c.get('D') == len(ids) and c.get('E') == sum(int(COUNT.search(v).group(1).replace(',', '')) for v in labs), f'Contents counts {s}')
    need(('†' in c['B']) == (s in set(spec.NEW) | set(spec.POSTHOC_OLD)), f'Contents dagger {s}')
sr = zn.read(npath['Sources']).decode()
ids = [v for _, _, v in re.findall(r'<hyperlink ref="A(\d+)" location="\'([^\']+)\'!A1" display="([^"]+)"', sr) if v != '← Contents']
need(sorted(ids) == sorted(all_labels), f'Sources blocks differ from labels: {sorted(set(ids) ^ set(all_labels))[:6]}')
need(ids == sorted(ids, key=lambda b: (spec.ORDER.index(re.match(r'S\d+[a-z]', b).group(0)), int((b.split('.') + ['1'])[1]))), 'Sources order')
checks['blocks'] = len(all_labels)

# 6. stale references in any text that is still referenced
text = zn.read('xl/sharedStrings.xml').decode() + ''.join(zn.read(npath[s]).decode() for s in ['Contents', 'Sources'] + list(spec.NEW) + [v for k, v in spec.MAP.items() if k.startswith('S9')])
for pat in ['Figure 2e panel', 'Figure 5f',   # Figure 2e is now the composition panel (2026-10-11 relettering); the stale text was the route hold-out
            'Figure S11c', 'Figure S20h', 'Tables S6–S8', 'Table S9. Analyses', 'Supplementary Tables S1–S9']:
    need(pat not in text, f'stale text: {pat}')
OUTV = HERE.parent / 'out' / 'Supplementary_Tables_S1-S10.v11.verify.json'
OUTV.write_text(json.dumps(dict(checks=checks, problems=problems), indent=1, ensure_ascii=False))
print(json.dumps(checks)); print('problems', len(problems)); [print(' ', p) for p in problems[:20]]
