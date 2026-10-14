"""Verify the S4j candidate against the production workbook and the frozen source tables. Writes verify.json."""
import json, math, os, re, sys, zipfile
from pathlib import Path
import pandas as pd, openpyxl
ROOT = Path(os.environ['REPRO_ROOT']); OUT = Path(os.environ['S4J_OUT'])
OLD = ROOT / 'table/Supplementary_Tables_S1-S10.xlsx'; NEW = OUT / 'Supplementary_Tables_S1-S10.xlsx'
build = json.load(open(OUT / 'build.json')); checks = {}; fail = []
def ok(name, cond, detail=None):
    checks[name] = dict(pass_=bool(cond), detail=detail)
    if not cond: fail.append(name)

zo, zn = zipfile.ZipFile(OLD), zipfile.ZipFile(NEW)
po, pn = set(zo.namelist()), set(zn.namelist())
added = sorted(pn - po); ok('no parts removed', po <= pn, sorted(po - pn))
expected_added = sorted([build['sheet_part'], re.sub(r'worksheets/(sheet\d+\.xml)$', r'worksheets/_rels/\1.rels', build['sheet_part'])] + build['tables'])
ok('added parts', added == expected_added, added)
changed = sorted(p for p in po if zo.read(p) != zn.read(p))
expected_changed = sorted(['[Content_Types].xml', 'docProps/app.xml', 'xl/workbook.xml', 'xl/_rels/workbook.xml.rels',
                           'xl/worksheets/sheet1.xml', 'xl/worksheets/sheet45.xml', 'xl/tables/table152.xml'])
ok('only expected existing parts changed', changed == expected_changed, changed)
ok('styles.xml unchanged', zo.read('xl/styles.xml') == zn.read('xl/styles.xml'))

wo = openpyxl.load_workbook(OLD, read_only=True, data_only=True); wn = openpyxl.load_workbook(NEW, read_only=True, data_only=True)
so, sn = wo.sheetnames, wn.sheetnames
i = so.index('S4i'); ok('sheet order', sn == so[:i+1] + ['S4j'] + so[i+1:], sn[i-1:i+3])
vals = lambda ws: [tuple(c for c in r) for r in ws.iter_rows(values_only=True)]

# S4j blocks equal the source tables cell by cell
def norm(v):
    if v is None or (isinstance(v, float) and math.isnan(v)) or v is pd.NA: return 'NA'
    if type(v).__name__ in ('bool', 'bool_'): return str(bool(v))
    if isinstance(v, float) or type(v).__name__.startswith('float'): return float(v)
    if type(v).__name__.startswith(('int', 'uint')): return int(v)
    return str(v)
grid = vals(wn['S4j'])
src1 = pd.read_csv(ROOT / build['blocks'][0]['file'], sep='\t')
cols1 = ['human_gene_id', 'mouse_gene_id', 'human_gene_symbol', 'mouse_gene_symbol', 'P_S_measurable', 'P_S_human_direction',
         'P_S_mouse_direction', 'P_M_measurable', 'P_M_human_direction', 'P_M_mouse_direction', 'case_first_three', 'threshold']
b1 = src1[src1['case_all_four'] == True][cols1].reset_index(drop=True)
b2 = pd.read_csv(ROOT / build['blocks'][1]['file'], sep='\t')
def find(title): return next(k for k, r in enumerate(grid) if r and isinstance(r[0], str) and r[0].startswith(title))
for tag, df in [('S4j.1', b1), ('S4j.2', b2)]:
    k = find(tag); head = [c for c in grid[k+3] if c is not None]
    ok(f'{tag} header', head == list(df.columns), head)
    body = [[c for c in r[:len(df.columns)]] for r in grid[k+4:k+4+len(df)]]
    exp = [[norm(v) for v in row] for row in df.itertuples(index=False)]
    got = [[(float(c) if isinstance(c, (int, float)) and not isinstance(c, bool) else c) for c in r] for r in body]
    expf = [[(float(c) if isinstance(c, (int, float)) else c) for c in r] for r in exp]
    ok(f'{tag} values equal source ({len(df)} rows)', got == expf)
    ok(f'{tag} title row count', grid[k][0].endswith(f'({len(df)} rows)'), grid[k][0])
n48 = b1[b1.threshold == 0.02]
cells = {k: int(((n48[a] == x) & (n48[b] == x)).sum()) for k, (a, b, x) in dict(
    PS_prox=('P_S_human_direction', 'P_S_mouse_direction', 'proximal'), PS_dist=('P_S_human_direction', 'P_S_mouse_direction', 'distal'),
    PM_prox=('P_M_human_direction', 'P_M_mouse_direction', 'proximal'), PM_dist=('P_M_human_direction', 'P_M_mouse_direction', 'distal')).items()}
ok('48 pairs and cells 7/14/27/3 at 0.02', len(n48) == 48 and cells == dict(PS_prox=7, PS_dist=14, PM_prox=27, PM_dist=3), cells)
ok('12 pairs at 0.05', int((b1.threshold == 0.05).sum()) == 12)
ok('main48 ids identical in both blocks', set(zip(n48.human_gene_id, n48.mouse_gene_id)) == set(zip(b2[b2.set_name == 'main48'].human_gene_id, b2[b2.set_name == 'main48'].mouse_gene_id)))
ok('CSTF3, DDX6, PAIP2 in main48', {'CSTF3', 'DDX6', 'PAIP2'} <= set(n48.human_gene_symbol))

# Contents and Sources: old rows reappear unchanged after the inserted rows
for name, key, k in [('Contents', 'S4i', 1), ('Sources', None, 2)]:
    o, n = vals(wo[name]), vals(wn[name])
    a = next(j for j, r in enumerate(o) if r and r[0] == 'S4i') if key else max(j for j, r in enumerate(o) if r and isinstance(r[0], str) and r[0].startswith('S4i.'))
    ok(f'{name} rows preserved around insertion', o[:a+1] == n[:a+1] and o[a+1:] == n[a+1+k:], (len(o), len(n)))
    ok(f'{name} inserted rows', [r[0] for r in n[a+1:a+1+k]] == (['S4j'] if k == 1 else ['S4j.1', 'S4j.2']), n[a+1:a+1+k])
c = vals(wn['Contents']); row = next(r for r in c if r and r[0] == 'S4j')
ok('Contents S4j row', row[:5] == ('S4j', 'Orthologs with same-direction responses in both species †', 'Figure 3a, S19c', 2, 187), row)

# hyperlinks: every link targets the sheet named in its row
def links(z, part, sheetxml):
    x = z.read(part).decode(); rows = {}
    for m in re.finditer(r'<c r="A(\d+)"[^>]*t="inlineStr"><is><t[^>]*>([^<]*)</t>', x): rows[int(m.group(1))] = m.group(2)
    return [(int(r), loc, rows.get(int(r))) for r, loc in re.findall(r'<hyperlink ref="A(\d+)" location="([^"]+)"', x)]
bad = [(r, l, v) for r, l, v in links(zn, 'xl/worksheets/sheet1.xml', None) if v and re.match(r'^S\d+[a-z]$', v) and l != f"'{v}'!A1"]
ok('Contents links point to their own sheet', not bad, bad[:5])
bad = [(r, l, v) for r, l, v in links(zn, 'xl/worksheets/sheet45.xml', None) if v and re.match(r'^S\d+[a-z]\.', v) and l != f"'{v.split('.')[0]}'!A1"]
ok('Sources links point to their block sheet', not bad, bad[:5])
nl_old = len(re.findall(r'<hyperlink ', zo.read('xl/worksheets/sheet45.xml').decode())); nl_new = len(re.findall(r'<hyperlink ', zn.read('xl/worksheets/sheet45.xml').decode()))
ok('Sources gained two links', nl_new == nl_old + 2, (nl_old, nl_new))
t = zn.read('xl/tables/table152.xml').decode(); nrows = len(vals(wn['Sources']))
ok('Sources table range covers new rows', 'ref="A5:F207"' in t and 'autoFilter ref="A5:F207"' in t, re.findall(r'ref="[^"]+"', t))
names = [re.search(r'displayName="([^"]+)"', zn.read(p).decode()).group(1) for p in pn if re.match(r'xl/tables/table\d+\.xml$', p)]
ok('table display names unique', len(names) == len(set(names)), len(names))
ids = [re.search(r'\bid="(\d+)"', zn.read(p).decode()).group(1) for p in pn if re.match(r'xl/tables/table\d+\.xml$', p)]
ok('table ids unique', len(ids) == len(set(ids)), len(ids))
app = zn.read('docProps/app.xml').decode()
ok('docProps sheet list', '<ns1:i4>60</ns1:i4>' in app and 'size="60"' in app and '<ns1:lpstr>S4i</ns1:lpstr><ns1:lpstr>S4j</ns1:lpstr>' in app)
for s in so:
    if s in ('Contents', 'Sources'): continue
    if vals(wo[s]) != vals(wn[s]): fail.append(f'values changed in {s}')
ok('values of all other sheets unchanged', not [f for f in fail if f.startswith('values changed')])
json.dump(dict(old_sha256=build['source_sha256'], new_sha256=build['output_sha256'], checks=checks, failures=fail), open(OUT / 'verify.json', 'w'), indent=1, ensure_ascii=False)
print('FAIL' if fail else 'PASS', fail)
for k, v in checks.items(): print(('ok  ' if v['pass_'] else 'FAIL'), k)
