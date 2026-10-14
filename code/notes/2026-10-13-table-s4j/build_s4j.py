"""Table S4j: orthologs with same-direction responses in both species (the main48 set) and the
functional-category membership behind Figure S19c.

Inserted at the XML level into a copy of the production workbook, after S4i. Every existing sheet part
stays byte-identical except Contents (one inserted row), Sources (two inserted rows and the table range),
workbook.xml, its rels, [Content_Types].xml and docProps/app.xml. Styles are the existing cell styles of
the S4 sheets (43 title, 44 caption, 16 link, 46 block, 47 note, 48 header, 49 text, 50 float, 51 int),
so styles.xml is unchanged. Values are copied from the frozen result tables without recomputation.

Env: REPRO_ROOT (project root holding the inputs), S4J_OUT (output directory).
"""
import hashlib, json, math, os, re, zipfile
from pathlib import Path
from xml.sax.saxutils import escape
import pandas as pd

ROOT = Path(os.environ['REPRO_ROOT']); OUT = Path(os.environ['S4J_OUT']); OUT.mkdir(parents=True, exist_ok=True)
SRC = ROOT / 'table/Supplementary_Tables_S1-S10.xlsx'
DST = OUT / 'Supplementary_Tables_S1-S10.xlsx'
D = 'notes/2026-10-10-conserved-response/function/'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
NAME = 'S4j'; ANCHOR = 'S4i'
TITLE = 'Orthologs with same-direction responses in both species'
CAPTION = ('Part of Table S4: Cross-species comparisons and cis-feature models  ·  Specified after the initial results  ·  '
           'Related figures: Figure 3a, S19c')
COLS1 = ['human_gene_id', 'mouse_gene_id', 'human_gene_symbol', 'mouse_gene_symbol', 'P_S_measurable', 'P_S_human_direction',
         'P_S_mouse_direction', 'P_M_measurable', 'P_M_human_direction', 'P_M_mouse_direction', 'case_first_three', 'threshold']
BLOCKS = [
    dict(title='Orthologs with same-direction responses in both species', file='all_measured_ortholog_classification.tsv.gz',
         filt=lambda d: d[d['case_all_four'] == True][COLS1],
         rule='case_all_four == True; columns human_gene_id to P_M_mouse_direction, case_first_three and threshold; no statistical recomputation.',
         note=('One-to-one high-confidence orthologs (16,092 pairs) whose recurrent responses agreed in direction in both species within one '
               'contrast: proximal in both or distal in both, from source to P (P_S) or in P relative to M (P_M). At |ΔDWUI| ≥ 0.02 the four '
               'cells hold 7 (P–S proximal), 14 (P–S distal), 27 (P–M proximal) and 3 (P–M distal) pairs, 48 unique pairs; at 0.05 they hold '
               '12 unique pairs. Directions are within-species recurrence calls (proximal, distal or nonrecurrent); NA, contrast not measurable '
               'in that species. case_first_three marks the 45 pairs in the first three cells, a sensitivity set. Each ortholog pair counts once '
               'however many cells it occupies.')),
    dict(title='Functional-category membership of the shared sets', file='primary_case_categories.tsv',
         filt=lambda d: d, rule='All rows; no statistical recomputation.',
         note=('main48, the 48 pairs of S4j.1 at |ΔDWUI| ≥ 0.02. hN79, the 79 unique ortholog pairs whose equal-weighted proximal response in '
               'the two naive human conditions relative to primed coincided with mouse P–M proximal recurrence. Human GO annotation. _known, '
               'the gene is annotated in the category namespace; _member, the gene belongs to the category: GO:0003723 RNA binding, '
               'GO:0006397 mRNA processing, GO:0006417 regulation of translation. Enrichment tests of these categories are in Table S7f and '
               'Figure S19c.')),
]

# ------------------------------------------------------------------ helpers
def col(j):
    s = ''
    while j: j, r = divmod(j - 1, 26); s = chr(65 + r) + s
    return s
def tcell(ref, text, style):
    return f'<c r="{ref}" s="{style}" t="inlineStr"><is><t xml:space="preserve">{escape(str(text))}</t></is></c>'
def ncell(ref, value, style):
    return f'<c r="{ref}" s="{style}" t="n"><v>{value}</v></c>'
def value_cell(ref, v):
    if v is None or (isinstance(v, float) and math.isnan(v)) or v is pd.NA: return tcell(ref, 'NA', 49)
    tn = type(v).__name__
    if tn in ('bool', 'bool_'): return tcell(ref, str(bool(v)), 49)
    if tn.startswith(('int', 'uint')): return ncell(ref, int(v), 51)
    if isinstance(v, float) or tn.startswith('float'): return ncell(ref, repr(float(v)), 50)
    return tcell(ref, v, 49)

# ------------------------------------------------------------------ 1. build the sheet and its tables
zin = zipfile.ZipFile(SRC)
read = lambda p: zin.read(p).decode('utf-8')
wbx = read('xl/workbook.xml'); wrels = read('xl/_rels/workbook.xml.rels')
rels = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', wrels)) | {a: b for b, a in re.findall(r'Target="([^"]+)"[^>]*Id="(rId\d+)"', wrels)}
sheet_rid = {n: r for n, r in re.findall(r'<sheet name="([^"]+)" sheetId="\d+" r:id="(rId\d+)"', wbx)}
path_of = lambda n: 'xl/' + rels[sheet_rid[n]]
anchor_xml = read(path_of(ANCHOR))
grab = lambda tag: re.search(r'<%s[ >].*?</%s>|<%s[^>]*/>' % (tag, tag, tag), anchor_xml, re.S).group(0)

next_sheet_id = max(int(x) for x in re.findall(r'sheetId="(\d+)"', wbx)) + 1
next_rid = max(int(x) for x in re.findall(r'Id="rId(\d+)"', wrels)) + 1
names = zin.namelist()
next_sheet_file = max(int(re.search(r'sheet(\d+)\.xml$', p).group(1)) for p in names if re.match(r'xl/worksheets/sheet\d+\.xml$', p)) + 1
table_ids = [int(re.search(r'\bid="(\d+)"', read(p)).group(1)) for p in names if re.match(r'xl/tables/table\d+\.xml$', p)]
table_files = [int(re.search(r'table(\d+)\.xml$', p).group(1)) for p in names if re.match(r'xl/tables/table\d+\.xml$', p)]
next_table = max(table_ids + table_files) + 1

rows = [f'<row r="1" ht="20" customHeight="1">{tcell("A1", f"Table {NAME}. {TITLE}", 43)}</row>',
        f'<row r="2" ht="20" customHeight="1">{tcell("A2", CAPTION, 44)}</row>',
        f'<row r="3" ht="20" customHeight="1">{tcell("A3", "← Contents", 16)}</row>']
tables = []; sources = []; widths = {}; r = 5; maxcol = 1
for i, b in enumerate(BLOCKS, 1):
    path = ROOT / D / b['file']; full = pd.read_csv(path, sep='\t'); df = b['filt'](full).reset_index(drop=True)
    label = '%s.%d  %s (%d rows)' % (NAME, i, b['title'], len(df))
    rows.append('<row r="%d">' % r + tcell('A%d' % r, label, 46) + '</row>')
    rows.append(f'<row r="{r+1}">{tcell(f"A{r+1}", b["note"], 47)}</row>')
    rows.append(f'<row r="{r+2}">{tcell(f"A{r+2}", "Rows shown: " + b["rule"], 47)}</row>')
    h = r + 3
    rows.append(f'<row r="{h}" ht="42" customHeight="1">' + ''.join(tcell(f'{col(j)}{h}', c, 48) for j, c in enumerate(df.columns, 1)) + '</row>')
    for k, rec in enumerate(df.itertuples(index=False), h + 1):
        rows.append(f'<row r="{k}">' + ''.join(value_cell(f'{col(j)}{k}', v) for j, v in enumerate(rec, 1)) + '</row>')
    last = h + len(df); ncol = len(df.columns); maxcol = max(maxcol, ncol)
    ref = f'A{h}:{col(ncol)}{last}'
    tid = next_table + i - 1
    tcols = ''.join(f'<tableColumn id="{j}" name="{escape(c)}"/>' for j, c in enumerate(df.columns, 1))
    tables.append((f'xl/tables/table{tid}.xml',
                   f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<table xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                   f'id="{tid}" name="T_{NAME}_{i}" displayName="T_{NAME}_{i}" ref="{ref}" headerRowCount="1"><autoFilter ref="{ref}"/>'
                   f'<tableColumns count="{ncol}">{tcols}</tableColumns><tableStyleInfo showRowStripes="0" showColumnStripes="0"/></table>'))
    for j, c in enumerate(df.columns, 1):
        sample = [str(x) for x in df.iloc[:300, j - 1]]
        widths[j] = max(widths.get(j, 7), min(40, max(len(c) * 1.1, max([0] + [len(x) for x in sample])) + 3))
    sources.append(dict(block=f'{NAME}.{i}', file=D + b['file'], shown=len(df), total=len(full), rule=b['rule'] + ' ' + b['note'], sha=sha(path)))
    r = last + 3
widths[1] = max(widths.get(1, 7), 28.0)
cols = '<cols>' + ''.join(f'<col min="{j}" max="{j}" width="{w:.1f}" customWidth="1"/>' for j, w in sorted(widths.items())) + '</cols>'
lastrow = r - 3
sheet_xml = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
             '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
             'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
             + grab('sheetPr') + f'<dimension ref="A1:{col(maxcol)}{lastrow}"/>' + grab('sheetViews') + grab('sheetFormatPr') + cols
             + '<sheetData>' + ''.join(rows) + '</sheetData>' + grab('hyperlinks') + grab('pageMargins')
             + f'<tableParts count="{len(tables)}">' + ''.join(f'<tablePart r:id="rId{j}"/>' for j in range(1, len(tables) + 1)) + '</tableParts></worksheet>')
sheet_path = f'xl/worksheets/sheet{next_sheet_file}.xml'
sheet_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
              + ''.join(f'<Relationship Id="rId{j}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/table" '
                        f'Target="../tables/{Path(p).name}"/>' for j, (p, _) in enumerate(tables, 1)) + '</Relationships>')

# ------------------------------------------------------------------ 2. package-level edits (string level)
parts = {}
anchor_sheet = re.search(r'<sheet name="%s"[^>]*/>' % ANCHOR, wbx).group(0)
parts['xl/workbook.xml'] = wbx.replace(anchor_sheet, anchor_sheet + f'<sheet name="{NAME}" sheetId="{next_sheet_id}" r:id="rId{next_rid}" />', 1)
parts['xl/_rels/workbook.xml.rels'] = wrels.replace('</Relationships>',
    f'<Relationship Id="rId{next_rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
    f'Target="worksheets/sheet{next_sheet_file}.xml" /></Relationships>')
ct = read('[Content_Types].xml')
add = f'<Override PartName="/{sheet_path}" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml" />'
add += ''.join(f'<Override PartName="/{p}" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.table+xml" />' for p, _ in tables)
parts['[Content_Types].xml'] = ct.replace('</Types>', add + '</Types>')
app = read('docProps/app.xml'); n_old = len(re.findall(r'<sheet ', wbx))
assert f'<ns1:i4>{n_old}</ns1:i4>' in app and f'size="{n_old}" baseType="lpstr"' in app
app = app.replace(f'<ns1:i4>{n_old}</ns1:i4>', f'<ns1:i4>{n_old+1}</ns1:i4>').replace(f'size="{n_old}" baseType="lpstr"', f'size="{n_old+1}" baseType="lpstr"')
app = app.replace(f'<ns1:lpstr>{ANCHOR}</ns1:lpstr>', f'<ns1:lpstr>{ANCHOR}</ns1:lpstr><ns1:lpstr>{NAME}</ns1:lpstr>', 1)
parts['docProps/app.xml'] = app

def insert_rows(xml, after, new_rows, ncols_letter):
    """Insert rows after row `after`, shifting later rows, cell refs, hyperlinks and the dimension."""
    k = len(new_rows)
    def shift_row(m):
        rn = int(m.group(1) or m.group(2))
        if rn <= after: return m.group(0)
        body = re.sub(r'(<c r="[A-Z]+)(\d+)"', lambda c: f'{c.group(1)}{int(c.group(2))+k}"', m.group(0))
        return re.sub(r'^<row r="\d+"', f'<row r="{rn+k}"', body)
    xml = re.sub(r'<row r="(\d+)"[^>]*?/>|<row r="(\d+)"[^>]*>.*?</row>', shift_row, xml, flags=re.S)
    xml = re.sub(r'(<hyperlink ref="[A-Z]+)(\d+)"', lambda m: f'{m.group(1)}{int(m.group(2))+k if int(m.group(2))>after else int(m.group(2))}"', xml)
    pos = re.search(r'<row r="%d"[^>]*>.*?</row>' % after, xml, re.S).end()
    xml = xml[:pos] + ''.join(new_rows) + xml[pos:]
    xml = re.sub(r'<dimension ref="A1:([A-Z]+)(\d+)"', lambda m: f'<dimension ref="A1:{m.group(1)}{int(m.group(2))+k}"', xml)
    return xml
def row_styles(xml, rn):
    row = re.search(r'<row r="%d"[^>]*>.*?</row>' % rn, xml, re.S).group(0)
    return [int(s) for s in re.findall(r'<c r="[A-Z]+\d+" s="(\d+)"', row)]
def a_value(xml, rn):
    m = re.search(r'<row r="%d"[^>]*><c r="A%d"[^>]*t="inlineStr"><is><t[^>]*>([^<]*)</t>' % (rn, rn), xml)
    return m.group(1) if m else None

# Contents: one row after S4i, styled like the second sheet row of a table section (S4b)
cp = path_of('Contents'); cx = read(cp)
crow = {a_value(cx, n): n for n in range(1, 200) if a_value(cx, n)}
after = crow[ANCHOR]; st = row_styles(cx, crow['S4b']); nr = after + 1
vals = [NAME, TITLE + ' †', 'Figure 3a, S19c', len(BLOCKS), sum(s['shown'] for s in sources)]
newrow = (f'<row r="{nr}" spans="1:5">' + tcell(f'A{nr}', vals[0], st[0]) + tcell(f'B{nr}', vals[1], st[1]) + tcell(f'C{nr}', vals[2], st[2])
          + ncell(f'D{nr}', vals[3], st[3]) + ncell(f'E{nr}', vals[4], st[4]) + '</row>')
cx = insert_rows(cx, after, [newrow], 'E')
cx = cx.replace(f'<hyperlink ref="A{after}" ', f'<hyperlink ref="A{nr}" location="\'{NAME}\'!A1" display="{NAME}" /><hyperlink ref="A{after}" ', 1) \
     if f'<hyperlink ref="A{after}" ' in cx else cx
assert f'location="\'{NAME}\'!A1"' in cx, 'Contents hyperlink not inserted'
parts[cp] = cx

# Sources: rows after the last S4i block, styles alternating by row parity as in the existing sheet
sp = path_of('Sources'); sx = read(sp)
srow = {a_value(sx, n): n for n in range(1, 400) if a_value(sx, n)}
after = max(n for v, n in srow.items() if v.startswith(ANCHOR + '.'))
even = row_styles(sx, next(n for n in sorted(srow.values()) if n % 2 == 0 and n > 5))
odd = row_styles(sx, next(n for n in sorted(srow.values()) if n % 2 == 1 and n > 5))
new = []
for i, s in enumerate(sources):
    rn = after + 1 + i; stl = even if rn % 2 == 0 else odd
    new.append(f'<row r="{rn}">' + tcell(f'A{rn}', s['block'], stl[0]) + tcell(f'B{rn}', s['file'], stl[1]) + ncell(f'C{rn}', s['shown'], stl[2])
               + ncell(f'D{rn}', s['total'], stl[3]) + tcell(f'E{rn}', s['rule'], stl[4]) + tcell(f'F{rn}', s['sha'], stl[5]) + '</row>')
sx = insert_rows(sx, after, new, 'F')
links = ''.join(f'<hyperlink ref="A{after+1+i}" location="\'{NAME}\'!A1" display="{s["block"]}" />' for i, s in enumerate(sources))
m = re.search(r'<hyperlink ref="A%d" [^>]*/>' % after, sx)
sx = sx[:m.end()] + links + sx[m.end():]
parts[sp] = sx
srels = read(re.sub(r'worksheets/(sheet\d+\.xml)$', r'worksheets/_rels/\1.rels', sp))
tpath = 'xl/tables/' + Path(re.search(r'Target="([^"]+)"', srels).group(1)).name
tx = read(tpath)
tx = re.sub(r'ref="A5:F(\d+)"', lambda m: f'ref="A5:F{int(m.group(1))+len(sources)}"', tx)
parts[tpath] = tx

# ------------------------------------------------------------------ 3. write
parts[sheet_path] = sheet_xml
parts[f'xl/worksheets/_rels/sheet{next_sheet_file}.xml.rels'] = sheet_rels
for p, x in tables: parts[p] = x
with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as z:
    for info in zin.infolist():
        z.writestr(info, parts.pop(info.filename).encode('utf-8') if info.filename in parts else zin.read(info.filename))
    for p, x in parts.items(): z.writestr(p, x.encode('utf-8'))
json.dump(dict(source_workbook='table/Supplementary_Tables_S1-S10.xlsx', source_sha256=sha(SRC), output_sha256=sha(DST),
               sheet=NAME, after=ANCHOR, sheet_part=sheet_path, tables=[p for p, _ in tables], blocks=sources, production_replaced=False),
          open(OUT / 'build.json', 'w'), indent=1, ensure_ascii=False)
print('ok', [(s['block'], s['shown']) for s in sources], sheet_path, [p for p, _ in tables])
