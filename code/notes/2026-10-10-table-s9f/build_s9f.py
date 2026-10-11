"""S9f candidate: (1) native one-sheet workbook with S9e styles and full float precision; (2) XML-level append to a copy of the
production workbook, so every existing sheet part stays byte-identical (same approach as reviews/2026-10-08-.../append_s9.py)."""
import copy, hashlib, json, math, os, re, sys, zipfile, xml.etree.ElementTree as E
from pathlib import Path
import pandas as pd, openpyxl
import openpyxl.compat.strings as _cs, openpyxl.cell._writer as _cw
def _safe(v, _old=_cs.safe_string):
    if isinstance(v, float) and not (math.isnan(v) or math.isinf(v)): return repr(v)
    return _old(v)
_cs.safe_string = _safe; _cw.safe_string = _safe          # keep full stored precision (openpyxl default is %.16g)
from openpyxl.worksheet.table import Table, TableStyleInfo
from openpyxl.worksheet.hyperlink import Hyperlink
from openpyxl.utils import get_column_letter as L
ROOT = Path(os.environ['REPRO_ROOT']); OUTDIR = Path(os.environ['S9F_OUT']); OUTDIR.mkdir(parents=True, exist_ok=True)
SRC = ROOT/'table/Supplementary_Tables_S1-S9.xlsx'; NATIVE = OUTDIR/'S9f-native.xlsx'; OUT = OUTDIR/'Supplementary_Tables_S1-S9.xlsx'
D = 'notes/2026-10-10-marker-reference-sensitivity/'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
COMMON = ('polyApipe primary catalog. D00, original P − original M; D01, original P − M*; D10, P* − original M; D11, P* − M*. '
          'Positive ΔDWUI is distal; fraction_proximal counts ΔDWUI ≤ −threshold. Rows with status not_estimable (N = 0) are not zero responses. '
          'Contributing libraries are not independent replicates.')
BLOCKS = [
    dict(title='Fixed original P with marker-defined M*, primary rule', file='unique_support_bridge_route_summary.tsv',
         q="background == 'common_M_bridge' and version == 'core_third_strict' and setting == 'supported50' and threshold == 0.02",
         rule='background == common_M_bridge, version == core_third_strict, setting == supported50, threshold == 0.02; no statistical recomputation.',
         note='D00 and D01 in each route share genes, contributing libraries and hierarchical weights. Liu2020 M* drew mainly on cells originally labelled S. Guan2022 has no same-library P–M baseline. '+COMMON),
    dict(title='Fixed original P with marker-defined M*, all marker versions, PAS settings and thresholds', file='unique_support_bridge_route_summary.tsv',
         q="background == 'common_M_bridge'", rule='background == common_M_bridge; all versions, settings and thresholds; no statistical recomputation.', note=COMMON),
    dict(title='Fixed original M with marker-defined P*', file='unique_support_bridge_route_summary.tsv',
         q="background == 'common_P_bridge'", rule='background == common_P_bridge; all rows; no statistical recomputation.',
         note='Its own common background; changes cannot be compared with the fixed-P block to apportion shares. '+COMMON),
    dict(title='Joint redefinition of P* and M*', file='unique_support_route_summary.tsv', q=None, rule='All rows; no statistical recomputation.',
         note='common_four: genes and libraries measurable in D00, D01, D10 and D11, where mean ΔDWUI satisfies D11 − D00 = (D01 − D00) + (D10 − D00); full: each contrast on its own background. Under the primary rule Guo2019–CIP is not estimable on common_four, and under the native panel its P* came almost entirely from cells originally labelled E. '+COMMON),
    dict(title='Pluripotency-score tertiles within original M', file='gradient_route_summary.tsv', q=None, rule='All rows; no statistical recomputation.',
         note='Original P against low, middle and high pluripotency-score thirds of original M within each library, on genes and libraries common to all three. Score tertiles are not calibrated purity levels. Liu2020 and Guan2022 are not estimable.'),
    dict(title='Original labels of marker-defined classes per library', file='classification_confusion.tsv', q=None, rule='All rows; no statistical recomputation.',
         note='Counts of candidate metacells by marker class (P*, M*, unassigned/mixed) and original label, per library and marker version.'),
]
# ---------------------------------------------------------------- 1. native sheet
base = openpyxl.load_workbook(SRC, read_only=False); t = base['S9e']
ST = {k: dict(font=copy.copy(t[a].font), fill=copy.copy(t[a].fill), alignment=copy.copy(t[a].alignment), border=copy.copy(t[a].border), number_format=t[a].number_format)
      for k, a in dict(title='A1', caption='A2', link='A3', block='A5', note='A6', shown='A7', head='A8', text='A9', num='I9').items()}
del base
nw = openpyxl.Workbook(); ws = nw.active; ws.title = 'S9f'; ws.sheet_view.showGridLines = False
def put(r, c, v, k, nf=None):
    cell = ws.cell(r, c); cell.value = v
    for a in ('font', 'fill', 'alignment', 'border'): setattr(cell, a, copy.copy(ST[k][a]))
    cell.number_format = nf or ST[k]['number_format']; return cell
put(1, 1, 'Table S9f. Marker-defined P and M references', 'title')
put(2, 1, 'Table S9. Analyses specified after the initial results. Reference-state sensitivity of P–M comparisons.', 'caption')
put(3, 1, '← Contents', 'link')
for r in (1, 2, 3): ws.row_dimensions[r].height = 20.0
r = 5; widths = {}; sources = []
for i, b in enumerate(BLOCKS, 1):
    path = ROOT/D/b['file']; df = pd.read_csv(path, sep='\t'); n_src = len(df)
    if b['q']: df = df.query(b['q'])
    put(r, 1, f"S9f.{i}  {b['title']} ({len(df)} rows)", 'block'); put(r+1, 1, b['note'], 'note'); put(r+2, 1, 'Rows shown: '+b['rule'], 'shown')
    h = r+3; ws.row_dimensions[h].height = 42.0
    for j, col in enumerate(df.columns, 1): put(h, j, col, 'head')
    for k, row in enumerate(df.itertuples(index=False), h+1):
        for j, v in enumerate(row, 1):
            if isinstance(v, float) and math.isnan(v): put(k, j, 'NA', 'text')
            elif type(v).__name__ in ('bool', 'bool_'): put(k, j, str(bool(v)), 'text')
            elif type(v).__name__.startswith(('int', 'uint')): put(k, j, int(v), 'num', '0')
            elif isinstance(v, float) or type(v).__name__.startswith('float'): put(k, j, float(v), 'num', '0.0000')
            else: put(k, j, str(v), 'text')
    last = h+len(df)
    tb = Table(displayName=f'T_S9f_{i}', ref=f'A{h}:{L(len(df.columns))}{max(last, h+1)}')
    tb.tableStyleInfo = TableStyleInfo(name=None, showRowStripes=False, showColumnStripes=False); ws.add_table(tb)
    for j, col in enumerate(df.columns, 1):
        sample = [str(x) for x in df.iloc[:300, j-1]]
        widths[j] = max(widths.get(j, 7), min(40, max(len(col)*1.1, max([0]+[len(x) for x in sample]))+3))
    sources.append([f'S9f.{i}', D+b['file'], len(df), n_src, b['rule']+' '+b['note'], sha(path)])
    r = last+3
for j, w in widths.items(): ws.column_dimensions[L(j)].width = w
ws.freeze_panes = 'B4'
nw.save(NATIVE)
# ---------------------------------------------------------------- 2. XML-level append
N = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'; R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
P = 'http://schemas.openxmlformats.org/package/2006/relationships'; C = 'http://schemas.openxmlformats.org/package/2006/content-types'
NS = {'s': N}; S = lambda x: '{'+N+'}'+x
def xml(tree, namespace=N):
    E.register_namespace('', namespace); E.register_namespace('r', R)
    return E.tostring(tree, encoding='utf-8', xml_declaration=True)
def cell(row, ref, value, style):
    c = E.SubElement(row, S('c'), {'r': ref, 's': str(style)})
    if isinstance(value, (int, float)): E.SubElement(c, S('v')).text = repr(value) if isinstance(value, float) else str(value)
    elif value is not None:
        c.set('t', 'inlineStr'); E.SubElement(E.SubElement(c, S('is')), S('t')).text = str(value)
def rowstyles(tree, number):
    row = tree.find('.//s:row[@r="'+str(number)+'"]', NS); return [int(c.get('s', '0')) for c in row]
def append_hyper(tree, ref, location, display):
    links = tree.find('s:hyperlinks', NS)
    if links is None:
        links = E.Element(S('hyperlinks'))
        ix = next((i for i, e in enumerate(tree) if e.tag in {S(t) for t in ['printOptions', 'pageMargins', 'pageSetup', 'headerFooter', 'tableParts', 'extLst']}), len(tree))
        tree.insert(ix, links)
    E.SubElement(links, S('hyperlink'), {'ref': ref, 'location': location, 'display': display})
parts = {}
with zipfile.ZipFile(SRC) as old, zipfile.ZipFile(NATIVE) as new:
    os_ = E.fromstring(old.read('xl/styles.xml')); ns_ = E.fromstring(new.read('xl/styles.xml'))
    oldnum = os_.find(S('numFmts')); newnum = ns_.find(S('numFmts'))
    codes = {el.get('formatCode'): int(el.get('numFmtId')) for el in oldnum}; nextnum = max([163]+list(codes.values()))+1; nummap = {}
    if newnum is not None:
        for el in newnum:
            code = el.get('formatCode'); idx = int(el.get('numFmtId'))
            if code not in codes:
                ce = copy.deepcopy(el); ce.set('numFmtId', str(nextnum)); oldnum.append(ce); codes[code] = nextnum; nextnum += 1
            nummap[idx] = codes[code]
    oldnum.set('count', str(len(oldnum))); offsets = {}
    for name in ['fonts', 'fills', 'borders', 'cellStyleXfs', 'cellXfs', 'dxfs']:
        dest = os_.find(S(name)); source = ns_.find(S(name))
        if source is None: continue
        if dest is None:
            dest = E.Element(S(name), {'count': '0'}); ix = next((i for i, e in enumerate(os_) if e.tag in {S('tableStyles'), S('extLst')}), len(os_)); os_.insert(ix, dest)
        offsets[name] = len(dest)
        for item in source:
            item = copy.deepcopy(item)
            if name in ['cellStyleXfs', 'cellXfs']:
                for attr, group in [('fontId', 'fonts'), ('fillId', 'fills'), ('borderId', 'borders')]:
                    if attr in item.attrib: item.set(attr, str(int(item.get(attr))+offsets[group]))
                if 'xfId' in item.attrib: item.set('xfId', str(int(item.get('xfId'))+offsets.get('cellStyleXfs', 0)))
                if int(item.get('numFmtId', '0')) in nummap: item.set('numFmtId', str(nummap[int(item.get('numFmtId'))]))
            dest.append(item)
        dest.set('count', str(len(dest)))
    shared = E.fromstring(new.read('xl/sharedStrings.xml')) if 'xl/sharedStrings.xml' in new.namelist() else None
    owb = E.fromstring(old.read('xl/workbook.xml')); orel = E.fromstring(old.read('xl/_rels/workbook.xml.rels'))
    oldrels = {e.get('Id'): e.get('Target') for e in orel}
    sheets = owb.find(S('sheets')); source_sheet = next(s for s in sheets if s.get('name') == 'Sources')
    source_path = 'xl/'+oldrels[source_sheet.get('{'+R+'}id')]
    linkstyle = E.fromstring(old.read(source_path)).find('.//s:c[@r="A3"]', NS).get('s')
    source_index = list(sheets).index(source_sheet)
    nextsheet = max(int(s.get('sheetId')) for s in sheets)+1
    nextrid = max(int(e.get('Id').removeprefix('rId')) for e in orel)+1
    nexttable = max(int(re.search(r'table(\d+)\.xml$', p).group(1)) for p in old.namelist() if re.match(r'xl/tables/table\d+\.xml$', p))+1
    sheetfile = existing = [p for p in old.namelist() if re.match(r'xl/worksheets/sheet\d+\.xml$', p)]
    while f'xl/worksheets/sheet{nextsheet}.xml' in existing: nextsheet += 1
    tree = E.fromstring(new.read('xl/worksheets/sheet1.xml'))
    for c in tree.iter(S('c')):
        if 's' in c.attrib: c.set('s', str(int(c.get('s'))+offsets['cellXfs']))
        if c.get('t') == 's':
            v = c.find(S('v')); si = shared[int(v.text)]; c.remove(v); c.set('t', 'inlineStr'); inline = E.SubElement(c, S('is'))
            for x in si: inline.append(copy.deepcopy(x))
    for e in list(tree.iter(S('row')))+list(tree.iter(S('col'))):
        attr = 's' if e.tag == S('row') else 'style'
        if attr in e.attrib: e.set(attr, str(int(e.get(attr))+offsets['cellXfs']))
    tree.find('.//s:c[@r="A3"]', NS).set('s', linkstyle); append_hyper(tree, 'A3', "'Contents'!A1", '← Contents')
    newpath = f'xl/worksheets/sheet{nextsheet}.xml'; newparts = [newpath]
    sr = E.fromstring(new.read('xl/worksheets/_rels/sheet1.xml.rels'))
    for rr in sr:
        assert rr.get('Type').endswith('/table'), rr.attrib
        tpath = 'xl/tables/'+Path(rr.get('Target')).name
        table = E.fromstring(new.read(tpath)); table.set('id', str(nexttable))
        tdest = f'xl/tables/table{nexttable}.xml'; parts[tdest] = xml(table); newparts.append(tdest)
        rr.set('Target', f'../tables/table{nexttable}.xml'); nexttable += 1
    parts[f'xl/worksheets/_rels/sheet{nextsheet}.xml.rels'] = xml(sr, P); parts[newpath] = xml(tree)
    s = E.Element(S('sheet'), {'name': 'S9f', 'sheetId': str(nextsheet), '{'+R+'}id': f'rId{nextrid}'})
    sheets.insert(list(sheets).index(source_sheet), s)
    E.SubElement(orel, '{'+P+'}Relationship', {'Id': f'rId{nextrid}', 'Type': R+'/worksheet', 'Target': f'worksheets/sheet{nextsheet}.xml'})
    for dn in owb.findall('s:definedNames/s:definedName', NS):
        if dn.get('localSheetId') is not None and int(dn.get('localSheetId')) >= source_index: dn.set('localSheetId', str(int(dn.get('localSheetId'))+1))
    os_.find(S('cellXfs')).set('count', str(len(os_.find(S('cellXfs')))))
    parts['xl/styles.xml'] = xml(os_); parts['xl/workbook.xml'] = xml(owb); parts['xl/_rels/workbook.xml.rels'] = xml(orel, P)
    ct = E.fromstring(old.read('[Content_Types].xml'))
    for p in newparts:
        kind = 'worksheet' if '/worksheets/' in p else 'table'
        E.SubElement(ct, '{'+C+'}Override', {'PartName': '/'+p, 'ContentType': f'application/vnd.openxmlformats-officedocument.spreadsheetml.{kind}+xml'})
    parts['[Content_Types].xml'] = xml(ct, C)
    contents_path = 'xl/'+oldrels[next(s for s in sheets if s.get('name') == 'Contents').get('{'+R+'}id')]
    contents = E.fromstring(old.read(contents_path)); sd = contents.find(S('sheetData'))
    rlast = max(int(row.get('r')) for row in sd); styles = rowstyles(contents, rlast-1)
    row = E.SubElement(sd, S('row'), {'r': str(rlast+1)})
    for j, v in enumerate(['S9f', 'Marker-defined P and M references', 'Reference-state sensitivity of P–M comparisons', len(BLOCKS), sum(x[2] for x in sources)]):
        cell(row, f'{chr(65+j)}{rlast+1}', v, styles[j])
    append_hyper(contents, f'A{rlast+1}', "'S9f'!A1", 'S9f'); contents.find(S('dimension')).set('ref', f'A1:E{rlast+1}'); parts[contents_path] = xml(contents)
    sources_x = E.fromstring(old.read(source_path)); sd = sources_x.find(S('sheetData')); rr_ = max(int(row.get('r')) for row in sd)
    for i, vals in enumerate(sources):
        st_ = rowstyles(sources_x, rr_-1 if i % 2 == 0 else rr_)
        row = E.SubElement(sd, S('row'), {'r': str(rr_+1+i)})
        for j, v in enumerate(vals): cell(row, f'{chr(65+j)}{rr_+1+i}', v, st_[j])
    dim = sources_x.find(S('dimension'))
    if dim is not None: dim.set('ref', re.sub(r'\d+$', str(rr_+len(sources)), dim.get('ref')))
    for t_ in sources_x.iter(S('autoFilter')): pass
    parts[source_path] = xml(sources_x)
    with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED) as z:
        for info in old.infolist():
            z.writestr(info, parts.pop(info.filename) if info.filename in parts else old.read(info.filename))
        for name, data in parts.items(): z.writestr(name, data)
json.dump(dict(source_workbook='table/Supplementary_Tables_S1-S9.xlsx', source_sha256=sha(SRC), output='notes/2026-10-10-table-s9f/Supplementary_Tables_S1-S9.xlsx',
               output_sha256=sha(OUT), native_sha256=sha(NATIVE), blocks=sources, production_replaced=False), open(OUTDIR/'build.json', 'w'), indent=1, ensure_ascii=False)
print('ok', [(x[0], x[2]) for x in sources])
