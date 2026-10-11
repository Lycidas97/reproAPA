"""Merge step only: append the prebuilt S9f-native.xlsx (made by build_s9f.py, which needs ~3 GB RAM to read the workbook styles)
to a copy of the production workbook at XML level. Existing sheet parts stay byte-identical. Writes only to this folder."""
import copy, hashlib, json, os, re, zipfile, xml.etree.ElementTree as E
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
SRC = ROOT/'table/Supplementary_Tables_S1-S9.xlsx'; NATIVE = HERE/'S9f-native.xlsx'; OUT = HERE/'Supplementary_Tables_S1-S9.xlsx'; OUTDIR = HERE
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
built = json.loads((HERE/'build.json').read_text()); sources = built['blocks']; BLOCKS = sources
assert sha(SRC) == built['source_sha256'] and sha(NATIVE) == built['native_sha256']
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
built['output_sha256_device'] = sha(OUT); (HERE/'build.json').write_text(json.dumps(built, indent=1, ensure_ascii=False))
print('ok', [(x[0], x[2]) for x in sources])
