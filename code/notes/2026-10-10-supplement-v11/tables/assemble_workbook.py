#!/usr/bin/env python3
"""Assemble the v11 supplementary-table workbook (candidate) from notes/2026-10-10-table-s9f/Supplementary_Tables_S1-S9.xlsx.
Production table/ is never written. XML-level, so every part not listed below is copied byte-for-byte.
  * old worksheets: renamed and reordered by spec_v11.MAP/ORDER. In each, only the title (A1), the caption (A2) and the
    block labels in column A are rewritten (inline strings); hyperlink targets and Excel table names follow the new names.
  * NOTE_FIXES: shared-string notes citing figure or table numbers that changed; each must match exactly once.
  * spec_v11.NEW worksheets are written with the cell styles of the S9f worksheet; every source row is kept.
  * spec_v11.APPEND blocks follow the last block of a carried-over worksheet in that worksheet's label, header and row
    styles, with one note row; the existing rows are not touched. Each source row is kept.
  * Contents and Sources are rebuilt in the new order with their existing row styles. A dagger marks worksheets specified
    after the initial results (spec_v11.NEW and POSTHOC_OLD)."""
import csv, hashlib, html, json, math, re, sys, zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(q for q in [HERE, *HERE.parents] if (q / 'notes').is_dir())
sys.path.insert(0, str(HERE))
import spec_v11 as spec

SRC = ROOT / 'notes/2026-10-10-table-s9f/Supplementary_Tables_S1-S9.xlsx'
OUT = HERE.parent / 'out' / 'Supplementary_Tables_S1-S10.v11.xlsx'
REPORT = HERE.parent / 'out' / 'Supplementary_Tables_S1-S10.v11.build.json'
N = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
PR = 'http://schemas.openxmlformats.org/package/2006/relationships'
sha = lambda b: hashlib.sha256(b).hexdigest()

NOTE_FIXES = [  # (old text, new text): figure/table numbers that changed in v11
    ('Copied from the Figure 2e panel values: correlation (r) and RMSE between held-out-route means of full-fit class members and either the centres refitted without that route or the full-fit centres.',
     'Correlation (r) and RMSE between held-out-route means of full-fit class members and either the centres refitted without that route or the full-fit centres, copied from the plotted values of the route hold-out panel (the panel itself is not shown).'),
    ('in the same_series_process_classes block (primary analysis); these are the counts drawn in Figure S11c.',
     'in the same_series_process_classes block (primary analysis).'),
    ('(primary analysis of Figure S11c; the 26 sensitivity analyses', '(primary analysis; the 26 sensitivity analyses'),
    ('the rug rule of Figure 5f; Trip12 is drawn in Figure S20h, and Figure 5f shows Riok3.',
     'the rug rule of Figure 4d; Trip12 is drawn in Figure S15h, and Figure 4d shows Riok3.'),
    ('interval analyses of Tables S6–S8,', 'interval analyses of Tables S5, S6 and S9,'),
]
CONTENTS_TITLE = 'Supplementary Tables S1–S10'
CONTENTS_NOTE = ('Worksheets carried over from the frozen polyApipe v2-primary workbook keep its values and formats. Worksheets marked † hold '
                 'analyses specified after the initial results; their blocks are copied from the source tables, or derived as stated in the block, '
                 'without statistical recomputation. Every block is filterable, with full stored precision; NA and missing cells remain distinct '
                 'from zero. Sources records exact filters and file SHA-256.')
SOURCES_NOTE = ('Paths beginning with notes/ or reviews/ are project-relative; other paths keep the bundle-relative convention of the frozen '
                'workbook, and entries beginning “derived from” name the inputs of a derived block. Filtered source files retain the original '
                'full-file SHA-256 here.')

# ------------------------------------------------------------------ helpers
CTRL = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f]')
def esc(s):
    s = html.escape(s, quote=False)
    return CTRL.sub(lambda m: '_x%04X_' % ord(m.group(0)), s)
def col(n):
    s = ''
    while n: n, r = divmod(n - 1, 26); s = chr(65 + r) + s
    return s
def c_str(ref, s, v): return f'<c r="{ref}" s="{s}" t="inlineStr"><is><t xml:space="preserve">{esc(v)}</t></is></c>'
def c_num(ref, s, v): return f'<c r="{ref}" s="{s}" t="n"><v>{v}</v></c>'
def c_empty(ref, s): return f'<c r="{ref}" s="{s}"/>'
CELL_A = re.compile(r'<c r="A(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', re.S)
CELL = re.compile(r'<c r="([A-Z]+)(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', re.S)
ROW = re.compile(r'<row r="(\d+)"[^>]*?(?:/>|>.*?</row>)', re.S)
def style_of(attrs):
    m = re.search(r'\bs="(\d+)"', attrs); return m.group(1) if m else '0'

z = zipfile.ZipFile(SRC)
ss_raw = z.read('xl/sharedStrings.xml').decode('utf-8')
SI = re.findall(r'<si>(.*?)</si>', ss_raw, re.S)
def si_text(i): return html.unescape(''.join(re.findall(r'<t[^>]*>(.*?)</t>', SI[i], re.S)))
def value(attrs, inner):
    inner = inner or ''
    t = re.search(r'\bt="(\w+)"', attrs); t = t.group(1) if t else 'n'
    if t == 's': return si_text(int(re.search(r'<v>(\d+)</v>', inner).group(1)))
    if t == 'inlineStr': return html.unescape(''.join(re.findall(r'<t[^>]*>(.*?)</t>', inner, re.S)))
    v = re.search(r'<v>(.*?)</v>', inner, re.S)
    if not v: return None
    if t in ('str', 'e'): return html.unescape(v.group(1))
    if t == 'b': return v.group(1) == '1'
    return int(v.group(1)) if re.fullmatch(r'-?\d+', v.group(1)) else float(v.group(1))
def cells_of_row(rowxml): return {c: (a, i) for c, _, a, i in CELL.findall(rowxml)}

# ------------------------------------------------------------------ workbook structure
wb = z.read('xl/workbook.xml').decode()
wrels = z.read('xl/_rels/workbook.xml.rels').decode()
RELS = {}
for m in re.finditer(r'<Relationship\b([^>]*?)/?>', wrels):
    a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1))); RELS[a['Id']] = a['Target']
def part(target): return 'xl/' + target.lstrip('/').removeprefix('xl/')
SHEETS = []
for m in re.finditer(r'<sheet\b([^>]*?)/>', wb):
    a = dict(re.findall(r'([\w:]+)="([^"]*)"', m.group(1)))
    SHEETS.append(dict(name=a['name'], sheetId=int(a['sheetId']), rid=a['r:id'], path=part(RELS[a['r:id']])))
BYNAME = {s['name']: s for s in SHEETS}
assert set(BYNAME) == set(spec.MAP) | set(spec.DROP_OLD) | {'Contents', 'Sources'}, set(BYNAME) ^ (set(spec.MAP) | set(spec.DROP_OLD) | {'Contents', 'Sources'})
INV = {v: k for k, v in spec.MAP.items()}
POSTHOC = set(spec.NEW) | set(spec.POSTHOC_OLD)
def rel_text(sheet): return spec.REL.get(sheet) or spec.POSTHOC_OLD.get(sheet, '')

s9f = z.read(BYNAME['S9f']['path']).decode()
def sty(ref): return style_of(re.search(rf'<c r="{ref}"([^>]*)', s9f).group(1))
ST = dict(title=sty('A1'), caption=sty('A2'), link=sty('A3'), block=sty('A5'), note=sty('A6'), head=sty('A8'), text=sty('A9'),
          float=sty('J9'), int=sty('L9'))
assert ST == dict(title='43', caption='44', link='16', block='46', note='47', head='48', text='49', float='50', int='51'), ST
parts = {}
report = dict(source=str(SRC.relative_to(ROOT)), source_sha256=sha(SRC.read_bytes()), production_replaced=False,
              renamed={}, note_fixes=[], new_sheets={})

# ------------------------------------------------------------------ 1. note fixes (shared strings; Contents A3 may be inline)
contents_xml = z.read(BYNAME['Contents']['path']).decode()
for a, b in NOTE_FIXES:
    ea, eb = esc(a), esc(b)
    n_ss, n_ct = ss_raw.count(ea), contents_xml.count(ea)
    assert n_ss + n_ct == 1, (a[:60], n_ss, n_ct)
    if n_ss: ss_raw = ss_raw.replace(ea, eb)
    else: contents_xml = contents_xml.replace(ea, eb)
    report['note_fixes'].append(dict(old=a, new=b, where='sharedStrings' if n_ss else 'Contents'))
parts['xl/sharedStrings.xml'] = ss_raw.encode('utf-8')

# ------------------------------------------------------------------ 2. old worksheets
LABEL = re.compile(r'^(S\d+[a-z])(\.\d+)?(\s{2}.*)$', re.S)
COUNT = re.compile(r'(\d[\d,]*) rows')
RELEL = re.compile(r'<Relationship\b[^>]*?(?:/>|>\s*</Relationship>)', re.S)
def relink(xml):
    return re.sub(r'location="(?:\'|&apos;)([^\'&]+)(?:\'|&apos;)!',
                  lambda m: f'location="\'{spec.MAP.get(m.group(1), m.group(1))}\'!', xml)
def label_rows(xml, old):
    """old block index -> (label row, label text), read from the source sheet"""
    out = {}
    for m in CELL_A.finditer(xml):
        r = int(m.group(1))
        if r >= 5 and re.search(r'\bt="(s|inlineStr|str)"', m.group(2)):
            v = value(m.group(2), m.group(3)); mm = LABEL.match(v)
            if mm and mm.group(1) == old:
                i = int(mm.group(2)[1:]) if mm.group(2) else 1
                assert i not in out, (old, i); out[i] = (r, v)
    return out
def shift_ref(ref, o): return re.sub(r'([A-Z]+)(\d+)', lambda q: f'{q.group(1)}{int(q.group(2)) - o}', ref)
SKIP = set()                                   # source parts not written (dropped sheets and tables)
for old in spec.DROP_OLD:
    p = BYNAME[old]['path']; SKIP.add(p)
    rp = p.replace('worksheets/', 'worksheets/_rels/') + '.rels'
    if rp in z.namelist():
        SKIP.add(rp); SKIP.update('xl/tables/' + Path(t).name for t in re.findall(r'Target="([^"]+)"', z.read(rp).decode()))
labels, BID, info_old = {}, {}, {}
for old, new in spec.MAP.items():
    p = BYNAME[old]['path']; xml = z.read(p).decode('utf-8')
    assert '<f>' not in xml and '<f ' not in xml
    lab = label_rows(xml, old); order_ = sorted(lab); assert order_, old
    drop = set(spec.DROP_BLOCKS.get(old, [])); assert drop <= set(order_), (old, drop, order_)
    kept = [i for i in order_ if i not in drop]
    newidx = {i: k for k, i in enumerate(kept, 1)}; single = len(kept) == 1
    bid = lambda i, new=new, newidx=newidx, single=single: new if single else f'{new}.{newidx[i]}'
    BID[old] = (drop, bid)
    starts = [lab[i][0] for i in order_]
    span = {i: (lab[i][0], starts[k + 1] - 1 if k + 1 < len(starts) else 10 ** 9) for k, i in enumerate(order_)}
    dspans = [span[i] for i in sorted(drop)]
    def in_drop(r, dspans=dspans): return any(a <= r <= b for a, b in dspans)
    def off(r, dspans=dspans): return sum(b - a + 1 for a, b in dspans if b < r)
    def block_of(r, span=span): return next(i for i, (a, b) in span.items() if a <= r <= b)
    done = dict(A1=0, A2=0)
    def fix(m):
        r, attrs, inner = int(m.group(1)), m.group(2), m.group(3)
        s = style_of(attrs)
        if r == 1:
            v = value(attrs, inner); assert isinstance(v, str) and v.startswith(f'Table {old}. '), (old, v)
            done['A1'] += 1; return c_str('A1', s, f'Table {new}. ' + spec.TITLE_OVERRIDE.get(new, v[len(f'Table {old}. '):]))
        if r == 2:
            v = value(attrs, inner); assert isinstance(v, str) and (v.startswith('Part of Table') or v.startswith('Table S9.')), (old, v)
            done['A2'] += 1; return c_str('A2', s, spec.caption(new))
        if r >= 5 and re.search(r'\bt="(s|inlineStr|str)"', attrs):
            v = value(attrs, inner); mm = LABEL.match(v)
            if mm and mm.group(1) == old:
                i = int(mm.group(2)[1:]) if mm.group(2) else 1
                return m.group(0) if i in drop else c_str(f'A{r}', s, bid(i) + mm.group(3))
        return m.group(0)
    xml = CELL_A.sub(fix, xml)
    assert done == dict(A1=1, A2=1), (old, done)
    labels[new] = [bid(i) for i in kept]
    nrows = sum(int(COUNT.search(lab[i][1]).group(1).replace(',', '')) for i in kept)
    rels_path = p.replace('worksheets/', 'worksheets/_rels/') + '.rels'
    rels_x = z.read(rels_path).decode() if rels_path in z.namelist() else None
    if drop:                                     # remove dropped blocks and move the following rows up
        head, rest = xml.split('<sheetData>', 1); body, tail = rest.split('</sheetData>', 1)
        rows_ = list(ROW.finditer(body)); assert ''.join(m.group(0) for m in rows_) == body, old
        out = []
        for m in rows_:
            r = int(m.group(1))
            if in_drop(r): continue
            o = off(r); rx = m.group(0)
            if o:
                rx = re.sub(r'^<row r="\d+"', f'<row r="{r - o}"', rx)
                rx = re.sub(r'<c r="([A-Z]+)' + str(r) + '"', lambda q: f'<c r="{q.group(1)}{r - o}"', rx)
            out.append(rx)
        body = ''.join(out); maxr = max(int(x) for x in re.findall(r'<row r="(\d+)"', body))
        head = re.sub(r'<dimension ref="([A-Z]+\d+):([A-Z]+)\d+"', lambda q: f'<dimension ref="{q.group(1)}:{q.group(2)}{maxr}"', head)
        def cf(m):
            el = m.group(0); r0 = int(re.search(r'sqref="[A-Z]+(\d+)', el).group(1))
            if in_drop(r0): return ''
            o = off(r0)
            el = re.sub(r'sqref="([^"]+)"', lambda q: f'sqref="{shift_ref(q.group(1), o)}"', el)
            el, n = re.subn(r'ROW\(\)-(\d+)', lambda q: f'ROW()-{int(q.group(1)) - o}', el); assert n == 1, el[:200]
            return el
        tail = re.sub(r'<conditionalFormatting\b.*?</conditionalFormatting>', cf, tail, flags=re.S)
        assert '<mergeCell' not in tail and '<dataValidation' not in tail, old
        assert all(int(h) < min(a for a, _ in dspans) for h in re.findall(r'<hyperlink ref="[A-Z]+(\d+)"', tail)), old
        xml = head + '<sheetData>' + body + '</sheetData>' + tail
    if rels_x is not None:                       # Excel tables: drop those of dropped blocks, shift and rename the rest
        tabs = []
        for el in RELEL.findall(rels_x):
            if not el.split('Type="')[1].startswith(R + '/table'): continue
            rid = re.search(r'\bId="([^"]+)"', el).group(1); tp = 'xl/tables/' + Path(re.search(r'Target="([^"]+)"', el).group(1)).name
            tx = z.read(tp).decode('utf-8'); h = int(re.search(r'\sref="[A-Z]+(\d+):', tx).group(1))
            tabs.append((h, rid, tp, tx, el))
        k = 0
        for h, rid, tp, tx, el in sorted(tabs):
            if drop and block_of(h) in drop:
                SKIP.add(tp); rels_x = rels_x.replace(el, '')
                xml, n = re.subn(r'<tablePart r:id="' + re.escape(rid) + r'"\s*/>', '', xml); assert n == 1, (old, rid)
                continue
            k += 1; o = off(h) if drop else 0
            a = tx.index('<table'); e = tx.index('>', a) + 1
            tag, n = re.subn(r'\b(name|displayName)="T_' + re.escape(old) + r'_\d+"', lambda q: f'{q.group(1)}="T_{new}_{k}"', tx[a:e]); assert n == 2, (tp, tx[a:e][:200])
            tx = tx[:a] + tag + tx[e:]
            if o: tx = re.sub(r'(\sref=")([^"]+)"', lambda q: f'{q.group(1)}{shift_ref(q.group(2), o)}"', tx)
            parts[tp] = tx.encode('utf-8')
        if drop:
            xml = re.sub(r'<tableParts count="\d+"', f'<tableParts count="{k}"', xml)
            parts[rels_path] = rels_x.encode('utf-8')
    parts[p] = relink(xml).encode('utf-8')
    info_old[new] = (len(kept), nrows)
    report['renamed'][old] = dict(new=new, blocks=labels[new], dropped_blocks=sorted(drop))
report['dropped_sheets'] = spec.DROP_OLD

# ------------------------------------------------------------------ 3. new worksheets
INT = re.compile(r'-?\d+'); FLOAT = re.compile(r'-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?')
def typed(ref, v):
    if v in ('', 'NA', 'nan', 'NaN'): return c_str(ref, ST['text'], 'NA')
    if re.fullmatch(r'-?0\d+.*', v): return c_str(ref, ST['text'], v)            # codes with leading zeros stay text
    if INT.fullmatch(v) and len(v.lstrip('-')) <= 15: return c_num(ref, ST['int'], int(v))
    if FLOAT.fullmatch(v) and math.isfinite(float(v)): return c_num(ref, ST['float'], repr(float(v)))
    return c_str(ref, ST['text'], v)
existing = z.namelist()
next_sheet = max(int(re.search(r'sheet(\d+)\.xml$', n).group(1)) for n in existing if re.match(r'xl/worksheets/sheet\d+\.xml$', n)) + 1
next_table = max(int(re.search(r'table(\d+)\.xml$', n).group(1)) for n in existing if re.match(r'xl/tables/table\d+\.xml$', n)) + 1
next_sid = max(s['sheetId'] for s in SHEETS) + 1
next_rid = max(int(k.removeprefix('rId')) for k in RELS) + 1
new_entries, sources_new, info_new = {}, [], {}

for key, sh in spec.NEW.items():
    title = sh['title']
    rows = [f'<row r="1" ht="20" customHeight="1">' + c_str('A1', ST['title'], f'Table {key}. {title}') + '</row>',
            f'<row r="2" ht="20" customHeight="1">' + c_str('A2', ST['caption'], spec.caption(key)) + '</row>',
            f'<row r="3" ht="20" customHeight="1">' + c_str('A3', ST['link'], '← Contents') + '</row>']
    r, maxc, widths, tables, ids, total = 5, 1, {1: 10.7}, [], [], 0
    for i, b in enumerate(sh['blocks'], 1):
        path = ROOT / b['f']; raw = path.read_bytes()
        with open(path, newline='', encoding='utf-8') as fh:
            data = list(csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE))
        hdr, body = data[0], data[1:]
        assert body and all(hdr) and len(set(hdr)) == len(hdr), (key, i, hdr)
        assert all(len(x) == len(hdr) for x in body), (key, i, 'ragged rows')
        bid = f'{key}.{i}'; ids.append(bid); total += len(body)
        label = f"{bid}  {b['t']} ({len(body)} rows)"
        rows.append(f'<row r="{r}">' + c_str(f'A{r}', ST['block'], label) + '</row>')
        rows.append(f'<row r="{r+1}">' + c_str(f'A{r+1}', ST['note'], b['note']) + '</row>')
        rows.append(f'<row r="{r+2}">' + c_str(f'A{r+2}', ST['note'], 'Rows shown: ' + b['rule']) + '</row>')
        h = r + 3
        rows.append(f'<row r="{h}" ht="42" customHeight="1">' + ''.join(c_str(f'{col(j)}{h}', ST['head'], x) for j, x in enumerate(hdr, 1)) + '</row>')
        for k, rec in enumerate(body, h + 1):
            rows.append(f'<row r="{k}">' + ''.join(typed(f'{col(j)}{k}', x) for j, x in enumerate(rec, 1)) + '</row>')
        last = h + len(body); maxc = max(maxc, len(hdr))
        for j, name in enumerate(hdr, 1):
            sample = [x[j - 1] for x in body[:300]]
            widths[j] = max(widths.get(j, 7), min(40, max(len(name) * 1.1, max([0] + [len(x) for x in sample])) + 3))
        ref = f'A{h}:{col(len(hdr))}{last}'
        tid = next_table; next_table += 1; tables.append(tid)
        parts[f'xl/tables/table{tid}.xml'] = (
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<table xmlns="{N}" id="{tid}" name="T_{key}_{i}" displayName="T_{key}_{i}" '
            f'ref="{ref}" headerRowCount="1"><autoFilter ref="{ref}"/><tableColumns count="{len(hdr)}">'
            + ''.join(f'<tableColumn id="{j}" name="{html.escape(x, quote=True)}"/>' for j, x in enumerate(hdr, 1))
            + '</tableColumns><tableStyleInfo showRowStripes="0" showColumnStripes="0"/></table>').encode('utf-8')
        sources_new.append((key, i, [bid, b['f'], len(body), len(body), b['rule'] + ' ' + b['note'], sha(raw)]))
        r = last + 3
    lastrow = r - 3
    xml = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<worksheet xmlns="{N}" xmlns:r="{R}">'
           '<sheetPr><outlinePr summaryBelow="1" summaryRight="1"/><pageSetUpPr/></sheetPr>'
           f'<dimension ref="A1:{col(maxc)}{lastrow}"/>'
           '<sheetViews><sheetView showGridLines="0" workbookViewId="0"><pane xSplit="1" ySplit="3" topLeftCell="B4" activePane="bottomRight" state="frozen"/>'
           '<selection pane="topRight"/><selection pane="bottomLeft"/><selection pane="bottomRight" activeCell="A1" sqref="A1"/></sheetView></sheetViews>'
           '<sheetFormatPr baseColWidth="8" defaultRowHeight="15"/><cols>'
           + ''.join(f'<col min="{j}" max="{j}" width="{w:.1f}" customWidth="1"/>' for j, w in sorted(widths.items()))
           + '</cols><sheetData>' + ''.join(rows) + '</sheetData>'
           '<hyperlinks><hyperlink ref="A3" location="\'Contents\'!A1" display="← Contents"/></hyperlinks>'
           '<pageMargins left="0.75" right="0.75" top="1" bottom="1" header="0.5" footer="0.5"/>'
           f'<tableParts count="{len(tables)}">' + ''.join(f'<tablePart r:id="rId{k}"/>' for k in range(1, len(tables) + 1)) + '</tableParts></worksheet>')
    sp = f'xl/worksheets/sheet{next_sheet}.xml'
    parts[sp] = xml.encode('utf-8')
    parts[f'xl/worksheets/_rels/sheet{next_sheet}.xml.rels'] = (
        f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="{PR}">'
        + ''.join(f'<Relationship Id="rId{k}" Type="{R}/table" Target="../tables/table{t}.xml"/>' for k, t in enumerate(tables, 1))
        + '</Relationships>').encode('utf-8')
    new_entries[key] = dict(sheetId=next_sid, rid=f'rId{next_rid}', path=sp, tables=tables)
    next_sheet += 1; next_sid += 1; next_rid += 1
    labels[key] = ids; info_new[key] = (title, len(sh['blocks']), total)
    report['new_sheets'][key] = dict(part=sp, rows=total, blocks=[x[2] for x in sources_new if x[0] == key])

# ------------------------------------------------------------------ 3b. blocks appended to carried-over worksheets (table ids after the new worksheets)
APP_HEAD = dict(text='3', num='4')                                    # header styles of the carried-over worksheets
APP_ROW = [dict(text='8', int='9', float='17'), dict(text='12', int='13', float='18')]   # alternating data-row styles
def kind_of(values):
    vals = [v for v in values if v != '']
    if vals and all(INT.fullmatch(v) and len(v.lstrip('-')) <= 15 and not re.fullmatch(r'-?0\d+', v) for v in vals): return 'int'
    if vals and all(FLOAT.fullmatch(v) and math.isfinite(float(v)) for v in vals): return 'float'
    return 'text'
def app_cell(ref, v, kind, st):
    if v == '': return c_empty(ref, st['text' if kind == 'text' else kind])
    if kind == 'int': return c_num(ref, st['int'], int(v))
    if kind == 'float': return c_num(ref, st['float'], repr(float(v)))
    return c_str(ref, st['text'], v)
info_app, app_tables = {}, []
report['appended'] = {}
for new, blocks in spec.APPEND.items():
    old = INV[new]; p = BYNAME[old]['path']; xml = parts[p].decode('utf-8')
    head, rest = xml.split('<sheetData>', 1); body, tail = rest.split('</sheetData>', 1)
    # the worksheet's own styles, read from its first block
    lab_r = next(int(m.group(1)) for m in CELL_A.finditer(body) if int(m.group(1)) >= 5 and LABEL.match(str(value(m.group(2), m.group(3)) or '')))
    rows_b = {int(m.group(1)): m.group(0) for m in ROW.finditer(body)}
    assert style_of(CELL_A.search(rows_b[lab_r]).group(2)) == '0', (new, 'label style')
    hs = {style_of(a) for _, _, a, _ in CELL.findall(rows_b[lab_r + 1])}; assert hs <= set(APP_HEAD.values()), (new, hs)
    for k, st in ((lab_r + 2, APP_ROW[0]), (lab_r + 3, APP_ROW[1])):
        ds = {style_of(a) for _, _, a, _ in CELL.findall(rows_b[k])}; assert ds <= set(st.values()), (new, k, ds)
    maxr = max(rows_b); r = maxr + 3; out, ids, total = [], [], 0
    rels_path = p.replace('worksheets/', 'worksheets/_rels/') + '.rels'
    rels_x = parts[rels_path].decode() if rels_path in parts else z.read(rels_path).decode()
    next_r = max(int(x) for x in re.findall(r'\bId="rId(\d+)"', rels_x)) + 1
    tparts, maxc = [], 0
    for i, b in enumerate(blocks, len(labels[new]) + 1):
        path = ROOT / b['f']; raw = path.read_bytes()
        with open(path, newline='', encoding='utf-8') as fh:
            data = list(csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE))
        hdr, recs = data[0], data[1:]
        assert recs and all(hdr) and len(set(hdr)) == len(hdr) and all(len(x) == len(hdr) for x in recs), (new, i)
        kinds = [kind_of([x[j] for x in recs]) for j in range(len(hdr))]
        bid = f'{new}.{i}'; ids.append(bid); total += len(recs); nc = len(hdr); maxc = max(maxc, nc); sp = f' spans="1:{nc}"'
        out.append(f'<row r="{r}"{sp}>' + c_str(f'A{r}', '0', f"{bid}  {b['name']}    {len(recs):,} rows") + '</row>')
        out.append(f'<row r="{r+1}"{sp}>' + c_str(f'A{r+1}', ST['note'], b['note']) + '</row>')
        h = r + 2
        out.append(f'<row r="{h}"{sp}>' + ''.join(c_str(f'{col(j)}{h}', APP_HEAD['text' if kd == 'text' else 'num'], x)
                                              for j, (x, kd) in enumerate(zip(hdr, kinds), 1)) + '</row>')
        for k, rec in enumerate(recs):
            rn = h + 1 + k; st = APP_ROW[k % 2]
            out.append(f'<row r="{rn}"{sp}>' + ''.join(app_cell(f'{col(j)}{rn}', v, kd, st) for j, (v, kd) in enumerate(zip(rec, kinds), 1)) + '</row>')
        last = h + len(recs); ref = f'A{h}:{col(nc)}{last}'
        tid = next_table; next_table += 1; app_tables.append(tid); rid = f'rId{next_r}'; next_r += 1
        parts[f'xl/tables/table{tid}.xml'] = (
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<table xmlns="{N}" id="{tid}" name="T_{new}_{i}" displayName="T_{new}_{i}" '
            f'ref="{ref}" totalsRowShown="0"><autoFilter ref="{ref}"/><tableColumns count="{nc}">'
            + ''.join(f'<tableColumn id="{j}" name="{html.escape(x, quote=True)}"/>' for j, x in enumerate(hdr, 1))
            + '</tableColumns><tableStyleInfo showFirstColumn="0" showLastColumn="0" showRowStripes="1" showColumnStripes="0"/></table>').encode('utf-8')
        rels_x = rels_x.replace('</Relationships>', f'<Relationship Id="{rid}" Type="{R}/table" Target="../tables/table{tid}.xml"/></Relationships>')
        tparts.append(f'<tablePart r:id="{rid}"/>')
        sources_new.append((new, i, [bid, b['f'], len(recs), len(recs), spec.ALL + ' ' + b['note'], sha(raw)]))
        r = last + 3
    lastrow = r - 3
    m = re.search(r'<dimension ref="A1:([A-Z]+)(\d+)"/>', head); assert m, new
    def cnum(c):
        n = 0
        for ch in c: n = n * 26 + ord(ch) - 64
        return n
    head = head.replace(m.group(0), f'<dimension ref="A1:{col(max(cnum(m.group(1)), maxc))}{lastrow}"/>')
    defined = max(int(x) for x in re.findall(r'<col min="\d+" max="(\d+)"', head))
    if maxc > defined:
        head = head.replace('</cols>', f'<col min="{defined + 1}" max="{maxc}" width="18" customWidth="1"/></cols>')
    tc = re.search(r'<tableParts count="(\d+)">(.*?)</tableParts>', tail, re.S); assert tc, new
    tail = tail.replace(tc.group(0), f'<tableParts count="{int(tc.group(1)) + len(tparts)}">' + tc.group(2) + ''.join(tparts) + '</tableParts>')
    parts[p] = (head + '<sheetData>' + body + ''.join(out) + '</sheetData>' + tail).encode('utf-8')
    parts[rels_path] = rels_x.encode('utf-8')
    labels[new] = labels[new] + ids; info_app[new] = (len(blocks), total)
    report['appended'][new] = dict(blocks=ids, rows=total, first_row=maxr + 3, last_row=lastrow, sources=[b['f'] for b in blocks])
# ------------------------------------------------------------------ 4. Contents
rows_ct = {int(m.group(1)): m.group(0) for m in ROW.finditer(contents_xml)}
old_info = {}
for rn, rx in rows_ct.items():
    if rn < 7: continue
    c = cells_of_row(rx)
    if 'A' in c and value(*c['A']) in spec.MAP:
        old_info[value(*c['A'])] = tuple(value(*c[L]) if L in c else None for L in 'BDE')
assert set(old_info) == set(spec.MAP), set(spec.MAP) - set(old_info)
assert value(*cells_of_row(rows_ct[1])['A']) == 'Supplementary Tables S1–S9'
assert value(*cells_of_row(rows_ct[2])['A']).startswith('S1–S8 values')
r1 = re.sub(r'<c r="A1".*?</c>', lambda m: c_str('A1', '1', CONTENTS_TITLE), rows_ct[1], count=1, flags=re.S)
r2 = re.sub(r'<c r="A2".*?</c>', lambda m: c_str('A2', '2', CONTENTS_NOTE), rows_ct[2], count=1, flags=re.S)
out_rows, links, cr = [r1, r2, rows_ct[3], rows_ct[5]], [], 6
sheets_only = [s for s in spec.ORDER if s not in ('Contents', 'Sources')]
for t in dict.fromkeys(spec.table_of(s) for s in sheets_only):
    out_rows.append(f'<row r="{cr}" spans="1:5">' + c_str(f'A{cr}', 5, f'Table {t}') + c_str(f'B{cr}', 5, spec.GROUPS[t])
                    + ''.join(c_empty(f'{L}{cr}', 6) for L in 'CDE') + '</row>'); cr += 1
    for k, s in enumerate([x for x in sheets_only if spec.table_of(x) == t]):
        st = (7, 8, 8, 9, 10) if k % 2 == 0 else (11, 12, 12, 13, 14)
        if s in spec.NEW: title, nb, nr = info_new[s]
        else:
            t0, nb0, nr0 = old_info[INV[s]]; nb, nr = info_old[s]
            if INV[s] not in spec.DROP_BLOCKS: assert (int(nb0), int(nr0)) == (nb, nr), (s, nb0, nr0, nb, nr)   # counts agree with old Contents
            if s in info_app: nb, nr = nb + info_app[s][0], nr + info_app[s][1]                         # appended blocks
            title = spec.TITLE_OVERRIDE.get(s, t0)
        rel = rel_text(s)
        out_rows.append(f'<row r="{cr}" spans="1:5">' + c_str(f'A{cr}', st[0], s) + c_str(f'B{cr}', st[1], title + (' †' if s in POSTHOC else ''))
                        + (c_str(f'C{cr}', st[2], rel) if rel else c_empty(f'C{cr}', st[2]))
                        + c_num(f'D{cr}', st[3], int(nb)) + c_num(f'E{cr}', st[4], int(nr)) + '</row>')
        links.append(f'<hyperlink ref="A{cr}" location="\'{s}\'!A1" display="{s}"/>'); cr += 1
ct = contents_xml
ct = re.sub(r'<sheetData>.*</sheetData>', lambda m: '<sheetData>' + ''.join(out_rows) + '</sheetData>', ct, count=1, flags=re.S)
ct = re.sub(r'<hyperlinks>.*?</hyperlinks>', lambda m: '<hyperlinks>' + ''.join(links) + '</hyperlinks>', ct, count=1, flags=re.S)
ct = re.sub(r'<dimension ref="[^"]*"', f'<dimension ref="A1:E{cr-1}"', ct, count=1)
parts[BYNAME['Contents']['path']] = ct.encode('utf-8')

# ------------------------------------------------------------------ 5. Sources
src_xml = z.read(BYNAME['Sources']['path']).decode('utf-8')
rows_sr = {int(m.group(1)): m.group(0) for m in ROW.finditer(src_xml)}
assert value(*cells_of_row(rows_sr[2])['A']).startswith('S1–S8 paths')
recs, src_fixes = [], {}
for rn in sorted(rows_sr):
    if rn < 6: continue
    c = cells_of_row(rows_sr[rn])
    v = [value(*c[L]) if L in c else None for L in 'ABCDEF']
    m = re.match(r'^(S\d+[a-z])(\.(\d+))?$', v[0]); assert m and (m.group(1) in spec.MAP or m.group(1) in spec.DROP_OLD), v[0]
    for j, x in enumerate(v):                       # Sources repeats block notes as inline text: apply the same fixes
        if isinstance(x, str):
            for a, b in NOTE_FIXES:
                if a in x: x = x.replace(a, b); src_fixes[a] = src_fixes.get(a, 0) + 1
            v[j] = x
    if m.group(1) in spec.DROP_OLD: continue
    drop, bid = BID[m.group(1)]; i = int(m.group(3) or 1)
    if i in drop: continue
    new = spec.MAP[m.group(1)]; v[0] = bid(i)
    project = str(v[1]).startswith(('notes/', 'reviews/'))
    assert project or new not in POSTHOC, (v[0], v[1])                      # dagger blocks are all project-relative
    assert not project or (ROOT / str(v[1])).exists(), (v[0], v[1])        # project-relative paths must exist
    recs.append((spec.ORDER.index(new), int(v[0].split('.')[1]) if '.' in v[0] else 1, v))
for key, i, v in sources_new: recs.append((spec.ORDER.index(key), i, v))
recs.sort(key=lambda x: (x[0], x[1]))
r2 = re.sub(r'<c r="A2".*?</c>', lambda m: c_str('A2', '15', SOURCES_NOTE), rows_sr[2], count=1, flags=re.S)
out_rows = [rows_sr[1], r2, rows_sr[3], rows_sr[5]]
links = ['<hyperlink ref="A3" location="\'Contents\'!A1" display="← Contents"/>']
for k, (_, _, v) in enumerate(recs):
    rn = 6 + k; st = (7, 8, 9, 9, 8, 8) if rn % 2 == 0 else (11, 12, 13, 13, 12, 12)
    cells = []
    for j, (L, x) in enumerate(zip('ABCDEF', v)):
        if x is None: cells.append(c_empty(f'{L}{rn}', st[j]))
        elif isinstance(x, (int, float)): cells.append(c_num(f'{L}{rn}', st[j], x))
        else: cells.append(c_str(f'{L}{rn}', st[j], str(x)))
    out_rows.append(f'<row r="{rn}">' + ''.join(cells) + '</row>')
    target = re.match(r'^(S\d+[a-z])', v[0]).group(1)
    links.append(f'<hyperlink ref="A{rn}" location="\'{target}\'!A1" display="{v[0]}"/>')
sx = src_xml
sx = re.sub(r'<sheetData>.*</sheetData>', lambda m: '<sheetData>' + ''.join(out_rows) + '</sheetData>', sx, count=1, flags=re.S)
sx = re.sub(r'<hyperlinks>.*?</hyperlinks>', lambda m: '<hyperlinks>' + ''.join(links) + '</hyperlinks>', sx, count=1, flags=re.S)
sx = re.sub(r'<dimension ref="([A-Z]+\d+):([A-Z]+)\d+"', lambda m: f'<dimension ref="{m.group(1)}:{m.group(2)}{5 + len(recs)}"', sx, count=1)
parts[BYNAME['Sources']['path']] = sx.encode('utf-8')
report['sources_note_fixes'] = src_fixes
assert all(src_fixes.get(a, 0) == 1 for a, _ in NOTE_FIXES[:4]), src_fixes   # the four block notes, once each
# the Sources Excel table must cover the rebuilt rows (header in row 5)
srels = BYNAME['Sources']['path'].replace('worksheets/', 'worksheets/_rels/') + '.rels'
for t in re.findall(r'Target="([^"]+)"', z.read(srels).decode()):
    tp = 'xl/tables/' + Path(t).name; tx = z.read(tp).decode('utf-8')
    tx, n = re.subn(r'ref="A5:F\d+"', f'ref="A5:F{5 + len(recs)}"', tx); assert n == 2, (tp, n)
    parts[tp] = tx.encode('utf-8')
in_sheets = sorted(b for ids in labels.values() for b in ids); in_sources = sorted(v[0] for _, _, v in recs)
assert in_sheets == in_sources, set(in_sheets) ^ set(in_sources)

# ------------------------------------------------------------------ 6. workbook, rels, content types, app.xml
assert spec.ORDER[0] == 'Contents' and spec.ORDER[-1] == 'Sources'
def entry(s): return BYNAME[s] if s in ('Contents', 'Sources') else new_entries[s] if s in spec.NEW else BYNAME[INV[s]]
sheet_xml = ''.join(f'<sheet name="{s}" sheetId="{entry(s)["sheetId"]}" r:id="{entry(s)["rid"]}" />' for s in spec.ORDER)
wb2, n = re.subn(r'<sheets>.*?</sheets>', lambda m: '<sheets>' + sheet_xml + '</sheets>', wb, count=1, flags=re.S); assert n == 1
parts['xl/workbook.xml'] = wb2.encode('utf-8')
for old in spec.DROP_OLD:
    rid = BYNAME[old]['rid']
    wrels, n = re.subn(r'<Relationship\b(?=[^>]*\bId="' + re.escape(rid) + r'")[^>]*?(?:/>|>\s*</Relationship>)', '', wrels); assert n == 1, (old, rid)
parts['xl/_rels/workbook.xml.rels'] = wrels.replace('</Relationships>', ''.join(
    f'<Relationship Id="{e["rid"]}" Type="{R}/worksheet" Target="worksheets/{Path(e["path"]).name}" />' for e in new_entries.values())
    + '</Relationships>').encode('utf-8')
ctypes = z.read('[Content_Types].xml').decode()
add = ''.join(f'<Override PartName="/{e["path"]}" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml" />'
              + ''.join(f'<Override PartName="/xl/tables/table{t}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.table+xml" />'
                        for t in e['tables']) for e in new_entries.values())
for t in app_tables: add += f'<Override PartName="/xl/tables/table{t}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.table+xml" />'
for sp in sorted(SKIP):
    if sp.endswith('.rels'): continue
    ctypes, n = re.subn(r'<Override PartName="/' + re.escape(sp) + r'"[^>]*/>', '', ctypes); assert n == 1, sp
parts['[Content_Types].xml'] = ctypes.replace('</Types>', add + '</Types>').encode('utf-8')
app = z.read('docProps/app.xml').decode()
pre = re.search(r'<(\w+):vector size="\d+" baseType="lpstr">', app).group(1)
app, n1 = re.subn(r'(<\w+:lpstr>Worksheets</\w+:lpstr></\w+:variant><\w+:variant><\w+:i4>)\d+', lambda m: m.group(1) + str(len(spec.ORDER)), app)
app, n2 = re.subn(r'<TitlesOfParts>.*?</TitlesOfParts>', lambda m: f'<TitlesOfParts><{pre}:vector size="{len(spec.ORDER)}" baseType="lpstr">'
                  + ''.join(f'<{pre}:lpstr>{s}</{pre}:lpstr>' for s in spec.ORDER) + f'</{pre}:vector></TitlesOfParts>', app, flags=re.S)
assert n1 == 1 and n2 == 1
parts['docProps/app.xml'] = app.encode('utf-8')

# ------------------------------------------------------------------ write
changed = sorted(parts)
with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED) as zo:
    for info in z.infolist():
        if info.filename in SKIP: assert info.filename not in parts, info.filename; continue
        zo.writestr(info, parts.pop(info.filename) if info.filename in parts else z.read(info.filename))
    for name, data in parts.items(): zo.writestr(name, data, compress_type=zipfile.ZIP_DEFLATED)
report.update(output=str(OUT.relative_to(ROOT)), output_sha256=sha(OUT.read_bytes()), order=spec.ORDER, blocks=len(in_sheets),
              sources_rows=len(recs), parts_written=changed)
REPORT.write_text(json.dumps(report, indent=1, ensure_ascii=False))
print('ok', OUT.name, '| sheets', len(spec.ORDER), '| blocks', len(in_sheets), '| new data rows', sum(x[2][2] for x in sources_new),
      '| parts written', len(changed))
