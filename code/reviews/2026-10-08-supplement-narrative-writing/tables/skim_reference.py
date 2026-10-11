"""Make a small, read-only rendering copy of the original native sheet layout."""
from pathlib import Path
import argparse, copy, re, zipfile, xml.etree.ElementTree as E
ROOT=Path(__file__).resolve().parents[3]; HERE=Path(__file__).resolve().parent
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
P='http://schemas.openxmlformats.org/package/2006/relationships'
C='http://schemas.openxmlformats.org/package/2006/content-types'
ns={'s':N,'r':R}; E.register_namespace('',N); E.register_namespace('r',R)
ap=argparse.ArgumentParser();ap.add_argument('--final',action='store_true');args=ap.parse_args()
src=(HERE/'Supplementary_Tables_S1-S9.xlsx') if args.final else ROOT/'notes/2026-10-06-supplement-tables/Supplementary_Tables_S1-S8.xlsx'
selected=['Contents','Sources','S9a','S9b','S9c','S9d','S9e'] if args.final else ['Contents','S4e']
with zipfile.ZipFile(src) as z:
    shared=E.fromstring(z.read('xl/sharedStrings.xml'))
    wb=E.fromstring(z.read('xl/workbook.xml')); sheets=wb.find('s:sheets',ns)
    rels=E.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    targets={r.get('Id'):r.get('Target') for r in rels}
    chosen=[s for s in list(sheets) if s.get('name') in selected]
    for s in list(sheets):
        if s not in chosen:sheets.remove(s)
    dn=wb.find('s:definedNames',ns)
    if dn is not None:wb.remove(dn)
    keep={s.get('{'+R+'}id') for s in chosen}
    for r in list(rels):
        if r.get('Type').endswith('/worksheet') and r.get('Id') not in keep:rels.remove(r)
        if r.get('Type').endswith('/sharedStrings'):rels.remove(r)
    parts={}
    for s in chosen:
        p='xl/'+targets[s.get('{'+R+'}id')]
        tree=E.fromstring(z.read(p)); sd=tree.find('s:sheetData',ns)
        cap=(107 if s.get('name')=='S9e' else 1000 if s.get('name') in ['Contents','Sources'] else 20) if args.final else 20
        for row in list(sd):
            if int(row.get('r'))>cap:sd.remove(row);continue
            for cell in row:
                if cell.get('t')=='s':
                    v=cell.find('s:v',ns);si=shared[int(v.text)];cell.remove(v);cell.set('t','inlineStr')
                    inline=E.SubElement(cell,'{'+N+'}is')
                    for t in si:inline.append(copy.deepcopy(t))
        for tag in ['tableParts','hyperlinks','drawing','legacyDrawing','extLst']:
            el=tree.find('s:'+tag,ns)
            if el is not None:tree.remove(el)
        last=max(int(row.get('r')) for row in sd)
        col=max((re.match('[A-Z]+',c.get('r')).group(0) for row in sd for c in row),key=lambda x:(len(x),x))
        dim=tree.find('s:dimension',ns)
        if dim is None:dim=E.Element('{'+N+'}dimension');tree.insert(0,dim)
        dim.set('ref',f'A1:{col}{last}')
        parts[p]=E.tostring(tree,encoding='utf-8',xml_declaration=True)
    ct=E.fromstring(z.read('[Content_Types].xml'))
    for el in list(ct):
        p=el.get('PartName','').lstrip('/')
        if p.startswith('xl/worksheets/') and p not in parts:ct.remove(el)
        if p.startswith(('xl/tables/','xl/sharedStrings')):ct.remove(el)
    parts['xl/workbook.xml']=E.tostring(wb,encoding='utf-8',xml_declaration=True)
    E.register_namespace('',P)
    parts['xl/_rels/workbook.xml.rels']=E.tostring(rels,encoding='utf-8',xml_declaration=True)
    E.register_namespace('',C)
    parts['[Content_Types].xml']=E.tostring(ct,encoding='utf-8',xml_declaration=True)
    for p in ['xl/styles.xml','xl/theme/theme1.xml','_rels/.rels','docProps/core.xml','docProps/app.xml']:
        parts[p]=z.read(p)
    target=HERE/('final-skim.xlsx' if args.final else 'reference-skim.xlsx')
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as out:
        for p,v in parts.items():out.writestr(p,v)
print(target)
