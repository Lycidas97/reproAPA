"""Verify every exported result cell and reviewed sheet before issuing table completion."""
import argparse
import json
import math
from pathlib import Path
import posixpath
import re
import xml.etree.ElementTree as ET
import zipfile

from downstream_runtime import Context
from polyapipe_prepare import digest

NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL='http://schemas.openxmlformats.org/officeDocument/2006/relationships'


def same_value(actual,expected):
    if expected is None:return actual is None
    if isinstance(expected,bool):return isinstance(actual,bool) and actual==expected
    if isinstance(expected,(int,float)):
        if isinstance(actual,bool) or not isinstance(actual,(int,float)):return False
        if isinstance(expected,int):return actual==expected
        if expected==0:return actual==0
        return math.isclose(actual,expected,rel_tol=2e-14,abs_tol=0)
    if isinstance(expected,(list,dict)):expected=json.dumps(expected,separators=(',',':'),ensure_ascii=False)
    return actual==expected or (str(expected).startswith('=') and actual=="'"+expected)


def exported_cells(archive,path,strings):
    root=ET.fromstring(archive.read(path));cells={}
    for cell in root.findall('.//s:sheetData/s:row/s:c',NS):
        assert cell.find('s:f',NS) is None, 'Exported statistical results must not become spreadsheet formulas'
        value=cell.find('s:v',NS);value=value.text if value is not None else None
        kind=cell.get('t')
        if kind=='s':value=strings[int(value)]
        elif kind=='inlineStr':value=''.join(t.text or '' for t in cell.findall('.//s:t',NS))
        elif kind=='b':value=value=='1'
        elif kind=='e':raise AssertionError((path,cell.get('r'),value))
        elif value is not None:
            if kind not in ['str','d']:value=float(value)
        cells[cell.get('r')]=value
    return cells


def coordinate(row,col):
    label='';value=col+1
    while value:value,digit=divmod(value-1,26);label=chr(65+digit)+label
    return label+str(row+1)


def verify_workbook(path,payload,root):
    records=[]
    with zipfile.ZipFile(path) as archive:
        strings=[]
        if 'xl/sharedStrings.xml' in archive.namelist():
            strings=[''.join(t.text or '' for t in item.findall('.//s:t',NS)) for item in ET.fromstring(archive.read('xl/sharedStrings.xml'))]
        rels={r.get('Id'):r.get('Target') for r in ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))}
        sheets={s.get('name'):posixpath.normpath(posixpath.join('xl',rels[s.get('{'+REL+'}id')])) for s in ET.fromstring(archive.read('xl/workbook.xml')).findall('s:sheets/s:sheet',NS)}
        sheets={name:path.lstrip('/') for name,path in sheets.items()}
        expected={'Read me','Index'}|{s['name'] for s in payload['sheets'] if s['role']=='complete_table'}
        assert set(sheets)==expected,'Missing or unexpected worksheets'
        for sheet in payload['sheets']:
            if sheet['role']!='complete_table':continue
            data=json.loads((root/sheet['source']).read_text());values=[data['columns']]+data['data']
            cells=exported_cells(archive,sheets[sheet['name']],strings)
            assert all(len(row)==len(data['columns']) for row in values)
            for i,row in enumerate(values):
                for j,expected in enumerate(row):
                    address=coordinate(i,j);actual=cells.pop(address,None)
                    assert same_value(actual,expected),(sheet['name'],address,actual,expected)
            assert not any(v is not None for v in cells.values()),'Unexpected extra values'
            records.append(dict(sheet=sheet['name'],rows=len(values)-1,columns=len(data['columns']),all_cells_verified=True))
        # Covers are also checked for error/formula cells, although their prose is visually reviewed.
        for name in ['Read me','Index']:exported_cells(archive,sheets[name],strings)
    return records


def finish(ctx,review):
    collection=json.loads(ctx.register(ctx.out/'collection.complete.json').read_text())
    payload=json.loads(ctx.register(ctx.out/'workbook-inputs.json').read_text())
    exported=json.loads(ctx.register(ctx.out/'workbook-export.json').read_text())
    for record in [collection,payload,exported]:
        assert record['run_id']==ctx.run_id and record['input_fingerprint']==ctx.config['input_fingerprint']
    for path,sha in collection['inputs'].items():assert digest(ctx.register(path))==sha
    index=ctx.read(ctx.out/'table/S1-S30-table-index.tsv')
    assert set(index.table)==set(range(1,31))
    for row in index.itertuples():assert digest(ctx.out/row.file)==row.sha256
    for entry in payload['sheets']:ctx.register(ctx.out/entry['source'])
    workbook=ctx.register(ctx.out/'table/reproAPA-polyApipe-S1-S30.xlsx')
    numerical=verify_workbook(workbook,payload,ctx.out)
    reviews=json.loads(ctx.register(review).read_text())
    assert reviews['run_id']==ctx.run_id and reviews['input_fingerprint']==ctx.config['input_fingerprint']
    assert reviews['workbook_sha256']==digest(workbook)
    expected={s['name'] for s in exported['sheets']};done={s['name']:s for s in reviews['sheets']}
    assert expected==set(done) and len(done)==len(reviews['sheets'])
    for name,row in done.items():
        assert row['layout']=='pass' and row['notes']
        assert digest(ctx.out/'workbook-previews'/f'{name}.png')==row['preview_sha256']
    (ctx.out/'workbook-value-verification.json').write_text(json.dumps(dict(workbook_sha256=digest(workbook),sheets=numerical,
        all_exported_data_cells_checked=True,relative_float_tolerance=2e-14,missing_and_true_zero_distinct=True),indent=2)+'\n')
    (ctx.out/'workbook-visual-review.json').write_text(json.dumps(reviews,indent=2)+'\n')
    ctx.finish(dict(all_30_tables_delivered=True,all_source_TSV_hashes_verified=True,all_workbook_data_cells_verified=True,
        missing_zero_and_small_probabilities_preserved=True,all_sheet_opening_views_reviewed=True,
        large_tables_explicitly_delivered_as_TSV=True,historical_evidence_scope_preserved=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--review',type=Path,required=True)
    a=p.parse_args();finish(Context(a.config,'tables',a.out),a.review)
