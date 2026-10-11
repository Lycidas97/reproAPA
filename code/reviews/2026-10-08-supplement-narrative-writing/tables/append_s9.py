"""Append artifact-authored S9 to a copy; never rebuild the original 42 data sheets."""
from pathlib import Path
import copy, hashlib, json, re, zipfile, xml.etree.ElementTree as E

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
P='http://schemas.openxmlformats.org/package/2006/relationships'
C='http://schemas.openxmlformats.org/package/2006/content-types'
NS={'s':N}; S=lambda x:'{'+N+'}'+x
def xml(tree, namespace=N):
    E.register_namespace('',namespace); E.register_namespace('r',R)
    return E.tostring(tree,encoding='utf-8',xml_declaration=True)
def sha(data):return hashlib.sha256(data).hexdigest()
def cell(row, ref, value, style):
    c=E.SubElement(row,S('c'),{'r':ref,'s':str(style)})
    if isinstance(value,(int,float)):E.SubElement(c,S('v')).text=str(value)
    elif value is not None:
        c.set('t','inlineStr');E.SubElement(E.SubElement(c,S('is')),S('t')).text=str(value)
    return c
def put(tree,ref,value):
    c=tree.find('.//s:c[@r="'+ref+'"]',NS)
    assert c is not None,ref
    for e in list(c):c.remove(e)
    c.set('t','inlineStr');E.SubElement(E.SubElement(c,S('is')),S('t')).text=value
def rowstyles(tree,number):
    row=tree.find('.//s:row[@r="'+str(number)+'"]',NS)
    return [int(c.get('s','0')) for c in row]
def append_hyper(tree,ref,location,display):
    links=tree.find('s:hyperlinks',NS)
    if links is None:
        links=E.Element(S('hyperlinks'))
        # Hyperlinks precede print/page settings in spreadsheet schema order.
        ix=next((i for i,e in enumerate(tree) if e.tag in {S(t) for t in ['printOptions','pageMargins','pageSetup','headerFooter','tableParts','extLst']}),len(tree))
        tree.insert(ix,links)
    E.SubElement(links,S('hyperlink'),{'ref':ref,'location':location,'display':display})

payload=json.loads((HERE/'s9-payload.json').read_text()); logs=json.loads((HERE/'s9-blocks.json').read_text())
baseline=ROOT/payload['baseline']; out=HERE/'Supplementary_Tables_S1-S9.xlsx'
assert sha(baseline.read_bytes())==payload['baseline_sha256']
with zipfile.ZipFile(baseline) as old,zipfile.ZipFile(HERE/'S9-native.xlsx') as new:
    parts={}
    os=E.fromstring(old.read('xl/styles.xml')); ns=E.fromstring(new.read('xl/styles.xml'))
    oldnum=os.find(S('numFmts')); newnum=ns.find(S('numFmts'))
    if oldnum is None:oldnum=E.Element(S('numFmts'),{'count':'0'});os.insert(0,oldnum)
    codes={el.get('formatCode'):int(el.get('numFmtId')) for el in oldnum}
    nextnum=max([163]+list(codes.values()))+1; nummap={}
    if newnum is not None:
        for el in newnum:
            code=el.get('formatCode'); idx=int(el.get('numFmtId'))
            if code=='[>=1000000]0.000E+00;[<=-1000000]0.000E+00;0.0000':code='0.0000'
            if code not in codes:
                copyel=copy.deepcopy(el);copyel.set('numFmtId',str(nextnum));oldnum.append(copyel)
                codes[code]=nextnum;nextnum+=1
            nummap[idx]=codes[code]
    oldnum.set('count',str(len(oldnum)))
    offsets={}
    for name in ['fonts','fills','borders','cellStyleXfs','cellXfs','dxfs']:
        dest=os.find(S(name));source=ns.find(S(name))
        if source is None:continue
        if dest is None:
            dest=E.Element(S(name),{'count':'0'})
            ix=next((i for i,e in enumerate(os) if e.tag in {S('tableStyles'),S('extLst')}),len(os));os.insert(ix,dest)
        offsets[name]=len(dest)
        for item in source:
            item=copy.deepcopy(item)
            if name in ['cellStyleXfs','cellXfs']:
                for attr,group in [('fontId','fonts'),('fillId','fills'),('borderId','borders')]:
                    if attr in item.attrib:item.set(attr,str(int(item.get(attr))+offsets[group]))
                if 'xfId' in item.attrib:item.set('xfId',str(int(item.get('xfId'))+offsets.get('cellStyleXfs',0)))
                if int(item.get('numFmtId','0')) in nummap:item.set('numFmtId',str(nummap[int(item.get('numFmtId'))]))
            if name=='dxfs':
                for nf in item.iter(S('numFmt')):
                    if int(nf.get('numFmtId','0')) in nummap:nf.set('numFmtId',str(nummap[int(nf.get('numFmtId'))]))
            dest.append(item)
        dest.set('count',str(len(dest)))
    parts['xl/styles.xml']=xml(os)
    shared=E.fromstring(new.read('xl/sharedStrings.xml')) if 'xl/sharedStrings.xml' in new.namelist() else None
    owb=E.fromstring(old.read('xl/workbook.xml')); orel=E.fromstring(old.read('xl/_rels/workbook.xml.rels'))
    nwb=E.fromstring(new.read('xl/workbook.xml')); nrel=E.fromstring(new.read('xl/_rels/workbook.xml.rels'))
    rels={e.get('Id'):e.get('Target') for e in nrel}
    oldrels={e.get('Id'):e.get('Target') for e in orel}
    sheets=owb.find(S('sheets')); source_sheet=next(s for s in sheets if s.get('name')=='Sources')
    source_path='xl/'+oldrels[source_sheet.get('{'+R+'}id')]
    linkstyle=E.fromstring(old.read(source_path)).find('.//s:c[@r="A3"]',NS).get('s')
    source_index=list(sheets).index(source_sheet)
    nextsheet=max(int(s.get('sheetId')) for s in sheets)+1
    nextrid=max(int(e.get('Id').removeprefix('rId')) for e in orel)+1
    oldtable=max(int(re.search(r'table(\d+)\.xml$',p).group(1)) for p in old.namelist() if re.match(r'xl/tables/table\d+\.xml$',p))
    nexttable=oldtable+1; newparts=[]
    special_styles={}
    for added in nwb.find(S('sheets')):
        name=added.get('name'); spath='xl/'+rels[added.get('{'+R+'}id')].lstrip('/')
        if spath.startswith('xl/xl/'):spath=spath[3:]
        tree=E.fromstring(new.read(spath))
        for c in tree.iter(S('c')):
            if 's' in c.attrib:c.set('s',str(int(c.get('s'))+offsets['cellXfs']))
            if c.get('t')=='s':
                v=c.find(S('v'));si=shared[int(v.text)];c.remove(v);c.set('t','inlineStr');inline=E.SubElement(c,S('is'))
                for x in si:inline.append(copy.deepcopy(x))
            elif c.get('t') in [None,'n']:
                v=c.find(S('v'))
                if v is not None:
                    number=float(v.text)
                    if number!=0 and not number.is_integer() and (abs(number)<.001 or abs(number)>=1e6):
                        original_style=int(c.get('s','0'));key=(original_style,codes['0.000E+00'])
                        if key not in special_styles:
                            xf=copy.deepcopy(os.find(S('cellXfs'))[original_style]);xf.set('numFmtId',str(key[1]));xf.set('applyNumberFormat','1')
                            special_styles[key]=len(os.find(S('cellXfs')));os.find(S('cellXfs')).append(xf)
                        c.set('s',str(special_styles[key]))
        # Excel doubles cannot retain every digit of a 64-bit random seed; preserve such integer identifiers as text.
        cells={c.get('r'):c for c in tree.iter(S('c'))}
        for block in payload['blocks']:
            if block['sheet']!=name:continue
            header=next(l['header_row'] for l in logs if l['block']==block['block'])
            for ri,row in enumerate(block['values'],header+1):
                for ci,value in enumerate(row):
                    if isinstance(value,str) and re.fullmatch(r'[-+]?\d+',value) and abs(int(value))>2**53-1:
                        letter='';n=ci+1
                        while n:n,q=divmod(n-1,26);letter=chr(65+q)+letter
                        put(tree,f'{letter}{ri}',value)
        for e in list(tree.iter(S('row')))+list(tree.iter(S('col'))):
            attr='s' if e.tag==S('row') else 'style'
            if attr in e.attrib:e.set(attr,str(int(e.get(attr))+offsets['cellXfs']))
        for rule in tree.iter(S('cfRule')):
            if 'dxfId' in rule.attrib:rule.set('dxfId',str(int(rule.get('dxfId'))+offsets.get('dxfs',0)))
        for view in tree.iter(S('sheetView')):view.set('showGridLines','0')
        # Match the original builder's maximum width across all blocks, including sample data.
        widths={}
        for block in payload['blocks']:
            if block['sheet']!=name:continue
            for ci,heading in enumerate(block['columns'],1):
                sample=[str(row[ci-1]) for row in block['values'][:300] if row[ci-1] is not None]
                widths[ci]=max(widths.get(ci,7),min(40,max(len(heading)*1.1,max([0]+[len(x) for x in sample]))+3))
        cols=tree.find(S('cols'))
        if cols is None:
            cols=E.Element(S('cols'));tree.insert(list(tree).index(tree.find(S('sheetData'))),cols)
        for c in list(cols):cols.remove(c)
        for ci,width in widths.items():E.SubElement(cols,S('col'),{'min':str(ci),'max':str(ci),'width':str(width),'customWidth':'1'})
        put(tree,'A3','← Contents');append_hyper(tree,'A3',"'Contents'!A1",'← Contents')
        tree.find('.//s:c[@r="A3"]',NS).set('s',linkstyle)
        newpath=f'xl/worksheets/sheet{nextsheet}.xml';parts[newpath]=xml(tree);newparts.append(newpath)
        oldnspath=Path(spath); relpath=str(oldnspath.parent/'_rels'/(oldnspath.name+'.rels'))
        if relpath in new.namelist():
            sr=E.fromstring(new.read(relpath))
            for rr in sr:
                assert rr.get('Type').endswith('/table'),rr.attrib
                target=rr.get('Target'); tpath=str((oldnspath.parent/target)).replace('worksheets/../','')
                if target.startswith('/'):tpath=target.lstrip('/')
                table=E.fromstring(new.read(tpath));table.set('id',str(nexttable))
                if 'headerRowDxfId' in table.attrib:table.set('headerRowDxfId',str(int(table.get('headerRowDxfId'))+offsets.get('dxfs',0)))
                for tc in table.iter(S('tableColumn')):
                    for attr in ['dataDxfId','headerRowDxfId','totalsRowDxfId']:
                        if attr in tc.attrib:tc.set(attr,str(int(tc.get(attr))+offsets.get('dxfs',0)))
                tsi=table.find(S('tableStyleInfo'))
                if tsi is not None:
                    tsi.attrib.pop('name',None);tsi.set('showRowStripes','0');tsi.set('showColumnStripes','0')
                tdest=f'xl/tables/table{nexttable}.xml';parts[tdest]=xml(table);newparts.append(tdest)
                rr.set('Target',f'../tables/table{nexttable}.xml');nexttable+=1
            parts[f'xl/worksheets/_rels/sheet{nextsheet}.xml.rels']=xml(sr,P)
        s=E.Element(S('sheet'),{'name':name,'sheetId':str(nextsheet),'{'+R+'}id':f'rId{nextrid}'})
        sheets.insert(list(sheets).index(source_sheet),s)
        E.SubElement(orel,'{'+P+'}Relationship',{'Id':f'rId{nextrid}','Type':R+'/worksheet','Target':f'worksheets/sheet{nextsheet}.xml'})
        nextsheet+=1;nextrid+=1
    for dn in owb.findall('s:definedNames/s:definedName',NS):
        if dn.get('localSheetId')==str(source_index):dn.set('localSheetId',str(source_index+5))
    os.find(S('cellXfs')).set('count',str(len(os.find(S('cellXfs')))))
    parts['xl/styles.xml']=xml(os)
    parts['xl/workbook.xml']=xml(owb);parts['xl/_rels/workbook.xml.rels']=xml(orel,P)
    ct=E.fromstring(old.read('[Content_Types].xml'))
    for p in newparts:
        kind='worksheet' if '/worksheets/' in p else 'table'
        E.SubElement(ct,'{'+C+'}Override',{'PartName':'/'+p,'ContentType':f'application/vnd.openxmlformats-officedocument.spreadsheetml.{kind}+xml'})
    parts['[Content_Types].xml']=xml(ct,C)

    contents_path='xl/'+oldrels[next(s for s in sheets if s.get('name')=='Contents').get('{'+R+'}id')]
    contents=E.fromstring(old.read(contents_path));sd=contents.find(S('sheetData'))
    put(contents,'A1','Supplementary Tables S1–S9')
    put(contents,'A2','S1–S8 values and formats are preserved from the frozen polyApipe v2-primary workbook. S9 copies later specified analyses from both PAS workflows without statistical recomputation. Every block is filterable, with full stored precision; NA and missing cells remain distinct from zero. Sources records exact filters and file SHA-256. S9 companion sources and code are archived separately from the existing S1–S8 source package.')
    put(contents,'C5','Related figures or analyses')
    r=max(int(row.get('r')) for row in sd)+1
    groupstyle=rowstyles(contents,6)
    row=E.SubElement(sd,S('row'),{'r':str(r)})
    for j,v in enumerate(['Table S9','Analyses specified after the initial results',None,None,None]):cell(row,f'{chr(65+j)}{r}',v,groupstyle[j])
    for i,(name,title) in enumerate(payload['titles'].items()):
        r+=1;row=E.SubElement(sd,S('row'),{'r':str(r)})
        styles=rowstyles(contents,7+(i%2));bs=[b for b in payload['blocks'] if b['sheet']==name]
        for j,v in enumerate([name,title,payload['analyses'][name],len(bs),sum(b['rows'] for b in bs)]):cell(row,f'{chr(65+j)}{r}',v,styles[j])
        append_hyper(contents,f'A{r}',f"'{name}'!A1",name)
    contents.find(S('dimension')).set('ref',f'A1:E{r}');parts[contents_path]=xml(contents)
    sources=E.fromstring(old.read(source_path));sd=sources.find(S('sheetData'));r=max(int(row.get('r')) for row in sd)
    firstnew=r+1
    put(sources,'A2','S1–S8 paths and hashes retain the original bundle-relative convention. S9 paths are project-relative; filtered source files retain original full-file SHA-256 here. The S9 companion ZIP manifest also records each archived filtered file, original hash, exact row filter, and code hash.')
    for i,b in enumerate(payload['blocks']):
        r+=1;row=E.SubElement(sd,S('row'),{'r':str(r)})
        styles=rowstyles(sources,6+((firstnew-6+i)%2))
        rule=b['filter']+(' '+b['notes'] if b['notes'] else '')
        hashes='; '.join(v['sha256'] for v in b['inputs'])
        vals=[b['block'],b['source'],b['rows'],b['source_rows'],rule,hashes]
        for j,v in enumerate(vals):cell(row,f'{chr(65+j)}{r}',v,styles[j])
        append_hyper(sources,f'A{r}',f"'{b['sheet']}'!A1",b['block'])
    sources.find(S('dimension')).set('ref',f'A1:F{r}');parts[source_path]=xml(sources)
    for p in old.namelist():
        if re.match(r'xl/tables/table\d+\.xml$',p):
            table=E.fromstring(old.read(p))
            if table.get('name')=='T_Sources':
                table.set('ref',f'A5:F{r}');table.find(S('autoFilter')).set('ref',f'A5:F{r}');parts[p]=xml(table)
    # Update only workbook metadata; original content sheets and relationships stay byte-identical.
    A='http://schemas.openxmlformats.org/officeDocument/2006/extended-properties';V='http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes'
    app=E.fromstring(old.read('docProps/app.xml'))
    titles=app.find('{'+A+'}TitlesOfParts/{'+V+'}vector')
    if titles is not None:
        for e in list(titles):titles.remove(e)
        for s in sheets:E.SubElement(titles,'{'+V+'}lpstr').text=s.get('name')
        titles.set('size',str(len(sheets)))
    pairs=app.find('{'+A+'}HeadingPairs/{'+V+'}vector')
    if pairs is not None:
        ints=list(pairs.iter('{'+V+'}i4'))
        if ints:ints[0].text=str(len(sheets))
    parts['docProps/app.xml']=xml(app,A)
    core=E.fromstring(old.read('docProps/core.xml'));title=core.find('{http://purl.org/dc/elements/1.1/}title')
    if title is not None:title.text='Supplementary Tables S1–S9'
    parts['docProps/core.xml']=xml(core,'http://schemas.openxmlformats.org/package/2006/metadata/core-properties')
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as dest:
        for item in old.infolist():dest.writestr(item,parts.pop(item.filename,old.read(item.filename)))
        for p,v in parts.items():dest.writestr(p,v)
print(json.dumps({'workbook':str(out),'sheets':len(sheets),'new_blocks':len(logs),'new_rows':sum(b['rows'] for b in payload['blocks']),'sha256':sha(out.read_bytes())},indent=2))
