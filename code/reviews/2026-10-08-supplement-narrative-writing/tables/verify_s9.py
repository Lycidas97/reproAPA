"""Independent readback against frozen files, plus original OOXML preservation checks."""
from pathlib import Path
import collections, csv, gzip, hashlib, json, math, re, zipfile, xml.etree.ElementTree as E
import openpyxl
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
payload=json.loads((HERE/'s9-payload.json').read_text()); logs=json.loads((HERE/'s9-blocks.json').read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(rel):
    p=ROOT/rel
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else p.open()) as f:
        r=csv.DictReader(f,delimiter='\t');rows=list(r);return r.fieldnames,rows
def typed(v):
    if v=='' or v is None:return None
    if isinstance(v,bool):return str(v).lower()
    if isinstance(v,(float,int)):return v
    if re.fullmatch(r'[-+]?\d+',v) and abs(int(v))>2**53-1:return v
    if re.fullmatch(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?',v):return float(v)
    return v
def equal(a,b):
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        return math.isclose(a,b,rel_tol=3e-15,abs_tol=1e-300)
    return a==b
orig='notes/2026-10-08-narrative-figures/opposite-sign/opposite-sign.json'
obj=json.loads((ROOT/orig).read_text())
ep='notes/2026-10-06-supplement-apadetect-p50/contract-r2/bundle/tables/endpoints/route_effects.tsv.gz'
_,endpoints=read(ep)
primary={(r['gene_id'],r['route']):r for r in endpoints if r['setting']=='supported50_p08' and r['design']=='same_library' and r['contrast']=='P_M' and r['route'] in ['CIP','WOT']}
def derived(name):
    out=[]
    for symbol,g in obj.items():
        if name=='opposite-sign-decomposition.tsv':
            for x in g['decomposition']:
                r=dict(gene_symbol=symbol,gene_id=g['gene_id'],pas_source='polyapipe',setting='supported50_p08',design='same_library',dwui_pm_recurrent=g['dwui_pm_recurrent'],**x)
                r['family_start0'],r['family_end0']=r.pop('family')
                for k in ['terminal_shift_nt','within','mass','outside']:r.setdefault(k,'NA')
                out.append(r)
        elif name=='opposite-sign-intervals.tsv':
            for iv in g['strict_intervals']:
                for route in ['CIP','WOT']:
                    for cutoff in ['20','50']:
                        x=iv[route][cutoff]
                        out.append(dict(gene_symbol=symbol,gene_id=g['gene_id'],pas_source='polyapipe',route=route,interval_start0=iv['interval'][0],interval_end0=iv['interval'][1],interval_length_nt=iv['length'],family_UMI_cutoff=int(cutoff),n_paired_libraries=x['n'],dI_pp=x['dI_pp']))
        elif name=='opposite-sign-p50-primary.tsv':
            for route in ['CIP','WOT']:
                x=primary.get((g['gene_id'],route));p=g['p50']
                out.append(dict(gene_symbol=symbol,gene_id=g['gene_id'],pas_source='apadetect-p50',route=route,setting='supported50_p08',design='same_library',delta_DWUI=float(x['delta_DWUI']) if x else 'NA',n_pairs=int(x['n_pairs']) if x else 'NA',dwui_pm_recurrent=p['dwui_pm_recurrent'],interval_pairs=p['interval_pairs'],strict_pairs=p['strict_pairs'],local_evaluable='yes' if p['interval_pairs'] else 'no_interval_pair'))
        else:raise AssertionError(name)
    return out
expected={}; hashes={}; na=collections.Counter(); blanks=collections.Counter(); totalcells=0
for b in payload['blocks']:
    for p in b['inputs']:
        assert sha(ROOT/p['path'])==p['sha256'],p['path'];hashes[p['path']]=p['sha256']
    if '/tables/data/' in b['source']:rows=derived(Path(b['source']).name)
    else:
        _,rows=read(b['source'])
        if b['source'].endswith('corrected-tests.tsv') and b['sheet'] in ['S9a','S9b','S9d']:
            families={'S9a':lambda r:r['family'].startswith('A1-'),'S9d':lambda r:r['family']=='A6-shape-endpoint'}
            if b['sheet']=='S9b':
                family='A2-cycle-scores' if 'plan-a-analysis/' in b['source'] else 'A2-cycle-alignment-redesign'
                rows=[r for r in rows if r['family']==family]
            else:rows=[r for r in rows if families[b['sheet']](r)]
    assert len(rows)==b['rows'],b['block']
    assert all(len(str(r.get(c,'')))<=32767 for r in rows for c in b['columns'])
    matrix=[[typed(r[c]) for c in b['columns']] for r in rows]
    # Independently compare the intermediate payload too, rather than trusting it as ground truth.
    for actual,wanted in zip(b['values'],matrix):
        assert all(equal(typed(a),w) for a,w in zip(actual,wanted)),b['block']
    expected[b['block']]=matrix
    na[b['sheet']]+=sum(v=='NA' for row in matrix for v in row)
    blanks[b['sheet']]+=sum(v is None for row in matrix for v in row)
    totalcells+=sum(map(len,matrix))
workbook=HERE/'Supplementary_Tables_S1-S9.xlsx'; baseline=ROOT/payload['baseline']
N='http://schemas.openxmlformats.org/spreadsheetml/2006/main';R='http://schemas.openxmlformats.org/officeDocument/2006/relationships';NS={'s':N}
with zipfile.ZipFile(baseline) as old,zipfile.ZipFile(workbook) as new:
    oldwb=E.fromstring(old.read('xl/workbook.xml'));newwb=E.fromstring(new.read('xl/workbook.xml'))
    oldrels={r.get('Id'):r.get('Target') for r in E.fromstring(old.read('xl/_rels/workbook.xml.rels'))}
    oldsheet={s.get('name'):'xl/'+oldrels[s.get('{'+R+'}id')] for s in oldwb.find('s:sheets',NS)}
    names=[s.get('name') for s in newwb.find('s:sheets',NS)]
    assert names==list(oldsheet)[:-1]+list(payload['titles'])+['Sources'],names
    preserved=[]
    for name,path in oldsheet.items():
        if name not in ['Contents','Sources']:
            assert old.read(path)==new.read(path),name
            rel=Path(path).parent/'_rels'/(Path(path).name+'.rels')
            if str(rel) in old.namelist():assert old.read(str(rel))==new.read(str(rel)),name+'rels'
            preserved.append(name)
    assert len(preserved)==43
    assert old.read('xl/sharedStrings.xml')==new.read('xl/sharedStrings.xml')
    assert old.read('xl/theme/theme1.xml')==new.read('xl/theme/theme1.xml')
    oldstyles=E.fromstring(old.read('xl/styles.xml'));newstyles=E.fromstring(new.read('xl/styles.xml'))
    style_checks={}
    for group in ['numFmts','fonts','fills','borders','cellStyleXfs','cellXfs','cellStyles','dxfs']:
        a=oldstyles.find('s:'+group,NS);b=newstyles.find('s:'+group,NS)
        if a is None:continue
        assert b is not None
        for i,item in enumerate(a):assert E.tostring(item)==E.tostring(b[i]),(group,i)
        style_checks[group]=len(a)
    # Every original table definition is retained except Sources, whose range is extended.
    oldtables=0
    for p in old.namelist():
        if re.fullmatch(r'xl/tables/table\d+\.xml',p):
            tree=E.fromstring(old.read(p))
            if tree.get('name')!='T_Sources':assert old.read(p)==new.read(p),p;oldtables+=1
    newrels={r.get('Id'):r.get('Target') for r in E.fromstring(new.read('xl/_rels/workbook.xml.rels'))}
    newpaths={'xl/'+newrels[s.get('{'+R+'}id')] for s in newwb.find('s:sheets',NS) if s.get('name') in payload['titles']}
    for p in newpaths:
            tree=E.fromstring(new.read(p));assert not list(tree.iter('{'+N+'}f')),p
            assert not any(c.get('t')=='e' for c in tree.iter('{'+N+'}c')),p
            assert not tree.find('.//s:c[@r="A3"]/s:f',NS)
    # S1-S8 provenance/contents data cells are unchanged; allowed top metadata cells only.
    preservation_metadata={}
    for name,allowed in [('Contents',{'A1','A2','C5'}),('Sources',{'A2'})]:
        a=E.fromstring(old.read(oldsheet[name]));b=E.fromstring(new.read(oldsheet[name]));bm={c.get('r'):c for c in b.iter('{'+N+'}c')}
        count=0
        for c in a.iter('{'+N+'}c'):
            if c.get('r') not in allowed:assert E.tostring(c)==E.tostring(bm[c.get('r')]),(name,c.get('r'));count+=1
        preservation_metadata[name]=count
    allnames=[]
    for p in new.namelist():
        if re.fullmatch(r'xl/tables/table\d+\.xml',p):allnames.append(E.fromstring(new.read(p)).get('name'))
    assert len(allnames)==len(set(allnames))
    formats={0:'General',1:'0',2:'0.00'}
    formats.update({int(e.get('numFmtId')):e.get('formatCode') for e in newstyles.find('s:numFmts',NS)})
    sheetpaths={s.get('name'):'xl/'+newrels[s.get('{'+R+'}id')] for s in newwb.find('s:sheets',NS)}
    format_samples=[]
    for name,ref,negative in [('S9e','I10',True),('S9d','E9',True),('S9e','G25',False),('S9e','U11',False)]:
        tree=E.fromstring(new.read(sheetpaths[name]));c=tree.find('.//s:c[@r="'+ref+'"]',NS)
        value=float(c.find('s:v',NS).text);xf=newstyles.find('s:cellXfs',NS)[int(c.get('s'))];code=formats[int(xf.get('numFmtId','0'))]
        assert value<0 if negative else value>0
        assert code in ['0','0.0000','0.000E+00'],(name,ref,code)
        if ref=='U11':assert code=='0'
        format_samples.append({'sheet':name,'cell':ref,'stored_value':value,'numFmt':code})
wb=openpyxl.load_workbook(workbook,read_only=True,data_only=False)
max_difference=0.0;readback_cells=0
for sheet in payload['titles']:
    sh=wb[sheet];checks={}
    for l in logs:
        if l['sheet']==sheet:
            b=next(b for b in payload['blocks'] if b['block']==l['block'])
            checks[l['header_row']]=(l,b,0)
            for i,row in enumerate(expected[l['block']],1):checks[l['header_row']+i]=(l,row,i)
    for rownum,row in enumerate(sh.iter_rows(values_only=True),1):
        if rownum not in checks:continue
        l,wanted,i=checks[rownum];wanted=wanted['columns'] if i==0 else wanted
        for j,value in enumerate(wanted):
            actual=row[j];assert equal(actual,value),(sheet,l['block'],rownum,j,actual,value)
            if i:readback_cells+=1
            if isinstance(value,(int,float)) and isinstance(actual,(int,float)):max_difference=max(max_difference,abs(actual-value))
    assert len(checks)==sum(l['rows']+1 for l in logs if l['sheet']==sheet)
wb.close()
assert readback_cells==totalcells
# Final result families and intended availability boundaries.
def source_rows(fragment):
    b=next(b for b in payload['blocks'] if fragment in b['source']);return [dict(zip(b['columns'],r)) for r in expected[b['block']]]
a1=source_rows('plan-a-redesign/corrected-tests.tsv');assert len(a1)==12
assert collections.Counter(r['family'] for r in a1)=={'A1-conditional-reference-redesign':8,'A1-direction-composition-redesign':4}
a6=source_rows('A6-measurable-endpoint-recurrence.tsv');assert len(a6)==6
score=next(expected[b['block']] for b in payload['blocks'] if b['sheet']=='S9b' and 'plan-a-analysis/corrected-tests.tsv' in b['source']);assert len(score)==102
cy=next((b,expected[b['block']]) for b in payload['blocks'] if b['sheet']=='S9b' and 'plan-a-redesign/corrected-tests.tsv' in b['source'])
cyrows=[dict(zip(cy[0]['columns'],r)) for r in cy[1]];assert len(cyrows)==16
assert all(r['raw_p']=='NA' and r['q_BH_family']=='NA' and r['n_libraries']==0 for r in cyrows)
assert all(not any(x in b['source'] for x in ['pm_route_means','pm_route_mean_tests','A5-','A7-']) for b in payload['blocks'])
p50=source_rows('opposite-sign-p50-primary.tsv')
assert all(r['setting']=='supported50_p08' and r['design']=='same_library' for r in p50)
key={(r['gene_symbol'],r['route']):r for r in p50}
assert equal(key['Riok3','CIP']['delta_DWUI'],0.000072688531097783379)
assert equal(key['Riok3','WOT']['delta_DWUI'],-0.0065759460357852229)
assert equal(key['Fgfr1op2','WOT']['delta_DWUI'],0.015410781814164706)
decomp=source_rows('opposite-sign-decomposition.tsv');fg=[r for r in decomp if r['gene_symbol']=='Fgfr1op2'];assert len(fg)==2
for r in fg:assert math.isclose(sum(r[k] for k in ['within','mass','outside']),r['dDWUI'],abs_tol=1e-12)
iv=source_rows('opposite-sign-intervals.tsv');assert {r['interval_length_nt'] for r in iv if r['gene_symbol']=='Fgfr1op2'}=={970,492}
assert sha(baseline)==payload['baseline_sha256']
result={'status':'PASS','workbook':str(workbook.relative_to(ROOT)),'sha256':sha(workbook),'baseline_sha256_unchanged':sha(baseline),'baseline_sheets':len(oldsheet),'baseline_lettered_sheets':len(preserved),'original_data_sheets_byte_identical':preserved,'original_table_definitions_byte_identical':oldtables,'original_style_prefixes_preserved':style_checks,'original_contents_sources_cells_preserved':preservation_metadata,'sheets':len(names),'new_blocks':len(logs),'new_rows':sum(b['rows'] for b in payload['blocks']),'data_cells_verified_against_independently_reread_sources':readback_cells,'maximum_absolute_numeric_export_difference':max_difference,'NA_text_cells':dict(na),'blank_cells':dict(blanks),'source_files_SHA_verified':len(hashes),'no_new_cell_formulas_or_error_cells':True,'targeted_primary_p50_corrected':True,'Fgfr1op2_two_intervals_and_decomposition_preserved':True,'cycle16_unavailable_not_zero':True,'A1_two_models_two_PAS_separate12tests':True,'A6_all6classes_and4tests':True,'unused_result_branches_excluded':True}
(HERE/'workbook-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
(HERE/'final-format-inspection.json').write_text(json.dumps({'status':'PASS','negative_signs_and_integer_coordinates':format_samples},indent=2)+'\n')
