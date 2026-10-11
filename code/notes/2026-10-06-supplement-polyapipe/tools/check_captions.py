# Check supplementary legends: panel letters match the receipt; every number appears in values.json or receipt notes (within rounding).
# Run from notes/2026-10-06-supplement-polyapipe: python3 tools/check_captions.py
import json,glob,re,sys
order=['D01','D02','D03','G1','G2','G3']+[f'D{i:02d}' for i in range(4,28)]
def nums(o,out):
    if isinstance(o,dict):
        for k,v in o.items(): nums(k,out); nums(v,out)
    elif isinstance(o,list):
        for v in o: nums(v,out)
    elif isinstance(o,bool): pass
    elif isinstance(o,(int,float)): out.add(o)
    elif isinstance(o,str):
        for t in re.findall(r'-?\d+(?:\.\d+)?',o): out.add(float(t))
ALT='../2026-10-06-alt-catalogue-suppfig/final'
for i,f in list(enumerate(order,1))+[(31,'A01')]:
    cdir = ALT if f=='A01' else f'candidates/{f}'
    cap=open(f'captions/S{i}.en.txt').read()
    assert cap.startswith(f'Figure S{i}. '),i
    r=json.load(open(glob.glob(f'{cdir}/*.receipt.json')[0]))
    pan=sorted(p.lower() for p in r['panels'])
    # letters referenced
    refs=set()
    for m in re.finditer(r'\(([a-h])(?:[,–-] ?([a-h]))?\)',cap):
        a,b=m.group(1),m.group(2)
        if b and '–' in m.group(0): refs|={chr(c) for c in range(ord(a),ord(b)+1)}
        else: refs|={a}|({b} if b else set())
    miss=set(pan)-refs; extra=refs-set(pan)
    vals=set(); nums(json.load(open(glob.glob(f'{cdir}/*.values.json')[0])),vals); nums(r.get('notes',[]),vals)
    body=re.sub(r'Figure S\d+\.','',cap)
    body=re.sub(r'Figure \d[a-z]?|Figure S\d+[a-z, ]*','',body)
    unmatched=[]
    for t in re.findall(r'(?<![A-Za-z0-9])[−-]?\d[\d,]*(?:\.\d+)?|(?<![\d])\.\d+',body):
        s=t.replace(',','').replace('−','-')
        try: x=float(s)
        except: continue
        ok=any(abs(abs(v)-abs(x))<=(0.5*10**-(len(s.split('.')[1]) if '.' in s else 0))+1e-9 for v in vals if isinstance(v,(int,float)))
        if not ok: unmatched.append(t)
    print(f'S{i} {f} panels={pan} missing={sorted(miss)} extra={sorted(extra)} words={len(cap.split())}')
    if unmatched: print('   unmatched:',' '.join(dict.fromkeys(unmatched)))
