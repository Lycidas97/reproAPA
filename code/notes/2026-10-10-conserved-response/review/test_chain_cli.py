"""Test standalone production chain CLI on local fake chains, no algorithm import.
Expected coordinates/statuses below are hand-defined independently of the parser.
"""
from pathlib import Path
import gzip,hashlib,json,subprocess,sys
O=Path(__file__).resolve().parent;S=O.parent/'sites'
chain='''chain 100 chrP 1000 + 100 130 chrQ 2000 + 200 230 1
10 0 0
20

chain 100 chrG 1000 + 100 125 chrQ 2000 + 200 225 2
10 5 5
10

chain 100 chrI 1000 + 100 120 chrQ 2000 + 200 225 3
10 0 5
10

chain 100 chrN 1000 + 100 130 chrR 2000 - 200 230 4
30

chain 100 chrD 1000 + 100 130 chrD1 2000 + 200 230 5
30

chain 90 chrD 1000 + 100 130 chrD2 2000 + 300 330 6
30
'''
p=O/'toy-chain.chain.gz'
with gzip.GzipFile(filename=str(p),mode='wb',mtime=0) as f:f.write(chain.encode())
tests=[('contiguous_8','chrP','+',108,116,'unique',('chrQ','+',208,216)),('plus_7','chrP','+',101,108,'unique',('chrQ','+',201,208)),('gene_minus_8','chrP','-',108,116,'unique',('chrQ','-',208,216)),('right_boundary_8','chrP','+',122,130,'unique',('chrQ','+',222,230)),('outside_boundary_8','chrP','+',123,131,'gap_or_partial',None),('source_gap_8','chrG','+',108,116,'gap_or_partial',None),('query_gap_8','chrI','+',108,116,'gap_or_partial',None),('inversion_plus_8','chrN','+',105,113,'unique',('chrR','-',1787,1795)),('inversion_minus_7','chrN','-',105,112,'unique',('chrR','+',1788,1795)),('multichain_7','chrD','+',106,113,'ambiguous',None),('unmapped_7','chrZ','+',100,107,'unmapped',None)]
queries=[dict(query_id=n,species='human',chrom=c,strand=s,start0=a,end0=b) for n,c,s,a,b,_,_ in tests]
resources={k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k in ['chain_hg38_mm39','chain_mm39_hg38']}
code='REMOTE_RESOURCES='+repr(resources)+'\nREMOTE_QUERIES='+repr(queries)+'\n'+(S/'chain_mapper.py').read_text()
r=subprocess.run([sys.executable,'-'],input=code,text=True,capture_output=True,check=True)
out={x['query_id']:x for x in map(json.loads,r.stdout.splitlines())}
for n,c,s,a,b,status,coord in tests:
 x=out[n];assert x['status']==status,(n,x,status)
 assert all(h['end0']-h['start0']==b-a for h in x['hits'])
 if coord:
  h=x['hits'][0];assert (h['chrom'],h['strand'],h['start0'],h['end0'])==coord,(n,h,coord)
 if status=='ambiguous':assert len(x['hits'])==2
report={'status':'PASS','tests':len(tests),'method':'Standalone production CLI only for toy tests; expected coordinates/status hand specified. No production algorithm imports or use to calculate independent real-data statistics.','mapper_sha256':hashlib.sha256((S/'chain_mapper.py').read_bytes()).hexdigest(),'checks':out}
(O/'chain_toy_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':'PASS','tests':len(tests)}))
