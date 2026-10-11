"""Verify every archived byte stream and the exact filtered-ledger exception."""
from pathlib import Path
import csv, hashlib, io, json, zipfile
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def digest(stream):
    h=hashlib.sha256()
    while data:=stream.read(2**20):h.update(data)
    return h.hexdigest()
receipt=json.loads((HERE/'source-archive-verification.json').read_text());path=ROOT/receipt['archive']
with path.open('rb') as f:assert digest(f)==receipt['sha256']
with zipfile.ZipFile(path) as z:
    rows=list(csv.DictReader(io.StringIO(z.read('FILES-SHA256.tsv').decode()),delimiter='\t'))
    assert set(z.namelist())=={r['path'] for r in rows}|{'README.txt','FILES-SHA256.tsv'}
    subsets={}
    for r in rows:
        with z.open(r['path']) as f:assert digest(f)==r['archive_sha256'],r['path']
        with (ROOT/r['path']).open('rb') as f:assert digest(f)==r['original_sha256'],r['path']
        assert z.getinfo(r['path']).file_size==int(r['bytes'])
        assert not any(x in r['path'] for x in ['/initial-','/A5-','/A7-','pm_route_means','pm_route_mean_tests'])
        if r['archive_sha256']!=r['original_sha256']:
            data=list(csv.DictReader(io.StringIO(z.read(r['path']).decode()),delimiter='\t'))
            allowed={'A1-conditional-reference-redesign','A1-direction-composition-redesign','A2-cycle-alignment-redesign'} if 'redesign/' in r['path'] else {'A2-cycle-scores','A6-shape-endpoint'}
            assert {x['family'] for x in data}==allowed
            orig=list(csv.DictReader((ROOT/r['path']).open(),delimiter='\t'))
            assert data==[x for x in orig if x['family'] in allowed]
            subsets[r['path']]={'rows':len(data),'families':sorted(allowed)}
    assert len(subsets)==2
result={'status':'PASS','archive':receipt['archive'],'sha256':receipt['sha256'],'member_files':len(rows)+2,'all_manifest_SHA_and_ZIP_CRC_verified':True,'all_original_source_SHA_verified':True,'exact_filtered_ledger_copies':subsets,'unused_A5_A7_A8_result_products_excluded':True}
(HERE/'archive-content-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
