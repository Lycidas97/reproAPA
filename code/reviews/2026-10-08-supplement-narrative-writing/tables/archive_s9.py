"""Package frozen S9 inputs and original code; exclude unused result branches."""
from pathlib import Path
import csv, gzip, hashlib, io, json, zipfile
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
payload=json.loads((HERE/'s9-payload.json').read_text())
files={}; rows=[]
def add(path,role='frozen input/output',expected=None,filter_note='All bytes preserved.'):
    p=Path(path); p=p if p.is_absolute() else ROOT/p
    rel=str(p.relative_to(ROOT)); data=p.read_bytes(); digest=hashlib.sha256(data).hexdigest()
    if expected is not None:assert digest==expected,(rel,digest,expected)
    if rel not in files:
        files[rel]=data;rows.append(dict(path=rel,role=role,original_sha256=digest,archive_sha256=digest,filter=filter_note,bytes=len(data)))
def manifest(rel,field=None):
    add(rel,'original provenance')
    obj=json.loads((ROOT/rel).read_text());obj=obj[field] if field else obj
    for p,h in obj.items():
        if isinstance(h,str) and len(h)==64 and Path(p).is_absolute():add(p,'frozen analysis input',h)
for b in payload['blocks']:
    for item in b['inputs']:add(item['path'],'S9 displayed source',item['sha256'])
for rel,field in [
 ('reviews/2026-10-08-plan-a-redesign/background/results/input-audit.json','input_sha256'),
 ('reviews/2026-10-08-plan-a-analysis/cycle/input_sha256.json',None),
 ('reviews/2026-10-08-plan-a-redesign/cycle/input_sha256.json',None),
 ('reviews/2026-10-08-plan-a-redesign/cross-species/input-audit.json','input_sha256'),
 ('reviews/2026-10-08-plan-a-analysis/bridges/input-hashes.json',None)]:manifest(rel,field)
for folder,stems in {
 'reviews/2026-10-08-plan-a-redesign':['PROTOCOL.md','CYCLE-FREEZE-ADDENDUM.md','integrate.py'],
 'reviews/2026-10-08-plan-a-redesign/background':['analyze.py','verify.py','test_model.py','PROTOCOL.zh.md','results/verification.json','results/raw-tests.tsv','freeze-receipt.json'],
 'reviews/2026-10-08-plan-a-redesign/cycle':['run_cycle.py','validate.py','test_numeric.py','diagnose_contexts.py','collect.py','extract_remote.py','PROTOCOL.md','validation.json','freeze-receipt.json','execution-receipt.json','completion-receipt.json','additional_pas_subsets.complete.json'],
 'reviews/2026-10-08-plan-a-redesign/cross-species':['run_cross_species.py','verify.py','PROTOCOL.zh.md','verification.json','result.json','seeds.tsv'],
 'reviews/2026-10-08-plan-a-analysis':['PROTOCOL.md','integrate.py'],
 'reviews/2026-10-08-plan-a-analysis/cycle':['run_cycle.py','validate.py','collect.py','extract_remote.py','collect_additional_scores.py','additional_scores_remote.py','PROTOCOL.md','cycle_score_tests.tsv','validation.json','score_coverage.json'],
 'reviews/2026-10-08-plan-a-analysis/bridges':['run_bridges.py','PROTOCOL.md','run-receipt.json'],
 'reviews/2026-10-08-narrative-logic-chain/netshift-local':['compute.py','README.md'],
 'notes/2026-10-08-narrative-figures/opposite-sign':['derive_opposite.py'],
}.items():
    for stem in stems:add(folder+'/'+stem,'original code/protocol/receipt')
for base in ['notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables','notes/2026-10-06-supplement-apadetect-p50/contract-r2/bundle/tables']:
    for stem in ['endpoints/route_effects.tsv.gz','union_pool/figure5_union_claims.tsv','candidate_screen/recurrent_region_pairs.tsv.gz','recurrence/recurrence_members.tsv.gz']:
        add(base+'/'+stem,'netshift frozen input')
base='notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/'
for stem in ['intervals/annotation/mouse/CIP/models.json','intervals/annotation/mouse/WOT/models.json',
             'state_vectors/route_effects.primary.tsv.gz','state_vectors/PAS_axes.primary.tsv.gz',
             'candidate_screen/recurrent_region_pairs.tsv.gz','recurrence/recurrence_members.tsv.gz']:
    add(base+stem,'opposite-sign frozen input')
# The root correction ledgers contain unused A5/A7/A8 branches. Archive only rows used by S9.
for rel,rule in [('reviews/2026-10-08-plan-a-redesign/corrected-tests.tsv',lambda r:r['family'].startswith('A1-') or r['family']=='A2-cycle-alignment-redesign'),
                 ('reviews/2026-10-08-plan-a-analysis/corrected-tests.tsv',lambda r:r['family'] in ['A2-cycle-scores','A6-shape-endpoint'])]:
    reader=csv.DictReader(io.StringIO(files[rel].decode()),delimiter='\t'); selected=[r for r in reader if rule(r)]
    buf=io.StringIO(newline=''); w=csv.DictWriter(buf,fieldnames=reader.fieldnames,delimiter='\t');w.writeheader();w.writerows(selected)
    files[rel]=buf.getvalue().encode()
    receipt=next(r for r in rows if r['path']==rel);receipt.update(archive_sha256=hashlib.sha256(files[rel]).hexdigest(),bytes=len(files[rel]),filter='Only S9-used families retained; original full-file SHA is separate. '+str(sorted({r['family'] for r in selected})))
for stem in ['prepare_s9.py','build_s9.mjs','append_s9.py','verify_s9.py','verify_archive.py','archive_s9.py','skim_reference.py','render_final.mjs','README.zh.md','baseline-inspection.json','workbook-verification.json','final-format-inspection.json','table-sources-S9.tsv','s9-blocks.json','data/opposite-sign.primary-corrected.json']:
    add(str((HERE/stem).relative_to(ROOT)),'S9 export/correction code or manifest')
readme='''S9 companion source and code archive

Scope: frozen inputs/outputs for S9a A1, S9b raw cycle scores and support/balance diagnostics,
S9c A3, S9d A6, S9e netshift and opposite-sign decomposition. No statistics were rerun.
Paths are relative to the reproAPA project root. Original scripts/protocols are included
for scientific provenance; multi-analysis protocols/scripts can mention unused branches,
but A5/A7/A8 result products and route-mean result tables are excluded.

The two root corrected-tests.tsv ledgers are filtered copies, not byte-identical originals:
only S9-used families remain. FILES-SHA256.tsv records both original full-file SHA and
archived subset SHA plus exact filters. All other frozen source files preserve their bytes.
Workbook Sources records original source hashes; table-sources-S9.tsv records row and
column projections. The original S1-S8 archive and workbook are separate and unchanged.

The old opposite-sign.json and derive_opposite.py are archived as provenance and contain
a p50 setting-selection defect (substring matched supported50_p08_mc5). Do not use their
p50 dDWUI values. The new data/opposite-sign.primary-corrected.json and p50-primary.tsv
use exact setting=supported50_p08, contrast=P_M, design=same_library, CIP/WOT.
prepare_s9.py performs this copy-only correction and exports the three flattened blocks.
The polyApipe native-axis decomposition is preserved from the frozen JSON.

Original cycle caches/expression inputs support reproduction downstream of the original
data extraction. Remote raw sequencing matrices are not redistributed here; extraction
scripts/manifests preserve their provenance. Running original scripts may require the
repository's original runtime and paths; building the workbook uses the bundled artifact
runtime and the separate preserved S1-S8 workbook. Do not execute legacy scripts over
the frozen outputs without choosing an isolated run directory.
'''
rows.sort(key=lambda r:r['path']); buf=io.StringIO(newline='');w=csv.DictWriter(buf,fieldnames=['path','role','original_sha256','archive_sha256','filter','bytes'],delimiter='\t');w.writeheader();w.writerows(rows)
files['FILES-SHA256.tsv']=buf.getvalue().encode();files['README.txt']=readme.encode()
out=HERE/'Supplementary_Source_Code_S9.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for name,data in sorted(files.items()):
        info=zipfile.ZipInfo(name,(2026,10,8,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;z.writestr(info,data)
receipt={'archive':str(out.relative_to(ROOT)),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'files':len(files),'bytes':out.stat().st_size,'uncompressed_bytes':sum(map(len,files.values())),'filtered_root_ledgers':2,'all_other_source_bytes_preserved':True}
(HERE/'source-archive-verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
(HERE/'FILES-SHA256.tsv').write_text(buf.getvalue())
print(json.dumps(receipt,indent=2))
