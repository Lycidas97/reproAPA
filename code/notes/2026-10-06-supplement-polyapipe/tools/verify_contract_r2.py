"""Verify contract r2: sealed r1 bundle untouched, r2 tables byte-identical to r1 (or to the pinned origin for added
tables), the r2 interface loads, and every candidate receipt read only bound r2 tables."""
import json, glob, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[2]/'analysis_code/pas_rerun'))
import supdraw as S
from supdraw import digest
from result_records import ResultBundle
R1, R2 = S.BUNDLE_R1, S.BUNDLE
out = dict()
m1 = json.loads((R1/'bundle-complete.json').read_text())
out['r1_seal_intact'] = all(digest(R1/n) == m1[k] for n, k in [('plot-contract.json', 'plot_contract_sha256'),
                            ('table-manifest.json', 'table_manifest_sha256'), ('result-records.json', 'result_records_sha256')])
b2 = ResultBundle(R2); out['r2_interface_loads'] = True; out['r2_components'] = len(b2.components); out['r2_panels'] = len(b2.allowed)
t1 = {r['path']: r for r in json.loads((R1/'table-manifest.json').read_text())['tables']}
t2 = {r['path']: r for r in json.loads((R2/'table-manifest.json').read_text())['tables']}
out['r1_tables'] = len(t1); out['r2_tables'] = len(t2)
out['r1_tables_identical_in_r2'] = all(p in t2 and t2[p]['sha256'] == r['sha256'] and digest(R2/p) == r['sha256'] for p, r in t1.items())
added = {p: r for p, r in t2.items() if p not in t1}
out['added_tables'] = {p: dict(sha256=r['sha256'], origin=r['origin'], origin_sha_ok=digest(S.ROOT/r['origin']) == r['sha256']) for p, r in added.items()}
reads = []
for f in sorted(glob.glob(str(HERE.parent/'candidates/*/*.receipt.json'))):
    r = json.loads(Path(f).read_text())
    for t, v in r['tables'].items():
        for pnl in v['panels']:
            reads.append(dict(figure=r['figure'], panel=pnl, table=t, mode=v['mode'],
                              bound=t in b2.allowed.get(pnl, set()), sha_ok=t2[t]['sha256'] == v['sha256']))
out['receipt_reads'] = len(reads)
out['unbound_or_changed_reads'] = [x for x in reads if not (x['bound'] and x['sha_ok'])]
out['receipt_contract_revisions'] = sorted({json.loads(Path(f).read_text()).get('contract_revision') for f in glob.glob(str(HERE.parent/'candidates/*/*.receipt.json'))})
ok = out['r1_seal_intact'] and out['r1_tables_identical_in_r2'] and all(v['origin_sha_ok'] for v in out['added_tables'].values()) and not out['unbound_or_changed_reads']
out['status'] = 'pass' if ok else 'FAIL'
(HERE.parent/'verification'/'contract-r2.json').write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(json.dumps({k: v for k, v in out.items() if k not in ('added_tables',)}, indent=1)[:1500])
