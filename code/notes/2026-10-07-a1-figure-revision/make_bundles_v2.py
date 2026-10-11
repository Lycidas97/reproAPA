"""v2 (2026-10-07, later): also replaces cross_controls/cross_species_power.tsv with the eligibility-masked ladder
rerun (notes/2026-10-07-cross-null-followup). Derive figure-input bundles with the A1-corrected cross-overlap tables (2026-10-07).

Every file is hard-linked from the sealed source bundle except the two replaced tables and the interface JSONs whose
SHA fields change. Replaced tables keep the original schema:
  overlap_threshold_layer_sensitivity.tsv  <- notes/2026-10-07-cross-overlap-null-a1/polyapipe-grid-results.tsv (verbatim)
  overlap_null_results.tsv                 <- joint rows: the grid row (supported50, 0.02) plus q over two contrasts;
                                              Quintile/Range rows: analytic/null mean, 95% range, p and q from
                                              polyapipe-null-results.tsv. Fields not in the receipt are left empty
                                              (analytic_null_sd, null_median, strata, singleton_strata, movable_null_sd,
                                              seed, movable fractions); no figure reads them.
"""
import csv, hashlib, io, json, os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
A1 = ROOT/'notes/2026-10-07-cross-overlap-null-a1'
HERE = Path(__file__).resolve().parent
SOURCES = {'supplement-v2': ROOT/'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle'}
POWER = 'tables/cross_controls/cross_species_power.tsv'
FOLLOW = ROOT/'notes/2026-10-07-cross-null-followup/ladder-results/polyapipe/ladder/cross_species_power.tsv'
NULL, GRID = 'tables/cross_overlap/overlap_null_results.tsv', 'tables/cross_overlap/overlap_threshold_layer_sensitivity.tsv'
INTERFACE = ['table-manifest.json', 'result-records.json', 'result-records.tsv', 'result.json', 'bundle-complete.json']

def sha(b): return hashlib.sha256(b).hexdigest()
def rows(p): return list(csv.DictReader(open(p), delimiter='\t'))
def tsv(cols, rs):
    s = io.StringIO(); w = csv.DictWriter(s, cols, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rs)
    return s.getvalue().encode()

grid_bytes = (A1/'polyapipe-grid-results.tsv').read_bytes()
grid = rows(A1/'polyapipe-grid-results.tsv'); fixed = {r['contrast']+r['model']: r for r in rows(A1/'polyapipe-null-results.tsv')}

def corrected_null(old_path):
    old = rows(old_path); cols = list(old[0].keys()); out = []
    for o in old:
        f = fixed[o['contrast']+o['model']]; n = {c: '' for c in cols}
        for c in ['contrast', 'model', 'observed', 'orthologs', 'permutations']: n[c] = o[c]
        assert int(f['observed']) == int(o['observed'])
        if o['model'] == 'joint_measurement_PAS':
            g = [r for r in grid if r['contrast'] == o['contrast'] and r['layer'] == 'supported50' and abs(float(r['threshold'])-.02) < 1e-9]
            assert len(g) == 1; g = g[0]
            for c in cols:
                if c in g and c not in ('q_16',): n[c] = g[c]
            assert abs(float(g['empirical_p'])-float(f['p'])) < 1e-12 and abs(float(g['analytic_null_mean'])-float(f['analytic_expected'])) < 1e-9
        else:
            n.update(analytic_null_mean=f['analytic_expected'], null_mean=f['null_mean'], null_q025=f['null_q025'],
                     null_q975=f['null_q975'], empirical_p=f['p'])
        n['q_within_model_two_contrasts'] = f['q_two_contrasts']
        out.append(n)
    return tsv(cols, out)

report = {}
for label, src in SOURCES.items():
    dst = HERE/f'bundle-{label}'
    assert not dst.exists(), dst
    new = {NULL: corrected_null(src/NULL), GRID: grid_bytes, POWER: FOLLOW.read_bytes()}
    for rel in new: assert open(src/rel).readline() == new[rel].decode().splitlines()[0]+'\n', rel
    swap = {sha((src/rel).read_bytes()): sha(b) for rel, b in new.items()}
    sizes = {rel: len(b) for rel, b in new.items()}
    for folder, dirs, files in os.walk(src):
        rel_dir = Path(folder).relative_to(src)
        for name in files:
            rel = str(rel_dir/name) if str(rel_dir) != '.' else name
            target = dst/rel; target.parent.mkdir(parents=True, exist_ok=True)
            if rel in new: target.write_bytes(new[rel])
            elif rel in INTERFACE: continue
            else: os.link(src/rel, target)
    for name in INTERFACE[:-1]:
        if not (src/name).exists(): continue
        text = (src/name).read_text()
        for a, b in swap.items(): text = text.replace(a, b)
        if name == 'table-manifest.json':
            d = json.loads(text)
            for t in d['tables']:
                if t['path'] in sizes: t['bytes'] = sizes[t['path']]
            text = json.dumps(d, indent=1, ensure_ascii=False)+'\n'
        (dst/name).write_text(text)
    marker = json.loads((src/'bundle-complete.json').read_text())
    marker['table_manifest_sha256'] = sha((dst/'table-manifest.json').read_bytes())
    rec = json.loads((dst/'result-records.json').read_text())
    rec['table_manifest_sha256'] = marker['table_manifest_sha256']
    (dst/'result-records.json').write_text(json.dumps(rec, indent=1, ensure_ascii=False)+'\n')
    marker['result_records_sha256'] = sha((dst/'result-records.json').read_bytes())
    marker['derived_from'] = dict(bundle=str(src.relative_to(ROOT)), reason='A1 eligibility-corrected cross-overlap nulls (2026-10-07)',
                                  replaced={rel: dict(old=o, new=swap[o]) for rel in new for o in swap if sha((src/rel).read_bytes()) == o})
    (dst/'bundle-complete.json').write_text(json.dumps(marker, indent=2)+'\n')
    assert rec['plot_contract_sha256'] == marker['plot_contract_sha256']
    report[label] = marker['derived_from']
    sys.path.insert(0, str(ROOT/'analysis_code/pas_rerun'))
    from result_records import ResultBundle
    ResultBundle(dst)  # interface loads
(HERE/'bundles-receipt-v2.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=1))
