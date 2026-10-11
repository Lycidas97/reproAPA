"""Independent raw-workbook, SciPy-null, and raw-metacell QKI audit."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import sys
import zipfile
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact, false_discovery_control

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE / 'source/xlrd-2.0.2-py2.py3-none-any.whl'))
import xlrd

checks, max_error = 0, 0.0

def check(condition):
    global checks
    assert condition
    checks += 1

def equal(left, right, tolerance=1e-9):
    global max_error
    if pd.isna(left) or pd.isna(right):
        check(pd.isna(left) and pd.isna(right))
    else:
        error = abs(float(left) - float(right))
        max_error = max(max_error, error)
        check(error < tolerance)

def workbook(name):
    sheet = xlrd.open_workbook(HERE / 'source' / name).sheet_by_index(0)
    header = sheet.row_values(1)
    records = []
    for i in range(2, sheet.nrows):
        records.append(dict(zip(header, sheet.row_values(i))))
    return records

def main():
    freeze = json.loads((HERE / 'protocol_freeze.json').read_text())
    for relative, digest in freeze['files'].items():
        check(hashlib.sha256((HERE / relative).read_bytes()).hexdigest() == digest)
    with zipfile.ZipFile(HERE / 'source/supplements.zip') as archive:
        check(archive.testzip() is None)
        for name in ['KRNB_A_2294222_SM0800.xls', 'KRNB_A_2294222_SM0790.xls']:
            check(archive.read(name) == (HERE / 'source' / name).read_bytes())
    roster = pd.read_csv(ROOT / 'notes/2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz', sep='\t')
    roster = roster[roster.contrast.eq('P_M')]
    features = pd.read_csv(ROOT / 'notes/2026-10-10-recurrence-utr-ortholog/ortholog_assay_features.tsv.gz', sep='\t')
    features = features[features.contrast.eq('P_M') & np.isclose(features.threshold, .02)]
    covariates = { (r.species, r.gene_id): r._asdict() for r in features.itertuples(index=False) }
    events = pd.read_csv(HERE.parent / 'event_inventory.tsv.gz', sep='\t')
    strict = {}
    for key, part in events.groupby(['species', 'gene_id']):
        strict[key] = (all(part.driver_class.eq('strict_tandem_UTR')), all(part.full_primary_axis_tandem))
    output = pd.read_csv(HERE / 'ortholog_external_roster.tsv.gz', sep='\t')
    tests = pd.read_csv(HERE / 'overlap_tests.tsv', sep='\t')
    check(len(output) == 1502 and len(tests) == 36)
    check(len(roster) == 751 and roster.foreground_002.sum() == 27)
    external = {}
    for label, filename in [('QKI_KD', 'KRNB_A_2294222_SM0790.xls'), ('EMT', 'KRNB_A_2294222_SM0800.xls')]:
        raw = workbook(filename)
        check(len(raw) == len(set(r['name'] for r in raw)))
        check(all(r['name'].startswith('ENSG') for r in raw))
        external[label] = {r['name']: r for r in raw if math.isfinite(r['effect']) and math.isfinite(r['fdr_zero'])}
        saved = output[output.comparison.eq(label)].set_index('human_gene_id')
        check(int(saved.external_evaluable.sum()) == 605)
        check(int((saved.foreground_002 & saved.external_evaluable).sum()) == 25)
        for member in roster.itertuples():
            found = external[label].get(member.human_gene_id)
            record = saved.loc[member.human_gene_id]
            check(bool(record.external_evaluable) == (found is not None))
            if found is None:
                check(pd.isna(record.effect) and pd.isna(record.fdr_zero) and bool(record.missing_reason))
            else:
                equal(record.effect, found['effect'])
                equal(record.fdr_zero, found['fdr_zero'])
            driver = strict.get(('human', member.human_gene_id), (False, False))[0] and strict.get(('mouse', member.mouse_gene_id), (False, False))[0]
            axis = strict.get(('human', member.human_gene_id), (False, False))[1] and strict.get(('mouse', member.mouse_gene_id), (False, False))[1]
            check(record.driver_strict == driver and record.axis_strict == axis)
        expected = {'QKI_KD': (94, 25), 'EMT': (27, 52)}[label]
        counts = tuple(sum(r['fdr_zero'] < .05 and (r['effect'] < 0 if negative else r['effect'] > 0) for r in raw) for negative in [True, False])
        check(counts == expected)
    coarse = ['measured_pattern', 'proximal_capacity', 'distal_capacity', 'PAS_number_bin']
    full = coarse + ['mean_GEX_bin', 'APA_coverage_bin', 'baseline_bin']
    recalculated_p = []
    for test in tests.itertuples():
        groups, observed, n_core, K_total, N_total, variable = {}, 0, 0, 0, 0, 0
        for member in roster.itertuples():
            found = external[test.comparison].get(member.human_gene_id)
            if found is None:
                continue
            driver = strict[('human', member.human_gene_id)][0] and strict[('mouse', member.mouse_gene_id)][0]
            axis = strict[('human', member.human_gene_id)][1] and strict[('mouse', member.mouse_gene_id)][1]
            if (test.panel == 'driver_strict' and not driver) or (test.panel == 'axis_strict' and not axis):
                continue
            positive = found['effect'] < 0 if test.direction == 'proximal' else found['effect'] > 0
            response = found['fdr_zero'] < .05 and positive
            foreground = bool(member.foreground_002)
            key = () if test.model == 'unmatched' else tuple(covariates[(s, getattr(member, s+'_gene_id'))][c]
                for s in ['human', 'mouse'] for c in (coarse if test.model == 'COARSE' else full))
            group = groups.setdefault(key, [0, 0, 0])
            group[0] += 1; group[1] += response; group[2] += foreground
            N_total += 1; K_total += response; n_core += foreground; observed += foreground and response
        pmf = np.array([1.0])
        for N, K, n in groups.values():
            possible = np.arange(n+1)
            component = hypergeom.pmf(possible, N, K, n)
            pmf = np.convolve(pmf, component)
            if min(n, K) > max(0, n+K-N):
                variable += n
        x = np.arange(len(pmf)); expected = float(x @ pmf); variance = float((x-expected)**2 @ pmf)
        p = min(1.0, float(pmf[observed:].sum())) if variance > 1e-14 else np.nan
        check((N_total, n_core, K_total, observed, variable) == (test.N, test.n_core, test.n_external_response, test.observed, test.variable_foreground))
        equal(expected, test.expected); equal(variance, test.null_variance); equal(p, test.p_value)
        check((test.status == 'zero_null_variance') == (variance <= 1e-14))
        if test.model == 'unmatched' and variance > 1e-14:
            table = [[observed, n_core-observed], [K_total-observed, N_total-n_core-K_total+observed]]
            equal(fisher_exact(table, alternative='greater').pvalue, p)
        recalculated_p.append(p if np.isfinite(p) else 1.0)
    q = false_discovery_control(recalculated_p, method='bh')
    for left, right in zip(q, tests.q_BH36):
        equal(left, right)
    old = ROOT / 'notes/2026-10-02-rbp-expression'
    archive = np.load(old / 'expression_inputs.npz', allow_pickle=False)
    produced = pd.read_csv(HERE / 'qki_expression_pairs.tsv', sep='\t').set_index(['study_id', 'library_key'])
    old_pairs = pd.read_csv(old / 'paired_library_effects.tsv', sep='\t')
    for study in ['GSE118258', 'GSE150311', 'GSE103221', 'GSE115943']:
        metadata = json.loads(str(archive[study + '_meta_json']))
        names = list(archive[study + '_factors']); factor = 'QKI' if study in ['GSE118258', 'GSE150311'] else 'Qki'
        index = names.index(factor); sums = {}; cells = {}; times = {}
        for i, item in enumerate(metadata):
            if item['node'] not in ['P', 'M']:
                continue
            key = item['library_key'], item['node']
            sums.setdefault(key, []).append(math.log2(1 + float(archive[study+'_counts'][i,index])/float(archive[study+'_total_gex'][i])*1e6))
            cells[key] = cells.get(key, 0) + item['n_cells']; times[key] = str(item['experimental_time'])
        pairs = []
        for library in sorted(set(k[0] for k in sums)):
            keys = [(library, 'P'), (library, 'M')]
            if not all(k in sums and len(sums[k]) >= 3 and cells[k] >= 30 for k in keys):
                continue
            delta = sum(sums[keys[0]])/len(sums[keys[0]]) - sum(sums[keys[1]])/len(sums[keys[1]])
            row = produced.loc[(study, library)]; equal(delta, row.delta_log2CPM1)
            baseline = old_pairs[old_pairs.study_id.eq(study) & old_pairs.library_key.eq(library) & old_pairs.factor.eq(factor)]
            check(len(baseline) == 1); equal(delta, baseline.delta_log2CPM1.iloc[0])
            pairs.append((times[keys[0]], delta))
        summary = pd.read_csv(HERE / 'qki_expression_summary.tsv', sep='\t').set_index('study_id').loc[study]
        equal(sum(v for _, v in pairs)/len(pairs), summary.delta_log2CPM1)
        dates = {}
        for time, value in pairs:
            dates.setdefault(time, []).append(value)
        equal(sum(sum(v)/len(v) for v in dates.values())/len(dates), summary.time_equal_delta_log2CPM1)
    report = dict(status='PASS', checks=checks, skipped=0, max_absolute_numeric_error=max_error,
                  source_zip_CRC_pass=True, raw_author_response_counts_reproduced=True,
                  independent_scipy_null_and_fisher_pass=True, raw_qki_expression_and_old_tables_pass=True,
                  timestamp_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (HERE / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
