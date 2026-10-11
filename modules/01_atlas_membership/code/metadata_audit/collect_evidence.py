"""Read existing alignment provenance on 213docker; never rerun an analysis."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
scope = json.loads((ROOT / "sources/alignment-scope.json").read_text())
remote = r'''
import csv, hashlib, json, re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

root = Path('/root/nfsdata/apadetect/data_collection/reprogramming_differentiation_apa')
jobs = list(csv.DictReader((root/'tables/main_raw_quantification_jobs.tsv').open(), delimiter='\t'))
libraries = list(csv.DictReader((root/'tables/main_raw_to_expression_bam_libraries.tsv').open(), delimiter='\t'))
job_by_id = {r['library_id']: r for r in jobs}
library_by_id = {r['library_id']: r for r in libraries}
records, sources = [], {}

def read_source(path):
    text = path.read_text()
    sources[str(path)] = {'sha256': hashlib.sha256(text.encode()).hexdigest(), 'text': text}
    return text

for row in SCOPE:
    bam = Path(row['bam_path']).resolve(strict=True)
    run = bam.parent.parent
    versions = json.loads(read_source(run/'_versions'))
    invocation = read_source(run/'_invocation')
    cmdline = read_source(run/'_cmdline').strip()
    rec = dict(row, run_dir=str(run), resolved_bam=str(bam), cellranger_version=versions['pipelines'])
    for key in ['reference_path', 'chemistry', 'r1_length', 'r2_length', 'include_introns', 'force_cells', 'recovered_cells', 'no_bam', 'no_secondary_analysis']:
        found = re.search(r'^\s*'+key+r'\s*=\s*(.*?),\s*$', invocation, re.M)
        rec[key] = json.loads(found.group(1)) if found else 'UNRESOLVED'
    meta = library_by_id.get(row['library_id'], {})
    rec['run_accessions'] = meta.get('run_list', '')
    if not rec['run_accessions']:
        rec['run_accessions'] = ';'.join(sorted(set(re.findall(r'SRR\d+', cmdline))))
    rec['sample_name'] = meta.get('sample_name', run.name)
    rec['table_chemistry'] = job_by_id.get(row['library_id'], {}).get('chemistry', '')
    rec['chemistry_table_disagrees'] = bool(rec['table_chemistry'] and rec['table_chemistry'] != rec['chemistry'])
    metrics = list(csv.DictReader(read_source(run/'outs/metrics_summary.csv').splitlines()))[0]
    for key, column in [('called_cells','Estimated Number of Cells'), ('reads','Number of Reads'), ('median_genes','Median Genes per Cell'), ('median_umis','Median UMI Counts per Cell')]:
        rec[key] = int(metrics[column].replace(',',''))
    for key in ['raw_feature_bc_matrix.h5', 'filtered_feature_bc_matrix.h5', 'possorted_genome_bam.bam.bai', 'web_summary.html']:
        f = run/'outs'/key
        rec[key+'_present'] = f.exists() and f.stat().st_size > 0
    read_path = re.search(r'"read_path":\s*"([^"]+)"', invocation)
    if read_path:
        q = Path(read_path.group(1))
        rec['fastq_input_dir'] = str(q)
        rec['fastq_source_paths'] = sorted({str(f.resolve()) for f in q.glob('*.fastq*')})
    records.append(rec)

references = {}
for path in sorted({r['reference_path'] for r in records}):
    references[path] = json.loads(read_source(Path(path)/'reference.json'))
for rel in ['workflow/Snakefile','workflow/config/raw_download.yaml','tables/reference_genome_registry.tsv','tables/sra_fastq_route_plan.tsv','tables/dataset_audit.tsv']:
    read_source(root/rel)
print(json.dumps({'checked_at_utc':datetime.now(timezone.utc).isoformat(), 'audit_scope':'existing log and output metadata; no BAM body scan or biological reanalysis', 'records':records, 'references':references, 'sources':sources},ensure_ascii=False))
'''
result = subprocess.run(
    ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "213docker", "python3 -"],
    input="SCOPE = " + repr(scope) + "\n" + remote, text=True, capture_output=True, check=True,
)
data = json.loads(result.stdout)
(ROOT / "sources/remote-alignment-evidence.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
fields = [k for k in data['records'][0] if k != 'fastq_source_paths']
with (ROOT / 'alignment-library-audit.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(data['records'])
summary = {}
for row in data['records']:
    study = row['dataset_id'].split('_')[0].upper()
    if study not in summary:
        summary[study] = dict(libraries=0, called_cells=0, chemistry={}, versions={}, include_introns={}, table_disagreements=[])
    s = summary[study]
    s['libraries'] += 1
    s['called_cells'] += row['called_cells']
    for key, target in [('chemistry','chemistry'), ('cellranger_version','versions'), ('include_introns','include_introns')]:
        value = str(row[key])
        s[target][value] = s[target].get(value,0) + 1
    if row['chemistry_table_disagrees']:
        s['table_disagreements'].append(row['sample_id'])
(ROOT / 'alignment-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
