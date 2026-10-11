#!/usr/bin/env python3
"""Prepare the added study with the installed QC release's existing builders."""
import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import h5py
import yaml
from snakemake.utils import validate

ROOT = Path('/root/nfsdata/apadetect/analysis_projects/reprogramming_unified_rna_qc_v3.2_release')
COLLECTION = Path('/root/nfsdata/apadetect/data_collection/reprogramming_differentiation_apa')
DATASET = 'gse221739_human_pluripotency_reprogramming'
NAME = 'GSE221739_production'
sys.path.insert(0, str(ROOT / 'production'))
from build_inputs import write_study, write_config
sys.path.insert(0, str(ROOT / 'scripts'))
from common import load_config, read_manifest

output = ROOT / 'runs/gse221739_production'
if output.exists() and any(output.iterdir()):
    raise RuntimeError('Refusing to replace an existing added-study QC run')
config_path = ROOT / 'configs' / (NAME + '.yaml')
if config_path.exists():
    raise RuntimeError('Added-study config already exists; inspect before replacing')

source = COLLECTION / 'metadata_raw/geo/GSE221739_family.soft.gz'
metadata = {}
sample = None
with gzip.open(source, 'rt') as handle:
    for line in handle:
        line = line.rstrip('\n')
        if line.startswith('^SAMPLE = '):
            sample = line.split(' = ', 1)[1]
            metadata[sample] = {}
        elif sample and line.startswith('!Sample_title = '):
            metadata[sample]['title'] = line.split(' = ', 1)[1]
        elif sample and line.startswith('!Sample_characteristics_ch1 = '):
            key, value = line.split(' = ', 1)[1].split(': ', 1)
            metadata[sample][key] = value

library_source = COLLECTION / 'tables/main_raw_to_expression_bam_libraries.tsv'
with library_source.open() as handle:
    libraries = [r for r in csv.DictReader(handle, delimiter='\t') if r['dataset_id'] == DATASET]
rows = []
for item in libraries:
    sample = item['sample_id']
    meta = metadata[sample]
    match = re.fullmatch(r'D(\d+), pool ([AB])(?:, replicate (\d+))?', meta['title'])
    if match is None or meta['time'] != 'D' + match.group(1) or meta['pool'] != match.group(2):
        raise ValueError('GEO title/characteristic mismatch: ' + sample)
    if meta['treatment'] != 'mRNA transfections':
        raise ValueError('Unexpected treatment: ' + sample)
    if not Path(item['matrix_bam_done']).is_file():
        raise FileNotFoundError(item['matrix_bam_done'])
    matrix = Path(item['cellranger_run_dir']) / 'outs/filtered_feature_bc_matrix.h5'
    with h5py.File(matrix) as handle:
        shape = tuple(int(v) for v in handle['matrix/shape'][:])
        n = len(handle['matrix/barcodes'])
        if shape[1] != n:
            raise ValueError('H5 barcode/shape mismatch: ' + sample)
    with (Path(item['archive_qc_dir']) / 'metrics_summary.csv').open() as handle:
        expected = int(next(csv.DictReader(handle))['Estimated Number of Cells'].replace(',', ''))
    if n != expected:
        raise ValueError('Cell Ranger metrics/H5 count mismatch: ' + sample)
    rows.append(dict(
        study_id='GSE221739', species='Homo sapiens', dataset_id=DATASET,
        library_id=item['library_id'], sample_id=sample, sample_accession=sample,
        sample_title=meta['title'], experimental_time=meta['time'], condition='mRNA_reprogramming',
        author_state='NOT_USED', metadata_derivation='GEO sample title and characteristics: collection day, mRNA transfections; pool and capture retained in sample title',
        source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + sample,
        input_adapter='cellranger_filtered_h5', input_path=str(matrix), expected_called_cells=n,
        basic_qc_evidence_adapter='policy_document_recompute', prior_basic_qc_path='contract_inputs/historical_qc/report.md',
        diagnostic_evidence_adapter='scrublet_policy_recompute', prior_diagnostic_path='contract_inputs/historical_qc/production_policy_recommendation.tsv'))
assert len(rows) == 17 and sum(r['expected_called_cells'] for r in rows) == 65686
write_study(ROOT, NAME, rows, [])
write_config(ROOT, NAME, 'GSE221739', 'Homo sapiens', 'GSE221739 complete 17-library QC', rows,
             r'^D(?:3|5|7|9|11|13|15)[AB](?:_rep[12])?_GSM[0-9]+$', 'gse221739_production')
config = yaml.safe_load(config_path.read_text())
config['output_mode'] = 'resume'
config_path.write_text(yaml.safe_dump(config, sort_keys=False))
validate(config, str(ROOT / 'config.schema.yaml'))
baseline = load_config(ROOT / 'configs/GSE118258_production.yaml')
keys = ['adapters', 'basic_qc', 'sample_gates', 'scrublet', 'global_cell_id_template']
assert all(config[k] == baseline[k] for k in keys)
read_manifest(ROOT, config)
report = {'status': 'PASS', 'libraries': len(rows), 'called_cells': 65686,
          'same_method_as': 'GSE118258_production',
          'method_contract_sha256': hashlib.sha256(json.dumps({k:config[k] for k in keys}, sort_keys=True).encode()).hexdigest(),
          'source_sha256': {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,library_source,config_path]},
          'old_six_study_registry_modified': False}
preflight = ROOT / 'production/gse221739_preflight.json'
preflight.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report))
