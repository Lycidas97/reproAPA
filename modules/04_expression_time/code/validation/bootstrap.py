from pathlib import Path
import json
import numpy as np
import pandas as pd
from common import *

RUN.mkdir(exist_ok=True)
(RUN / 'shared').mkdir(exist_ok=True)
assert not (RUN / 'shared/ready.json').exists(), 'Already initialized; do not replace a running baseline.'
old = json.loads((ROOT / 'P00/baseline_input_manifest.json').read_text())['files']
for r in old:
    if r['path'].startswith(str(BASE)):
        assert sha(r['path']) == r['sha256'], r['path']
x = pd.read_csv(BASE / 'raw_cell_nodes.tsv.gz', sep='\t', low_memory=False)
x['state_key'] = np.where(x.species_key.eq('human'), x.original_state_id,
                         x.dataset.fillna('') + '|' + x.local_state_id.fillna(''))
x['sample_key'] = np.where(x.species_key.eq('human'), x.sample_accession, x.library_key)
x['raw_library_id'] = x.cell_id.str.split('|').str[1]
x['barcode'] = x.cell_id.str.rsplit('|', n=1).str[-1]
labels = pd.read_csv(ROOT / 'inputs/S5-local-state-annotations.tsv', sep='\t').rename(
    columns={'original_state_id':'state_key', 'study_id':'true_study'})
x = x.merge(labels[['state_key','true_study','display_group','marker_evidence',
                    'source_annotation','source_annotation_status','source_caveat']],
            on='state_key', how='left', validate='many_to_one', sort=False)
assert len(x) == 551201 and x.cell_id.is_unique and x.display_group.notna().all()
assert not x.duplicated(['sample_key','barcode']).any()
x['atlas_index'] = np.arange(len(x))
x.to_parquet(RUN / 'shared/atlas.parquet', index=False)
manifests = []
for name in ['GSE118258_production','GSE150311_production','GSE178325_regression',
             'GSE221739_production','GSE103221_production','GSE171127_production','GSE115943_production']:
    f = QC / 'manifests' / (name + '.tsv')
    d = pd.read_csv(f, sep='\t'); d['manifest_source'] = str(f)
    manifests.append(d)
m = pd.concat(manifests, ignore_index=True)
assert m.library_id.is_unique and m.sample_accession.is_unique and len(m) == 210
m = m[m.library_id.isin(x.raw_library_id)].copy()
assert len(m) == 161 and set(m.library_id) == set(x.raw_library_id)
assert m.input_adapter.eq('cellranger_filtered_h5').all()
assert all(Path(p).is_file() for p in m.input_path)
m.to_csv(RUN / 'shared/raw_count_manifest.tsv', sep='\t', index=False)
receipt([BASE/'raw_cell_nodes.tsv.gz', ROOT/'inputs/S5-local-state-annotations.tsv',
         *[QC/'manifests'/(n+'.tsv') for n in ['GSE118258_production','GSE150311_production',
         'GSE178325_regression','GSE221739_production','GSE103221_production','GSE171127_production','GSE115943_production']]],
        RUN/'shared/input_receipt.json')
write_json(RUN/'shared/ready.json', dict(cells=len(x),libraries=len(m),all_raw_paths_exist=True))
print('Frozen execution membership and 161 raw count paths verified.', flush=True)
