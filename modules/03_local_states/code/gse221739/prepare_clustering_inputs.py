#!/usr/bin/env python3
"""Add GSE221739 to the frozen clustering workflow without changing its old release."""

import copy
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
import yaml

ROOT = Path('/root/nfsdata/apadetect/analysis_projects')
QC = ROOT / 'reprogramming_unified_rna_qc_v3.2_release'
CLUSTER = ROOT / 'reprogramming_qc_v3_2_clustering_release_v1'
RUN = QC / 'runs/gse221739_production'
TARGET = CLUSTER / 'gse221739_generated'
sys.path.insert(0, str(CLUSTER / 'scripts'))
from build_qc_clustering_inputs import ordered_sha256, set_sha256, tier_contract


def main():
    acceptance = json.loads((RUN / 'acceptance.json').read_text())
    assert acceptance['status'] == 'PASS'
    assert acceptance['libraries'] == 17 and acceptance['called_cells'] == 65686
    assert json.loads((RUN / 'manual_visual_review.json').read_text())['status'] == 'PASS'
    tiers, memberships = tier_contract(RUN)
    all_cells = pd.concat(memberships.values(), ignore_index=True)
    assert len(all_cells) == 65686 and all_cells.global_cell_id.is_unique
    assert all_cells.study_id.eq('GSE221739').all()
    manifest_path = QC / 'manifests/GSE221739_production.tsv'
    manifest = pd.read_csv(manifest_path, sep='\t', dtype=str)
    assert len(manifest) == 17 and manifest.library_id.is_unique
    assert set(manifest.library_id) == set(all_cells.library_id)
    assert manifest.input_adapter.eq('cellranger_filtered_h5').all()
    assert all(Path(p).is_file() for p in manifest.input_path)
    for name in ('primary_pass', 'sensitivity_only'):
        assert set(memberships[name].library_id) == set(manifest.library_id), name

    config = copy.deepcopy(yaml.safe_load((CLUSTER / 'production_generated/configs/gse118258.yaml').read_text()))
    config['study_id'] = 'GSE221739'
    config['species'] = sorted(manifest.species.unique().tolist())
    config['membership_path'] = str(RUN / 'primary_membership.tsv.gz')
    config['sensitivity_membership_path'] = str(RUN / 'sensitivity_membership.tsv.gz')
    config['technical_membership_path'] = str(RUN / 'technical_membership.tsv.gz')
    config['library_metadata_path'] = str(manifest_path)
    preflight = json.loads((QC / 'production/gse221739_preflight.json').read_text())
    config['qc_contract'] = {
        'method_contract_sha256': preflight['method_contract_sha256'],
        'production_acceptance_path': str(RUN / 'acceptance.json'),
        'independent_validation_path': str(RUN / 'independent_validation.tsv'),
        'tiers': tiers,
    }
    primary = memberships['primary_pass']
    sensitivity = memberships['sensitivity_only']
    scope = config['scopes']['studywide']
    scope.update({
        'n_cells': len(primary),
        'source_membership_order_sha256': ordered_sha256(primary.global_cell_id),
        'cell_set_sha256': set_sha256(primary.global_cell_id),
        'sensitivity_n_cells': len(sensitivity),
        'sensitivity_source_membership_order_sha256': ordered_sha256(sensitivity.global_cell_id),
        'sensitivity_cell_set_sha256': set_sha256(sensitivity.global_cell_id),
        'expression_sources': [
            {'selector': 'library_id', 'value': r.library_id, 'mode': '10x_h5', 'path': r.input_path}
            for r in manifest.sort_values('library_id', kind='mergesort').itertuples(index=False)
        ],
    })
    project = copy.deepcopy(yaml.safe_load((CLUSTER / 'production_generated/config.production.yaml').read_text()))
    project['project_config_path'] = str(TARGET / 'config.production.yaml')
    project['output_root'] = 'results_qc_v3_2_gse221739_v1'
    project['dataset_configs'] = {'gse221739': 'gse221739_generated/configs/gse221739.yaml'}
    assert not (TARGET / 'config.production.yaml').exists(), 'Refusing to overwrite existing config'
    output = CLUSTER / project['output_root']
    assert not output.exists() or not any(output.iterdir()), 'Refusing to reuse nonempty output'
    (TARGET / 'configs').mkdir(parents=True, exist_ok=True)
    config_path = TARGET / 'configs/gse221739.yaml'
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    project_path = TARGET / 'config.production.yaml'
    project_path.write_text(yaml.safe_dump(project, sort_keys=False))
    report = {
        'status': 'PASS', 'study_id': 'GSE221739', 'libraries': 17,
        'tiers': tiers, 'same_method_as': 'production_generated/configs/gse118258.yaml',
        'n_resamples': project['selection']['n_resamples'],
        'n_resolutions': len(project['selection']['resolutions']),
        'old_release_modified': False,
        'config_sha256': hashlib.sha256(config_path.read_bytes()).hexdigest(),
        'project_config_sha256': hashlib.sha256(project_path.read_bytes()).hexdigest(),
    }
    (TARGET / 'input_acceptance.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
