#!/usr/bin/env python3
"""Freeze unique human cells, exact raw features, reference sampling and time semantics."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return pd.read_csv(path, sep='\t', dtype=str)


def decode(values):
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def metadata(row):
    study, time, title = row.study_id, row.experimental_time, row.sample_title.strip()
    result = dict(series=study, time_kind='unordered_reference', actual_day=np.nan,
                  stage='', stage_order=np.nan, within_stage_day=np.nan, donor='',
                  pool='', capture='', numeric_time_eligible=False)
    if study in ('GSE118258', 'GSE221739'):
        match = re.fullmatch(r'D(\d+(?:\.\d+)?)', time)
        if match is None:
            raise ValueError(f'Unexpected numerical time: {time}')
        result.update(time_kind='actual_day', actual_day=float(match[1]), numeric_time_eligible=True)
        if study == 'GSE221739':
            match = re.search(r'pool ([AB])', title)
            if match is None:
                raise ValueError(f'Missing GSE221739 pool: {title}')
            result['pool'] = match[1]
            capture = re.search(r'replicate (\d+)', title)
            result['capture'] = capture[1] if capture else 'unsplit'
    elif study == 'GSE150311':
        if time == 'D0_D7_pool':
            result['time_kind'] = 'pooled_days_unresolved'
        elif time == 'endpoint':
            result['time_kind'] = 'endpoint_culture_reference'
        else:
            raise ValueError(f'Unexpected GSE150311 time: {time}')
    elif study == 'GSE178325':
        match = re.search(r'-(0330|0605|0618|0809|1117|1230)$', title)
        result['donor'] = match[1] if match else ''
        result['series'] = study + ('::donor_' + result['donor'] if match else '::reference')
        stage = re.fullmatch(r'Stage_(IV|III|II|I)(?:_day_(\d+(?:\.\d+)?))?', time)
        if stage:
            result.update(time_kind='ordered_stage', stage=stage[1],
                          stage_order={'I': 1, 'II': 2, 'III': 3, 'IV': 4}[stage[1]])
            if stage[2] is not None:
                result['within_stage_day'] = float(stage[2])
        elif time in ('source', 'endpoint', 'pluripotent_reference'):
            result['time_kind'] = {'source': 'source_anchor', 'endpoint': 'endpoint_reference',
                                   'pluripotent_reference': 'external_reference'}[time]
        else:
            raise ValueError(f'Unexpected chemical time: {time}')
    else:
        raise ValueError(study)
    # Conditions are retained in sampling to keep perturbation/selection branches.
    result['reference_stratum'] = '|'.join([result['series'], time, row.condition, title if '0330' in title else ''])
    return result


def balanced_sample(frame, cap, seed):
    selected = []
    for stratum, group in frame.groupby('reference_stratum', sort=True):
        rng = np.random.default_rng(seed + int(hashlib.sha256(stratum.encode()).hexdigest()[:8], 16))
        queues = [rng.permutation(g.sort_values('cell_id').index.to_numpy())
                  for _, g in group.groupby('library_key', sort=True)]
        # Round-robin quotas retain small libraries before redistributing unused quota.
        choice = []
        for position in range(max(map(len, queues))):
            for queue in queues:
                if position < len(queue):
                    choice.append(queue[position])
                    if len(choice) == min(cap, len(group)):
                        break
            if len(choice) == min(cap, len(group)):
                break
        selected.extend(choice)
    return frame.loc[selected].sort_values('cell_id').reset_index(drop=True)


def temporal_plan(manifest, route):
    edges = []
    groups = [(study, manifest[manifest.study_id.eq(study)], 'actual_day')
              for study in ['GSE118258', 'GSE221739']]
    for name, rows in route.groupby('dataset', sort=True):
        selected = manifest[manifest.study_id.eq('GSE178325') & manifest.library_id.isin(rows.library_id)]
        groups.append((name, selected, 'ordered_stage_without_elapsed_days'))
    for name, rows, mode in groups:
        nodes = []
        for label, node in rows.groupby('experimental_time', sort=False):
            first = node.iloc[0]
            if mode == 'actual_day':
                order = (float(first.actual_day), 0)
            elif first.time_kind == 'source_anchor':
                order = (0, 0)
            elif first.time_kind == 'endpoint_reference':
                order = (5, 0)
            elif first.time_kind == 'ordered_stage':
                order = (int(first.stage_order), float(first.within_stage_day) if pd.notna(first.within_stage_day) else 0)
            else:
                raise ValueError(f'Unsupported temporal node: {name} {label}')
            nodes.append((order, label, node))
        nodes.sort(key=lambda x: x[0])
        for left, right in zip(nodes[:-1], nodes[1:]):
            edges.append({'route_id': name, 'time_from': left[1], 'time_to': right[1],
                          'relation': mode, 'elapsed_days': right[0][0] - left[0][0] if mode == 'actual_day' else np.nan,
                          'source_libraries': json.dumps(left[2].library_key.tolist()),
                          'target_libraries': json.dumps(right[2].library_key.tolist()),
                          'status': 'PLANNED_NOT_SOLVED', 'shared_source_counts_as_independent_replicate': False})
    return pd.DataFrame(edges)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    target = Path(config['output_root']) / 'prepared'
    target.mkdir(parents=True, exist_ok=True)
    if (target / 'preflight.json').exists():
        raise FileExistsError('Prepared inputs already exist; inspect before rebuilding')
    primary_parts, sensitivity_parts, technical_parts, manifests, provenance = [], [], [], [], []
    study_rows = []
    for study, spec in config['studies'].items():
        qc = Path(config['qc_root']) / 'runs' / spec['qc_run']
        scope = Path(config['clustering_root']) / spec['clustering_release']
        manifest_path = Path(config['qc_root']) / 'manifests' / spec['qc_manifest']
        manifest = read(manifest_path)
        assert len(manifest) == spec['expected_libraries'] and manifest.library_id.is_unique
        assert manifest.study_id.eq(study).all() and manifest.species.eq('Homo sapiens').all()
        validation = read(qc / 'independent_validation.tsv')
        assert len(validation) and validation.status.eq('PASS').all()
        tiers = {tier: read(qc / f'{tier}_membership.tsv.gz') for tier in ['primary', 'sensitivity', 'technical']}
        merged = pd.concat(tiers.values())
        assert merged.global_cell_id.is_unique
        assert set(merged.library_id) == set(manifest.library_id)
        assert len(merged) == pd.to_numeric(manifest.expected_called_cells).sum()
        for tier, cells in tiers.items():
            assert len(cells) == spec['expected_' + tier]
            assert cells.study_id.eq(study).all()
        partition = read(scope / 'frozen/cell_clusters.tsv.gz')
        projection = read(scope / 'sensitivity/projection/sensitivity_projection.tsv.gz')
        assert set(partition.cell_id) == set(tiers['primary'].global_cell_id)
        assert set(projection.cell_id) == set(tiers['sensitivity'].global_cell_id)
        assert partition.cell_id.is_unique and projection.cell_id.is_unique
        assert partition.cluster_id.notna().all() and partition.cluster_id.nunique() == spec['expected_states']
        transfer = read(scope / 'full_partition_audit/cluster_transfer_audit.tsv')
        programs = read(scope / 'program_audit/cluster_program_scores.tsv')
        labels = partition[['cell_id', 'cluster_id']].merge(
            transfer[['cluster_id', 'cluster_release_class', 'review_reasons']], on='cluster_id', validate='many_to_one'
        ).merge(programs[['cluster_id', 'requires_program_review', 'program_review_reasons']], on='cluster_id', validate='many_to_one')
        labels['original_state_id'] = study + '::studywide::' + labels.cluster_id
        primary = tiers['primary'].rename(columns={'global_cell_id': 'cell_id'}).merge(labels, on='cell_id', validate='one_to_one')
        sensitivity = tiers['sensitivity'].rename(columns={'global_cell_id': 'cell_id'}).merge(
            projection[['cell_id', 'assigned_cluster', 'maximum_cluster_probability', 'normalized_cluster_entropy', 'projection_status']],
            on='cell_id', validate='one_to_one')
        sensitivity['projected_original_state_id'] = study + '::studywide::' + sensitivity.assigned_cluster
        primary_parts.append(primary)
        sensitivity_parts.append(sensitivity)
        technical_parts.append(tiers['technical'].rename(columns={'global_cell_id': 'cell_id'}))
        manifests.append(manifest)
        release = json.loads((scope / 'release/scope_release.json').read_text())
        study_rows.append(dict(study_id=study, libraries=len(manifest), primary=len(primary), sensitivity=len(sensitivity),
                               technical=len(tiers['technical']), states=partition.cluster_id.nunique(), clustering_status=release['status']))
        for p in [manifest_path, qc / 'independent_validation.tsv', qc / 'tier_counts_and_hashes.tsv',
                  scope / 'frozen/cell_clusters.tsv.gz', scope / 'sensitivity/projection/sensitivity_projection.tsv.gz',
                  scope / 'release/scope_release.json']:
            provenance.append({'path': str(p), 'sha256': digest(p)})

    manifest = pd.concat(manifests, ignore_index=True)
    manifest['library_key'] = manifest.study_id + '::' + manifest.library_id
    assert manifest.library_key.is_unique
    derived = pd.DataFrame([metadata(r) for r in manifest.itertuples(index=False)])
    manifest = pd.concat([manifest, derived], axis=1)
    # A source library can support several route comparisons without duplicating expression rows.
    route = pd.read_csv(config['route_registry'], sep='\t', dtype=str,
                        usecols=['route_cell_id', 'dataset', 'study_id', 'cell_id', 'library_id', 'qc_tier', 'cluster_id'])
    route = route[route.study_id.eq('GSE178325')].copy()
    assert route.route_cell_id.is_unique
    primary = pd.concat(primary_parts, ignore_index=True)
    sensitivity = pd.concat(sensitivity_parts, ignore_index=True)
    technical = pd.concat(technical_parts, ignore_index=True)
    assert set(route.cell_id) <= set(pd.concat([primary, sensitivity]).cell_id)
    route['route_local_state_id'] = route.dataset + '::' + route.cluster_id
    route.to_csv(target / 'chemical_route_memberships.tsv.gz', sep='\t', index=False)
    route[['dataset', 'study_id', 'library_id']].drop_duplicates().to_csv(target / 'chemical_route_libraries.tsv', sep='\t', index=False)
    edges = temporal_plan(manifest, route)
    assert not edges.route_id.eq('GSE150311').any()
    assert edges.loc[edges.relation.ne('actual_day'), 'elapsed_days'].isna().all()
    edges.to_csv(target / 'temporal_edge_plan.tsv', sep='\t', index=False)
    safe = ['study_id', 'library_id', 'library_key', 'sample_accession', 'sample_title', 'experimental_time', 'condition',
            'species', 'series', 'time_kind', 'actual_day', 'stage', 'stage_order', 'within_stage_day', 'donor', 'pool', 'capture', 'numeric_time_eligible', 'reference_stratum']
    primary = primary.merge(manifest[safe], on=['study_id', 'library_id'], validate='many_to_one')
    sensitivity = sensitivity.merge(manifest[safe], on=['study_id', 'library_id'], validate='many_to_one')
    all_ids = pd.concat([primary.cell_id, sensitivity.cell_id, technical.cell_id])
    assert all_ids.is_unique
    assert len(primary) == config['expected_primary'] and len(sensitivity) == config['expected_sensitivity']
    assert len(technical) == config['expected_technical'] and len(manifest) == config['expected_libraries']
    assert primary.original_state_id.nunique() == config['expected_studywide_states']
    all_memberships = pd.concat([primary, sensitivity, technical], ignore_index=True)
    memberships_by_library = all_memberships.groupby(['study_id', 'library_id'], sort=False)
    common, feature_rows, gene_symbols, matrix_provenance = None, [], {}, []
    for row in manifest.itertuples(index=False):
        path = Path(row.input_path)
        assert path.is_file() and row.input_adapter == 'cellranger_filtered_h5'
        with h5py.File(path, 'r') as handle:
            g = handle['matrix']
            keep = decode(g['features/feature_type'][:]) == 'Gene Expression'
            genes = decode(g['features/id'][:])[keep]
            symbols = decode(g['features/name'][:])[keep]
            barcodes = decode(g['barcodes'][:])
            assert pd.Index(genes).is_unique and pd.Index(barcodes).is_unique
            assert all(str(v).startswith('ENSG') for v in genes)
            selected = memberships_by_library.get_group((row.study_id, row.library_id))
            assert set(selected.barcode) == set(barcodes)
            common = set(genes) if common is None else common & set(genes)
            for gene, symbol in zip(genes, symbols):
                gene_symbols.setdefault(gene, set()).add(symbol)
            feature_rows.append({'library_key': row.library_key, 'n_gex_features': len(genes), 'n_barcodes': len(barcodes),
                                 'gene_order_sha256': hashlib.sha256(('\n'.join(genes) + '\n').encode()).hexdigest()})
        matrix_provenance.append({'path': str(path), 'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
    assert len(common) >= config['model']['n_hvg']
    features = pd.DataFrame([{'gene_id': g, 'gene_symbol': sorted(gene_symbols[g])[0],
                              'observed_symbols': '|'.join(sorted(gene_symbols[g]))} for g in sorted(common)])
    reference = balanced_sample(primary, config['reference_cells_per_stratum'], config['seed'])
    assert reference.cell_id.is_unique and set(reference.library_key) == set(primary.library_key)
    primary['is_reference_cell'] = primary.cell_id.isin(reference.cell_id)
    primary['qc_tier'] = 'primary_pass'
    sensitivity['qc_tier'] = 'sensitivity_only'
    technical['qc_tier'] = 'technical_exclude'
    for name, frame in [('primary_cells', primary), ('sensitivity_cells', sensitivity), ('technical_exclude_cells', technical), ('reference_cells', reference)]:
        frame.sort_values('cell_id').to_csv(target / f'{name}.tsv.gz', sep='\t', index=False, compression={'method': 'gzip', 'mtime': 0})
    manifest['primary_cells'] = manifest.library_key.map(primary.groupby('library_key').size())
    manifest['sensitivity_cells'] = manifest.library_key.map(sensitivity.groupby('library_key').size()).fillna(0).astype(int)
    manifest['reference_cells'] = manifest.library_key.map(reference.groupby('library_key').size()).fillna(0).astype(int)
    manifest.to_csv(target / 'libraries.tsv', sep='\t', index=False)
    features.to_csv(target / 'shared_features.tsv', sep='\t', index=False)
    pd.DataFrame(feature_rows).to_csv(target / 'raw_feature_audit.tsv', sep='\t', index=False)
    pd.DataFrame(study_rows).to_csv(target / 'study_summary.tsv', sep='\t', index=False)
    reference.groupby(['study_id', 'reference_stratum']).size().rename('reference_cells').to_csv(target / 'reference_strata.tsv', sep='\t')
    summary = dict(status='PASS', phase='INPUT_CONTRACT_READY', primary_cells=len(primary), sensitivity_cells=len(sensitivity),
                   technical_exclude_cells=len(technical), libraries=len(manifest), studies=len(study_rows), studywide_states=primary.original_state_id.nunique(),
                   shared_genes=len(features), reference_cells=len(reference), reference_by_study=reference.groupby('study_id').size().to_dict(),
                   time_kinds=manifest.time_kind.value_counts().to_dict(),
                   chemical_route_rows=len(route), chemical_route_unique_cells=route.cell_id.nunique(),
                   planned_temporal_edges=len(edges), temporal_plans_solved=False,
                   chemical_route_duplicates_in_primary_input=0, all_primary_ids_unique=True, all_raw_barcodes_reconciled=True,
                   config_sha256=digest(config_path), joint_model_trained=False)
    (target / 'source_provenance.json').write_text(json.dumps({'files': provenance, 'raw_matrix_stats': matrix_provenance}, indent=2) + '\n')
    (target / 'preflight.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
