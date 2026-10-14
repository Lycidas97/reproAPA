"""Validate the current code inventory without loading scientific inputs."""

import csv
from hashlib import sha256
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN = {'.bam', '.bai', '.fastq', '.fq', '.h5', '.h5ad', '.loom', '.mtx',
             '.npy', '.npz', '.parquet', '.pickle', '.pkl', '.rds', '.xlsx',
             '.pdf', '.png', '.svg', '.gz', '.zip', '.tar'}
METADATA_ROLES = {'configuration', 'contract', 'workbook_contract',
                  'display_configuration', 'interface_fixture', 'model_feature_list',
                  'resource_registry'}
FIELDS = ['module', 'repository_path', 'audited_source_path', 'bytes', 'sha256']


def rows(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def safe_path(value):
    relative = Path(value)
    assert not relative.is_absolute() and '..' not in relative.parts, value
    assert relative.parts, value
    return ROOT / relative


def validate_current_edition():
    edition = json.loads((ROOT / 'provenance/current-edition.json').read_text())
    scope = edition['scope']
    assert scope == dict(main_figures=4, supplementary_figures=21,
                         supplementary_tables=10, workbook_data_sheets=58,
                         workbook_total_sheets=60, source_blocks=202)
    package = json.loads(safe_path(edition['data_package']).read_text())
    assert package['version'] == '2026-10-13 v11f'
    assert package['reserved_doi'] == '10.5281/zenodo.23347745'
    files = {r['filename']: r for r in package['files']}
    assert len(files) == len(package['files']) == 87
    assert len([name for name in files if name.endswith('.tar.gz')]) == package['archive_count'] == 81
    assert package['companion_count'] == 6
    for name, field in [('FILES.tsv', 'files_manifest'), ('ARCHIVES.tsv', 'archives_manifest')]:
        path = safe_path(package[field])
        assert path.stat().st_size == files[name]['bytes']
        assert sha256(path.read_bytes()).hexdigest() == files[name]['sha256']
    frozen = {r['path']: r for r in rows(safe_path(package['files_manifest']))}
    assert len(frozen) == 2964
    verification = json.loads(safe_path(package['local_verification']).read_text())
    assert verification['status'] == 'PASS' and verification['all_restored_content_sha256_verified']
    assert verification['verified_member_count'] == len(frozen)
    outputs = rows(safe_path(edition['outputs']))
    assert len(outputs) == len({r['output'] for r in outputs}) == 35
    expected = {f'Figure {i}' for i in range(1, 5)} | {f'Figure S{i}' for i in range(1, 22)}
    expected |= {f'Table S{i}' for i in range(1, 11)}
    assert {r['output'] for r in outputs} == expected
    for row in outputs:
        assert safe_path(row['producer']).is_file(), row
        assert row['sha256'] == frozen['project/' + row['workspace_path']]['sha256'], row
    workbook = edition['workbook']
    assert len(workbook['sheets']) == len(set(workbook['sheets'])) == 60 and 'S4j' in workbook['sheets']
    assert len([n for n in workbook['sheets'] if n not in {'Contents', 'Sources'}]) == 58
    assert workbook['sha256'] == frozen['project/' + workbook['path']]['sha256']
    assert workbook['base_sha256'] == frozen['project/' + workbook['base_workbook']]['sha256']
    assert workbook['blocks'] == 202 and safe_path(workbook['builder']).is_file()
    assert all(r['sha256'] == workbook['sha256'] for r in outputs if r['role'] == 'supplementary_table')
    for item in edition['updated_sources'] + [dict(path=p) for p in edition['added_sources']]:
        path = safe_path(item['path'])
        record = frozen['project/' + item['path'].removeprefix('code/')]
        assert path.stat().st_size == int(record['bytes'])
        assert sha256(path.read_bytes()).hexdigest() == record['sha256'], item
    current_hashes = rows(safe_path(edition['manuscript_hashes']))
    assert len({r['path'] for r in current_hashes}) == len(current_hashes)
    assert all(re.fullmatch(r'[0-9a-f]{64}', r['sha256']) and int(r['bytes']) > 0 for r in current_hashes)
    return scope


def main():
    modules = sorted((ROOT / 'modules').glob('[0-9][0-9]_*'))
    assert len(modules) == 10, 'Expected ten module entry guides'
    for module in modules:
        assert (module / 'README.md').is_file(), module
    for module in modules[:4]:
        assert (module / 'code').is_dir(), module
    assert not (ROOT / 'revisions').exists(), 'Obsolete release snapshots remain'

    manifest = rows(ROOT / 'provenance/source_manifest.tsv')
    recorded = {row['repository_path']: row for row in manifest}
    assert len(recorded) == len(manifest), 'Duplicate source paths'
    selection = json.loads((ROOT / 'provenance/source-selection.json').read_text())['files']
    assert len(selection) == len(manifest), 'Selection and manifest differ'
    for row in selection:
        assert {key: row[key] for key in FIELDS} == recorded[row['repository_path']], row
    for row in manifest:
        path = safe_path(row['repository_path'])
        assert path.is_file() and not path.is_symlink(), path
        payload = path.read_bytes()
        assert len(payload) == int(row['bytes']), path
        assert sha256(payload).hexdigest() == row['sha256'], path

    files = [p for p in ROOT.rglob('*') if '.git' not in p.relative_to(ROOT).parts
             and '__pycache__' not in p.parts and p.is_file()]
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        assert not path.is_symlink(), path
        assert path.suffix.lower() not in FORBIDDEN, path
        if path.suffix.lower() in {'.tsv', '.csv'}:
            assert (relative.startswith(('docs/', 'provenance/')) or
                    (relative in recorded and recorded[relative]['module'] in METADATA_ROLES)), path
        if relative.startswith(('code/', 'modules/', 'environments/')) and path.name != 'README.md':
            assert relative in recorded, f'Unregistered source: {relative}'
    for path in files:
        if path.suffix == '.py':
            compile(path.read_bytes(), str(path), 'exec')

    overrides = json.loads((ROOT / 'config/source-overrides.json').read_text())
    for item in overrides.values():
        path = safe_path(item['path'])
        assert sha256(path.read_bytes()).hexdigest() == item['sha256'], path
    release = json.loads((ROOT / 'provenance/release.json').read_text())
    scope = release['scope']
    assert scope == dict(main_figures=4, supplementary_figures=21,
                         supplementary_tables=10, workbook_data_sheets=57,
                         workbook_total_sheets=59)
    outputs = rows(safe_path(release['outputs']))
    assert len(outputs) == len({row['output'] for row in outputs}) == 35
    expected = {f'Figure {i}' for i in range(1, 5)}
    expected |= {f'Figure S{i}' for i in range(1, 22)}
    expected |= {f'Table S{i}' for i in range(1, 11)}
    assert {row['output'] for row in outputs} == expected
    for row in outputs:
        assert safe_path(row['producer']).is_file(), row
    for role, key in [('main_figure', 'main_figures'),
                      ('supplementary_figure', 'supplementary_figures'),
                      ('supplementary_table', 'supplementary_tables')]:
        assert sum(row['role'] == role for row in outputs) == scope[key], role
    assert len([p for p in (ROOT / 'code/analysis_code/pas_rerun').glob('test_*.py')]) == release['pas_unit_suite_files'] == 35
    registry = json.loads((ROOT / 'environments/registry.json').read_text())
    assert {n for component in registry['components'] for n in component['module_guides']} == set(range(1, 11))
    assert len({component['component'] for component in registry['components']}) == len(registry['components'])
    for component in registry['components']:
        assert component['status'] and component['limits'], component
        for item in component['evidence']:
            path = safe_path(item['path'])
            assert sha256(path.read_bytes()).hexdigest() == item['sha256'], item
    environments = json.loads((ROOT / 'environments/reconstruction.json').read_text())
    definitions = environments['environments']
    for component in registry['components']:
        assert component['reproduction']['environment'] in definitions, component
        target = definitions[component['reproduction']['environment']]
        assert component['reproduction']['python'] == target['python'], component['component']
        assert component['reproduction']['versions'] == target['versions'], component['component']
        for snapshot in component['current_snapshots']:
            assert snapshot in environments['current_snapshots'], snapshot
    for name, definition in definitions.items():
        assert definition['python'] and definition['versions'] is not None, name
        if definition['kind'] != 'resolved-reconstruction':
            continue
        assert definition['resolution']['status'] == 'PASS', name
        lock = safe_path(definition['requirements_lock'])
        assert sha256(lock.read_bytes()).hexdigest() == definition['lock_sha256'], name
        pins = dict(re.findall(r'^([A-Za-z0-9_.-]+)==([^\s\\]+)', lock.read_text(), re.M))
        assert pins == definition['versions'], name
        assert all(pins[key] == value for key, value in definition['known_pins'].items()), name
        for line in lock.read_text().replace('\\\n', '').splitlines():
            if line.strip() and not line.startswith('#'):
                assert re.search(r'--hash=sha256:[0-9a-f]{64}', line), (name, line)
    releases = json.loads(safe_path(environments['public_release_evidence']).read_text())
    for package in releases['packages']:
        assert any(not artifact['yanked'] for artifact in package['files']), package['name']
        assert all(re.fullmatch(r'[0-9a-f]{64}', artifact['sha256']) for artifact in package['files'])
    pas = json.loads((ROOT / 'environments/evidence/pas-runtime-and-structural-cache.json').read_text())
    assert pas['cache_lock_hashes_match'] and pas['receipt_and_controller_hashes_match_current_config']
    assert pas['legacy_policy_constants_identical']
    assert all(all(row.values()) for row in pas['policy_constants_and_gtf_signatures']['gtf_signature_comparison'].values())
    index = json.loads(safe_path(release['source_data']).read_text())
    data = json.loads((ROOT / 'provenance/source-data-selection.json').read_text())
    assert index['status'] == 'PASS' and index['member_crc_and_sha256_verified']
    assert index['source_members'] == len(data['files'])
    assert index['members'] == len(data['files']) + 4
    assert index['coverage'] == data['coverage_counts']
    for data_key, scope_key in [('main_figures', 'main_figures'),
                                ('supplementary_figures', 'supplementary_figures'),
                                ('tables', 'supplementary_tables'),
                                ('data_sheets', 'workbook_data_sheets')]:
        assert index['coverage'][data_key] == scope[scope_key], data_key
    inputs = json.loads(safe_path(data['input_archives']).read_text())['archives']
    assert index['input_archives'] == data['input_archives']
    archives = {row['path']: row for row in inputs}
    assert len(archives) == len(inputs)
    assert index['companion'] == {k: data['companion'][k] for k in index['companion']}
    assert index['additional_companions'] == data['additional_companions']
    for companion in [data['companion'], *data['additional_companions']]:
        recorded = archives[companion['path']]
        assert all(recorded[k] == companion[k] for k in ['bytes', 'sha256', 'members']), companion
    input_files = {row['path']: row for archive in inputs for row in archive.get('files', [])}
    sources = {row['path']: row for row in data['files']}
    assert len(sources) == len(data['files'])
    assert len({row['archive_path'] for row in data['files']}) == len(data['files'])
    assert data['coverage_counts']['required_input_sources'] == len(data['required_input_sources'])
    for row in data['required_input_sources']:
        assert row['sha256'] == input_files[row['path']]['sha256'], row
        assert row['bytes'] == input_files[row['path']]['bytes'], row
        sources[row['path']] = row
    for row in outputs:
        assert row['sha256'] == sources[row['workspace_path']]['sha256'], row
    workbook = release['workbook']
    sheets = workbook['sheets']
    assert len(sheets) == len(set(sheets)) == scope['workbook_total_sheets']
    assert len([name for name in sheets if name not in {'Contents', 'Sources'}]) == scope['workbook_data_sheets']
    assert workbook['sha256'] == sources[workbook['path']]['sha256']
    assert workbook['source_base_sha256'] == sources[workbook['source_base']]['sha256']
    assert all(row['sha256'] == workbook['sha256'] for row in outputs if row['role'] == 'supplementary_table')
    assert release['recorded_checks']['workbook']['problems'] == []
    assert release['recorded_checks']['workbook']['checks']['sheets'] == len(sheets)
    assert all(release['recorded_checks']['main_figures']['set_copies'].values())
    checks = json.loads((ROOT / 'provenance/verification.json').read_text())
    assert checks['unit_tests']['failures'] == checks['unit_tests']['errors'] == checks['unit_tests']['skipped'] == 0
    assert checks['followup_unit_tests']['failures'] == checks['followup_unit_tests']['errors'] == checks['followup_unit_tests']['skipped'] == 0
    assert checks['workbook']['data_sheets'] == scope['workbook_data_sheets']
    assert checks['workbook']['total_sheets'] == scope['workbook_total_sheets']
    assert checks['workbook']['sha256'] == workbook['sha256']
    with safe_path(release['manuscript_hashes']).open(newline='') as stream:
        manuscript = list(csv.DictReader(stream, delimiter='\t'))
    assert len({row['path'] for row in manuscript}) == len(manuscript)
    assert all(re.fullmatch(r'[0-9a-f]{64}', row['sha256']) and int(row['bytes']) > 0 for row in manuscript)
    for path in [ROOT / 'README.md', ROOT / 'code/README.md',
                 *(ROOT / 'docs').glob('*.md'), *(m / 'README.md' for m in modules)]:
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            if re.match(r'^[a-z]+://', target) or target.startswith('#'):
                continue
            assert (path.parent / target.split('#')[0]).exists(), (path, target)
    current_scope = validate_current_edition()
    print(f'PASS: {len(manifest)} source records, '
          f'{sum(p.suffix == ".py" for p in files)} Python files; '
          f'{current_scope["main_figures"]} main figures, '
          f'{current_scope["supplementary_figures"]} supplementary figures, '
          f'{current_scope["supplementary_tables"]} tables / '
          f'{current_scope["workbook_data_sheets"]} data sheets; current v11f and historical source/environment hashes')



if __name__ == '__main__':
    main()
