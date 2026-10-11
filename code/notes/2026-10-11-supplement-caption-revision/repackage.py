"""Refresh the existing v11b selections with current files after caption revision."""
import importlib.util
import json
from pathlib import Path
import shutil
import zipfile

RUN = Path(__file__).resolve().parent
ROOT = RUN.parents[1]
spec = importlib.util.spec_from_file_location(
    'packager', ROOT / 'github/reproapa-analysis/scripts/package_source_data.py')
pack = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pack)
NAMES = ('Supplementary_Source_Code_v11b.zip',
         'Supplementary-Source-Data-20261011-v11b.zip')
STAGE = RUN / 'repacked'
BACKUP = RUN / 'before-archives'


def artifacts(archive):
    return (archive, archive.with_suffix('.manifest.tsv'),
            archive.with_suffix('.verification.json'),
            archive.with_suffix('.zip.sha256'))


def selection(archive):
    with zipfile.ZipFile(archive) as z:
        return json.loads(z.read('SOURCE-SELECTION.json'))


def bind(relative):
    path = ROOT / relative
    with zipfile.ZipFile(path) as z:
        count = len(z.namelist())
    return dict(path=relative, bytes=path.stat().st_size,
                sha256=pack.digest(path), members=count)


def refresh(selected, additions):
    rows = {row['path']: row for row in selected['files']}
    for relative, role in additions:
        rows.setdefault(relative, dict(path=relative, archive_path='project/' + relative,
                                       roles=[role], owners=['supplement caption revision']))
    for relative, row in rows.items():
        path = ROOT / relative
        row.update(bytes=path.stat().st_size, sha256=pack.digest(path))
    selected['files'] = sorted(rows.values(), key=lambda row: row['path'])
    if 'source_files' in selected['coverage_counts']:
        selected['coverage_counts']['source_files'] = len(rows)


def build(selected, name):
    pin = RUN / 'work' / (name + '.selection.json')
    pin.write_text(json.dumps(selected, indent=2) + '\n')
    output = STAGE / name
    pack.package(ROOT, pin, output)
    for path in artifacts(output):
        shutil.copy2(path, ROOT / 'table' / path.name)


if __name__ == '__main__':
    # Never silently overwrite an earlier staged run or its pre-revision archives.
    assert not BACKUP.exists(), BACKUP
    assert not any((STAGE / name).exists() for name in NAMES), STAGE
    code, data = [selection(ROOT / 'table' / name) for name in NAMES]
    BACKUP.mkdir()
    for name in NAMES:
        for path in artifacts(ROOT / 'table' / name):
            shutil.copy2(path, BACKUP / path.name)
    original_workbook = pack.digest(ROOT / 'table/Supplementary_Tables_S1-S10.xlsx')
    historical = bind('table/Supplementary-Source-Data-20261009.zip')
    s9 = bind('table/Supplementary_Source_Code_S9.zip')
    prefix = RUN.relative_to(ROOT).as_posix() + '/'
    try:
        code['scope'] += ' Includes the final supplementary caption and native-state ID revision.'
        code['companion'] = historical
        code['additional_companions'] = [s9]
        code['archive_readme'] = code['archive_readme'].replace(
            'The unchanged Supplementary_Source_Code_S9.zip is a required companion, bound\nby COMPANION.json.',
            'The unchanged 9 October source-data archive is bound by COMPANION.json.\n'
            'Supplementary_Source_Code_S9.zip is bound by additional_companions in\n'
            'SOURCE-SELECTION.json and remains required.')
        code['archive_readme'] += (
            '\nThe final supplementary caption revision defines states and statistics and\n'
            'restores the sensitivity pointer to current Table S5d.6/S5d.7. S1/S2 now\n'
            'show native state IDs. Their current compose.py supersedes its earlier\n'
            'captured version; notes/2026-10-11-supplement-caption-revision/render.py\n'
            'rebuilds only these two figures without promoting them to Fig/.\n')
        refresh(code, [(prefix + 'render.py', 'figure_rendering_entrypoint'),
                       (prefix + 'repackage.py', 'packaging_entrypoint')])
        build(code, NAMES[0])

        data['scope'] += ' Final supplementary captions and S1/S2 state-ID assets.'
        data['companion'] = historical
        data['additional_companions'] = [bind('table/' + NAMES[0]), s9]
        data['archive_readme'] = data['archive_readme'].replace(
            'No supplementary figure changed.',
            'S1/S2 subsequently gained native state column IDs, with all numerical\n'
            'values preserved. The other 19 supplementary figure assets are unchanged.')
        data['archive_readme'] += (
            '\nThe final caption revision is recorded in\n'
            'notes/2026-10-11-supplement-caption-revision/README.md. It defines state\n'
            'codes, statistics and evidence tiers, and restores the sensitivity\n'
            'pointer to Table S5d.6/S5d.7. The already verified 200-block workbook\n'
            'is unchanged by this revision. New S1/S2 rendering receipts and values\n'
            'are under that revision directory. Older candidate receipts retain\n'
            'their original identity; Fig/figure-manifest.tsv pins production assets.\n'
            'Packaging follows this rendering and does not rerun scientific analyses.\n')
        additions = [(prefix + 'README.md', 'caption_revision_record'),
                     ('table/README.md', 'current_source_archive_guide'),
                     ('Fig/README.md', 'figure_provenance_guide'),
                     (prefix + 'work/figure-qa.json', 'figure_revision_verification'),
                     (prefix + 'work/state-column-map.tsv', 'native_state_column_mapping')]
        for key, number in (('A', 1), ('B', 2)):
            for suffix in ('pdf', 'receipt.json', 'values.json'):
                additions.append((prefix + f'rendered/{key}/SuppFig.S{number}.{suffix}',
                                  'current_S1_S2_rendering_record'))
        refresh(data, additions)
        assert pack.digest(ROOT / 'table/Supplementary_Tables_S1-S10.xlsx') == original_workbook
        build(data, NAMES[1])
    except BaseException:
        for path in BACKUP.iterdir():
            shutil.copy2(path, ROOT / 'table' / path.name)
        raise

    # Both ZIPs have passed streamed member CRC/SHA checks before promotion.
    assert pack.digest(ROOT / 'table/Supplementary_Tables_S1-S10.xlsx') == original_workbook
    records = []
    for name in NAMES:
        receipt = json.loads((ROOT / 'table' / name).with_suffix('.verification.json').read_text())
        records.append(dict(archive=name, sha256=receipt['sha256'], members=receipt['members']))
    (RUN / 'work/repackage-result.json').write_text(json.dumps(dict(
        status='PASS', workbook_sha256=original_workbook, archives=records,
        pre_revision_archives=str(BACKUP.relative_to(ROOT))), indent=2) + '\n')
