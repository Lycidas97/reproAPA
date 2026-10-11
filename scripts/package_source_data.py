"""Package pinned source files and bind the required companion."""

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
import zipfile


def digest(path):
    value = sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def safe_path(root, value):
    relative = Path(value)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError(f'Unsafe source path: {value}')
    path = root / relative
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Source escapes workspace: {value}')
    return path


def package(workspace, selection_path, output):
    selection = json.loads(selection_path.read_text())
    files = selection['files']
    assert len({x['path'] for x in files}) == len(files)
    assert len({x['archive_path'] for x in files}) == len(files)
    if output.exists():
        raise FileExistsError(output)
    companion = selection['companion']
    companions = [companion, *selection.get('additional_companions', [])]
    for required in companions:
        companion_path = safe_path(workspace, required['path'])
        assert companion_path.stat().st_size == required['bytes']
        assert digest(companion_path) == required['sha256'], 'Required companion changed'
    # Validate every pinned source before creating an output.
    for row in files:
        path = safe_path(workspace, row['path'])
        safe_path(workspace, row['archive_path'])
        assert path.is_file() and not path.is_symlink(), path
        assert path.stat().st_size == row['bytes'] and digest(path) == row['sha256'], path
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.partial')
    if temporary.exists():
        raise FileExistsError(temporary)
    manifest = io.StringIO(newline='')
    fields = ['path', 'archive_path', 'bytes', 'sha256', 'roles', 'owners']
    writer = csv.DictWriter(manifest, fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows({**x, 'roles': '; '.join(x['roles']),
                     'owners': '; '.join(x['owners'])} for x in files)
    readme = selection['archive_readme']
    archive_date = tuple(int(v) for v in selection['date'].split('-')) + (0, 0, 0)
    try:
        with zipfile.ZipFile(temporary, 'w', allowZip64=True) as archive:
            for row in files:
                path = safe_path(workspace, row['path'])
                info = zipfile.ZipInfo(row['archive_path'], archive_date)
                info.compress_type = (zipfile.ZIP_STORED if path.suffix.lower() in
                                      {'.gz', '.npz', '.xlsx', '.pdf', '.png', '.zip'}
                                      else zipfile.ZIP_DEFLATED)
                info._compresslevel = 1
                info.external_attr = 0o100644 << 16
                value = sha256()
                size = 0
                with path.open('rb') as source, archive.open(info, 'w', force_zip64=True) as target:
                    for chunk in iter(lambda: source.read(8 * 1024 * 1024), b''):
                        value.update(chunk)
                        size += len(chunk)
                        target.write(chunk)
                assert size == row['bytes'] and value.hexdigest() == row['sha256'], path
            archive.writestr('SOURCE-MANIFEST.tsv', manifest.getvalue())
            archive.writestr('SOURCE-SELECTION.json', json.dumps(selection, indent=2))
            archive.writestr('COMPANION.json', json.dumps(companion, indent=2))
            archive.writestr('README.md', readme)
        # Stream every archived member: CRC is checked by zipfile and SHA by this ledger.
        with zipfile.ZipFile(temporary) as archive:
            assert len(archive.namelist()) == len(files) + 4
            assert len(set(archive.namelist())) == len(archive.namelist())
            for row in files:
                value = sha256()
                size = 0
                with archive.open(row['archive_path']) as stream:
                    for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                        value.update(chunk)
                        size += len(chunk)
                assert size == row['bytes'] and value.hexdigest() == row['sha256'], row
        for required in companions:
            companion_path = safe_path(workspace, required['path'])
            assert digest(companion_path) == required['sha256'], 'Required companion changed during packaging'
        temporary.rename(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    receipt = {'status': 'PASS', 'date': selection['date'], 'archive': output.name,
               'bytes': output.stat().st_size, 'sha256': digest(output),
               'members': len(files) + 4, 'source_members': len(files),
               'source_uncompressed_bytes': sum(x['bytes'] for x in files),
               'member_crc_and_sha256_verified': True,
               'coverage': selection['coverage_counts'],
               'companion': {k: companion[k] for k in ['path', 'bytes', 'sha256', 'members']},
               'additional_companions': selection.get('additional_companions', []),
               'full_scientific_analysis_rerun': False,
               'figure_rendering_performed': False,
               'figure_provenance_limits': selection['figure_provenance_limits']}
    output.with_suffix('.verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
    output.with_suffix(output.suffix + '.sha256').write_text(receipt['sha256'] + '  ' + output.name + '\n')
    output.with_suffix('.manifest.tsv').write_text(manifest.getvalue())
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--selection', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'provenance/source-data-selection.json')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    package(args.workspace.resolve(), args.selection.resolve(), args.out.resolve())
