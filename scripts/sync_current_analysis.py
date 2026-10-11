"""Synchronize the pinned current source inventory; reject drift before writes."""

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import shutil


HERE = Path(__file__).resolve()
DEFAULT_REPOSITORY = (HERE.parents[1] if HERE.parent.name == 'scripts'
                      else HERE.parents[1] / 'github/reproapa-analysis')
FIELDS = ['module', 'repository_path', 'audited_source_path', 'bytes', 'sha256']


def safe_path(root, value):
    relative = Path(value)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError(f'Unsafe repository path: {value}')
    path = root / relative
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes repository: {value}')
    return path


def synchronize(workspace, repository):
    inventory = json.loads((repository / 'provenance/source-selection.json').read_text())
    with (repository / 'provenance/source_manifest.tsv').open(newline='') as stream:
        manifest = list(csv.DictReader(stream, delimiter='\t'))
    entries = inventory['files']
    projected = [{key: row[key] for key in FIELDS} for row in entries]
    if projected != manifest or len({row['repository_path'] for row in entries}) != len(entries):
        raise ValueError('Source selection and manifest differ')
    pending = []
    # Complete all source checks before copying even the first file.
    for row in entries:
        destination = safe_path(repository, row['repository_path'])
        alias = row['source_alias']
        if alias == 'PAPER':
            source = safe_path(workspace, row['source_path'])
        elif alias == 'CAPTURED':
            source = safe_path(repository, row['source_path'])
        elif alias == 'EXTERNAL' and Path(row['source_path']).is_absolute():
            source = Path(row['source_path'])
        else:
            raise ValueError(f'Unknown source alias: {alias}')
        payload = source.read_bytes()
        if len(payload) != int(row['bytes']) or sha256(payload).hexdigest() != row['sha256']:
            raise ValueError(f'Pinned source changed: {source}; revise the inventory explicitly')
        if destination.is_symlink():
            raise ValueError(f'Symlink destination: {destination}')
        pending.append((source, destination))
    release = json.loads((repository / 'provenance/release.json').read_text())
    with safe_path(repository, release['manuscript_hashes']).open(newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            source = safe_path(workspace, row['path'])
            if sha256(source.read_bytes()).hexdigest() != row['sha256']:
                raise ValueError(f'Manuscript changed: {source}; revise the edition metadata')
    copied = 0
    for source, destination in pending:
        if source.resolve() == destination.resolve():
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        copied += 1
    print(f'PASS: {len(entries)} pinned sources checked; {copied} files synchronized')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--repository', type=Path, default=DEFAULT_REPOSITORY)
    args = parser.parse_args()
    synchronize(args.workspace.resolve(), args.repository.resolve())
