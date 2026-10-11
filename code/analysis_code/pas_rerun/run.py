#!/usr/bin/env python3
"""Fail-closed orchestration for the PAS rerun. No implicit legacy commands."""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from intake import HERE, artifact, digest, require, validate

ROOT = HERE.parents[1]


def load_workflow(path):
    data = json.loads(Path(path).read_text())
    seen = set()
    for n in data['nodes']:
        require(re.fullmatch(r'[a-z][a-z0-9_]*', n['id']) and n['id'] not in seen, 'Invalid/duplicate node')
        require(set(n['depends_on']) <= seen, f'Workflow not topologically ordered: {n["id"]}')
        seen.add(n['id'])
    return data


def binding_errors(workflow, bindings, base):
    problems = []
    nodes = bindings.get('nodes', {})
    expected = {n['id'] for n in workflow['nodes']}
    if set(nodes) != expected:
        problems.append('Binding IDs do not exactly cover the full workflow')
    for n in workflow['nodes']:
        b = nodes.get(n['id'], {})
        if not isinstance(b.get('argv'), list) or not b.get('argv') or not all(isinstance(a,str) and a for a in b['argv']):
            problems.append(n['id']+': no adapted command')
        if not b.get('code_files'):
            problems.append(n['id']+': no pinned adapter code')
        for code in b.get('code_files', []):
            try:
                artifact(code, base)
            except (ValueError, OSError) as exc:
                problems.append(n['id']+': '+str(exc))
    return problems


def verify_result(path, node, run_id, fingerprint):
    r = json.loads(path.read_text())
    require(r.get('node') == node and r.get('run_id') == run_id and r.get('input_fingerprint') == fingerprint,
            f'{node}: result belongs to different node/run/inputs')
    require(r.get('status') == 'complete', f'{node}: failed, skipped or incomplete result')
    require(isinstance(r.get('artifacts'), list) and r['artifacts'], f'{node}: empty artifact list')
    require(isinstance(r.get('checks'), dict) and r['checks'] and all(v is True for v in r['checks'].values()), f'{node}: verification missing/failed')
    require({'source_data', 'audit'} <= {a.get('role') for a in r['artifacts']}, f'{node}: source_data and audit artifacts are required')
    for a in r['artifacts']:
        require(a.get('role'), f'{node}: missing artifact role')
        p = artifact(a, path.parent)
        require(p.is_relative_to(path.parent.resolve()), f'{node}: output must be inside its own run directory')
    return r


def execute(workflow, bindings, manifest_path, run_dir, *, scope_path=None):
    report = validate(manifest_path, scope_path=scope_path)
    p = Path(run_dir).resolve()
    require(not p.exists(), f'Run directory already exists: {p}; use a new directory (no stale-output reuse)')
    run_id = p.name
    # Content fingerprint includes all bound data, workflow and adapter definitions.
    import hashlib
    fingerprint = hashlib.sha256(json.dumps({'inputs':report['files'], 'manifest':digest(manifest_path), 'workflow':workflow, 'bindings':bindings}, sort_keys=True).encode()).hexdigest()
    p.mkdir(parents=True)
    (p/'inputs.validated.json').write_text(json.dumps(report, indent=2)+'\n')
    receipts = {}
    status = dict(status='running', run_id=run_id, input_fingerprint=fingerprint, completed=[])

    def checkpoint():
        (p/'run.json').write_text(json.dumps(status, indent=2)+'\n')

    checkpoint()
    try:
        for n in workflow['nodes']:
            name = n['id']; dest = p/name; dest.mkdir()
            deps = {key: str(receipts[key]) for key in n['depends_on']}
            for key in deps:
                verify_result(receipts[key], key, run_id, fingerprint)
            (dest/'dependencies.json').write_text(json.dumps(deps, indent=2)+'\n')
            env = dict(os.environ, PAS_INPUT_MANIFEST=str(Path(manifest_path).resolve()), PAS_RUN_ID=run_id,
                       PAS_INPUT_FINGERPRINT=fingerprint, PAS_NODE_ID=name, PAS_NODE_DIR=str(dest),
                       PAS_RUN_DIR=str(p), PAS_DEPENDENCIES_JSON=str(dest/'dependencies.json'), PAS_PROJECT_ROOT=str(ROOT))
            status['active_node'] = name; checkpoint()
            print('RUN', name, flush=True)
            with (dest/'command.log').open('w') as log:
                subprocess.run(bindings['nodes'][name]['argv'], cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            receipt = dest/'result.json'
            verify_result(receipt, name, run_id, fingerprint)
            receipts[name] = receipt; status['completed'].append(name); checkpoint()
        # Recheck every output and input at the end, not just process exit codes.
        for name, receipt in receipts.items():
            verify_result(receipt, name, run_id, fingerprint)
        for path, expected in report['files'].items():
            require(digest(path) == expected, f'Input changed during run: {path}')
        status['status'] = 'complete'; status.pop('active_node', None); checkpoint()
    except Exception as exc:
        status.update(status='failed', error=str(exc)); checkpoint()
        raise
    return status


def audit_baseline():
    with (HERE/'baseline-files.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    mismatches = [r['file'] for r in rows if not (ROOT/r['file']).is_file() or digest(ROOT/r['file']) != r['sha256']]
    require(not mismatches, 'Baseline changed; reconcile inventory before rerun: '+', '.join(mismatches))
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['plan', 'validate', 'preflight', 'run', 'audit'])
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--bindings', type=Path, default=HERE/'bindings.template.json')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    workflow = load_workflow(HERE/'workflow.json')
    if args.command == 'audit':
        print(json.dumps({'status':'passed','baseline_files':audit_baseline()})); return
    if args.command == 'plan':
        for n in workflow['nodes']:
            print(n['id']+'\t'+','.join(n['depends_on'])+'\t'+n['title'])
        return
    require(args.manifest is not None, '--manifest is required')
    if args.command == 'validate':
        print(json.dumps(validate(args.manifest), ensure_ascii=False, indent=2)); return
    bindings = json.loads(args.bindings.read_text())
    errors = binding_errors(workflow, bindings, args.bindings.resolve().parent)
    try:
        audit_baseline()
    except ValueError as exc:
        errors.append(str(exc))
    if errors:
        raise ValueError('Full rerun is NOT ready:\n'+'\n'.join(errors))
    if args.command == 'preflight':
        report = validate(args.manifest)
        print(json.dumps({'status':'ready', 'libraries':len(report['libraries']), 'nodes':len(workflow['nodes'])})); return
    require(args.out is not None, '--out is required')
    result = execute(workflow, bindings, args.manifest, args.out)
    errors = binding_errors(workflow, bindings, args.bindings.resolve().parent)
    if errors:
        result.update(status='failed', error='Adapter code changed during run: '+ '; '.join(errors))
        (args.out/'run.json').write_text(json.dumps(result, indent=2)+'\n')
        raise ValueError(result['error'])
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, ImportError, subprocess.CalledProcessError) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
