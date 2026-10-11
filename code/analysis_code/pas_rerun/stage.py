#!/usr/bin/env python3
"""Assemble a reviewable manuscript from a single verified rerun, never in place."""
import argparse
import csv
import json
import shutil
from pathlib import Path

from intake import HERE, artifact, digest, require
from run import ROOT, audit_baseline, load_workflow, verify_result


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def check_package(package_path, run_dir):
    package_path=Path(package_path).resolve(); run_dir=Path(run_dir).resolve()
    pkg=json.loads(package_path.read_text()); state=json.loads((run_dir/'run.json').read_text())
    require(state['status']=='complete' or (state['status']=='running' and state.get('active_node')=='release'), 'Run has not reached the release step')
    require(pkg.get('run_id')==state['run_id'] and pkg.get('input_fingerprint')==state['input_fingerprint'], 'Package/run mismatch')
    nodes=[n['id'] for n in load_workflow(HERE/'workflow.json')['nodes'] if n['id']!='release']
    require(set(nodes)<=set(state['completed']), 'Not all analysis/render/table/text nodes completed')
    produced=set()
    for n in nodes:
        receipt=run_dir/n/'result.json'
        result=verify_result(receipt,n,state['run_id'],state['input_fingerprint'])
        produced.update((receipt.parent/a['path']).resolve() for a in result['artifacts'])
    files=pkg.get('files',[]);targets=[r['target'] for r in files]
    require(len(set(targets))==len(targets), 'Duplicate replacement targets')
    expected_figs={r['target'] for r in read_tsv(HERE/'figures.tsv') if r['action']=='rebuild'}
    require({r['target'] for r in files if r['kind']=='figure'}==expected_figs, 'Missing/extra PAS figure replacements')
    require(set(range(8,31))=={int(v) for r in files if r['kind']=='table' for v in r.get('table_ids',[])}, 'Tables S8–S30 must all have new data assets')
    allowed_sources={r['file'] for r in read_tsv(HERE/'baseline-files.tsv') if r['file'].startswith(('main/','supplement/')) and r['file'].endswith('.tex')}
    require({r['target'] for r in files if r['kind']=='split_source'}<=allowed_sources, 'Only existing split sources may be replaced; structural edits require tex_sections separately')
    for r in files:
        target=Path(r['target'])
        require(not target.is_absolute() and '..' not in target.parts, 'Unsafe replacement target')
        require(r['kind'] in {'figure','table','split_source','bibliography'}, 'Unknown replacement kind')
        if r['kind']=='bibliography':
            require(r['target']=='references.bib', 'Candidate bibliography must replace references.bib')
        if r['kind']=='table':
            require(r['target'].startswith('table/'), 'Table assets must be under table/')
        p=artifact(r,package_path.parent)
        require(p in produced, f'Unregistered or historical result in package: {p}')
        if r['kind']=='figure':
            with p.open('rb') as f:
                require(f.read(5)==b'%PDF-',f'Not a PDF: {p}')
    review=artifact(pkg['text_review'],package_path.parent)
    require(review in produced, 'Text-review ledger not in this run')
    expected={r['id']:r for r in read_tsv(HERE/'text-targets.tsv')}
    rows=read_tsv(review)
    require(len(rows)==len(expected) and {r['id'] for r in rows}==set(expected),'Text review has missing/duplicate entries')
    for row in rows:
        require(row['baseline_sha256']==expected[row['id']]['sha256'],'Text review baseline mismatch')
        require(row['status'] in {'rewritten','verified_unchanged','not_result_verified'} and row.get('evidence'), 'Unresolved text/claim review')
        if row['status']=='rewritten':
            require(expected[row['id']]['file'] in targets, 'Rewritten paragraph lacks new split source')
    return pkg


def stage(package_path, run_dir, destination):
    audit_baseline()
    pkg=check_package(package_path,run_dir)
    dest=Path(destination).resolve()
    require(not dest.exists(), 'Candidate destination already exists')
    require(dest!=ROOT and not ROOT.is_relative_to(dest),'Candidate must not replace the project')
    dest.mkdir(parents=True)
    for folder in ['main','supplement','Definitions']:
        shutil.copytree(ROOT/folder,dest/folder)
    (dest/'scripts').mkdir();shutil.copy2(ROOT/'scripts/tex_sections.py',dest/'scripts/tex_sections.py')
    for name in ['Makefile','main.tex','supplement.tex','references.bib']:
        shutil.copy2(ROOT/name,dest/name)
    (dest/'Fig').mkdir();(dest/'table').mkdir()
    registry=read_tsv(ROOT/'Fig/figure-manifest.tsv')
    for row in registry:
        shutil.copy2(ROOT/'Fig'/row['file'],dest/'Fig'/row['file'])
    # Reusable metadata TSVs only: never carry a combined workbook containing old S8–S15.
    for p in (ROOT/'table').glob('S[1-7]-*.tsv'):
        shutil.copy2(p,dest/'table'/p.name)
    for row in pkg['files']:
        source=artifact(row,Path(package_path).resolve().parent)
        target=dest/row['target'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        if row['kind']=='figure':
            entry=next(r for r in registry if 'Fig/'+r['file']==row['target'])
            entry['sha256']=row['sha256'];entry['source_candidate']=str(source)
            entry['title']=row.get('title',entry['title'])
    with (dest/'Fig/figure-manifest.tsv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(registry[0]),delimiter='\t');writer.writeheader();writer.writerows(registry)
    (dest/'PAS-RERUN.json').write_text(json.dumps({'run_id':pkg['run_id'],'input_fingerprint':pkg['input_fingerprint'],'status':'staged_not_compiled','package_sha256':digest(package_path),'required_checks':['make check','make pdf','main and supplement page visual review']},indent=2)+'\n')
    return dest


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    print(stage(args.package,args.run,args.out))
