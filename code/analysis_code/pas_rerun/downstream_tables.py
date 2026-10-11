"""Collect verified source tables; workbook creation is a separate local rendering step."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import zipfile

import numpy as np
import pandas as pd

from downstream_runtime import Context
from polyapipe_prepare import digest


HISTORICAL={'recount','consensus'}
HISTORICAL_FILES={'cross_trace/historical_published_support.tsv',
    'caller_comparison/legacy_recurrence_input.tsv.gz','caller_comparison/legacy_members_input.tsv.gz'}
MIXED_FILES={'cross_trace/member_migration.tsv','external_pas/historical_roster_migration.tsv',
    'coordinate_states/family_migration.tsv','examples/ten_locus_migration.tsv'}


def json_scalar(value):
    if isinstance(value,np.generic):return value.item()
    raise TypeError(f'Unsupported workbook value: {type(value).__name__}')


def workbook_input(ctx,frame,label):
    typed=ctx.out/'workbook-inputs'/f'{label}.json';typed.parent.mkdir(exist_ok=True)
    values=frame.astype(object).where(frame.notna(),None).to_numpy().tolist()
    typed.write_text(json.dumps(dict(columns=frame.columns.tolist(),data=values),allow_nan=False,
        default=json_scalar)+'\n')
    return str(typed.relative_to(ctx.out))


def verified_sources(ctx,nodes):
    artifacts={};index=[]
    for node in nodes:
        receipt=ctx.analysis/node/'result.json';record=json.loads(ctx.register(receipt).read_text())
        assert record['node']==node and record['status']=='complete' and all(record['checks'].values())
        assert record['run_id']==ctx.run_id and record['input_fingerprint']==ctx.config['input_fingerprint']
        for item in record['artifacts']:
            relative=str(Path(node)/item['path']);path=ctx.analysis/relative
            assert path.is_file() and digest(path)==item['sha256'],relative
            artifacts[relative]=(path,item['sha256'])
            scope=item.get('evidence_scope','historical_original_pipeline' if node in HISTORICAL else 'current_native_rerun')
            if relative in HISTORICAL_FILES:scope='historical_original_pipeline'
            elif relative in MIXED_FILES:scope='mixed_explicitly_labeled'
            elif node=='caller_comparison':scope='native_vs_historical_comparison'
            index.append(dict(node=node,source=relative,sha256=item['sha256'],bytes=path.stat().st_size,
                evidence_scope=scope,native_polyapipe_evidence=scope=='current_native_rerun'))
    return artifacts,pd.DataFrame(index)


def collect(ctx):
    if ctx.config.get('reference_variant')=='v2_primary':
        assert 'reference_control_source_archive' in ctx.config['resources'],'Full v3 control archive is required for this integrated-reference delivery'
    bindings=json.loads(ctx.resource('table_bindings').read_text())
    workflow=json.loads(ctx.resource('workflow').read_text())
    nodes=[r['id'] for r in workflow['nodes'] if r['id'] not in ['render','tables','text','release']]
    artifacts,source_index=verified_sources(ctx,nodes+['figure_scores','render']+ctx.config.get('extra_source_nodes',[]))
    scopes=source_index.set_index('source').evidence_scope.to_dict()
    assert set(bindings)=={str(i) for i in range(8,31) if i!=29}
    manuscript=Path(ctx.config['manuscript_root']);inventory=[];sheets=[]
    for number in range(1,8):
        files=list((manuscript/'table').glob(f'S{number}-*.tsv'));assert len(files)==1
        p=ctx.register(files[0]);target=ctx.out/'table'/p.name;target.parent.mkdir(exist_ok=True);shutil.copy2(p,target)
        inventory.append(dict(table=number,part='metadata',file=str(target.relative_to(ctx.out)),source=str(p),sha256=digest(target),evidence_scope='frozen_GEX_metadata',rows=len(pd.read_csv(p,sep='\t'))))
        sheets.append(dict(table=number,name=f'S{number}',source=workbook_input(ctx,ctx.read(target),f'S{number}'),role='complete_table',
            evidence_scope='frozen_GEX_metadata',origin='table/'+p.name,data_file=str(target.relative_to(ctx.out))))
    for number,patterns in bindings.items():
        chosen=[]
        for pattern in patterns:
            matches=sorted(ctx.analysis.glob(pattern));assert matches,f'Missing Table S{number} source: {pattern}'
            for path in matches:
                relative=str(path.relative_to(ctx.analysis));assert relative in artifacts,f'Unreceipted source: {relative}'
                if relative not in chosen:chosen.append(relative)
        for i,relative in enumerate(chosen,1):
            source,expected=artifacts[relative];node=Path(relative).parts[0]
            label=f'S{number}_{i:02d}';target=ctx.out/'table'/f'{label}-{source.name}'
            shutil.copy2(source,target);assert digest(target)==expected
            frame=ctx.read(source) if '.tsv' in source.name else pd.DataFrame([json.loads(source.read_text())])
            inventory.append(dict(table=int(number),part=label,file=str(target.relative_to(ctx.out)),source=relative,sha256=expected,
                evidence_scope=scopes[relative],rows=len(frame)))
            # Full precision machine-readable tables are always delivered. Very large gene-level tables
            # remain in those TSVs; the workbook index explicitly identifies this delivery format.
            if len(frame)*len(frame.columns)<=200000 and len(frame)<=100000:
                sheets.append(dict(table=int(number),name=label,source=workbook_input(ctx,frame,label),role='complete_table',
                    evidence_scope=inventory[-1]['evidence_scope'],origin=relative,data_file=str(target.relative_to(ctx.out))))
            else:sheets.append(dict(table=int(number),name=label,source=str(target.relative_to(ctx.out)),role='complete_TSV_only',
                evidence_scope=inventory[-1]['evidence_scope'],origin=relative,data_file=str(target.relative_to(ctx.out)),reason='large gene-level table supplied losslessly as TSV'))
    source_index['included_in_archive']=source_index.source.map(lambda v:'.tsv' in v or Path(v).suffix in ['.json','.npz','.py'])
    source_index['archive_member']=source_index.source.where(source_index.included_in_archive,'')
    ctx.table(source_index,'table/S29-source-data-index.tsv')
    inventory.append(dict(table=29,part='source_index',file='table/S29-source-data-index.tsv',source='45 verified statistical nodes, whole-genome figure scores and verified rendering source data',
        sha256=digest(ctx.out/'table/S29-source-data-index.tsv'),evidence_scope='mixed_explicitly_labeled',rows=len(source_index)))
    sheets.append(dict(table=29,name='S29',source=workbook_input(ctx,source_index,'S29'),role='complete_table',
        evidence_scope='mixed_explicitly_labeled',origin='table/S29-source-data-index.tsv',data_file='table/S29-source-data-index.tsv'))
    # Archive source tables and numerical/audit records, with historical scope retained in their path.
    with zipfile.ZipFile(ctx.out/'table/S29-statistical-source-data.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as archive:
        for relative,(path,expected) in artifacts.items():
            if '.tsv' in path.name or path.suffix in ['.json','.npz','.py']:
                archive.write(path,relative)
        archive.write(ctx.out/'table/S29-source-data-index.tsv','SOURCE-INDEX.tsv')
    if 'execution_provenance_archive' in ctx.config['resources']:
        source=ctx.resource('execution_provenance_archive')
        target=ctx.out/'table/S29-execution-provenance.zip';shutil.copy2(source,target)
        with zipfile.ZipFile(target) as archive:
            manifest=json.loads(archive.read('MANIFEST.json'))
            assert manifest['primary_run_id']==ctx.run_id
            assert manifest['primary_input_fingerprint']==ctx.config['input_fingerprint']
            import hashlib
            for path,sha in manifest['files'].items():assert hashlib.sha256(archive.read(path)).hexdigest()==sha
        inventory.append(dict(table=29,part='execution_provenance',file=str(target.relative_to(ctx.out)),source=str(source),
            sha256=digest(target),evidence_scope='analysis_code_and_configuration',rows=0))
    # Keep the final archive binding outside that archive to avoid a self-hash cycle.
    target=ctx.out/'table/S29-delivery-config.json';shutil.copy2(ctx.config_path,target)
    inventory.append(dict(table=29,part='delivery_config',file=str(target.relative_to(ctx.out)),source=str(ctx.config_path),
        sha256=digest(target),evidence_scope='analysis_code_and_configuration',rows=0))
    if 'reference_control_source_archive' in ctx.config['resources']:
        source=ctx.resource('reference_control_source_archive')
        target=ctx.out/'table/S18-v3-reference-control-source-data.zip';shutil.copy2(source,target)
        identities=ctx.read(ctx.analysis/'reference_comparison/reference_runs.tsv')
        control=identities[identities.reference.eq('v3')].iloc[0]
        with zipfile.ZipFile(target) as archive:
            manifest=json.loads(archive.read('MANIFEST.json'))
            assert manifest['run_id']==control.run_id and manifest['input_fingerprint']==control.input_fingerprint
            assert manifest['reference_variant']=='v3_ge62_comparison'
            assert set(manifest['completed_nodes'])==set(nodes+['figure_scores','cross_sharing_range10']+ctx.config.get('control_extra_source_nodes',[]))
            import hashlib
            for path,sha in manifest['files'].items():assert hashlib.sha256(archive.read(path)).hexdigest()==sha
        inventory.append(dict(table=18,part='full_reference_control',file=str(target.relative_to(ctx.out)),source=str(source),
            sha256=digest(target),evidence_scope='v3_reference_sensitivity_with_labeled_historical_audits',rows=0))
    ctx.table(inventory,'table/S1-S30-table-index.tsv')
    assert set(r['table'] for r in inventory)==set(range(1,31))
    payload=dict(run_id=ctx.run_id,input_fingerprint=ctx.config['input_fingerprint'],sheets=sheets,
        full_precision='per-table TSV and S29 source archive; XLSX displays IEEE-754 values up to 15 significant decimal digits',
        historical_scope='S18 and S29 include explicitly historical original-pipeline audit, never native polyApipe evidence',
        workbook_status='awaiting local artifact-tool export and visual verification')
    (ctx.out/'workbook-inputs.json').write_text(json.dumps(payload,indent=2)+'\n')
    # Completion is intentionally issued only after workbook export/verification adds its files.
    (ctx.out/'collection.complete.json').write_text(json.dumps(dict(run_id=ctx.run_id,input_fingerprint=ctx.config['input_fingerprint'],
        inputs=ctx.inputs,table_ids=list(range(1,31)),artifacts=len(artifacts),workbook_complete=False),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();collect(Context(a.config,'tables',a.out))
