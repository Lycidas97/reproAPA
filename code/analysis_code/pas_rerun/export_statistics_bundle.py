"""Freeze verified result tables for independent plotting and manuscript lookup."""
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
from downstream_runtime import Context
from pipeline import digest,read,write,require,verified_artifacts
from manuscript_contract import EVIDENCE,MAIN,SCOPES
from result_records import export_records


def read_contract(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))


def panel_contracts(root,available):
    main=read_contract(root/'main-panel-contract.tsv');supp=read_contract(root/'supplement-panel-contract.tsv')
    evidence=read_contract(root/'evidence-rerun-contract.tsv')
    require({r['evidence_id'] for r in evidence}==set(EVIDENCE),'Unbound manuscript evidence')
    require({f"F{r['figure']}{r['panel']}" for r in main}==set(MAIN),'Unbound main panel')
    panels=[]
    for r in main:
        key=f"F{r['figure']}{r['panel']}";panels.append(dict(panel=key,figure='F'+r['figure'],title=r['title'],nodes=MAIN[key],scope=r['action']))
    for r in supp:
        tags=re.findall(r'E\d\d[a-z]?',r['evidence'])
        ids={k for t in tags for k in EVIDENCE if k==t or (len(t)==3 and k.startswith(t))}
        require(bool(ids),f"No evidence binding for {r['new_panel']}")
        nodes=sorted({n for key in ids for n in EVIDENCE[key]})
        panels.append(dict(panel=r['new_panel'],figure=r['production_figure'],title=r.get('claims',''),
                           evidence_ids=sorted(ids),nodes=nodes,scope=r['action']))
    bindings=root/'supplement-producer-bindings.tsv'
    if bindings.exists():
        by_panel={p['panel']:p for p in panels}
        for r in read_contract(bindings):
            require(r['panel'] in by_panel,'Unknown supplementary binding: '+r['panel'])
            p=by_panel[r['panel']];p['nodes']=sorted(set(p['nodes'])|set(r['producers'].split(';')))
    for p in panels:require(set(p['nodes'])<=set(available),f"Missing producers: {p['panel']}")
    return panels,evidence


def require_numerical_products(contracts,rows):
    required=read_contract(contracts/'required-statistics-tables.tsv')
    extra=contracts/'supplement-statistics-tables.tsv'
    if extra.exists():required+=read_contract(extra)
    available={r['path'].removeprefix('tables/'):r for r in rows}
    for rule in required:
        path=rule['path'];require(path in available,f'Missing required statistical product: {path}')
        missing=set(rule['required_columns'].split(';'))-set(available[path]['columns'])
        require(not missing,f'Missing statistical columns in {path}: {sorted(missing)}')
    return required


def run(ctx):
    deps=ctx.config['jobs']['bundle']['depends_on'];contracts=Path(ctx.config['contracts_root'])
    panels,evidence=panel_contracts(contracts,deps);rows=[];receipts={}
    for node in deps:
        p=ctx.analysis/node/'result.json';record=read(ctx.register(p))
        require(record['node']==node and record['run_id']==ctx.run_id and record['input_fingerprint']==ctx.config['input_fingerprint'],f'Cross-run result: {node}')
        verified_artifacts(ctx.analysis/node,record);receipts[node]=digest(p)
        if node=='core':continue  # Consumers never receive a raw count-matrix interface.
        for a in record['artifacts']:
            if a['path'].endswith(('.png','.pdf','.svg')):raise ValueError(f'Graphics inside statistics node: {node}/{a["path"]}')
            source=ctx.analysis/node/a['path'];dest=ctx.out/'tables'/node/a['path']
            dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
            require(digest(dest)==a['sha256'],'Bundle copy differs')
            columns=[]
            if a['path'].endswith(('.tsv','.tsv.gz')):
                import gzip
                opener=gzip.open if a['path'].endswith('.gz') else open
                with opener(dest,'rt') as f:columns=next(csv.reader(f,delimiter='\t'),[])
            rows.append(dict(node=node,path=str(dest.relative_to(ctx.out)),sha256=a['sha256'],bytes=dest.stat().st_size,
                columns=columns,role=a.get('role','source_data'),historical_fields_present=any('legacy' in c or c.startswith('old_') for c in columns)))
    # Frozen expression coordinates/labels remain available to the presentation
    # layer without giving it access to an expression or PAS count matrix.
    frozen=['frozen_cell_embedding','frozen_metacell_states','expression_graph_nodes']
    frozen+=sorted(k for k in ctx.config['resources'] if k.startswith('frozen_figure_'))
    for key in frozen:
        source=ctx.resource(key);dest=ctx.out/'tables/frozen'/source.name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,dest)
        import gzip
        opener=gzip.open if dest.name.endswith('.gz') else open
        columns=[]
        if dest.name.endswith(('.tsv','.tsv.gz')):
            with opener(dest,'rt') as f:columns=next(csv.reader(f,delimiter='\t'),[])
        rows.append(dict(node='frozen',path=str(dest.relative_to(ctx.out)),sha256=digest(dest),bytes=dest.stat().st_size,
            columns=columns,role='frozen_expression_framework',historical_fields_present=False))
    required=require_numerical_products(contracts,rows)
    for p in panels:
        p['tables']=[r['path'] for r in rows if (r['node'] in p['nodes'] or r['node']=='frozen') and r['columns']]
        p['required_tables']=['tables/'+r['path'] for r in required if r['path'].split('/')[0] in p['nodes']]
        require(bool(p['tables']),f"Panel has no result tables: {p['panel']}")
    write(ctx.out/'plot-contract.json',dict(schema='pas-statistics-bundle/1',run_id=ctx.run_id,
        input_fingerprint=ctx.config['input_fingerprint'],caller=ctx.config['caller'],analysis_route=ctx.config['analysis_route'],
        reference_variant=ctx.config['reference_variant'],reference_policy='v2_only_v3_excluded_by_user',panels=panels,
        rules=['Read bundled tables only; do not read raw counts or old outputs.',
               'No fitting, testing, bootstrap, candidate selection or multiplicity correction in plotting.',
               'Use recorded denominators, uncertainty types and test-family labels.',
               'Missing estimates stay unavailable. Empty candidates never trigger an old-gene fallback.',
               'Render into a new candidate directory; promotion is a separate reviewed action.']))
    write(ctx.out/'table-manifest.json',dict(tables=rows,node_receipts=receipts))
    records=export_records(ctx.out,ctx.config,panels,rows,contracts)
    ctx.table([dict(evidence_id=r['evidence_id'],production_table=r['production_table'],title=r['title'],
        producers=';'.join(EVIDENCE[r['evidence_id']]),scope=SCOPES.get(r['evidence_id'],'Full current-run cohort and its documented statistical unit/correction family'),
        table_count=sum(x['node'] in EVIDENCE[r['evidence_id']] for x in rows)) for r in evidence],'evidence-index.tsv')
    shutil.copytree(contracts,ctx.out/'manuscript-contracts')
    write(ctx.out/'bundle-complete.json',dict(status='statistics_bundle_complete',run_id=ctx.run_id,
        source_nodes=len(deps),main_panels=len(MAIN),supplement_panels=len(panels)-len(MAIN),evidence_bundles=len(EVIDENCE),
        plot_contract_sha256=digest(ctx.out/'plot-contract.json'),table_manifest_sha256=digest(ctx.out/'table-manifest.json'),
        result_records_sha256=digest(ctx.out/'result-records.json'),statistics_schema_id=records['statistics_schema_id'],
        required_statistical_products=len(required),new_figures_rendered=False,production_replaced=False))
    ctx.finish(dict(all_node_receipts_verified=True,all_panels_bound=True,all_evidence_bound=True,
                    current_run_only=True,no_render_or_promotion=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'bundle',a.out))
