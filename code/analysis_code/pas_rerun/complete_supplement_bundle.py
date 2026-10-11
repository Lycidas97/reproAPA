"""Append supplementary statistics to a completed run in a new immutable bundle.

Used after an already running release completes; it neither changes its code nor
restarts its count/statistics jobs. New preparations include the node in their DAG.
"""
import argparse
import csv
import gzip
from pathlib import Path
import shutil

from downstream_runtime import Context
from downstream_supplement_statistics import run as supplement
from export_statistics_bundle import panel_contracts,require_numerical_products
from pipeline import digest,read,write,require,verified_artifacts
from result_records import ResultBundle,export_records,source_identity

SUPPLEMENT_CONTRACTS=['supplement-producer-bindings.tsv','supplement-statistics-tables.tsv','supplement-data-interface.tsv']


def complete(bundle,config,contracts,out):
    bundle=Path(bundle).resolve();config=Path(config).resolve();contracts=Path(contracts).resolve();out=Path(out).resolve()
    require(not out.exists(),'Use a new supplement completion directory')
    require(not out.is_relative_to(bundle),'Keep the source bundle immutable')
    require(not any(p in out.parts for p in ['Fig','main','supplement','production-v15']),'Never write into production')
    original=ResultBundle(bundle);settings=read(config)
    require(source_identity(settings)==original.identity,'Config belongs to a different PAS input')
    manifest=read(bundle/'table-manifest.json')
    require('supplement_statistics' not in manifest['node_receipts'],'Supplement node already delivered')
    if 'count_inference' in manifest['node_receipts']:
        require(read(bundle/'tables/count_inference/audit.json')['checks'].get('robust_probability_tie_comparison') is True,
            'Count inference requires the exact probability-tie repair before supplementary delivery')
    for row in manifest['tables']:
        from result_records import contained
        require(digest(contained(bundle,row['path']))==row['sha256'],'Source table changed: '+row['path'])
    for name in SUPPLEMENT_CONTRACTS:require((contracts/name).is_file(),'Missing supplement contract: '+name)
    out.mkdir(parents=True)
    new_contracts=out/'contracts';shutil.copytree(bundle/'manuscript-contracts',new_contracts)
    # Preserve the main-figure contract used by this source run (v9 has 81
    # components), even while another task edits the working-tree main interface.
    for name in SUPPLEMENT_CONTRACTS:shutil.copy2(contracts/name,new_contracts/name)
    node_config=out/'supplement.config.json'
    write(node_config,dict(settings,analysis_root=str(bundle/'tables'),contracts_root=str(new_contracts),
        supplement_parent_bundle=str(bundle),supplement_code_sha256={str(p):digest(p) for p in Path(__file__).parent.glob('*.py') if not p.name.startswith('test_')}))
    node=out/'supplement_statistics';supplement(Context(node_config,'supplement_statistics',node))
    result=read(node/'result.json');verified_artifacts(node,result)
    require(result['run_id']==original.identity['run_id'] and result['input_fingerprint']==original.identity['input_fingerprint'],'Cross-source supplement result')
    dest=out/'bundle';dest.mkdir();shutil.copytree(bundle/'tables',dest/'tables')
    tables=list(manifest['tables'])
    for row in tables:require(digest(dest/row['path'])==row['sha256'],'Copied source table changed: '+row['path'])
    for item in result['artifacts']:
        source=node/item['path'];target=dest/'tables/supplement_statistics'/item['path'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        require(digest(target)==item['sha256'],'Supplement copy changed')
        columns=[]
        if target.name.endswith(('.tsv','.tsv.gz')):
            with (gzip.open(target,'rt') if target.suffix=='.gz' else target.open()) as f:columns=next(csv.reader(f,delimiter='\t'),[])
        tables.append(dict(node='supplement_statistics',path=str(target.relative_to(dest)),sha256=item['sha256'],bytes=target.stat().st_size,
            columns=columns,role=item.get('role','source_data'),historical_fields_present=False))
    receipts=dict(manifest['node_receipts'],supplement_statistics=digest(node/'result.json'))
    required=require_numerical_products(new_contracts,tables)
    panels,_=panel_contracts(new_contracts,receipts)
    for p in panels:
        p['tables']=[r['path'] for r in tables if r['columns'] and (r['node'] in p['nodes'] or r['node']=='frozen')]
        p['required_tables']=['tables/'+r['path'] for r in required if r['path'].split('/')[0] in p['nodes']]
    write(dest/'table-manifest.json',dict(manifest,tables=tables,node_receipts=receipts))
    write(dest/'plot-contract.json',dict(read(bundle/'plot-contract.json'),panels=panels))
    shutil.copytree(new_contracts,dest/'manuscript-contracts')
    records=export_records(dest,original.identity,panels,tables,new_contracts)
    write(dest/'supplement-parent.json',dict(source_bundle=str(bundle),source_marker_sha256=digest(bundle/'bundle-complete.json'),
        source_config_sha256=digest(config),supplement_result_sha256=digest(node/'result.json'),
        original_statistics_reused=True,BAM_recount=False,production_replaced=False))
    marker=read(bundle/'bundle-complete.json')
    write(dest/'bundle-complete.json',dict(marker,source_nodes=len(receipts),
        plot_contract_sha256=digest(dest/'plot-contract.json'),table_manifest_sha256=digest(dest/'table-manifest.json'),
        result_records_sha256=digest(dest/'result-records.json'),statistics_schema_id=records['statistics_schema_id'],
        required_statistical_products=len(required),supplement_parent_sha256=digest(dest/'supplement-parent.json'),
        supplement_statistics_complete=True,new_figures_rendered=False,production_replaced=False))
    ResultBundle(dest,original.identity)
    return dest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['bundle','config','contracts','out']:p.add_argument('--'+name,required=True)
    a=p.parse_args();print(complete(a.bundle,a.config,a.contracts,a.out))
