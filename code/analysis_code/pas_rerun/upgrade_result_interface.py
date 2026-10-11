"""Add the new result-record interface to a completed v8 bundle in a fresh folder.

Only verified, already computed tables are copied. This never starts a statistical
node, edits a running release, or changes the original bundle.
"""
import argparse
import csv
import gzip
from pathlib import Path
import shutil
from pipeline import digest,read,write,require
from export_statistics_bundle import panel_contracts,require_numerical_products
from result_records import contained,export_records,source_identity


def upgrade(bundle,contracts,out):
    bundle=Path(bundle).resolve();contracts=Path(contracts).resolve();out=Path(out).resolve()
    require(not out.exists(),'Use a new interface package directory')
    require(not out.is_relative_to(bundle),'Keep the source bundle immutable')
    require(not any(p in out.parts for p in ['Fig','main','supplement','production-v15']),'Never write into production')
    marker=read(bundle/'bundle-complete.json');require(marker['status']=='statistics_bundle_complete','Incomplete source bundle')
    for name,key in [('plot-contract.json','plot_contract_sha256'),('table-manifest.json','table_manifest_sha256')]:
        require(digest(bundle/name)==marker[key],'Unverified source bundle: '+name)
    original_contract=read(bundle/'plot-contract.json');identity=source_identity(original_contract)
    require(marker['run_id']==identity['run_id'],'Source run mismatch')
    manifest=read(bundle/'table-manifest.json');tables=manifest['tables']
    require(len({r['path'] for r in tables})==len(tables),'Duplicate source table')
    for r in tables:
        path=contained(bundle,r['path']);require(digest(path)==r['sha256'],'Changed source table: '+r['path'])
        if r['columns']:
            with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:
                require(next(csv.reader(f,delimiter='\t'),[])==r['columns'],'Source table header differs from manifest')
    required=require_numerical_products(contracts,tables)
    panels,_=panel_contracts(contracts,manifest['node_receipts'])
    for p in panels:
        p['tables']=[r['path'] for r in tables if r['columns'] and (r['node'] in p['nodes'] or r['node']=='frozen')]
        p['required_tables']=['tables/'+r['path'] for r in required if r['path'].split('/')[0] in p['nodes']]
    out.mkdir(parents=True)
    for r in tables:
        src=contained(bundle,r['path']);dest=contained(out,r['path']);dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dest);require(digest(dest)==r['sha256'],'Result copy differs')
    write(out/'table-manifest.json',manifest)
    write(out/'plot-contract.json',dict(original_contract,panels=panels))
    shutil.copytree(contracts,out/'manuscript-contracts')
    record=export_records(out,identity,panels,tables,contracts)
    parent=dict(source_bundle=str(bundle),source_marker_sha256=digest(bundle/'bundle-complete.json'),
                source_table_manifest_sha256=digest(bundle/'table-manifest.json'),
                source_plot_contract_sha256=digest(bundle/'plot-contract.json'),
                source_statistics_recomputed=False,production_replaced=False)
    write(out/'interface-parent.json',parent)
    write(out/'bundle-complete.json',dict(marker,
        plot_contract_sha256=digest(out/'plot-contract.json'),table_manifest_sha256=digest(out/'table-manifest.json'),
        result_records_sha256=digest(out/'result-records.json'),statistics_schema_id=record['statistics_schema_id'],
        required_statistical_products=len(required),interface_parent_sha256=digest(out/'interface-parent.json'),
        interface_only_upgrade=True,new_figures_rendered=False,production_replaced=False))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',required=True)
    p.add_argument('--contracts',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    print(upgrade(a.bundle,a.contracts,a.out))
