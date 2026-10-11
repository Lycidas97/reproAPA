"""One result-record interface for every PAS caller and processing source.

Records address immutable result tables, not raw counts. A processing source is
identified outside those tables: their biological `route` column remains intact.
"""
import csv
import gzip
import hashlib
import json
from pathlib import Path
from pipeline import digest,read,write,require

IDENTITY=['run_id','input_fingerprint','caller','analysis_route','reference_variant']
SCHEMA='pas-panel-result-records/1'


def content_id(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def source_identity(value):
    require(all(isinstance(value.get(k),str) and value[k] for k in IDENTITY),'Incomplete result source identity')
    require(value['reference_variant']=='v2_primary','Only v2_primary is in the current analysis contract')
    identity={k:value[k] for k in IDENTITY}
    return dict(identity,source_key=content_id(identity))


def contained(root,name):
    root=Path(root).resolve();relative=Path(name)
    require(not relative.is_absolute(),'Result paths must be relative')
    path=(root/relative).resolve()
    require(path.is_relative_to(root),'Result path escapes its bundle')
    return path


def export_records(out,config,panels,tables,contracts):
    """Materialize stable component IDs without copying or flattening large tables."""
    out=Path(out);contracts=Path(contracts)
    with (contracts/'panel-data-interface.tsv').open() as f:interface=list(csv.DictReader(f,delimiter='\t'))
    extra=contracts/'supplement-data-interface.tsv'
    if extra.exists():
        with extra.open() as f:interface+=list(csv.DictReader(f,delimiter='\t'))
    registry=read(contracts/'figure-renderers.json')
    identity=source_identity(config);by_path={r['path']:r for r in tables}
    allowed={p['panel']:set(p['tables']) for p in panels};components=[];seen=set()
    for item in interface:
        key=item['component_id'];panel=item['panel'];name='tables/'+item['table']
        require(key not in seen,'Duplicate component ID: '+key);seen.add(key)
        require(name in allowed.get(panel,set()),f'Unbound component {panel}: {name}')
        require(name in by_path,'Missing component table: '+name)
        meta=by_path[name];columns=item['read_columns'].split(';')
        require(set(columns)<=set(meta['columns']),'Component schema mismatch: '+key)
        components.append(dict(component_id=key,panel=panel,figure=panel[:-1],table=name,
            table_sha256=meta['sha256'],columns=columns,producer=meta['node'],
            selection=item['selection'],unit_or_uncertainty=item['unit_or_uncertainty'],
            source_key=identity['source_key'],record_id=identity['source_key']+':'+key))
    require(set(allowed)>=set(r['panel'] for r in components),'Unknown component panel')
    schema_columns={r['table']:r['columns'] for r in components}
    # A schema fingerprint must not depend on caller, values, table hashes or IDs.
    schema_id=content_id(dict(interface=interface,columns=schema_columns))
    record=dict(schema=SCHEMA,statistics_schema_id=schema_id,source=identity,
        plot_contract_sha256=digest(out/'plot-contract.json'),table_manifest_sha256=digest(out/'table-manifest.json'),
        components=components,figures=registry['figures'],
        tables=[dict(table=r['path'],table_sha256=r['sha256'],columns=r['columns'],producer=r['node'])
                for r in tables if r['columns']],
        rules=dict(biological_route_is_not_analysis_route=True,statistics_recomputed_by_renderer=False,
                   unavailable_is_not_zero=True,old_source_fallback=False,production_replacement=False))
    write(out/'result-records.json',record)
    with (out/'result-records.tsv').open('w') as f:
        names=['source_key',*IDENTITY,'statistics_schema_id','record_id','component_id','panel','figure','table',
               'table_sha256','producer','columns','selection','unit_or_uncertainty']
        writer=csv.DictWriter(f,fieldnames=names,delimiter='\t');writer.writeheader()
        for r in components:
            row={k:r[k] for k in names if k in r};row.update({k:identity[k] for k in ['source_key',*IDENTITY]})
            row.update(statistics_schema_id=schema_id,columns=';'.join(r['columns']));writer.writerow(row)
    return record


class ResultBundle:
    """The same table reader is used by every figure, independent of PAS input."""
    def __init__(self,bundle,expected=None):
        self.root=Path(bundle).resolve();marker=read(self.root/'bundle-complete.json')
        require(marker['status']=='statistics_bundle_complete','Result bundle is incomplete')
        for name,key in [('plot-contract.json','plot_contract_sha256'),('table-manifest.json','table_manifest_sha256'),
                         ('result-records.json','result_records_sha256')]:
            require(digest(self.root/name)==marker.get(key),'Unverified result interface: '+name)
        self.record=read(self.root/'result-records.json');require(self.record['schema']==SCHEMA,'Unsupported result interface')
        self.identity=source_identity(self.record['source'])
        require(self.identity==self.record['source'],'Source key mismatch')
        require(marker['run_id']==self.identity['run_id'],'Cross-run result marker')
        contract=read(self.root/'plot-contract.json')
        require(source_identity(contract)==self.identity,'Cross-source plot contract')
        for key in ['plot_contract_sha256','table_manifest_sha256']:
            require(self.record[key]==marker[key],'Result index and bundle disagree: '+key)
        if expected is not None:
            require(source_identity(expected)==self.identity,'Requested source differs from result bundle')
        self.components={r['component_id']:r for r in self.record['components']}
        require(len(self.components)==len(self.record['components']),'Duplicate component ID')
        self.manifest={r['path']:r for r in read(self.root/'table-manifest.json')['tables']}
        self.allowed={p['panel']:set(p['tables']) for p in contract['panels']}
        self.used={}

    def table(self,panel,component_id):
        require(component_id in self.components,'Unknown component: '+component_id)
        r=self.components[component_id]
        require(r['panel']==panel,'Component belongs to another panel')
        require(r['source_key']==self.identity['source_key'],'Component belongs to another PAS source')
        require(r['record_id']==self.identity['source_key']+':'+component_id,'Invalid component record identity')
        require(r['table'] in self.allowed.get(panel,set()),'Panel is not bound to this table')
        meta=self.manifest.get(r['table'])
        require(meta is not None and meta['sha256']==r['table_sha256'],'Component table fingerprint mismatch')
        path=contained(self.root,r['table'])
        require(digest(path)==r['table_sha256'],'Result table changed: '+r['table'])
        with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:
            columns=next(csv.reader(f,delimiter='\t'),[])
        require(set(r['columns'])<=set(columns),'Missing result columns: '+component_id)
        self.used[r['table']]=r['table_sha256']
        return path

    def provenance(self):
        return dict(source=self.identity,result_records_sha256=digest(self.root/'result-records.json'),
                    statistics_schema_id=self.record['statistics_schema_id'],tables=dict(self.used),
                    statistics_recomputed=False,production_replaced=False)


def collect_sources(items,out):
    """Collect completed bundles, including caller/processing sources from different runs."""
    sources=[];names=set();keys=set()
    import re
    for item in items:
        name=item['source_id'];require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',name) is not None,'Unsafe source_id')
        require(name not in names,'Duplicate source_id');names.add(name)
        bundle=ResultBundle(item['bundle'],item.get('identity'));key=bundle.identity['source_key']
        require(key not in keys,'The same source bundle is listed twice');keys.add(key)
        sources.append(dict(source_id=name,bundle=str(bundle.root),identity=bundle.identity,
            statistics_schema_id=bundle.record['statistics_schema_id'],
            result_records_sha256=digest(bundle.root/'result-records.json')))
    require(bool(sources),'No completed sources')
    require(len({s['statistics_schema_id'] for s in sources})==1,'Sources use different statistical interfaces')
    result=dict(schema='pas-result-series/1',status='complete_result_index',sources=sources,
        production_replaced=False,figures_rendered=False)
    require(not Path(out).exists(),'Use a new result-series path')
    write(out,result);return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description='Collect completed PAS result bundles without recalculating statistics')
    p.add_argument('--bundle',action='append',required=True,metavar='SOURCE_ID=PATH');p.add_argument('--out',required=True)
    a=p.parse_args();items=[]
    for value in a.bundle:
        name,sep,path=value.partition('=');require(bool(sep and path),'Expected SOURCE_ID=PATH')
        items.append(dict(source_id=name,bundle=path))
    collect_sources(items,a.out)
