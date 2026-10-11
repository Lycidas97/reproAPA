"""Freeze catalogue QC, example tracing and independent external feature audits."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import PROJECT,REMOTE
from downstream_external_structure import TRACKS


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);dest=REMOTE/out.name
    for name in ['code','legacy','resources']:(out/name).mkdir()
    config=json.loads(args.base.read_text())
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    for relative in ['notes/2026-10-02-GSE282319-long-read/analyze.py','notes/2026-10-02-rad23b-ago2-gse25310/check_overlap.py']:
        path=PROJECT/relative;name=digest(path)[:12]+'_'+path.name;shutil.copy2(path,out/'legacy'/name)
        config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(path))
    lr=PROJECT/'notes/2026-10-02-GSE282319-long-read';ago=PROJECT/'notes/2026-10-02-rad23b-ago2-gse25310/inputs'
    mapping=dict(legacy_cross_members=PROJECT/'notes/2026-09-26-endpoint-counting/inputs/figure4_frozen_members.tsv',
        legacy_cross_support=PROJECT/'notes/2026-09-26-endpoint-counting/inputs/figure4_published_support.tsv',
        longread_gtf=lr/'raw/GSE282319_Repro_filter_tst.gtf.gz',longread_classification=lr/'raw/GSE282319_Repro_unfiltered_tst_classification.txt.gz',
        longread_tpm=lr/'raw/GSE282319_Isoquant.transcript_name_id.tpm.txt.gz',longread_sample_inventory=lr/'sample_inventory_with_passage.tsv',
        longread_published_classification=lr/'filtered_supplement_classification.tsv',longread_raw_sample_metadata=lr/'raw/GSE282319_family.soft.gz',
        longread_raw_supplement1=lr/'raw/supplement-data1.xlsx',longread_raw_supplement2=lr/'raw/supplement-data2.xlsx')
    for source,destbuild in [('mm39','mm10'),('mm10','mm39'),('mm9','mm10'),('mm10','mm9')]:
        mapping['chain_'+source+'_'+destbuild]=ago/(source+'To'+destbuild[0].upper()+destbuild[1:]+'.over.chain.gz')
    for prefix,_,_ in TRACKS:
        for ext in ['json','bed']:mapping['ago2_'+prefix+'_'+ext]=ago/(prefix+'.'+ext)
    for key,path in mapping.items():
        name=key+'__'+path.name;shutil.copy2(path,out/'resources'/name)
        config['resources'][key]=dict(path=str(dest/'resources'/name),sha256=digest(path))
    config['resources'].update(json.loads(args.genomic_resources.read_text()))
    config['jobs']={node:dict(script=str(dest/'code'/script),arguments=['--node',node] if node in ['global','tandem','longread','ago2'] else [],
        depends_on=deps,requires_core=node!='pas_qc') for node,script,deps in [
        ('global','downstream_global.py',['state_vectors']),('tandem','downstream_global.py',['global']),
        ('pas_qc','downstream_pas_qc.py',[]),('cross_trace','downstream_trace.py',['cross_overlap','global']),
        ('examples','downstream_examples.py',['matched_pms','cross_trace','site_statistics','well_inference','union_pool']),
        ('longread','downstream_external_structure.py',['examples']),('ago2','downstream_external_structure.py',['examples'])]}
    config['upstream_dispatch_status']=[str(Path(config['analysis_root'])/('dispatch-'+v+'-status.json')) for v in ['v1-recovered','v2','v3']]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='native-count analyses; legacy six-rule recount paused and not a dependency',nodes=list(config['jobs']),
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),remote_destination=str(dest),nodes=len(config['jobs']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--genomic-resources',type=Path,required=True);main(p.parse_args())
