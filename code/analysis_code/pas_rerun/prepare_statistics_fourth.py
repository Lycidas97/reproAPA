"""Freeze native-window recount and its direction/consensus analyses."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import PROJECT,PROJECTS,REMOTE


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);dest=REMOTE/out.name
    for name in ['code','legacy','resources']:(out/name).mkdir()
    config=json.loads(args.base.read_text())
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    for relative in ['notes/2026-09-26-endpoint-counting/common.py','notes/2026-09-26-global-direction-audit/analyze.py',
                     'notes/2026-09-26-global-direction-audit/tandem_utr.py']:
        path=PROJECT/relative;name=digest(path)[:12]+'_'+path.name;shutil.copy2(path,out/'legacy'/name)
        config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(path))
    mapping=dict(native_libraries=PROJECT/'notes/2026-10-02-scape-polyapipe-comparison/atlas/input_libraries.tsv',
                 legacy_endpoint_ledger=PROJECT/'notes/2026-09-26-endpoint-counting/inputs/ledger.json')
    for sp,species,build in [('human','Homo_sapiens','GRCh38'),('mouse','Mus_musculus','GRCm39')]:
        mapping[sp+'_ensembl_gtf']=PROJECTS/'2026-09-13_human-mouse-ortholog-mapping/sources/web/ensembl116/genes'/species.lower()/(species+'.'+build+'.116.gtf.gz')
    for key,path in mapping.items():
        name=key+'__'+path.name;shutil.copy2(path,out/'resources'/name)
        config['resources'][key]=dict(path=str(dest/'resources'/name),sha256=digest(path))
    config['native_caller_root']='/root/nfsdata/apadetect/analysis_projects/reproapa_scape_polyapipe_20261002_v1/atlas'
    config['recount_workers']=2
    config['jobs']={node:dict(script=str(dest/'code'/script),arguments=['--node',node] if node in ['global','tandem'] else [],
        depends_on=deps,requires_core=core) for node,script,deps,core in [
        ('recount','downstream_recount.py',[],False),('consensus','downstream_consensus.py',['recount'],False),
        ('global','downstream_global.py',['state_vectors','recount'],True),('tandem','downstream_global.py',['global'],True)]}
    config['upstream_dispatch_status']=[str(Path(config['analysis_root'])/'dispatch-v1-recovered-status.json')]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='four additional adapters; full manuscript rerun incomplete',nodes=list(config['jobs']),
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),remote_destination=str(dest),nodes=len(config['jobs']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
