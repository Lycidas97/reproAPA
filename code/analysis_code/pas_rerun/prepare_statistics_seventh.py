"""Freeze independent PAS-seq, cis models, and user-selected historical audits."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import PROJECT,REMOTE


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);dest=REMOTE/out.name
    for name in ['code','legacy','resources','historical']:(out/name).mkdir()
    config=json.loads(args.base.read_text())
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    prefix='notes/2026-10-02-gse99922-pilot/'
    for relative in [prefix+'expanded_measurements.py',prefix+'matched_analysis.py',prefix+'vector_analysis.py',
        'github/reproapa-analysis/modules/09_cross_species/code/cis_supplement/scripts/recompute_models.py']:
        path=PROJECT/relative;name=digest(path)[:12]+'_'+path.name;shutil.copy2(path,out/'legacy'/name)
        config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(path))
    source=PROJECT/prefix;mapping=dict(external_pas_counts=source/'inputs/GSE99922_readcounts-table.txt.gz',
        external_mouse_gene_info=source/'inputs/Mus_musculus.gene_info.gz',external_vector_protocol=source/'vector-protocol.json',
        external_legacy_candidates=source/'inputs/frozen_candidates.tsv',external_sample_manifest=source/'inputs/sample_manifest.tsv',
        chain_mm9_mm10=source/'inputs/mm9ToMm10.over.chain.gz',chain_mm10_mm39=source/'inputs/mm10ToMm39.over.chain.gz',
        legacy_six_rule_scope=PROJECT/'notes/2026-10-03-polyapipe-rerun/legacy-six-rule-scope.json')
    for key,path in mapping.items():
        name=key+'__'+path.name;shutil.copy2(path,out/'resources'/name)
        config['resources'][key]=dict(path=str(dest/'resources'/name),sha256=digest(path))
    config['resources'].update(json.loads(args.cis_resources.read_text()))
    registry={'recount':[],'consensus':[]}
    for node,directories in [('recount',['fixed-catalog-results']),('consensus',['figure3-recount-recurrence','figure3-joint-consensus-null'])]:
        for directory in directories:
            source=PROJECT/'notes/2026-09-29-core-sensitivity'/directory
            for path in sorted(source.iterdir()):
                if not path.is_file() or not ('.tsv' in path.name or path.suffix in ['.json','.npz']):continue
                relative=Path(directory)/path.name;target=out/'historical'/relative;target.parent.mkdir(exist_ok=True);shutil.copy2(path,target)
                registry[node].append(dict(path=str(dest/'historical'/relative),relative=str(relative),sha256=digest(path)))
    path=out/'resources/legacy_six_rule_files.json';path.write_text(json.dumps(registry,indent=2)+'\n')
    config['resources']['legacy_six_rule_files']=dict(path=str(dest/'resources'/path.name),sha256=digest(path))
    config['jobs']={node:dict(script=str(dest/'code'/script),arguments=['--node',node] if node in ['recount','consensus'] else [],depends_on=deps) for node,script,deps in [
        ('recount','downstream_legacy_audit.py',[]),('consensus','downstream_legacy_audit.py',['recount']),
        ('external_pas','downstream_external_pas.py',['recurrence']),('cis','downstream_cis.py',['cross_continuous','cross_overlap'])]}
    config['upstream_dispatch_status']=[str(Path(config['analysis_root'])/('dispatch-'+v+'-status.json')) for v in ['v1-recovered','v2','v3']]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='native external and cis analyses; six-rule nodes explicitly historical per user decision',nodes=list(config['jobs']),
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),remote_destination=str(dest),nodes=len(config['jobs']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cis-resources',type=Path,required=True);main(p.parse_args())
