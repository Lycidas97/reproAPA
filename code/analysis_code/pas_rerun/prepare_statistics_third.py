"""Freeze shared-state and continuous-comparison jobs as an independent phase."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import PROJECT,REMOTE


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);dest=REMOTE/out.name
    for name in ['code','legacy']:(out/name).mkdir()
    config=json.loads(args.base.read_text())
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    prefix='github/reproapa-analysis/modules/09_cross_species/code/figure4/'
    sources=[prefix+'recurrence_predicts_cross_species_v26.py',prefix+'audit_figure4c_reused_controls_v53.py',
        'github/reproapa-analysis/modules/08_recurrence/code/apa_gex_state/analyze.py']
    for relative in sources:
        path=PROJECT/relative;name=digest(path)[:12]+'_'+path.name;shutil.copy2(path,out/'legacy'/name)
        config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(path))
    config['jobs']={node:dict(script=str(dest/'code'/script),arguments=[],depends_on=deps) for node,script,deps in [
        ('metacell','downstream_metacell.py',[]),('intermediate','downstream_intermediate.py',['recurrence']),
        ('apa_gex','downstream_apa_gex.py',['intermediate','recurrence']),('cross_continuous','downstream_continuous.py',['cross_overlap'])]}
    config['upstream_dispatch_status']=[str(Path(config['analysis_root'])/n) for n in ['dispatch-v1-recovered-status.json','dispatch-v2-status.json']]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='four additional adapters; full manuscript rerun still incomplete',nodes=list(config['jobs']),
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),remote_destination=str(dest),nodes=len(config['jobs']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
