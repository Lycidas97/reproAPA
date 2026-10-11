"""Freeze corrected backgrounds and trajectory coverage before pending jobs start."""
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
    for name in ['analyze_figure2_dynamic_architecture_v1.py','analyze_figure2_cde_dynamics_v1.py']:
        relative='github/reproapa-analysis/modules/07_apa_dynamics/code/figure2/'+name
        source=PROJECT/relative;target=digest(source)[:12]+'_'+name;shutil.copy2(source,out/'legacy'/target)
        config['source_overrides'][relative]=dict(path=str(dest/'legacy'/target),sha256=digest(source))
    order=['endpoints','metrics','recurrence','dwui_pool','union_pool','profiles','well_inference','matched_pms',
        'pool_inference','cross_sharing','functional','moran','profile_gex','recurrence_controls','measurement','profile_controls','count_inference']
    assert set(order)==set(config['jobs'])-{'state_vectors'}
    config['jobs']={node:config['jobs'][node] for node in order}
    for job in config['jobs'].values():job['script']=str(dest/'code'/Path(job['script']).name)
    # Profiles read completed core library states, not the endpoint-only vector bundle.
    config['jobs']['profiles']['depends_on']=[]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='17 pending statistics; running state_vectors remains unchanged',nodes=order,
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),remote_destination=str(dest),nodes=len(order))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
