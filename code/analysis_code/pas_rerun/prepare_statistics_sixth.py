"""Freeze coordinate and intermediate-state audits against native count outputs."""
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
    for relative in ['notes/2026-09-29-core-sensitivity/figure5_coordinate_audit.py',
        'analysis_code/release_2026-09-24/sources/F5/site-recurrence-v1/supplement-compressed-v1/prepare.py']:
        path=PROJECT/relative;name=digest(path)[:12]+'_'+path.name;shutil.copy2(path,out/'legacy'/name)
        config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(path))
    config['jobs']={'coordinate_states':dict(script=str(dest/'code/downstream_coordinate_states.py'),arguments=[],
        depends_on=['intervals','site_scan','dwui_pool','original_pool_sites'])}
    config['upstream_dispatch_status']=[str(Path(config['analysis_root'])/('dispatch-'+v+'-status.json')) for v in ['v1-recovered','v2']]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='native matrix coordinate and R/E audits; no BAM counting',nodes=list(config['jobs']),
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),remote_destination=str(dest),nodes=len(config['jobs']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
