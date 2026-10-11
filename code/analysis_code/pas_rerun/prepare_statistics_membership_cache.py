"""Freeze the verified metadata-index optimization for explicitly selected pending nodes."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import REMOTE


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);(out/'code').mkdir();dest=REMOTE/out.name
    config=json.loads(args.base.read_text());nodes=args.nodes.split(',')
    assert set(nodes)<=set(config['jobs'])
    config['jobs']={node:config['jobs'][node] for node in nodes}
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    for job in config['jobs'].values():job['script']=str(dest/'code'/Path(job['script']).name)
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(nodes=nodes,scope='same statistical rules; frozen metadata row indices cached',
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),destination=str(dest),nodes=nodes)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--nodes',required=True);main(p.parse_args())
