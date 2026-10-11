"""Freeze pending example/structure jobs with native inclusion influence tables."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import REMOTE,PROJECT


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);(out/'code').mkdir();(out/'resources').mkdir();dest=REMOTE/out.name
    config=json.loads(args.base.read_text());nodes=['tandem','cross_trace','examples','longread','ago2']
    config['jobs']={node:config['jobs'][node] for node in nodes}
    config['jobs']['examples']['depends_on'].append('intermediate')
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    for job in config['jobs'].values():job['script']=str(dest/'code'/Path(job['script']).name)
    source=PROJECT/'analysis_code/release_2026-09-24/sources/F5/site-recurrence-v1/supplement-compressed-v1/runs/annotation-pas-liu-E-r4-member-sensitivity/tables/candidates.tsv.gz'
    shutil.copy2(source,out/'resources/legacy_ten_loci.tsv.gz')
    config['resources']['legacy_ten_loci']=dict(path=str(dest/'resources/legacy_ten_loci.tsv.gz'),sha256=digest(source))
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(nodes=nodes,scope='New library-level interval inclusion, depth/time/leave-library sensitivity, fixed-source DWUI influence and all ten historical locus identities. Running global job preserved.',
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),destination=str(dest),nodes=nodes)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
