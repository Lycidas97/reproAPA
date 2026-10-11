"""Freeze complete process-state sensitivities before the shared-state phase starts."""
import argparse
import json
from pathlib import Path
import shutil

from polyapipe_prepare import digest
from prepare_statistics_package import PROJECT,REMOTE


def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=False);dest=REMOTE/out.name
    for name in ['code','legacy','resources']:(out/name).mkdir()
    config=json.loads(args.base.read_text())
    for source in Path(__file__).parent.glob('*.py'):
        if source.name!='downstream_common.py':shutil.copy2(source,out/'code'/source.name)
    relative='analysis_code/recurrence_correction_2026-09-25/sources/MOUSE/analysis/figure3-recurrence-common-reference-2026-09-25/analyze.py'
    source=PROJECT/relative;name=digest(source)[:12]+'_process.py';shutil.copy2(source,out/'legacy'/name)
    config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(source))
    source=PROJECT/'notes/2026-09-26-state-comparability-audit/metacell_expression_support.tsv.gz';shutil.copy2(source,out/'resources'/source.name)
    config['resources']['state_expression_support']=dict(path=str(dest/'resources'/source.name),sha256=digest(source))
    config['jobs'].pop('metacell')
    for job in config['jobs'].values():job['script']=str(dest/'code'/Path(job['script']).name)
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='remaining shared-state, APA/GEX, continuous jobs; metacell receipt preserved',
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),nodes=list(config['jobs']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
