"""Freeze the native caller comparison using existing integer matrices."""
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
    prefix='notes/2026-10-02-scape-polyapipe-comparison/'
    for name in ['atlas_effect_math.py','atlas_effect_engine.py']:
        relative=prefix+'scripts/'+name;source=PROJECT/relative;name=digest(source)[:12]+'_'+name
        shutil.copy2(source,out/'legacy'/name);config['source_overrides'][relative]=dict(path=str(dest/'legacy'/name),sha256=digest(source))
    inputs=PROJECT/prefix/'atlas/comparison_inputs'
    mapping=dict(caller_original_manifest=inputs/'preparation.complete.json',caller_source_hashes=inputs/'verified_source_sha256.json',
        caller_frozen_pairs=inputs/'frozen_pair_effects.tsv.gz',caller_frozen_routes=inputs/'frozen_route_effects.tsv.gz',
        caller_legacy_recurrence=PROJECT/'notes/2026-09-30-multimetric-statistics/nulls/study_event_covariates.tsv.gz',
        caller_legacy_members=PROJECT/'notes/2026-09-30-multimetric-statistics/nulls/recurrence_members.tsv.gz')
    for key,source in mapping.items():
        name=key+'__'+source.name;shutil.copy2(source,out/'resources'/name)
        config['resources'][key]=dict(path=str(dest/'resources'/name),sha256=digest(source))
    config['caller_comparison_root']='/root/nfsdata/apadetect/analysis_projects/reproapa_scape_polyapipe_20261002_v1/atlas'
    config['jobs']={'caller_comparison':dict(script=str(dest/'code/downstream_caller_comparison.py'),arguments=[],depends_on=['state_vectors'])}
    config['upstream_dispatch_status']=[str(Path(config['analysis_root'])/'dispatch-v1-recovered-status.json')]
    config['code_sha256']={str(dest/'code'/p.name):digest(p) for p in (out/'code').glob('*.py')}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'package.json').write_text(json.dumps(dict(scope='new native polyApipe against historical counts; SCAPE availability explicit',
        files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}),indent=2)+'\n')
    print(json.dumps(dict(package=str(out),nodes=list(config['jobs']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
