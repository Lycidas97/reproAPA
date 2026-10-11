"""Compute fresh core effects in a new run-owned native directory."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

from pipeline import digest, read, write, require


def run(config_path, out):
    cfg=read(config_path); out=Path(out); native=Path(cfg['native_run'])
    require(native==out/'native','Core may only write its own native directory')
    native.mkdir()
    source=Path(cfg['prepared_native_source'])
    manifest=read(source/'intake.native.json')
    require(len(manifest['studies'])==7 and sum(r['libraries'] for r in manifest['studies'])==161,
            'Production core requires all seven studies, never a pilot')
    for name in ['intake.native.json','annotation.complete.json','frozen_native_membership_audit.json','frozen_native_membership_audit.tsv']:
        (native/name).symlink_to(source/name)
    folders=[]
    for record in manifest['studies']:
        name=record['study_id']; src=source/'studies'/name; dest=native/'studies'/name;dest.mkdir(parents=True)
        for p in src.iterdir():
            if p.is_file() and p.name not in ['effects.log'] and p.suffix in ['.gz','.tsv','.npz']:
                (dest/p.name).symlink_to(p.resolve())
        folders.append(dest)
    script=Path(__file__).with_name('polyapipe_effects.py')
    def one(folder):
        with (folder/'effects.log').open('w') as log:
            subprocess.run([sys.executable,str(script),'--study',str(folder)],stdout=log,stderr=subprocess.STDOUT,check=True)
        r=read(folder/'effects/complete.json')
        require(r['status']=='core_effects_complete_not_full_rerun','Incomplete core effects')
        return dict(study=folder.name,summaries=r['summaries'])
    with ThreadPoolExecutor(max_workers=cfg.get('core_workers',2)) as pool:
        studies=list(pool.map(one,folders))
    checks=dict(full_seven_study_scope=True,new_effects_no_old_cache=True,frozen_counts_read_only=True)
    write(out/'audit.json',dict(studies=studies,checks=checks,source=str(source)))
    artifacts=[dict(path=str(p.relative_to(out)),sha256=digest(p),role='audit' if p.name=='audit.json' else 'source_data')
               for p in sorted(out.rglob('*')) if p.is_file() and not p.is_symlink() and p.name not in ['command.log','result.json']]
    write(out/'result.json',dict(node='core',run_id=cfg['run_id'],input_fingerprint=cfg['input_fingerprint'],
                               status='complete',checks=checks,artifacts=artifacts))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.config,a.out)
