"""Run independent study comparisons concurrently, then verify and collect all seven."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

from downstream_runtime import Context,STUDIES
from downstream_caller_comparison import migration
from polyapipe_prepare import digest

TABLES=['route_effects_with_full_missing.tsv.gz','effect_agreement_summary.tsv',
    'bidirectional_PAS_distances.tsv','complete_gene_background.tsv.gz',
    'caller_availability.tsv','integer_count_and_membership_audit.tsv']


def run(ctx):
    def study(sid):
        out=ctx.out/'shards'/sid;out.mkdir(parents=True,exist_ok=True)
        receipt=out/'result.json'
        if not receipt.exists():
            with (out/'command.log').open('w') as log:
                subprocess.run([sys.executable,'-u',str(Path(__file__).with_name('downstream_caller_comparison.py')),
                    '--config',str(ctx.config_path),'--out',str(out),'--study',sid],stdout=log,stderr=subprocess.STDOUT,check=True)
        return sid,receipt
    with ThreadPoolExecutor(max_workers=7) as pool:
        completed=list(pool.map(study,STUDIES))
    frames={name:[] for name in TABLES}
    for sid,path in completed:
        record=json.loads(ctx.register(path).read_text())
        assert record['node']=='caller_comparison_'+sid and record['status']=='complete' and all(record['checks'].values())
        assert record['run_id']==ctx.run_id and record['input_fingerprint']==ctx.config['input_fingerprint']
        artifacts={a['path']:a['sha256'] for a in record['artifacts']}
        for relative,expected in artifacts.items():assert digest(ctx.register(path.parent/relative))==expected
        for name in TABLES:
            assert name in artifacts
            f=ctx.read(path.parent/name);assert set(f.study_id)<={sid},(name,sid)
            if not len(f):assert name in TABLES[:2],(name,sid)
            frames[name].append(f)
    for name,parts in frames.items():ctx.table(pd.concat(parts,ignore_index=True),name)
    routes=pd.concat(frames[TABLES[0]],ignore_index=True)
    assert not routes.duplicated(['study_id','route','gene_id','contrast','design','view','contributors']).any()
    migration(ctx,routes)
    ctx.finish(dict(all_seven_study_shards_verified=True,no_BAM_recount=True,
        unchanged_per_study_comparison_function=True,full_missing_denominators_retained=True,
        unchanged_frozen_cells_pairs_and_weights=True,complete_legacy_roster_migration=True),
        comparison_scope='native polyApipe versus frozen original pipeline; SCAPE remains a separate reserved comparison',
        execution='Seven independent studies; identical per-study numerical functions and seeds; concatenation in frozen study order')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'caller_comparison',a.out))
