"""Retain an explicitly historical audit without treating it as native evidence."""
import argparse
import json
import shutil

from downstream_runtime import Context
from polyapipe_prepare import digest


def run(ctx):
    scope=json.loads(ctx.resource('legacy_six_rule_scope').read_text())
    assert scope['user_instruction']=='保留为明确标注的旧流程审计，不作为本轮结果'
    registry=json.loads(ctx.resource('legacy_six_rule_files').read_text())
    rows=[]
    for item in registry[ctx.node]:
        source=ctx.register(item['path']);assert digest(source)==item['sha256']
        target=ctx.out/'historical'/item['relative'];target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        if source.name=='audit.json':assert json.loads(source.read_text())['status']=='PASS'
        rows.append(dict(file=str(target.relative_to(ctx.out)),sha256=digest(target),evidence_scope='historical_original_pipeline',
            native_polyapipe_evidence=False,caption_required='Historical audit of the original PAS pipeline; not a polyApipe rerun result'))
    assert rows
    ctx.table(rows,'historical_evidence_manifest.tsv')
    (ctx.out/'SCOPE.json').write_text(json.dumps(dict(user_instruction=scope['user_instruction'],
        scientific_scope='historical_original_pipeline',native_result=False,new_BAM_counting=False,
        use_restriction='Never use these effects, cohorts, p values or q values as native polyApipe evidence'),ensure_ascii=False,indent=2)+'\n')
    ctx.finish(dict(user_selected_historical_audit=True,archived_files_hash_verified=True,
        no_new_BAM_counting=True,separate_from_native_inference=True),result_scope='historical_audit',native_result=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--node',choices=['recount','consensus'],required=True)
    a=p.parse_args();run(Context(a.config,a.node,a.out))
