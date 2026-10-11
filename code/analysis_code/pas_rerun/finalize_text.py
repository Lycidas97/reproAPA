"""Register reviewed candidate prose only after every original claim is resolved."""
import argparse
import csv
import json
from pathlib import Path
import re

from downstream_runtime import Context
from downstream_tables import verified_sources
from intake import HERE
from polyapipe_prepare import digest
from run import audit_baseline


def check_review(rows,targets,split,baseline):
    expected={r['id']:r for r in targets}
    assert len(rows)==len(expected) and {r['id'] for r in rows}==set(expected),'Incomplete text review'
    for row in rows:
        original=expected[row['id']]
        assert row['baseline_sha256']==original['sha256'] and row['file']==original['file']
        assert row['status'] in ['rewritten','verified_unchanged','not_result_verified'] and row['evidence'].strip(),row['id']
        candidate=split/original['file'];source=candidate if candidate.exists() else baseline/original['file']
        paragraphs=[' '.join(m.group().strip().split()) for m in re.finditer(r'\S[\s\S]*?(?=\n\s*\n|\Z)',source.read_text())]
        if row['status']=='rewritten':
            assert candidate.is_file() and digest(candidate)!=digest(baseline/original['file']),row['id']
        else:assert original['text'] in paragraphs,f"Unchanged paragraph was edited: {row['id']}"


def finish(ctx,review):
    audit_baseline()
    targets=list(csv.DictReader(ctx.register(HERE/'text-targets.tsv').open(),delimiter='\t'))
    rows=list(csv.DictReader(ctx.register(ctx.out/'text-review.tsv').open(),delimiter='\t'))
    baseline=Path(ctx.config['manuscript_root']);split=ctx.out/'split_sources'
    check_review(rows,targets,split,baseline)
    workflow=json.loads(ctx.resource('workflow').read_text())
    nodes=[n['id'] for n in workflow['nodes'] if n['id'] not in ['text','release']]
    verified_sources(ctx,nodes+['figure_scores']+ctx.config.get('extra_source_nodes',[]))
    decision=json.loads(ctx.register(review).read_text())
    assert decision['run_id']==ctx.run_id and decision['input_fingerprint']==ctx.config['input_fingerprint']
    assert decision['text_review_sha256']==digest(ctx.out/'text-review.tsv')
    required=['claims_checked_against_current_results','all_captions_match_current_panels',
        'empty_and_negative_results_reported','historical_audits_distinguished','original_manuscript_preserved']
    assert all(decision['checks'].get(k) is True for k in required) and decision['notes']
    files={str(p.relative_to(split)):digest(p) for p in split.rglob('*.tex')}
    assert files==decision['split_sources'],'Reviewed prose changed after review'
    assert digest(ctx.out/'references.bib')==decision['bibliography_sha256']
    allowed={r['file'] for r in csv.DictReader((HERE/'baseline-files.tsv').open(),delimiter='\t')}
    assert set(files)<=allowed,'Structural additions require the section-management workflow'
    (ctx.out/'text-release-review.json').write_text(json.dumps(decision,indent=2)+'\n')
    ctx.finish(dict(all_249_text_targets_resolved=True,unchanged_paragraphs_verified=True,
        current_statistics_figures_and_tables_verified=True,all_final_split_sources_reviewed=True,
        original_manuscript_preserved=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--review',required=True)
    a=p.parse_args();finish(Context(a.config,'text',a.out),Path(a.review))
