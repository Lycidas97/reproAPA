"""Accept a separate candidate only after compilation and complete page review."""
import argparse
import json
from pathlib import Path
import subprocess

from downstream_runtime import Context
from polyapipe_prepare import digest
from run import audit_baseline
from stage import check_package


def check_pages(candidate,review,register):
    import pymupdf
    assert set(review['documents'])=={'main','supplement'}
    for name,entry in review['documents'].items():
        pdf=register(candidate/'output/pdf'/f'{name}.pdf')
        assert digest(pdf)==entry['sha256'],'PDF changed after page review'
        with pymupdf.open(pdf) as document:count=len(document)
        pages=entry['pages']
        assert len(pages)==count and {p['page'] for p in pages}==set(range(1,count+1)),name
        for page in pages:
            assert page['layout']=='pass' and page['notes'],(name,page['page'])
            assert digest(register(Path(page['image'])))==page['image_sha256'],'Reviewed page image changed'


def finish(ctx,candidate,package,review):
    audit_baseline()
    payload=check_package(package,ctx.analysis)
    assert payload['run_id']==ctx.run_id and payload['input_fingerprint']==ctx.config['input_fingerprint']
    ctx.register(package)
    for row in payload['files']:
        assert digest(ctx.register(candidate/row['target']))==row['sha256'],'Candidate differs from its verified package'
    record=json.loads(ctx.register(review).read_text())
    assert record['run_id']==ctx.run_id and record['input_fingerprint']==ctx.config['input_fingerprint']
    for target in ['check','pdf']:
        with (ctx.out/f'make-{target}.log').open('w') as log:
            subprocess.run(['make',target],cwd=candidate,stdout=log,stderr=subprocess.STDOUT,check=True)
    check_pages(candidate,record,ctx.register)
    # Builds must not alter any split source or packaged asset.
    for row in payload['files']:assert digest(candidate/row['target'])==row['sha256']
    audit_baseline()
    (ctx.out/'page-review.json').write_text(json.dumps(record,indent=2)+'\n')
    ctx.finish(dict(package_unchanged=True,make_check_passed=True,make_pdf_passed=True,
        every_main_and_supplement_page_reviewed=True,original_manuscript_preserved=True))
    p=candidate/'PAS-RERUN.json';state=json.loads(p.read_text())
    assert state['run_id']==ctx.run_id and state['input_fingerprint']==ctx.config['input_fingerprint']
    state.update(status='complete_candidate_reviewed',release_receipt_sha256=digest(ctx.out/'result.json'))
    p.write_text(json.dumps(state,indent=2)+'\n')
    p=ctx.analysis/'run.json';state=json.loads(p.read_text());state['completed'].append('release')
    state.update(status='complete',active_node=None,note='Separate candidate compiled and all pages reviewed; original manuscript preserved.')
    state['receipt_sha256']['release']=digest(ctx.out/'result.json');p.write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['config','out','candidate','package','review']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();finish(Context(a.config,'release',a.out),a.candidate.resolve(),a.package.resolve(),a.review.resolve())
