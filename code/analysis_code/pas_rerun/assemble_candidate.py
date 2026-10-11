"""Stage a separate manuscript only after all analysis and review receipts verify."""
import argparse
import csv
import json
from pathlib import Path
import re

from intake import HERE,digest
from run import load_workflow,verify_result
from stage import stage


def assemble(config,package,destination):
    config=Path(config).resolve();c=json.loads(config.read_text());root=Path(c['analysis_root'])
    assert not Path(destination).exists()
    nodes=[n['id'] for n in load_workflow(HERE/'workflow.json')['nodes'] if n['id']!='release']
    receipts={n:verify_result(root/n/'result.json',n,c['run_id'],c['input_fingerprint']) for n in nodes}
    for n in ['figure_scores']+c.get('extra_source_nodes',[]):
        verify_result(root/n/'result.json',n,c['run_id'],c['input_fingerprint'])
    files=[]
    def add(path,target,kind,**extra):
        files.append(dict(path=str(path.resolve()),sha256=digest(path),target=target,kind=kind,**extra))
    for row in csv.DictReader((HERE/'figures.tsv').open(),delimiter='\t'):
        if row['action']=='rebuild':add(root/'render'/row['target'],row['target'],'figure',title=row['title'])
    for node in ['tables','text']:
        for p in sorted((root/node/'table').glob('*')):
            if not p.is_file():continue
            match=re.match(r'S(\d+)(?:_|-|\.)',p.name)
            ids=[int(match[1])] if match else list(range(8,31))
            ids=[i for i in ids if 8<=i<=30]
            if ids:add(p,'table/'+p.name,'table',table_ids=ids)
    for p in sorted((root/'text/split_sources').rglob('*.tex')):
        add(p,str(p.relative_to(root/'text/split_sources')),'split_source')
    add(root/'text/references.bib','references.bib','bibliography')
    review=root/'text/text-review.tsv'
    payload=dict(run_id=c['run_id'],input_fingerprint=c['input_fingerprint'],files=files,
        text_review=dict(path=str(review),sha256=digest(review)))
    package=Path(package).resolve();assert not package.exists();package.write_text(json.dumps(payload,indent=2)+'\n')
    state=dict(status='running',active_node='release',run_id=c['run_id'],input_fingerprint=c['input_fingerprint'],
        completed=nodes,receipt_sha256={n:digest(root/n/'result.json') for n in nodes},
        note='All nodes independently verified; separate candidate still requires compilation and final page review.')
    (root/'run.json').write_text(json.dumps(state,indent=2)+'\n')
    return stage(package,root,destination)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--package',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(assemble(a.config,a.package,a.out))
