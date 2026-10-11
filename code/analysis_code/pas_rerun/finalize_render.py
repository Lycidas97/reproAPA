"""Issue rendering completion only after all sources and figure-specific reviews pass."""
import argparse
import csv
import json
from pathlib import Path
import shutil

from downstream_runtime import Context
from downstream_tables import verified_sources
from polyapipe_prepare import digest


def finish(ctx,review):
    inventory=list(csv.DictReader(ctx.resource('figure_inventory').open(),delimiter='\t'))
    expected={r['target']:r for r in inventory if r['action']=='rebuild'}
    workflow=json.loads(ctx.resource('workflow').read_text())
    nodes=[r['id'] for r in workflow['nodes'] if r['id'] not in ['render','tables','text','release']]
    artifacts,_=verified_sources(ctx,nodes+['figure_scores']+ctx.config.get('extra_source_nodes',[]))
    reviews=json.loads(ctx.register(review).read_text())
    assert reviews['run_id']==ctx.run_id and reviews['input_fingerprint']==ctx.config['input_fingerprint']
    entries=reviews['figures'];assert len(entries)==len(expected)
    checks={r['target']:r for r in entries};assert set(checks)==set(expected)
    receipts=list(ctx.out.glob('Fig/*.sources.json'));assert len(receipts)==len(expected)
    render_config=json.loads(ctx.resource('render_config').read_text());protocol=[]
    for path in receipts:
        record=json.loads(path.read_text());target=record['target'];assert target in expected
        assert record['run_id']==ctx.run_id and record['input_fingerprint']==ctx.config['input_fingerprint']
        assert digest(ctx.out/target)==record['sha256']==checks[target]['sha256']
        if record['figure']=='Figure 1':
            retained=record['frozen_GEX_pixel_checks']
            assert len(retained)==2 and all(v['passed'] for v in retained)
        review_record=checks[target]
        assert review_record['layout']=='pass' and review_record['source_data']=='pass' and review_record['notes']
        assert (ctx.out/target).with_suffix('.png').is_file()
        for artifact in record['numerical_panel_exports']+record['render_code_snapshots']:
            assert digest(ctx.out/artifact['path'])==artifact['sha256']
        for source in record['source_files']:
            if source.startswith('resource:'):
                name=source.split(':',1)[1];r=render_config['resources'][name]
                copied=ctx.out/'resources'/name/Path(r['path']).name
                assert digest(copied)==r['sha256']
            else:assert source in artifacts
        record['visual_review_complete']=True;record['review']=review_record
        path.write_text(json.dumps(record,indent=2)+'\n')
        protocol.append(dict(figure=record['figure'],target=target,sha256=record['sha256'],panels=record['panels'],
            sources=record['source_files'],protocol=record['description'],review=review_record))
    source=ctx.config.get('figure_score_code_snapshot')
    if source:shutil.copytree(source,ctx.out/'resources/figure_scores_code',dirs_exist_ok=True)
    ctx.table([dict(figure=r['figure'],target=r['target'],panel=p,protocol=r['protocol'],source_files=json.dumps(r['sources']))
        for r in protocol for p in r['panels']],'panel-protocol-index.tsv')
    (ctx.out/'visual-review.json').write_text(json.dumps(reviews,indent=2)+'\n')
    (ctx.out/'render-progress.json').write_text(json.dumps(dict(run_id=ctx.run_id,input_fingerprint=ctx.config['input_fingerprint'],
        status='all_figures_and_visual_review_complete',figures=protocol),indent=2)+'\n')
    ctx.finish(dict(all_25_registered_PAS_figures=True,all_required_current_source_receipts_verified=True,
        lossless_plotted_arrays_archived=True,executed_render_code_snapshots=True,individual_visual_reviews=True,
        frozen_GEX_panels_retained=True,historical_audits_separately_labeled=True))
    # Rendering contains current results, retained GEX and explicitly historical
    # audit panels. Preserve those distinctions in the S29 archive inventory.
    receipt=ctx.out/'result.json';record=json.loads(receipt.read_text())
    for a in record['artifacts']:
        p=a['path']
        if 'historical_counting_steps/' in p:a['evidence_scope']='historical_original_pipeline'
        elif 'frozen_cell_embedding/' in p:a['evidence_scope']='frozen_GEX_metadata'
        elif p.startswith(('Fig/','panel-data/Figure_S4/','panel-data/Figure_S9/','resources/original_figure1/')):
            a['evidence_scope']='mixed_explicitly_labeled'
        elif p.startswith(('render-code/','resources/figure_scores_code','resources/figure_scores_code_audit')):
            a['evidence_scope']='analysis_code_and_provenance'
        else:a['evidence_scope']='current_native_rerun'
    receipt.write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--review',required=True)
    a=p.parse_args();finish(Context(a.config,'render',a.out),Path(a.review))
