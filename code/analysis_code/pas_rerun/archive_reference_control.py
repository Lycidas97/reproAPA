"""Archive the completed control statistics separately from primary inference."""
import argparse
import json
from pathlib import Path
import zipfile

from polyapipe_prepare import digest
from run import verify_result


def archive(config_path,workflow_path,out):
    config=json.loads(config_path.read_text());workflow=json.loads(workflow_path.read_text())
    assert config['reference_variant']=='v3_ge62_comparison'
    nodes=[n['id'] for n in workflow['nodes'] if n['id'] not in ['render','tables','text','release']]
    nodes+=['figure_scores','cross_sharing_range10']+config.get('extra_control_nodes',[])
    assert len(nodes)==len(set(nodes))
    root=Path(config['analysis_root']);files={}
    for node in nodes:
        receipt=root/node/'result.json'
        result=verify_result(receipt,node,config['run_id'],config['input_fingerprint'])
        files[str(receipt.relative_to(root))]=digest(receipt)
        files.update({str(Path(node)/a['path']):a['sha256'] for a in result['artifacts']})
    assert not out.exists(),'Archive destination already exists'
    manifest=dict(run_id=config['run_id'],input_fingerprint=config['input_fingerprint'],
        reference_variant=config['reference_variant'],completed_nodes=nodes,files=files,
        workflow_sha256=digest(workflow_path),
        evidence_scope='v3 reference sensitivity; retained historical audits remain labeled in their receipts')
    temporary=out.with_suffix('.partial.zip')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as bundle:
        for name,sha in files.items():
            path=root/name;assert digest(path)==sha
            bundle.write(path,name)
        bundle.writestr('MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
        bundle.write(config_path,'EXECUTION-CONFIG.json')
    temporary.replace(out)
    print(json.dumps(dict(path=str(out),sha256=digest(out),nodes=len(nodes),files=len(files))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True)
    p.add_argument('--workflow',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();archive(a.config,a.workflow,a.out)
