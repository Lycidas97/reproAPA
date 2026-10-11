"""Launch both pinned reference runs after reference conversion has completed."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from polyapipe_prepare import digest


def main(args):
    root=args.root.resolve();package=args.package.resolve() if args.package else Path(__file__).resolve().parent.parent;code=package/'code'
    plan=json.loads((package/'plan.json').read_text())
    for relative,expected in plan['package_files'].items():
        assert digest(package/relative)==expected,relative
    assert not (root/'statistics-launch.json').exists(),'A statistics launch already exists'
    def wait_for(marker):
        while not marker.exists():
            # A zombie or missing conversion worker is failure, never an infinite wait.
            state=Path('/proc')/str(args.reference_pid)/'stat'
            assert state.exists() and state.read_text().split()[2]!='Z','Reference conversion failed; inspect reference-replacement.log'
            time.sleep(15)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    launches=[]
    for variant in ['v2_primary','v3_ge62_comparison']:
        native=root/'runs'/variant
        wait_for(native/'reference_audit.json')
        audit=json.loads((native/'reference_audit.json').read_text())
        for path,expected in audit['inputs'].items():assert digest(path)==expected,path
        annotation=json.loads((native/'annotation.complete.json').read_text())
        assert annotation['status']=='complete' and len(annotation['studies'])==7
        for study in annotation['studies']:
            assert digest(native/'studies'/study['study_id']/'features.tsv.gz')==study['features_sha256']
        config=dict(plan['config'])
        config.update(native_run=str(native),analysis_root=str(root/'statistics'/variant),
            run_id='polyapipe_integrated_20261003_'+variant,
            input_fingerprint=digest(native/'intake.native.json'),reference_variant=variant)
        config['resources']=dict(config['resources'])
        reference=native/'reference/representative_support_sources.tsv.gz'
        config['resources']['external_pas_source_records']=dict(path=str(reference),sha256=digest(reference))
        config['resources']['support_reference_audit']=dict(path=str(native/'reference_audit.json'),sha256=digest(native/'reference_audit.json'))
        config_path=package/(variant+'.config.json')
        assert not config_path.exists()
        config_path.write_text(json.dumps(config,indent=2)+'\n')
        core_command=[plan['core_python'],'-u',str(code/'polyapipe_core_run.py'),'--run',str(native),
                      '--annotation-pid',str(args.reference_pid),'--workers','7']
        with (root/(variant+'.core.log')).open('w') as log:
            core=subprocess.Popen(core_command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        command=[plan['core_python'],'-u',str(code/'dispatch_ready_statistics.py'),'--config',str(config_path),
                 '--phase-id','full','--core-pid',str(core.pid),'--workers','12']
        with (root/(variant+'.statistics.log')).open('w') as log:
            dispatch=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        launches.append(dict(variant=variant,core_pid=core.pid,dispatch_pid=dispatch.pid,core_argv=core_command,
            dispatch_argv=command,config_sha256=digest(config_path),run_id=config['run_id'],input_fingerprint=config['input_fingerprint']))
        (root/'statistics-launch.json').write_text(json.dumps(dict(launches=launches,complete=False),indent=2)+'\n')
        print(json.dumps(launches[-1]),flush=True)
    marker=root/'runs/reference-replacement.complete.json';wait_for(marker)
    assert json.loads(marker.read_text())['status']=='reference_rebound_core_and_statistics_pending'
    (root/'statistics-launch.json').write_text(json.dumps(dict(launches=launches,complete=True,
        interpretation='launch complete; analyses and manuscript delivery still pending'),indent=2)+'\n')
    print(json.dumps(launches),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--reference-pid',type=int,required=True)
    p.add_argument('--package',type=Path)
    main(p.parse_args())
