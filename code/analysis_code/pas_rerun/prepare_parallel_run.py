"""Fork a current APADETECT input into a bounded, immutable parallel run.

Only the verified core from this exact input is imported. All downstream jobs
execute anew. The source scheduler, contracts and outputs remain untouched.
"""
import argparse
from pathlib import Path
import shutil

from pipeline import digest,read,write,require,validate_inputs,validate_cached,job_fingerprint

REPLACEMENTS=['pipeline.py','downstream_parallel.py','downstream_bundle.py','downstream_profiles.py',
              'downstream_coordinate_states.py','downstream_intermediate.py','downstream_measurement.py','downstream_metacell.py','downstream_catalog_qc.py']
NODE_WORKERS=dict(state_vectors=28,profiles=14,intermediate=21,measurement=14,metacell=7,pas_qc=7,coordinate_states=4,count_inference=32)


def prepare(source_config,out,cpu_budget=96):
    source_config=Path(source_config).resolve();old=read(source_config);out=Path(out).resolve()
    require(cpu_budget>=33,'Nested parallelism requires at least 33 CPU slots')
    require(not out.exists() and not out.is_relative_to(source_config.parent),'Use a fresh independent output directory')
    require(old['caller']=='APADETECT' and old['reference_variant']=='v2_primary','Only the current APADETECT v2 input can be imported')
    require(not any(p in out.parts for p in ['Fig','main','supplement','production-v15']),'Never write into production')
    validate_inputs(old)
    old_root=Path(old['analysis_root'])
    require(validate_cached(old_root,old,'core',{}),'Verified current-input core required')
    manifest=read(Path(old['prepared_native_source'])/'intake.native.json')
    require(manifest['run_id']==old['run_id'] and manifest.get('analysis_route','base')==old['analysis_route'],
            'Prepared input identity differs from the completed core')
    out.mkdir();code=out/'code';shutil.copytree(source_config.parent/'code',code)
    for name in REPLACEMENTS:shutil.copy2(Path(__file__).with_name(name),code/name)
    contracts=out/'contracts';shutil.copytree(old['contracts_root'],contracts)
    root=out/'statistics';root.mkdir()
    shutil.copytree(old_root/'core',root/'core',symlinks=True)
    cfg=dict(old,analysis_root=str(root),native_run=str(root/'core/native'),contracts_root=str(contracts),
             cpu_budget=cpu_budget,node_workers=NODE_WORKERS,count_workers=32,execution_revision='single-route-parallel-v2',
             jobs={n:dict(j,script=str(code/Path(j['script']).name)) for n,j in old['jobs'].items()},
             code_sha256={str(p):digest(p) for p in sorted(code.glob('*.py'))})
    for path in sorted(contracts.rglob('*')):
        if path.is_file():cfg['code_sha256'][str(path)]=digest(path)
    new_config=out/'pipeline.config.json';write(new_config,cfg)
    write(root/'core/pipeline-receipt.json',dict(fingerprint=job_fingerprint(cfg,'core',{}),
        result_sha256=digest(root/'core/result.json'),dependencies={},
        verified_current_input_core_import=str(old_root/'core'),source_config_sha256=digest(source_config)))
    validate_inputs(cfg);require(validate_cached(root,cfg,'core',{}),'Imported core did not validate')
    write(out/'core-import.json',dict(source_config=str(source_config),source_config_sha256=digest(source_config),
        source_receipt_sha256=digest(old_root/'core/pipeline-receipt.json'),
        source_result_sha256=digest(old_root/'core/result.json'),copied_result_sha256=digest(root/'core/result.json'),
        input_fingerprint=cfg['input_fingerprint'],run_id=cfg['run_id'],
        core_numeric_code_unchanged=True,downstream_results_imported=False,BAM_recount=False,source_modified=False))
    return new_config


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--from-config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--cpu-budget',type=int,default=96)
    a=p.parse_args();print(prepare(a.from_config,a.out,a.cpu_budget))
