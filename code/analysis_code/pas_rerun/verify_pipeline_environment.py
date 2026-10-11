"""Read-only resource and runtime preflight; does not certify scientific outputs."""
import argparse
import importlib
from pathlib import Path
import sys
from pipeline import read,write,digest
from prepare_pipeline import make_jobs


def run(environment,out):
    env=read(environment);results=[];modules={};node_modules={}
    for name in ['numpy','pandas','scipy','sklearn','pysam','pyBigWig','matplotlib']:
        try:
            module=importlib.import_module(name);modules[name]=getattr(module,'__version__','available')
        except ImportError as error:modules[name]={'error':str(error)}
    jobs=make_jobs(env,Path(__file__).parent,sys.executable)
    for name in sorted({Path(j['script']).stem for j in jobs.values()}|{'render_from_bundle'}):
        try:
            importlib.import_module(name);node_modules[name]='imported'
        except Exception as error:node_modules[name]={'error':repr(error)}
    for group in ['resources','source_overrides']:
        for name,item in env[group].items():
            p=Path(item['path'])
            actual=digest(p) if p.is_file() else None
            results.append(dict(group=group,name=name,path=str(p),expected=item['sha256'],actual=actual,match=actual==item['sha256']))
            if len(results)%20==0:print('verified',len(results),flush=True)
    ok=(all(r['match'] for r in results) and all(isinstance(v,str) for v in modules.values())
        and all(v=='imported' for v in node_modules.values()))
    write(out,dict(status='PASS' if ok else 'FAIL',scope='static resources and runtime only; no scientific run',
                   environment_sha256=digest(environment),python=sys.executable,python_version=sys.version,
                   modules=modules,node_modules=node_modules,files=results))
    print('PASS' if ok else 'FAIL',len(results),flush=True)
    if not ok:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--environment',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.environment,a.out)
