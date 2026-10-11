"""Upload small requests/script; retrieve only scalar group aggregates, never matrices."""
import subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
REMOTE='/root/nfsdata/apadetect/analysis_projects/reproapa_marker_ortholog_sharing_20261010'
subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker','mkdir -p '+REMOTE],check=True)
for n in ['gex_requests.tsv.gz','server_source_config.json','extract_gex_remote.py']:
 subprocess.run(['scp','-q',str(HERE/n),'213docker:'+REMOTE+'/'+n],check=True)
command="bash -ic 'PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/apasim_analysis/bin/python "+REMOTE+'/extract_gex_remote.py '+REMOTE+"'"
p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker',command],capture_output=True,text=True)
(HERE/'extract_gex.log').write_text(p.stdout+'\n'+p.stderr)
print(p.stdout,flush=True)
if p.returncode:raise SystemExit(p.returncode)
for n in ['group_GEX_aggregates.tsv.gz','GEX_source_manifest.json']:
 subprocess.run(['scp','-q','213docker:'+REMOTE+'/'+n,str(HERE/n)],check=True)
