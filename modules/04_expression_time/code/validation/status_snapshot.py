"""Read-only process verification and a dated snapshot; never infer completion from PID exit."""
from common import *

registry=json.loads((RUN/'launch_registry.json').read_text())
records=[]
for job,launch in sorted(registry.items()):
    path=RUN/job/'status.json'
    status=json.loads(path.read_text()) if path.exists() else {'phase':'STARTING'}
    proc=Path('/proc')/str(launch['pid'])
    command=(proc/'cmdline').read_bytes() if (proc/'cmdline').exists() else b''
    live=Path(launch['command'][-1]).name.encode() in command
    record=dict(job=job,launch_pid=launch['pid'],process_live=live,status=status,
                attempts=1+len(launch.get('prior_attempts',[])),log=str(RUN/job/'run.log'))
    if live and (proc/'status').exists():
        record['resources']={line.split(':',1)[0]:line.split(':',1)[1].strip()
            for line in (proc/'status').read_text().splitlines() if line.startswith(('VmRSS:','Threads:'))}
    if not live and status['phase'] not in ['STAGE_COMPLETE','FAILED']:
        record['attention']='Worker exited without a terminal receipt; investigate before claiming completion.'
    records.append(record)
snapshot=dict(captured_unix=time.time(),submitted_jobs=len(records),jobs=records,
    all_seven_submitted=len(records)==7,all_full_analyses_complete=False,
    backend='SSH-launched detached processes; per-job dependency and GPU waits are in the analysis scripts')
write_json(RUN/'status_snapshot.json',snapshot)
print(json.dumps(snapshot,ensure_ascii=False,indent=2))
