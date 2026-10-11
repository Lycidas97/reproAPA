from pathlib import Path
import json,subprocess
import pandas as pd
from prepare_events import HERE,CONFIG,save
config=json.loads(CONFIG.read_text());q=pd.read_csv(HERE/'sequence_queries.tsv',sep='\t')
code='CONFIG='+repr(config)+'\nQUERIES='+repr(q.to_dict('records'))+'\n'+(HERE/'extract_sequences.py').read_text()
partial=HERE/'full_aUTR_sequences.partial.tsv'
with partial.open('wb') as out,(HERE/'extract.log').open('wb') as log:
 subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker','PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/scanpy/bin/python -'],input=code.encode(),stdout=out,stderr=log,check=True)
partial.rename(HERE/'full_aUTR_sequences.tsv')
s=pd.read_csv(HERE/'full_aUTR_sequences.tsv',sep='\t');assert len(s)==len(q) and not s.duplicated(['species','gene_id']).any()
s['complete_ACGT']=s.sequence.str.len().eq(s.end0-s.start0)&s.sequence.str.fullmatch('[ACGT]+');save(s,'full_aUTR_sequences.tsv')
i=pd.read_csv(HERE/'gene_interval_inventory.before_sequence.tsv.gz',sep='\t');seq_lookup=s.set_index(['species','gene_id'])
for k,r in i[i.status.eq('sequence_pending')].iterrows():
 z=seq_lookup.loc[(r.species,r.gene_id)];i.loc[k,'status']='eligible' if z.complete_ACGT else 'ambiguous_or_incomplete_sequence'
 i.loc[k,'GC']=sum(z.sequence.count(x) for x in 'GC')/len(z.sequence)
save(i,'gene_interval_inventory.tsv.gz');save(i[i.foreground_002],'core_event_summary.tsv');print(i.groupby(['species','status']).size().to_string());print((HERE/'extract.log').read_text())
