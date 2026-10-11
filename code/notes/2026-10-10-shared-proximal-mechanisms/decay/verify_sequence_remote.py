"""Independent reference extraction, line-by-line; compare sequence hashes, not motifs."""
from pathlib import Path
import hashlib,json,subprocess
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
config=json.loads((ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json').read_text())
s=pd.read_csv(HERE/'full_aUTR_sequences.tsv',sep='\t')
rows=s[['species','gene_id','chrom','strand','start0','end0']].to_dict('records')
remote='''
import csv,hashlib,sys
from pathlib import Path
out=csv.writer(sys.stdout,delimiter='\\t');out.writerow(['species','gene_id','sequence_sha256','length','first12','last12'])
for sp in ('human','mouse'):
 path=Path(CONFIG['resources'][sp+'_genome']['path'])
 with path.open('rb') as file:assert hashlib.file_digest(file,'sha256').hexdigest()==CONFIG['resources'][sp+'_genome']['sha256']
 idx={}
 for raw in Path(str(path)+'.fai').read_text().splitlines():
  c,l,o,b,w=raw.split('\\t')[:5];idx[c.removeprefix('chr')]=(int(l),int(o),int(b),int(w))
 with path.open('rb') as file:
  for q in ROWS:
   if q['species']!=sp:continue
   l,o,b,w=idx[str(q['chrom']).removeprefix('chr')];start,end=int(q['start0']),int(q['end0']);parts=[];pos=start
   while pos<end:
    rowno=pos//b;col=pos%b;take=min(b-col,end-pos);file.seek(o+rowno*w);line=file.readline().strip();chunk=line[col:col+take];assert len(chunk)==take;parts.append(chunk);pos+=take
   text=b''.join(parts).decode().upper()
   if q['strand']=='-':text=''.join({'A':'T','T':'A','G':'C','C':'G','N':'N'}[x] for x in text[::-1])
   assert len(text)==end-start;out.writerow([sp,q['gene_id'],hashlib.sha256(text.encode()).hexdigest(),len(text),text[:12],text[-12:]])
'''
code='CONFIG='+repr(config)+'\nROWS='+repr(rows)+'\n'+remote
with (HERE/'independent_reference_sequences.tsv').open('wb') as out,(HERE/'independent_reference.log').open('wb') as log:
 subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker','PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/scanpy/bin/python -'],input=code.encode(),stdout=out,stderr=log,check=True)
f=pd.read_csv(HERE/'independent_reference_sequences.tsv',sep='\t').set_index(['species','gene_id']);assert len(f)==len(s)
for row in s.to_dict('records'):
 z=f.loc[(row['species'],row['gene_id'])];assert z.sequence_sha256==hashlib.sha256(row['sequence'].encode()).hexdigest();assert int(z.length)==len(row['sequence']);assert z.first12==row['sequence'][:12] and z.last12==row['sequence'][-12:]
(HERE/'independent_reference_validation.json').write_text(json.dumps({'status':'passed','full_sequences_compared':len(s),'independent_read_method':'read each indexed FASTA line and slice each line separately','hash_pinned_whole_FASTA':True,'mismatches':0},indent=2)+'\n');print('Independent reference sequence validation PASS:',len(s))
