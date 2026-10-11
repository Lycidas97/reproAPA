"""Read-only extraction from hash-pinned whole genome FASTA. Full aUTR only."""
from pathlib import Path
import csv,hashlib,json,sys

def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def reverse_complement(seq):return seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
writer=None
for species in ('human','mouse'):
 spec=CONFIG['resources'][species+'_genome'];path=Path(spec['path']);observed=sha(path);assert observed==spec['sha256']
 print(json.dumps(dict(species=species,path=str(path),sha256=observed,fai_sha256=sha(str(path)+'.fai'))),file=sys.stderr,flush=True)
 index={}
 with Path(str(path)+'.fai').open() as fi:
  for line in fi:
   chrom,length,offset,nb,nw=line.rstrip().split('\t')[:5];index[chrom]=tuple(map(int,(length,offset,nb,nw)))
 aliases={c.removeprefix('chr'):c for c in index};assert len(aliases)==len(index)
 with path.open('rb') as fa:
  for q in QUERIES:
   if q['species']!=species:continue
   c=aliases[str(q['chrom']).removeprefix('chr')];lo,hi=int(q['start0']),int(q['end0']);length,offset,nb,nw=index[c];assert 0<=lo<hi<=length
   begin=offset+(lo//nb)*nw+lo%nb;end=offset+((hi-1)//nb)*nw+(hi-1)%nb+1
   fa.seek(begin);seq=fa.read(end-begin).replace(b'\n',b'').replace(b'\r',b'').decode().upper();assert len(seq)==hi-lo
   if q['strand']=='-':seq=reverse_complement(seq)
   row=q|dict(reference_contig=c,sequence=seq,reference_sha256=observed)
   if writer is None:writer=csv.DictWriter(sys.stdout,fieldnames=list(row),delimiter='\t');writer.writeheader()
   writer.writerow(row)
