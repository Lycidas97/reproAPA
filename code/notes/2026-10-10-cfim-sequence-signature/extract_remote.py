"""Read pinned reference sequences; no remote writes, no labels and no liftover."""
from pathlib import Path
import csv
import hashlib
import json
import sys



def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN','TGCAN'))[::-1]


writer=None
for species in ['human','mouse']:
    spec=CONFIG['resources'][species+'_genome'];path=Path(spec['path'])
    observed=sha(path);assert observed==spec['sha256'],(path,observed,spec['sha256'])
    provenance=dict(species=species,path=str(path),sha256=observed,fai_sha256=sha(str(path)+'.fai'))
    print(json.dumps(provenance),file=sys.stderr,flush=True)
    index={}
    with Path(str(path)+'.fai').open() as fi:
        for line in fi:
            c,length,offset,line_bases,line_width=line.rstrip().split('\t')[:5]
            index[c]=tuple(map(int,[length,offset,line_bases,line_width]))
    with path.open('rb') as fa:
        aliases={c.removeprefix('chr'):c for c in index};assert len(aliases)==len(index)
        for q in QUERIES:
            if q['species']!=species:continue
            c=aliases[q['chrom'].removeprefix('chr')];lo,hi=q['start0'],q['end0']
            length,offset,line_bases,line_width=index[c]
            assert 0<=lo<hi<=length and hi-lo==100
            begin=offset+(lo//line_bases)*line_width+lo%line_bases
            end=offset+((hi-1)//line_bases)*line_width+(hi-1)%line_bases+1
            fa.seek(begin);sequence=fa.read(end-begin).replace(b'\n',b'').replace(b'\r',b'').decode('ascii').upper()
            assert len(sequence)==100
            if q['strand']=='-':sequence=reverse_complement(sequence)
            row=q|dict(reference_contig=c,sequence=sequence,reference_sha256=observed)
            if writer is None:writer=csv.DictWriter(sys.stdout,fieldnames=list(row),delimiter='\t');writer.writeheader()
            writer.writerow(row)
