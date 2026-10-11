"""Read-only replay of fixed mm10 windows using the already pinned resources."""
from pathlib import Path
import hashlib
import json
import subprocess

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CONFIG = ROOT / 'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json'

REMOTE = '''
from pathlib import Path
import csv,hashlib,json,sys
import numpy as np
import pyBigWig,pysam
for key in ['cis_mm10_fasta','cis_phastcons']:
    path=Path(RESOURCES[key]['path'])
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):digest.update(block)
    actual=digest.hexdigest()
    assert actual==RESOURCES[key]['sha256'],(key,actual)
    print(json.dumps(dict(key=key,path=str(path),sha256=actual,bytes=path.stat().st_size)),file=sys.stderr,flush=True)
bw=pyBigWig.open(RESOURCES['cis_phastcons']['path'])
fa=pysam.FastaFile(RESOURCES['cis_mm10_fasta']['path'])
sizes=bw.chroms()
columns=['id','mapped_chrom','mapped_cut','start0','end0','window_bases','finite_bases','finite_fraction','phastcons_replayed','gc_fraction','ACGT_bases','sequence_length']
writer=csv.DictWriter(sys.stdout,fieldnames=columns,delimiter='\\t');writer.writeheader()
for query in QUERIES:
    assert query['mapping_status']=='unique'
    chrom=query['mapped_chrom'];cut=int(query['mapped_cut'])
    requested_lo,requested_hi=cut-100,cut+50
    lo=max(0,requested_lo);hi=min(sizes[chrom],requested_hi)
    values=np.asarray(bw.values(chrom,lo,hi),float)
    finite=np.isfinite(values)
    sequence=fa.fetch(chrom,lo,hi).upper()
    acgt=sum(sequence.count(base) for base in 'ACGT')
    writer.writerow(dict(id=query['id'],mapped_chrom=chrom,mapped_cut=cut,start0=requested_lo,end0=requested_hi,
        window_bases=150,finite_bases=int(finite.sum()),finite_fraction=float(finite.sum()/150),
        phastcons_replayed=float(np.mean(values[finite])) if finite.any() else 'NA',
        gc_fraction=(sequence.count('G')+sequence.count('C'))/acgt if acgt else 'NA',
        ACGT_bases=acgt,sequence_length=len(sequence)))
bw.close();fa.close()
'''


def main():
    config = json.loads(CONFIG.read_text())
    query = pd.read_csv(HERE / 'fixed_window_queries.tsv', sep='\t', float_precision='round_trip')
    assert query.id.is_unique and query.mapping_status.eq('unique').all()
    code = 'RESOURCES=' + repr(config['resources']) + '\nQUERIES=' + repr(query.to_dict('records')) + '\n' + REMOTE
    command = 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/apasim_analysis/bin/python -'
    partial = HERE / 'fixed_window_metadata.partial.tsv'
    with partial.open('wb') as output, (HERE / 'fixed_window_remote.log').open('wb') as log:
        subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', '213docker', command],
                       input=code.encode(), stdout=output, stderr=log, check=True)
    metadata = pd.read_csv(partial, sep='\t', float_precision='round_trip')
    assert len(metadata) == len(query) and metadata.id.is_unique
    merged = query.merge(metadata, on='id', validate='one_to_one')
    assert merged.phastcons.isna().eq(merged.phastcons_replayed.isna()).all()
    error = float(np.max(np.abs(merged.phastcons.dropna() - merged.phastcons_replayed.dropna()))) if merged.phastcons.notna().any() else 0.
    assert error < 1e-12, error
    partial.rename(HERE / 'fixed_window_metadata.tsv')
    report = dict(status='PASS', fixed_window='mm10 genomic [cut-100,cut+50)',
                  n_sites=len(metadata), max_phastcons_replay_error=error,
                  minimum_finite_fraction=float(metadata.finite_fraction.min()),
                  full_150bp_sequences=int(metadata.sequence_length.eq(150).sum()),
                  all_150bp_ACGT=int(metadata.ACGT_bases.eq(150).sum()),
                  resources=[json.loads(line) for line in (HERE / 'fixed_window_remote.log').read_text().splitlines() if line.startswith('{')],
                  query_sha256=hashlib.sha256((HERE / 'fixed_window_queries.tsv').read_bytes()).hexdigest(),
                  config_sha256=hashlib.sha256(CONFIG.read_bytes()).hexdigest())
    (HERE / 'window_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
