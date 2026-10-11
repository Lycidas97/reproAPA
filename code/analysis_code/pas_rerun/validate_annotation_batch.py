"""Compare the batch backend with frozen legacy outputs; never write production."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import sys
import time

import numpy as np
import pandas as pd

from polyapipe_annotation_fast import BedtoolsCandidateIndex
import polyapipe_annotation_fast as fast
import polyapipe_assignment_batch as batch


def main(args):
    sys.path.insert(0, str(args.code))
    from postprocessing.core import annotation as ann
    from postprocessing.utils.sequence_qc import compute_sequence_qc
    spec=importlib.util.spec_from_file_location('driver',args.code/'polyapipe_annotate.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    args.out.mkdir(parents=True,exist_ok=True)
    batch_hash=hashlib.sha256(Path(batch.__file__).read_bytes()).hexdigest()
    fast_hash=hashlib.sha256(Path(fast.__file__).read_bytes()).hexdigest()
    reports=[]
    for species,study in [('human','GSE178325'),('mouse','GSE103221')]:
        cache=args.out/f'{species}-legacy.pickle'
        if cache.exists():
            with cache.open('rb') as f: saved=pickle.load(f)
        else:
            frame=pd.read_csv(args.run/'studies'/study/'coordinate_features.tsv.gz',sep='\t')
            chrom='chr1' if 'chr1' in set(frame.chrom) else '1'
            frame=frame[frame.chrom.eq(chrom)]
            frame=frame.iloc[np.linspace(0,len(frame)-1,1000,dtype=int)].reset_index(drop=True)
            with driver.REFERENCES[species]['cache'].open('rb') as f:
                intervals=pickle.load(f)['intervals']
            intervals={k:v[v.chrom.eq(chrom)].copy() for k,v in intervals.items()}
            peaks=pd.DataFrame(dict(chrom=frame.chrom,start=np.maximum(0,frame.cut_site-1),
                end=np.maximum(1,frame.cut_site),name=frame.atlas_pas_id,score=frame.cut_site,strand=frame.strand))
            pmap=ann.normalize_pas_features(peaks,species,frame.genome_build.iloc[0])
            start=time.perf_counter();candidates=ann.generate_annotation_candidates(pmap,intervals)
            candidate_seconds=time.perf_counter()-start
            start=time.perf_counter();expected=ann.choose_primary_assignments(pmap,candidates)
            assignment_seconds=time.perf_counter()-start
            bed=BedtoolsCandidateIndex(ann,intervals,args.bedtools)
            start=time.perf_counter();actual_candidates=bed.generate(pmap)
            bed_seconds=time.perf_counter()-start
            pd.testing.assert_frame_equal(candidates,actual_candidates,check_exact=True)
            seq=compute_sequence_qc(expected[['pas_id','chrom','cut_site','strand']],
                fasta_path=str(driver.REFERENCES[species]['fasta']),
                fasta_index_path=str(driver.REFERENCES[species]['fasta'])+'.fai')
            saved=dict(pmap=pmap,candidates=candidates,expected=expected,seq=seq,
                candidate_seconds=candidate_seconds,assignment_seconds=assignment_seconds,bed_seconds=bed_seconds)
            with cache.open('wb') as f:pickle.dump(saved,f)
        start=time.perf_counter();actual=batch.choose(ann,saved['pmap'],saved['candidates'])
        seconds=time.perf_counter()-start
        expected=saved['expected']
        pd.testing.assert_frame_equal(expected,actual,check_exact=True)
        pd.testing.assert_frame_equal(
            ann.build_confidence_tiers(ann._attach_sequence_qc(expected,saved['seq'])),
            ann.build_confidence_tiers(ann._attach_sequence_qc(actual,saved['seq'])),check_exact=True)
        report=dict(species=species,sites=len(actual),candidate_rows=len(saved['candidates']),
            exact_candidates=True,exact_assignments=True,exact_primary_tiers=True,batch_assignments_exact=True,
            seconds={k:saved[k] for k in ['candidate_seconds','assignment_seconds','bed_seconds']},
            batch_assignment_seconds=seconds,
            batch_source_sha256=batch_hash,fast_source_sha256=fast_hash)
        assert hashlib.sha256(Path(batch.__file__).read_bytes()).hexdigest()==batch_hash
        assert hashlib.sha256(Path(fast.__file__).read_bytes()).hexdigest()==fast_hash
        reports.append(report);print(json.dumps(report),flush=True)
        (args.out/'benchmark.json').write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--code',type=Path,required=True)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--bedtools',required=True);main(p.parse_args())
