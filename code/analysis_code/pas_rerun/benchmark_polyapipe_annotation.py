"""Profile and compare complete real annotation blocks without touching production."""
import argparse
import cProfile
import importlib.util
import json
from pathlib import Path
import pickle
import pstats
import sys
import time

import numpy as np
import pandas as pd

from polyapipe_annotation_fast import CandidateIndex, BedtoolsCandidateIndex, choose_with_cached_context


def main(args):
    sys.path.insert(0, str(args.code))
    from postprocessing.core import annotation as ann
    from postprocessing.utils.sequence_qc import compute_sequence_qc
    spec=importlib.util.spec_from_file_location('driver',args.code/'polyapipe_annotate.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    args.out.mkdir(parents=True,exist_ok=True)
    reports=[]
    for species, study in [('human','GSE178325'),('mouse','GSE103221')]:
        frame=pd.read_csv(args.run/'studies'/study/'coordinate_features.tsv.gz',sep='\t')
        chrom='chr1' if 'chr1' in set(frame.chrom) else '1'
        frame=frame[frame.chrom.eq(chrom)]
        # Spread queries over the chromosome: both strands, overlaps, gaps, TES.
        frame=frame.iloc[np.linspace(0,len(frame)-1,args.sites,dtype=int)].reset_index(drop=True)
        with driver.REFERENCES[species]['cache'].open('rb') as f:
            all_intervals=pickle.load(f)['intervals']
        intervals={k:v[v.chrom.eq(chrom)].copy() for k,v in all_intervals.items()}
        peaks=pd.DataFrame(dict(chrom=frame.chrom,start=np.maximum(0,frame.cut_site-1),
            end=np.maximum(1,frame.cut_site),name=frame.atlas_pas_id,score=frame.cut_site,strand=frame.strand))
        pmap=ann.normalize_pas_features(peaks,species,frame.genome_build.iloc[0])
        times={}
        profile=cProfile.Profile();profile.enable();start=time.perf_counter()
        expected=ann.generate_annotation_candidates(pmap,intervals)
        times['original_candidates_profiled']=time.perf_counter()-start
        start=time.perf_counter();features=ann.choose_primary_assignments(pmap,expected)
        times['original_assignments_profiled']=time.perf_counter()-start
        profile.disable();profile.dump_stats(str(args.out/f'{species}.prof'))
        with (args.out/f'{species}-profile.txt').open('w') as f:
            pstats.Stats(profile,stream=f).sort_stats('cumulative').print_stats(35)
        start=time.perf_counter();index=CandidateIndex(ann,intervals)
        times['fast_index_build']=time.perf_counter()-start
        start=time.perf_counter();actual=index.generate(pmap)
        times['fast_candidates']=time.perf_counter()-start
        # Exact row order, columns, dtypes and values; approximate floats forbidden.
        pd.testing.assert_frame_equal(expected,actual,check_exact=True)
        start=time.perf_counter();new_features=ann.choose_primary_assignments(pmap,actual)
        times['assignments_unprofiled']=time.perf_counter()-start
        pd.testing.assert_frame_equal(features,new_features,check_exact=True)
        start=time.perf_counter()
        unprofiled=ann.generate_annotation_candidates(pmap,intervals)
        times['original_candidates_unprofiled']=time.perf_counter()-start
        pd.testing.assert_frame_equal(unprofiled,actual,check_exact=True)
        if args.bedtools:
            start=time.perf_counter();bed=BedtoolsCandidateIndex(ann,intervals,args.bedtools)
            times['bedtools_reference_build']=time.perf_counter()-start
            start=time.perf_counter();bed_candidates=bed.generate(pmap)
            times['bedtools_candidates']=time.perf_counter()-start
            pd.testing.assert_frame_equal(expected,bed_candidates,check_exact=True)
            start=time.perf_counter();bed_features=choose_with_cached_context(ann,pmap,bed_candidates)
            times['cached_context_assignments']=time.perf_counter()-start
            pd.testing.assert_frame_equal(features,bed_features,check_exact=True)
        seq=compute_sequence_qc(features[['pas_id','chrom','cut_site','strand']],
            fasta_path=str(driver.REFERENCES[species]['fasta']),
            fasta_index_path=str(driver.REFERENCES[species]['fasta'])+'.fai')
        tiers=ann.build_confidence_tiers(ann._attach_sequence_qc(features,seq))
        new_tiers=ann.build_confidence_tiers(ann._attach_sequence_qc(new_features,seq))
        pd.testing.assert_frame_equal(tiers,new_tiers,check_exact=True)
        report=dict(species=species,sites=len(frame),candidate_rows=len(actual),
            assignments=features.primary_assignment_status.value_counts().to_dict(),seconds=times,
            exact_candidates=True,exact_assignments=True,exact_primary_tiers=True,
            bedtools_and_context_exact=bool(args.bedtools),
            candidate_speedup=times['original_candidates_unprofiled']/times['fast_candidates'])
        reports.append(report);print(json.dumps(report),flush=True)
        (args.out/'benchmark.json').write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    p.add_argument('--code',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--sites',type=int,default=1000)
    p.add_argument('--bedtools');main(p.parse_args())
