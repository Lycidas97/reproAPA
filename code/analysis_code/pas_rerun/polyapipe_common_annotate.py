"""Fresh structural/sequence annotation and v2 support for common-catalogue PAS."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import pickle
import sys
import time

import numpy as np
import pandas as pd
from polyapipe_annotate import REFERENCES, save, sha
from replace_pas_reference import support_fields


def initialize(code, root, reference, bedtools):
    global ANN, SEQUENCE, ROOT, REFS, BEDTOOLS, INTERVALS, INDEX
    sys.path.insert(0, code)
    from postprocessing.core import annotation
    from postprocessing.utils.sequence_qc import compute_sequence_qc
    ANN, SEQUENCE, ROOT, BEDTOOLS = annotation, compute_sequence_qc, Path(root), bedtools
    REFS = pd.read_csv(reference, sep='\t', low_memory=False)
    assert set(REFS.reference_variant) == {'v2_primary'}
    INTERVALS, INDEX = {}, {}


def annotate(task):
    from polyapipe_annotation_fast import BedtoolsCandidateIndex
    from polyapipe_assignment_batch import choose
    species, block, frame = task
    dest = ROOT/'annotation_chunks'/species/f'{block:06d}'; dest.mkdir(parents=True, exist_ok=False)
    chrom = frame.chrom.iloc[0]; key = (species, chrom)
    if species not in INTERVALS:
        with REFERENCES[species]['cache'].open('rb') as f: INTERVALS[species] = pickle.load(f)['intervals']
    if key not in INDEX:
        INDEX.clear()
        intervals = {k:v[v.chrom.eq(chrom)].copy() for k,v in INTERVALS[species].items()}
        INDEX[key] = BedtoolsCandidateIndex(ANN, intervals, BEDTOOLS)
    # The annotation API takes a genomic base; support and analysis take interbase.
    base = frame.annotation_cut_site_0based_base.to_numpy(np.int64)
    peaks = pd.DataFrame(dict(chrom=frame.chrom,start=base,end=base+1,name=frame.atlas_pas_id,score=base,strand=frame.strand))
    pmap = ANN.normalize_pas_features(peaks, species, frame.genome_build.iloc[0])
    np.testing.assert_array_equal(pmap.cut_site, base)
    candidates = INDEX[key].generate(pmap)
    feature = choose(ANN, pmap, candidates)
    seq = SEQUENCE(feature[['pas_id','chrom','cut_site','strand']],fasta_path=str(REFERENCES[species]['fasta']),
        fasta_index_path=str(REFERENCES[species]['fasta'])+'.fai')
    assert seq.sequence_qc_status.eq('computed').all()
    feature = ANN._attach_sequence_qc(feature, seq)
    result = feature.merge(ANN.build_confidence_tiers(feature),on='pas_id',validate='one_to_one')
    result = result.rename(columns={'raw_peak_id':'atlas_pas_id'})
    assert result.atlas_pas_id.tolist() == frame.atlas_pas_id.tolist()
    result['annotation_cut_site_0based_base'] = result.cut_site
    result['cut_site'] = frame.cut_site.to_numpy(np.int64)
    result = support_fields(result, REFS[REFS.species.eq(species)])
    result['support_reference_variant'] = 'v2_primary'
    result['support_distance_basis'] = 'integrated_representative_interbase'
    result['quantifiable'] = True; result['coordinate_supported'] = True
    result['analysis_view'] = 'polyapipe_common'; result['measurement_exclusion'] = ''
    save(result, dest/'features.tsv.gz'); save(seq,dest/'sequence_qc.tsv.gz')
    save(candidates,dest/'annotation_candidates.tsv.gz')
    record = dict(species=species,block=block,rows=len(frame),features_sha256=sha(dest/'features.tsv.gz'))
    (dest/'complete.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def main(root, code, reference, bedtools, equivalence, workers=32):
    from polyapipe_finalize_annotation import finalize
    from polyapipe_annotation_fast import __file__ as fast
    from polyapipe_assignment_batch import __file__ as batch
    root = Path(root); report = json.loads(Path(equivalence).read_text())
    assert {r['species'] for r in report} == {'human','mouse'}
    assert all(r['exact_candidates'] and r['exact_assignments'] and r['exact_primary_tiers'] and
        r['batch_assignments_exact'] and r['sites']>=1000 and r['fast_source_sha256']==sha(fast) and
        r['batch_source_sha256']==sha(batch) for r in report)
    paths = [Path(reference),Path(bedtools),Path(equivalence),Path(__file__),Path(fast),Path(batch)]
    paths.extend(sorted((Path(code)/'postprocessing').rglob('*.py')))
    paths.extend(p for r in REFERENCES.values() for p in r.values())
    coordinates = sorted((root/'studies').glob('*/coordinate_features.tsv.gz')); paths.extend(coordinates)
    inputs = {str(p.resolve()):sha(p) for p in paths}
    lock = root/'annotation_inputs_bedtools.json'; lock.write_text(json.dumps(inputs,indent=2)+'\n')
    all_sites = pd.concat([pd.read_csv(p,sep='\t') for p in coordinates],ignore_index=True)
    tasks = []
    for species, part in all_sites.groupby('species'):
        block = 0
        for _, chrom in part.groupby('chrom',sort=True):
            for start in range(0,len(chrom),1000):
                tasks.append((species,block,chrom.iloc[start:start+1000].reset_index(drop=True))); block+=1
    print(json.dumps(dict(phase='annotation',tasks=len(tasks),rows=len(all_sites),workers=workers)),flush=True)
    with ProcessPoolExecutor(max_workers=workers,initializer=initialize,initargs=(str(code),str(root),str(reference),str(bedtools))) as pool:
        futures = [pool.submit(annotate,task) for task in tasks]
        for n,future in enumerate(as_completed(futures),1):
            future.result()
            if n%100==0: print(json.dumps(dict(phase='annotation',completed=n,total=len(tasks))),flush=True)
    r = finalize(root,verify_lock=True,backend='bedtools',resumed_legacy_chunks=False)
    r['support_reference_variant'] = 'v2_primary'
    r['checks'].pop('external_original_source_coordinates')
    r['checks']['integrated_v2_representative_interbase_support'] = True
    (root/'annotation.complete.json').write_text(json.dumps(r,indent=2)+'\n')
    (root/'intake-complete.json').write_text(json.dumps(dict(status='READY_FOR_CORE',
        source_id='polyapipe_common',annotation_receipt_sha256=sha(root/'annotation.complete.json'),
        full_main_statistics_complete=False,time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['root','code','reference','bedtools','equivalence']:p.add_argument('--'+name,required=True)
    p.add_argument('--workers',type=int,default=32);a=p.parse_args()
    main(a.root,a.code,a.reference,a.bedtools,a.equivalence,a.workers)
