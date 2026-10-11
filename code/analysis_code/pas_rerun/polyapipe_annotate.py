"""Apply the pinned manuscript structural annotation to new polyApipe coordinates.

Run on the analysis host, where the original genomes and reference caches reside.
The native count policy remains unchanged; sequence QC affects primary eligibility.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
import numpy as np
import pandas as pd

SOURCES = ['polyadb_main', 'polyasite_broad', 'gencode_polya']
BASE = Path('/root/nfsdata/apadetect/data')
CACHE = BASE / 'reprogramming_differentiation_apa/analyses/apadetect_reproapa_study_atlas_v5_structural_priority_6aabc22_2026-08-14/reference_cache'
REF = Path('/root/nfsdata/multica/worktrees/apadetect/DOC-145-reproapa-minimal-validation/analyses/DOC-209-reproapa-external-pas-validation')
REFERENCES = dict(
    mouse=dict(cache=CACHE / 'reference_intervals_0d8f4d3b9c5a16bd599f.pickle',
               fasta=BASE / 'genome/mouse/refdata-gex-GRCm39-2024-A/fasta/genome.fa'),
    human=dict(cache=CACHE / 'reference_intervals_a8e01ce68c1efa471385.pickle',
               fasta=BASE / 'genome/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(frame, path):
    frame.to_csv(path, sep='\t', index=False,
                 compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)


def reference_index(out):
    manifest = json.loads((REF / 'run_manifest.json').read_text())
    source = REF / 'external_pas_catalog_sources.tsv.gz'
    expected = {x['path']: x['sha256'] for x in manifest['outputs']}[source.name]
    assert sha(source) == expected, 'External original-source coordinates changed'
    index, counts = {}, {}
    cols = ['species', 'target_build', 'source_family', 'target_chrom', 'target_strand',
            'target_cut_site', 'primary_contig']
    for part in pd.read_csv(source, sep='\t', usecols=cols, chunksize=500000):
        for key, group in part.groupby(['species', 'source_family', 'target_chrom', 'target_strand'], sort=False):
            assert group.target_build.eq('GRCh38_2024_A' if key[0] == 'human' else 'GRCm39_2024_A').all()
            assert group.primary_contig.astype(str).str.lower().eq('true').all()
            index.setdefault(key, []).append(group.target_cut_site.to_numpy(np.int64))
            counts[key[:2]] = counts.get(key[:2], 0) + len(group)
    result = {key: np.unique(np.concatenate(value)) for key, value in index.items()}
    wanted = pd.read_csv(REF / 'external_reference_manifest.tsv', sep='\t')
    for row in wanted.itertuples(index=False):
        assert counts[(row.species, row.source_family)] == row.retained_target_rows
    with (out / 'reference_index.pickle').open('wb') as f:
        pickle.dump(result, f)
    (out / 'reference_audit.json').write_text(json.dumps(dict(source=str(source), sha256=expected,
        records={':'.join(k): v for k, v in counts.items()}, coordinate_system='0-based-interbase'), indent=2)+'\n')


def initialize(code, root, backend='legacy', bedtools=None):
    global ann, sequence, INTERVALS, EXTERNAL, ROOT, BACKEND, BEDTOOLS, CANDIDATE_INDEX
    sys.path.insert(0, code)
    from postprocessing.core import annotation as ann
    from postprocessing.utils.sequence_qc import compute_sequence_qc as sequence
    ROOT = Path(root)
    INTERVALS = {}
    BACKEND, BEDTOOLS, CANDIDATE_INDEX = backend, bedtools, {}
    with (ROOT / 'reference_index.pickle').open('rb') as f:
        EXTERNAL = pickle.load(f)


def annotate(task):
    species, block_id, frame = task
    dest = ROOT / 'annotation_chunks' / species / f'{block_id:06d}'
    dest.mkdir(parents=True, exist_ok=True)
    marker = dest / 'complete.json'
    if marker.exists():
        record = json.loads(marker.read_text())
        assert record['features_sha256'] == sha(dest / 'features.tsv.gz')
        assert record['rows'] == len(frame), 'Chunk boundaries changed'
        assert all((dest/n).is_file() for n in ['sequence_qc.tsv.gz','annotation_candidates.tsv.gz'])
        return record
    started = time.time()
    chrom = frame.chrom.iloc[0]
    if species not in INTERVALS:
        with REFERENCES[species]['cache'].open('rb') as f:
            INTERVALS[species] = pickle.load(f)['intervals']
    timings = {}
    indexed = time.perf_counter()
    if BACKEND == 'bedtools':
        from polyapipe_annotation_fast import BedtoolsCandidateIndex
        from polyapipe_assignment_batch import choose
        key = (species, chrom)
        if key not in CANDIDATE_INDEX:
            CANDIDATE_INDEX.clear()  # Keep one chromosome per worker.
            intervals = {k: v[v.chrom.eq(chrom)].copy() for k, v in INTERVALS[species].items()}
            CANDIDATE_INDEX[key] = BedtoolsCandidateIndex(ann, intervals, BEDTOOLS)
    else:
        intervals = {k: v[v.chrom.eq(chrom)].copy() for k, v in INTERVALS[species].items()}
    timings['reference_index'] = time.perf_counter()-indexed
    peaks = pd.DataFrame(dict(chrom=frame.chrom, start=np.maximum(0, frame.cut_site-1),
        end=np.maximum(1, frame.cut_site), name=frame.atlas_pas_id, score=frame.cut_site, strand=frame.strand))
    pmap = ann.normalize_pas_features(peaks, species, frame.genome_build.iloc[0])
    tick = time.perf_counter()
    candidates = (CANDIDATE_INDEX[key].generate(pmap) if BACKEND == 'bedtools'
                  else ann.generate_annotation_candidates(pmap, intervals, cut_site_window_bp=25))
    timings['candidates'] = time.perf_counter()-tick
    tick = time.perf_counter()
    feature = (choose(ann, pmap, candidates) if BACKEND == 'bedtools'
               else ann.choose_primary_assignments(pmap, candidates))
    timings['assignments'] = time.perf_counter()-tick
    tick = time.perf_counter()
    seq = sequence(feature[['pas_id', 'chrom', 'cut_site', 'strand']],
        fasta_path=str(REFERENCES[species]['fasta']),
        fasta_index_path=str(REFERENCES[species]['fasta'])+'.fai')
    assert seq.sequence_qc_status.eq('computed').all(), seq.sequence_qc_status.value_counts().to_dict()
    feature = ann._attach_sequence_qc(feature, seq)
    tiers = ann.build_confidence_tiers(feature)
    result = feature.merge(tiers, on='pas_id', validate='one_to_one')
    result = result.rename(columns={'raw_peak_id': 'atlas_pas_id'})
    timings['sequence_and_tiers'] = time.perf_counter()-tick
    tick = time.perf_counter()
    for source in SOURCES:
        distances = np.full(len(result), -1, np.int64)
        cuts = np.full(len(result), -1, np.int64)
        for (c, strand), ids in result.groupby(['chrom', 'strand']).groups.items():
            alias = 'chr'+c if c.isdigit() or c in ['X', 'Y'] else ('chrM' if c == 'MT' else c)
            values = EXTERNAL.get((species, source, alias, strand))
            if values is None:
                continue
            query = result.loc[ids, 'cut_site'].to_numpy(np.int64)
            right = np.clip(np.searchsorted(values, query), 0, len(values)-1)
            left = np.maximum(right-1, 0)
            ix = np.where(abs(query-values[left]) <= abs(query-values[right]), left, right)
            distances[ids] = abs(query-values[ix]); cuts[ids] = values[ix]
            for k in np.linspace(0, len(ids)-1, min(3, len(ids)), dtype=int):
                assert abs(values-query[k]).min() == distances[ids[k]]
        result[source+'_distance_nt'] = distances
        result[source+'_nearest_cut'] = cuts
    for window in [25, 50, 100]:
        result[f'matched_source_count_{window}nt'] = sum(result[s+'_distance_nt'].between(0, window).astype(int) for s in SOURCES)
    timings['external_support'] = time.perf_counter()-tick
    tick = time.perf_counter()
    save(result, dest / 'features.tsv.gz')
    save(seq, dest / 'sequence_qc.tsv.gz')
    # Candidate-level provenance is retained, including excluded and ambiguous genes.
    save(candidates, dest / 'annotation_candidates.tsv.gz')
    timings['write_outputs'] = time.perf_counter()-tick
    record = dict(species=species, block=block_id, rows=len(frame), seconds=time.time()-started,
                  features_sha256=sha(dest / 'features.tsv.gz'), backend=BACKEND, timings=timings)
    temporary = dest/'complete.tmp.json'
    temporary.write_text(json.dumps(record, indent=2)+'\n');temporary.replace(marker)
    return record


def main(args):
    root = args.run.resolve()
    if args.cpu_manifest:
        allowed=json.loads(args.cpu_manifest.read_text())['allowed_cpu_ids']
        os.sched_setaffinity(0, set(allowed))
        assert args.workers <= len(allowed)
    # An interrupted run can resume only with byte-identical code and coordinates.
    inputs = {str(p): sha(p) for p in sorted(Path(args.code).rglob('*.py'))}
    inputs[str(Path(__file__).resolve())] = sha(__file__)
    from polyapipe_finalize_annotation import finalize, __file__ as assembly_source
    inputs[str(Path(assembly_source).resolve())] = sha(assembly_source)
    if args.backend == 'bedtools':
        from polyapipe_annotation_fast import __file__ as fast_source
        from polyapipe_assignment_batch import __file__ as batch_source
        assert args.equivalence_report and args.bedtools
        report=json.loads(args.equivalence_report.read_text())
        assert {r['species'] for r in report} == {'human','mouse'}
        assert all(r['batch_assignments_exact'] and r['exact_candidates'] and
                   r['exact_assignments'] and r['exact_primary_tiers'] and r['sites']>=1000 and
                   r['batch_source_sha256']==sha(batch_source) and
                   r['fast_source_sha256']==sha(fast_source) for r in report)
        for p in [Path(fast_source),Path(batch_source),Path(args.bedtools),args.equivalence_report]:
            inputs[str(p.resolve())]=sha(p)
    for p in sorted((root / 'studies').glob('*/coordinate_features.tsv.gz')):
        inputs[str(p)] = sha(p)
    for info in REFERENCES.values():
        for p in info.values():
            inputs[str(p)] = sha(p)
    lock = root / ('annotation_inputs_bedtools.json' if args.backend == 'bedtools' else 'annotation_inputs.json')
    if args.resume_legacy_lock:
        assert args.backend == 'bedtools' and args.chunk_size == 1000
        legacy=json.loads((root/'annotation_inputs.json').read_text())
        assert all(sha(p)==value for p,value in legacy.items()), 'Legacy inputs or code changed'
        inputs[str(root/'annotation_inputs.json')]=sha(root/'annotation_inputs.json')
    if lock.exists():
        assert json.loads(lock.read_text()) == inputs, 'Changed code/reference/input; use a new run'
    else:
        lock.write_text(json.dumps(inputs, indent=2)+'\n')
    if not (root / 'reference_audit.json').exists():
        reference_index(root)
    frames = [pd.read_csv(p, sep='\t') for p in sorted((root / 'studies').glob('*/coordinate_features.tsv.gz'))]
    all_sites = pd.concat(frames, ignore_index=True)
    tasks = []
    for species, part in all_sites.groupby('species'):
        ordinal = 0
        for _, chrom in part.groupby('chrom', sort=True):
            for start in range(0, len(chrom), args.chunk_size):
                tasks.append((species, ordinal, chrom.iloc[start:start+args.chunk_size].reset_index(drop=True)))
                ordinal += 1
    if args.pilot:
        tasks = [(s, i, f.iloc[:20].copy()) for s, i, f in tasks[:1]]
        # Pilot output uses a separate root, so it cannot be mistaken for a full block.
        pilot_root = root / 'pilot'; pilot_root.mkdir(exist_ok=True)
        (pilot_root / 'reference_index.pickle').symlink_to(root / 'reference_index.pickle')
        initialize(args.code, str(pilot_root), args.backend, args.bedtools)
        print(json.dumps(annotate(tasks[0])), flush=True)
        return
    print(json.dumps(dict(phase='annotation_running',backend=args.backend,workers=args.workers,
                          tasks=len(tasks),rows=len(all_sites))),flush=True)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize,
                             initargs=(args.code, str(root),args.backend,args.bedtools)) as pool:
        futures = [pool.submit(annotate, task) for task in tasks]
        for future in as_completed(futures):
            print(json.dumps(future.result()), flush=True)
    finalize(root, verify_lock=False, backend=args.backend, resumed_legacy_chunks=args.resume_legacy_lock)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--code', required=True)
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--chunk-size', type=int, default=1000)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--backend',choices=['legacy','bedtools'],default='legacy')
    parser.add_argument('--bedtools')
    parser.add_argument('--equivalence-report',type=Path)
    parser.add_argument('--resume-legacy-lock',action='store_true')
    parser.add_argument('--cpu-manifest',type=Path)
    main(parser.parse_args())
