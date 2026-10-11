"""Run native SCAPE over every frozen atlas library and annotated terminal UTR."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
import fcntl
import hashlib
import json
import os
import pickle
import signal
import subprocess
import threading
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'atlas/scape'
os.environ.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
    OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', TI_ARCH='x64',
    TI_CPU_MAX_NUM_THREADS='2', PYTHONUNBUFFERED='1')
if hasattr(os, 'sched_setaffinity'):
    os.sched_setaffinity(0, json.loads((ROOT / 'protocol/cpu_allocation.json').read_text())['allowed_cpu_ids'])
import numpy as np
import pandas as pd
import pysam
from scipy import sparse

SCAPE = '/root/miniconda3/envs/scape-apa/bin/scape'
EVENT_LOCK = threading.Lock()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def dump(record, path):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    temporary.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def events(kind, **fields):
    record = dict(time=now(), event=kind, **fields)
    with EVENT_LOCK:
        with (OUT / 'events.jsonl').open('a') as f:
            f.write(json.dumps(record) + '\n')
        print(json.dumps(record), flush=True)


def run(command, log, cwd):
    # Each native invocation has its own group: a failed infer may leave a watchdog.
    assert not log.exists(), f'Preserve the previous attempt before retrying: {log}'
    with log.open('w') as f:
        p = subprocess.Popen(command, stdout=f, stderr=subprocess.STDOUT,
            cwd=cwd, env=os.environ.copy(), start_new_session=True)
        dump(dict(time=now(), pid=p.pid, command=command), log.with_suffix('.process.json'))
        events('command_start', pid=p.pid, command=command, log=str(log))
        while p.poll() is None:
            time.sleep(5)
            with log.open('rb') as source:
                source.seek(max(0, log.stat().st_size - 131072))
                fatal = b'Traceback (most recent call last):' in source.read()
            if fatal:
                os.killpg(p.pid, signal.SIGTERM)
                try:
                    p.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                    p.wait()
                events('native_exception_group_stopped', pid=p.pid, log=str(log))
                raise RuntimeError(f'Native traceback (watchdog group stopped): {log}')
    events('command_exit', pid=p.pid, returncode=p.returncode, log=str(log))
    assert p.returncode == 0, f'Exit {p.returncode}: {log}'


def objects(path):
    with path.open('rb') as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                return


def regions(catalog):
    """Coalesce overlaps within a gene/contig/strand, without filling gaps."""
    records = []
    for (gene, chrom, strand), q in catalog.groupby(['gene_id', 'chrom', 'strand'], sort=True):
        intervals = []
        for row in q.sort_values(['annotated_start0', 'annotated_end0']).itertuples():
            start, end = int(row.annotated_start0), int(row.annotated_end0)
            assert 0 <= start < end and strand in ['+', '-']
            if intervals and start <= intervals[-1][1]:
                intervals[-1][1] = max(intervals[-1][1], end)
            else:
                intervals.append([start, end])
        for utr, (start, end) in enumerate(intervals, 1):
            records.append(dict(chrom=chrom, start=start, end=end, strand=strand,
                gene_id=gene, utr_id=utr))
    return pd.DataFrame(records).sort_values(['chrom', 'start', 'end', 'gene_id']).reset_index(drop=True)


def round_robin(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row['study_id'], []).append(row)
    return [group[i] for i in range(max(map(len, groups.values())))
        for group in groups.values() if i < len(group)]


def integer_assignments(models, barcode_index, members):
    mc = pd.Index(sorted(members.metacell_id.unique()))
    lookup = dict(zip(members.barcode, mc.get_indexer(members.metacell_id)))
    assert len(lookup) == len(members) and min(lookup.values()) >= 0
    assert set(barcode_index.values()) == set(lookup) and len(barcode_index) == len(lookup)
    features, rows, cols, values = [], [], [], []
    assigned = background = models_count = 0
    counted_cells = set()
    for model in models:
        chrom, gene, utr, span, strand = model.gene_info_str.split(':')
        start, end = map(int, span.split('-'))
        chrom = chrom if chrom.startswith('chr') else 'chr' + chrom
        offset = len(features)
        assert model.K == len(model.alpha_arr) == len(model.beta_arr)
        assert np.isfinite(model.alpha_arr).all() and strand in ['+', '-']
        for k, alpha in enumerate(model.alpha_arr):
            cut = start + float(alpha) + 1 if strand == '+' else end - float(alpha) + 1
            features.append(dict(site_id=f'SCAPE:{chrom}:{gene}:{utr}:{k}',
                chrom=chrom, strand=strand, gene_id=gene, cut_site=cut,
                model_beta=float(model.beta_arr[k]), misprime='not_available'))
        raw_labels, raw_cb = np.asarray(model.label_arr), np.asarray(model.cb_id_arr)
        assert len(raw_labels) == len(raw_cb) and np.isfinite(raw_labels).all() and np.isfinite(raw_cb).all()
        labels, cb = raw_labels.astype(np.int64), raw_cb.astype(np.int64)
        assert np.array_equal(labels, raw_labels) and np.array_equal(cb, raw_cb)
        assert np.isin(cb, list(barcode_index)).all(), 'Unknown native CB index'
        assert ((labels >= -1) & (labels <= model.K)).all(), 'Unknown native assignment label'
        valid = (labels >= 0) & (labels < model.K)
        assigned += int(valid.sum()); background += int((~valid).sum()); models_count += 1
        barcodes = [barcode_index[x] for x in cb[valid]]
        counted_cells.update(barcodes)
        indices = np.array([lookup[x] for x in barcodes], dtype=np.int64)
        q = pd.DataFrame(dict(mc=indices, site=labels[valid] + offset)).value_counts().reset_index(name='count')
        rows.extend(q.mc); cols.extend(q.site); values.extend(q['count'])
    counts = sparse.coo_matrix((np.asarray(values, dtype=np.int64),
        (np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64))),
        shape=(len(mc), len(features)), dtype=np.int64).tocsr()
    assert int(counts.sum()) == assigned
    f = pd.DataFrame(features, columns=['site_id', 'chrom', 'strand', 'gene_id',
        'cut_site', 'model_beta', 'misprime'])
    assert f.site_id.is_unique
    return f, mc, counts, dict(native_merged_models=models_count,
        native_assigned_molecules=assigned, native_background_assignments=background,
        observed_counted_cells=len(counted_cells), cells_without_native_counts=len(members)-len(counted_cells),
        zero_count_metacells=int((counts.sum(1).A1 == 0).sum()),
        integer_assignment_conservation=True, unknown_barcode_or_site_entries=0)


def library(row, manifest, wait_for_existing=False):
    out = OUT / row['study_id'] / row['library_id']
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'job.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | (0 if wait_for_existing else fcntl.LOCK_NB))
        if (out / 'complete.json').exists():
            return json.loads((out / 'complete.json').read_text())
        assert not (out / 'failure.traceback.txt').exists(), 'Preserve and resolve the failed native attempt before retrying'
        def phase(name, **fields):
            dump(dict(time=now(), phase=name, library_key=row['library_key'], **fields), out / 'status.json')
            events('library_phase', library_key=row['library_key'], phase=name, **fields)
        try:
            assert digest(Path(row['members_path'])) == row['members_sha256']
            assert digest(Path(row['barcodes_path'])) == row['barcodes_sha256']
            members = pd.read_csv(row['members_path'], sep='\t')
            barcodes = pd.read_csv(row['barcodes_path'], names=['barcode']).barcode
            assert len(members) == row['cells'] and members.metacell_id.nunique() == row['metacells']
            assert members.barcode.is_unique and barcodes.is_unique and set(barcodes) == set(members.barcode)
            source = Path(row['bam_path'])
            assert source.stat().st_size == row['bam_bytes'] and source.stat().st_mtime_ns == row['bam_mtime_ns']
            assert digest(Path(row['bai_path'])) == row['bai_sha256']
            poly = ROOT / 'atlas/polyapipe' / row['study_id'] / row['library_id']
            bam = poly / 'filtered.bam'
            phase('waiting_for_exact_member_filtered_bam')
            while not ((poly / 'filter.complete.json').exists() and (poly / 'quickcheck.log').exists()
                       and (poly / 'index.log').exists() and Path(str(bam) + '.bai').exists()):
                time.sleep(30)
            filt = json.loads((poly / 'filter.complete.json').read_text())
            assert bam.stat().st_size == filt['bytes']
            run(['/usr/bin/samtools', 'quickcheck', '-v', str(bam)], out / 'reused_bam.quickcheck.log', out)
            utr_path = Path(manifest['annotations'][row['species']]['path'])
            assert digest(utr_path) == manifest['annotations'][row['species']]['sha256']
            utrs = pd.read_csv(utr_path)
            with pysam.AlignmentFile(bam) as f:
                contigs = set(f.references)
            present = utrs.chrom.map(lambda c: c in contigs or c.removeprefix('chr') in contigs)
            utrs[~present].to_csv(out / 'annotated_utrs_missing_from_bam.csv', index=False)
            dump(dict(time=now(), filtered_bam=str(bam), filtered_bytes=filt['bytes'],
                members_sha256=row['members_sha256'], barcodes_sha256=row['barcodes_sha256'],
                annotated_utrs=len(utrs), available_utrs=int(present.sum()),
                missing_utrs=int((~present).sum()), available_chromosomes=sorted(utrs.loc[present, 'chrom'].unique())),
                out / 'input_audit.json')
            directories = []
            chromosomes = sorted(utrs.loc[present, 'chrom'].unique(),
                key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c))
            for chrom in chromosomes:
                dest = out / chrom; dest.mkdir(exist_ok=True); directories.append(dest)
                if (dest / 'complete.json').exists():
                    continue
                q = utrs[utrs.chrom.eq(chrom)]
                csv = dest / 'utrs.csv'; q.to_csv(csv, index=False)
                if not (dest / 'prepare.complete.json').exists():
                    phase('native_prepare_running', chromosome=chrom, annotated_utrs=len(q))
                    run([SCAPE, 'prepare_input', '--utr_file', str(csv), '--cb_file', row['barcodes_path'],
                        '--bam_file', str(bam), '--output_dir', str(dest), '--chunksize', '100'], dest / 'prepare.log', dest)
                    pkls = sorted((dest / 'pkl_input').glob('*.input.pkl'))
                    assert pkls, 'Native prepare did not produce its declared input shards'
                    model_regions = [info for p in pkls for info, data in objects(p)]
                    assert len(model_regions) == len(set(model_regions))
                    pd.DataFrame(dict(native_gene_info=model_regions)).to_csv(dest / 'prepared_regions.tsv', sep='\t', index=False)
                    dump(dict(time=now(), annotated_utrs=len(q), prepared_utrs=len(model_regions),
                        regions_without_native_input=len(q)-len(model_regions), native_read_threshold='>100',
                        input_shards=len(pkls)), dest / 'prepare.complete.json')
                pkls = sorted((dest / 'pkl_input').glob('*.input.pkl'))
                prepared = json.loads((dest / 'prepare.complete.json').read_text())
                if prepared['prepared_utrs'] == 0:
                    dump(dict(time=now(), status='no_regions_passed_native_input_threshold',
                        native_models=0, sites=0, assigned_molecules=0), dest / 'complete.json')
                    continue
                for pkl in pkls:
                    mark = dest / (pkl.name + '.complete.json')
                    if mark.exists():
                        continue
                    phase('native_infer_running', chromosome=chrom, shard=pkl.name)
                    run([SCAPE, 'infer_pa', '--pkl_input_file', str(pkl), '--output_dir', str(dest)],
                        dest / (pkl.name + '.infer.log'), dest)
                    result = dest / 'pkl_output' / (pkl.name[:-10] + '.res.pkl')
                    assert result.exists()
                    dump(dict(time=now(), output=str(result), bytes=result.stat().st_size), mark)
                phase('native_merge_running', chromosome=chrom)
                run([SCAPE, 'merge_pa', '--output_dir', str(dest), '--utr_merge', 'True'], dest / 'merge.log', dest)
                assert (dest / 'res.gene.pkl').exists()
                dump(dict(time=now(), input_shards=len(pkls), status='native_prepare_infer_merge_completed'), dest / 'complete.json')
            dump(dict(time=now(), chromosomes=len(directories), scientific_comparison_status='pending'), out / 'native.complete.json')
            phase('integer_assignment_audit_running')
            # Import the installed native class only after CPU/env restrictions are active.
            import scape.apa_core
            def models():
                for dest in directories:
                    if (dest / 'res.gene.pkl').exists():
                        index = pd.read_csv(dest / 'barcode_index.csv').set_index('index').CB.to_dict()
                        assert index == dict(enumerate(barcodes))
                        yield from objects(dest / 'res.gene.pkl')
            features, mc, counts, audit = integer_assignments(models(), dict(enumerate(barcodes)), members)
            features.to_csv(out / 'site_inventory.tsv.gz', sep='\t', index=False,
                compression={'method': 'gzip', 'mtime': 0})
            pd.DataFrame(dict(metacell_id=mc)).to_csv(out / 'metacell_index.tsv', sep='\t', index=False)
            sparse.save_npz(out / 'metacell_by_site_counts.npz', counts)
            audit.update(time=now(), library_key=row['library_key'], input_cells=len(members),
                input_metacells=len(mc), discovered_sites=len(features),
                counted_sites=int((counts.sum(0).A1 > 0).sum()), integer_molecules=int(counts.sum()),
                members_sha256=row['members_sha256'], scientific_comparison_status='pending',
                count_unit='native merged PAS assignments; not globally deduplicated across genes')
            dump(audit, out / 'count_audit.json')
            phase('waiting_for_original_bam_full_sha256')
            while not (poly / 'complete.json').exists():
                time.sleep(30)
            original = json.loads((poly / 'complete.json').read_text())
            assert original['input_cells'] == len(members) and original['input_metacells'] == len(mc)
            assert source.stat().st_size == row['bam_bytes'] and source.stat().st_mtime_ns == row['bam_mtime_ns']
            assert len(original['bam_sha256']) == 64
            audit.update(original_bam_sha256=original['bam_sha256'],
                original_bam_sha256_provenance=str(poly / 'complete.json'))
            dump(audit, out / 'complete.json'); phase('native_and_integer_assignment_audit_completed')
            return audit
        except Exception as error:
            (out / 'failure.traceback.txt').write_text(traceback.format_exc())
            phase('failed', error=repr(error)); raise


def scaleup_plan(manifest, total_workers):
    original_workers = manifest['max_parallel_libraries']
    assert original_workers < total_workers
    assert total_workers * manifest['native_threads_per_library'] <= len(manifest['allowed_cpu_ids'])
    rows = round_robin(manifest['input_libraries'])
    return rows[original_workers:], total_workers - original_workers


def original_alive(saved):
    folder = Path('/proc') / str(saved['pid'])
    try:
        fields = (folder / 'stat').read_text().rsplit(') ', 1)[1].split()
        if int(fields[19]) != saved['start_ticks'] or fields[0] == 'Z':
            return False
        assert (folder / 'cwd').resolve() == ROOT
        command = (folder / 'cmdline').read_bytes().split(b'\0')
        assert str(Path(__file__)).encode() in command and b'--total-parallel-libraries' not in command
        return True
    except FileNotFoundError:
        return False


def scaleup(total_workers, original_pid, original_start_ticks):
    manifest = json.loads((OUT / 'input_manifest.json').read_text())
    rows, added_workers = scaleup_plan(manifest, total_workers)
    saved = dict(pid=original_pid, start_ticks=original_start_ticks)
    with (OUT / 'scaleup.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert original_alive(saved), 'Original runner identity must be verified before adding workers'
        assert int((OUT / 'runner.pid').read_text()) == original_pid
        assert not (OUT / 'scaleup_launch.json').exists(), 'Do not replace a previous scale-up attempt'
        launch = dict(time=now(), pid=os.getpid(), original_runner=saved,
            total_parallel_libraries=total_workers, original_parallel_libraries=manifest['max_parallel_libraries'],
            added_parallel_libraries=added_workers, added_queue_libraries=len(rows),
            input_manifest_sha256=digest(OUT / 'input_manifest.json'), scheduler_sha256=digest(Path(__file__)),
            allowed_cpu_ids=sorted(os.sched_getaffinity(0)), native_threads_per_library=2,
            original_native_jobs_restarted=False, scientific_comparison_status='pending')
        dump(launch, OUT / 'scaleup_launch.json')
        (OUT / 'runner_before_scaleup.pid').write_text(str(original_pid) + '\n')
        (OUT / 'scaleup.pid').write_text(str(os.getpid()) + '\n')
        (OUT / 'runner.pid').write_text(str(os.getpid()) + '\n')
        events('scape_atlas_scaleup_start', **{k: v for k, v in launch.items() if k != 'time'})
        completed, failures = [], []
        def status(phase):
            dump(dict(time=now(), pid=os.getpid(), phase=phase, original_runner=saved,
                total_parallel_libraries=total_workers, added_parallel_libraries=added_workers,
                completed_added_queue_libraries=len(completed), failed_added_queue_libraries=len(failures),
                failures=failures, scientific_comparison_status='pending'), OUT / 'scaleup_status.json')
        try:
            status('running')
            with ThreadPoolExecutor(max_workers=added_workers) as pool:
                tasks = {pool.submit(library, row, manifest, True): row for row in rows}
                for task in as_completed(tasks):
                    try:
                        completed.append(task.result())
                    except Exception as error:
                        failures.append(dict(library_key=tasks[task]['library_key'], error=repr(error)))
                    status('running')
            status('waiting_for_original_runner_to_close')
            while original_alive(saved):
                time.sleep(30)
            previous = json.loads((OUT / 'stage_status.json').read_text())
            dump(previous, OUT / 'original_runner_final_stage_status.json')
            complete = []
            for row in manifest['input_libraries']:
                path = OUT / row['study_id'] / row['library_id'] / 'complete.json'
                assert path.exists(), f'Full library completion missing: {row["library_key"]}'
                record = json.loads(path.read_text())
                assert record['input_cells'] == row['cells'] and record['input_metacells'] == row['metacells']
                assert record['integer_assignment_conservation'] and len(record['original_bam_sha256']) == 64
                complete.append(record)
            assert not failures, failures
            dump(dict(time=now(), phase='completed', total_libraries=len(complete),
                completed_libraries=len(complete), failed_libraries=0, failures=[],
                original_runner_attempt_records=previous.get('failures', []),
                original_runner_stage_provenance=str(OUT / 'original_runner_final_stage_status.json'),
                scientific_comparison_status='pending'), OUT / 'stage_status.json')
            dump(dict(time=now(), libraries=len(complete), status='all_native_jobs_counts_and_source_hashes_completed',
                scientific_comparison_status='pending'), OUT / 'native.complete.json')
            status('all_native_libraries_completed; science_pending')
        except Exception:
            (OUT / 'scaleup_failure.traceback.txt').write_text(traceback.format_exc())
            status('failed_requires_handling')
            raise


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((OUT / 'input_manifest.json').read_text())
    with (OUT / 'runner.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        (OUT / 'runner.pid').write_text(str(os.getpid()) + '\n')
        rows = round_robin(manifest['input_libraries'])
        assert len(rows) == 161 and sum(r['cells'] for r in rows) == 551201
        completed, failures = [], []
        dump(dict(time=now(), phase='running', total_libraries=len(rows), completed_libraries=0,
            failed_libraries=0, failures=[], scientific_comparison_status='pending'), OUT / 'stage_status.json')
        events('scape_atlas_start', pid=os.getpid(), libraries=len(rows), cells=manifest['cells'],
            scope=manifest['scope'], parallel_libraries=manifest['max_parallel_libraries'],
            allowed_logical_cpus=len(os.sched_getaffinity(0)))
        with ThreadPoolExecutor(max_workers=manifest['max_parallel_libraries']) as pool:
            tasks = {pool.submit(library, row, manifest): row for row in rows}
            for task in as_completed(tasks):
                row = tasks[task]
                try:
                    completed.append(task.result())
                except Exception as error:
                    failures.append(dict(library_key=row['library_key'], error=repr(error)))
                dump(dict(time=now(), phase='running', total_libraries=len(rows),
                    completed_libraries=len(completed), failed_libraries=len(failures), failures=failures,
                    scientific_comparison_status='pending'), OUT / 'stage_status.json')
        if not failures:
            dump(dict(time=now(), libraries=len(completed), status='all_native_jobs_counts_and_source_hashes_completed',
                scientific_comparison_status='pending'), OUT / 'native.complete.json')
        events('scape_atlas_end', completed=len(completed), failed=len(failures))
        raise SystemExit(bool(failures))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--total-parallel-libraries', type=int)
    parser.add_argument('--original-runner-pid', type=int)
    parser.add_argument('--original-runner-start-ticks', type=int)
    args = parser.parse_args()
    if args.total_parallel_libraries is None:
        main()
    else:
        assert args.original_runner_pid is not None and args.original_runner_start_ticks is not None
        scaleup(args.total_parallel_libraries, args.original_runner_pid, args.original_runner_start_ticks)
