"""Shared identities, receipts and durable progress for the seven validation jobs."""
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECTS = Path('/root/nfsdata/apadetect/analysis_projects')
ROOT = PROJECTS / 'reproapa_integration_validation_20260913'
RUN = ROOT / 'execution_v1'
JOINT = PROJECTS / 'human_mouse_time_aware_umap_v1'
BASE = JOINT / 'results_raw_cell_species_integrated_v3_faiss2'
MOUSE = PROJECTS / 'mouse_time_aware_integration_pilot_v1'
HUMAN = PROJECTS / 'human_time_aware_integration_v1'
QC = PROJECTS / 'reprogramming_unified_rna_qc_v3.2_release'
SEED = 20260913


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    os.replace(tmp, path)


def progress(job, phase, **info):
    record = dict(job=job, phase=phase, pid=os.getpid(), updated_unix=time.time(),
                  full_analysis_complete=False, **info)
    write_json(RUN / job / 'status.json', record)
    print(json.dumps(record, ensure_ascii=False), flush=True)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def receipt(paths, output):
    write_json(output, [dict(path=str(p), bytes=Path(p).stat().st_size, sha256=sha(p)) for p in paths])


def obs():
    x = pd.read_parquet(RUN / 'shared/atlas.parquet')
    assert len(x) == 551201 and x.cell_id.is_unique
    return x


def wait_file(job, path, dependency):
    last = 0
    while not Path(path).exists():
        if time.time() - last > 300:
            progress(job, 'WAITING_INPUT', dependency=dependency, required=str(path))
            last = time.time()
        status = RUN / dependency / 'status.json'
        if status.exists() and json.loads(status.read_text())['phase'] == 'FAILED':
            raise RuntimeError(f'Upstream {dependency} failed; see {status}')
        time.sleep(20)


def finish(job, stage, **info):
    progress(job, 'STAGE_COMPLETE', stage=stage, **info)
    write_json(RUN / job / (stage + '.complete.json'), dict(stage=stage, **info))


def execute(job, fn):
    import traceback
    (RUN / job).mkdir(parents=True, exist_ok=True)
    progress(job, 'RUNNING', stage='initializing')
    try:
        fn()
    except BaseException as exc:
        progress(job, 'FAILED', error=str(exc), traceback=traceback.format_exc())
        raise
