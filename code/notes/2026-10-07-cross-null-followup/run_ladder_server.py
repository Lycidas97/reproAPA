"""Rerun masked ladder calibration from frozen effects, without native matrix extraction."""
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
available = os.sched_getaffinity(0)
os.sched_setaffinity(0, {2 if 2 in available else min(available)})
root = Path(__file__).resolve().parent
sys.path.insert(0, str(root / 'code'))
import numpy as np
import pandas as pd
import scipy
from downstream_runtime import Context
from downstream_cross import controls, feature_frame

prior = root.parent / 'reproapa_cross_overlap_null_a1_20261007_v1'
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
status = dict(status='running', sources={})
(root / 'environment.json').write_text(json.dumps(dict(python=platform.python_version(),
    numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__), indent=2) + '\n')

try:
    for source in ['polyapipe', 'apadetect-p50']:
        print('START', source, flush=True)
        folder = root / source
        folder.mkdir(exist_ok=True)
        config = json.loads((prior / source / 'config.json').read_text())
        (folder / 'source.config.json').write_text(json.dumps(config, indent=2) + '\n')
        feature_ctx = Context(folder / 'source.config.json', 'cross_overlap_features', folder / 'analysis/cross_overlap')
        feature_ctx.register(root / 'code/downstream_cross.py')
        original = feature_ctx.read(feature_ctx.analysis / 'cross_overlap/ortholog_matched_features.tsv.gz')
        orth = feature_ctx.read(feature_ctx.resource('ortholog'))
        nodes = feature_ctx.read(feature_ctx.analysis / 'cross_overlap/all_mapped_nodes.tsv.gz')
        paired = feature_ctx.read(feature_ctx.analysis / 'cross_overlap/paired_pm_study_effects.tsv.gz')
        frame = feature_frame(orth, nodes, paired, 'supported50')
        pd.testing.assert_frame_equal(frame[original.columns], original, check_dtype=False, rtol=0, atol=1e-12)
        feature_ctx.table(frame, 'ortholog_matched_features.tsv.gz')
        feature_ctx.finish(dict(original_features_and_labels_preserved=True,
            actual_effect_counts_rebuilt_from_frozen_effects=True))
        config['analysis_root'] = str(folder / 'analysis')
        for name in ['common.py', 'cross_species.py']:
            path = root / 'code/numerics' / name
            config['source_overrides']['notes/2026-10-02-null-model-controls/' + name] = dict(path=str(path), sha256=digest(path))
        config['code_sha256'].update({str(path): digest(path) for path in (root / 'code').glob('*.py')})
        config['ladder_rerun'] = dict(source_config_sha256=digest(folder / 'source.config.json'),
            correction='Actual measurable effect count added at every ladder step and in L3 injection')
        (folder / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
        status['sources'][source] = 'running'
        (root / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
        ctx = Context(folder / 'config.json', 'cross_controls', folder / 'ladder')
        for path in [Path(__file__), root / 'code/downstream_cross.py', root / 'code/downstream_runtime.py']:
            ctx.register(path)
        controls(ctx)
        status['sources'][source] = 'complete'
        (root / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
        print('DONE', source, flush=True)
    status['status'] = 'complete'
except Exception as error:
    status.update(status='failed', error=repr(error))
    raise
finally:
    (root / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
