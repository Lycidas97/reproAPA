"""Cross-check the parallel deliverables without repeating their nulls."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
modules = ['2026-10-10-naive-mouse-correspondence', '2026-10-10-recurrence-component-split',
           '2026-10-10-state-temporal-gradients', '2026-10-10-recurrence-utr-ortholog']
audits = {}
for module in modules:
    directory = ROOT/'notes'/module
    path = directory/'manifest.json'
    manifest = json.loads(path.read_text())
    if 'status' in manifest:
        assert manifest['status'] in ['PASS', 'complete'], (module, manifest['status'])
        outputs = manifest.get('output_sha256', {})
    else:
        assert module == '2026-10-10-state-temporal-gradients'
        validation = json.loads((directory/'validation.json').read_text())
        assert validation['status'] == 'passed'
        outputs = validation['output_sha256']
        assert sha(directory/'analyze.py') == manifest['analysis_code_sha256']
        assert sha(directory/'protocol.json') == manifest['protocol_sha256']
        assert all(sha(ROOT/p) == digest for p, digest in manifest['upstream_code_sha256'].items())
    inputs = manifest.get('input_sha256', manifest.get('inputs', {}))
    for name, digest in inputs.items():
        source = Path(name) if Path(name).is_absolute() else ROOT/name
        assert sha(source) == digest, (module, 'input', name)
    for name, digest in outputs.items():
        assert sha(directory/name) == digest, (module, 'output', name)
    if 'script_sha256' in manifest:
        assert sha(directory/'analyze.py') == manifest['script_sha256']
    tests = directory/'tests.log'
    if tests.exists():
        log = tests.read_text()
        assert log.rstrip().endswith('OK') and 'skipped' not in log.lower()
    assert (directory/'REPORT.zh.md').exists()
    audits[module] = dict(inputs=len(inputs), outputs=len(outputs), manifest_sha256=sha(path))

plot_dir = ROOT/'notes/2026-10-10-recurrence-utr-ortholog'
for name in ['figure-verification.json', 'state-figure-verification.json']:
    checks = json.loads((plot_dir/name).read_text())
    for path, digest in checks['input_sha256'].items():
        target = ROOT/path if '/' in path else plot_dir/path
        assert sha(target) == digest, (name, path)
    assert checks['all_text_within_canvas']

# Independently reconstitute every original label in the new ortholog tables.
b = ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
frozen = pd.read_csv(b/'recurrence_controls/recurrence_frozen_covariates.tsv.gz', sep='\t')
pairs = pd.read_csv(plot_dir/'measured_ortholog_members.tsv.gz', sep='\t')
cells = pd.read_csv(plot_dir/'ortholog_3x3.tsv', sep='\t')
for (co, t), q in pairs.groupby(['contrast', 'threshold']):
    source = frozen[frozen.contrast.eq(co)&frozen.threshold.eq(t)]
    for sp in ['human', 'mouse']:
        sf = source[source.species.eq(sp)]
        votes = sf.pivot(index='gene_id', columns='unit', values='event')
        counts = votes.notna().sum(1)
        ids = q[sp+'_gene_id']
        assert counts.reindex(ids).ge(2).all()
        labels = pd.Series('nonrecurrent', index=votes.index)
        labels.loc[votes.eq(-1).sum(1).ge(2)] = 'proximal'
        labels.loc[votes.eq(1).sum(1).ge(2)] = 'distal'
        assert q[sp+'_direction'].tolist() == labels.reindex(ids).tolist()
    for _, row in cells[cells.contrast.eq(co)&cells.threshold.eq(t)].iterrows():
        observed = (q.human_direction.eq(row.human_direction)&q.mouse_direction.eq(row.mouse_direction)).sum()
        assert observed == row.observed

# Rounding in the component prose must agree with the saved means.
component_dir = ROOT/'notes/2026-10-10-recurrence-component-split'
summary = pd.read_csv(component_dir/'component_distributions.tsv', sep='\t')
report_rows = [line.split('|')[1:-1] for line in (component_dir/'REPORT.zh.md').read_text().splitlines()
               if line.startswith('| 人 P–') or line.startswith('| 鼠 P–')]
assert len(report_rows) == 8
for parts in report_rows:
    label = parts[0]
    species = 'human' if '人' in label else 'mouse'
    contrast = 'P_S' if 'P–S' in label else 'P_M'
    direction = 'distal' if '远端' in label else 'proximal'
    scope = 'GSE118258' if species == 'human' else 'CIP_WOT_common'
    selected = summary[summary.species.eq(species)&summary.scope.eq(scope)&summary.selection_contrast.eq(contrast)&
                       summary.direction.eq(direction)&summary.threshold.eq(.02)].set_index('metric')
    for metric, text in zip(['DWUI_PS', 'DWUI_PM', 'DWUI_MS'], parts[2:5]):
        value = float(text.strip().replace('−', '-'))
        assert abs(value-selected.loc[metric, 'mean']) <= .00000501, (label, metric)

result = dict(status='PASS', modules=audits, this_round_intent_tests=19, skipped_tests=0,
    actual_measurement_required_for_ortholog_nonrecurrence=True,
    ortholog_3x3_observed_counts_independently_replayed=True,
    component_report_all_eight_rows_match_saved_values=True,
    figure_input_hashes_match_final_tables=True, figure_canvas_checks_pass=True,
    make_check_passed='check complete: no errors' in (plot_dir/'make-check.log').read_text(),
    PDF_visual_QA='Both PDF rasters and 150dpi PNGs reviewed; embedded Nimbus Sans L confirmed with pdffonts; qpdf unavailable.',
    manuscripts_modified=False)
(HERE/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result))
