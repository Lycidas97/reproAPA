"""Pin and package the promoted v14 manuscript and its code using the existing packager."""
import importlib.util
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[3]
RUN = Path(__file__).resolve().parent
REPO = ROOT / 'github/reproapa-analysis'
spec = importlib.util.spec_from_file_location('packager', REPO / 'scripts/package_source_data.py')
pack = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pack)


def add(rows, path, role, archive_path=None):
    path = Path(path)
    relative = path.relative_to(ROOT).as_posix()
    row = rows.setdefault(relative, dict(path=relative,
        archive_path=archive_path or 'project/' + relative,
        bytes=path.stat().st_size, sha256=pack.digest(path), roles=[], owners=['v14 production']))
    if role not in row['roles']:
        row['roles'].append(role)


def binding(path):
    with zipfile.ZipFile(path) as z:
        members = len(z.namelist())
    return dict(path=path.relative_to(ROOT).as_posix(), bytes=path.stat().st_size,
                sha256=pack.digest(path), members=members)


def save(selection, name, output):
    selection['files'] = sorted(selection['files'], key=lambda x: x['path'])
    file = RUN / name
    file.write_text(json.dumps(selection, indent=2) + '\n')
    pack.package(ROOT, file, output)


base = json.loads((REPO / 'provenance/source-data-selection.json').read_text())
candidate = json.loads((REPO / 'provenance/draft-source-data-selection.json').read_text())
inventory = json.loads((REPO / 'provenance/source-selection.json').read_text())
limits = ['Scientific analyses and model fitting were not rerun during packaging.',
          'Raw reads, full training matrices and scientific runtime resources remain external.',
          'Historical rendering receipts preserve the identities of their original runs.']

# Current workspace implementations supersede the earlier captured code where available.
code = {}
for row in inventory['files']:
    if row['source_alias'] == 'PAPER':
        source = ROOT / row['source_path']
        archive_path = 'project/' + row['source_path']
    else:
        source = REPO / row['repository_path']
        archive_path = 'captured/' + row['repository_path']
    add(code, source, 'analysis_code_or_configuration', archive_path)
for p in (ROOT / 'revision-backups/20261011-before-v11-apply/supplement').rglob('*.tex'):
    add(code, p, 'pinned_pre_revision_generator_input')
for p in (ROOT / 'revision-backups/20261011-before-v11-apply/Fig').glob('*.pdf'):
    add(code, p, 'pinned_pre_revision_generator_input')
for p in [ROOT / 'reviews/2026-10-10-divergence-writing/work/supp_methods_v11.py',
          ROOT / 'github/reproapa-analysis/scripts/package_source_data.py', Path(__file__).resolve()]:
    add(code, p, 'packaging_or_text_generator')
code_output = ROOT / 'table/Supplementary_Source_Code_v11.zip'
code_selection = dict(schema='reproapa-source-data/1', date='2026-10-11',
    scope='Current implementations and pinned inputs for the v14 production revision.',
    files=list(code.values()), companion=base['companion'],
    coverage_counts=dict(source_files=len(code)), figure_provenance_limits=limits,
    archive_readme='''# Source code for the 11 October 2026 production revision

SOURCE-MANIFEST.tsv pins every implementation and configuration. Restore project/
into the manuscript workspace. captured/ retains the separately captured upstream
implementations and environment records under their repository layout.

Current workspace implementations supersede the earlier captured version where
available. The frozen pre-revision figure PDFs and supplementary sources allow
the renumbering generators to run after production promotion.

The unchanged Supplementary_Source_Code_S9.zip is a required companion, bound
by COMPANION.json. This extension does not replace its historical inputs.
Full scientific input matrices and runtime resources remain external.
No scientific analyses were rerun while packaging.
''')
if not code_output.exists():
    save(code_selection, 'code-selection.json', code_output)
else:
    receipt = json.loads(code_output.with_suffix('.verification.json').read_text())
    assert receipt['sha256'] == pack.digest(code_output), 'Code archive changed'

data = {}
for row in candidate['files']:
    add(data, ROOT / row['path'], 'revision_source_or_render_record')
for row in candidate['required_primary_sources']:
    add(data, ROOT / row['path'], 'reused_renderer_input')
for name in ('unique_support_bridge_route_summary.tsv', 'unique_support_route_summary.tsv',
             'gradient_route_summary.tsv', 'classification_confusion.tsv'):
    add(data, ROOT / 'notes/2026-10-10-marker-reference-sensitivity' / name, 'S10a_source')
for root in ('main', 'supplement'):
    for p in (ROOT / root).rglob('*.tex'):
        add(data, p, 'production_manuscript')
    add(data, ROOT / root / '.texsplit.json', 'production_manifest')
for p in (ROOT / 'Fig').glob('*.pdf'):
    add(data, p, 'production_figure')
for p in [ROOT / 'Fig/figure-manifest.tsv', ROOT / 'references.bib',
          ROOT / 'table/Supplementary_Tables_S1-S10.xlsx',
          ROOT / 'reviews/2026-10-10-divergence-writing/spec/logic.json',
          ROOT / 'reviews/2026-10-10-divergence-writing/spec/evidence.json',
          code_output.with_suffix('.manifest.tsv'), code_output.with_suffix('.verification.json')]:
    add(data, p, 'production_output_or_provenance')
data_selection = dict(schema='reproapa-source-data/1', date='2026-10-11',
    scope='Promoted logic-v14 production: four main figures, 21 supplementary figures and Tables S1-S10.',
    files=list(data.values()), companion=binding(ROOT / 'table/Supplementary-Source-Data-20261009.zip'),
    additional_companions=[binding(code_output), base['companion']],
    coverage_counts=dict(main_figures=4, supplementary_figures=21, tables=10,
                         data_sheets=57, total_sheets=59, blocks=197, new_table_source_links=28),
    figure_provenance_limits=limits,
    archive_readme='''# Source data for the 11 October 2026 production revision

This archive records the promoted v14 manuscript, four main figures, 21
supplementary figures and Tables S1-S10 (57 data sheets, 197 blocks).
SOURCE-MANIFEST.tsv pins every file; restore project/ into the workspace.
The Sources worksheet records each block's input, filter and SHA-256 identity.

The unchanged 9 October source-data archive supplies retained historical inputs
and is required. COMPANION.json binds that identity. The additional_companions
in SOURCE-SELECTION.json bind Supplementary_Source_Code_v11.zip and the original
Supplementary_Source_Code_S9.zip. Supply these archives together for reuse.

Candidate rendering receipts are retained as historical records. Production
file identities are given separately by Fig/figure-manifest.tsv. No scientific
analyses or historical figure rendering were rerun while packaging. Raw reads,
full training matrices and scientific runtime resources remain external.
''')
save(data_selection, 'source-selection.json', ROOT / 'table/Supplementary-Source-Data-20261011-v11.zip')
