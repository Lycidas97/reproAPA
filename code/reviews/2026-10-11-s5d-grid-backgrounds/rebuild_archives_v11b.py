"""Package the v11b code and source-data archives after Table S5d gained blocks S5d.6-S5d.8 (issue I138).
Same packager and selection rules as production-20261011/rebuild_archives.py; differences are marked 'v11b'."""
import importlib.util
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[2]
RUN = Path(__file__).resolve().parent
REPO = ROOT / 'github/reproapa-analysis'
spec = importlib.util.spec_from_file_location('packager', REPO / 'scripts/package_source_data.py')
pack = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pack)
PRIOR = ROOT / 'reviews/2026-10-10-divergence-writing/production-20261011'


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
# v11b: provenance/draft-source-data-selection.json is no longer in the repository; the revision records and reused
# renderer inputs are taken from the v11a manifest, which pins the same selection, and are re-hashed here.
import csv
with open(ROOT / 'table/Supplementary-Source-Data-20261011-v11a.manifest.tsv', newline='') as fh:
    PREV = list(csv.DictReader(fh, delimiter='\t'))
inventory = json.loads((REPO / 'provenance/source-selection.json').read_text())
limits = ['Scientific analyses and model fitting were not rerun during packaging.',
          'Raw reads, full training matrices and scientific runtime resources remain external.',
          'Historical rendering receipts preserve the identities of their original runs.']
S5D = [ROOT / 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables/motif_controls' / n
       for n in ('grid_summary.tsv', 'background_ladder.tsv', 'motif_power.tsv')]

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
          REPO / 'scripts/package_source_data.py', PRIOR / 'rebuild_archives.py', Path(__file__).resolve()]:
    add(code, p, 'packaging_or_text_generator')
code_output = ROOT / 'table/Supplementary_Source_Code_v11b.zip'                          # v11b
code_selection = dict(schema='reproapa-source-data/1', date='2026-10-11',
    scope='Current implementations and pinned inputs for the v14 production revision (v11b).',
    files=list(code.values()), companion=base['companion'],
    coverage_counts=dict(source_files=len(code)), figure_provenance_limits=limits,
    archive_readme='''# Source code for the 11 October 2026 production revision (v11b)

SOURCE-MANIFEST.tsv pins every implementation and configuration. Restore project/
into the manuscript workspace. captured/ retains the separately captured upstream
implementations and environment records under their repository layout.

Current workspace implementations supersede the earlier captured version where
available. The frozen pre-revision figure PDFs and supplementary sources allow
the renumbering generators to run after production promotion.

This v11b snapshot supersedes Supplementary_Source_Code_v11.zip. The workbook
generator (notes/2026-10-10-supplement-v11/tables) now appends the matching grid,
random-priority backgrounds and detection power to Table S5d (blocks S5d.6-S5d.8),
and the related-figure letters of Tables S2f and S3h follow the relettered Figure 2.

The unchanged Supplementary_Source_Code_S9.zip is a required companion, bound
by COMPANION.json. This extension does not replace its historical inputs.
Full scientific input matrices and runtime resources remain external.
No scientific analyses were rerun while packaging.
''')
save(code_selection, 'code-selection.json', code_output)

data = {}
for role in ('revision_source_or_render_record', 'reused_renderer_input'):
    for row in PREV:
        if role in row['roles'].split('; '):
            add(data, ROOT / row['path'], role)
for name in ('unique_support_bridge_route_summary.tsv', 'unique_support_route_summary.tsv',
             'gradient_route_summary.tsv', 'classification_confusion.tsv'):
    add(data, ROOT / 'notes/2026-10-10-marker-reference-sensitivity' / name, 'S10a_source')
for p in S5D:                                                                             # v11b
    add(data, p, 'S5d_source')
for root in ('main', 'supplement'):
    for p in (ROOT / root).rglob('*.tex'):
        add(data, p, 'production_manuscript')
    add(data, ROOT / root / '.texsplit.json', 'production_manifest')
for p in (ROOT / 'Fig').glob('*.pdf'):
    add(data, p, 'production_figure')
for p in [ROOT / 'Fig/figure-manifest.tsv', ROOT / 'references.bib',
          ROOT / 'table/Supplementary_Tables_S1-S10.xlsx',
          ROOT / 'notes/2026-10-10-supplement-v11/out/Supplementary_Tables_S1-S10.v11.build.json',
          ROOT / 'notes/2026-10-10-supplement-v11/out/Supplementary_Tables_S1-S10.v11.verify.json',
          ROOT / 'reviews/2026-10-10-divergence-writing/spec/logic.json',
          ROOT / 'reviews/2026-10-10-divergence-writing/spec/evidence.json',
          code_output.with_suffix('.manifest.tsv'), code_output.with_suffix('.verification.json')]:
    add(data, p, 'production_output_or_provenance')
data_selection = dict(schema='reproapa-source-data/1', date='2026-10-11',
    scope='Promoted logic-v14 production: four main figures, 21 supplementary figures and Tables S1-S10 (v11b).',
    files=list(data.values()), companion=binding(ROOT / 'table/Supplementary-Source-Data-20261009.zip'),
    additional_companions=[binding(code_output), base['companion']],
    coverage_counts=dict(main_figures=4, supplementary_figures=21, tables=10,
                         data_sheets=57, total_sheets=59, blocks=200, new_table_source_links=28, appended_S5d_block_sources=3),
    figure_provenance_limits=limits,
    archive_readme='''# Source data for the 11 October 2026 production revision (v11b)

This archive records the promoted v14 manuscript, four main figures, 21
supplementary figures and Tables S1-S10 (57 data sheets, 200 blocks).
SOURCE-MANIFEST.tsv pins every file; restore project/ into the workspace.
The Sources worksheet records each block's input, filter and SHA-256 identity.

The unchanged 9 October source-data archive supplies retained historical inputs
and is required. COMPANION.json binds that identity. The additional_companions
in SOURCE-SELECTION.json bind Supplementary_Source_Code_v11b.zip and the original
Supplementary_Source_Code_S9.zip. Supply these archives together for reuse.

Candidate rendering receipts are retained as historical records. Production
file identities are given separately by Fig/figure-manifest.tsv. No scientific
analyses or historical figure rendering were rerun while packaging. Raw reads,
full training matrices and scientific runtime resources remain external.

This v11b snapshot supersedes v11a. Table S5d gains blocks S5d.6-S5d.8: the
12-setting matching grid, the random-priority same-gene and cross-gene
backgrounds, and the detection power of the matched-window test. Their three
source tables are archived here and are unchanged from the 9 October archive.
The related-figure letters of Tables S2f and S3h follow the relettered Figure 2.
Since v11a, the manuscript sources also carry the 11 October sentence-level
revisions, the revised main-figure captions and the updated archive names, and
Figures 1-4 carry the 11 October in-place label, lettering and color edits
recorded in Fig/figure-manifest.tsv. No supplementary figure changed.
''')
save(data_selection, 'source-selection.json', ROOT / 'table/Supplementary-Source-Data-20261011-v11b.zip')
