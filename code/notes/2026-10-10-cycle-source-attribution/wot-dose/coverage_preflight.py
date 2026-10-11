"""Inspect availability only: never load APA outcomes or cycle predictors."""
from pathlib import Path
import hashlib
import json

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TABLE = ROOT / 'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/state_vectors/pair_effects.primary.tsv.gz'
ORTHO = ROOT / 'notes/2026-10-10-recurrence-utr-ortholog/measured_ortholog_members.tsv.gz'
RECURRENT = ROOT / 'notes/2026-10-10-directional-recurrence/observed_members.tsv.gz'
UNITS = ROOT / 'notes/2026-10-10-cpa-cycle-core-dose/table_audit_wot_units.tsv'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


units = pd.read_csv(UNITS, sep='\t')
assert len(units) == 24 and units.library_key.is_unique
orthologs = pd.read_csv(ORTHO, sep='\t').query('contrast == "P_M" and threshold == .02')
assert len(orthologs) == 751 and orthologs.mouse_gene_id.is_unique
core = orthologs.query('human_direction == "proximal" and mouse_direction == "proximal"')
assert len(core) == 27
proximal = pd.read_csv(RECURRENT, sep='\t').query('species == "mouse" and contrast == "P_M" and threshold == .02 and direction == "proximal"')
assert len(proximal) == 372 and proximal.gene_id.is_unique
columns = ['study_id', 'route', 'contrast', 'design', 'setting', 'library_key_a', 'library_key_b', 'time_a', 'time_b', 'gene_id', 'n_mc_a', 'n_mc_b', 'n_cells_a', 'n_cells_b', 'n_joint_PAS']
pieces = []
for part in pd.read_csv(TABLE, sep='\t', usecols=columns, chunksize=50000):
    keep = part.study_id.eq('GSE115943') & part.route.eq('WOT') & part.contrast.eq('P_M') & part.design.eq('same_library') & part.setting.eq('supported50_p08')
    pieces.append(part[keep])
presence = pd.concat(pieces, ignore_index=True)
assert set(presence.library_key_a) == set(units.library_key)
assert presence.library_key_a.eq(presence.library_key_b).all()
assert presence.time_a.eq(presence.time_b).all()
assert not presence.duplicated(['library_key_a', 'gene_id']).any()
assert presence[['n_mc_a', 'n_mc_b']].ge(3).all().all()
assert presence[['n_cells_a', 'n_cells_b']].ge(30).all().all()
assert presence.n_joint_PAS.ge(2).all()
presence = presence.rename(columns={'library_key_a': 'library_key'})
panels = {
    'shared27': set(core.mouse_gene_id),
    'mouse_proximal372': set(proximal.gene_id),
    'ortholog_background751': set(orthologs.mouse_gene_id),
    'native_measured_union': set(presence.gene_id),
}
rows, common_rows, summaries = [], [], {}
for name, genes in panels.items():
    gene_by_lib = {lib: set(presence[presence.library_key.eq(lib)].gene_id) & genes for lib in units.library_key}
    shared = sorted(set.intersection(*gene_by_lib.values()))
    for lib, subset in gene_by_lib.items():
        unit = units[units.library_key.eq(lib)].iloc[0]
        rows.append(dict(panel=name, n_fixed_members=len(genes), library_key=lib, time=unit.time, day=unit.day, C_label=unit.C_label, n_native_measured=len(subset)))
    for time, group in units.groupby('time', sort=False):
        selected = [gene_by_lib[lib] for lib in group.library_key]
        paired = sorted(set.intersection(*selected))
        common_rows.append(dict(panel=name, time=time, day=float(group.day.iloc[0]), n_wells=len(group), n_within_date_common=len(paired), gene_ids=json.dumps(paired)))
    summaries[name] = dict(n_fixed_members=len(genes), all24_n_common=len(shared), all24_common_gene_ids=shared, min_library_measured=min(map(len, gene_by_lib.values())), max_library_measured=max(map(len, gene_by_lib.values())))

pd.DataFrame(rows).to_csv(HERE / 'coverage_by_library.tsv', sep='\t', index=False)
pd.DataFrame(common_rows).sort_values(['panel', 'day']).to_csv(HERE / 'coverage_by_date.tsv', sep='\t', index=False)
presence.to_csv(HERE / 'native_coverage_presence.tsv.gz', sep='\t', index=False, compression={'method':'gzip', 'mtime':0})
out = dict(stage='coverage_only', new_APA_outcomes_loaded=False, new_cycle_associations_calculated=False, panels=summaries, n_libraries=len(units), n_dates=units.time.nunique(), n_two_well_dates=int(units.groupby('time').size().eq(2).sum()), n_native_gene_library_rows=len(presence), inputs={str(p.relative_to(ROOT)): digest(p) for p in [TABLE, ORTHO, RECURRENT, UNITS]})
(HERE / 'coverage_summary.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
