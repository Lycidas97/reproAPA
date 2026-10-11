"""Supplement plot-contract revision r2 (2026-10-06), written as a NEW interface package next to the sealed bundle.

Follows analysis_code/pas_rerun/upgrade_result_interface.py: the source bundle is verified and never modified; its tables
are hard-linked (read-only, byte-identical, SHA re-checked), new tables are added under their own prefixes, and the
contract, manuscript-contract tables, result records and completion marker are regenerated. No statistic is recomputed.

Changes (also written to contract-r2/CHANGES.tsv):
  1. bind the tables that the house-style supplement read outside the r1 interface to the panels that read them
  2. add the frozen GEX state-grouping tables (S4-S6) and the Trip12 external-evidence tables (S29) to the bundle
  3. add panels G1A, G2A, G3A (new S4-S6) and D26H (S29h, per-site evidence tiers)
  4. register a component for every newly bound table that has columns
  5. record the printed supplement numbering (S1-S30) and register the supdraw renderers
"""
import csv, gzip, json, os, shutil, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'Fig').is_dir() and (p/'main.tex').exists())
sys.path.insert(0, str(ROOT/'analysis_code/pas_rerun'))
from pipeline import digest, read, write, require             # noqa
from result_records import contained, export_records, source_identity   # noqa

SRC = ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/supplement'
OUT = HERE.parent/'contract-r2'/'bundle'
TASK = HERE.parent
REV = 'r2-20261006'

# ------------------------------------------------------------------ what changes
BIND = [  # (panel, table, reason)
    ('D01A', 'tables/frozen/marker_selection.json', 'S1 frozen marker list and gene groups (JSON)'),
    ('D02A', 'tables/frozen/marker_selection.json', 'S2 frozen marker list and gene groups (JSON)'),
    ('D01A', 'tables/frozen/palette.json', 'S1 state-family colours (Figure 1 palette)'),
    ('D02A', 'tables/frozen/palette.json', 'S2 state-family colours (Figure 1 palette)'),
    ('D03A', 'tables/frozen/palette.json', 'S3 state-family colours (Figure 1 palette)'),
    ('D07F', 'tables/panel_statistics/dominant-stage-marginal-expectation.tsv', 'S10 mouse dominant-stage agreement expectation'),
    ('D11E', 'tables/panel_statistics/example-cumulative-PAS-curves.tsv.gz', 'S14 PDLIM7 same-contributor cumulative PAS usage'),
    ('D17A', 'tables/state_vectors/route_effects.primary.tsv.gz', 'S20 per-route P-S DWUI of the shared P-S members'),
    ('D25B', 'tables/site_statistics/site_factor_inventory.tsv.gz', 'S28b Trip12 transferred-motif RBP sites (Figure 5f third rug row)'),
    ('D25B', 'tables/candidate_screen/recurrent_region_pairs.tsv.gz', 'S28b Trip12 shared strict P-M intervals'),
]
EXT_SG = json.loads((HERE/'external-sources.json').read_text())
EXT_T12 = json.loads((HERE/'external-sources-trip12.json').read_text())
ADD = []   # (key, origin path, new bundle path, node, panels, sha)
for ext, node in ((EXT_SG, 'state_grouping'), (EXT_T12, 'trip12_evidence')):
    for key, f in ext['files'].items():
        if 'bundle_table' in f:
            ADD.append((key, f['path'], f['bundle_table'], node, f['panels'], f['sha256']))
NEW_PANELS = {
    'G1A': ('S4', 'Human frozen states: marker evidence by study and partition stability', ['tables/frozen/states.tsv', 'tables/frozen/marker_by_state.tsv', 'tables/frozen/palette.json']),
    'G2A': ('S5', 'Mouse frozen states: marker evidence by study and partition stability', ['tables/frozen/states.tsv', 'tables/frozen/marker_by_state.tsv', 'tables/frozen/palette.json']),
    'G3A': ('S6', 'Partition stability of 11 original partitions and boundary-state programme contrasts', ['tables/frozen/states.tsv']),
    'D26H': ('S29', 'Trip12: every predicted site in the two strict intervals by evidence tier (Figure 5f rug source)', []),
}
ORDER = ['D01', 'D02', 'D03', 'G1', 'G2', 'G3']+[f'D{i:02d}' for i in range(4, 28)]
SNUM = {fid: f'S{i+1}' for i, fid in enumerate(ORDER)}
SCRIPT = {fid: ('figures/D_markers.py' if fid in ('D01', 'D02') else 'figures/G_markers.py' if fid in ('G1', 'G2') else f'figures/{fid}.py') for fid in ORDER}


def header(path):
    if path.suffix not in ('.tsv', '.gz'): return []
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return next(csv.reader(f, delimiter='\t'), [])


def main():
    require(not OUT.exists(), 'Use a new interface package directory')
    marker = read(SRC/'bundle-complete.json'); require(marker['status'] == 'statistics_bundle_complete', 'Incomplete source bundle')
    for name, key in [('plot-contract.json', 'plot_contract_sha256'), ('table-manifest.json', 'table_manifest_sha256'),
                      ('result-records.json', 'result_records_sha256')]:
        require(digest(SRC/name) == marker[key], 'Unverified source bundle: '+name)
    contract = read(SRC/'plot-contract.json'); identity = source_identity(contract)
    manifest = read(SRC/'table-manifest.json'); tables = manifest['tables']
    by_path = {r['path']: r for r in tables}
    OUT.mkdir(parents=True)
    # 1. hard-link every source table after re-checking its SHA
    for r in tables:
        src = contained(SRC, r['path']); require(digest(src) == r['sha256'], 'Changed source table: '+r['path'])
        dst = contained(OUT, r['path']); dst.parent.mkdir(parents=True, exist_ok=True)
        os.link(src, dst)
    # 2. new tables
    changes = []
    new_rows = []
    for key, origin, newp, node, panels, sha in ADD:
        src = ROOT/origin; require(digest(src) == sha, f'{key}: origin SHA differs from the pinned value')
        require(newp not in by_path, 'New table path collides: '+newp)
        dst = contained(OUT, newp); dst.parent.mkdir(parents=True, exist_ok=True); os.link(src, dst)
        row = dict(node=node, path=newp, sha256=sha, bytes=src.stat().st_size, columns=header(src), role='source_data',
                   historical_fields_present=False, origin=origin, added_in=REV)
        tables.append(row); by_path[newp] = row; new_rows.append(row)
        changes.append(dict(change='add_table', panel=';'.join(panels), table=newp, detail=f'{key} from {origin} (SHA pinned)'))
    # 3. panels
    panels = contract['panels']; P = {p['panel']: p for p in panels}
    for pid, (fig, title, base) in NEW_PANELS.items():
        require(pid not in P, 'Panel exists: '+pid)
        p = dict(panel=pid, figure=fig, title=title, nodes=[], scope='render_bound_tables_only', tables=list(base),
                 required_tables=[], added_in=REV)
        panels.append(p); P[pid] = p
        changes.append(dict(change='add_panel', panel=pid, table='', detail=f'{fig}: {title}'))
    for pid, t, why in BIND:
        require(t in by_path, 'Unknown table '+t)
        if t not in P[pid]['tables']:
            P[pid]['tables'].append(t); changes.append(dict(change='bind_table', panel=pid, table=t, detail=why))
    for key, origin, newp, node, pids, sha in ADD:
        for pid in pids:
            if newp not in P[pid]['tables']: P[pid]['tables'].append(newp); changes.append(dict(change='bind_table', panel=pid, table=newp, detail=key))
    # printed numbering: legacy figure label kept, printed number added
    for p in panels:
        d = p['panel'][:-1]
        if d in SNUM:
            p['printed_figure'] = SNUM[d]
    contract_r2 = dict(contract, panels=panels, contract_revision=REV,
                       supplement_numbering=dict(order=ORDER, printed={k: v for k, v in SNUM.items()},
                                                 note='S1-S3 = D01-D03, S4-S6 = G1-G3 (state grouping, added 2026-10-06), S7-S30 = D04-D27; '
                                                      'panel field "figure" keeps the r1 label, "printed_figure" is the current number'))
    write(OUT/'plot-contract.json', contract_r2)
    write(OUT/'table-manifest.json', dict(manifest, tables=tables))
    # 4. manuscript contracts + components
    shutil.copytree(SRC/'manuscript-contracts', OUT/'manuscript-contracts')
    mc = OUT/'manuscript-contracts'
    os.chmod(mc/'supplement-data-interface.tsv', 0o600); os.chmod(mc/'figure-renderers.json', 0o600)
    with (mc/'supplement-data-interface.tsv').open() as f: iface = list(csv.DictReader(f, delimiter='\t'))
    have = {r['component_id'] for r in iface}
    with (mc/'panel-data-interface.tsv').open() as f: have |= {r['component_id'] for r in csv.DictReader(f, delimiter='\t')}
    adds = [(pid, t) for pid, t, _ in BIND]+[(pid, newp) for _, _, newp, _, pids, _ in ADD for pid in pids]
    adds += [(pid, t) for pid, (_, _, base) in NEW_PANELS.items() for t in base]
    for pid, t in adds:
        cols = by_path[t]['columns']
        cid = f'{pid}:{t[7:]}'
        if not cols or cid in have: continue
        iface.append(dict(component_id=cid, panel=pid, table=t[7:], read_columns=';'.join(cols),
                          selection=f'Contract {REV}: rows as filtered by the supdraw renderer; denominators retained',
                          unit_or_uncertainty='Use recorded denominators, measurement status and original uncertainty; no old-run fallback'))
        have.add(cid); changes.append(dict(change='add_component', panel=pid, table=t, detail=cid))
    with (mc/'supplement-data-interface.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(iface[0]), delimiter='\t'); w.writeheader(); w.writerows(iface)
    reg = read(mc/'figure-renderers.json')
    for fid in ORDER:
        pids = sorted(p['panel'] for p in panels if p['panel'][:-1] == fid)
        reg['figures'].append(dict(figure=SNUM[fid], design_figure=fid, script='notes/2026-10-06-supplement-polyapipe/'+SCRIPT[fid],
                                   sha256=digest(TASK/SCRIPT[fid]), panels=pids, interface='result_records_supdraw',
                                   status='candidate_rendered', added_in=REV))
    write(mc/'figure-renderers.json', reg)
    changes.append(dict(change='register_renderers', panel='', table='', detail=f'{len(ORDER)} supplement renderers (S1-S30)'))
    record = export_records(OUT, identity, panels, tables, mc)
    parent = dict(source_bundle=str(SRC.relative_to(ROOT)), source_marker_sha256=digest(SRC/'bundle-complete.json'),
                  source_table_manifest_sha256=digest(SRC/'table-manifest.json'), source_plot_contract_sha256=digest(SRC/'plot-contract.json'),
                  source_result_records_sha256=digest(SRC/'result-records.json'), tables_hard_linked=True,
                  source_statistics_recomputed=False, production_replaced=False, contract_revision=REV)
    write(OUT/'interface-parent.json', parent)
    write(OUT/'bundle-complete.json', dict(marker, plot_contract_sha256=digest(OUT/'plot-contract.json'),
          table_manifest_sha256=digest(OUT/'table-manifest.json'), result_records_sha256=digest(OUT/'result-records.json'),
          statistics_schema_id=record['statistics_schema_id'], interface_parent_sha256=digest(OUT/'interface-parent.json'),
          interface_only_upgrade=True, contract_revision=REV, supplement_panels=marker['supplement_panels']+len(NEW_PANELS),
          new_figures_rendered=False, production_replaced=False))
    with (OUT.parent/'CHANGES.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=['change', 'panel', 'table', 'detail'], delimiter='\t'); w.writeheader(); w.writerows(changes)
    print(OUT, len(tables), 'tables;', len(record['components']), 'components;', len(changes), 'changes')


if __name__ == '__main__':
    main()
