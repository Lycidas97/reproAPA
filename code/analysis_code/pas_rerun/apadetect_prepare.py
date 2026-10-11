"""Receive completed APADETECT counts on the frozen expression framework.

No BAM access, discovery, statistical inference, or production-figure replacement.
The full PAS axis is retained; unavailable windows never become measured zeros.
"""
import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from aggregate import aggregate_counts
from polyapipe_prepare import digest, save
from replace_pas_reference import canonical, support_fields


def require(value, message):
    if not value:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def checked(path, inputs, expected=None):
    path = Path(path).resolve()
    actual = digest(path)
    require(expected is None or actual == expected, f'Changed input: {path}')
    require(str(path) not in inputs or inputs[str(path)] == actual, f'Input changed during intake: {path}')
    inputs[str(path)] = actual
    return path


def receipt(folder, inputs):
    folder = Path(folder)
    record = read_json(checked(folder / 'complete.json', inputs))
    require(record['status'] == 'COMPLETE' and bool(record['artifacts']), f'Incomplete producer: {folder}')
    for name, sha in record['artifacts'].items():
        p = (folder / name).resolve()
        require(p.is_relative_to(folder.resolve()), f'Artifact escapes producer directory: {name}')
        checked(p, inputs, sha)
    return record


def truth(series):
    require(series.notna().all(), f'Missing Boolean: {series.name}')
    s = series.astype(str).str.lower()
    require(s.isin(['true', 'false']).all(), f'Invalid Boolean: {series.name}')
    return s.eq('true')


def apply_view(features, view):
    """Independent analysis route: same measured counts, fresh eligible PAS axis."""
    f = features.copy()
    supported = truth(f.coordinate_supported)
    quantified = truth(f.quantifiable)
    selected = np.ones(len(f), dtype=bool)
    if view != 'base':
        require(view.startswith('reject_original_') and view in f, 'Unknown ASNR analysis route')
        selected = f[view].notna().to_numpy() & ~truth(f[view].fillna(False)).to_numpy()
    f['analysis_view'] = view
    f['measurement_eligible'] = supported & quantified & selected
    f['measurement_exclusion'] = np.select(
        [~quantified, ~supported, ~selected],
        ['unavailable_count_window', 'unavailable_coordinate', 'ASNR_route_excluded_or_unknown'], default='')
    for col in ['keep_for_apa_primary', 'keep_for_sensitivity', 'keep_for_pas_quant']:
        structural = 'structural_' + col
        if structural not in f:
            f[structural] = truth(f[col])
        f[col] = truth(f[structural]) & f.measurement_eligible
    return f


def convert_features(raw, references, tiers, view='base'):
    """Use upstream genomic-base annotations, then export interbase cuts."""
    f = raw.copy().reset_index(drop=True)
    require(f.name.is_unique and f.strand.isin(['+', '-']).all(), 'Invalid PAS identity/strand')
    supported = truth(f.coordinate_supported)
    quantified = truth(f.quantifiable)
    require(np.array_equal(f.loc[supported, 'cut_site'], f.loc[supported, 'score']), 'Upstream coordinate mismatch')
    require(set(tiers.pas_id) == set(f.pas_id) and tiers.pas_id.is_unique, 'Confidence-tier axis mismatch')
    # Annotation and sequence QC were performed at the original base coordinate.
    # Keep those values distinct from the interbase analysis interface.
    f['native_cut_site_0based_base'] = f.score.astype(np.int64)
    f['annotation_cut_site_0based_base'] = f.cut_site
    f['cut_site'] = f.score.astype(np.int64) + f.strand.eq('+').astype(np.int64)
    require(np.array_equal(f.loc[supported, 'terminal_cut_boundary0'], f.loc[supported, 'cut_site']),
            'Independent terminal-protection coordinate conversion disagrees')
    f['atlas_pas_id'] = f.name
    f['atlas_column_0based'] = np.arange(len(f))
    f['feature_index_1based'] = np.arange(1, len(f) + 1)
    for source, target in [('min_cut', 'member_cut_site_min'), ('max_cut', 'member_cut_site_max')]:
        f[target] = f[source].astype(np.int64) + f.strand.eq('+').astype(np.int64)
    f['cluster_width_bp'] = f.max_cut - f.min_cut
    f['native_peaks'] = f.member_count
    f['n_libraries'] = f.library_count
    f['atlas_status'] = np.where(f.library_count.ge(2), 'core', 'single_library')
    f = f.drop(columns=[c for c in tiers if c != 'pas_id' and c in f])
    f = f.merge(tiers, on='pas_id', how='left', sort=False, validate='one_to_one')
    require(f.name.tolist() == raw.name.tolist(), 'Tier join reordered PAS axis')
    f = apply_view(f, view)
    # Never inherit the old external reference used by the upstream annotator.
    aliases = {canonical(c): c for c in f.chrom.unique()}
    require(len(aliases) == f.chrom.nunique(), 'Ambiguous chromosome aliases')
    refs = references.copy()
    refs['target_chrom'] = refs.target_chrom.map(lambda c: aliases.get(canonical(c)))
    refs = refs[refs.target_chrom.notna()]
    f = support_fields(f, refs)
    f['support_reference_variant'] = 'v2_primary'
    f['support_distance_basis'] = 'integrated_representative_interbase'
    f['count_origin'] = 'apadetect_0.4.1_study_uniform_recount'
    return f


def aggregate_library(matrix, pas_ids, expected_pas, barcodes, members, metacells):
    require(pas_ids == expected_pas, 'Library PAS row order differs from study catalog')
    require(len(set(barcodes)) == len(barcodes), 'Duplicated raw barcode')
    require(matrix.shape == (len(pas_ids), len(barcodes)), 'PAS-by-cell matrix dimensions differ')
    require(np.issubdtype(matrix.dtype, np.integer) and (matrix.data >= 0).all(), 'Counts must be nonnegative integer UMI')
    require(members.barcode.is_unique and members.cell_id.is_unique, 'Repeated frozen cell')
    require(set(members.barcode) <= set(barcodes), 'Frozen cells missing; suffix guessing is forbidden')
    require(set(members.metacell_id) == set(metacells.metacell_id), 'Metacell membership differs')
    return aggregate_counts(matrix.T.tocsr(), [{'original_barcode': b} for b in barcodes],
                            members.to_dict('records'), metacells.to_dict('records'))


def library_metacells(meta, library_key):
    """Select the frozen identity; library_id can be absent or a GEX alias."""
    require('library_key' in meta and meta.library_key.notna().all(), 'Frozen library keys missing')
    selected = meta[meta.library_key.eq(library_key)]
    require(not selected.empty, f'No frozen metacells for exact library key: {library_key}')
    return selected


def tier_function(path):
    tree = ast.parse(Path(path).read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'build_confidence_tiers']
    require(len(functions) == 1, 'Confidence-tier implementation missing')
    namespace = {'pd': pd}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['build_confidence_tiers']


def validate_handoff(root, inputs):
    plan = read_json(checked(root / 'plan.json', inputs))
    active = read_json(checked(root / 'active_run.json', inputs))
    state = read_json(checked(root / active['status_file'], inputs))
    require(state['status'] == 'COMPLETE' and not state['failures'], 'Upstream still running or failed')
    require([state[k] for k in ['study_catalogs_complete', 'recount_complete', 'studies_final_complete']] == [7, 161, 7],
            'Incomplete seven-study release')
    handoff = read_json(checked(root / 'handoff.json', inputs))
    require(handoff['status'] == 'STUDY_CATALOG_RECOUNT_ASNR_TERMINAL_COMPLETE', 'Final handoff missing')
    require(handoff['plan_sha256'] == digest(root / 'plan.json'), 'Handoff plan changed')
    require(handoff['execution_revision_sha256'] == digest(root / 'resume_plan.json'), 'Execution revision changed')
    require(plan['caller_version'] == '0.4.1' and len(plan['jobs']) == 210, 'Unexpected discovery release')
    require(len({j['library_id'] for j in plan['jobs']}) == 210, 'Duplicate discovery library')
    require(sum(j['in_fixed_cell_atlas'] for j in plan['jobs']) == 161, 'Fixed atlas scope changed')
    require({s['study'] for s in handoff['studies']} == set(plan['studies']) and len(handoff['studies']) == 7,
            'Study handoff differs from plan')
    return plan, handoff


def prepare(config, out, study=None):
    config = Path(config).resolve(); cfg = read_json(config); inputs = {}
    checked(config, inputs)
    checked(Path(__file__), inputs)
    require(cfg['analysis_route_policy'] == 'parallel_primary_sources' and cfg['reference_variant'] == 'v2_primary', 'Analysis route contract changed')
    root = Path(cfg['upstream_root']); baseline = Path(cfg['frozen_native_run'])
    if study is None:
        plan, handoff = validate_handoff(root, inputs)
    else:
        # A completed study can validate the adapter while another study finishes.
        # This branch cannot publish the full annotation or core-ready marker.
        plan = read_json(checked(root / 'plan.json', inputs))
        require(plan['caller_version'] == '0.4.1' and len(plan['jobs']) == 210, 'Unexpected pilot release')
        require(study in plan['studies'], 'Unknown pilot study')
        source = root / 'studies' / study / 'final'
        handoff = {'studies': [dict(study=study, output_dir=str(source),
                                   receipt=read_json(checked(source / 'complete.json', inputs)))]}
    registry = pd.read_csv(checked(cfg['native_libraries']['path'], inputs, cfg['native_libraries']['sha256']), sep='\t')
    require(len(registry) == 161 and registry.library_id.is_unique, 'Frozen library registry changed')
    upstream = {j['library_id']: j for j in plan['jobs'] if j['in_fixed_cell_atlas']}
    require(set(upstream) == set(registry.library_id), 'Upstream/frozen library sets differ')
    refs = pd.read_csv(checked(cfg['reference']['path'], inputs, cfg['reference']['sha256']), sep='\t', low_memory=False)
    require(set(refs.reference_variant) == {'v2_primary'}, 'Mixed reference versions')
    ann = checked(cfg['annotation_code']['path'], inputs, cfg['annotation_code']['sha256'])
    build_tiers = tier_function(ann)
    out = Path(out).resolve()
    require(out.is_relative_to(Path(cfg['output_root']).resolve()), 'Intake output must stay under the configured candidate root')
    require(not out.exists(), 'Use a fresh candidate directory; existing results are immutable')
    out.mkdir(parents=True)
    summaries = []; libraries = []; products = {}
    selected = [r for r in handoff['studies'] if study is None or r['study'] == study]
    require(bool(selected), 'Unknown pilot study')
    for rec in selected:
        source = root / 'studies' / rec['study'] / 'final'
        require(Path(rec['output_dir']) == source, 'Handoff output path differs')
        final = receipt(source, inputs)
        require(final == rec['receipt'], 'Handoff receipt changed')
        raw = pd.read_csv(source / 'features_asnr_terminal.tsv.gz', sep='\t', low_memory=False)
        require(raw.name.tolist() == (source / 'pas.tsv').read_text().splitlines(), 'Study PAS axis changed')
        lib = registry[registry.library_id.isin([k for k, j in upstream.items() if j['study'] == rec['study']])]
        require(lib.study_id.nunique() == 1 and lib.species.nunique() == 1, 'Ambiguous study mapping')
        require(lib.library_key.is_unique, 'Ambiguous frozen library keys')
        sid = lib.study_id.iloc[0]; sp = lib.species.iloc[0]
        require(raw.species.dropna().eq(sp).all(), 'Species mismatch')
        raw['pas_id'] = raw.pas_id.fillna('unannotated:' + raw.name)
        raw['primary_assignment_status'] = raw.primary_assignment_status.fillna('unavailable_coordinate')
        f = convert_features(raw, refs[refs.species.eq(sp)], build_tiers(raw), 'base')
        dest = out / 'studies' / sid; dest.mkdir(parents=True)
        meta = pd.read_csv(checked(baseline / 'studies' / sid / 'metacells.tsv.gz', inputs), sep='\t')
        require(meta.metacell_id.is_unique, 'Repeated metacell in frozen order')
        # Reuse expression/state metadata, never old PAS coverage or gene eligibility.
        meta = meta[[c for c in meta if 'pas' not in c.lower() and not c.startswith('genes_with_ge2')
                     and c != 'primary_detected_genes']]
        blocks = []; block_ids = []
        for row in lib.itertuples(index=False):
            d = root / 'recount' / row.library_id
            lr = receipt(d, inputs)
            require(lr['study'] == rec['study'] and lr['library_id'] == row.library_id, 'Library identity mismatch')
            members = pd.read_csv(checked(row.members_path, inputs, row.members_sha256), sep='\t')
            mc = library_metacells(meta, row.library_key)
            require(len(members) == row.cells and len(mc) == row.metacells, 'Frozen cells/metacells changed')
            matrix = sparse.load_npz(d / 'pas_by_cell.npz').tocsr()
            require(int(matrix.sum()) == lr['umi'], 'Recount molecule total differs from receipt')
            require(matrix[~truth(raw.quantifiable).to_numpy()].nnz == 0, 'Unavailable windows contain counts')
            block = aggregate_library(matrix, (d / 'pas.tsv').read_text().splitlines(), raw.name.tolist(),
                                      (d / 'barcodes.tsv').read_text().splitlines(), members, mc)
            blocks.append(block); block_ids.extend(mc.metacell_id)
            libraries.append(dict(study_id=sid, library_id=row.library_id, library_key=row.library_key,
                                  cells=len(members), metacells=len(mc), quantified_umi=lr['umi'],
                                  selected_cell_umi=int(block.sum()), exact_barcode_mapping=True))
        order = pd.Index(block_ids).get_indexer(meta.metacell_id)
        require((order >= 0).all() and len(order) == len(block_ids), 'Metacell reordering lost members')
        matrix = sparse.vstack(blocks, format='csr')[order]
        f['current_target_umi'] = np.asarray(matrix.sum(axis=0)).ravel()
        save(f, dest / 'features.tsv.gz'); save(meta, dest / 'metacells.tsv.gz')
        sparse.save_npz(dest / 'metacell_by_pas_counts.npz', matrix)
        require((sparse.load_npz(dest / 'metacell_by_pas_counts.npz') != matrix).nnz == 0, 'Matrix readback differs')
        for name in ['selected_metacells.tsv.gz', 'sampling_design.tsv', 'endpoint_library_design.tsv']:
            p = checked(baseline / 'studies' / sid / name, inputs)
            (dest / name).write_bytes(p.read_bytes())
        m = pd.read_csv(dest / 'selected_metacells.tsv.gz', sep='\t')
        require(np.array_equal(meta.metacell_id.iloc[m.matrix_row], m.metacell_id), 'Selected MC indices changed')
        summaries.append(dict(study_id=sid, species=sp, libraries=len(lib), cells=int(lib.cells.sum()),
                              metacells=len(meta), atlas_pas=len(f), molecules=int(matrix.sum()),
                              unavailable_windows=int((~truth(raw.quantifiable)).sum()),
                              primary_pas=int(f.keep_for_apa_primary.sum()),
                              supported50_primary_pas=int((f.keep_for_apa_primary & f.matched_source_count_50nt.ge(1)).sum())))
        products[sid] = {p.name: digest(p) for p in sorted(dest.iterdir()) if p.is_file()}
        print(json.dumps(summaries[-1]), flush=True)
    for path, sha in inputs.items():
        require(digest(path) == sha, f'Input changed during intake: {path}')
    full = study is None
    save(pd.DataFrame(libraries), out / 'frozen_native_membership_audit.tsv')
    audit = dict(status='PASS', libraries=len(libraries), cells=sum(r['cells'] for r in summaries),
                 metacells=sum(r['metacells'] for r in summaries),
                 checks=dict(exact_frozen_barcode_membership=True, integer_counts_conserved=True,
                             unquantifiable_not_measured_zero=True, old_PAS_columns_discarded=True))
    write_json(out / 'frozen_native_membership_audit.json', audit)
    write_json(out / 'intake.native.json', dict(status='prepared_full' if full else 'prepared_pilot',
        caller='APADETECT', caller_version='0.4.1', run_id=cfg['run_id']+'__base', analysis_route='base',
        analysis_route_policy=cfg['analysis_route_policy'],
        cut_coordinate_system='0-based-interbase', studies=summaries, inputs=inputs,
        outputs=products, full_manuscript_rerun_complete=False, new_BAM_recount=False))
    if full:
        require(len(summaries) == 7 and len(libraries) == 161, 'Cannot publish partial full-release intake')
        write_json(out / 'annotation.complete.json', dict(status='complete', rows=sum(s['atlas_pas'] for s in summaries),
                   studies=products, support_reference_variant='v2_primary', source='apadetect structural annotation plus verified tiers'))
    write_json(out / 'intake-complete.json', dict(status='READY_FOR_CORE' if full else 'PILOT_ONLY',
        inputs=inputs, studies=summaries, outputs=products, core_or_statistics_run=False,
        figures_replaced=False, full_manuscript_rerun_complete=False))


def derive_route(source, view, out):
    """Rebind a parallel primary route without aggregating or recounting again."""
    source, out = Path(source).resolve(), Path(out).resolve()
    require(view != 'base', 'Use the prepared base directory for the base route')
    require(not out.exists() and not out.is_relative_to(source), 'Use a fresh sibling route directory')
    old = read_json(source / 'intake.native.json')
    require(old.get('analysis_route', old.get('primary_view')) == 'base', 'Derive every route from base; never chain filters')
    inputs = {str(source / 'intake.native.json'): digest(source / 'intake.native.json'), str(Path(__file__).resolve()): digest(__file__)}
    summaries, products = [], {}
    out.mkdir(parents=True)
    for row in old['studies']:
        sid = row['study_id']; folder = source / 'studies' / sid
        dest = out / 'studies' / sid; dest.mkdir(parents=True)
        for name, h in old['outputs'][sid].items():
            checked(folder / name, inputs, h)
        f = pd.read_csv(folder / 'features.tsv.gz', sep='\t', low_memory=False)
        result = apply_view(f, view)
        require(result.atlas_pas_id.equals(f.atlas_pas_id), 'Parallel route changed stored PAS order')
        save(result, dest / 'features.tsv.gz')
        for name in old['outputs'][sid]:
            if name != 'features.tsv.gz':
                (dest / name).symlink_to(folder / name)
        summary = dict(row, primary_pas=int(result.keep_for_apa_primary.sum()),
                       supported50_primary_pas=int((result.keep_for_apa_primary & result.matched_source_count_50nt.ge(1)).sum()))
        summaries.append(summary)
        products[sid] = {p.name:digest(p) for p in dest.iterdir()}
    for name in ['frozen_native_membership_audit.json', 'frozen_native_membership_audit.tsv']:
        checked(source / name, inputs); (out / name).symlink_to(source / name)
    for path,h in inputs.items():
        require(digest(path)==h, f'Base route changed: {path}')
    run_base = old['run_id'].removesuffix('__base')
    record = dict(old, run_id=run_base+'__'+view, analysis_route=view,
                  analysis_route_policy='parallel_primary_sources', studies=summaries,
                  inputs=inputs, outputs=products, base_intake_sha256=digest(source / 'intake.native.json'))
    record.pop('primary_view', None)
    write_json(out / 'intake.native.json', record)
    full = (source / 'annotation.complete.json').exists()
    if full:
        write_json(out / 'annotation.complete.json', dict(status='complete', rows=sum(s['atlas_pas'] for s in summaries),
                   studies=products, support_reference_variant='v2_primary', analysis_route=view))
    write_json(out / 'intake-complete.json', dict(status='READY_FOR_CORE' if full else 'PILOT_ONLY',
        inputs=inputs, studies=summaries, outputs=products, analysis_route=view,
        interpretation='parallel primary data source; all effects/cohorts/statistics require fresh computation',
        shared_immutable_counts=True, new_BAM_recount=False, core_or_statistics_run=False,
        figures_replaced=False, full_manuscript_rerun_complete=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--pilot-study', help='Prepare one study; never publishes the full-run annotation marker')
    p.add_argument('--from-prepared', type=Path, help='Reuse one immutable prepared base for a parallel route')
    p.add_argument('--view', default='base')
    a = p.parse_args()
    if a.from_prepared:
        require(a.config is None and a.pilot_study is None, 'Route derivation uses the prepared base provenance')
        derive_route(a.from_prepared, a.view, a.out)
    else:
        require(a.config is not None and a.view == 'base', 'Prepare base first; derive other routes from it')
        prepare(a.config, a.out, a.pilot_study)
