"""Cohort and integer-count eligibility only; no usage, DWUI or cycle effects."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = ROOT / 'notes/2026-10-10-marker-reference-sensitivity'
STUDIES = {'GSE103221': ('mouse', ['CIP']), 'GSE115943': ('mouse', ['WOT']),
           'GSE118258': ('human', ['GSE118258']), 'GSE150311': ('human', ['GSE150311']),
           'GSE178325': ('human', ['178325_0618', '178325_0809']),
           'GSE221739': ('human', ['GSE221739'])}


def save(frame, name):
    pd.DataFrame(frame).to_csv(HERE / name, sep='\t', index=False, na_rep='NA',
                             compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def frame(z, key):
    return pd.DataFrame(json.loads(str(z[key])))


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    files = [OLD / 'inputs.npz',
             ROOT / 'notes/2026-10-10-proliferation-apa/metacell_cycle_scores.tsv.gz',
             ROOT / 'notes/2026-10-10-recurrence-utr-ortholog/measured_ortholog_members.tsv.gz',
             ROOT / 'notes/2026-10-10-directional-recurrence/observed_members.tsv.gz']
    hashes = {str(p): digest(p) for p in files}
    z = np.load(files[0], allow_pickle=False)
    cycle_ids = pd.read_csv(files[1], sep='\t', usecols=['metacell_id'])
    assert cycle_ids.metacell_id.is_unique
    cycles = set(cycle_ids.metacell_id)
    measured = pd.read_csv(files[2], sep='\t')
    core = measured[measured.contrast.eq('P_M') & measured.threshold.eq(.02) &
                    measured.human_direction.eq('proximal') & measured.mouse_direction.eq('proximal')]
    assert len(core) == 27 and core.homology_type.eq('ortholog_one2one').all()
    recurrent = pd.read_csv(files[3], sep='\t')
    recurrent = recurrent[recurrent.contrast.eq('P_M') & recurrent.threshold.eq(.02) & recurrent.direction.eq('proximal')]
    sets, roster = {}, []
    for species, expected in [('human', 144), ('mouse', 372)]:
        sets[species] = {'core27': set(core[species + '_gene_id']),
                         'complete_PM_proximal': set(recurrent.loc[recurrent.species.eq(species), 'gene_id'])}
        assert len(sets[species]['complete_PM_proximal']) == expected
        assert sets[species]['core27'] <= sets[species]['complete_PM_proximal']
        for name, genes in sets[species].items():
            roster.extend(dict(species=species, gene_set=name, gene_id=gene) for gene in sorted(genes))
    save(roster, 'fixed_gene_sets.tsv')
    cohort_rows, pair_rows, eligibility, coverage = [], [], [], []
    for sid, (species, routes) in STUDIES.items():
        meta = frame(z, sid + '_meta_json')
        selected = frame(z, sid + '_selected_json')
        selected = selected[selected.route.isin(routes) & selected.purity.ge(.8) & selected.node.isin(['S', 'M', 'P'])]
        assert selected.metacell_id.is_unique
        assert selected.metacell_id.isin(cycles).all()
        assert np.array_equal(meta.metacell_id.iloc[selected.matrix_row], selected.metacell_id)
        sampling = frame(z, sid + '_sampling_json')
        groups = {}
        for (route, lib, state), q in selected.groupby(['route', 'library_key', 'node'], sort=True):
            d = sampling[sampling.route.eq(route) & sampling.library_key.eq(lib) & sampling.state.eq(state)]
            assert len(d) == 1
            info = d.iloc[0].to_dict()
            info.update(n_mc=len(q), n_cells=int(q.n_cells.sum()), qualified=len(q) >= 3 and q.n_cells.sum() >= 30,
                        documented_source=state == 'S' and str(info['time']) in ['D0', 'source', 'Stage_I_day_0'],
                        pooled_time=str(info['time']) == 'D0_D7_pool')
            cohort_rows.append(info)
            if info['qualified']:
                groups[(route, lib, state)] = (q.matrix_row.to_numpy(int), info)
        candidates = []
        for route in routes:
            slibs = [(lib, item) for (r, lib, state), item in groups.items() if r == route and state == 'S']
            mlibs = [(lib, item) for (r, lib, state), item in groups.items() if r == route and state == 'M']
            for slib, (_, sinfo) in slibs:
                for mlib, (_, minfo) in mlibs:
                    schemes = ['cross_observed_S']
                    if sinfo['documented_source']:
                        schemes.append('cross_documented_source')
                    if slib == mlib:
                        schemes.append('same_library')
                    elif sinfo['time'] == minfo['time']:
                        schemes.append('same_time_cross_library')
                    for scheme in schemes:
                        pair = dict(species=species, study_id=sid, route=route, design=scheme,
                                    library_key_S=slib, library_key_M=mlib,
                                    time_S=sinfo['time'], time_M=minfo['time'],
                                    sample_unit_S=sinfo['sample_unit'], sample_unit_M=minfo['sample_unit'],
                                    replicate_status_S=sinfo['replicate_status'], replicate_status_M=minfo['replicate_status'],
                                    source_documented=sinfo['documented_source'], same_time=sinfo['time'] == minfo['time'],
                                    same_library=slib == mlib, n_mc_S=sinfo['n_mc'], n_mc_M=minfo['n_mc'],
                                    n_cells_S=sinfo['n_cells'], n_cells_M=minfo['n_cells'])
                        candidates.append(pair)
                        pair_rows.append(pair)
        features = frame(z, sid + '_features_json')
        features['column'] = np.arange(len(features))
        wanted = sets[species]['complete_PM_proximal']
        # This checks only pre-specified genes, not direction or usage.
        features = features[features.coordinate_eligible & features.matched_source_count_50nt.ge(1) &
                            features.primary_gene_id.isin(wanted)]
        data, indices, indptr = (z[sid + suffix] for suffix in ['_pas_data', '_pas_indices', '_pas_indptr'])
        nrows = int(z[sid + '_pas_shape'][0])
        gene_profiles = {}
        for gene, part in features.groupby('primary_gene_id', sort=True):
            if len(part) < 2:
                continue
            c = np.zeros((nrows, len(part)), dtype=data.dtype)
            for j, col in enumerate(part.column):
                start, end = indptr[col:col + 2]
                c[indices[start:end], j] = data[start:end]
            profiles = {}
            for key, (rr, info) in groups.items():
                keep = rr[c[rr].sum(1) >= 20]
                q = c[keep]
                profiles[key] = (len(keep), int(meta.n_cells.iloc[keep].sum()), q.sum(0), (q > 0).sum(0))
            gene_profiles[gene] = profiles
        for pair in candidates:
            for name, genes in sets[species].items():
                n_eligible = 0
                for gene in sorted(genes):
                    why = 'fewer_than_two_supported_coordinate_PAS'
                    if gene in gene_profiles:
                        a = gene_profiles[gene][(pair['route'], pair['library_key_S'], 'S')]
                        b = gene_profiles[gene][(pair['route'], pair['library_key_M'], 'M')]
                        why = ('insufficient_gene_qualified_MC' if min(a[0], b[0]) < 3 else
                               'insufficient_gene_qualified_cells' if min(a[1], b[1]) < 30 else
                               'fewer_than_two_joint_supported_PAS' if (((a[2] + b[2]) >= 10) & ((a[3] + b[3]) >= 2)).sum() < 2 else '')
                    okay = not why
                    n_eligible += okay
                    eligibility.append(pair | dict(gene_set=name, gene_id=gene, eligible=okay, reason=why))
                coverage.append(pair | dict(gene_set=name, n_fixed_genes=len(genes), n_measurable=n_eligible))
        print(sid, 'cohorts', sum(row['study_id'] == sid for row in cohort_rows),
              'SM_candidate_designs', len(candidates), 'target_genes_axis_eligible', len(gene_profiles), flush=True)
    save(cohort_rows, 'cohort_preflight.tsv')
    save(pair_rows, 'candidate_pair_design.tsv')
    save(eligibility, 'target_gene_eligibility.tsv.gz')
    save(coverage, 'target_pair_coverage.tsv')
    c = pd.DataFrame(coverage)
    route = []
    for (species, sid, routes_) in [(v[0], k, v[1]) for k, v in STUDIES.items()]:
        for route_ in routes_:
            for design in ['same_library', 'same_time_cross_library', 'cross_documented_source', 'cross_observed_S']:
                for name, genes in sets[species].items():
                    sub = c[c.study_id.eq(sid) & c.route.eq(route_) & c.design.eq(design) & c.gene_set.eq(name)]
                    e = pd.DataFrame(eligibility)
                    sube = e[e.study_id.eq(sid) & e.route.eq(route_) & e.design.eq(design) & e.gene_set.eq(name) & e.eligible]
                    route.append(dict(species=species, study_id=sid, route=route_, design=design, gene_set=name,
                                      candidate_pairs=len(sub), n_fixed_genes=len(genes),
                                      n_measurable_any_pair=sube.gene_id.nunique(),
                                      n_S_libraries=sub.library_key_S.nunique(), n_M_libraries=sub.library_key_M.nunique(),
                                      min_pair_genes=sub.n_measurable.min() if len(sub) else np.nan,
                                      max_pair_genes=sub.n_measurable.max() if len(sub) else np.nan))
    save(route, 'route_preflight.tsv')
    assert all(digest(p) == expected for p, expected in hashes.items())
    validation = dict(status='PASS_ELIGIBILITY_ONLY', input_sha256=hashes,
                      no_APA_usage_or_DWUI_or_cycle_effect_computed=True,
                      roster_sizes={species: {name: len(genes) for name, genes in current.items()} for species, current in sets.items()},
                      candidate_design_rows=len(pair_rows), target_pair_gene_qualification_rows=len(eligibility),
                      preflight_script_sha256=digest(Path(__file__)))
    (HERE / 'preflight_validation.json').write_text(json.dumps(validation, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    main()
