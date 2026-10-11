"""Independent integer-count reconstruction of every fixed-core S-M pair."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FOLDER = HERE / 'sm-core'
ROOT = HERE.parents[1]


def read(path, **kwargs):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kwargs)


def equal(a, b):
    assert np.allclose(a, b, rtol=0, atol=1e-11, equal_nan=True)


def main():
    archive = np.load(ROOT / 'notes/2026-10-10-marker-reference-sensitivity/inputs.npz', allow_pickle=False)
    scores = read(HERE / 'independent_mc_scores.tsv.gz').set_index('metacell_id')
    roster = read(FOLDER / 'fixed_gene_sets.tsv')
    roster = roster[roster.gene_set.eq('core27')]
    eligibility = read(FOLDER / 'SM_pair_eligibility.tsv.gz')
    pairs = read(FOLDER / 'SM_pair_effects.tsv.gz')
    key = ['study_id', 'route', 'design', 'gene_id', 'library_key_a', 'library_key_b']
    assert not pairs.duplicated(key).any()
    indexed = pairs.set_index(key)
    n_qualification, n_measured = 0, 0
    for study, pe in eligibility.groupby('study_id', sort=True):
        species = pe.species.iloc[0]
        genes = set(roster.loc[roster.species.eq(species), 'gene_id'])
        pe = pe[pe.gene_id.isin(genes)]
        meta = pd.DataFrame(json.loads(str(archive[study + '_meta_json'])))
        selected = pd.DataFrame(json.loads(str(archive[study + '_selected_json'])))
        selected = selected[selected.purity.ge(.8) & selected.node.isin(['S', 'M'])]
        groups = {k: part.matrix_row.to_numpy(int) for k, part in selected.groupby(['route', 'library_key', 'node'])
                  if len(part) >= 3 and part.n_cells.sum() >= 30}
        features = pd.DataFrame(json.loads(str(archive[study + '_features_json'])))
        features['column'] = np.arange(len(features))
        data, indices, indptr = (archive[study + suffix] for suffix in ['_pas_data', '_pas_indices', '_pas_indptr'])
        assert np.issubdtype(data.dtype, np.integer) and (data >= 0).all()
        for gene, ge in pe.groupby('gene_id', sort=True):
            full = features[features.primary_gene_id.eq(gene) & features.coordinate_eligible]
            f = full[full.matched_source_count_50nt.ge(1)]
            if len(f) < 2:
                assert not ge.eligible.any()
                n_qualification += len(ge)
                continue
            signed = full.cut_site.to_numpy() * np.where(full.strand.eq('+'), 1, -1)
            positions = (f.cut_site.to_numpy() * np.where(f.strand.eq('+'), 1, -1) - signed.min()) / np.ptp(signed)
            equal(f.position, positions)
            counts = np.zeros((len(meta), len(f)), dtype=data.dtype)
            for j, col in enumerate(f.column):
                start, stop = indptr[col:col + 2]
                counts[indices[start:stop], j] = data[start:stop]
            profiles = {}
            for k, rows in groups.items():
                depth = counts[rows].sum(1)
                keep = rows[depth >= 20]
                c = counts[keep]
                profiles[k] = dict(rows=keep, c=c, cells=meta.n_cells.iloc[keep].sum())
            for row in ge.itertuples():
                a = profiles[row.route, row.library_key_a, 'M']
                b = profiles[row.route, row.library_key_b, 'S']
                joint = ((a['c'].sum(0) + b['c'].sum(0)) >= 10) & (((a['c'] > 0).sum(0) + (b['c'] > 0).sum(0)) >= 2)
                accepted = min(len(a['rows']), len(b['rows'])) >= 3 and min(a['cells'], b['cells']) >= 30 and joint.sum() >= 2
                assert accepted == row.eligible
                n_qualification += 1
                if not accepted:
                    continue
                found = indexed.loc[tuple(getattr(row, field) for field in key)]
                ua = np.mean(a['c'] / a['c'].sum(1)[:, None], axis=0)
                ub = np.mean(b['c'] / b['c'].sum(1)[:, None], axis=0)
                assert json.loads(found.PAS_ids) == f.atlas_pas_id.tolist()
                equal(positions, json.loads(found.position))
                equal(ua, json.loads(found.usage_a))
                equal(ub, json.loads(found.usage_b))
                equal(a['c'].sum(0), json.loads(found.pooled_molecules_a))
                equal(b['c'].sum(0), json.loads(found.pooled_molecules_b))
                equal(-(ua - ub) @ positions, found.delta_DWUI)
                assert json.loads(found.qualified_MC_ids_a) == meta.metacell_id.iloc[a['rows']].tolist()
                assert json.loads(found.qualified_MC_ids_b) == meta.metacell_id.iloc[b['rows']].tolist()
                for name in ['Cycle', 'Adjusted_Cycle']:
                    sa = scores.loc[meta.metacell_id.iloc[a['rows']], name].mean()
                    sb = scores.loc[meta.metacell_id.iloc[b['rows']], name].mean()
                    equal(sb - sa, found['delta_' + name + '_exact'])
                    ra = groups[row.route, row.library_key_a, 'M']
                    rb = groups[row.route, row.library_key_b, 'S']
                    equal(scores.loc[meta.metacell_id.iloc[rb], name].mean() - scores.loc[meta.metacell_id.iloc[ra], name].mean(),
                          found['delta_' + name + '_fullcohort'])
                n_measured += 1
        print(study, 'core raw reconstruction complete', flush=True)
    result = dict(status='PASS', raw_fixed_core_qualification_rows=n_qualification, raw_fixed_core_measured_pairs=n_measured,
                  raw_full_axis_reconstructed=True, low_depth_MC_reexcluded=True, raw_cycle_recovered=True,
                  missing_not_zero_imputed=True, tests_skipped=0)
    (HERE / 'sm_raw_independent_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
