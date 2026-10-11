"""Requalify PAS states and endpoint effects on the new study-native axes.

Only the expression scaffold and physical endpoint design are frozen. No old
gene roster, PAS coordinate, PAS-qualified cell list, or effect is used here.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from polyapipe_prepare import digest, save

METRICS = ['DWUI', 'RWUI', 'PDUI', 'PPUI']
SETTINGS = [('supported50_p08', True, .8, 3), ('supported50_p09', True, .9, 3),
            ('supported50_p08_mc5', True, .8, 5), ('expanded_primary_p08', False, .8, 3)]


def weights(frame):
    """Equal source date/unit then target date/unit; equal captures last."""
    levels = [['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b'],
              ['time_a', 'sample_unit_a', 'time_b'], ['time_a', 'sample_unit_a'], ['time_a'], []]
    q = frame.reset_index(drop=True).copy()
    result = 1 / q.groupby(levels[0], dropna=False).time_a.transform('size').to_numpy(float)
    for child, parent in zip(levels[:-1], levels[1:]):
        distinct = q[child].drop_duplicates()
        if not parent:
            result /= len(distinct)
        else:
            n = distinct.groupby(parent, dropna=False).size().rename('_n').reset_index()
            result /= q[parent].merge(n, on=parent, how='left', validate='many_to_one')._n.to_numpy()
    assert abs(result.sum()-1) < 1e-12
    return result


def metric_weights(sites):
    signed = sites.cut_site.to_numpy() * np.where(sites.strand.eq('+'), 1, -1)
    order = np.argsort(signed, kind='stable')
    rank = np.empty(len(sites), int); rank[order] = np.arange(1, len(sites)+1)
    return np.column_stack([sites.position, rank / len(sites), rank == len(sites), rank == 1])


def state(data, rows, meta):
    c = data[rows]
    good = c.sum(1) >= 20
    rows = rows[good]; c = c[good]
    usage = (c / c.sum(1)[:, None]).mean(0) if len(c) else np.zeros(data.shape[1])
    return dict(rows=rows, n_mc=len(rows), n_cells=int(meta.n_cells.iloc[rows].sum()), usage=usage,
                molecules=c.sum(0), nonzero_mc=(c > 0).sum(0))


def reason(a, b, minimum_mc):
    if min(a['n_mc'], b['n_mc']) < minimum_mc:
        return 'insufficient_gene_qualified_metacells'
    if min(a['n_cells'], b['n_cells']) < 30:
        return 'insufficient_gene_qualified_cells'
    joint = (a['molecules'] + b['molecules'] >= 10) & (a['nonzero_mc'] + b['nonzero_mc'] >= 2)
    if joint.sum() < 2:
        return 'fewer_than_two_joint_supported_PAS'
    return ''


def study(folder, settings=SETTINGS):
    required = ['features.tsv.gz', 'metacell_by_pas_counts.npz', 'metacells.tsv.gz',
                'selected_metacells.tsv.gz', 'endpoint_library_design.tsv', 'sampling_design.tsv']
    inputs = {str(folder/n): digest(folder/n) for n in required}
    f = pd.read_csv(folder/'features.tsv.gz', sep='\t', low_memory=False)
    f = f[f.keep_for_apa_primary.astype(str).str.lower().eq('true') & f.primary_gene_id.notna()].copy()
    signed = f.cut_site * np.where(f.strand.eq('+'), 1, -1)
    shape = f.groupby('primary_gene_id').agg(n=('cut_site','nunique'), c=('chrom','nunique'), s=('strand','nunique'))
    f['position'] = (signed - signed.groupby(f.primary_gene_id).transform('min')) / (
        signed.groupby(f.primary_gene_id).transform('max') - signed.groupby(f.primary_gene_id).transform('min'))
    f['coordinate_eligible'] = f.primary_gene_id.isin(shape.query('n>=2 and c==1 and s==1').index)
    original = pd.read_csv(folder/'metacells.tsv.gz', sep='\t')
    meta = pd.read_csv(folder/'selected_metacells.tsv.gz', sep='\t').reset_index(drop=True)
    assert meta.metacell_id.is_unique
    assert np.array_equal(original.metacell_id.iloc[meta.matrix_row], meta.metacell_id)
    counts = sparse.load_npz(folder/'metacell_by_pas_counts.npz').tocsr()[meta.matrix_row.to_numpy()].tocsc()
    assert np.issubdtype(counts.dtype, np.integer) and (counts.data >= 0).all()
    design = pd.read_csv(folder/'endpoint_library_design.tsv', sep='\t')
    out = folder/'effects'; out.mkdir(exist_ok=False)
    save(f, out/'primary_axes.tsv.gz')
    all_summaries = []
    for setting, supported, purity, minimum_mc in settings:
        if supported:
            sf = f[f.coordinate_eligible & f.matched_source_count_50nt.ge(1)].copy()
        else:
            sf = f[f.coordinate_eligible].copy()
        groups = {key: np.asarray(ids) for key, ids in meta[meta.purity.ge(purity)].groupby(
                  ['route', 'library_key', 'node'], sort=False).groups.items()}
        groups = {key: rows for key, rows in groups.items()
                  if len(rows) >= minimum_mc and meta.n_cells.iloc[rows].sum() >= 30}
        pairs, routes, eligibility, libstates = [], [], [], []
        background = shape.copy()
        background['retained_PAS'] = sf.groupby('primary_gene_id').size().reindex(background.index, fill_value=0)
        background['status'] = np.select([background.c.ne(1) | background.s.ne(1) | background.n.lt(2),
            background.retained_PAS.lt(2)], ['unusable_primary_axis', 'fewer_than_two_retained_PAS'], default='endpoint_requalification_required')
        save(background.reset_index(), out/f'{setting}_gene_background.tsv.gz')
        study_id = str(meta.study_id.iloc[0]).replace('WOT_GSE115943', 'GSE115943')
        species = f.species.iloc[0]
        for gene, sites in sf.groupby('primary_gene_id', sort=True):
            base = dict(species=species, study_id=study_id, gene_id=gene,
                        gene_symbol=sites.primary_gene_symbol.iloc[0], setting=setting)
            if len(sites) < 2:
                for d in design.to_dict('records'):
                    eligibility.append(dict(d, **base, status='not_measurable', reason='fewer_than_two_retained_PAS'))
                continue
            w = metric_weights(sites)
            data = counts[:, sites.feature_index_1based.to_numpy()-1].toarray()
            states = {key: state(data, rows, meta) for key, rows in groups.items()}
            site_ids = json.dumps(sites.atlas_pas_id.tolist())
            for route in sorted(meta.route.dropna().unique()):
                # Figure 2 uses all gene-qualified MCs from prequalified groups
                # for route PAS support; its 30-cell gate is at group level.
                selected = [s for (r, _, _), s in states.items() if r == route]
                if not selected:
                    continue
                joint = (sum(s['molecules'] for s in selected) >= 10) & (sum(s['nonzero_mc'] for s in selected) >= 2)
                if joint.sum() < 2:
                    continue
                for (r, lib, node), s in states.items():
                    if r != route or s['n_mc'] < minimum_mc:
                        continue
                    mu = s['usage'] @ w
                    libstates.append(dict(base, route=r, library_key=lib, node=node,
                        position_index=mu[0], qualified_metacells=s['n_mc'], qualified_cells=s['n_cells'],
                        pas_counts=int(s['molecules'].sum()), PAS_number=len(sites),
                        **{m: float(v) for m, v in zip(METRICS, mu)}))
            for (route, contrast, comparison), planned in design.groupby(['route','contrast','design'], sort=False):
                eligible = []
                for d in planned.to_dict('records'):
                    a = states.get((route, d['library_key_a'], 'P'))
                    b = states.get((route, d['library_key_b'], contrast[-1]))
                    why = 'state_absent_at_purity_threshold' if a is None or b is None else reason(a, b, minimum_mc)
                    eligibility.append(dict(d, **base, status='not_measurable' if why else 'eligible', reason=why,
                        n_mc_a=0 if a is None else a['n_mc'], n_mc_b=0 if b is None else b['n_mc'],
                        n_cells_a=0 if a is None else a['n_cells'], n_cells_b=0 if b is None else b['n_cells']))
                    if why:
                        continue
                    delta = a['usage'] - b['usage']
                    assert abs(delta.sum()) < 1e-12
                    eligible.append(dict(d, **base, delta_usage=delta, target_usage=a['usage'], reference_usage=b['usage'],
                        metacell_ids_a=meta.metacell_id.iloc[a['rows']].tolist(), metacell_ids_b=meta.metacell_id.iloc[b['rows']].tolist(),
                        n_mc_a=a['n_mc'], n_mc_b=b['n_mc'], n_cells_a=a['n_cells'], n_cells_b=b['n_cells'],
                        pas_counts_a=int(a['molecules'].sum()), pas_counts_b=int(b['molecules'].sum()),
                        molecules_a=a['molecules'].tolist(), molecules_b=b['molecules'].tolist()))
                if not eligible:
                    continue
                eq = pd.DataFrame(eligible); ww = weights(eq)
                a = ww @ np.stack(eq.target_usage); b = ww @ np.stack(eq.reference_usage); delta = a-b
                assert np.isclose(a.sum(),1) and np.isclose(b.sum(),1) and abs(delta.sum())<1e-12
                effect = delta @ w
                assert abs(effect).max() <= 1+1e-12 and np.ptp(delta) <= 2+1e-12
                record = dict(base, route=route, contrast=contrast, design=comparison, PAS_number=len(sites),
                    n_pairs=len(eq), PAS_ids=site_ids, target_usage=json.dumps(a.tolist()), reference_usage=json.dumps(b.tolist()),
                    delta_usage=json.dumps(delta.tolist()), MPRO=float(np.ptp(delta)),
                    APA_coverage=float(ww @ ((eq.pas_counts_a/eq.n_mc_a+eq.pas_counts_b/eq.n_mc_b)/2)),
                    reference_concentration=float(np.square(b).sum()), reference_min_usage=float(b.min()),
                    reference_max_usage=float(b.max()), rwui_min_weight=float(w[:,1].min()), rwui_max_weight=1.)
                record.update({'delta_'+m: float(v) for m,v in zip(METRICS,effect)})
                record.update({'reference_'+m: float(v) for m,v in zip(METRICS,b@w)})
                record.update({'target_'+m: float(v) for m,v in zip(METRICS,a@w)})
                routes.append(record)
                for rec, weight in zip(eligible, ww):
                    rec['weight']=float(weight); rec['PAS_ids']=site_ids
                    rec.update({'delta_'+m: float(v) for m,v in zip(METRICS,rec['delta_usage']@w)})
                    rec['MPRO']=float(np.ptp(rec['delta_usage']))
                    for key in ['delta_usage','target_usage','reference_usage','metacell_ids_a','metacell_ids_b','molecules_a','molecules_b']:
                        rec[key]=json.dumps(np.asarray(rec[key]).tolist())
                    pairs.append(rec)
            if len(routes) and len(routes) % 5000 == 0:
                print(folder.name, setting, gene, len(routes), flush=True)
        for name, records in [('pair_effects',pairs),('route_effects',routes),('eligibility',eligibility),('library_states',libstates)]:
            table = pd.DataFrame(records)
            if table.empty:
                columns = ['species','study_id','setting','route','gene_id','gene_symbol']
                columns += (['node','library_key','position_index','qualified_metacells','qualified_cells']
                            if name == 'library_states' else ['contrast','design','status','reason'])
                table = pd.DataFrame(columns=columns)
            save(table, out/f'{setting}_{name}.tsv.gz')
        states_frame = pd.DataFrame(libstates)
        if len(states_frame):
            nodes = states_frame.groupby(['species','study_id','setting','route','node','gene_id','gene_symbol']).agg(
                position_index=('position_index','mean'), n_libraries=('library_key','nunique'),
                qualified_metacells=('qualified_metacells','sum')).reset_index()
            save(nodes, out/f'{setting}_node_positions.tsv.gz')
        else:
            save(pd.DataFrame(columns=['species','study_id','setting','route','node','gene_id','gene_symbol',
                 'position_index','n_libraries','qualified_metacells']), out/f'{setting}_node_positions.tsv.gz')
        summary = dict(study_id=study_id, setting=setting, genes_with_primary_axis=int(f[f.coordinate_eligible].primary_gene_id.nunique()),
            retained_genes=sf.primary_gene_id.nunique(), endpoint_tests=len(eligibility), qualified_pairs=len(pairs),
            route_effects=len(routes), library_states=len(libstates))
        all_summaries.append(summary); print(json.dumps(summary), flush=True)
    assert all(digest(p)==value for p,value in inputs.items())
    (out/'complete.json').write_text(json.dumps(dict(status='core_effects_complete_not_full_rerun',
        inputs=inputs, summaries=all_summaries, script_sha256=digest(__file__)),indent=2)+'\n')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True)
    p.add_argument('--primary-only', action='store_true');a=p.parse_args()
    study(a.study, SETTINGS[:1] if a.primary_only else SETTINGS)
