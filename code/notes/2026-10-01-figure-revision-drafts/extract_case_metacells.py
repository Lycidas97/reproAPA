"""Recover display distributions on the exact frozen libraries and PAS axes."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'display-data'
PROJECTS = Path('/Users/lycidas/work/科研/projects')
MOUSE = PROJECTS / '2026-09-02_mouse-time-aware-integration-pilot'
HUMAN = PROJECTS / '2026-09-09_human-time-aware-integration'
INPUTS = {}


def register(path):
    with path.open('rb') as stream:
        INPUTS[str(path)] = hashlib.file_digest(stream, 'sha256').hexdigest()
    return path


def read(path, **kwargs):
    return pd.read_csv(register(path), sep='\t', **kwargs)


def main():
    triplets = read(ROOT/'notes/2026-10-01-figure3-matched-metrics/weighted_triplets.tsv.gz')
    triplets = triplets[triplets.gene_symbol.eq('Nudt21') & triplets.setting.eq('supported50_p08')]
    pairs = read(ROOT/'notes/2026-09-30-multimetric-statistics/pair_effects.tsv.gz')
    pairs = pairs[pairs.gene_symbol.str.upper().isin(['SZRD1', 'RABIF']) & pairs.route.isin(['CIP', 'GSE118258'])]
    matched = read(ROOT/'table/Supplementary-Source-Data-20261001-matched-metrics/matched_state_vectors.tsv.gz')
    matched = matched[matched.gene_symbol.eq('Nudt21') & matched.setting.eq('supported50_p08')]
    effects = read(ROOT/'notes/2026-10-01-figure4-state-pas-v6/data/case_effects.tsv')
    wanted, axes, targets = [], {}, {}
    for (study, route, gene), g in triplets.groupby(['study_id', 'route', 'gene_id']):
        r = matched[matched.study_id.eq(study) & matched.route.eq(route)].iloc[0]
        axes[(study, gene)] = (json.loads(r.PAS_ids), np.asarray(json.loads(r.position_weights)))
        for state, side in [('S', 'b'), ('P', 'a'), ('M', 'a')]:
            targets[(study, route, gene, state)] = float(r[state+'_DWUI'])
            for lib, q in g.groupby('library_key_'+side):
                cols = [state+'_n_mc', state+'_n_cells', state+'_pas_counts', state+'_WPP']
                assert all(q[col].nunique() == 1 for col in cols)
                wanted.append(dict(study_id=study, route=route, gene_id=gene, gene_symbol='Nudt21',
                    library_key=lib, state=state, n_mc=int(q[cols[0]].iloc[0]), n_cells=int(q[cols[1]].iloc[0]),
                    pas_counts=int(q[cols[2]].iloc[0]), expected_DWUI=q[cols[3]].iloc[0], weight=q.weight.sum()))
    for (study, route, gene), g in pairs.groupby(['study_id', 'route', 'gene_id']):
        gene_symbol = g.gene_symbol.iloc[0].upper()
        f = effects[effects.study_id.eq(study) & effects.route.eq(route) & effects.gene_id.eq(gene)].set_index('contrast')
        ids = json.loads(f.loc['P_S'].PAS_ids)
        positions = np.asarray([int(p.split(':')[3]) for p in ids])
        distance = positions-positions.min() if ids[0].split(':')[4] == '+' else positions.max()-positions
        axes[(study, gene)] = (ids, distance/distance.max())
        for state, contrast, side, metric in [('S', 'P_S', 'b', 'reference_WPP'),
                                              ('P', 'P_S', 'a', 'target_WPP'),
                                              ('M', 'P_M', 'b', 'reference_WPP')]:
            targets[(study, route, gene, state)] = float(f.loc[contrast, metric])
            for lib, q in g[g.contrast.eq(contrast)].groupby('library_key_'+side):
                cols = ['n_mc_'+side, 'n_cells_'+side, 'pas_counts_'+side, metric]
                assert all(q[col].nunique() == 1 for col in cols)
                wanted.append(dict(study_id=study, route=route, gene_id=gene, gene_symbol=gene_symbol,
                    library_key=lib, state=state, n_mc=int(q[cols[0]].iloc[0]), n_cells=int(q[cols[1]].iloc[0]),
                    pas_counts=int(q[cols[2]].iloc[0]), expected_DWUI=q[cols[3]].iloc[0], weight=q.weight.sum()))
    wanted = pd.DataFrame(wanted)
    keys = ['study_id', 'route', 'gene_id', 'state']
    assert not wanted.duplicated(keys+['library_key']).any()
    np.testing.assert_allclose(wanted.groupby(keys).weight.sum(), 1, atol=1e-12)
    points, libraries = [], []
    for study, g in wanted.groupby('study_id'):
        stem = study.lower()+('_human_pluripotency_reprogramming' if study == 'GSE118258'
                              else '_pluripotency_reprogramming_core')
        folder = ((HUMAN/'data/results_metacell_pas_v2/studies') if study == 'GSE118258'
                  else (MOUSE/'data/metacell-pas-v7-v1/studies')) / stem
        meta = read(MOUSE/'data/apa-majority-shared-v1/profiles'/stem/'selected_metacells.tsv.gz')
        original = read(folder/'metacells.tsv.gz', usecols=['metacell_id'])
        assert np.array_equal(original.metacell_id.iloc[meta.matrix_row], meta.metacell_id)
        meta = meta[meta.library_key.isin(g.library_key) & meta.node.isin(['S', 'P', 'M']) & meta.purity.ge(.8)]
        meta = meta.sort_values('matrix_row').reset_index(drop=True)
        pids = {pid for gene in g.gene_id.unique() for pid in axes[(study, gene)][0]}
        feature_path = register(folder/'features.tsv.gz')
        pieces = []
        for chunk in pd.read_csv(feature_path, sep='\t', usecols=['atlas_pas_id', 'feature_index_1based', 'primary_gene_id'], chunksize=100000):
            pieces.append(chunk[chunk.atlas_pas_id.isin(pids)])
        features = pd.concat(pieces).set_index('atlas_pas_id')
        assert set(features.index) == pids and features.index.is_unique
        matrix_path = register(folder/'metacell_by_pas_counts.npz')
        audit = json.loads((ROOT/f'notes/2026-10-01-figure3-matched-metrics/{study}_extraction_audit.json').read_text())
        for path in [feature_path, matrix_path]:
            assert INPUTS[str(path)] == audit['inputs'][str(path)]['sha256']
        print('Loading frozen counts', study, flush=True)
        matrix = sparse.load_npz(matrix_path).tocsr()
        selected = matrix[meta.matrix_row.to_numpy()][:, features.feature_index_1based.to_numpy()-1].toarray()
        del matrix
        for gene, q in g.groupby('gene_id'):
            ids, weights = axes[(study, gene)]
            assert features.loc[ids].primary_gene_id.eq(gene).all()
            counts = selected[:, features.index.get_indexer(ids)]
            total = counts.sum(axis=1)
            values = np.divide(counts@weights, total, out=np.full(len(meta), np.nan), where=total > 0)
            for r in q.itertuples():
                mask = (meta.route.eq(r.route) & meta.library_key.eq(r.library_key) & meta.node.eq(r.state)).to_numpy() & (total >= 20)
                block = meta.loc[mask, ['study_id', 'route', 'library_key', 'node', 'metacell_id', 'n_cells', 'purity']].copy()
                block = block.rename(columns={'node': 'state'})
                block['gene_id'], block['gene_symbol'] = gene, r.gene_symbol
                block['PAS_UMIs'], block['DWUI'] = total[mask], values[mask]
                assert len(block) == r.n_mc and block.n_cells.sum() == r.n_cells
                assert block.PAS_UMIs.sum() == r.pas_counts
                error = abs(block.DWUI.mean()-r.expected_DWUI)
                assert error < 1e-8, (study, gene, r.state, error)
                points.append(block)
                libraries.append(dict(r._asdict(), DWUI=block.DWUI.mean(), error=error))
        print('Recovered', study, 'libraries', len(g), flush=True)
    points = pd.concat(points, ignore_index=True)
    libraries = pd.DataFrame(libraries).drop(columns='Index')
    comparisons = []
    for key, q in libraries.groupby(keys):
        observed = q.DWUI.to_numpy()@q.weight.to_numpy()
        assert abs(observed-targets[key]) < 1e-8
        p = points.set_index(keys).loc[[key]]
        comparisons.append(dict(zip(keys, key), gene_symbol=q.gene_symbol.iloc[0], n_metacells=len(p),
            n_libraries=len(q), DWUI=observed, expected_DWUI=targets[key], min=p.DWUI.min(),
            q25=p.DWUI.quantile(.25), median=p.DWUI.median(), q75=p.DWUI.quantile(.75),
            max=p.DWUI.max(), zero_fraction=p.DWUI.eq(0).mean()))
    for data, name in [(points, 'case-metacells.tsv'), (libraries, 'case-library-means.tsv'),
                       (pd.DataFrame(comparisons), 'case-distribution-summary.tsv')]:
        data.to_csv(OUT/name, sep='\t', index=False, float_format='%.17g')
    report = dict(status='PASS', metacells=len(points), library_state_records=len(libraries),
        max_library_error=libraries.error.max(), route_state_summaries=len(comparisons),
        inputs_sha256=INPUTS, extractor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        outputs_sha256={name: hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in
                       ['case-metacells.tsv', 'case-library-means.tsv', 'case-distribution-summary.tsv']})
    (HERE/'case-metacells-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(pd.DataFrame(comparisons).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
