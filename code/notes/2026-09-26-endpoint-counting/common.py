from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
BASE = Path('/root/nfsdata/apadetect/analysis_projects')
OLD = BASE / 'reproapa_pas_reliability_2026_09_13'
MOUSE = BASE / 'mouse_time_aware_integration_pilot_v1'
HUMAN = BASE / 'human_time_aware_integration_v1'
PAS = {'mouse': MOUSE / 'results_metacell_pas_v7_v1', 'human': HUMAN / 'results_metacell_pas_v2'}
MODES = ['original', 'atlas_dedup', 'transferred_observed', 'transferred', 'primary_clipped', 'all_clipped']
STUDIES = {
    'GSE118258': 'gse118258_human_pluripotency_reprogramming',
    'GSE103221': 'gse103221_pluripotency_reprogramming_core',
    'GSE178325': 'gse178325_human_chemical_pluripotency_reprogramming',
    'GSE115943': 'gse115943_pluripotency_reprogramming_core',
}


def read(path, **kwargs):
    return pd.read_csv(path, sep='\t', low_memory=False, **kwargs)


def save(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, float_format='%.12g',
                 compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dump(value, path):
    Path(path).write_text(json.dumps(value, indent=2, default=lambda x: x.item(), allow_nan=False) + '\n')


def canon(chrom):
    return chrom if str(chrom).startswith('chr') else ('chrM' if chrom == 'MT' else 'chr' + str(chrom))


def collapse(pairs):
    keys = ['species', 'study_id', 'route', 'view', 'mode', 'gene_id']
    x = pairs.groupby(keys + ['time_a', 'sample_unit_a', 'time_b', 'sample_unit_b']).delta.mean()
    x = x.groupby(keys + ['time_a', 'sample_unit_a', 'time_b']).mean()
    x = x.groupby(keys + ['time_a', 'sample_unit_a']).mean()
    dates = x.groupby(keys + ['time_a']).mean().rename('delta').reset_index()
    effects = dates.groupby(keys).delta.mean().reset_index()
    coverage = pairs.groupby(keys).agg(n_pairs=('delta', 'size'), n_targets=('library_key_a', 'nunique'),
                                      n_sources=('library_key_b', 'nunique')).reset_index()
    return effects.merge(coverage, on=keys, validate='one_to_one'), dates


def library_stats(counts, features, meta, allowed, mode, common_good=None):
    genes = pd.Index(features.primary_gene_id.unique())
    gi = genes.get_indexer(features.primary_gene_id)
    mapping = sparse.csr_matrix((np.ones(len(features)), (np.arange(len(features)), gi)),
                                shape=(len(features), len(genes)))
    weights = mapping.multiply(features.relative.to_numpy()[:, None]).tocsr()
    totals = (counts @ mapping).toarray()
    good = (totals >= 20) & genes.isin(allowed)[None, :]
    if common_good is not None:
        good &= common_good
    n_mc = good.sum(0)
    n_cells = (good * meta.n_cells.to_numpy()[:, None]).sum(0)
    values = np.divide((counts @ weights).toarray(), totals, out=np.zeros_like(totals), where=totals > 0)
    valid = (n_mc >= 3) & (n_cells >= 30)
    if mode + '_available' in features:
        bad = features.loc[~features[mode + '_available'], 'primary_gene_id']
        valid &= ~genes.isin(bad)
    apa = np.divide(np.where(good, values, 0).sum(0), n_mc,
                    out=np.full(len(genes), np.nan), where=valid)
    coo = counts.tocoo()
    take = good[coo.row, gi[coo.col]]
    site_sum = np.bincount(coo.col[take], weights=coo.data[take], minlength=len(features))
    site_mc = np.bincount(coo.col[take], minlength=len(features))
    return dict(APA=apa, n_mc=n_mc, n_cells=n_cells, site_sum=site_sum, site_mc=site_mc), good


def paired_stats(a, b, features):
    genes = pd.Index(features.primary_gene_id.unique())
    gi = genes.get_indexer(features.primary_gene_id)
    supported = ((a['site_sum'] + b['site_sum']) >= 10) & ((a['site_mc'] + b['site_mc']) >= 2)
    n_sites = np.bincount(gi, weights=supported, minlength=len(genes))
    valid = np.isfinite(a['APA']) & np.isfinite(b['APA']) & (n_sites >= 2)
    return a['APA'] - b['APA'], valid
