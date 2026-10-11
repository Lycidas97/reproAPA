"""Build a count-conserving study atlas from the completed native polyApipe run."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(frame, path):
    frame.to_csv(path, sep='\t', index=False, float_format='%.17g',
                 compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)


def clusters(cuts):
    """The manuscript's adjacent <=25, total width <=50 rule, without chaining."""
    labels = np.zeros(len(cuts), dtype=np.int64)
    first = previous = int(cuts[0])
    label = 0
    for i, cut in enumerate(cuts):
        if cut - previous > 25 or cut - first > 50:
            label += 1
            first = int(cut)
        labels[i] = label
        previous = int(cut)
    return labels


def build(source, out, frozen):
    out.mkdir(parents=True, exist_ok=False)
    libraries = pd.read_csv(source / 'input_libraries.tsv', sep='\t')
    assert len(libraries) == 161 and libraries.library_id.is_unique
    members = pd.read_csv(frozen, sep='\t')
    inputs, summaries = {}, []
    inputs[str(source / 'input_libraries.tsv')] = digest(source / 'input_libraries.tsv')
    inputs[str(frozen)] = digest(frozen)
    for study, lib in libraries.groupby('study_id', sort=False):
        dest = out / 'studies' / study
        dest.mkdir(parents=True)
        frames, matrices, metas, audits = [], [], [], []
        for row in lib.itertuples(index=False):
            folder = source / 'polyapipe' / study / row.library_id
            files = [folder / n for n in ['site_inventory.tsv.gz', 'metacell_index.tsv',
                     'metacell_by_site_counts.npz', 'count_audit.json', 'complete.json']]
            inputs.update({str(p): digest(p) for p in files})
            f = pd.read_csv(files[0], sep='\t')
            m = pd.read_csv(files[1], sep='\t')
            c = sparse.load_npz(files[2]).tocsr()
            audit = json.loads(files[3].read_text())
            assert c.shape == (len(m), len(f))
            assert np.issubdtype(c.dtype, np.integer) and (c.data >= 0).all()
            assert int(c.sum()) == audit['integer_molecules']
            assert len(m) == row.metacells and f.site_id.is_unique and m.metacell_id.is_unique
            mp = members[(members.library_id == row.library_id) & (members.role == 'membership')].iloc[0]
            assert digest(mp.path) == mp.sha256
            inputs[mp.path] = mp.sha256
            membership = pd.read_csv(mp.path, sep='\t', compression='gzip')
            assert set(m.metacell_id) == set(membership.metacell_id)
            assert len(membership) == row.cells
            assert membership.barcode.str.split('-').str[0].is_unique
            m['library_key'] = row.library_key
            m['n_cells'] = m.metacell_id.map(membership.groupby('metacell_id').size())
            assert m.n_cells.notna().all() and int(m.n_cells.sum()) == row.cells
            f['native_column_0based'] = np.arange(len(f))
            f['library_id'] = row.library_id
            f['native_cut_site_1based'] = f.cut_site
            # GFF is inclusive: plus cleavage follows end; minus precedes start.
            f['cut_site'] = f.cut_site - f.strand.eq('-').astype(int)
            f['native_molecules'] = np.asarray(c.sum(0)).ravel()
            frames.append(f); matrices.append(c); metas.append(m)
            audits.append(dict(library_id=row.library_id, **audit))
        origins = pd.concat(frames, ignore_index=True)
        origins = origins.sort_values(['chrom', 'strand', 'cut_site', 'library_id', 'site_id'])
        labels, offset = [], 0
        for _, part in origins.groupby(['chrom', 'strand'], sort=False):
            value = clusters(part.cut_site.to_numpy()) + offset
            labels.extend(value)
            offset = int(value[-1]) + 1
        origins['atlas_column_0based'] = labels
        f = origins.groupby('atlas_column_0based', sort=True).agg(
            chrom=('chrom', 'first'), strand=('strand', 'first'),
            member_cut_site_min=('cut_site', 'min'), member_cut_site_max=('cut_site', 'max'),
            n_libraries=('library_id', 'nunique'), native_peaks=('site_id', 'size'),
            current_target_umi=('native_molecules', 'sum')).reset_index()
        f['cut_site'] = (f.member_cut_site_min + f.member_cut_site_max) // 2
        f['cluster_width_bp'] = f.member_cut_site_max - f.member_cut_site_min
        f['atlas_pas_id'] = [f'polyApipe:{study}:{c}:{s}:{p}' for c, s, p in zip(f.chrom, f.strand, f.cut_site)]
        f['feature_index_1based'] = f.atlas_column_0based + 1
        f['species'] = lib.species.iloc[0]
        f['genome_build'] = 'GRCh38_2024_A' if lib.species.iloc[0] == 'human' else 'GRCm39_2024_A'
        f['atlas_status'] = np.where(f.n_libraries >= 2, 'core', 'single_library')
        assert f.atlas_pas_id.is_unique and f.cluster_width_bp.le(50).all()
        blocks = []
        for row, c in zip(lib.itertuples(index=False), matrices):
            mapping = origins[origins.library_id.eq(row.library_id)].sort_values('native_column_0based')
            assert np.array_equal(mapping.native_column_0based, np.arange(c.shape[1]))
            projection = sparse.csr_matrix((np.ones(c.shape[1], dtype=np.int64),
                (np.arange(c.shape[1]), mapping.atlas_column_0based)), shape=(c.shape[1], len(f)))
            block = c @ projection
            assert int(block.sum()) == int(c.sum())
            blocks.append(block)
        counts = sparse.vstack(blocks, format='csr')
        meta = pd.concat(metas, ignore_index=True)
        original_meta = pd.read_csv(source / 'comparison_inputs' / study / 'all_metacells.tsv.gz', sep='\t')
        assert set(original_meta.metacell_id) == set(meta.metacell_id)
        order = pd.Index(meta.metacell_id).get_indexer(original_meta.metacell_id)
        counts = counts[order].tocsr()
        assert np.array_equal(np.asarray(counts.sum(0)).ravel(), f.current_target_umi)
        assert len(meta) == lib.metacells.sum() and meta.n_cells.sum() == lib.cells.sum()
        save(f, dest / 'coordinate_features.tsv.gz')
        save(origins, dest / 'native_peak_to_atlas.tsv.gz')
        save(original_meta, dest / 'metacells.tsv.gz')
        sparse.save_npz(dest / 'metacell_by_pas_counts.npz', counts)
        reread = sparse.load_npz(dest / 'metacell_by_pas_counts.npz')
        assert (reread != counts).nnz == 0
        for name in ['all_metacells.tsv.gz', 'selected_metacells.tsv.gz', 'endpoint_library_design.tsv', 'sampling_design.tsv']:
            p = source / 'comparison_inputs' / study / name
            inputs[str(p)] = digest(p)
            (dest / name).write_bytes(p.read_bytes())
        save(pd.DataFrame(audits), dest / 'native_library_audit.tsv')
        record = dict(study_id=study, species=lib.species.iloc[0], libraries=len(lib),
                      cells=int(lib.cells.sum()), metacells=len(meta), native_peaks=len(origins),
                      atlas_pas=len(f), molecules=int(counts.sum()), count_conservation=True)
        summaries.append(record)
        print(json.dumps(record), flush=True)
    save(pd.DataFrame(summaries), out / 'atlas_summary.tsv')
    manifest = dict(status='native_atlas_prepared_not_full_rerun', caller='polyApipe',
        caller_commit='43d866e15b22e80974b827a871c3c4e1d7266677',
        cut_coordinate_system='0-based-interbase', scope='161 fixed GEX atlas libraries',
        native_filtering='native polyApipe; no ASNR/junction postfilter',
        clustering=dict(adjacent_bp=25, maximum_width_bp=50, representative='integer midpoint'),
        counts='sum native UR-deduplicated peak counts; no new read-level deduplication',
        inputs=inputs, studies=summaries)
    (out / 'intake.native.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--frozen', type=Path, required=True)
    a = p.parse_args()
    build(a.source.resolve(), a.out.resolve(), a.frozen.resolve())
